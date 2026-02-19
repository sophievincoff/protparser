# protparser/parallel_inventory.py
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Sequence, Set, Tuple, Literal, Dict, Any

import concurrent.futures as cf
import pandas as pd
import traceback

from protparser.rcsb import RCSBStructure

PDBID = str
ExecutorKind = Literal["thread", "process"]


@dataclass
class InventoryJobResult:
    structures: List[RCSBStructure]

    # Combined inventories (all PDBs)
    chain_inventory_df: pd.DataFrame
    entity_inventory_df: pd.DataFrame
    entry_inventory_df: pd.DataFrame

    failed: List[Tuple[PDBID, str]]  # (pdb_id, error_string)


def build_inventories(
    pdb_ids: Sequence[str],
    *,
    cache_dir: str | Path,
    include_rcsb_rest: bool = True,
    valid_extra_tokens: Optional[Set[str]] = None,
    parallel: bool = True,
    n_workers: Optional[int] = None,
    executor: ExecutorKind = "thread",
) -> InventoryJobResult:
    """
    Orchestrator:
      - create RCSBStructure objects (optionally in parallel)
      - build chain, entity, and entry inventories on each
      - concat each inventory type across all PDB IDs

    Notes:
      - Default executor="thread" to be safe to call from anywhere (no __main__ guard).
      - If you choose executor="process", you may need a __main__ guard in some contexts.
    """
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)

    pdb_ids_norm = _normalize_pdb_ids(pdb_ids)
    if not pdb_ids_norm:
        return InventoryJobResult(
            structures=[],
            chain_inventory_df=pd.DataFrame(),
            entity_inventory_df=pd.DataFrame(),
            entry_inventory_df=pd.DataFrame(),
            failed=[],
        )

    if n_workers is None:
        # threads: allow higher; processes: keep modest
        n_workers = min(32, len(pdb_ids_norm)) if executor == "thread" else min(8, len(pdb_ids_norm))

    Executor = cf.ThreadPoolExecutor if executor == "thread" else cf.ProcessPoolExecutor

    structures: List[RCSBStructure] = []
    failures: List[Tuple[str, str]] = []

    def _build_one(pid: str) -> Dict[str, Any]:
        """
        Return dict so we can safely aggregate without relying on side effects.
        We still include the RCSBStructure object for report writing.
        """
        try:
            s = RCSBStructure(
                pid,
                cache_dir=cache_dir,
                download=True,
                valid_extra_tokens=valid_extra_tokens,
            )

            # 1) CHAIN inventory (also usually prepares other cached state)
            s.build_inventory(include_rcsb_rest=include_rcsb_rest)
            chain_df = getattr(s, "inventory_df", None)
            if chain_df is None:
                raise RuntimeError("RCSBStructure.build_inventory did not set inventory_df")

            # 2) ENTITY inventory
            # Make sure you have this method on your class (you do).
            ent_df = s.build_entity_inventory(include_rcsb_rest=include_rcsb_rest)
            if ent_df is None:
                raise RuntimeError("RCSBStructure.build_entity_inventory returned None")

            # 3) ENTRY inventory
            entry_df = s.build_entry_inventory()
            if entry_df is None:
                raise RuntimeError("RCSBStructure.build_entry_inventory returned None")

            return {
                "ok": True,
                "pdb_id": pid,
                "structure": s,
                "chain_df": chain_df,
                "entity_df": ent_df,
                "entry_df": entry_df,
            }
        except Exception as e:
            return {
                "ok": False,
                "pdb_id": pid,
                "error": f"{type(e).__name__}: {e}\n{traceback.format_exc()}",
            }

    results: List[Dict[str, Any]] = []
    if parallel and n_workers > 1:
        if executor == "process":
            raise RuntimeError(
                "executor='process' is not supported with this 'call anywhere' version. "
                "Use executor='thread' (default) or I can give you a process-safe variant."
            )

        with Executor(max_workers=n_workers) as ex:
            futs = [ex.submit(_build_one, pid) for pid in pdb_ids_norm]
            for f in cf.as_completed(futs):
                results.append(f.result())
    else:
        for pid in pdb_ids_norm:
            results.append(_build_one(pid))

    # Aggregate
    chain_dfs: List[pd.DataFrame] = []
    entity_dfs: List[pd.DataFrame] = []
    entry_dfs: List[pd.DataFrame] = []

    for r in results:
        if r.get("ok"):
            structures.append(r["structure"])

            cdf = r.get("chain_df")
            if isinstance(cdf, pd.DataFrame) and not cdf.empty:
                chain_dfs.append(cdf)

            edf = r.get("entity_df")
            if isinstance(edf, pd.DataFrame) and not edf.empty:
                entity_dfs.append(edf)

            idf = r.get("entry_df")
            if isinstance(idf, pd.DataFrame) and not idf.empty:
                entry_dfs.append(idf)
        else:
            failures.append((r.get("pdb_id", ""), r.get("error", "unknown_error")))

    chain_inventory_df = pd.concat(chain_dfs, ignore_index=True, sort=False) if chain_dfs else pd.DataFrame()
    entity_inventory_df = pd.concat(entity_dfs, ignore_index=True, sort=False) if entity_dfs else pd.DataFrame()
    entry_inventory_df = pd.concat(entry_dfs, ignore_index=True, sort=False) if entry_dfs else pd.DataFrame()

    return InventoryJobResult(
        structures=structures,
        chain_inventory_df=chain_inventory_df,
        entity_inventory_df=entity_inventory_df,
        entry_inventory_df=entry_inventory_df,
        failed=failures,
    )


def combine_inventories(
    structures: Sequence[RCSBStructure],
) -> Dict[str, pd.DataFrame]:
    """
    Combine inventories from already-built RCSBStructure objects.

    Returns dict with keys: chain/entity/entry.
    """
    chain_dfs: List[pd.DataFrame] = []
    entity_dfs: List[pd.DataFrame] = []
    entry_dfs: List[pd.DataFrame] = []

    for s in structures:
        cdf = getattr(s, "inventory_df", None)
        if isinstance(cdf, pd.DataFrame) and not cdf.empty:
            chain_dfs.append(cdf)

        edf = getattr(s, "entity_inventory_df", None)
        if isinstance(edf, pd.DataFrame) and not edf.empty:
            entity_dfs.append(edf)

        idf = getattr(s, "entry_inventory_df", None)
        if isinstance(idf, pd.DataFrame) and not idf.empty:
            entry_dfs.append(idf)

    return {
        "chain": pd.concat(chain_dfs, ignore_index=True, sort=False) if chain_dfs else pd.DataFrame(),
        "entity": pd.concat(entity_dfs, ignore_index=True, sort=False) if entity_dfs else pd.DataFrame(),
        "entry": pd.concat(entry_dfs, ignore_index=True, sort=False) if entry_dfs else pd.DataFrame(),
    }


def write_inventory_reports(
    result: InventoryJobResult,
    *,
    mode: Literal["group", "per_pdb"] = "per_pdb",
    include_rcsb_rest: bool = True,
    output: Literal["string", "file"] = "file",
    out_path: Optional[str | Path] = None,
) -> str | Path | dict[str, Path]:
    """
    Write reports for an InventoryJobResult.

    mode:
      - "group": one combined report for all PDBs (requires out_path if output='file')
      - "per_pdb": one report per PDB ID (out_path treated as output directory if output='file')

    output:
      - "string": return text (group) or dict[pdb_id -> text] (per_pdb)
      - "file": return Path (group) or dict[pdb_id -> Path] (per_pdb)
    """
    structures = result.structures or []
    if not structures:
        empty = "No structures available in result.structures.\n"
        if output == "string":
            return empty if mode == "group" else {}
        if out_path is None:
            raise ValueError("out_path required for output='file'")
        p = Path(out_path)
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(empty, encoding="utf-8")
        return p

    if mode == "group":
        text_parts: list[str] = []
        for s in sorted(structures, key=lambda x: x.pdb_id):
            text_parts.append(
                s.generate_report(include_rcsb_rest=include_rcsb_rest, output="string")
            )
            text_parts.append("\n" + ("-" * 60) + "\n")
        group_text = "".join(text_parts).rstrip() + "\n"

        if output == "string":
            return group_text

        if out_path is None:
            raise ValueError("out_path required for output='file' in mode='group'")
        out_path = Path(out_path)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(group_text, encoding="utf-8")
        return out_path

    # mode == "per_pdb"
    if output == "string":
        out: dict[str, str] = {}
        for s in structures:
            out[s.pdb_id.lower()] = s.generate_report(include_rcsb_rest=include_rcsb_rest, output="string")
        return out

    # output == "file"
    if out_path is None:
        raise ValueError("out_path required for output='file' in mode='per_pdb' (directory)")
    out_dir = Path(out_path)
    out_dir.mkdir(parents=True, exist_ok=True)

    written: dict[str, Path] = {}
    for s in structures:
        pid = s.pdb_id.upper()
        p = out_dir / f"{pid}_report.txt"
        s.generate_report(include_rcsb_rest=include_rcsb_rest, output="file", out_path=p)
        written[s.pdb_id.lower()] = p
    return written


def _normalize_pdb_ids(pdb_ids: Sequence[str]) -> List[str]:
    seen = set()
    out: List[str] = []
    for pid in pdb_ids:
        if pid is None:
            continue
        p = str(pid).strip().lower()
        if not p or p in seen:
            continue
        seen.add(p)
        out.append(p)
    return out