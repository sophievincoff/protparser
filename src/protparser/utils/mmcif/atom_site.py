# Utilities to extract resolved residue positions and non-canonical residue tags
# from the _atom_site category of an mmCIF file.
#
# get_resolved_positions_from_atom_site: Returns sequence positions with ATOM records for a given chain.
# _scan_chain_resnames_from_atom_site: Collects all residue names observed for a chain in ATOM records.
# _extract_position_tags_from_atom_site_noncanon: Identifies non-canonical residues and tags their positions per chain.


from __future__ import annotations
from pathlib import Path
from typing import Dict, Set
from collections import defaultdict
import gemmi

from .cif_columns import safe_find_column, cif_get_row_count, clean_str, safe_int
from .seq_cleaning import ALLOWED_POLYMER_RESNAMES

def get_resolved_positions_from_atom_site(cif_path: Path, chain_id: str) -> set[int]:
    doc = gemmi.cif.read_file(str(cif_path))
    block = doc.sole_block()
    atom = block.find_mmcif_category("_atom_site")
    if atom is None:
        return set()

    col_group = safe_find_column(atom, "group_PDB")
    col_asym  = safe_find_column(atom, "label_asym_id")
    col_seq   = safe_find_column(atom, "label_seq_id")

    if col_asym is None or col_seq is None:
        return set()

    resolved = set()
    n = len(atom) if hasattr(atom, "__len__") else atom.get_row_count()
    for i in range(n):
        if col_group is not None:
            g = str(col_group[i]).strip()
            if g != "ATOM":
                continue

        if str(col_asym[i]).strip() != str(chain_id).strip():
            continue

        s = str(col_seq[i]).strip()
        if not s or s in {".", "?"}:
            continue
        try:
            resolved.add(int(s))
        except ValueError:
            pass

    return resolved

def _scan_chain_resnames_from_atom_site(block: gemmi.cif.Block, chain_id_label: str) -> Set[str]:
    """
    Return set of 3-letter residue names (label_comp_id) observed for this chain
    in ATOM records. Uses label_asym_id matching.
    """
    atom = block.find_mmcif_category("_atom_site")
    if atom is None:
        return set()

    col_group = safe_find_column(atom, "group_PDB")       # optional
    col_asym  = safe_find_column(atom, "label_asym_id")   # required-ish
    col_comp  = safe_find_column(atom, "label_comp_id")   # required-ish

    if col_asym is None or col_comp is None:
        return set()

    out: Set[str] = set()
    n = len(atom) if hasattr(atom, "__len__") else atom.get_row_count()
    chain_id_label = str(chain_id_label).strip()

    for i in range(n):
        if col_group is not None:
            g = str(col_group[i]).strip()
            if g != "ATOM":
                continue

        if str(col_asym[i]).strip() != chain_id_label:
            continue

        comp = str(col_comp[i]).strip().upper()
        if not comp or comp in {".", "?"}:
            continue

        out.add(comp)

    return out

def _extract_position_tags_from_atom_site_noncanon(
    block: gemmi.cif.Block,
    *,
    allowed_resnames: Optional[Set[str]] = None,
) -> Dict[str, Dict[int, Set[str]]]:
    """
    From _atom_site, find residues that are non-natural (label_comp_id not in ALLOWED_POLYMER_RESNAMES)
    and tag that position with the observed residue name (e.g. PTR, SEP, MSE).
    Output format matches modification_feature extractor:
      {chain: {pos: set(tags)}}
    """
    allowed = allowed_resnames if allowed_resnames is not None else ALLOWED_POLYMER_RESNAMES
    
    atom = block.find_mmcif_category("_atom_site")
    if atom is None:
        return {}

    col_group = safe_find_column(atom, "group_PDB")
    col_chain = safe_find_column(atom, "label_asym_id")
    col_seq   = safe_find_column(atom, "label_seq_id")
    col_comp  = safe_find_column(atom, "label_comp_id")

    if col_chain is None or col_seq is None or col_comp is None:
        return {}

    out: Dict[str, Dict[int, Set[str]]] = defaultdict(lambda: defaultdict(set))
    n = cif_get_row_count(atom)

    for i in range(n):
        if col_group is not None:
            g = clean_str(col_group[i])
            if g and g != "ATOM":
                continue

        chain = clean_str(col_chain[i])
        pos = safe_int(col_seq[i])
        comp = clean_str(col_comp[i]).upper()

        if not chain or pos is None:
            continue

        # Only tag "non-natural" polymer residue names
        if comp and (comp not in allowed):
            out[chain][pos].add(comp)

    return out