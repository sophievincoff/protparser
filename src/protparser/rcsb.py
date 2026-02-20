from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Literal  # <-- CHANGED (added Literal)
import json
import re
import os

import numpy as np
import math
import pandas as pd
import requests
import gemmi

# --- Your existing mmCIF utilities ---
from protparser.utils.mmcif.cif_columns import (
    safe_find_column,
    col_required,
    cif_get_row_count,
    safe_int,
    clean_str,
)

from protparser.utils.mmcif.seq_cleaning import (
    clean_mmcif_seq_text,
    _seq_is_clean,
    expand_allowed_tokens
)

from protparser.utils.mmcif.seq_build import (
    build_resolved_sequences,
    build_resolved_index_map,
    map_canon_pos_to_resolved_only_indices,
    _positions_to_csv,
)

from protparser.utils.mmcif.atom_site import (
    get_resolved_positions_from_atom_site,
    _extract_position_tags_from_atom_site_noncanon,
)

from protparser.utils.mmcif.modification_feature import (
    entry_has_protein_modification_flag,
    extract_mod_features_by_chain,
    extract_position_tags_from_modification_feature,
    merge_pos_tags,
    _compute_chain_flags,
    annotate_seq_with_position_tags,
    annotate_masked_seq_with_position_tags,
    _parse_yes_no_flag,
    build_entity_source_lookup
)

from protparser.utils.mmcif.seq_cleaning import (
    VALID_AAS, ALLOWED_POLYMER_RESNAMES
)
from protparser.utils.pandas_utils import harmonize_nulls_to_nan

from protparser.utils.mmcif.entity import build_entity_lookup

from protparser.utils.data_utils import (
    _norm_str, 
    _csv_sorted_unique,
    _pipe_sorted_unique
)

class RCSBStructure:
    """
    Single PDB entry object.

    Responsibilities:
      1) Ensure a cached mmCIF exists (reuse silently if present).
      2) Build a full chain inventory (one row per label_asym_id chain) from CIF.
      3) Build a PDB-entry inventory (one row per PDB) from CIF (structure-level flags).
      4) Optional REST enrichment into the chain inventory.
      5) Add per-PDB binary interaction flags (computed from chain inventory).

    Public usage:
      r = RCSBStructure("1a2b", cache_dir="rcsb_cache")
      chain_df = r.build_inventory(include_rcsb_rest=True)
      entry_df = r.build_entry_inventory()
    """
    
    _FASTA_TAXID_RE = re.compile(r"^(?P<species>.*?)(?:\s*\((?P<taxid>\d+)\))?\s*$")
    _CHAIN_PREFIX_RE = re.compile(r"^\s*Chains?\s+", flags=re.IGNORECASE)
    _CHAIN_TOKEN_RE = re.compile(
        r"""^\s*
            (?P<pdb>[A-Za-z0-9]+)
            \s*
            (?:\[\s*auth\s*[:=]?\s*(?P<auth>[A-Za-z0-9]+)\s*\])?
            \s*$""",
        flags=re.IGNORECASE | re.VERBOSE,
    )

    def __init__(
        self,
        pdb_id: str,
        *,
        cache_dir: Optional[str | Path] = None,
        cif_path: Optional[str | Path] = None,
        download: bool = True,
        valid_extra_tokens: Optional[Set[str]] = None,
    ):
        self.pdb_id = str(pdb_id).strip().lower()
        # Set cache dir and create it if it doesn't exist
        self.cache_dir = Path(cache_dir) if cache_dir is not None else None
        if self.cache_dir is not None:
            self.cache_dir.mkdir(parents=True, exist_ok=True)

        self.cif_path: Optional[Path] = Path(cif_path) if cif_path is not None else None

        # inventory_df is explicitly "chain-level inventory"
        self.inventory_df: Optional[pd.DataFrame] = None

        # entry-level inventory (one row per pdb_id)
        self.entry_inventory_df: Optional[pd.DataFrame] = None
        
        # entity-level inventory (one row per entity)
        self.entity_inventory_df: Optional[pd.DataFrame] = None

        if self.cif_path is None and download:
            self.cif_path = self._ensure_cif()

        extra_valid, extra_res = expand_allowed_tokens(valid_extra_tokens)
        self.valid_aas = VALID_AAS | extra_valid
        self.allowed_resnames = ALLOWED_POLYMER_RESNAMES | extra_res
        
        # storing the FASTA 
        self.fasta_text: Optional[str] = None
        self.fasta_records_df: Optional[pd.DataFrame] = None

    # -------------------------
    # Public API
    # -------------------------
    def build_inventory(self, *, include_rcsb_rest: bool = True) -> pd.DataFrame:
        """
        Build and store the CHAIN-LEVEL inventory for this PDB entry.

        Also ensures entry_inventory_df exists (built from CIF) so structure-level
        annotations are available separately.
        """
        cif_path = self.cif_path or self._ensure_cif()

        # Build chain-level inventory
        inv = self._build_inventory_from_cif(cif_path)

        if include_rcsb_rest:
            rest_df = self._build_rest_chain_table()
            if rest_df is not None and not rest_df.empty:
                inv = inv.merge(
                    rest_df,
                    on=["pdb_id", "chain_id", "entity_id"],
                    how="left",
                    suffixes=("", "_rest"),
                )
                
        if "chain_id_auth_rest" in inv.columns:
            inv["chain_id_auth_cif"] = inv["chain_id_auth"]  # keep original
            m = inv["chain_id_auth_rest"].notna() & inv["chain_id_auth_rest"].astype(str).str.strip().ne("")
            inv.loc[m, "chain_id_auth"] = inv.loc[m, "chain_id_auth_rest"]
        else:
            inv = inv.drop(columns=["chain_id_auth"])

        # compute binary info for entry inventory, but do NOT keep these cols on chain inventory
        inv_with_binary = self._annotate_binary_interactions(inv)

        # drop entry-level binary cols from chain inventory
        binary_cols = ["binary_n_unique_protein_seqs", "binary_stoichiometry", "is_binary_interaction"]
        inv = inv_with_binary.drop(columns=[c for c in binary_cols if c in inv_with_binary.columns], errors="ignore")

        inv = harmonize_nulls_to_nan(inv)
        self.inventory_df = inv
        
        if self.entity_inventory_df is None:
            self.entity_inventory_df = self._build_entity_inventory_from_chain_inventory(self.inventory_df)

        # ensure entry inventory uses the version that *does* have binary cols available
        if self.entry_inventory_df is None:
            self.entry_inventory_df = self._build_entry_inventory_from_cif(cif_path, chain_inv=inv_with_binary)

        return inv

    def build_entry_inventory(self) -> pd.DataFrame:
        """
        Build and store the ENTRY-LEVEL inventory (one row per pdb_id).
        This is where structure-level flags belong.
        """
        cif_path = self.cif_path or self._ensure_cif()
            
        chain_inv = self.inventory_df
        if chain_inv is None:
            chain_inv = self.build_inventory(include_rcsb_rest=False)

        # --- some redundancy: ensure binary cols exist for entry aggregation ---
        chain_inv = self._annotate_binary_interactions(chain_inv)

        self.entry_inventory_df = self._build_entry_inventory_from_cif(cif_path, chain_inv=chain_inv)
        return self.entry_inventory_df

    def get_inventory_df(self) -> pd.DataFrame:
        """Back-compat name: returns CHAIN-LEVEL inventory."""
        if self.inventory_df is None:
            raise RuntimeError("inventory_df is not built yet. Call build_inventory() first.")
        return self.inventory_df

    def get_entry_inventory_df(self) -> pd.DataFrame:
        if self.entry_inventory_df is None:
            raise RuntimeError("entry_inventory_df is not built yet. Call build_entry_inventory() first.")
        return self.entry_inventory_df
    
    def load_fasta(self, *, force: bool = False) -> str:
        """
        Ensure FASTA is downloaded and loaded into self.fasta_text.
        Returns the FASTA text.
        """
        if force:
            # delete cached and re-download
            out_dir = self.cache_dir if self.cache_dir is not None else Path(".")
            target = out_dir / f"{self.pdb_id.upper()}.fasta"
            if target.exists():
                target.unlink()

        fasta_path = self._ensure_fasta()
        self.fasta_text = fasta_path.read_text(encoding="utf-8")
        return self.fasta_text

    def write_fasta(self, out_path: str | Path) -> Path:
        """
        Write the downloaded FASTA to a file.
        """
        if self.fasta_text is None:
            self.load_fasta()
        out_path = Path(out_path)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(self.fasta_text or "", encoding="utf-8")
        return out_path
    
    # -------------------------
    # CIF management (silent reuse)
    # -------------------------
    def _ensure_cif(self) -> Path:
        """
        Ensure CIF exists in cache_dir. If it already exists, reuse silently.
        """
        out_dir = self.cache_dir if self.cache_dir is not None else Path(".")
        out_dir.mkdir(parents=True, exist_ok=True)

        target = out_dir / f"{self.pdb_id.upper()}.cif"
        if target.exists():
            return target

        p = download_rcsb(self.pdb_id, struct_format="cif", output_dir=str(out_dir))
        return Path(p)
    
    # -------------------------
    # Ensure fasta exists
    # -------------------------
    def _ensure_fasta(self) -> Path:
        out_dir = self.cache_dir if self.cache_dir is not None else Path(".")
        out_dir.mkdir(parents=True, exist_ok=True)

        target = out_dir / f"{self.pdb_id.upper()}.fasta"
        if target.exists():
            return target

        url = f"https://www.rcsb.org/fasta/entry/{self.pdb_id.upper()}"
        r = requests.get(url, timeout=30)
        if r.status_code != 200:
            raise RuntimeError(f"GET {url} failed: {r.status_code}")
        target.write_text(r.text, encoding="utf-8")
        return target

    # -------------------------
    # Chain inventory from CIF
    # -------------------------
    def _build_inventory_from_cif(self, cif_path: Path) -> pd.DataFrame:
        doc = gemmi.cif.read_file(str(cif_path))
        block = doc.sole_block()
        entity_lookup = build_entity_lookup(block)
        
        entity_source_lookup = build_entity_source_lookup(block)

        asym = block.find_mmcif_category("_struct_asym")
        if asym is None:
            raise RuntimeError("Missing _struct_asym category")

        chain_id_col = col_required(asym, "id")          # label_asym_id
        entity_id_col = col_required(asym, "entity_id")
        auth_chain_col = safe_find_column(asym, "pdbx_auth_asym_id")  # optional

        ent_poly = block.find_mmcif_category("_entity_poly")
        if ent_poly is None:
            raise RuntimeError("Missing _entity_poly category")

        ent_poly_ent = col_required(ent_poly, "entity_id")
        seq_can_col  = col_required(ent_poly, "pdbx_seq_one_letter_code_can")
        poly_type_col = safe_find_column(ent_poly, "type")  # optional

        # entity_id -> canonical polymer sequence + type
        entity_to_seq: Dict[str, str] = {}
        entity_to_type: Dict[str, str] = {}

        for i in range(cif_get_row_count(ent_poly_ent)):
            eid = clean_str(ent_poly_ent[i])
            seq = clean_mmcif_seq_text(seq_can_col[i])
            entity_to_seq[eid] = seq
            if poly_type_col is not None:
                entity_to_type[eid] = clean_str(poly_type_col[i])

        # CIF-wide PTM feature summaries + position tags
        feat_by_chain = extract_mod_features_by_chain(block)

        mf_tags = extract_position_tags_from_modification_feature(block)
        atom_noncanon_tags = _extract_position_tags_from_atom_site_noncanon(
            block,
            allowed_resnames=self.allowed_resnames,
        )
        pos_tags_by_chain = merge_pos_tags(mf_tags, atom_noncanon_tags)

        rows: List[Dict[str, Any]] = []
        pdb_id = self.pdb_id

        for i in range(cif_get_row_count(chain_id_col)):
            chain_id = clean_str(chain_id_col[i])     # label chain
            entity_id = clean_str(entity_id_col[i])
            
            src = entity_source_lookup.get(entity_id, {}) or {}
            species_name = src.get("species_name", np.nan)
            taxid = src.get("taxid", np.nan)
                        
            entity_meta = entity_lookup.get(entity_id, {})

            chain_id_auth = chain_id
            if auth_chain_col is not None:
                v = clean_str(auth_chain_col[i])
                if v:
                    chain_id_auth = v

            seq = entity_to_seq.get(entity_id, "")
            ok, why = _seq_is_clean(seq, allowed_aas=self.valid_aas) if seq else (False, "no_sequence")

            invalid_letters = ""
            if not ok:
                s = re.sub(r"\s+", "", seq).upper() if seq else ""
                bad = sorted(set(s) - set("ACDEFGHIKLMNPQRSTVWYU"))
                invalid_letters = "".join(bad) if bad else (why or "")

            # resolved positions + resolved sequences
            resolved_pos = get_resolved_positions_from_atom_site(cif_path, chain_id)
            seq_res_masked, seq_res_only, n_res, frac_res = build_resolved_sequences(seq, resolved_pos)

            # PTM / noncanon position sets (canonical indexing)
            ptm_pos_can = set(mf_tags.get(chain_id, {}).keys()) if mf_tags.get(chain_id) else set()
            noncanon_pos_can = set(atom_noncanon_tags.get(chain_id, {}).keys()) if atom_noncanon_tags.get(chain_id) else set()

            canon_to_res = build_resolved_index_map(len(seq), resolved_pos)
            ptm_pos_res = map_canon_pos_to_resolved_only_indices(ptm_pos_can, canon_to_res)
            noncanon_pos_res = map_canon_pos_to_resolved_only_indices(noncanon_pos_can, canon_to_res)

            # token-annotated sequences
            chain_tags = pos_tags_by_chain.get(chain_id, {})
            seq_can_annot = annotate_seq_with_position_tags(seq, chain_tags)
            seq_res_masked_annot = annotate_masked_seq_with_position_tags(seq_res_masked, seq, chain_tags)

            # chain flags (PTM/noncanonical + reasons)
            flags = _compute_chain_flags(block, chain_id_label=chain_id, chain_id_auth=chain_id_auth)

            # detailed mod feature summary
            feat = feat_by_chain.get(chain_id, None)
            ptm_categories = sorted(feat["categories"]) if feat is not None else []
            ptm_categories = [clean_str(x) for x in ptm_categories]
            ptm_types      = sorted(feat["types"]) if feat is not None else []
            ptm_comp_ids   = sorted(feat["comp_ids"]) if feat is not None else []
            ptm_examples   = feat["examples"] if feat is not None else []
            ptm_n_features = int(feat["n"]) if feat is not None else 0

            rows.append({
                "pdb_id": pdb_id,
                "chain_id": chain_id,
                "chain_id_auth": chain_id_auth,
                "entity_id": entity_id,
    
                
                # --- entity-level fields
                "entity_type": entity_meta.get("entity_type", ""),
                "entity_description": entity_meta.get("entity_description", ""),
                "entity_src_method": entity_meta.get("entity_src_method", ""),
                #"entity_formula_weight": entity_meta.get("entity_formula_weight", np.nan),
                "entity_number_of_molecules": entity_meta.get("entity_number_of_molecules", np.nan),
                "entity_mutation": entity_meta.get("entity_mutation", np.nan),
                
                "species_name": species_name,
                "species_taxid": taxid,

                # optional extras
                "entity_ec": entity_meta.get("entity_ec", np.nan),
                "entity_fragment": entity_meta.get("entity_fragment", np.nan),
                "entity_details": entity_meta.get("entity_details", np.nan),

                # if we care about ligands/water names when chains are non-polymers:
                "nonpoly_name": entity_meta.get("nonpoly_name", np.nan),
                #"comp_id": entity_meta.get("comp_id", np.nan),
                "chem_comp_name": entity_meta.get("chem_comp_name", np.nan),
                #"chem_comp_formula": entity_meta.get("chem_comp_formula", np.nan),
                #"chem_comp_formula_weight": entity_meta.get("chem_comp_formula_weight", np.nan),

                "seq_can": seq,
                "seq_len": len(seq) if seq else 0,
                "is_valid_aa": bool(ok),
                "invalid_letters": invalid_letters,
                "entity_poly_type": entity_to_type.get(entity_id, ""),

                # PTM / noncanonical flags (CHAIN-LEVEL ONLY)
                "has_ptm": flags["has_ptm"],
                "has_ptm_reason": flags["has_ptm_reason"],
                "ptm_detect_sources": flags["ptm_detect_sources"],
                "ptm_resnames": flags["ptm_resnames"],
                "ptm_mf_chain_hit": flags["ptm_mf_chain_hit"],

                "has_noncanonical_resname": flags["has_noncanonical_resname"],
                "noncanonical_resnames": flags["noncanonical_resnames"],

                # detailed PTM features (CSV-friendly)
                "ptm_categories": "|".join(ptm_categories) if ptm_categories else np.nan,
                "ptm_types": "|".join(ptm_types) if ptm_types else np.nan,
                "ptm_comp_ids": "|".join(ptm_comp_ids) if ptm_comp_ids else np.nan,
                "ptm_n_features": ptm_n_features,
                "ptm_examples_json": self._json_or_nan(ptm_examples),

                # resolved sequence info
                "seq_resolved_masked": seq_res_masked,
                "seq_resolved_only": seq_res_only,
                "n_resolved": n_res,
                "frac_resolved": frac_res,

                # annotated sequences
                "seq_can_annot": seq_can_annot,
                "seq_resolved_masked_annot": seq_res_masked_annot,

                # positions
                "ptm_pos_can": _positions_to_csv(ptm_pos_can),
                "noncanon_pos_can": _positions_to_csv(noncanon_pos_can),
                "ptm_pos_resolved_only": _positions_to_csv(ptm_pos_res),
                "noncanon_pos_resolved_only": _positions_to_csv(noncanon_pos_res),
            })

        return pd.DataFrame(rows)

    # -------------------------
    # Entry-level inventory (structure-level flags live here)
    # -------------------------
    def _build_entry_inventory_from_cif(self, cif_path: Path, *, chain_inv: Optional[pd.DataFrame] = None) -> pd.DataFrame:
        """
        Build a one-row DataFrame with ENTRY/STRUCTURE-level annotations.

        We also include a few derived rollups computed from the chain inventory
        (e.g., whether any chain has PTM/noncanon evidence, and which chains).
        """
        doc = gemmi.cif.read_file(str(cif_path))
        block = doc.sole_block()
        
        # --- PDB title ---
        pdb_title = np.nan
        st = block.find_mmcif_category("_struct")
        if st is not None:
            col_title = safe_find_column(st, "title")
            if col_title is not None:
                try:
                    t = clean_str(col_title[0])
                    pdb_title = t if t not in ("", "?", ".") else np.nan
                except Exception:
                    pdb_title = np.nan

        entry_flag = entry_has_protein_modification_flag(block)  # likely "Y"/"N"/None

        # Normalize to bool/NA for cleaner downstream usage
        if entry_flag is None or (isinstance(entry_flag, str) and not entry_flag.strip()):
            entry_has_ptm = np.nan
        else:
            if isinstance(entry_flag, bool):
                entry_has_ptm = entry_flag
            else:
                entry_has_ptm = _parse_yes_no_flag(entry_flag)

        # If chain inventory is not provided, use what we have (or build it)
        if chain_inv is None:
            chain_inv = self.inventory_df

        if chain_inv is None:
            # last resort: build chain inventory (no REST)
            chain_inv = self.build_inventory(include_rcsb_rest=False)

        c = chain_inv.copy() if chain_inv is not None else pd.DataFrame()
        # --- ensure binary columns exist (compute them if missing) ---
        binary_cols = {"binary_n_unique_protein_seqs", "binary_stoichiometry", "is_binary_interaction"}
        if not c.empty and not binary_cols.issubset(set(c.columns)):
            c = self._annotate_binary_interactions(c)
        if not c.empty:
            # Determine protein chains count (same logic as your binary detection)
            is_protein_chain = c.get("entity_poly_type", "").astype(str).str.lower().str.startswith("polypeptide")
            
            # --- basic normalized fields ---
            ptype = c.get("entity_poly_type", "").astype(str).str.strip().str.lower()
            edesc = c.get("entity_description", "").astype(str).str.strip().str.lower()

            chain_id = c.get("chain_id", pd.Series([], dtype=str)).astype(str).str.strip()
            chain_auth = c.get("chain_id_auth", pd.Series([], dtype=str)).astype(str).str.strip()

            # --- masks ---
            is_protein = ptype.str.startswith("polypeptide")
            is_dna = ptype.eq("polydeoxyribonucleotide")
            is_rna = ptype.eq("polyribonucleotide")

            # Water is a nonpoly entity description of "water"
            is_water = edesc.eq("water")

            # Nonpolymer other: NOT polymer-type present (or ptype empty/"nan") AND NOT water
            # (Your chain table labels non-polymers by ptype == "nan" after stringify; handle robustly)
            is_polymer = ptype.ne("") & ~ptype.eq("nan")
            is_nonpoly_other = (~is_polymer) & (~is_water)
            
            # --- protein species rollups (unique across protein entities) ---
            species_lookup = build_entity_source_lookup(block)

            prot_entity_ids = (
                c.loc[is_protein, "entity_id"].dropna().astype(str).map(lambda x: x.strip()).tolist()
                if "entity_id" in c.columns else []
            )

            species_names = []
            species_taxids = []

            for eid in prot_entity_ids:
                info = species_lookup.get(eid, {})
                n = info.get("species_name", np.nan)
                tx = info.get("taxid", np.nan)
                if isinstance(n, str) and n.strip():
                    species_names.append(n.strip())
                if isinstance(tx, (int, np.integer)) and not (isinstance(tx, float) and np.isnan(tx)):
                    species_taxids.append(str(int(tx)))

            protein_species_names = "|".join(sorted(set(species_names))) if species_names else np.nan
            protein_species_taxids = "|".join(sorted(set(species_taxids))) if species_taxids else np.nan


            # --- chain lists (label and auth) ---
            protein_chains = _csv_sorted_unique(c.loc[is_protein, "chain_id"].astype(str).tolist()) if "chain_id" in c.columns else np.nan
            protein_chains_auth = _csv_sorted_unique(c.loc[is_protein, "chain_id_auth"].astype(str).tolist()) if "chain_id_auth" in c.columns else np.nan

            dna_chains = _csv_sorted_unique(c.loc[is_dna, "chain_id"].astype(str).tolist()) if "chain_id" in c.columns else np.nan
            dna_chains_auth = _csv_sorted_unique(c.loc[is_dna, "chain_id_auth"].astype(str).tolist()) if "chain_id_auth" in c.columns else np.nan

            rna_chains = _csv_sorted_unique(c.loc[is_rna, "chain_id"].astype(str).tolist()) if "chain_id" in c.columns else np.nan
            rna_chains_auth = _csv_sorted_unique(c.loc[is_rna, "chain_id_auth"].astype(str).tolist()) if "chain_id_auth" in c.columns else np.nan

            water_chains = _csv_sorted_unique(c.loc[is_water, "chain_id"].astype(str).tolist()) if "chain_id" in c.columns else np.nan
            water_chains_auth = _csv_sorted_unique(c.loc[is_water, "chain_id_auth"].astype(str).tolist()) if "chain_id_auth" in c.columns else np.nan

            nonpoly_other_chains = _csv_sorted_unique(c.loc[is_nonpoly_other, "chain_id"].astype(str).tolist()) if "chain_id" in c.columns else np.nan
            nonpoly_other_chains_auth = _csv_sorted_unique(c.loc[is_nonpoly_other, "chain_id_auth"].astype(str).tolist()) if "chain_id_auth" in c.columns else np.nan

            # --- PTM chains (auth) ---
            ptm_chains_auth = (
                _csv_sorted_unique(c.loc[c["has_ptm"].astype(bool), "chain_id_auth"].astype(str).tolist())
                if "has_ptm" in c.columns and "chain_id_auth" in c.columns
                else np.nan
            )

            # --- name lists ---
            # Use entity_description primarily; for nonpolymers you might prefer nonpoly_name/chem_comp_name if present.
            protein_names = _pipe_sorted_unique(c.loc[is_protein, "entity_description"].astype(str).tolist()) if "entity_description" in c.columns else np.nan
            dna_names = _pipe_sorted_unique(c.loc[is_dna, "entity_description"].astype(str).tolist()) if "entity_description" in c.columns else np.nan
            rna_names = _pipe_sorted_unique(c.loc[is_rna, "entity_description"].astype(str).tolist()) if "entity_description" in c.columns else np.nan

            # All nonpolymers (including water): prefer nonpoly_name, else chem_comp_name, else entity_description
            if any(col in c.columns for col in ["nonpoly_name", "chem_comp_name", "entity_description"]):
                nonpoly_name_source = (
                    c.get("nonpoly_name", pd.Series([np.nan] * len(c)))
                    .fillna(c.get("chem_comp_name", pd.Series([np.nan] * len(c))))
                    .fillna(c.get("entity_description", pd.Series([np.nan] * len(c))))
                )
                nonpolymer_names = _pipe_sorted_unique(nonpoly_name_source.loc[~is_polymer].astype(str).tolist())
            else:
                nonpolymer_names = np.nan


            n_chains_total = int(c["chain_id"].nunique()) if "chain_id" in c.columns else int(len(c))
            n_protein_chains = int(is_protein_chain.sum()) if "entity_poly_type" in c.columns else 0

            # rollups for chain-level PTM/noncanon evidence
            has_any_chain_ptm = bool(c.get("has_ptm", False).astype(bool).any()) if "has_ptm" in c.columns else False
            has_any_chain_noncanon = bool(c.get("has_noncanonical_resname", False).astype(bool).any()) if "has_noncanonical_resname" in c.columns else False

            ptm_chains = (
                ",".join(sorted(c.loc[c["has_ptm"].astype(bool), "chain_id"].astype(str).unique().tolist()))
                if "has_ptm" in c.columns and "chain_id" in c.columns
                else ""
            )
            noncanon_chains = (
                ",".join(sorted(c.loc[c["has_noncanonical_resname"].astype(bool), "chain_id"].astype(str).unique().tolist()))
                if "has_noncanonical_resname" in c.columns and "chain_id" in c.columns
                else ""
            )
            
            # --- invalid AA rollups (polypeptide chains only) ---
            is_polyL = ptype.eq("polypeptide(l)")  # only these count for invalid-aa summary

            if "is_valid_aa" in c.columns and "chain_id" in c.columns:
                invalid_mask = is_polyL & (~c["is_valid_aa"].astype(bool))
                has_any_chain_invalid_aa = bool(invalid_mask.any())

                invalid_aa_chains = _csv_sorted_unique(
                    c.loc[invalid_mask, "chain_id"].astype(str).tolist()
                )

                invalid_aa_chains_auth = _csv_sorted_unique(
                    c.loc[invalid_mask, "chain_id_auth"].astype(str).tolist()
                ) if "chain_id_auth" in c.columns else np.nan
            else:
                has_any_chain_invalid_aa = False
                invalid_aa_chains = np.nan
                invalid_aa_chains_auth = np.nan

            
            # --- protein UniProt IDs (pipe-joined unique across protein chains) ---
            if "uniprot_ids" in c.columns:
                prot_uniprot_series = c.loc[is_protein, "uniprot_ids"].dropna().astype(str)

                all_ids = []
                for val in prot_uniprot_series:
                    parts = [p.strip() for p in val.split(",") if p.strip()]
                    all_ids.extend(parts)

                protein_uniprot_ids = (
                    "|".join(sorted(set(all_ids))) if all_ids else np.nan
                )
            else:
                protein_uniprot_ids = np.nan

            # binary interaction fields (if present on chain table)
            is_binary = bool(c["is_binary_interaction"].iloc[0]) if "is_binary_interaction" in c.columns and len(c) else False
            stoich = c["binary_stoichiometry"].iloc[0] if "binary_stoichiometry" in c.columns and len(c) else np.nan
            # get reduced_stoich
            reduced_stoich = np.nan
            if isinstance(stoich, str):
                stoich_a, stoich_b = [int(x) for x in stoich.split(":")]
                stoich_gcd = math.gcd(stoich_a,stoich_b)
                stoich_a = int(stoich_a/stoich_gcd)
                stoich_b = int(stoich_b/stoich_gcd)
                reduced_stoich = f"{stoich_a}:{stoich_b}"
            nuniq = int(c["binary_n_unique_protein_seqs"].iloc[0]) if "binary_n_unique_protein_seqs" in c.columns and len(c) else 0
        else:
            pdb_title = ""
            n_chains_total = 0
            n_protein_chains = 0
            has_any_chain_ptm = False
            has_any_chain_noncanon = False
            ptm_chains = ""
            noncanon_chains = ""
            is_binary = False
            stoich = np.nan
            reduced_stoich = np.nan
            nuniq = 0
            protein_chains = protein_chains_auth = ptm_chains_auth = np.nan
            dna_chains = dna_chains_auth = np.nan
            rna_chains = rna_chains_auth = np.nan
            water_chains = water_chains_auth = np.nan
            nonpoly_other_chains = nonpoly_other_chains_auth = np.nan
            protein_names = dna_names = rna_names = nonpolymer_names = np.nan
            protein_uniprot_ids = np.nan
            protein_species_names = np.nan
            protein_species_taxids = np.nan
            has_any_chain_invalid_aa = False
            invalid_aa_chains = np.nan
            invalid_aa_chains_auth = np.nan


        row = {
            "pdb_id": self.pdb_id,
            "pdb_title": pdb_title,
            
            "protein_chains": protein_chains,
            "protein_chains_auth": protein_chains_auth,
            
            "protein_names": protein_names,
            "protein_uniprot_ids": protein_uniprot_ids,
            
            "protein_species_names": protein_species_names,
            "protein_species_taxids": protein_species_taxids,

            "dna_names": dna_names,
            "rna_names": rna_names,
            "nonpolymer_names": nonpolymer_names,
            
            "entry_has_protein_modification": entry_has_ptm,
            "entry_has_protein_modification_raw": entry_flag if entry_flag is not None else np.nan,

            "n_chains_total": n_chains_total,
            "n_protein_chains": n_protein_chains,

            "has_any_chain_ptm_evidence": has_any_chain_ptm,
            "ptm_chains": ptm_chains if ptm_chains else np.nan,
            "ptm_chains_auth": ptm_chains_auth,

            "has_any_chain_noncanonical_evidence": has_any_chain_noncanon,
            "noncanonical_chains": noncanon_chains if noncanon_chains else np.nan,
            
            "has_any_chain_invalid_aa": has_any_chain_invalid_aa,
            "invalid_aa_chains": invalid_aa_chains,
            "invalid_aa_chains_auth": invalid_aa_chains_auth,

            "is_binary_interaction": is_binary,
            "binary_n_unique_protein_seqs": nuniq,
            "binary_stoichiometry": stoich,
            "reduced_binary_stoichiometry": reduced_stoich,

            "dna_chains": dna_chains,
            "dna_chains_auth": dna_chains_auth,

            "rna_chains": rna_chains,
            "rna_chains_auth": rna_chains_auth,

            "water_chains": water_chains,
            "water_chains_auth": water_chains_auth,

            "nonpolymer_other_chains": nonpoly_other_chains,
            "nonpolymer_other_chains_auth": nonpoly_other_chains_auth,
        }

        return pd.DataFrame([row])

    # -------------------------
    # Per-PDB binary interaction flags (computed from chain inventory)
    # -------------------------
    @staticmethod
    def _annotate_binary_interactions(inv: pd.DataFrame) -> pd.DataFrame:
        """
        Adds:
          - binary_n_unique_protein_seqs
          - is_binary_interaction
          - binary_stoichiometry
        """
        if inv is None or inv.empty:
            return inv

        out = inv.copy()
        if "entity_poly_type" not in out.columns:
            out["entity_poly_type"] = ""

        def _is_protein_poly_type(x: Any) -> bool:
            if x is None:
                return False
            return str(x).strip().lower().startswith("polypeptide")

        out["pdb_id"] = out["pdb_id"].astype(str).str.strip().str.lower()
        out["seq_can"] = out["seq_can"].astype(str)

        prot_mask = out["entity_poly_type"].apply(_is_protein_poly_type)
        prot_mask &= out["seq_can"].notna() & out["seq_can"].ne("") & out["seq_can"].str.lower().ne("nan")

        prot = out.loc[prot_mask, ["pdb_id", "seq_can"]].copy()
        if prot.empty:
            out["binary_n_unique_protein_seqs"] = 0
            out["is_binary_interaction"] = False
            out["binary_stoichiometry"] = np.nan
            return out

        counts = (
            prot.groupby(["pdb_id", "seq_can"], as_index=False)
                .size()
                .rename(columns={"size": "n_chains"})
        )

        n_unique = (
            counts.groupby("pdb_id", as_index=False)["seq_can"]
                  .nunique()
                  .rename(columns={"seq_can": "binary_n_unique_protein_seqs"})
        )

        top2 = (
            counts.sort_values(["pdb_id", "n_chains"], ascending=[True, False])
                  .groupby("pdb_id", as_index=False)
                  .head(2)
        )

        stoich_rows = []
        for pdb_id, g in top2.groupby("pdb_id"):
            vals = g["n_chains"].tolist()
            if len(vals) == 2:
                vals = sorted(vals, reverse=True)
                stoich = f"{vals[0]}:{vals[1]}"
            else:
                stoich = np.nan
            stoich_rows.append({"pdb_id": pdb_id, "binary_stoichiometry": stoich})

        stoich = pd.DataFrame(stoich_rows)

        out = out.merge(n_unique, on="pdb_id", how="left")
        out = out.merge(stoich[["pdb_id", "binary_stoichiometry"]], on="pdb_id", how="left")

        out["binary_n_unique_protein_seqs"] = out["binary_n_unique_protein_seqs"].fillna(0).astype(int)
        out["is_binary_interaction"] = out["binary_n_unique_protein_seqs"].eq(2)
        out.loc[~out["is_binary_interaction"], "binary_stoichiometry"] = np.nan
        return out

    # -------------------------
    # REST enrichment (kept inside class; JSON strings)
    # -------------------------
    def _build_rest_chain_table(self) -> pd.DataFrame:
        """
        Per-chain REST annotations:
          - UniProt IDs + alignments
          - Pfam / disorder / disorder_binding
          - mutation positions
        """
        pdb = self.pdb_id.upper()

        entry = self._get_json(f"https://data.rcsb.org/rest/v1/core/entry/{pdb}")
        poly_entity_ids = entry.get("rcsb_entry_container_identifiers", {}).get("polymer_entity_ids", [])
        if not poly_entity_ids:
            return pd.DataFrame(columns=["pdb_id", "chain_id", "entity_id"])

        rows: List[Dict[str, Any]] = []

        for polymer_entity in poly_entity_ids:
            pol = self._get_json(f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb}/{polymer_entity}")
            
            entity_poly = pol.get("entity_poly", {}) or {}
            conf = {
                "rcsb_conflict_count": entity_poly.get("rcsb_conflict_count", np.nan),
                "rcsb_deletion_count": entity_poly.get("rcsb_deletion_count", np.nan),
                "rcsb_insertion_count": entity_poly.get("rcsb_insertion_count", np.nan),
                "rcsb_mutation_count": entity_poly.get("rcsb_mutation_count", np.nan),
                "rcsb_non_std_monomer_count": entity_poly.get("rcsb_non_std_monomer_count", np.nan),
                "rcsb_sample_sequence_length": entity_poly.get("rcsb_sample_sequence_length", np.nan),
            }
            conf_json = self._json_or_nan(conf)

            ids = pol.get("rcsb_polymer_entity_container_identifiers", {}) or {}
            chains = ids.get("asym_ids", []) or []
            auths  = ids.get("auth_asym_ids", []) or []
            entity_id = clean_str(ids.get("entity_id", ""))

            # UniProt alignments
            alignment_info: Dict[str, Any] = {}
            for al in pol.get("rcsb_polymer_entity_align", []) or []:
                if al.get("reference_database_name") != "UniProt":
                    continue
                uniprot_id = al.get("reference_database_accession")
                if not uniprot_id:
                    continue
                alignment_info.setdefault(uniprot_id, [])
                for reg in al.get("aligned_regions", []) or []:
                    alignment_info[uniprot_id].append({
                        "entity_beg_seq_id": reg.get("entity_beg_seq_id"),
                        "length": reg.get("length"),
                        "ref_beg_seq_id": reg.get("ref_beg_seq_id"),
                    })

            uniprot_ids = sorted(list(alignment_info.keys()))
            uniprot_ids_str = ",".join(uniprot_ids) if uniprot_ids else np.nan
            uniprot_align_json = self._json_or_nan(alignment_info if alignment_info else None)

            # Features
            pfam = []
            disorder = []
            disorder_binding = []
            mutation_pos: Set[int] = set()

            for feat in pol.get("rcsb_polymer_entity_feature", []) or []:
                ftype = feat.get("type", "")
                if ftype == "Pfam":
                    pf = {"feature_id": feat.get("feature_id"), "name": feat.get("name"), "feature_positions": []}
                    for fp in feat.get("feature_positions", []) or []:
                        pf["feature_positions"].append({
                            "beg_seq_id": fp.get("beg_seq_id"),
                            "end_seq_id": fp.get("end_seq_id"),
                        })
                    pfam.append(pf)

                elif ftype == "disorder":
                    d = {"source": feat.get("provenance_source"), "name": feat.get("name"), "feature_positions": []}
                    for fp in feat.get("feature_positions", []) or []:
                        d["feature_positions"].append({
                            "beg_seq_id": fp.get("beg_seq_id"),
                            "values": fp.get("values"),
                        })
                    disorder.append(d)

                elif ftype == "disorder_binding":
                    d = {"source": feat.get("provenance_source"), "name": feat.get("name"), "feature_positions": []}
                    for fp in feat.get("feature_positions", []) or []:
                        d["feature_positions"].append({
                            "beg_seq_id": fp.get("beg_seq_id"),
                            "values": fp.get("values"),
                        })
                    disorder_binding.append(d)

                elif ftype == "mutation":
                    for fp in feat.get("feature_positions", []) or []:
                        beg = fp.get("beg_seq_id")
                        end = fp.get("end_seq_id", beg)
                        if isinstance(beg, int) and isinstance(end, int):
                            for k in range(beg, end + 1):
                                mutation_pos.add(k)

            # chain/auth pairing (best-effort)
            if len(auths) == len(chains) and len(chains) > 0:
                chain_auth_pairs = list(zip(chains, auths))
            else:
                chain_auth_pairs = [(c, auths[0] if auths else c) for c in chains]

            for chain_id, auth_id in chain_auth_pairs:
                rows.append({
                    "pdb_id": self.pdb_id,
                    "chain_id": str(chain_id).strip(),
                    "entity_id": entity_id,

                    "chain_id_auth_rest": str(auth_id).strip(),

                    "uniprot_ids": uniprot_ids_str,
                    "uniprot_alignments_json": uniprot_align_json,

                    "pfam_json": self._json_or_nan(pfam if pfam else None),
                    "disordered_regions_json": self._json_or_nan(disorder if disorder else None),
                    "disordered_binding_sites_json": self._json_or_nan(disorder_binding if disorder_binding else None),
                    "mutation_pos_can": _positions_to_csv(mutation_pos),
                    
                    # Sequence conflicts
                    "rcsb_conflicts_json": conf_json,
                    "rcsb_conflict_count": conf["rcsb_conflict_count"],
                    "rcsb_deletion_count": conf["rcsb_deletion_count"],
                    "rcsb_insertion_count": conf["rcsb_insertion_count"],
                    "rcsb_mutation_count": conf["rcsb_mutation_count"],
                    "rcsb_non_std_monomer_count": conf["rcsb_non_std_monomer_count"],
                    "rcsb_sample_sequence_length": conf["rcsb_sample_sequence_length"],
                })

        return pd.DataFrame(rows)
    
    def build_entity_inventory(self, *, include_rcsb_rest: bool = True) -> pd.DataFrame:
        """
        Build and store ENTITY-LEVEL inventory (one row per entity_id).
        """
        # ensure chain inventory exists (and has REST if requested)
        if self.inventory_df is None:
            self.build_inventory(include_rcsb_rest=include_rcsb_rest)
        else:
            if include_rcsb_rest and ("uniprot_ids" not in self.inventory_df.columns):
                self.build_inventory(include_rcsb_rest=True)

        self.entity_inventory_df = self._build_entity_inventory_from_chain_inventory(self.inventory_df)
        return self.entity_inventory_df

    def get_entity_inventory_df(self) -> pd.DataFrame:
        if self.entity_inventory_df is None:
            raise RuntimeError("entity_inventory_df is not built yet. Call build_entity_inventory() first.")
        return self.entity_inventory_df

    def generate_report(
        self,
        *,
        include_rcsb_rest: bool = True,
        output: Literal["string", "file"] = "string",
        out_path: Optional[str | Path] = None,
    ) -> str | Path:
        """
        Generate a concise bullet-style report for this PDB.

        Includes:
          - ENTRY-LEVEL structure flags (e.g. entry_has_protein_modification)
          - chain count, protein chain count
          - binary interaction summary
          - PER-CHAIN summaries as before
        """
        # Ensure chain inventory exists with the requested REST setting
        if self.inventory_df is None:
            self.build_inventory(include_rcsb_rest=include_rcsb_rest)
        else:
            if include_rcsb_rest and ("uniprot_ids" not in self.inventory_df.columns):
                self.build_inventory(include_rcsb_rest=True)

        # Ensure entry inventory exists (separate table)
        if self.entry_inventory_df is None:
            self.build_entry_inventory()
            
        # Ensure entity inventory exists
        if self.entity_inventory_df is None:
            self.build_entity_inventory(include_rcsb_rest=include_rcsb_rest)

        inv = self.inventory_df
        entry = self.entry_inventory_df
        ent = self.entity_inventory_df

        if inv is None or inv.empty:
            text = f"{self.pdb_id.upper()} — no inventory rows.\n"
            return self._report_output(text, output=output, out_path=out_path)

        inv = inv.copy()
        inv["pdb_id"] = inv["pdb_id"].astype(str).str.lower()

        lines: list[str] = []
        pdb = self.pdb_id.upper()

        # ---- Entry-level header
        lines.append(f"{pdb}")

        if entry is not None and not entry.empty:
            e0 = entry.iloc[0]
            entry_ptm = e0.get("entry_has_protein_modification", np.nan)
            entry_ptm_str = "NA" if pd.isna(entry_ptm) else str(bool(entry_ptm))
            lines.append(f"- Entry flag has_protein_modification: {entry_ptm_str}")

        # ---- Basic counts + binary summary
        is_binary = bool(entry.iloc[0].get("is_binary_interaction", False)) if entry is not None and len(entry) else False
        stoich = entry.iloc[0].get("binary_stoichiometry", np.nan) if entry is not None and len(entry) else np.nan
        nuniq = int(entry.iloc[0].get("binary_n_unique_protein_seqs", 0) or 0) if entry is not None and len(entry) else 0

        if "entity_poly_type" in inv.columns:
            n_prot = int((inv["entity_poly_type"].astype(str).str.lower().str.startswith("polypeptide")).sum())
        else:
            n_prot = 0

        lines.append(f"- Chains: {inv['chain_id'].nunique()} (protein chains: {n_prot})")
        lines.append(f"- Binary interaction (protein seqs): {is_binary} (n_unique={nuniq}, stoich={stoich if pd.notna(stoich) else 'NA'})")
        
        # ---- Entity summary
        if ent is not None and not ent.empty:
            lines.append("- Entities:")

            ent2 = ent.copy()
            if "entity_id" in ent2.columns:
                ent2 = ent2.sort_values(["entity_id"])

            for _, erow in ent2.iterrows():
                eid = str(erow.get("entity_id", "")).strip()
                etype = erow.get("entity_type", np.nan)
                edesc = erow.get("entity_description", np.nan)

                etype_str = str(etype).strip() if pd.notna(etype) and str(etype).strip() else ""
                edesc_str = str(edesc).strip() if pd.notna(edesc) and str(edesc).strip() else ""

                chains_disp = erow.get("chain_display", np.nan)
                if pd.isna(chains_disp) or not str(chains_disp).strip():
                    # fallback if you removed chain_display from entity inventory
                    cids = str(erow.get("chain_ids", "")).strip()
                    auths = str(erow.get("chain_id_auths", "")).strip()
                    if cids and auths:
                        c_list = cids.split("|")
                        a_list = auths.split("|")
                        # best-effort pairing: if sets, just print both; if equal len, zip
                        if len(c_list) == len(a_list):
                            chains_disp = ",".join([self._fmt_chain(c, a) for c, a in zip(c_list, a_list)])
                        else:
                            chains_disp = f"chains={cids}; auths={auths}"
                    else:
                        chains_disp = "NA"
                else:
                    chains_disp = str(chains_disp).strip()

                varying = erow.get("entity_varying_columns", np.nan)
                varying_str = str(varying).strip() if pd.notna(varying) and str(varying).strip() else "none"

                header_bits = [f"entity {eid}"]
                if etype_str:
                    header_bits.append(etype_str)
                if edesc_str:
                    header_bits.append(edesc_str)

                lines.append(f"  - " + " | ".join(header_bits))
                lines.append(f"    - Chains: {chains_disp}")
                lines.append(f"    - Varying columns: {varying_str}")

        # ---- Per-chain bullets (unchanged logic except no entry-flag mention)
        inv = inv.sort_values(["chain_id"])
        lines.append("- Chains:")
        for _, row in inv.iterrows():
            chain = str(row.get("chain_id", "")).strip()
            auth = str(row.get("chain_id_auth", "")).strip()
            eid = str(row.get("entity_id", "")).strip()
            ptype = str(row.get("entity_poly_type", "")).strip()
            ptype = "non-polymer" if ptype=="nan" else ptype

            seq_len = int(row.get("seq_len", 0) or 0)
            n_res = int(row.get("n_resolved", 0) or 0)
            frac_res = row.get("frac_resolved", np.nan)
            frac_str = f"{float(frac_res):.3f}" if pd.notna(frac_res) else "NA"

            uniprot = row.get("uniprot_ids", np.nan)
            uniprot_str = str(uniprot) if pd.notna(uniprot) and str(uniprot).strip() else "NA"

            is_valid = row.get("is_valid_aa", True)
            invalid_letters = row.get("invalid_letters", "")

            has_ptm = bool(row.get("has_ptm", False))
            ptm_pos = row.get("ptm_pos_can", np.nan)
            ptm_comp = row.get("ptm_comp_ids", np.nan)
            ptm_types = row.get("ptm_types", np.nan)
            ptm_cats = row.get("ptm_categories", np.nan)

            has_noncanon = bool(row.get("has_noncanonical_resname", False))
            noncanon_pos = row.get("noncanon_pos_can", np.nan)
            noncanon_resnames = row.get("noncanonical_resnames", np.nan)

            chain_disp = self._fmt_chain(chain, auth)
            lines.append(
                f"  - Chain {chain_disp}"
                f" | entity {eid}"
                f"{f' | {ptype}' if ptype else ''}"
            )

            # we only care about invalid if it's a polymer
            if not(ptype=="non-polymer"):
                lines.append(
                    f"    - Lengths: canonical={seq_len}, resolved={n_res} (frac={frac_str})"
                    f" | UniProt={uniprot_str}"
                )
                
                if "polypeptide" in ptype:
                    if not bool(is_valid):
                        inv_str = str(invalid_letters) if pd.notna(invalid_letters) else "unknown"
                        lines.append(f"    - Invalids: invalid={inv_str}")
                    else:
                        lines.append(f"    - Invalids: none")
                        
                    if has_noncanon:
                        lines.append(
                            "      - Noncanonicals: "
                            f"pos={noncanon_pos if pd.notna(noncanon_pos) and str(noncanon_pos).strip() else 'NA'}"
                            f"; resnames={noncanon_resnames if pd.notna(noncanon_resnames) else 'NA'}"
                        )
                    else:
                        lines.append("    - Noncanonicals: none detected")

                # PTMs are valid for other molecule types too, but I don't think we care if ligand or water
                if has_ptm:
                    lines.append(
                        "    - PTMs: "
                        f"pos={ptm_pos if pd.notna(ptm_pos) and str(ptm_pos).strip() else 'NA'}"
                        f"; comp_ids={ptm_comp if pd.notna(ptm_comp) else 'NA'}"
                        f"; types={ptm_types if pd.notna(ptm_types) else 'NA'}"
                        f"; categories={ptm_cats if pd.notna(ptm_cats) else 'NA'}"
                    )
                else:
                    lines.append("    - PTMs: none detected")

        text = "\n".join(lines).rstrip() + "\n"
        return self._report_output(text, output=output, out_path=out_path)

    @classmethod
    def _parse_fasta_header(cls, header: str) -> dict[str, Any]:
        """
        header without leading '>'
        """
        parts = header.split("|")
        rec_id = parts[0].strip() if len(parts) > 0 else ""
        chain_field = parts[1].strip() if len(parts) > 1 else ""
        pdb_name = parts[2].strip() if len(parts) > 2 else ""
        org_field = parts[3].strip() if len(parts) > 3 else ""

        pdb_chains, auth_chains = cls._parse_fasta_chain_field(chain_field)

        species = np.nan
        taxid = np.nan
        if org_field:
            m = cls._FASTA_TAXID_RE.match(org_field.strip())
            if m:
                sp = (m.group("species") or "").strip()
                tx = (m.group("taxid") or "").strip()
                species = sp if sp else np.nan
                taxid = int(tx) if tx.isdigit() else np.nan

        return {
            "record_id": rec_id,
            "chains_pdb": ",".join(pdb_chains) if pdb_chains else np.nan,
            "chains_auth": ",".join(auth_chains) if auth_chains else np.nan,
            "pdb_name": pdb_name if pdb_name else np.nan,
            "species": species,
            "species_id": taxid,
        }

    @staticmethod
    def _parse_fasta_chain_field(chain_field: str) -> tuple[list[str], list[str]]:
        s = (chain_field or "").strip()
        s = RCSBStructure._CHAIN_PREFIX_RE.sub("", s).strip()
        if not s:
            return [], []

        parts = [p.strip() for p in s.split(",") if p.strip()]
        pdb_chains: list[str] = []
        auth_chains: list[str] = []

        for p in parts:
            m = RCSBStructure._CHAIN_TOKEN_RE.match(p)
            if not m:
                # last-resort cleanup for weird spacing inside brackets
                p2 = re.sub(r"\[\s+", "[", p)
                p2 = re.sub(r"\s+\]", "]", p2)
                p2 = re.sub(r"\s+", " ", p2).strip()
                m = RCSBStructure._CHAIN_TOKEN_RE.match(p2)

            if m:
                pdb = m.group("pdb")
                auth = m.group("auth") or pdb
                pdb_chains.append(pdb)
                auth_chains.append(auth)
            else:
                # conservative fallback
                pdb = re.split(r"\s|\[", p, maxsplit=1)[0].strip()
                pdb_chains.append(pdb if pdb else p)
                auth_chains.append(pdb if pdb else p)

        return pdb_chains, auth_chains
    
    @staticmethod
    def _fmt_chain(chain_id: Any, auth_id: Any) -> str:
        c = "" if chain_id is None else str(chain_id).strip()
        a = "" if auth_id is None else str(auth_id).strip()

        # normalize common "missing" sentinels
        if c.lower() in ("", "nan", "none", "?", "."):
            c = ""
        if a.lower() in ("", "nan", "none", "?", "."):
            a = ""

        if not c:
            return "NA"
        if a and a != c:
            return f"{c}[auth: {a}]"
        return c

    @staticmethod
    def _build_entity_inventory_from_chain_inventory(chain_inv: pd.DataFrame) -> pd.DataFrame:
        if chain_inv is None or chain_inv.empty:
            return pd.DataFrame()

        df = chain_inv.copy()

        # cols we ALWAYS aggregate as pipe-joined
        always_join = {"chain_id", "chain_id_auth"}

        # ensure stable types
        for c in ["chain_id", "chain_id_auth", "entity_id", "pdb_id"]:
            if c in df.columns:
                df[c] = df[c].astype(str)

        def _norm(x: Any) -> str:
            if x is None:
                return ""
            if isinstance(x, float) and np.isnan(x):
                return ""
            s = str(x).strip()
            if s in ("", "nan", "None", "?", "."):
                return ""
            return s

        def _pipe_unique(series: pd.Series) -> Any:
            vals = [_norm(v) for v in series.tolist()]
            vals = [v for v in vals if v]  # drop empties
            if not vals:
                return np.nan
            uniq = sorted(set(vals))
            return "|".join(uniq)

        def _collapse_group(g: pd.DataFrame) -> dict[str, Any]:
            out: dict[str, Any] = {
                "pdb_id": g["pdb_id"].iloc[0] if "pdb_id" in g.columns else np.nan,
                "entity_id": g["entity_id"].iloc[0],
            }

            varying_cols: set[str] = set()

            # Convenience: precompute chain display strings for summary/reporting
            # chain_ids and chain_id_auths are aligned row-wise, not set-wise
            if "chain_id" in g.columns and "chain_id_auth" in g.columns:
                chain_disp = []
                for cid, auth in zip(g["chain_id"].tolist(), g["chain_id_auth"].tolist()):
                    s = RCSBStructure._fmt_chain(cid, auth)
                    if s != "NA":
                        chain_disp.append(s)
                out["chain_display"] = ",".join(chain_disp) if chain_disp else np.nan

            cols = [c for c in g.columns if c not in ("entity_id",)]

            for c in cols:
                if c in always_join:
                    # keep these as pipe-joined lists, but do NOT count them as "varying"
                    out[c + "s"] = _pipe_unique(g[c])
                    continue

                nn = g[c].dropna()
                if nn.empty:
                    out[c] = np.nan
                    continue

                normed = [_norm(v) for v in nn.tolist()]
                normed = [v for v in normed if v]
                if not normed:
                    out[c] = np.nan
                    continue

                uniq = sorted(set(normed))
                if len(uniq) == 1:
                    out[c] = uniq[0]
                else:
                    out[c] = "|".join(uniq)
                    varying_cols.add(c)

            # NEW: record which columns varied (pipe-joined because >1 unique value)
            out["entity_varying_columns"] = ",".join(sorted(varying_cols)) if varying_cols else np.nan

            return out

        grouped = df.groupby(["pdb_id", "entity_id"], as_index=False)
        rows = [_collapse_group(g) for _, g in grouped]
        return pd.DataFrame(rows)

    def build_fasta_inventory(self, *, force_download: bool = False) -> pd.DataFrame:
        """
        One row per FASTA record (i.e., per polymer entity sequence in the RCSB entry FASTA).
        """
        fasta = self.load_fasta(force=force_download)
        rows: list[dict[str, Any]] = []

        header = None
        seq_lines: list[str] = []

        def _flush():
            nonlocal header, seq_lines
            if header is None:
                return
            meta = self._parse_fasta_header(header)
            seq = "".join(seq_lines).replace(" ", "").replace("\t", "").strip()
            rows.append({
                "pdb_id": self.pdb_id,
                **meta,
                "sequence": seq,
                "seq_len": len(seq),
            })
            header = None
            seq_lines = []

        for line in fasta.splitlines():
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                _flush()
                header = line[1:].strip()
            else:
                seq_lines.append(line)

        _flush()

        df = pd.DataFrame(rows)
        df = harmonize_nulls_to_nan(df)
        self.fasta_records_df = df
        return df

    @staticmethod
    def _report_output(text: str, *, output: str, out_path: Optional[str | Path]) -> str | Path:
        if output == "string":
            return text
        if out_path is None:
            raise ValueError("out_path is required when output='file'")
        out_path = Path(out_path)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(text, encoding="utf-8")
        return out_path

    # -------------------------
    # Tiny class-local JSON/REST utilities
    # -------------------------
    @staticmethod
    def _get_json(url: str, *, timeout: int = 30) -> Dict[str, Any]:
        r = requests.get(url, timeout=timeout)
        if r.status_code != 200:
            raise RuntimeError(f"GET {url} failed: {r.status_code}")
        return r.json()

    @staticmethod
    def _json_or_nan(x: Any) -> Any:
        if x is None:
            return np.nan
        if isinstance(x, float) and np.isnan(x):
            return np.nan
        try:
            return json.dumps(x, ensure_ascii=False)
        except Exception:
            return np.nan


def download_rcsb_apicall(pdb_id, struct_format='cif', output_dir=None):
    full_file_name = f"{pdb_id}.{struct_format}"
    full_file_name_uppercase = f"{pdb_id.upper()}.{struct_format}"
    output_path = full_file_name_uppercase
    if output_dir is not None:
        os.makedirs(output_dir, exist_ok=True)
        output_path = f"{output_dir}/{full_file_name_uppercase}"

    url = f"https://files.rcsb.org/download/{full_file_name}"
    response = requests.get(url)
    return response, output_path


def download_rcsb(pdb_id, struct_format='cif', convert_if_fail=False, output_dir=None):
    """
    Download mmCIF file with provided pdb_id and optional output_dir.
    Return: path to downloaded file if successful.
    """
    uppercase_id = pdb_id.upper()
    lowercase_id = pdb_id.lower()

    for formatted_pdb_id in [uppercase_id, lowercase_id]:
        response, output_path = download_rcsb_apicall(
            formatted_pdb_id, struct_format=struct_format, output_dir=output_dir
        )

        if response.status_code == 200:
            with open(output_path, 'wb') as file:
                file.write(response.content)
            return output_path
        else:
            if struct_format == 'pdb' and convert_if_fail:
                response, output_path = download_rcsb_apicall(
                    formatted_pdb_id, struct_format='cif', output_dir=output_dir
                )
                if response.status_code == 200:
                    with open(output_path, 'wb') as file:
                        file.write(response.content)
                    convert_cif_to_pdb(output_path, output_path.replace('.cif', '.pdb'))
                    os.remove(output_path)
                    raise Exception(
                        f"Deleted original .cif file download. See {output_path.replace('.cif','.pdb')} for the PDB"
                    )
                else:
                    raise Exception(
                        f"Failed to download {pdb_id} file (searched as: {formatted_pdb_id}). Status code: {response.status_code}"
                    )
            else:
                raise Exception(
                    f"Failed to download {pdb_id} file (searched as: {formatted_pdb_id}). Status code: {response.status_code}"
                )

    return None