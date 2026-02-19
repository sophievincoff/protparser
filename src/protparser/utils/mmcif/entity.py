from __future__ import annotations

from typing import Any, Dict

import numpy as np
import gemmi

from protparser.utils.mmcif.cif_columns import (
    safe_find_column,
    col_required,
    cif_get_row_count,
    safe_int,
    clean_str,
)


def build_entity_lookup(block: gemmi.cif.Block) -> dict[str, dict[str, Any]]:
    """
    Build an entity lookup from an mmCIF block.

    Returns:
        entity_id -> dict with keys like:
          - entity_type
          - entity_description
          - entity_src_method
          - entity_formula_weight
          - entity_number_of_molecules
          - entity_mutation
          - entity_ec
          - entity_fragment
          - entity_details
          - comp_id / nonpoly_name (from _pdbx_entity_nonpoly)
          - chem_comp_name / chem_comp_formula / chem_comp_formula_weight (from _chem_comp)
    """
    out: dict[str, dict[str, Any]] = {}

    # -------------------------
    # _entity: id/type/description + extra fields you want
    # -------------------------
    ent = block.find_mmcif_category("_entity")
    if ent is not None:
        ent_id = col_required(ent, "id")

        ent_type = safe_find_column(ent, "type")
        ent_desc = safe_find_column(ent, "pdbx_description")
        ent_src  = safe_find_column(ent, "src_method")
        ent_fw   = safe_find_column(ent, "formula_weight")

        ent_num_mol = safe_find_column(ent, "pdbx_number_of_molecules")
        ent_ec      = safe_find_column(ent, "pdbx_ec")
        ent_mut     = safe_find_column(ent, "pdbx_mutation")
        ent_frag    = safe_find_column(ent, "pdbx_fragment")
        ent_details = safe_find_column(ent, "details")

        src_method_dict = {
            "man": "man (entity isolated from a genetically manipulated source)",
            "nat": "nat (entity isolated from a natural source)",
            "syn": "syn (entity obtained synthetically)"                           
        }
        for i in range(cif_get_row_count(ent_id)):
            eid = clean_str(ent_id[i])
            d = out.setdefault(eid, {})

            d["entity_type"] = clean_str(ent_type[i]) if ent_type is not None else ""
            d["entity_description"] = clean_str(ent_desc[i]) if ent_desc is not None else ""
            d["entity_src_method"] = src_method_dict.get(clean_str(ent_src[i])) if ent_src is not None else ""

            # formula_weight is sometimes numeric-ish but often a string in CIF; keep string + numeric
            fw_raw = clean_str(ent_fw[i]) if ent_fw is not None else ""
            d["entity_formula_weight_raw"] = fw_raw if fw_raw else np.nan
            try:
                d["entity_formula_weight"] = float(fw_raw) if fw_raw not in ("", "?", ".") else np.nan
            except Exception:
                d["entity_formula_weight"] = np.nan

            # NEW: number of molecules (stoichiometry at entity level)
            n_raw = clean_str(ent_num_mol[i]) if ent_num_mol is not None else ""
            d["entity_number_of_molecules_raw"] = n_raw if n_raw else np.nan
            n_int = safe_int(n_raw) if n_raw else None
            d["entity_number_of_molecules"] = int(n_int) if isinstance(n_int, int) else np.nan

            # NEW: mutation string (often comma-separated like "C34S, N63Q")
            mut_raw = clean_str(ent_mut[i]) if ent_mut is not None else ""
            d["entity_mutation"] = mut_raw if mut_raw not in ("", "?", ".") else np.nan

            # other optional goodies
            ec_raw = clean_str(ent_ec[i]) if ent_ec is not None else ""
            d["entity_ec"] = ec_raw if ec_raw not in ("", "?", ".") else np.nan

            frag_raw = clean_str(ent_frag[i]) if ent_frag is not None else ""
            d["entity_fragment"] = frag_raw if frag_raw not in ("", "?", ".") else np.nan

            details_raw = clean_str(ent_details[i]) if ent_details is not None else ""
            d["entity_details"] = details_raw if details_raw not in ("", "?", ".") else np.nan

    # -------------------------
    # _pdbx_entity_nonpoly: entity_id -> (name, comp_id) for non-polymers/water
    # -------------------------
    nonpoly = block.find_mmcif_category("_pdbx_entity_nonpoly")
    if nonpoly is not None:
        # Some CIFs have this category but omit entity_id (or are malformed/empty).
        # Be permissive: if we can't find an entity id column, skip gracefully.
        np_eid = safe_find_column(nonpoly, "entity_id")
        if np_eid is None:
            # (optional) try a couple of known oddball alternatives
            np_eid = safe_find_column(nonpoly, "entity") or safe_find_column(nonpoly, "entityid")

        if np_eid is not None:
            np_name = safe_find_column(nonpoly, "name")
            np_comp = safe_find_column(nonpoly, "comp_id")

            for i in range(cif_get_row_count(np_eid)):
                eid = clean_str(np_eid[i])
                if not eid or eid in ("?", "."):
                    continue

                d = out.setdefault(eid, {})

                if np_comp is not None:
                    comp_id = clean_str(np_comp[i])
                    d["comp_id"] = comp_id if comp_id not in ("", "?", ".") else ""

                if np_name is not None:
                    name = clean_str(np_name[i])
                    d["nonpoly_name"] = name if name not in ("", "?", ".") else np.nan
        # else: skip silently (or log)

    # -------------------------
    # _chem_comp: comp_id -> name/formula (only for components present in entry)
    # -------------------------
    chem = block.find_mmcif_category("_chem_comp")
    chem_by_id: dict[str, dict[str, Any]] = {}
    if chem is not None:
        cc_id = col_required(chem, "id")
        cc_name = safe_find_column(chem, "name")
        cc_formula = safe_find_column(chem, "formula")
        cc_fw = safe_find_column(chem, "formula_weight")

        for i in range(cif_get_row_count(cc_id)):
            cid = clean_str(cc_id[i])
            chem_by_id[cid] = {
                "chem_comp_name": clean_str(cc_name[i]) if cc_name is not None else "",
                "chem_comp_formula": clean_str(cc_formula[i]) if cc_formula is not None else "",
                "chem_comp_formula_weight": clean_str(cc_fw[i]) if cc_fw is not None else "",
            }

    # attach chem_comp info to entities that have a comp_id
    for eid, d in out.items():
        comp = d.get("comp_id", "")
        if comp and comp in chem_by_id:
            d.update(chem_by_id[comp])

    return out