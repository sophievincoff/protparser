# Parse mmCIF modification annotations (_pdbx_modification_feature / entry flags) and
# convert them into per-chain, per-residue tag maps for sequence annotation (e.g., disulfide, ACE). 
# 
# entry_has_protein_modification_flag: Reads the entry-level “has_protein_modification” flag (Y/N) if present.
# chain_has_modification_feature: Checks whether a specific chain is referenced in _pdbx_modification_feature.
# extract_mod_features_by_chain: Summarizes modification features per chain (categories/types/residues + examples).
# extract_position_tags_from_modification_feature: Converts modification features into {chain: {pos: {tags}}}.
# mk_mod_label: Normalizes labels used in annotation tokens (e.g., disulfide).
# annotate_seq_with_position_tags: Replaces residues in a sequence with modification tokens at tagged positions.
# annotate_masked_seq_with_position_tags: Same as above, but preserves masked (‘-’) positions and only annotates resolved residues.
# add_pos_tag: Adds a tag to a position map with basic input validation.
# merge_pos_tags: Merges multiple {chain: {pos: tags}} maps into one.

from __future__ import annotations
from typing import Any, Dict, Optional, Set
from collections import defaultdict
import gemmi
import numpy as np

from .cif_columns import safe_find_column, cif_get_row_count, clean_str, safe_int, truthy_cif_flag
from .seq_cleaning import ALLOWED_POLYMER_RESNAMES, COMMON_PTM_RESNAMES, AA3_TO_1, PTM_RESNAME_TO_BASE

def entry_has_protein_modification_flag(block: gemmi.cif.Block) -> Optional[bool]:
    """
    If present, reads _pdbx_entry_details.has_protein_modification (Y/N).
    Returns True/False or None if missing/unparseable.
    """
    det = block.find_mmcif_category("_pdbx_entry_details")
    if det is None:
        return None
    col = safe_find_column(det, "has_protein_modification")
    if col is None:
        return None
    # Usually one row; take first non-missing
    try:
        v = col[0]
    except Exception:
        return None
    return truthy_cif_flag(v)


def chain_has_modification_feature(block: gemmi.cif.Block, chain_id_label: str, chain_id_auth: Optional[str]) -> Optional[bool]:
    """
    If present, uses _pdbx_modification_feature to see if this chain is listed.
    Returns True/False or None if category/columns absent.
    """
    mf = block.find_mmcif_category("_pdbx_modification_feature")
    if mf is None:
        return None

    # These vary; try a few likely column names
    col_label = safe_find_column(mf, "label_asym_id")
    col_auth  = safe_find_column(mf, "auth_asym_id")

    if col_label is None and col_auth is None:
        return None

    n = len(mf) if hasattr(mf, "__len__") else mf.get_row_count()
    chain_id_label = str(chain_id_label).strip()
    chain_id_auth = str(chain_id_auth).strip() if chain_id_auth is not None else None

    for i in range(n):
        if col_label is not None:
            v = str(col_label[i]).strip()
            if v and v not in {".", "?"} and v == chain_id_label:
                return True
        if col_auth is not None and chain_id_auth:
            v = str(col_auth[i]).strip()
            if v and v not in {".", "?"} and v == chain_id_auth:
                return True

    return False

def extract_mod_features_by_chain(block: gemmi.cif.Block) -> Dict[str, Dict[str, Any]]:
    """
    Returns dict keyed by label_asym_id (chain) with:
      categories, types, comp_ids, n, examples(list)

    Robust to missing category / missing columns / gemmi RuntimeError on find_column.
    """
    cat = block.find_mmcif_category("_pdbx_modification_feature")
    if cat is None:
        return {}

    # Use safe wrapper everywhere (gemmi may raise on missing columns)
    col_chain = safe_find_column(cat, "label_asym_id")
    col_comp  = safe_find_column(cat, "label_comp_id")
    col_seq   = safe_find_column(cat, "label_seq_id")
    col_type  = safe_find_column(cat, "type")
    col_cat   = safe_find_column(cat, "category")

    # partner residue fields (for disulfides/crosslinks)
    col_chain2 = safe_find_column(cat, "modified_residue_label_asym_id")
    col_comp2  = safe_find_column(cat, "modified_residue_label_comp_id")
    col_seq2   = safe_find_column(cat, "modified_residue_label_seq_id")

    # If we can't even attribute a feature to a chain, we can't build per-chain features.
    # (Some CIFs may only have auth_* columns, or different schema.)
    if col_chain is None:
        # OPTIONAL fallback: try auth_asym_id if present
        col_chain = safe_find_column(cat, "auth_asym_id")
        col_chain2 = safe_find_column(cat, "modified_residue_auth_asym_id") if col_chain2 is None else col_chain2

    if col_chain is None:
        return {}

    out: Dict[str, Dict[str, Any]] = defaultdict(lambda: {
        "categories": set(),
        "types": set(),
        "comp_ids": set(),
        "n": 0,
        "examples": []
    })

    n = cif_get_row_count(cat)

    def _norm(x: Any) -> str:
        s = str(x).strip() if x is not None else ""
        return "" if (not s or s in {".", "?"}) else s

    def _add(chain: Any, comp: Any, seq: Any, typ: Any, catg: Any,
             chain_other: Any, comp_other: Any, seq_other: Any):
        c = _norm(chain)
        if not c:
            return

        d = out[c]
        d["n"] += 1

        cc = _norm(comp)
        ss = _norm(seq)
        tt = _norm(typ)
        gg = _norm(catg)

        if cc:
            d["comp_ids"].add(cc)
        if tt:
            d["types"].add(tt)
        if gg:
            d["categories"].add(gg)

        if len(d["examples"]) < 5:
            other = ""
            c2 = _norm(chain_other)
            if c2:
                other = f"-> {c2}:{_norm(comp_other)}{_norm(seq_other)}"
            d["examples"].append(f"{c}:{cc}{ss} {gg} {other}".strip())

    for i in range(n):
        _add(
            col_chain[i],
            col_comp[i] if col_comp is not None else None,
            col_seq[i]  if col_seq  is not None else None,
            col_type[i] if col_type is not None else None,
            col_cat[i]  if col_cat  is not None else None,
            col_chain2[i] if col_chain2 is not None else None,
            col_comp2[i]  if col_comp2  is not None else None,
            col_seq2[i]   if col_seq2   is not None else None,
        )
        if col_chain2 is not None:
            _add(
                col_chain2[i],
                col_comp2[i] if col_comp2 is not None else None,
                col_seq2[i]  if col_seq2  is not None else None,
                col_type[i]  if col_type  is not None else None,
                col_cat[i]   if col_cat   is not None else None,
                col_chain[i],
                col_comp[i] if col_comp is not None else None,
                col_seq[i]  if col_seq  is not None else None,
            )

    return out

def extract_position_tags_from_modification_feature(block: gemmi.cif.Block) -> Dict[str, Dict[int, Set[str]]]:
    """
    Parse _pdbx_modification_feature into:
      { chain_label_asym_id : { position(int) : set(tags) } }

    Heuristics:
      # For terminal mods like ACE/NH2: use modified_residue_label_seq_id (polymer residue)
      # For disulfides: tag BOTH label_seq_id and modified_residue_label_seq_id on their chains
      # Otherwise: tag the polymer residue position we can identify (prefer modified_residue_*, else label_*)
    """
    cat = block.find_mmcif_category("_pdbx_modification_feature")
    if cat is None:
        return {}

    # columns we might use (some optional)
    col_chain1 = safe_find_column(cat, "label_asym_id")
    col_seq1   = safe_find_column(cat, "label_seq_id")
    col_comp1  = safe_find_column(cat, "label_comp_id")

    col_chain2 = safe_find_column(cat, "modified_residue_label_asym_id")
    col_seq2   = safe_find_column(cat, "modified_residue_label_seq_id")
    col_comp2  = safe_find_column(cat, "modified_residue_label_comp_id")

    col_cat    = safe_find_column(cat, "category")
    col_type   = safe_find_column(cat, "type")

    if col_chain1 is None:
        return {}

    out: Dict[str, Dict[int, Set[str]]] = defaultdict(lambda: defaultdict(set))

    n = cif_get_row_count(cat)
    for i in range(n):
        chain1 = clean_str(col_chain1[i] if col_chain1 is not None else None)
        chain2 = clean_str(col_chain2[i] if col_chain2 is not None else None)

        seq1 = safe_int(col_seq1[i] if col_seq1 is not None else None)
        seq2 = safe_int(col_seq2[i] if col_seq2 is not None else None)

        comp1 = clean_str(col_comp1[i] if col_comp1 is not None else None)  # e.g. CYS, ACE, NH2
        comp2 = clean_str(col_comp2[i] if col_comp2 is not None else None)  # e.g. CYS, THR, PRO
        catv  = clean_str(col_cat[i]  if col_cat is not None  else None)    # e.g. 'Disulfide bridge'
        typev = clean_str(col_type[i] if col_type is not None else None)

        # decide "what to call the modification" (right side of token)
        # # disulfide is best labeled by category
        # # terminal mods often show comp1=ACE/NH2 and category describes it too; use comp1 ("ACE", "NH2")
        # # otherwise use comp1 if it looks like a special residue, else fallback to category/type
        label = ""
        if catv:
            label = mk_mod_label(catv)  # normalizes disulfide -> "disulfide"
        if label == "disulfide":
            # tag both sides at their polymer positions (seq1/seq2)
            if chain1:
                add_pos_tag(out[chain1], seq1, "disulfide")
            if chain2:
                add_pos_tag(out[chain2], seq2, "disulfide")
            continue

        # terminal mods: comp1 often ACE or NH2
        if comp1 in {"ACE", "NH2"}:
            # tag the polymer residue position (usually seq2), not the "ACE/NH2 pseudo-residue"
            if chain2:
                add_pos_tag(out[chain2], seq2, comp1)
            elif chain1:
                # fallback if modified_residue_* absent
                add_pos_tag(out[chain1], seq1, comp1)
            continue

        # general case: attach something meaningful
        label = mk_mod_label(comp1) or mk_mod_label(catv) or mk_mod_label(typev)
        # prefer modified residue position if present
        if chain2 and seq2 is not None:
            add_pos_tag(out[chain2], seq2, label)
        elif chain1 and seq1 is not None:
            add_pos_tag(out[chain1], seq1, label)

    return out

def mk_mod_label(comp_or_cat: str) -> str:
    """
    Normalize the PTM label that goes on the right of '+'
    """
    s = clean_str(comp_or_cat)
    if not s:
        return ""
    # standardize common ones
    low = s.lower()
    if "disulfide" in low:
        return "disulfide"
    # keep ACE/NH2 as-is, and otherwise keep original category/comp_id
    return s

def add_pos_tag(tags_by_pos: Dict[int, Set[str]], pos: Optional[int], tag: str):
    if pos is None or pos <= 0:
        return
    if not tag:
        return
    tags_by_pos[pos].add(tag)
    
def merge_pos_tags(*maps: Dict[str, Dict[int, Set[str]]]) -> Dict[str, Dict[int, Set[str]]]:
    out: Dict[str, Dict[int, Set[str]]] = defaultdict(lambda: defaultdict(set))
    for m in maps:
        for chain, posmap in (m or {}).items():
            for pos, tags in posmap.items():
                out[chain][pos].update(tags)
    return out

def _compute_chain_flags(block: gemmi.cif.Block, chain_id_label: str, chain_id_auth: Optional[str]) -> Dict[str, Any]:
    """
    Returns a dict of chain-level flags + debug metadata.
    """
    resnames = _scan_chain_resnames_from_atom_site(block, chain_id_label)

    # ---- noncanonical residue name detection ----
    noncanon = sorted({r for r in resnames if r not in ALLOWED_POLYMER_RESNAMES})
    has_noncanonical = len(noncanon) > 0

    # ---- PTM detection components ----
    mf = _chain_has_modification_feature(block, chain_id_label, chain_id_auth)  # True/False/None
    entry_flag = _entry_has_protein_modification_flag(block)                   # True/False/None

    # residue-name heuristic
    ptm_resnames = sorted({r for r in resnames if r in COMMON_PTM_RESNAMES})
    resname_heur = len(ptm_resnames) > 0

    sources = []
    if mf is True:
        sources.append("pdbx_modification_feature")
    if entry_flag is True:
        sources.append("entry_flag")
    if resname_heur:
        sources.append("resname_heuristic")

    # remove entry_flag if it's the ONLY source
    if set(sources)==set(["entry_flag"]):
        sources = []
        
    has_ptm = len(sources) > 0

    # ---- human-readable reason ----
    if not has_ptm:
        reason = ""
    else:
        parts = []
        if mf is True:
            parts.append("chain listed in _pdbx_modification_feature")
        # Commented out below because we do not want entry_flag contributing to chain features
        #if entry_flag is True:
            #parts.append("_pdbx_entry_details.has_protein_modification=Y")
        if resname_heur:
            parts.append(f"atom_site residue names contain {','.join(ptm_resnames)}")
        reason = "; ".join(parts)
        
    return {
        "has_ptm": bool(has_ptm),
        "ptm_detect_sources": ",".join(sources),
        "ptm_resnames": ",".join(ptm_resnames),
        "ptm_entry_flag": entry_flag,     # may be None
        "ptm_mf_chain_hit": mf,           # may be None

        "has_noncanonical_resname": bool(has_noncanonical),
        "noncanonical_resnames": ",".join(noncanon),
        "has_ptm_reason": reason,
    }
    
def _chain_has_modification_feature(block: gemmi.cif.Block, chain_id_label: str, chain_id_auth: Optional[str]) -> Optional[bool]:
    """
    If present, uses _pdbx_modification_feature to see if this chain is listed.
    Returns True/False or None if category/columns absent.
    """
    mf = block.find_mmcif_category("_pdbx_modification_feature")
    if mf is None:
        return None

    # These vary; try a few likely column names
    col_label = safe_find_column(mf, "label_asym_id")
    col_auth  = safe_find_column(mf, "auth_asym_id")

    if col_label is None and col_auth is None:
        return None

    n = len(mf) if hasattr(mf, "__len__") else mf.get_row_count()
    chain_id_label = str(chain_id_label).strip()
    chain_id_auth = str(chain_id_auth).strip() if chain_id_auth is not None else None

    for i in range(n):
        if col_label is not None:
            v = str(col_label[i]).strip()
            if v and v not in {".", "?"} and v == chain_id_label:
                return True
        if col_auth is not None and chain_id_auth:
            v = str(col_auth[i]).strip()
            if v and v not in {".", "?"} and v == chain_id_auth:
                return True

    return

def _entry_has_protein_modification_flag(block: gemmi.cif.Block) -> Optional[bool]:
    """
    If present, reads _pdbx_entry_details.has_protein_modification (Y/N).
    Returns True/False or None if missing/unparseable.
    """
    det = block.find_mmcif_category("_pdbx_entry_details")
    if det is None:
        return None
    col = safe_find_column(det, "has_protein_modification")
    if col is None:
        return None
    # Usually one row; take first non-missing
    try:
        v = col[0]
    except Exception:
        return None
    return truthy_cif_flag(v)

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

def _parse_yes_no_flag(x: Any) -> Any:
    """
    Parse mmCIF-ish yes/no flags robustly.
    Returns True/False, or np.nan if unknown/missing.
    Accepts: Y/N, YES/NO, TRUE/FALSE, 1/0, with or without quotes/prefixes like '=Y'.
    """
    if x is None:
        return np.nan
    if isinstance(x, float) and np.isnan(x):
        return np.nan
    if isinstance(x, (bool, np.bool_)):
        return bool(x)
    if isinstance(x, (int, np.integer)):
        if x == 1:
            return True
        if x == 0:
            return False
        return np.nan

    s = str(x).strip().upper()

    # strip common quoting
    if len(s) >= 2 and ((s[0] == s[-1] == "'") or (s[0] == s[-1] == '"')):
        s = s[1:-1].strip().upper()

    # handle patterns like "=Y", "HAS_PROTEIN_MODIFICATION=Y"
    if "=" in s:
        s = s.split("=")[-1].strip().upper()

    # final cleanup
    s = s.strip(" ;,")

    if s in {"Y", "YES", "TRUE", "T", "1"}:
        return True
    if s in {"N", "NO", "FALSE", "F", "0"}:
        return False

    return np.nan

def annotate_seq_with_position_tags(seq: str, pos_tags: Dict[int, Set[str]]) -> str:
    """
    Replace positions in a 1-letter sequence with tokens like:
      <Y+PTR> or <C+disulfide> or <T+ACE>
    If multiple tags at same position, uses <A+TAG1|TAG2>.

    NOTE: Output string length will differ from original.
    """
    if not seq:
        return seq

    # Fast path: no tags
    if not pos_tags:
        return seq

    out_parts: List[str] = []
    L = len(seq)

    for pos1 in range(1, L + 1):
        aa = seq[pos1 - 1]
        tags = pos_tags.get(pos1, None)
        if not tags:
            out_parts.append(aa)
            continue

        # normalize + dedup tags
        norm = []
        for t in sorted(tags):
            tt = mk_mod_label(t)
            if tt:
                norm.append(tt)

        if not norm:
            out_parts.append(aa)
            continue

        # if the tag itself is a PTM residue name like PTR/SEP/TPO/MSE, prefer base letter mapping
        # (but still keep the observed letter as "original" if you want — here we override only if aa looks wrong)
        orig = aa
        # If aa is '-' (masked), preserve '-' as orig in token? Usually you want the real base. We'll infer base from tag.
        if orig == "-":
            # choose a base from the first tag if possible
            orig = base_letter_for_position("", 0, fallback_comp=norm[0])

        tag_str = "|".join(norm)
        out_parts.append(f"<{orig}+{tag_str}>")

    return "".join(out_parts)

def annotate_masked_seq_with_position_tags(masked_seq: str, full_seq: str, pos_tags: Dict[int, Set[str]]) -> str:
    """
    Annotate ONLY resolved positions (non '-' in masked_seq).
    If masked_seq[pos] == '-', keep '-' even if there is a tag at that canonical position.
    Tokens use base letter from full_seq.
    """
    if not masked_seq:
        return masked_seq
    if not pos_tags:
        return masked_seq

    out_parts: List[str] = []
    L = len(masked_seq)

    for pos1 in range(1, L + 1):
        ch = masked_seq[pos1 - 1]

        # Keep unresolved as dash, always
        if ch == "-":
            out_parts.append("-")
            continue

        tags = pos_tags.get(pos1, None)
        if not tags:
            out_parts.append(ch)
            continue

        base = base_letter_for_position(full_seq, pos1, fallback_comp=None)

        norm = []
        for t in sorted(tags):
            tt = mk_mod_label(t)
            if tt:
                norm.append(tt)

        if not norm:
            out_parts.append(ch)
            continue

        out_parts.append(f"<{base}+{'|'.join(norm)}>")
    return "".join(out_parts)

def base_letter_for_position(seq: str, pos1: int, fallback_comp: Optional[str] = None) -> str:
    """
    Given 1-indexed pos, return base letter:
      - from seq if in-range
      - else from PTM_RESNAME_TO_BASE/AA3_TO_1 if fallback_comp provided
      - else 'X'
    """
    if seq and 1 <= pos1 <= len(seq):
        return seq[pos1 - 1]
    if fallback_comp:
        c = fallback_comp.upper()
        if c in PTM_RESNAME_TO_BASE:
            return PTM_RESNAME_TO_BASE[c]
        if c in AA3_TO_1:
            return AA3_TO_1[c]
    return "X"

def build_entity_source_lookup(block: gemmi.cif.Block) -> dict[str, dict[str, Any]]:
    """
    Best-effort entity_id -> {'species_name': str/np.nan, 'taxid': int/np.nan}
    Uses common mmCIF source categories.
    """
    out: dict[str, dict[str, Any]] = {}

    def _set(eid: str, name: Any, taxid: Any):
        eid = clean_str(eid)
        if not eid or eid in ("?", "."):
            return
        d = out.setdefault(eid, {"species_name": np.nan, "taxid": np.nan})

        n = clean_str(name) if name is not None else ""
        if n and n not in ("?", "."):
            d["species_name"] = n

        tx = clean_str(taxid) if taxid is not None else ""
        # taxid sometimes stored as string/int-ish
        if tx and tx not in ("?", "."):
            try:
                d["taxid"] = int(tx)
            except Exception:
                pass

    # ---- Genetically engineered source ----
    cat = block.find_mmcif_category("_entity_src_gen")
    if cat is not None:
        eid = safe_find_column(cat, "entity_id")
        name = safe_find_column(cat, "pdbx_gene_src_scientific_name")
        tax  = safe_find_column(cat, "pdbx_gene_src_ncbi_taxonomy_id")
        if eid is not None:
            for i in range(cif_get_row_count(eid)):
                _set(eid[i], name[i] if name is not None else None, tax[i] if tax is not None else None)

    # ---- Natural source ----
    cat = block.find_mmcif_category("_entity_src_nat")
    if cat is not None:
        eid = safe_find_column(cat, "entity_id")
        name = safe_find_column(cat, "pdbx_organism_scientific")
        tax  = safe_find_column(cat, "pdbx_ncbi_taxonomy_id")
        if eid is not None:
            for i in range(cif_get_row_count(eid)):
                _set(eid[i], name[i] if name is not None else None, tax[i] if tax is not None else None)

    # ---- Synthetic source (schema varies; try common tag names) ----
    cat = block.find_mmcif_category("_pdbx_entity_src_syn")
    if cat is not None:
        eid = safe_find_column(cat, "entity_id")
        name = (
            safe_find_column(cat, "organism_scientific") or
            safe_find_column(cat, "pdbx_organism_scientific") or
            safe_find_column(cat, "organism_scientific_name")
        )
        tax = (
            safe_find_column(cat, "ncbi_taxonomy_id") or
            safe_find_column(cat, "pdbx_ncbi_taxonomy_id")
        )
        if eid is not None:
            for i in range(cif_get_row_count(eid)):
                _set(eid[i], name[i] if name is not None else None, tax[i] if tax is not None else None)

    return out
