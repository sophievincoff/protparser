# Sequence-cleaning utilities to normalize mmCIF polymer sequences and
# validate that they contain only allowed amino-acid characters.

# _seq_is_clean: Checks whether a sequence contains only allowed amino-acid letters and returns a failure reason if not.
# clean_mmcif_seq_text: Cleans raw mmCIF sequence text, removing semicolon blocks, whitespace, and newlines.

from __future__ import annotations
import re
from typing import Optional, Tuple, Set

# Natural 20 + selenocysteine U
VALID_AAS = set("ACDEFGHIKLMNPQRSTVWYU")
DISALLOWED_AAS_DEFAULT = set(["X", "O", "U","B", "Z", "J"])  # X unknown, O pyrrolysine, other characters I tend not to like

# -----------------------------
# Residue-name / PTM detection
# -----------------------------

STANDARD_RESNAMES_20 = {
    "ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE",
    "LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL",
}

ALLOWED_POLYMER_RESNAMES = set(STANDARD_RESNAMES_20)

# Fallback heuristic PTM residue names (used when explicit categories are absent)
COMMON_PTM_RESNAMES = {
    "MSE",                 # selenomethionine
    "SEP", "TPO", "PTR",   # phospho Ser/Thr/Tyr
    "CSO", "CME", "CSD", "OCS",
    "HYP",
    "KCX",
    "MLY", "M2L", "M3L",
}

# If you want “PTM” only (exclude PCMs), you can later filter by ptm_category values
# (e.g., keep only rows where ptm_categories contains "PTM").
# For now, "has_ptm" means “has protein modification (PCM/PTM) observed”.
# ------------------------------------------------------------
# 3-letter -> 1-letter (for common residues)
AA3_TO_1 = {
    "ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E","GLY":"G","HIS":"H","ILE":"I",
    "LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P","SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V",
    "SEC":"U",
}

# PTM residues you may see in atom_site; what the "original" letter should be
PTM_RESNAME_TO_BASE = {
    "PTR": "Y",
    "TPO": "T",
    "SEP": "S",
    "MSE": "M",   # selenomethionine; you might or might not want to treat as PTM
}

def _seq_is_clean(
    seq: str,
    *,
    allowed_aas: Optional[Set[str]] = None,
) -> Tuple[bool, Optional[str]]:
    """
    Returns (ok, reason_if_not_ok)
    """
    allowed = allowed_aas if allowed_aas is not None else VALID_AAS

    seq = re.sub(r"\s+", "", seq).upper()
    if not seq:
        return False, "no_sequence"

    bad = sorted(set(seq) - VALID_AAS)
    if bad:
        # includes any letters besides allowed set (including B, Z, J, etc.)
        return False, f"contains_non_allowed_letters:{''.join(bad)}"

    return True, None

def clean_mmcif_seq_text(x: Optional[str]) -> str:
    """
    Clean an mmCIF polymer sequence field (e.g. _entity_poly.pdbx_seq_one_letter_code_can).

    Handles:
      - semicolon-delimited multi-line text blocks
      - embedded newlines/spaces
      - stray leading/trailing semicolons produced by some parsers
    """
    if x is None:
        return ""
    s = str(x)

    # Trim outer whitespace first
    s = s.strip()

    # If the *value* includes mmCIF semicolon block markers, remove them.
    # In practice (e.g., via gemmi), the returned string often literally begins/ends with ';'.
    if s.startswith(";"):
        s = s[1:]  # drop the opening ';'
    # drop a trailing ';' if present after stripping
    s = s.strip()
    if s.endswith(";"):
        s = s[:-1]

    # Remove all whitespace/newlines inside the sequence
    s = re.sub(r"\s+", "", s)

    return s

def expand_allowed_tokens(
    extra_tokens: Optional[Set[str]]
) -> tuple[Set[str], Set[str]]:
    """
    Given user-specified tokens (1-letter or 3-letter),
    return:
      (extra_valid_aas, extra_allowed_resnames)

    Example:
      {"U"}   -> ({"U"}, {"SEC"})
      {"SEC"} -> ({"U"}, {"SEC"})
    """
    if not extra_tokens:
        return set(), set()

    extra_valid = set()
    extra_resnames = set()

    for t in extra_tokens:
        if not t:
            continue
        s = str(t).strip().upper()

        if len(s) == 1:
            # one-letter code
            extra_valid.add(s)
            # try to infer 3-letter
            for k, v in AA3_TO_1.items():
                if v == s:
                    extra_resnames.add(k)

        elif len(s) == 3:
            # three-letter code
            extra_resnames.add(s)
            if s in AA3_TO_1:
                extra_valid.add(AA3_TO_1[s])

    return extra_valid, extra_resnames