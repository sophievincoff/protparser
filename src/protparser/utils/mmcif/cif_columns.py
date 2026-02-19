# Safe helpers for accessing and normalizing mmCIF table columns and values.
#
# safe_find_column: Safely retrieves a column from a CIF category without raising errors.
# col_required: Retrieves a required column, raising a clear error if missing.
# cif_get_row_count: Returns the number of rows in a CIF table.
# truthy_cif_flag: Normalizes CIF yes/no-style flags to True/False/None.
# safe_int: Safely converts a CIF value to int, handling missing values.
# clean_str: Cleans CIF strings, removing placeholders like . and ?.

from __future__ import annotations
from typing import Any, Optional
import gemmi

def safe_find_column(cat: gemmi.cif.Table, name: str):
    """
    gemmi's find_column *sometimes raises* RuntimeError instead of returning None
    (depending on gemmi version/build). Make it safe.
    """
    if cat is None:
        return None
    try:
        return cat.find_column(name)
    except Exception:
        return None

def col_required(cat: gemmi.cif.Table, name: str):
    """
    Required column lookup (consistent error).
    """
    c = safe_find_column(cat, name)
    if c is None:
        raise KeyError(f"Missing mmCIF column '{name}' in category '{getattr(cat, 'name', '<unknown>')}'")
    return c

def cif_get_row_count(cat: gemmi.cif.Table) -> int:
    return len(cat) if hasattr(cat, "__len__") else cat.get_row_count()

def truthy_cif_flag(x: Any) -> Optional[bool]:
    """
    Normalize CIF yes/no-ish fields.
    Returns True/False, or None if can't interpret.
    """
    if x is None:
        return None
    s = str(x).strip().strip('"').strip("'").upper()
    if not s or s in {".", "?"}:
        return None
    if s in {"Y", "YES", "TRUE", "T", "1"}:
        return True
    if s in {"N", "NO", "FALSE", "F", "0"}:
        return False
    return None

def safe_int(x: Any) -> Optional[int]:
    if x is None:
        return None
    s = str(x).strip()
    if not s or s in {".", "?"}:
        return None
    try:
        return int(s)
    except Exception:
        return None
    
def clean_str(x: Any) -> str:
    """
    Normalize entity_id to a canonical string form.
    Examples:
      1   -> "1"
      "1" -> "1"
      " 1 " -> "1"
      None -> ""
    """
    if x is None:
        return ""
    s = str(x).strip()
    if not s or s in {".", "?"}:
        return ""
    # if numeric-like, normalize
    if s.isdigit():
        s = str(int(s))
    
    if s and ((s.startswith("'") and s.endswith("'")) or (s.startswith('"') and s.endswith('"'))):
        s = s[1:-1].strip()
    return s