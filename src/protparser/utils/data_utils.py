import numpy as np
from typing import Any

def _norm_str(x: Any) -> str:
    if x is None:
        return ""
    if isinstance(x, float) and np.isnan(x):
        return ""
    s = str(x).strip()
    if s.lower() in ("", "nan", "none", "?", "."):
        return ""
    return s

def _pipe_sorted_unique(vals: list[str]) -> Any:
    vals = [_norm_str(v) for v in vals]
    vals = [v for v in vals if v]
    if not vals:
        return np.nan
    return "|".join(sorted(set(vals)))

def _csv_sorted_unique(vals: list[str]) -> Any:
    vals = [_norm_str(v) for v in vals]
    vals = [v for v in vals if v]
    if not vals:
        return np.nan
    return ",".join(sorted(set(vals)))