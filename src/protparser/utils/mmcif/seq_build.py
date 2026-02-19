# Helpers for turning a canonical sequence + resolved residue positions into masked/resolved-only sequences
# and index mappings between canonical and resolved-only coordinates.

# build_resolved_sequences: Builds a masked sequence (unresolved = '-') plus the resolved-only sequence and coverage stats.
# build_resolved_index_map: Creates a mapping from canonical residue positions to resolved-only indices.
# map_canon_pos_to_resolved_only_indices: Converts canonical positions to resolved-only indices, dropping unresolved.
# _positions_to_csv: Serializes a set of positions as a CSV string (or NaN if empty).

from __future__ import annotations
from typing import Dict, Set, Any
import numpy as np

def build_resolved_sequences(full_seq: str, resolved_positions: set[int]) -> tuple[str, str, int, float]:
    """
    full_seq is the canonical polymer sequence (length L).
    resolved_positions are 1..L indices that have coordinates.
    """
    L = len(full_seq)
    masked = []
    only = []
    for i in range(1, L + 1):
        aa = full_seq[i - 1]
        if i in resolved_positions:
            masked.append(aa)
            only.append(aa)
        else:
            masked.append("-")
    n_res = len(resolved_positions)
    frac = (n_res / L) if L else 0.0
    return "".join(masked), "".join(only), n_res, frac

def build_resolved_index_map(L: int, resolved_positions: Set[int]) -> Dict[int, int]:
    """
    Map canonical 1-indexed position -> resolved_only 1-indexed position,
    for positions that are resolved.
    """
    m: Dict[int, int] = {}
    j = 0
    for i in range(1, L + 1):
        if i in resolved_positions:
            j += 1
            m[i] = j
    return m

def map_canon_pos_to_resolved_only_indices(pos: Set[int], canon_to_res: Dict[int, int]) -> Set[int]:
    """
    Convert canonical positions -> resolved_only indices (1-indexed), dropping positions not resolved.
    """
    return {canon_to_res[p] for p in pos if p in canon_to_res}

def _positions_to_csv(pos: Set[int]) -> Any:
    """
    Convert set of 1-indexed positions -> "1,2,30" or np.nan if empty.
    """
    if not pos:
        return np.nan
    return ",".join(str(i) for i in sorted(pos))