# Helper to split a list of items into roughly equal-sized chunks for parallel processing.

# chunkify: Divides a list into n_chunks contiguous sublists of approximately equal size.

import math
from typing import List

def chunk_by_size(items: List[str], chunk_size: int) -> List[List[str]]:
    """
    Convert a desired chunk_size (items per chunk) into n_chunks for chunkify().
    We keep chunkify() untouched (it expects n_chunks, not chunk_size).
    """
    chunk_size = max(1, int(chunk_size))
    n = len(items)
    if n == 0:
        return []
    n_chunks = max(1, math.ceil(n / chunk_size))
    return chunkify(items, n_chunks)


def chunkify(items: List[str], n_chunks: int) -> List[List[str]]:
    """Split items into n_chunks as evenly as possible."""
    n_chunks = max(1, int(n_chunks))
    n = len(items)
    if n == 0:
        return [[] for _ in range(n_chunks)]
    # ceil(n / n_chunks) sized blocks
    k = math.ceil(n / n_chunks)
    return [items[i:i+k] for i in range(0, n, k)]