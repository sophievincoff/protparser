# Pandas helper to standardize missing values by converting null-like entries to np.nan
# and normalizing column dtypes.

# harmonize_nulls_to_nan: Converts various null representations to np.nan and normalizes column dtypes.

from __future__ import annotations
import numpy as np
import pandas as pd

def harmonize_nulls_to_nan(df: pd.DataFrame, *, also_blank_strings=True, keep_datetime=False) -> pd.DataFrame:
    out = df.copy()
    if also_blank_strings:
        out = out.replace({"": pd.NA, "None": pd.NA, "nan": pd.NA})
    out = out.convert_dtypes()
    for c in out.columns:
        dt = out[c].dtype
        is_ext = isinstance(dt, pd.api.extensions.ExtensionDtype)
        if keep_datetime and pd.api.types.is_datetime64_any_dtype(dt):
            continue
        if is_ext:
            out[c] = out[c].astype(object)
    out = out.where(~out.isna(), np.nan)
    return out