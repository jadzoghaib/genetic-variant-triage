"""Small shared helpers for the core layer."""

from __future__ import annotations

import pandas as pd


def is_missing(value) -> bool:
    """True for None, float NaN, and pandas NA.

    The `value != value` NaN idiom is not enough here: a LEFT JOIN on an integer
    column yields pandas' nullable NA, whose comparisons return NA rather than
    True, and calling bool() on that raises.
    """
    try:
        return bool(pd.isna(value))
    except (TypeError, ValueError):     # arrays, unhashable objects
        return False
