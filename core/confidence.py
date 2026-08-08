"""Structural confidence — how much weight a structural interpretation can bear.

Phase 2 replaced a binary pLDDT gate with three ordered tiers, because the
binary version was wrong in a way that mattered: EGFR L858R, the most important
activating mutation in non-small-cell lung cancer, has a pLDDT of 51.2.
AlphaFold models the activation loop poorly because it is genuinely flexible —
but the residue is covered by crystal structures. A pLDDT-only gate would have
down-weighted the single most clinically consequential variant in the dataset.

Experimental coverage therefore OVERRIDES pLDDT. It is not averaged with it.
"""

from __future__ import annotations

EXPERIMENTAL = "experimental"
PREDICTED_CONFIDENT = "predicted_confident"
PREDICTED_WEAK = "predicted_weak"

TIER_ORDER = (EXPERIMENTAL, PREDICTED_CONFIDENT, PREDICTED_WEAK)

# AlphaFold's own banding thresholds.
PLDDT_VERY_HIGH = 90.0
PLDDT_CONFIDENT = 70.0
PLDDT_LOW = 50.0


def plddt_band(plddt: float) -> str:
    """AlphaFold's four-band vocabulary for a per-residue confidence score."""
    if plddt < PLDDT_LOW:
        return "very_low"
    if plddt < PLDDT_CONFIDENT:
        return "low"
    if plddt < PLDDT_VERY_HIGH:
        return "confident"
    return "very_high"


def evidence_tier(is_solved: bool, plddt: float | None) -> str:
    """Ordered structural evidence tier for one residue.

    `is_solved` means some experimental structure covers this position. It wins
    outright — a solved residue is solved regardless of what the prediction
    thought of it.
    """
    if is_solved:
        return EXPERIMENTAL
    if plddt is not None and plddt >= PLDDT_CONFIDENT:
        return PREDICTED_CONFIDENT
    return PREDICTED_WEAK


def supports_structural_claim(tier: str) -> bool:
    """Whether a structural interpretation of a variant at this residue is
    defensible at all. Weak predictions are shown, but never counted in
    headline figures."""
    return tier in (EXPERIMENTAL, PREDICTED_CONFIDENT)


def assign(df, solved_col: str = "is_solved", plddt_col: str = "plddt"):
    """Add `evidence_tier` and `plddt_band` columns to a dataframe."""
    from core._util import is_missing

    out = df.copy()
    out["evidence_tier"] = [
        evidence_tier(False if is_missing(s) else bool(s),
                      None if is_missing(p) else float(p))
        for s, p in zip(out[solved_col], out[plddt_col])
    ]
    out["plddt_band"] = [
        None if is_missing(p) else plddt_band(float(p)) for p in out[plddt_col]
    ]
    return out
