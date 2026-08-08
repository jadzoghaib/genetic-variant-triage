"""Unit tests for the structural confidence tiers. No database required."""

import pandas as pd
import pytest

from core import confidence as C


@pytest.mark.parametrize("plddt,expected", [
    (0.0, "very_low"), (49.9, "very_low"),
    (50.0, "low"), (69.9, "low"),
    (70.0, "confident"), (89.9, "confident"),
    (90.0, "very_high"), (100.0, "very_high"),
])
def test_plddt_band_boundaries(plddt, expected):
    assert C.plddt_band(plddt) == expected


@pytest.mark.parametrize("plddt,expected", [
    (69.9, C.PREDICTED_WEAK),
    (70.0, C.PREDICTED_CONFIDENT),
    (98.6, C.PREDICTED_CONFIDENT),
])
def test_unsolved_residues_fall_back_to_plddt(plddt, expected):
    assert C.evidence_tier(False, plddt) == expected


def test_experimental_coverage_overrides_low_plddt():
    """The EGFR L858R case, which is why this tier exists.

    L858R is the most important activating mutation in non-small-cell lung
    cancer and sits at pLDDT 51.2 — AlphaFold models the activation loop poorly
    because it is genuinely flexible. Crystal structures cover it. A
    pLDDT-only gate would have down-weighted the most clinically consequential
    variant in the dataset.
    """
    assert C.evidence_tier(True, 51.2) == C.EXPERIMENTAL
    assert C.evidence_tier(True, 0.0) == C.EXPERIMENTAL
    assert C.evidence_tier(False, 51.2) == C.PREDICTED_WEAK


def test_missing_plddt_is_not_treated_as_confident():
    assert C.evidence_tier(False, None) == C.PREDICTED_WEAK


def test_only_solved_and_confident_support_a_structural_claim():
    assert C.supports_structural_claim(C.EXPERIMENTAL)
    assert C.supports_structural_claim(C.PREDICTED_CONFIDENT)
    assert not C.supports_structural_claim(C.PREDICTED_WEAK)


def test_tiers_are_ordered_best_first():
    assert C.TIER_ORDER == (C.EXPERIMENTAL, C.PREDICTED_CONFIDENT, C.PREDICTED_WEAK)


def test_assign_over_a_dataframe():
    df = pd.DataFrame({"is_solved": [True, False, False], "plddt": [51.2, 95.0, 30.0]})
    out = C.assign(df)
    assert list(out["evidence_tier"]) == [
        C.EXPERIMENTAL, C.PREDICTED_CONFIDENT, C.PREDICTED_WEAK]
    assert list(out["plddt_band"]) == ["low", "very_high", "very_low"]
