"""Unit tests for the triage matrix. No database required."""

from itertools import product

import pandas as pd
import pytest

from core import triage as T


# ── totality ──────────────────────────────────────────────────────────────

def test_matrix_covers_every_combination():
    """No (significance, prediction) pair may fall through to a default."""
    expected = set(product(T.SIGNIFICANCES, T.AM_CLASSES))
    assert set(T.MATRIX) == expected


def test_every_matrix_value_is_a_declared_class():
    declared = {
        T.UPGRADE, T.DOWNGRADE, T.REMAINS_UNCERTAIN, T.CONCORDANT, T.DISCORDANT,
        T.MODEL_UNINFORMATIVE, T.NOVEL_CANDIDATE, T.UNASSERTED, T.NOT_TRIAGED,
    }
    assert set(T.MATRIX.values()) <= declared


# ── the documented matrix, cell by cell ───────────────────────────────────

@pytest.mark.parametrize("significance,am_class,expected", [
    # uncertain: the product's reason for existing
    ("VUS",         "LPath", T.UPGRADE),
    ("VUS",         "LBen",  T.DOWNGRADE),
    ("VUS",         "Amb",   T.REMAINS_UNCERTAIN),
    ("CONFLICTING", "LPath", T.UPGRADE),
    ("CONFLICTING", "LBen",  T.DOWNGRADE),
    ("CONFLICTING", "Amb",   T.REMAINS_UNCERTAIN),
    # asserted: agreement and contradiction
    ("PATH",   "LPath", T.CONCORDANT),
    ("PATH",   "LBen",  T.DISCORDANT),
    ("PATH",   "Amb",   T.MODEL_UNINFORMATIVE),
    ("BENIGN", "LBen",  T.CONCORDANT),
    ("BENIGN", "LPath", T.DISCORDANT),
    ("BENIGN", "Amb",   T.MODEL_UNINFORMATIVE),
    # no pathogenicity claim
    (None,     "LPath", T.NOVEL_CANDIDATE),
    (None,     "LBen",  T.UNASSERTED),
    (None,     "Amb",   T.UNASSERTED),
    # an assertion on a different axis is not a pathogenicity call
    ("OTHER",  "LPath", T.NOT_TRIAGED),
    ("OTHER",  "LBen",  T.NOT_TRIAGED),
])
def test_classify(significance, am_class, expected):
    assert T.classify(significance, am_class) == expected


def test_only_actionable_classes_reach_a_worklist():
    assert T.ACTIONABLE == {T.UPGRADE, T.DOWNGRADE, T.DISCORDANT}
    assert T.CONCORDANT not in T.ACTIONABLE
    assert T.NOT_TRIAGED not in T.ACTIONABLE


# ── absence, in all its spellings ─────────────────────────────────────────

@pytest.mark.parametrize("value", [None, "", "ABSENT", "NO_ASSERTION", float("nan")])
def test_all_absence_spellings_normalise_identically(value):
    assert T.normalise_significance(value) == T.ABSENT


def test_no_assertion_triages_as_absent_not_as_other():
    """An empty CLNSIG means ClinVar made no claim. Treating it as OTHER would
    silently drop 1,114 real variants out of triage."""
    assert T.classify("NO_ASSERTION", "LPath") == T.NOVEL_CANDIDATE
    assert T.classify("OTHER", "LPath") == T.NOT_TRIAGED


# ── input validation ──────────────────────────────────────────────────────

def test_unknown_significance_raises():
    with pytest.raises(ValueError, match="unknown clinical significance"):
        T.classify("Likely_pathogenic", "LPath")   # raw ClinVar label, not collapsed


def test_missing_prediction_raises():
    """No prediction means the variant is not missense on the canonical
    transcript; asking for its triage class is a caller bug."""
    for bad in (None, "", "likely_pathogenic"):
        with pytest.raises(ValueError, match="unknown AlphaMissense class"):
            T.classify("VUS", bad)


# ── model strength ────────────────────────────────────────────────────────

@pytest.mark.parametrize("score,expected", [
    (1.0, "strong"), (0.9, "strong"), (0.1, "strong"), (0.0, "strong"),
    (0.89, "moderate"), (0.564, "moderate"), (0.34, "moderate"), (0.11, "moderate"),
    (0.5, "ambiguous"), (0.4, "ambiguous"),
])
def test_model_strength_boundaries(score, expected):
    assert T.model_strength(score) == expected


# ── review priority ───────────────────────────────────────────────────────

def test_non_actionable_classes_get_no_priority():
    priority, reasons = T.review_priority(T.CONCORDANT, "experimental", 4, 0.99)
    assert priority == "none" and reasons == []


def test_high_priority_requires_structure_stars_and_a_strong_call():
    assert T.review_priority(T.UPGRADE, "experimental", 2, 0.99)[0] == "high"
    # drop any one requirement and it is no longer high
    assert T.review_priority(T.UPGRADE, "predicted_confident", 2, 0.99)[0] == "medium"
    assert T.review_priority(T.UPGRADE, "experimental", 1, 0.99)[0] == "medium"
    assert T.review_priority(T.UPGRADE, "experimental", 2, 0.80)[0] == "medium"


def test_weak_structure_caps_priority_at_low():
    assert T.review_priority(T.UPGRADE, "predicted_weak", 4, 0.99)[0] == "low"


def test_priority_reasons_are_reported():
    _, reasons = T.review_priority(T.DISCORDANT, "experimental", 3, 0.97)
    assert any("class=" in r for r in reasons)
    assert any("structure=" in r for r in reasons)
    assert any("clinvar=" in r for r in reasons)


# ── dataframe application ─────────────────────────────────────────────────

def test_assign_over_a_dataframe():
    df = pd.DataFrame({
        "significance": ["VUS", "PATH", None, "OTHER"],
        "am_class": ["LPath", "LBen", "LPath", "Amb"],
        "evidence_tier": ["experimental", "experimental", "predicted_weak", "experimental"],
        "stars": [2, 3, float("nan"), 1],
        "am_pathogenicity": [0.97, 0.05, 0.95, 0.5],
    })
    out = T.assign(df)
    assert list(out["triage_class"]) == [
        T.UPGRADE, T.DISCORDANT, T.NOVEL_CANDIDATE, T.NOT_TRIAGED]
    assert list(out["review_priority"]) == ["high", "high", "none", "none"]
    assert out["triage_class"].notna().all()
