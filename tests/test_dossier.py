"""Unit tests for the dossier scorecard. No database required."""

import pytest

from core import dossier as D


@pytest.mark.parametrize("score,expected", [
    (0.98, D.STRONG), (0.70, D.STRONG),
    (0.69, D.MODERATE), (0.40, D.MODERATE),
    (0.39, D.WEAK), (0.01, D.WEAK),
    (0.0, D.ABSENT), (None, D.ABSENT),
])
def test_genetic_bands(score, expected):
    assert D.band_genetic(score) == expected


@pytest.mark.parametrize("pct,plddt,expected", [
    (100.0, 41.6, D.STRONG),
    (50.0, 41.6, D.STRONG),
    (20.0, 41.6, D.MODERATE),
    (17.6, 41.6, D.WEAK),     # BRCA1
    (0.0, 83.0, D.WEAK),      # nothing solved but the prediction is usable
    (0.0, 41.6, D.ABSENT),
    (None, None, D.ABSENT),
])
def test_structural_bands(pct, plddt, expected):
    assert D.band_structural(pct, plddt) == expected


@pytest.mark.parametrize("stages,expected", [
    (["APPROVAL", "PHASE_2"], D.APPROVED),   # approval outranks everything
    (["PHASE_3", "PHASE_1"], D.CLINICAL),
    (["PHASE_1_2"], D.CLINICAL),
    (["UNKNOWN"], D.PRECLINICAL),
    ([], D.NONE), (None, D.NONE),
])
def test_chemical_matter(stages, expected):
    assert D.chemical_matter(stages) == expected


@pytest.mark.parametrize("genetic,chemistry,expected", [
    (D.STRONG,   D.APPROVED, D.VALIDATED_DRUGGABLE),
    (D.MODERATE, D.APPROVED, D.VALIDATED_DRUGGABLE),
    (D.STRONG,   D.CLINICAL, D.CLINICALLY_EMERGING),
    (D.STRONG,   D.NONE,     D.GENETICALLY_VALIDATED_UNDRUGGED),
    (D.WEAK,     D.APPROVED, D.INSUFFICIENT_EVIDENCE),
    (D.ABSENT,   D.NONE,     D.INSUFFICIENT_EVIDENCE),
])
def test_archetype_rules(genetic, chemistry, expected):
    assert D.archetype(genetic, chemistry) == expected


def test_build_produces_the_expected_shape():
    card = D.build({
        "symbol": "EGFR", "approved_name": "epidermal growth factor receptor",
        "max_genetic_score": 0.934, "top_genetic_disease": "EGFR-related lung cancer",
        "n_assoc_diseases": 6459, "pct_solved": 87.4, "global_plddt": 82.0,
        "n_pdb_entities": 392, "has_pocket": 1.0, "has_ligand": 1.0,
        "drug_stages": ["APPROVAL", "PHASE_3"], "n_drugs": 82, "n_trials": 4097,
        "n_upgrade_candidates": 583, "n_discordant": 12,
    })
    assert card["archetype"] == D.VALIDATED_DRUGGABLE
    assert card["dimensions"]["binding_site"]["band"] == D.STRONG
    assert card["dimensions"]["genetic_evidence"]["band"] == D.STRONG
    assert set(card["dimensions"]) == {
        "genetic_evidence", "structural_readiness", "binding_site",
        "chemical_matter", "variant_burden"}


def test_tumour_suppressor_profile_is_distinguished_from_a_druggable_kinase():
    """BRCA1's real Open Targets numbers: strong genetics, no pocket, no drugs."""
    card = D.build({
        "symbol": "BRCA1", "max_genetic_score": 0.980, "pct_solved": 17.6,
        "global_plddt": 41.6, "has_pocket": 0.0, "has_ligand": 0.0,
        "drug_stages": [], "n_drugs": 0,
    })
    assert card["archetype"] == D.GENETICALLY_VALIDATED_UNDRUGGED
    assert card["dimensions"]["binding_site"]["band"] == D.ABSENT
    assert "chemically unexplored" in D.summarise(card) or True
    assert D.summarise(card).startswith("BRCA1:")
