"""Integration tests — the core rules applied to the real store.

These are the credibility tests. The unit tests prove the rules are internally
consistent; these prove they produce the right answer on variants whose status
is established independently of this pipeline.
"""

import pytest

import queries
from core import confidence as C
from core import dossier as D
from core import triage as T

# Phase 0 spike baseline, measured 2026-08-07 before the schema existed, scoped
# to the three genes it ran on. Treated as a reference point rather than a
# fixture: the pipeline now rebuilds weekly against live ClinVar, so these
# numbers drift as variants are reclassified. See
# test_phase0_baseline_has_not_drifted_implausibly.
PHASE0_GENES = ("BRCA1", "TP53", "PTEN")
PHASE0 = {"uncertain": 6488, "upgrades": 1308, "downgrades": 4807}

# BRCA1 variants whose classification is established outside this pipeline.
# The pathogenic set are the canonical RING (C61G) and BRCT domain missense
# variants; the benign set are the common population polymorphisms.
BRCA1_PATHOGENIC = ["C61G", "T1691K", "R1699Q", "R1699W", "V1736A", "M1775R"]
BRCA1_BENIGN = ["P871L", "E1038G", "K1183R", "S1613G"]


@pytest.fixture(scope="module")
def triaged(con):
    df = queries.worklist(con)
    return T.assign(C.assign(df))


# ── totality on real data ─────────────────────────────────────────────────

def test_every_prediction_receives_a_triage_class(triaged):
    assert len(triaged) > 0
    assert triaged["triage_class"].notna().all()
    assert triaged["evidence_tier"].notna().all()
    assert set(triaged["triage_class"]) <= set(T.MATRIX.values())


def test_class_counts_partition_the_dataset(triaged):
    assert triaged["triage_class"].value_counts().sum() == len(triaged)


def test_priority_is_assigned_exactly_to_actionable_classes(triaged):
    actionable = triaged["triage_class"].isin(T.ACTIONABLE)
    has_priority = triaged["review_priority"] != "none"
    assert (actionable == has_priority).all()


# ── external ground truth ─────────────────────────────────────────────────

def test_known_pathogenic_brca1_variants_are_concordant(triaged):
    """All six are expert-reviewed pathogenic and all should be confirmed by the
    model — if any came back DISCORDANT the premise would be in doubt."""
    rows = triaged[(triaged.symbol == "BRCA1")
                   & triaged.protein_variant.isin(BRCA1_PATHOGENIC)]
    assert set(rows.protein_variant) == set(BRCA1_PATHOGENIC)
    assert (rows.significance == "PATH").all()
    assert (rows.am_class == "LPath").all()
    assert (rows.triage_class == T.CONCORDANT).all()


def test_known_benign_brca1_polymorphisms_are_concordant(triaged):
    rows = triaged[(triaged.symbol == "BRCA1")
                   & triaged.protein_variant.isin(BRCA1_BENIGN)]
    assert set(rows.protein_variant) == set(BRCA1_BENIGN)
    assert (rows.significance == "BENIGN").all()
    assert (rows.am_class == "LBen").all()
    assert (rows.triage_class == T.CONCORDANT).all()


def test_brca1_pathogenic_sit_in_solved_domains_and_benign_do_not(triaged):
    """The structural signal is real, not decorative: BRCA1's pathogenic missense
    variants cluster in the folded RING/BRCT domains that crystallise, while the
    common polymorphisms sit in the disordered central linker."""
    brca1 = triaged[triaged.symbol == "BRCA1"]
    path = brca1[brca1.protein_variant.isin(BRCA1_PATHOGENIC)]
    benign = brca1[brca1.protein_variant.isin(BRCA1_BENIGN)]
    assert path.is_solved.all()
    assert not benign.is_solved.any()
    assert path.plddt.min() > benign.plddt.max()


def test_egfr_driver_mutations_are_all_predicted_pathogenic(triaged):
    rows = triaged[(triaged.symbol == "EGFR")
                   & triaged.protein_variant.isin(["L858R", "T790M", "G719S", "L861Q"])]
    assert len(rows) >= 4
    assert (rows.am_class == "LPath").all()


def test_l858r_is_rescued_by_experimental_coverage(triaged):
    """The case that forced the three-tier model into existence."""
    row = triaged[(triaged.symbol == "EGFR")
                  & (triaged.protein_variant == "L858R")].iloc[0]
    assert row.plddt < C.PLDDT_CONFIDENT          # prediction alone is unreliable
    assert row.is_solved                          # but crystal structures cover it
    assert row.evidence_tier == C.EXPERIMENTAL


# ── regression against the Phase 0 spike ──────────────────────────────────

def _phase0_slice(triaged):
    return triaged[triaged.symbol.isin(PHASE0_GENES)
                   & triaged.significance.isin(["VUS", "CONFLICTING"])]


def test_uncertain_variants_are_partitioned_exactly(triaged):
    """Every uncertain variant lands in exactly one of three outcomes.

    This is the property the Phase 0 spike actually established — that the
    normalised schema reproduces the flat spike's logic — and unlike a count it
    is true no matter what ClinVar does next.
    """
    df = _phase0_slice(triaged)
    up = int((df.triage_class == T.UPGRADE).sum())
    down = int((df.triage_class == T.DOWNGRADE).sum())
    same = int((df.triage_class == T.REMAINS_UNCERTAIN).sum())
    assert up + down + same == len(df), "an uncertain variant fell outside the matrix"
    assert len(df) > 0


def test_phase0_baseline_has_not_drifted_implausibly(triaged):
    """Guard the pipeline, not the upstream data.

    The Phase 0 figures were measured against a fixed 2026-08-07 snapshot, and
    asserting them exactly used to be the regression check. That stopped being
    tenable once the site rebuilt weekly against live ClinVar: on 2026-08-10 the
    build failed because ClinVar had reclassified two variants out of
    "uncertain", one of them an upgrade candidate — which is precisely the event
    this tool exists to surface. A test that breaks when the data does its job
    is testing the wrong thing.

    So this checks the shape instead: a few percent of movement is ClinVar
    working, while a large swing means the join, the transcript mapping or a
    connector has broken.
    """
    df = _phase0_slice(triaged)
    got = {
        "uncertain": len(df),
        "upgrades": int((df.triage_class == T.UPGRADE).sum()),
        "downgrades": int((df.triage_class == T.DOWNGRADE).sum()),
    }
    for key, baseline in PHASE0.items():
        drift = abs(got[key] - baseline) / baseline
        assert drift < 0.10, (
            f"{key} moved {drift:.1%} from the 2026-08-07 baseline "
            f"({baseline} -> {got[key]}). Reclassification is normal; a swing "
            f"this large is a pipeline fault."
        )


def test_no_assertion_variants_are_triaged_not_discarded(triaged):
    """Empty CLNSIG previously collapsed to OTHER and dropped out of triage."""
    absent = triaged[triaged.significance.isna()
                     | (triaged.significance == "NO_ASSERTION")]
    assert len(absent) > 0
    assert set(absent.triage_class) <= {T.NOVEL_CANDIDATE, T.UNASSERTED}


# ── dossier ───────────────────────────────────────────────────────────────

@pytest.fixture(scope="module")
def cards(con):
    facts = queries.target_facts(con)
    out = {}
    for _, row in facts.iterrows():
        f = row.to_dict()
        f["drug_stages"] = queries.drug_stages(con, f["uniprot_acc"])
        out[f["symbol"]] = D.build(f)
    return out


def test_every_loaded_target_has_strong_genetic_evidence(cards):
    """The gene set is chosen from established disease genes, so a weak band
    here means the Open Targets join went wrong, not that the biology is
    marginal."""
    assert len(cards) >= 4
    for symbol, card in cards.items():
        assert card["dimensions"]["genetic_evidence"]["band"] == D.STRONG, symbol


def test_dossier_separates_a_druggable_kinase_from_tumour_suppressors(cards):
    assert cards["EGFR"]["archetype"] == D.VALIDATED_DRUGGABLE
    for symbol in ("BRCA1", "PTEN"):
        card = cards[symbol]
        assert card["archetype"] == D.GENETICALLY_VALIDATED_UNDRUGGED
        assert card["dimensions"]["binding_site"]["band"] == D.ABSENT
        assert card["dimensions"]["chemical_matter"]["n_drugs"] == 0


def test_structural_readiness_reflects_experimental_coverage(cards):
    assert cards["TP53"]["dimensions"]["structural_readiness"]["band"] == D.STRONG
    assert cards["BRCA1"]["dimensions"]["structural_readiness"]["band"] == D.WEAK
