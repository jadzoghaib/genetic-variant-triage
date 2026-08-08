"""The triage matrix — reconciling a clinical assertion with a model prediction.

This is the product. Everything else is plumbing around it.

The matrix is total by construction: it is built over the full cross product of
significance x prediction, so no combination can silently fall through to a
default. `test_triage.py` asserts that totality.

A deliberate constraint: Locus surfaces *candidates for expert review*. It never
asserts a clinical classification. The class names say what was observed
(a model and a curator disagree), not what should be done about it.
"""

from __future__ import annotations

from itertools import product

from core._util import is_missing

# ── triage classes ────────────────────────────────────────────────────────
UPGRADE = "reclass_upgrade"            # uncertain, model says pathogenic
DOWNGRADE = "reclass_downgrade"        # uncertain, model says benign
REMAINS_UNCERTAIN = "remains_uncertain"
CONCORDANT = "concordant"
DISCORDANT = "discordant"              # model contradicts a confident curator
MODEL_UNINFORMATIVE = "model_uninformative"
NOVEL_CANDIDATE = "novel_candidate"    # no assertion, model says pathogenic
UNASSERTED = "unasserted"
NOT_TRIAGED = "not_triaged"            # assertion is on a non-pathogenicity axis

#: Classes that put a variant on someone's worklist.
ACTIONABLE = frozenset({UPGRADE, DOWNGRADE, DISCORDANT})

# ── input vocabularies ────────────────────────────────────────────────────
ABSENT = "ABSENT"
#: Values that mean "ClinVar makes no pathogenicity claim" and must triage
#: identically to a variant ClinVar has never seen.
_ABSENT_ALIASES = frozenset({None, "", "ABSENT", "NO_ASSERTION"})

SIGNIFICANCES = ("PATH", "BENIGN", "VUS", "CONFLICTING", "OTHER", ABSENT)
AM_CLASSES = ("LPath", "Amb", "LBen")

_UNCERTAIN = {"LPath": UPGRADE, "Amb": REMAINS_UNCERTAIN, "LBen": DOWNGRADE}
_ASSERTED_PATH = {"LPath": CONCORDANT, "Amb": MODEL_UNINFORMATIVE, "LBen": DISCORDANT}
_ASSERTED_BENIGN = {"LPath": DISCORDANT, "Amb": MODEL_UNINFORMATIVE, "LBen": CONCORDANT}
_ABSENT_ROW = {"LPath": NOVEL_CANDIDATE, "Amb": UNASSERTED, "LBen": UNASSERTED}

_ROWS = {
    "VUS": _UNCERTAIN,
    "CONFLICTING": _UNCERTAIN,
    "PATH": _ASSERTED_PATH,
    "BENIGN": _ASSERTED_BENIGN,
    "OTHER": dict.fromkeys(AM_CLASSES, NOT_TRIAGED),
    ABSENT: _ABSENT_ROW,
}

#: The full matrix, materialised so gaps are impossible rather than unlikely.
MATRIX: dict[tuple[str, str], str] = {
    (sig, am): _ROWS[sig][am] for sig, am in product(SIGNIFICANCES, AM_CLASSES)
}


def normalise_significance(significance: str | None) -> str:
    """Map every 'ClinVar has no pathogenicity opinion' spelling onto ABSENT.

    Missing values are included: a LEFT JOIN that found no assertion surfaces
    as NaN or pandas NA once the frame round-trips through pandas.
    """
    if is_missing(significance) or significance in _ABSENT_ALIASES:
        return ABSENT
    if significance not in SIGNIFICANCES:
        raise ValueError(f"unknown clinical significance: {significance!r}")
    return significance


def classify(significance: str | None, am_class: str | None) -> str:
    """Triage class for one variant.

    A missing prediction is not a triage outcome — it means the variant is not
    missense on the canonical transcript, and the caller should not have asked.
    """
    if am_class not in AM_CLASSES:
        raise ValueError(f"unknown AlphaMissense class: {am_class!r}")
    return MATRIX[(normalise_significance(significance), am_class)]


def model_strength(am_pathogenicity: float) -> str:
    """How far the model's score sits from its own ambiguous band.

    The published thresholds are 0.34 and 0.564; 'strong' additionally requires
    the score to be near a rail, where the model is least likely to be
    borderline by accident.
    """
    p = am_pathogenicity
    if p >= 0.9 or p <= 0.1:
        return "strong"
    if p >= 0.564 or p <= 0.34:
        return "moderate"
    return "ambiguous"


def review_priority(
    triage_class: str,
    evidence_tier: str,
    stars: int | None,
    am_pathogenicity: float,
) -> tuple[str, list[str]]:
    """Worklist priority, with the reasons that produced it.

    Rule-based rather than a composite score: an analyst must be able to see
    why a variant reached the top of the queue, and a weighted sum cannot be
    argued with.
    """
    if triage_class not in ACTIONABLE:
        return "none", []

    from core.confidence import EXPERIMENTAL, supports_structural_claim

    stars = stars or 0
    strength = model_strength(am_pathogenicity)
    reasons = [f"class={triage_class}", f"structure={evidence_tier}",
               f"clinvar={stars}star", f"model={strength}"]

    if evidence_tier == EXPERIMENTAL and stars >= 2 and strength == "strong":
        return "high", reasons
    if supports_structural_claim(evidence_tier) and stars >= 1 and strength != "ambiguous":
        return "medium", reasons
    return "low", reasons


def assign(df):
    """Add `triage_class`, `model_strength`, `review_priority`, `priority_reasons`.

    Expects `significance`, `am_class`, `evidence_tier`, `stars`,
    `am_pathogenicity`.
    """
    out = df.copy()
    out["triage_class"] = [
        classify(s, a) for s, a in zip(out["significance"], out["am_class"])
    ]
    out["model_strength"] = [model_strength(p) for p in out["am_pathogenicity"]]

    priorities, reasons = [], []
    for tc, tier, stars, p in zip(out["triage_class"], out["evidence_tier"],
                                  out["stars"], out["am_pathogenicity"]):
        pr, rs = review_priority(tc, tier, None if is_missing(stars) else int(stars), p)
        priorities.append(pr)
        reasons.append("; ".join(rs))
    out["review_priority"] = priorities
    out["priority_reasons"] = reasons
    return out
