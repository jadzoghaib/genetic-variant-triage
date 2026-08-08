"""Target dossier scorecard.

Deliberately NOT a weighted composite score. A single number would be
unarguable and would hide which evidence carried it — the opposite of what a
target-selection decision needs. Instead each dimension is banded against a
stated threshold, and an archetype is derived from a small set of readable
rules. The raw evidence travels alongside every band.
"""

from __future__ import annotations

# ── dimension bands ───────────────────────────────────────────────────────
STRONG, MODERATE, WEAK, ABSENT = "strong", "moderate", "weak", "absent"

# ── chemical matter states ────────────────────────────────────────────────
APPROVED, CLINICAL, PRECLINICAL, NONE = "approved", "clinical", "preclinical", "none"

# ── archetypes ────────────────────────────────────────────────────────────
VALIDATED_DRUGGABLE = "validated druggable target"
CLINICALLY_EMERGING = "genetically validated, clinically emerging"
GENETICALLY_VALIDATED_UNDRUGGED = "genetically validated, chemically unexplored"
INSUFFICIENT_EVIDENCE = "insufficient genetic evidence"

CLINICAL_STAGES = {"PHASE_1", "PHASE_1_2", "PHASE_2", "PHASE_2_3", "PHASE_3"}


def band_genetic(max_genetic_score: float | None) -> str:
    if not max_genetic_score:
        return ABSENT
    if max_genetic_score >= 0.7:
        return STRONG
    if max_genetic_score >= 0.4:
        return MODERATE
    return WEAK


def band_structural(pct_solved: float | None, global_plddt: float | None) -> str:
    """Structural readiness: experimental coverage first, prediction as fallback."""
    pct = pct_solved or 0.0
    if pct >= 50:
        return STRONG
    if pct >= 20:
        return MODERATE
    if pct > 0:
        return WEAK
    if (global_plddt or 0) >= 70:
        return WEAK       # nothing solved, but the prediction is usable
    return ABSENT


def chemical_matter(stages: list[str] | None) -> str:
    stages = [s for s in (stages or []) if s]
    if not stages:
        return NONE
    if "APPROVAL" in stages:
        return APPROVED
    if any(s in CLINICAL_STAGES for s in stages):
        return CLINICAL
    return PRECLINICAL


def archetype(genetic: str, chemistry: str) -> str:
    if genetic in (WEAK, ABSENT):
        return INSUFFICIENT_EVIDENCE
    if chemistry == APPROVED:
        return VALIDATED_DRUGGABLE
    if chemistry in (CLINICAL, PRECLINICAL):
        return CLINICALLY_EMERGING
    return GENETICALLY_VALIDATED_UNDRUGGED


def build(facts: dict) -> dict:
    """Assemble one target's scorecard from already-fetched facts.

    `facts` is a plain dict so this stays free of any database dependency.
    """
    genetic = band_genetic(facts.get("max_genetic_score"))
    structural = band_structural(facts.get("pct_solved"), facts.get("global_plddt"))
    chemistry = chemical_matter(facts.get("drug_stages"))
    has_pocket = bool(facts.get("has_pocket"))
    has_ligand = bool(facts.get("has_ligand"))

    return {
        "symbol": facts.get("symbol"),
        "approved_name": facts.get("approved_name"),
        "archetype": archetype(genetic, chemistry),
        "dimensions": {
            "genetic_evidence": {
                "band": genetic,
                "max_genetic_score": facts.get("max_genetic_score"),
                "top_disease": facts.get("top_genetic_disease"),
                "n_associated_diseases": facts.get("n_assoc_diseases"),
            },
            "structural_readiness": {
                "band": structural,
                "pct_residues_solved": facts.get("pct_solved"),
                "n_pdb_entities": facts.get("n_pdb_entities"),
                "global_plddt": facts.get("global_plddt"),
            },
            "binding_site": {
                "band": STRONG if (has_pocket and has_ligand)
                        else MODERATE if (has_pocket or has_ligand) else ABSENT,
                "has_pocket": has_pocket,
                "has_ligand": has_ligand,
            },
            "chemical_matter": {
                "band": chemistry,
                "n_drugs": facts.get("n_drugs"),
                "n_trials": facts.get("n_trials"),
            },
            "variant_burden": {
                "upgrade_candidates": facts.get("n_upgrade_candidates"),
                "discordant": facts.get("n_discordant"),
            },
        },
    }


def summarise(card: dict) -> str:
    """One-line human summary. The archetype plus the evidence that decided it."""
    d = card["dimensions"]
    return (
        f"{card['symbol']}: {card['archetype']} "
        f"(genetic={d['genetic_evidence']['band']}, "
        f"structure={d['structural_readiness']['band']}, "
        f"pocket={d['binding_site']['band']}, "
        f"chemistry={d['chemical_matter']['band']})"
    )
