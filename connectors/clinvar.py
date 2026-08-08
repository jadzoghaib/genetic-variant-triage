"""ClinVar connector — domain parsing on top of the tabix range fetcher.

Mechanics (index parsing, BGZF ranges) live in clinvar_regions.py; this module
knows what a ClinVar record means.
"""

from __future__ import annotations

import pandas as pd

import clinvar_regions as cr

SOURCE = "clinvar"

# ClinVar review status -> star rating.
STARS = {
    "practice_guideline": 4,
    "reviewed_by_expert_panel": 3,
    "criteria_provided,_multiple_submitters,_no_conflicts": 2,
    "criteria_provided,_single_submitter": 1,
    "criteria_provided,_conflicting_classifications": 1,
    "criteria_provided,_conflicting_interpretations": 1,
}


def parse_info(info: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for field in info.split(";"):
        if "=" in field:
            k, _, v = field.partition("=")
            out[k] = v
    return out


def collapse_significance(s: str) -> str:
    """Collapse ClinVar significance into decision-relevant buckets.

    Order matters: 'Conflicting' and 'Uncertain' are tested before the
    Pathogenic/Benign checks so compound labels land correctly.

    Two distinctions this draws that an earlier version did not, each of which
    was silently corrupting the triage input:

      * An absent or `not_provided` CLNSIG is NO_ASSERTION — the ClinVar record
        exists but makes no pathogenicity claim. That is not the same as OTHER,
        and it must triage alongside variants ClinVar has never seen.
      * Membership tests miss compound labels. `Pathogenic/Likely_pathogenic/
        Pathogenic,_low_penetrance` is pathogenic, but it is in no fixed set, so
        prefix matching is used instead.

    OTHER is reserved for real assertions on a different axis — drug_response,
    risk_factor, protective — which are not pathogenicity calls at all.
    """
    if not s or s == "not_provided":
        return "NO_ASSERTION"
    if "Conflicting" in s:
        return "CONFLICTING"
    if "Uncertain" in s or "no_classification" in s:
        return "VUS"
    if s.startswith(("Pathogenic", "Likely_pathogenic")):
        return "PATH"
    if s.startswith(("Benign", "Likely_benign")):
        return "BENIGN"
    return "OTHER"


def load_index() -> cr.TabixIndex:
    return cr.load_index()


def fetch_gene(
    index: cr.TabixIndex, symbol: str, chrom: str, beg: int, end: int
) -> tuple[pd.DataFrame, bytes, str]:
    """ClinVar records attributed to `symbol` within a genomic span.

    The span comes from the target's own AlphaMissense file, so both sides of
    the join describe the same region by construction. The span can still
    overlap neighbouring genes (NBR2 overlaps BRCA1), so GENEINFO is checked
    per record.

    Duplicate genomic coordinates do occur — the same substitution can appear
    under more than one ClinVar record. They are collapsed to the best-reviewed
    assertion, since clinical_assertion is keyed one-per-variant.
    """
    lines = cr.fetch_region(index, chrom, beg, end)
    raw = "\n".join(lines).encode("utf-8")
    url = f"{cr.CLINVAR_VCF}#{chrom}:{beg}-{end}"

    rows = []
    for line in lines:
        f = line.split("\t")
        if len(f) < 8:
            continue
        c, pos, vid, ref, alt, _qual, _filt, info = f[:8]
        d = parse_info(info)
        genes = {g.split(":")[0] for g in d.get("GENEINFO", "").split("|") if g}
        if symbol not in genes:
            continue
        sig_raw = d.get("CLNSIG", "")
        revstat = d.get("CLNREVSTAT", "")
        for a in alt.split(","):
            rows.append({
                "chrom": c,
                "pos": int(pos),
                "ref": ref,
                "alt": a,
                "variant_id": f"{c}-{pos}-{ref}-{a}",
                "clinvar_id": vid,
                "significance_raw": sig_raw,
                "significance": collapse_significance(sig_raw),
                "review_status": revstat,
                "stars": STARS.get(revstat, 0),
                # Aggregated across ALL transcripts — kept for lineage only.
                # Never predicate on it; see SPEC.md Phase 0b.
                "molecular_consequence": d.get("MC", ""),
                "is_snv": len(ref) == 1 and len(a) == 1,
            })

    df = pd.DataFrame(rows)
    if df.empty:
        return df, raw, url

    before = len(df)
    df = (df.sort_values("stars", ascending=False)
            .drop_duplicates(subset="variant_id", keep="first")
            .reset_index(drop=True))
    df.attrs["duplicates_collapsed"] = before - len(df)
    return df, raw, url
