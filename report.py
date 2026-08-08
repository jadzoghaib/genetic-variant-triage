"""Locus analytical report — the core rules applied to the whole store.

This is what Phase 4's UI will render. Keeping it as a text report first proves
the analytics stand on their own before any pixels are involved.
"""

from __future__ import annotations

import pandas as pd

import queries
from core import confidence as C
from core import dossier as D
from core import triage as T

pd.set_option("display.width", 200)


def hr(title: str) -> None:
    print("\n" + "=" * 92)
    print(title)
    print("=" * 92)


def main() -> None:
    con = queries.connect()
    df = T.assign(C.assign(queries.worklist(con)))

    hr("1. TRIAGE CLASS DISTRIBUTION")
    pivot = (df.pivot_table(index="triage_class", columns="symbol",
                            values="variant_id", aggfunc="count", fill_value=0))
    pivot["TOTAL"] = pivot.sum(axis=1)
    print(pivot.sort_values("TOTAL", ascending=False).to_string())

    hr("2. REVIEW WORKLIST — actionable variants by priority")
    act = df[df.triage_class.isin(T.ACTIONABLE)]
    print(act.pivot_table(index="review_priority", columns="triage_class",
                          values="variant_id", aggfunc="count", fill_value=0)
             .reindex(["high", "medium", "low"]).to_string())
    print(f"\n  {len(act):,} actionable of {len(df):,} predictions "
          f"({100*len(act)/len(df):.1f}%)")
    print(f"  {(act.review_priority == 'high').sum():,} at high priority — "
          f"the queue a curator would actually work through")

    hr("3. HIGH-PRIORITY RECLASSIFICATION CANDIDATES (top 12 by model score)")
    top = (act[(act.review_priority == "high") & (act.triage_class == T.UPGRADE)]
           .sort_values("am_pathogenicity", ascending=False)
           .head(12))
    # variant_id is shown because one protein change can arise from more than
    # one genomic substitution (EGFR G901R, C797S) — without it those look like
    # duplicate rows rather than distinct variants with distinct assertions.
    print(top[["symbol", "protein_variant", "variant_id", "significance", "stars",
               "am_pathogenicity", "am_class", "plddt", "evidence_tier"]]
          .to_string(index=False))

    hr("4. DISCORDANT — the model contradicts a well-reviewed curator")
    disc = (df[(df.triage_class == T.DISCORDANT) & (df.stars >= 2)]
            .sort_values(["stars", "am_pathogenicity"], ascending=[False, False]))
    print(f"  {len(disc)} cases at >=2 stars\n")
    print(disc[["symbol", "protein_variant", "significance", "stars",
                "am_pathogenicity", "am_class", "evidence_tier"]]
          .head(10).to_string(index=False))

    hr("5. EVIDENCE TIER x TRIAGE CLASS (where structural claims are supportable)")
    print(df[df.triage_class.isin(T.ACTIONABLE)]
          .pivot_table(index="evidence_tier", columns="triage_class",
                       values="variant_id", aggfunc="count", fill_value=0)
          .reindex(list(C.TIER_ORDER)).to_string())

    hr("6. TARGET DOSSIERS")
    facts = queries.target_facts(con)
    burden = (df[df.triage_class == T.UPGRADE].groupby("symbol").size()
                .rename("n_upgrade_candidates"))
    for _, row in facts.iterrows():
        f = row.to_dict()
        f["drug_stages"] = queries.drug_stages(con, f["uniprot_acc"])
        f["n_upgrade_candidates"] = int(burden.get(f["symbol"], 0))
        card = D.build(f)
        d = card["dimensions"]
        print(f"\n  {card['symbol']} — {card['approved_name']}")
        print(f"    ARCHETYPE: {card['archetype']}")
        print(f"    genetic      {d['genetic_evidence']['band']:9s} "
              f"max={d['genetic_evidence']['max_genetic_score']:.3f} "
              f"({d['genetic_evidence']['top_disease']}) "
              f"of {d['genetic_evidence']['n_associated_diseases']:,} diseases")
        print(f"    structure    {d['structural_readiness']['band']:9s} "
              f"{d['structural_readiness']['pct_residues_solved']:.1f}% solved, "
              f"{d['structural_readiness']['n_pdb_entities']} PDB entities, "
              f"global pLDDT {d['structural_readiness']['global_plddt']:.1f}")
        print(f"    binding site {d['binding_site']['band']:9s} "
              f"pocket={d['binding_site']['has_pocket']} "
              f"ligand={d['binding_site']['has_ligand']}")
        print(f"    chemistry    {d['chemical_matter']['band']:9s} "
              f"{d['chemical_matter']['n_drugs']} drugs, "
              f"{d['chemical_matter']['n_trials']:,} trial reports")
        print(f"    variants     {d['variant_burden']['upgrade_candidates']:,} "
              f"reclassification-upgrade candidates")

    con.close()


if __name__ == "__main__":
    main()
