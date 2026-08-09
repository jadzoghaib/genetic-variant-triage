"""Ingest — populate the ontology.

Phase 1 (variant layer): AlphaFold DB + AlphaMissense + ClinVar.
Phase 2 (dossier layer): Open Targets + RCSB PDB.

Load order follows the ontology's dependency structure, and every insert
carries the evidence_id of the retrieval it came from.
"""

from __future__ import annotations

import pandas as pd

import locus_db as db
from connectors import alphafold as af
from connectors import clinvar as cv
from connectors import opentargets as ot
from connectors import rcsb

# symbol -> (uniprot accession, ensembl gene id)
#
# The set is chosen to span the axes the tool has to handle, not to be a list of
# famous genes. Every archetype the dossier can emit appears at least once, and
# the structural range runs from a protein that is mostly disordered to one that
# is almost entirely solved.
#
# Adding a gene is this one edit plus a re-export — but anything over ~2,700
# residues will trip the variant_effect -> residue foreign key, because
# AlphaFold DB splits those into F1/F2 fragments and this schema assumes one.
# That constraint was written to fail loudly rather than drop rows quietly.
GENES = {
    "BRCA1": ("P38398", "ENSG00000012048"),   # pLDDT 41.6 — largely disordered
    "TP53":  ("P04637", "ENSG00000141510"),   # pLDDT 75.1 — mixed
    "PTEN":  ("P60484", "ENSG00000171862"),   # pLDDT 83.0 — well folded
    "EGFR":  ("P00533", "ENSG00000146648"),   # druggable kinase, 392 PDB entities
    "CFTR":  ("P13569", "ENSG00000001626"),   # recessive; a transporter with
                                              # approved modulator drugs
    "KRAS":  ("P01116", "ENSG00000133703"),   # 189 aa, pLDDT 91.5 — the classic
                                              # "undruggable" oncogene that only
                                              # recently acquired inhibitors
    "MLH1":  ("P40692", "ENSG00000076242"),   # mismatch repair, Lynch syndrome
    "SCN1A": ("P35498", "ENSG00000144285"),   # 2,009 aa ion channel, and the one
                                              # non-cancer gene in the set
}

ASSEMBLY = "GRCh38"

STRUCTURE_COLS = ["structure_id", "uniprot_acc", "source", "model_version", "method",
                  "resolution", "coverage_start", "coverage_end", "global_plddt",
                  "file_url", "title", "evidence_id"]


def _band(plddt: float) -> str:
    if plddt < 50:
        return "very_low"
    if plddt < 70:
        return "low"
    if plddt < 90:
        return "confident"
    return "very_high"


# ─────────────────────────────────────────────── Phase 1: variant layer ──

def ingest_variants(con, symbol: str, acc: str) -> dict:
    stats: dict[str, int] = {}

    entry, raw, url = af.entry(acc)
    ev_af = db.record_evidence(con, af.SOURCE, url, raw,
                               source_version=f"model_v{entry['latestVersion']}")
    seq = af.sequence(entry)

    con.execute(
        """INSERT INTO target (uniprot_acc, symbol, protein_length, global_plddt,
                               frac_very_low, status, evidence_id)
           VALUES (?, ?, ?, ?, ?, ?, ?)""",
        [acc, symbol, len(seq), entry["globalMetricValue"],
         entry["fractionPlddtVeryLow"], "screened", ev_af],
    )
    con.execute(
        f"INSERT INTO structure ({', '.join(STRUCTURE_COLS)}) "
        f"VALUES ({', '.join('?' * len(STRUCTURE_COLS))})",
        [entry["modelEntityId"], acc, "alphafold", str(entry["latestVersion"]),
         entry.get("toolUsed"), None, entry.get("uniprotStart", 1),
         entry.get("uniprotEnd", len(seq)), entry["globalMetricValue"],
         entry["pdbUrl"], f"AlphaFold predicted structure of {symbol}", ev_af],
    )
    stats["protein_length"] = len(seq)

    plddt_df, raw, url = af.plddt(entry, symbol)
    ev_plddt = db.record_evidence(con, af.SOURCE, url, raw,
                                  source_version=f"model_v{entry['latestVersion']}",
                                  record_count=len(plddt_df))
    if len(plddt_df) != len(seq):
        raise RuntimeError(
            f"{symbol}: pLDDT covers {len(plddt_df)} residues but sequence is {len(seq)}"
        )
    plddt_df = plddt_df.assign(
        uniprot_acc=acc,
        wt_aa=[seq[p - 1] for p in plddt_df["position"]],
        confidence_band=plddt_df["plddt"].map(_band),
        evidence_id=ev_plddt,
    )
    stats["residues"] = db.insert_df(
        con, "residue", plddt_df,
        ["uniprot_acc", "position", "wt_aa", "plddt", "confidence_band", "evidence_id"],
    )

    am, raw, url = af.alphamissense(entry, symbol)
    ev_am = db.record_evidence(con, "alphamissense", url, raw, record_count=len(am))
    am["variant_id"] = (am["chrom"] + "-" + am["pos"].astype(str) + "-"
                        + am["ref"] + "-" + am["alt"])

    am_variants = am.drop_duplicates(subset="variant_id").assign(
        assembly=ASSEMBLY, is_snv=True, discovered_via="alphamissense", evidence_id=ev_am)
    stats["variants_from_am"] = db.insert_df(
        con, "variant", am_variants,
        ["variant_id", "assembly", "chrom", "pos", "ref", "alt",
         "is_snv", "discovered_via", "evidence_id"])

    stats["variant_effects"] = db.insert_df(
        con, "variant_effect", am.assign(uniprot_acc=acc, is_mane_select=True,
                                         evidence_id=ev_am),
        ["variant_id", "transcript_id", "uniprot_acc", "aa_pos", "ref_aa", "alt_aa",
         "protein_variant", "am_pathogenicity", "am_class", "is_mane_select",
         "evidence_id"])

    chrom, beg, end = am["chrom"].iloc[0], int(am["pos"].min()), int(am["pos"].max())
    cvdf, raw, url = cv.fetch_gene(_CV_INDEX, symbol, chrom, beg, end)
    ev_cv = db.record_evidence(con, cv.SOURCE, url, raw, record_count=len(cvdf),
                               source_version=_CV_INDEX.last_modified)
    stats["clinvar_records"] = len(cvdf)

    known = set(am_variants["variant_id"])
    stats["variants_clinvar_only"] = db.insert_df(
        con, "variant",
        cvdf[~cvdf["variant_id"].isin(known)].assign(
            assembly=ASSEMBLY, discovered_via="clinvar", evidence_id=ev_cv),
        ["variant_id", "assembly", "chrom", "pos", "ref", "alt",
         "is_snv", "discovered_via", "evidence_id"])

    stats["clinical_assertions"] = db.insert_df(
        con, "clinical_assertion", cvdf.assign(evidence_id=ev_cv),
        ["variant_id", "clinvar_id", "significance_raw", "significance",
         "review_status", "stars", "molecular_consequence", "evidence_id"])

    edges = pd.concat([
        pd.DataFrame({"variant_id": am_variants["variant_id"],
                      "attribution": "alphamissense_transcript"}),
        pd.DataFrame({"variant_id": cvdf["variant_id"], "attribution": "clinvar_geneinfo"}),
    ], ignore_index=True)
    edges = (edges.groupby("variant_id", as_index=False)["attribution"]
                  .agg(lambda s: "+".join(sorted(set(s)))).assign(uniprot_acc=acc))
    stats["target_variant_edges"] = db.insert_df(
        con, "target_variant", edges, ["uniprot_acc", "variant_id", "attribution"])

    db.emit_event(con, "target.screened", "target", acc,
                  payload={"symbol": symbol, "protein_length": len(seq)})
    db.emit_event(con, "variant.scored", "target", acc,
                  payload={"effects": stats["variant_effects"]})
    return stats


# ─────────────────────────────────────────────── Phase 2: dossier layer ──

def ingest_dossier(con, symbol: str, acc: str, ensembl_id: str) -> dict:
    stats: dict[str, int] = {}

    t, raw, url = ot.fetch(ensembl_id, acc)
    ev_ot = db.record_evidence(con, ot.SOURCE, url, raw, source_version=_OT_RELEASE)
    rec = ot.flatten(t, acc)

    con.execute(
        """UPDATE target SET ensembl_id = ?, approved_name = ?, biotype = ?,
                             n_assoc_diseases = ?, n_drugs = ?
           WHERE uniprot_acc = ?""",
        [t["id"], t.get("approvedName"), t.get("biotype"),
         (t.get("associatedDiseases") or {}).get("count"),
         (t.get("drugAndClinicalCandidates") or {}).get("count"), acc],
    )
    stats["assoc_diseases_total"] = (t.get("associatedDiseases") or {}).get("count", 0)
    stats["drugs_total"] = (t.get("drugAndClinicalCandidates") or {}).get("count", 0)

    ins = db.insert_records
    stats["diseases"] = ins(con, "disease", rec["disease"],
                            ["efo_id", "name", "evidence_id"], ["efo_id"], ev_ot, True)
    ins(con, "therapeutic_area", rec["therapeutic_area"], ["ta_id", "name"], ["ta_id"],
        or_ignore=True)
    ins(con, "disease_therapeutic_area", rec["disease_therapeutic_area"],
        ["efo_id", "ta_id"], ["efo_id", "ta_id"], or_ignore=True)

    stats["associations"] = ins(
        con, "association", rec["association"],
        ["uniprot_acc", "efo_id", "overall_score", "novelty", "evidence_id"],
        ["uniprot_acc", "efo_id"], ev_ot)
    stats["datatype_scores"] = ins(
        con, "association_datatype_score", rec["association_datatype_score"],
        ["uniprot_acc", "efo_id", "datatype_id", "score"],
        ["uniprot_acc", "efo_id", "datatype_id"])

    stats["tractability"] = ins(
        con, "tractability", rec["tractability"],
        ["uniprot_acc", "modality", "label", "value", "evidence_id"],
        ["uniprot_acc", "modality", "label"], ev_ot)
    stats["prioritisation"] = ins(
        con, "target_prioritisation", rec["target_prioritisation"],
        ["uniprot_acc", "metric_key", "metric_value", "evidence_id"],
        ["uniprot_acc", "metric_key"], ev_ot)

    stats["drugs"] = ins(con, "drug", rec["drug"],
                         ["chembl_id", "name", "drug_type", "evidence_id"],
                         ["chembl_id"], ev_ot, True)
    stats["mechanisms"] = ins(
        con, "drug_mechanism", rec["drug_mechanism"],
        ["chembl_id", "mechanism_of_action", "action_type", "target_name", "evidence_id"],
        ["chembl_id", "mechanism_of_action"], ev_ot, True)
    stats["target_drug"] = ins(
        con, "target_drug", rec["target_drug"],
        ["uniprot_acc", "chembl_id", "max_clinical_stage", "evidence_id"],
        ["uniprot_acc", "chembl_id"], ev_ot)
    stats["clinical_reports"] = ins(
        con, "clinical_report", rec["clinical_report"],
        ["report_id", "source", "url", "trial_phase", "trial_status", "evidence_id"],
        ["report_id"], ev_ot, True)
    stats["report_links"] = ins(
        con, "target_drug_report", rec["target_drug_report"],
        ["uniprot_acc", "chembl_id", "report_id"],
        ["uniprot_acc", "chembl_id", "report_id"])

    # ── RCSB experimental structures ────────────────────────────────────
    ids, total = rcsb.search_entities(acc)
    stats["pdb_entities_total"] = total
    if ids:
        structures, coverage, raw = rcsb.fetch_details(ids, acc)
        ev_pdb = db.record_evidence(con, rcsb.SOURCE, rcsb.SEARCH, raw,
                                    record_count=len(structures))
        cov_df = pd.DataFrame(coverage).drop_duplicates(
            subset=["structure_id", "ref_beg", "ref_end"]) if coverage else pd.DataFrame()

        # Span a structure covers, derived from its own aligned regions.
        bounds = (cov_df.groupby("structure_id")
                        .agg(coverage_start=("ref_beg", "min"),
                             coverage_end=("ref_end", "max"))
                        .reset_index()) if len(cov_df) else pd.DataFrame(
                            columns=["structure_id", "coverage_start", "coverage_end"])

        sdf = (pd.DataFrame(structures)
                 .drop_duplicates(subset="structure_id")
                 .merge(bounds, on="structure_id", how="left")
                 .assign(global_plddt=None, evidence_id=ev_pdb))
        sdf["file_url"] = "https://www.rcsb.org/structure/" + \
                          sdf["structure_id"].str.split("_").str[0]
        stats["pdb_structures"] = db.insert_df(con, "structure", sdf, STRUCTURE_COLS)
        stats["coverage_spans"] = db.insert_df(
            con, "structure_coverage", cov_df,
            ["structure_id", "uniprot_acc", "ref_beg", "ref_end"])
        db.emit_event(con, "structure.linked", "target", acc,
                      payload={"source": "pdb", "entities": len(sdf)})
    else:
        stats["pdb_structures"] = 0
        stats["coverage_spans"] = 0

    con.execute("UPDATE target SET n_pdb_entities = ? WHERE uniprot_acc = ?",
                [total, acc])
    db.emit_event(con, "target.dossier_loaded", "target", acc,
                  payload={"diseases": stats["associations"], "drugs": stats["drugs"]})
    return stats


# ────────────────────────────────────────────────────────────────── main ──

_CV_INDEX = None
_OT_RELEASE = None


def main() -> None:
    global _CV_INDEX, _OT_RELEASE
    con = db.connect()
    db.init_schema(con)

    print("loading ClinVar tabix index...")
    _CV_INDEX = cv.load_index()
    _OT_RELEASE = ot.release()
    print(f"  ClinVar archive     {_CV_INDEX.last_modified}")
    print(f"  Open Targets release {_OT_RELEASE}")

    rows = []
    for symbol, (acc, ensembl_id) in GENES.items():
        print(f"\n=== {symbol} ({acc} / {ensembl_id}) ===")
        db.clear_target(con, acc)
        s = {"symbol": symbol}
        s.update(ingest_variants(con, symbol, acc))
        print("  variant layer:", {k: v for k, v in s.items() if k != "symbol"})
        s.update(ingest_dossier(con, symbol, acc, ensembl_id))
        print("  dossier layer:", {k: s[k] for k in
                                   ("assoc_diseases_total", "associations", "drugs_total",
                                    "drugs", "clinical_reports", "pdb_entities_total",
                                    "pdb_structures", "coverage_spans") if k in s})
        rows.append(s)

    con.execute("CHECKPOINT")
    print("\n" + "=" * 100)
    print(pd.DataFrame(rows).to_string(index=False))
    print(f"\nwritten to {db.DB_PATH}")
    con.close()


if __name__ == "__main__":
    main()
