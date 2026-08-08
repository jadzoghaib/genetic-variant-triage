"""Phase 2 exit criteria.

Beyond integrity, this checks the thing Phase 2 was actually for: whether the
dossier layer discriminates between targets, and whether experimental structure
coverage improves the confidence model over pLDDT alone.
"""

from __future__ import annotations

import locus_db as db

con = db.connect()
failures: list[str] = []


def check(label: str, sql: str) -> int:
    n = con.execute(sql).fetchone()[0]
    print(f"  [{'PASS' if n == 0 else 'FAIL'}] {label}: {n:,}")
    if n:
        failures.append(label)
    return n


def show(title: str, sql: str) -> None:
    print(f"\n{title}")
    print(con.execute(sql).df().to_string(index=False))


print("=" * 84)
print("1. REFERENTIAL + PROVENANCE INTEGRITY (dossier layer)")
print("=" * 84)
check("associations with no disease",
      "SELECT COUNT(*) FROM association a LEFT JOIN disease d USING (efo_id) WHERE d.efo_id IS NULL")
check("datatype scores with no association",
      """SELECT COUNT(*) FROM association_datatype_score s LEFT JOIN association a
         USING (uniprot_acc, efo_id) WHERE a.efo_id IS NULL""")
check("target_drug rows with no drug",
      "SELECT COUNT(*) FROM target_drug td LEFT JOIN drug d USING (chembl_id) WHERE d.chembl_id IS NULL")
check("report links with no clinical_report",
      """SELECT COUNT(*) FROM target_drug_report r LEFT JOIN clinical_report c
         USING (report_id) WHERE c.report_id IS NULL""")
check("coverage spans with no structure",
      """SELECT COUNT(*) FROM structure_coverage sc LEFT JOIN structure s
         USING (structure_id) WHERE s.structure_id IS NULL""")
check("coverage spans outside the protein",
      """SELECT COUNT(*) FROM structure_coverage sc JOIN target t USING (uniprot_acc)
         WHERE sc.ref_beg < 1 OR sc.ref_end > t.protein_length""")
check("disease_therapeutic_area with no therapeutic_area",
      """SELECT COUNT(*) FROM disease_therapeutic_area dta LEFT JOIN therapeutic_area ta
         USING (ta_id) WHERE ta.ta_id IS NULL""")
for tbl in ("disease", "association", "tractability", "target_prioritisation",
            "drug", "drug_mechanism", "target_drug", "clinical_report", "structure"):
    check(f"{tbl} rows with unresolvable evidence_id",
          f"""SELECT COUNT(*) FROM {tbl} t LEFT JOIN evidence e USING (evidence_id)
              WHERE e.evidence_id IS NULL""")
check("targets missing dossier enrichment",
      "SELECT COUNT(*) FROM target WHERE ensembl_id IS NULL OR approved_name IS NULL")

print("\n" + "=" * 84)
print("2. DOES THE DOSSIER DISCRIMINATE BETWEEN TARGETS?")
print("=" * 84)
show("Druggability signal (Open Targets prioritisation + drugs + structures):", """
    SELECT t.symbol,
           t.n_assoc_diseases AS diseases,
           t.n_drugs          AS drugs,
           t.n_pdb_entities   AS pdb,
           MAX(CASE WHEN p.metric_key = 'hasPocket'   THEN p.metric_value END) AS has_pocket,
           MAX(CASE WHEN p.metric_key = 'hasLigand'   THEN p.metric_value END) AS has_ligand,
           MAX(CASE WHEN p.metric_key = 'isInMembrane' THEN p.metric_value END) AS membrane,
           MAX(CASE WHEN tr.label = 'Approved Drug' AND tr.modality = 'SM'
                    THEN tr.value END) AS sm_approved_drug
    FROM target t
    LEFT JOIN target_prioritisation p ON p.uniprot_acc = t.uniprot_acc
    LEFT JOIN tractability tr         ON tr.uniprot_acc = t.uniprot_acc
    GROUP BY t.symbol, t.n_assoc_diseases, t.n_drugs, t.n_pdb_entities
    ORDER BY t.n_drugs DESC
""")

show("Top genetic-evidence associations (datatype = genetic_association):", """
    SELECT t.symbol, d.name AS disease,
           ROUND(a.overall_score, 3) AS overall,
           ROUND(s.score, 3) AS genetic
    FROM association a
    JOIN target t USING (uniprot_acc)
    JOIN disease d USING (efo_id)
    JOIN association_datatype_score s
      ON s.uniprot_acc = a.uniprot_acc AND s.efo_id = a.efo_id
     AND s.datatype_id = 'genetic_association'
    QUALIFY ROW_NUMBER() OVER (PARTITION BY t.symbol ORDER BY s.score DESC) <= 2
    ORDER BY t.symbol, genetic DESC
""")

show("Clinical report sources (AACT = ClinicalTrials.gov, rest are regulatory):", """
    SELECT source, COUNT(*) AS reports FROM clinical_report
    GROUP BY source ORDER BY reports DESC LIMIT 6
""")

print("\n" + "=" * 84)
print("3. THE THREE-TIER CONFIDENCE MODEL")
print("   Phase 1 could only ask 'is pLDDT high?'. With experimental coverage")
print("   the question becomes 'is this residue actually solved?'")
print("=" * 84)
show("All residues by evidence tier:", """
    SELECT t.symbol,
           COUNT(*) FILTER (WHERE cov.covered) AS experimental,
           COUNT(*) FILTER (WHERE NOT cov.covered AND r.plddt >= 70) AS predicted_confident,
           COUNT(*) FILTER (WHERE NOT cov.covered AND r.plddt <  70) AS predicted_weak,
           ROUND(100.0 * COUNT(*) FILTER (WHERE cov.covered) / COUNT(*), 1) AS pct_solved
    FROM residue r
    JOIN target t USING (uniprot_acc)
    JOIN LATERAL (SELECT EXISTS (
            SELECT 1 FROM structure_coverage sc
            WHERE sc.uniprot_acc = r.uniprot_acc
              AND r.position BETWEEN sc.ref_beg AND sc.ref_end) AS covered) cov ON TRUE
    GROUP BY t.symbol ORDER BY pct_solved DESC
""")

show("Reclassification-upgrade candidates by evidence tier (the payoff):", """
    SELECT t.symbol,
           COUNT(*) AS upgrade_candidates,
           COUNT(*) FILTER (WHERE cov.covered) AS on_solved_residue,
           COUNT(*) FILTER (WHERE NOT cov.covered AND r.plddt >= 70) AS predicted_confident,
           COUNT(*) FILTER (WHERE NOT cov.covered AND r.plddt <  70) AS predicted_weak
    FROM variant_effect ve
    JOIN clinical_assertion ca USING (variant_id)
    JOIN residue r ON r.uniprot_acc = ve.uniprot_acc AND r.position = ve.aa_pos
    JOIN target  t ON t.uniprot_acc = ve.uniprot_acc
    JOIN LATERAL (SELECT EXISTS (
            SELECT 1 FROM structure_coverage sc
            WHERE sc.uniprot_acc = ve.uniprot_acc
              AND ve.aa_pos BETWEEN sc.ref_beg AND sc.ref_end) AS covered) cov ON TRUE
    WHERE ca.significance IN ('VUS','CONFLICTING') AND ve.am_class = 'LPath'
    GROUP BY t.symbol ORDER BY upgrade_candidates DESC
""")

print("\n" + "=" * 84)
print("4. SPOT-CHECK — do the canonical EGFR oncogenic variants land correctly?")
print("=" * 84)
show("Known EGFR driver / resistance mutations:", """
    SELECT ve.protein_variant,
           ROUND(ve.am_pathogenicity, 3) AS am_score, ve.am_class,
           COALESCE(CAST(ca.significance AS VARCHAR), '(not in ClinVar)') AS clinvar,
           ca.stars, ROUND(r.plddt, 1) AS plddt,
           EXISTS (SELECT 1 FROM structure_coverage sc
                   WHERE sc.uniprot_acc = ve.uniprot_acc
                     AND ve.aa_pos BETWEEN sc.ref_beg AND sc.ref_end) AS solved
    FROM variant_effect ve
    JOIN target t ON t.uniprot_acc = ve.uniprot_acc AND t.symbol = 'EGFR'
    JOIN residue r ON r.uniprot_acc = ve.uniprot_acc AND r.position = ve.aa_pos
    LEFT JOIN clinical_assertion ca USING (variant_id)
    WHERE ve.protein_variant IN ('L858R','T790M','G719S','G719A','L861Q','C797S')
    ORDER BY ve.aa_pos
""")

print("\n" + "=" * 84)
print("5. TRAVERSAL — target to drug to mechanism to trial, one query")
print("=" * 84)
show("Approved drugs acting on a target, with mechanism and trial count:", """
    SELECT t.symbol, d.name AS drug, d.drug_type,
           SUBSTR(MIN(m.mechanism_of_action), 1, 44) AS mechanism,
           COUNT(DISTINCT tdr.report_id) AS trials
    FROM target t
    JOIN target_drug td USING (uniprot_acc)
    JOIN drug d USING (chembl_id)
    LEFT JOIN drug_mechanism m ON m.chembl_id = d.chembl_id
    LEFT JOIN target_drug_report tdr
           ON tdr.uniprot_acc = t.uniprot_acc AND tdr.chembl_id = d.chembl_id
    WHERE td.max_clinical_stage = 'APPROVAL'
    GROUP BY t.symbol, d.name, d.drug_type
    ORDER BY trials DESC LIMIT 6
""")

print("\n" + "=" * 84)
print("6. STORE CONTENTS")
print("=" * 84)
tables = [r[0] for r in con.execute(
    "SELECT table_name FROM information_schema.tables WHERE table_schema='main' ORDER BY 1"
).fetchall()]
for tbl in tables:
    n = con.execute(f"SELECT COUNT(*) FROM {tbl}").fetchone()[0]
    print(f"  {tbl:28s} {n:>9,}")

print("\n" + "=" * 84)
print("VERDICT: " + ("PASS — Phase 3 unblocked" if not failures
                     else f"FAIL — {len(failures)} check(s): {failures}"))
print("=" * 84)
con.close()
raise SystemExit(1 if failures else 0)
