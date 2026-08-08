"""The only module that reads the database.

Core logic receives dataframes and dicts from here and knows nothing about
DuckDB. Keeping every query in one place is what lets the UI layer be swapped
without touching the rules.
"""

from __future__ import annotations

import pandas as pd

import locus_db as db

# One row per (variant, transcript) that has both a prediction and a residue —
# i.e. every missense variant on the canonical transcript. The LEFT JOIN to
# clinical_assertion is deliberate: variants ClinVar has never seen are a triage
# outcome (novel_candidate), not a missing row.
WORKLIST_SQL = """
SELECT
    t.symbol,
    ve.variant_id,
    v.chrom, v.pos, v.ref, v.alt,
    ve.uniprot_acc, ve.transcript_id, ve.protein_variant,
    ve.aa_pos, ve.ref_aa, ve.alt_aa,
    ve.am_pathogenicity, ve.am_class,
    CAST(ca.significance AS VARCHAR) AS significance,
    ca.significance_raw, ca.stars, ca.clinvar_id,
    r.plddt,
    EXISTS (
        SELECT 1 FROM structure_coverage sc
        WHERE sc.uniprot_acc = ve.uniprot_acc
          AND ve.aa_pos BETWEEN sc.ref_beg AND sc.ref_end
    ) AS is_solved
FROM variant_effect ve
JOIN variant v USING (variant_id)
JOIN target  t ON t.uniprot_acc = ve.uniprot_acc
JOIN residue r ON r.uniprot_acc = ve.uniprot_acc AND r.position = ve.aa_pos
LEFT JOIN clinical_assertion ca USING (variant_id)
"""


def worklist(con, symbol: str | None = None) -> pd.DataFrame:
    sql = WORKLIST_SQL + (" WHERE t.symbol = ?" if symbol else "")
    return con.execute(sql, [symbol] if symbol else []).df()


TARGET_FACTS_SQL = """
WITH solved AS (
    SELECT r.uniprot_acc,
           100.0 * COUNT(*) FILTER (WHERE EXISTS (
               SELECT 1 FROM structure_coverage sc
               WHERE sc.uniprot_acc = r.uniprot_acc
                 AND r.position BETWEEN sc.ref_beg AND sc.ref_end)) / COUNT(*) AS pct_solved
    FROM residue r GROUP BY r.uniprot_acc
),
genetic AS (
    SELECT s.uniprot_acc, MAX(s.score) AS max_genetic_score,
           ARG_MAX(d.name, s.score) AS top_genetic_disease
    FROM association_datatype_score s
    JOIN disease d USING (efo_id)
    WHERE s.datatype_id = 'genetic_association'
    GROUP BY s.uniprot_acc
),
prio AS (
    SELECT uniprot_acc,
           MAX(CASE WHEN metric_key = 'hasPocket' THEN metric_value END) AS has_pocket,
           MAX(CASE WHEN metric_key = 'hasLigand' THEN metric_value END) AS has_ligand
    FROM target_prioritisation GROUP BY uniprot_acc
),
trials AS (
    SELECT uniprot_acc, COUNT(DISTINCT report_id) AS n_trials
    FROM target_drug_report GROUP BY uniprot_acc
)
SELECT t.symbol, t.uniprot_acc, t.approved_name, t.global_plddt,
       t.n_assoc_diseases, t.n_drugs, t.n_pdb_entities,
       solved.pct_solved,
       genetic.max_genetic_score, genetic.top_genetic_disease,
       prio.has_pocket, prio.has_ligand,
       COALESCE(trials.n_trials, 0) AS n_trials
FROM target t
LEFT JOIN solved  ON solved.uniprot_acc  = t.uniprot_acc
LEFT JOIN genetic ON genetic.uniprot_acc = t.uniprot_acc
LEFT JOIN prio    ON prio.uniprot_acc    = t.uniprot_acc
LEFT JOIN trials  ON trials.uniprot_acc  = t.uniprot_acc
ORDER BY t.symbol
"""


def target_facts(con) -> pd.DataFrame:
    return con.execute(TARGET_FACTS_SQL).df()


def drug_stages(con, uniprot_acc: str) -> list[str]:
    rows = con.execute(
        "SELECT DISTINCT max_clinical_stage FROM target_drug WHERE uniprot_acc = ?",
        [uniprot_acc],
    ).fetchall()
    return [r[0] for r in rows if r[0]]


def evidence_for(con, table: str, key_col: str, key: str) -> pd.DataFrame:
    """Lineage lookup: which retrieval produced this row."""
    return con.execute(
        f"""SELECT e.source_system, e.resource_url, e.source_version,
                   e.retrieved_at, e.payload_hash
            FROM {table} x JOIN evidence e USING (evidence_id)
            WHERE x.{key_col} = ?""",
        [key],
    ).df()


RESIDUE_PROFILE_SQL = """
SELECT
    r.position, r.wt_aa, r.plddt, r.confidence_band,
    EXISTS (
        SELECT 1 FROM structure_coverage sc
        WHERE sc.uniprot_acc = r.uniprot_acc
          AND r.position BETWEEN sc.ref_beg AND sc.ref_end
    ) AS is_solved,
    MAX(ve.am_pathogenicity)   AS max_am,
    -- AlphaMissense's own class for the most damaging substitution at this
    -- residue. Shipping it means no consumer re-derives the 0.34/0.564 cuts;
    -- a threshold duplicated in the browser is a threshold that can drift.
    ARG_MAX(ve.am_class, ve.am_pathogenicity) AS max_am_class,
    AVG(ve.am_pathogenicity)   AS mean_am,
    COUNT(ve.variant_id)       AS n_predictions,
    COUNT(ca.variant_id)       AS n_clinvar
FROM residue r
LEFT JOIN variant_effect ve
       ON ve.uniprot_acc = r.uniprot_acc AND ve.aa_pos = r.position
LEFT JOIN clinical_assertion ca ON ca.variant_id = ve.variant_id
WHERE r.uniprot_acc = ?
GROUP BY r.uniprot_acc, r.position, r.wt_aa, r.plddt, r.confidence_band
ORDER BY r.position
"""


def residue_profile(con, uniprot_acc: str) -> pd.DataFrame:
    """One row per residue: structural confidence plus variant burden.

    Drives both the 3D colouring and the sequence profile, so the two views
    cannot disagree about a residue.
    """
    return con.execute(RESIDUE_PROFILE_SQL, [uniprot_acc]).df()


def predicted_structure(con, uniprot_acc: str) -> dict | None:
    row = con.execute(
        """SELECT structure_id, file_url, global_plddt, coverage_start, coverage_end
           FROM structure WHERE uniprot_acc = ? AND source = 'alphafold' LIMIT 1""",
        [uniprot_acc],
    ).fetchone()
    if not row:
        return None
    return dict(zip(
        ["structure_id", "file_url", "global_plddt", "coverage_start", "coverage_end"],
        row))


def variant_evidence(con, variant_id: str) -> pd.DataFrame:
    """Full lineage for one variant: every source that contributed to its row."""
    return con.execute("""
        SELECT 'variant' AS object, e.source_system, e.resource_url,
               e.source_version, e.retrieved_at
        FROM variant x JOIN evidence e USING (evidence_id) WHERE x.variant_id = ?
        UNION ALL
        SELECT 'variant_effect', e.source_system, e.resource_url,
               e.source_version, e.retrieved_at
        FROM variant_effect x JOIN evidence e USING (evidence_id) WHERE x.variant_id = ?
        UNION ALL
        SELECT 'clinical_assertion', e.source_system, e.resource_url,
               e.source_version, e.retrieved_at
        FROM clinical_assertion x JOIN evidence e USING (evidence_id) WHERE x.variant_id = ?
    """, [variant_id, variant_id, variant_id]).df()


def targets(con) -> pd.DataFrame:
    return con.execute(
        "SELECT symbol, uniprot_acc, approved_name FROM target ORDER BY symbol").df()


def associations(con, uniprot_acc: str, limit: int = 25) -> pd.DataFrame:
    """Top disease associations with the datatype breakdown that decides whether
    the evidence is genetic or merely literature co-mention."""
    return con.execute("""
        SELECT d.name AS disease, a.efo_id,
               ROUND(a.overall_score, 3) AS overall,
               ROUND(MAX(CASE WHEN s.datatype_id = 'genetic_association'
                              THEN s.score END), 3) AS genetic,
               ROUND(MAX(CASE WHEN s.datatype_id = 'literature'
                              THEN s.score END), 3) AS literature,
               ROUND(MAX(CASE WHEN s.datatype_id = 'known_drug'
                              THEN s.score END), 3) AS known_drug,
               ROUND(MAX(CASE WHEN s.datatype_id = 'somatic_mutation'
                              THEN s.score END), 3) AS somatic,
               string_agg(DISTINCT ta.name, ', ') AS therapeutic_areas
        FROM association a
        JOIN disease d USING (efo_id)
        LEFT JOIN association_datatype_score s
               ON s.uniprot_acc = a.uniprot_acc AND s.efo_id = a.efo_id
        LEFT JOIN disease_therapeutic_area dta ON dta.efo_id = a.efo_id
        LEFT JOIN therapeutic_area ta ON ta.ta_id = dta.ta_id
        WHERE a.uniprot_acc = ?
        GROUP BY d.name, a.efo_id, a.overall_score
        ORDER BY a.overall_score DESC
        LIMIT ?
    """, [uniprot_acc, limit]).df()


def drugs(con, uniprot_acc: str) -> pd.DataFrame:
    return con.execute("""
        SELECT d.name AS drug, d.chembl_id, d.drug_type,
               td.max_clinical_stage AS stage,
               MIN(m.mechanism_of_action) AS mechanism,
               COUNT(DISTINCT tdr.report_id) AS trials
        FROM target_drug td
        JOIN drug d USING (chembl_id)
        LEFT JOIN drug_mechanism m ON m.chembl_id = d.chembl_id
        LEFT JOIN target_drug_report tdr
               ON tdr.uniprot_acc = td.uniprot_acc AND tdr.chembl_id = d.chembl_id
        WHERE td.uniprot_acc = ?
        GROUP BY d.name, d.chembl_id, d.drug_type, td.max_clinical_stage
        ORDER BY trials DESC, d.name
    """, [uniprot_acc]).df()


def experimental_structures(con, uniprot_acc: str, limit: int = 20) -> pd.DataFrame:
    return con.execute("""
        SELECT s.structure_id, s.method, s.resolution,
               s.coverage_start, s.coverage_end, s.title, s.file_url
        FROM structure s
        WHERE s.uniprot_acc = ? AND s.source = 'pdb'
        ORDER BY s.resolution NULLS LAST
        LIMIT ?
    """, [uniprot_acc, limit]).df()


def tractability(con, uniprot_acc: str) -> pd.DataFrame:
    return con.execute("""
        SELECT modality, label, value FROM tractability
        WHERE uniprot_acc = ? AND value ORDER BY modality, label
    """, [uniprot_acc]).df()


def prioritisation(con, uniprot_acc: str) -> pd.DataFrame:
    return con.execute("""
        SELECT metric_key, metric_value FROM target_prioritisation
        WHERE uniprot_acc = ? ORDER BY metric_key
    """, [uniprot_acc]).df()


def connect(read_only: bool = True):
    return db.connect(read_only=read_only)
