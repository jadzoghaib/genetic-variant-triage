-- Locus ontology — Phase 1.
--
-- Modelling rules this encodes (see SPEC.md §1):
--   * Object-centric, not table-centric. A genomic variant, a transcript-level
--     prediction, and a clinical assertion are three different things with
--     three different sources; they are never flattened into one row.
--   * Every object carries provenance. `evidence` holds one row per retrieval;
--     objects point at the retrieval they came from.
--   * `event` and `triage_decision` are append-only. Analyst decisions supersede
--     rather than overwrite.
--
-- Phase 0b evidence for the variant / variant_effect split: ClinVar reports
-- molecular consequences aggregated across all transcripts (48 for BRCA1, 39
-- for TP53), while AlphaMissense predicts against MANE Select only. A variant
-- that is synonymous on MANE can be missense on another transcript. Those are
-- different facts about different objects.
--
-- Phase 2 will add disease, association, compound and trial for the target
-- dossier. They are deliberately NOT declared here: their shape should be cut
-- against real Open Targets / ChEMBL payloads, not guessed.

CREATE TYPE am_class              AS ENUM ('LPath', 'Amb', 'LBen');
-- NO_ASSERTION is distinct from OTHER: the ClinVar record exists but makes no
-- pathogenicity claim (empty or not_provided CLNSIG). OTHER means a real
-- assertion on a different axis — drug_response, risk_factor, protective.
CREATE TYPE clinical_significance AS ENUM
    ('PATH', 'BENIGN', 'VUS', 'CONFLICTING', 'OTHER', 'NO_ASSERTION');
CREATE TYPE confidence_band       AS ENUM ('very_low', 'low', 'confident', 'very_high');
CREATE TYPE target_status         AS ENUM ('unscreened', 'screened', 'shortlisted', 'deprioritized');
CREATE TYPE structure_source      AS ENUM ('alphafold', 'pdb');

CREATE SEQUENCE IF NOT EXISTS event_seq START 1;


-- ─────────────────────────────────────────────────────────── provenance ──

-- One row per retrieval, not per field. Every variant_effect row produced by a
-- single AlphaMissense fetch genuinely shares that fetch's lineage, so pointing
-- them all at one evidence row is both accurate and compact.
CREATE TABLE IF NOT EXISTS evidence (
    evidence_id    VARCHAR PRIMARY KEY,
    source_system  VARCHAR   NOT NULL,
    resource_url   VARCHAR   NOT NULL,
    source_version VARCHAR,
    retrieved_at   TIMESTAMP NOT NULL,
    payload_hash   VARCHAR   NOT NULL,
    record_count   INTEGER
);

CREATE TABLE IF NOT EXISTS event (
    event_id    BIGINT    PRIMARY KEY,
    event_type  VARCHAR   NOT NULL,
    object_type VARCHAR   NOT NULL,
    object_key  VARCHAR   NOT NULL,
    occurred_at TIMESTAMP NOT NULL,
    actor       VARCHAR   NOT NULL,
    payload     VARCHAR
);


-- ────────────────────────────────────────────────────────────── objects ──

CREATE TABLE IF NOT EXISTS target (
    uniprot_acc      VARCHAR PRIMARY KEY,
    symbol           VARCHAR NOT NULL UNIQUE,
    ensembl_id       VARCHAR,
    approved_name    VARCHAR,
    biotype          VARCHAR,
    protein_length   INTEGER NOT NULL,
    global_plddt     DOUBLE,
    frac_very_low    DOUBLE,
    -- Totals reported by the source, which exceed what is stored: the dossier
    -- shows the top N and must be able to say "of how many".
    n_assoc_diseases INTEGER,
    n_drugs          INTEGER,
    n_pdb_entities   INTEGER,
    status           target_status NOT NULL,
    evidence_id      VARCHAR NOT NULL REFERENCES evidence (evidence_id)
);

CREATE TABLE IF NOT EXISTS structure (
    structure_id   VARCHAR PRIMARY KEY,        -- 'AF-P38398-F1', or a PDB id
    uniprot_acc    VARCHAR NOT NULL REFERENCES target (uniprot_acc),
    source         structure_source NOT NULL,
    model_version  VARCHAR,
    method         VARCHAR,
    resolution     DOUBLE,                     -- NULL for predicted models
    coverage_start INTEGER,
    coverage_end   INTEGER,
    global_plddt   DOUBLE,                     -- predicted models only
    file_url       VARCHAR,
    title          VARCHAR,
    evidence_id    VARCHAR NOT NULL REFERENCES evidence (evidence_id)
);

-- The join object between sequence, structure and variant. Without it the two
-- halves of the product are adjacent tables rather than one graph.
CREATE TABLE IF NOT EXISTS residue (
    uniprot_acc     VARCHAR NOT NULL REFERENCES target (uniprot_acc),
    position        INTEGER NOT NULL,
    wt_aa           VARCHAR NOT NULL,
    plddt           DOUBLE,
    confidence_band confidence_band,
    evidence_id     VARCHAR NOT NULL REFERENCES evidence (evidence_id),
    PRIMARY KEY (uniprot_acc, position)
);

-- Genomic level. Assembly-scoped, transcript-agnostic, gene-agnostic. Holds the
-- union of everything AlphaMissense predicts on and everything ClinVar asserts
-- about — including variants with no prediction and predictions with no
-- assertion.
CREATE TABLE IF NOT EXISTS variant (
    variant_id     VARCHAR PRIMARY KEY,        -- '17-43045682-T-C'
    assembly       VARCHAR NOT NULL,
    chrom          VARCHAR NOT NULL,
    pos            INTEGER NOT NULL,
    ref            VARCHAR NOT NULL,
    alt            VARCHAR NOT NULL,
    is_snv         BOOLEAN NOT NULL,
    discovered_via VARCHAR NOT NULL,           -- which source first introduced it
    evidence_id    VARCHAR NOT NULL REFERENCES evidence (evidence_id)
);

-- Transcript level. One row per (variant, transcript). AlphaMissense publishes
-- MANE Select only, so is_mane_select is TRUE throughout today — but the grain
-- is correct, so adding VEP consequences later needs no migration.
--
-- The composite FK to residue is deliberate: it guarantees every predicted
-- effect lands on a modelled residue. It will fail loudly for proteins over
-- 2,700 residues, where AlphaFold DB splits models into F1/F2 fragments and
-- this schema's single-fragment assumption breaks. That is the intended
-- behaviour — surface it rather than silently drop rows.
CREATE TABLE IF NOT EXISTS variant_effect (
    variant_id       VARCHAR NOT NULL REFERENCES variant (variant_id),
    transcript_id    VARCHAR NOT NULL,
    uniprot_acc      VARCHAR NOT NULL,
    aa_pos           INTEGER NOT NULL,
    ref_aa           VARCHAR NOT NULL,
    alt_aa           VARCHAR NOT NULL,
    protein_variant  VARCHAR NOT NULL,
    am_pathogenicity DOUBLE  NOT NULL,
    am_class         am_class NOT NULL,
    is_mane_select   BOOLEAN NOT NULL,
    evidence_id      VARCHAR NOT NULL REFERENCES evidence (evidence_id),
    PRIMARY KEY (variant_id, transcript_id),
    FOREIGN KEY (uniprot_acc, aa_pos) REFERENCES residue (uniprot_acc, position)
);

-- ClinVar's opinion about a genomic variant. Separate from the variant itself
-- because it is a different source with a different refresh cadence, and
-- because most variants have no assertion at all.
--
-- molecular_consequence is stored for lineage ONLY. It aggregates across every
-- transcript, so predicating on it over-counts missense (Phase 0b). The
-- authoritative test for "missense on the canonical transcript" is the
-- existence of a variant_effect row.
CREATE TABLE IF NOT EXISTS clinical_assertion (
    variant_id            VARCHAR PRIMARY KEY REFERENCES variant (variant_id),
    clinvar_id            VARCHAR,
    significance_raw      VARCHAR,
    significance          clinical_significance NOT NULL,
    review_status         VARCHAR,
    stars                 INTEGER NOT NULL,
    molecular_consequence VARCHAR,
    evidence_id           VARCHAR NOT NULL REFERENCES evidence (evidence_id)
);

-- Target ↔ variant is genuinely many-to-many: overlapping genes (NBR2 overlaps
-- BRCA1) mean one variant can be attributed to several targets. This edge is
-- what keeps ClinVar-only variants — those with no prediction — reachable from
-- their target, so the traversal test still passes for them.
CREATE TABLE IF NOT EXISTS target_variant (
    uniprot_acc VARCHAR NOT NULL REFERENCES target (uniprot_acc),
    variant_id  VARCHAR NOT NULL REFERENCES variant (variant_id),
    attribution VARCHAR NOT NULL,   -- 'alphamissense_transcript' | 'clinvar_geneinfo'
    PRIMARY KEY (uniprot_acc, variant_id)
);


-- ──────────────────────────────────────── dossier layer (Phase 2) ──
--
-- Cut against real payloads, not guessed. Probing collapsed the four planned
-- connectors to two: Open Targets already aggregates ChEMBL mechanisms AND
-- clinical trials (its clinicalReports are dominated by AACT, the Aggregate
-- Analysis of ClinicalTrials.gov), so neither needs its own connector and the
-- flaky ClinicalTrials.gov dependency leaves the critical path entirely.

CREATE TABLE IF NOT EXISTS disease (
    efo_id      VARCHAR PRIMARY KEY,
    name        VARCHAR NOT NULL,
    evidence_id VARCHAR NOT NULL REFERENCES evidence (evidence_id)
);

CREATE TABLE IF NOT EXISTS therapeutic_area (
    ta_id VARCHAR PRIMARY KEY,
    name  VARCHAR NOT NULL
);

CREATE TABLE IF NOT EXISTS disease_therapeutic_area (
    efo_id VARCHAR NOT NULL REFERENCES disease (efo_id),
    ta_id  VARCHAR NOT NULL REFERENCES therapeutic_area (ta_id),
    PRIMARY KEY (efo_id, ta_id)
);

-- Link object: the association IS the evidence-bearing thing, not a join row.
CREATE TABLE IF NOT EXISTS association (
    uniprot_acc   VARCHAR NOT NULL REFERENCES target (uniprot_acc),
    efo_id        VARCHAR NOT NULL REFERENCES disease (efo_id),
    overall_score DOUBLE  NOT NULL,
    novelty       DOUBLE,
    evidence_id   VARCHAR NOT NULL REFERENCES evidence (evidence_id),
    PRIMARY KEY (uniprot_acc, efo_id)
);

-- The breakdown that decides whether an association is real genetic evidence or
-- just literature co-mention. Collapsing it into the overall score would throw
-- away the only thing that makes the dossier trustworthy.
CREATE TABLE IF NOT EXISTS association_datatype_score (
    uniprot_acc VARCHAR NOT NULL,
    efo_id      VARCHAR NOT NULL,
    datatype_id VARCHAR NOT NULL,
    score       DOUBLE  NOT NULL,
    PRIMARY KEY (uniprot_acc, efo_id, datatype_id),
    FOREIGN KEY (uniprot_acc, efo_id) REFERENCES association (uniprot_acc, efo_id)
);

CREATE TABLE IF NOT EXISTS tractability (
    uniprot_acc VARCHAR NOT NULL REFERENCES target (uniprot_acc),
    modality    VARCHAR NOT NULL,
    label       VARCHAR NOT NULL,
    value       BOOLEAN NOT NULL,
    evidence_id VARCHAR NOT NULL REFERENCES evidence (evidence_id),
    PRIMARY KEY (uniprot_acc, modality, label)
);

CREATE TABLE IF NOT EXISTS target_prioritisation (
    uniprot_acc  VARCHAR NOT NULL REFERENCES target (uniprot_acc),
    metric_key   VARCHAR NOT NULL,
    metric_value DOUBLE,
    evidence_id  VARCHAR NOT NULL REFERENCES evidence (evidence_id),
    PRIMARY KEY (uniprot_acc, metric_key)
);

CREATE TABLE IF NOT EXISTS drug (
    chembl_id   VARCHAR PRIMARY KEY,
    name        VARCHAR,
    drug_type   VARCHAR,
    evidence_id VARCHAR NOT NULL REFERENCES evidence (evidence_id)
);

CREATE TABLE IF NOT EXISTS drug_mechanism (
    chembl_id           VARCHAR NOT NULL REFERENCES drug (chembl_id),
    mechanism_of_action VARCHAR NOT NULL,
    action_type         VARCHAR,
    target_name         VARCHAR,
    evidence_id         VARCHAR NOT NULL REFERENCES evidence (evidence_id),
    PRIMARY KEY (chembl_id, mechanism_of_action)
);

CREATE TABLE IF NOT EXISTS target_drug (
    uniprot_acc        VARCHAR NOT NULL REFERENCES target (uniprot_acc),
    chembl_id          VARCHAR NOT NULL REFERENCES drug (chembl_id),
    max_clinical_stage VARCHAR,
    evidence_id        VARCHAR NOT NULL REFERENCES evidence (evidence_id),
    PRIMARY KEY (uniprot_acc, chembl_id)
);

-- Trials and regulatory records. `source` distinguishes AACT (ClinicalTrials.gov)
-- from FDA / EMA / ATC / DailyMed labels, which are not trials at all.
CREATE TABLE IF NOT EXISTS clinical_report (
    report_id    VARCHAR PRIMARY KEY,
    source       VARCHAR,
    url          VARCHAR,
    trial_phase  VARCHAR,
    trial_status VARCHAR,
    evidence_id  VARCHAR NOT NULL REFERENCES evidence (evidence_id)
);

CREATE TABLE IF NOT EXISTS target_drug_report (
    uniprot_acc VARCHAR NOT NULL,
    chembl_id   VARCHAR NOT NULL,
    report_id   VARCHAR NOT NULL REFERENCES clinical_report (report_id),
    PRIMARY KEY (uniprot_acc, chembl_id, report_id),
    FOREIGN KEY (uniprot_acc, chembl_id) REFERENCES target_drug (uniprot_acc, chembl_id)
);

-- UniProt-coordinate spans covered by an experimental structure. Kept as spans
-- rather than a per-residue flag so coverage is a query, not a denormalisation:
-- a residue is experimentally covered iff its position falls in some span.
-- This upgrades the confidence model from pLDDT alone to three tiers —
-- experimentally solved > confidently predicted > predicted but unreliable.
CREATE TABLE IF NOT EXISTS structure_coverage (
    structure_id VARCHAR NOT NULL REFERENCES structure (structure_id),
    uniprot_acc  VARCHAR NOT NULL REFERENCES target (uniprot_acc),
    ref_beg      INTEGER NOT NULL,
    ref_end      INTEGER NOT NULL,
    PRIMARY KEY (structure_id, ref_beg, ref_end)
);


-- ────────────────────────────────────────────────── governed actions ──

-- Append and supersede, never update in place. Current state is the set of
-- rows where superseded_by IS NULL.
CREATE TABLE IF NOT EXISTS triage_decision (
    decision_id   VARCHAR PRIMARY KEY,
    variant_id    VARCHAR   NOT NULL REFERENCES variant (variant_id),
    analyst       VARCHAR   NOT NULL,
    prior_class   VARCHAR,
    new_class     VARCHAR   NOT NULL,
    rationale     VARCHAR   NOT NULL,
    decided_at    TIMESTAMP NOT NULL,
    superseded_by VARCHAR
);
