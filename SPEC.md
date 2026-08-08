# Locus — Target & Variant Intelligence Console

A coordination surface for two linked questions in early drug discovery and clinical genetics:

1. **Which variants in this gene are misclassified?** — AlphaMissense predictions reconciled against ClinVar's clinical assertions, mapped onto AlphaFold structure, gated by model confidence.
2. **Is this target worth pursuing?** — genetic evidence, structural coverage, chemical matter, and competitive activity assembled into one evidenced dossier.

Both views share one ontology. A variant traverses to its residue, its residue to its structure, its structure to its target, its target to the diseases and trials that justify working on it.

**Non-commercial.** AlphaMissense data is CC BY-NC-SA 4.0. Portfolio and research use only.

---

## 0. Verified data sources

All endpoints below were live-checked on 2026-08-07. No API keys required for any of them.

| Source | Endpoint | Verified | Notes |
|---|---|---|---|
| AlphaFold DB — metadata | `GET https://alphafold.ebi.ac.uk/api/prediction/{uniprot_acc}` | 200 | Returns a **list** (one entry per isoform). Filter `isUniProt: true` for canonical. |
| AlphaFold DB — structure | `pdbUrl` / `cifUrl` field from the API response | 200, 1.19 MB | **Do not hand-build this URL.** Current version is `v6`; hardcoded `v4` 404s. Always read the URL from the API. |
| AlphaFold DB — pLDDT | `plddtDocUrl` (JSON), or B-factor column of the PDB | 200 | Per-residue confidence. |
| **AlphaMissense — per protein** | `https://alphafold.ebi.ac.uk/files/AF-{acc}-F1-aa-substitutions.csv` | 200, 668 KB, 35,397 rows | `protein_variant,am_pathogenicity,am_class` |
| **AlphaMissense — genomic** | `https://alphafold.ebi.ac.uk/files/AF-{acc}-F1-hg38.csv` | 200 | `CHROM,POS,REF,ALT,genome,uniprot_id,transcript_id,protein_variant,am_pathogenicity,am_class` |
| **ClinVar — by region** | `clinvar.vcf.gz` + `.tbi` over HTTP range requests | 200 / 206 | **Preferred.** See "Targeted region fetch" below. ~0.3–0.9 MB per gene. |
| ClinVar — bulk (not used) | `.../vcf_GRCh38/clinvar.vcf.gz` (193 MB) or `variant_summary.txt.gz` (442 MB) | 200 | NCBI serves at **~101 KB/s**, and parallel connections only reach ~177 KB/s aggregate — per-IP throttling. 18–32 min. Avoided. |
| ClinVar — per gene | `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=clinvar&term={gene}[gene]` | 200 | 16,053 BRCA1 records. Rate-limited to 3 req/s without a key. |
| Open Targets | `POST https://api.platform.opentargets.org/api/v4/graphql` | 200 | **The dossier's single source.** Associations + datatype breakdown, tractability, prioritisation, drugs, ChEMBL mechanisms, and clinical trials. |
| gnomAD | `POST https://gnomad.broadinstitute.org/api` | 200 | GraphQL. Population allele frequencies. |
| RCSB PDB — search | `POST https://search.rcsb.org/rcsbsearch/v2/query` | 200 | Accession → polymer entity ids. Returns **204** when there are no hits. |
| RCSB PDB — detail | `POST https://data.rcsb.org/graphql` | 200 | Method, resolution, and UniProt-coordinate aligned regions. |
| ~~ChEMBL~~ | — | — | **Not used.** `Drug.mechanismsOfAction` in Open Targets carries the same mechanism data. |
| ~~ClinicalTrials.gov v2~~ | — | — | **Not used.** Open Targets `clinicalReports` are sourced predominantly from AACT (Aggregate Analysis of ClinicalTrials.gov) — 3,987 of 4,205 stored reports. This removes the dependency that returned HTTP 500 on every endpoint during verification. |

### The architectural find

AlphaFold DB **serves AlphaMissense annotations itself**, per UniProt accession, in both protein and hg38 genomic coordinate form. This removes the 5 GB Zenodo bulk download from the critical path, and the hg38 file arrives with `CHROM,POS,REF,ALT` — which is the exact join key to ClinVar. The two halves of the product were designed to be joined; the join key is given.

### Targeted region fetch — no bulk download anywhere

The ClinVar VCF is bgzip-compressed with a tabix index, and NCBI honours HTTP range requests. `clinvar_regions.py` parses the `.tbi` **linear index** (the binning index is skipped — unnecessary for region lookup), resolves a genomic region to a compressed byte range, fetches only those bytes, and decompresses the BGZF members.

Measured: BRCA1 0.91 MB / 19.5 s, TP53 0.34 MB / 5.1 s, PTEN 0.28 MB / 5.9 s — **~1.5 MB and ~30 s total, versus 193 MB and 18–32 minutes.** The index itself is 610 KB, fetched once.

Consequence for the architecture: **Locus never downloads a bulk file.** ClinVar becomes an on-demand per-target fetch exactly like AlphaFold, the app stays portable, and adding a target costs seconds rather than a re-download. The region is derived from `min(POS)`/`max(POS)` of the target's own AlphaMissense file, so both sides of the join describe the same span by construction.

### The two AlphaMissense files use different class vocabularies

`AF-{acc}-F1-aa-substitutions.csv` emits `LPath` / `LBen` / `Amb`. `AF-{acc}-F1-hg38.csv` emits `likely_pathogenic` / `likely_benign` / `ambiguous`. Same source, same release, different spelling. Any filter written against one silently returns zero rows against the other. Normalise on ingest and raise on unrecognised values — never `.map()` without checking for unmapped leftovers.

### AlphaMissense class thresholds — derived empirically from the BRCA1 file

| Class | `am_pathogenicity` range | BRCA1 count |
|---|---|---|
| `LBen` (likely benign) | ≤ 0.3399 | 24,611 |
| `Amb` (ambiguous) | 0.34 – 0.564 | 5,170 |
| `LPath` (likely pathogenic) | ≥ 0.564 | 5,616 |

Use the `am_class` column directly. Do not re-derive from thresholds — the column is authoritative and the cutoffs may shift between releases.

---

## 1. Ontology

### Object types

**Target**
| Property | Type | Source |
|---|---|---|
| `ensembl_id` (PK) | str | Open Targets |
| `approved_symbol` | str | Open Targets |
| `uniprot_acc` | str | Open Targets → UniProt |
| `biotype` | str | Open Targets |
| `protein_length` | int | AlphaFold |
| `global_plddt` | float | AlphaFold `globalMetricValue` |
| `disordered_fraction` | float | AlphaFold `fractionPlddtVeryLow` |
| `status` | enum | derived |

Links: → `Structure` (1:n), → `Residue` (1:n), → `Association` (1:n), → `Compound` (n:m), → `Trial` (n:m)
Lifecycle: `unscreened → screened → shortlisted | deprioritized`

**Residue** — the join object between sequence, structure, and variant. Without it the app is two tables next to each other; with it, it is one graph.
| Property | Type | Source |
|---|---|---|
| `uniprot_acc, position` (PK) | str, int | AlphaFold |
| `wt_aa` | char | AlphaFold sequence |
| `plddt` | float | AlphaFold B-factor |
| `confidence_band` | enum | derived: `very_low <50 / low 50-70 / confident 70-90 / very_high >90` |
| `mean_am_pathogenicity` | float | derived over 19 substitutions |
| `in_domain` | str? | UniProt features |

Links: → `Target` (n:1), → `Variant` (1:19)

**Variant / VariantEffect / ClinicalAssertion — split three ways (revised after Phase 0b)**

The original design put genomic coordinates, the AlphaMissense prediction, and the ClinVar assertion on one `Variant` row. Phase 0b showed that is wrong: ClinVar aggregates molecular consequences across **all** transcripts (48 for BRCA1, 39 for TP53) while AlphaMissense predicts against MANE Select only, so one genomic variant can be synonymous on MANE and missense elsewhere. Three sources, three refresh cadences, three grains — three tables.

| Table | Grain | Source | Key |
|---|---|---|---|
| `variant` | one genomic substitution | union of AlphaMissense ∪ ClinVar | `variant_id` = `chrom-pos-ref-alt` |
| `variant_effect` | one (variant, transcript) prediction | AlphaMissense | `(variant_id, transcript_id)` |
| `clinical_assertion` | ClinVar's opinion on a variant | ClinVar | `variant_id` |
| `target_variant` | target ↔ variant edge | both | `(uniprot_acc, variant_id)` |

`target_variant` exists because overlapping genes (NBR2/BRCA1) make the relationship many-to-many, and because it keeps ClinVar-only variants — those with no prediction — reachable from their target so the traversal test still passes for them.

`variant_effect` carries a **composite foreign key to `residue (uniprot_acc, position)`**, guaranteeing every prediction lands on a modelled residue. It will fail loudly for proteins over ~2,700 residues where AlphaFold DB splits models into F1/F2 fragments — intended, so the fragment assumption surfaces rather than silently dropping rows.

> **Rule: never predicate on ClinVar's `MC` field.** It aggregates across transcripts and over-counts missense. The authoritative test for "missense on the canonical transcript" is the existence of a `variant_effect` row. `MC` is stored for lineage only.

> **Target PK is `uniprot_acc`, not `ensembl_id`** (deviation from the original ontology). Everything structural — AlphaFold, AlphaMissense, residues — is keyed on UniProt, so making the schema wait on Open Targets would invert the dependency. `ensembl_id` is a nullable attribute filled in Phase 2.

**Structure** — `source` (`alphafold`|`pdb`), `identifier`, `method`, `resolution`, `coverage_start/end`, `model_version`, `file_url`. Links → `Target` (n:1).

**Disease** — `efo_id`, `name`, `therapeutic_area`.

**Association** (link object carrying evidence) — `overall_score`, `genetic_association_score`, `known_drug_score`, `literature_score`, `datasource_breakdown`. Links `Target` ↔ `Disease`.

**Compound** — `chembl_id`, `pref_name`, `max_phase`, `mechanism_of_action`, `molecule_type`.

**Trial** — `nct_id`, `phase`, `overall_status`, `sponsor`, `condition`, `start_date`.

**Evidence** — first-class provenance. `evidence_id`, `claim_ref` (object type + id + field), `source_system`, `source_url`, `retrieved_at`, `source_version`, `payload_hash`. Polymorphic link to any object. Every displayed number resolves to one of these.

**TriageDecision** — append-and-supersede, never updated in place. `decision_id`, `variant_id`, `analyst`, `prior_class`, `new_class`, `rationale`, `decided_at`, `superseded_by`.

### Events

| Event | Emitted by | Touches |
|---|---|---|
| `target.screened` | ingest completion | Target |
| `structure.linked` | AlphaFold/PDB ingest | Structure, Target |
| `variant.scored` | AlphaMissense ingest | Variant |
| `variant.clinvar_matched` | ClinVar join | Variant |
| `variant.flagged` | triage rule engine | Variant |
| `variant.reviewed` | analyst action | Variant, TriageDecision |
| `target.shortlisted` / `target.deprioritized` | analyst action | Target |
| `dossier.exported` | analyst action | Target |

Every event carries a real timestamp and a source system. Ordering is derivable.

### Actions (preconditions → state change → log)

| Action | Precondition | Effect |
|---|---|---|
| `flag_for_review(variant)` | `triage_class ∈ {reclassification_candidate, discordant}` | status → `under_review`; emits `variant.flagged` |
| `record_triage_decision(variant, new_class, rationale)` | variant under review; rationale non-empty | appends `TriageDecision`, supersedes prior; emits `variant.reviewed` |
| `shortlist_target(target, rationale)` | dossier complete (≥4 of 5 sources returned) | status → `shortlisted`; emits `target.shortlisted` |
| `export_dossier(target)` | `status ≠ unscreened` | emits evidence-stamped HTML; logs `dossier.exported` |

---

## 2. The derived logic — the actual product

### Triage matrix

| ClinVar assertion | AlphaMissense | `triage_class` |
|---|---|---|
| VUS / conflicting | `LPath` | **`reclass_upgrade`** ← the headline |
| VUS / conflicting | `LBen` | `reclass_downgrade` |
| VUS / conflicting | `Amb` | `remains_uncertain` |
| Pathogenic / Likely path. | `LPath` | `concordant` |
| Benign / Likely benign | `LBen` | `concordant` |
| Pathogenic / Likely path. | `LBen` | **`discordant`** ← model contradicts clinic |
| Benign / Likely benign | `LPath` | **`discordant`** |
| absent from ClinVar | `LPath` | `novel_candidate` |
| absent from ClinVar | `LBen` / `Amb` | `unasserted` |

### The confidence gate — the caveat as a feature

Every variant carries `residue_plddt` from its residue. Where pLDDT < 50 the local structure is effectively unmodelled, so structural interpretation of that variant is not supportable — the UI must mark it and exclude it from headline counts, with the exclusion visible rather than silent.

This is not defensive decoration. **BRCA1 has a global pLDDT of 41.59 and 80.4% of residues in the very-low band** — it is largely intrinsically disordered. A tool that renders BRCA1 confidently is lying.

**Corrected by Phase 0:** the pre-build assumption was that BRCA1's disorder would gut the structural layer. It does not. Reclassification candidates concentrate in the *folded* domains — 76.7% of BRCA1's sit at pLDDT ≥ 70, rising to 93.8% (PTEN) and 97.9% (TP53). Pathogenic missense variation clusters where function lives (RING, BRCT), and those regions model well even in a protein that is disordered on average. **Protein-average pLDDT is the wrong statistic; per-residue pLDDT at the variant positions is the right one.** The gate is a filter on a real minority (65 BRCA1 candidates below 50), not a wholesale exclusion — and it must be applied per residue, never per protein.

Two-panel money shot: the same structure colored by **AlphaMissense pathogenicity** and by **pLDDT confidence**, side by side. Where they disagree is where the science actually is.

### Three tiers, not two (added in Phase 2)

Experimental structure coverage from RCSB replaces the binary pLDDT gate with an ordered evidence tier:

**experimentally solved** > **confidently predicted** (pLDDT ≥ 70, no structure) > **predicted but unreliable** (pLDDT < 70, no structure)

Phase 2 proved this is a correction, not a refinement. **EGFR L858R — the single most important activating mutation in non-small-cell lung cancer — sits at pLDDT 51.2**, and L861Q at 43.4: AlphaFold models the activation loop poorly because it is genuinely flexible. A pLDDT-only gate would have quietly down-weighted the most clinically consequential variant in the dataset. Both residues are covered by crystal structures, so the experimental tier rescues them.

Coverage is domain-specific, never whole-protein — EGFR 1M14 covers residues 695–1022 (kinase domain), 1IVO covers 25–646 (ectodomain). Stored as spans in `structure_coverage`, so "is this residue solved?" is a query rather than a denormalised flag.

| Target | % residues solved | Upgrade candidates on a solved residue |
|---|---|---|
| TP53 | 100.0% | 427 / 427 |
| PTEN | 100.0% | 529 / 529 |
| EGFR | 87.4% | 573 / 583 |
| BRCA1 | 17.6% | 258 / 352 |

BRCA1 is the informative row: only 17.6% of the protein is experimentally solved, yet **73% of its reclassification candidates fall on solved residues** — pathogenic variation concentrates in the RING and BRCT domains, which are exactly the parts that crystallise.

---

## 3. Architecture

```
connectors/          one module per source, each returns (records, Evidence[])
  alphafold.py         metadata, structure file, pLDDT, AM annotations
  clinvar.py           bulk gz → DuckDB, filtered to demo genes
  opentargets.py       GraphQL: associations, tractability, known drugs
  gnomad.py            GraphQL: allele frequencies
  rcsb.py              experimental structure coverage
  chembl.py            compounds, mechanisms
  ctgov.py             trials — DEGRADABLE, dossier renders without it
        │
        ▼
  ingest.py          per-target orchestration, evidence stamping, HTTP cache
        │
        ▼
  DuckDB (locus.duckdb)   one table per object type + evidence + events + decisions
        │
        ▼
  core/                pure functions, no I/O, fully testable
    triage.py            the matrix above
    confidence.py        pLDDT banding and gating
    dossier.py           target scorecard assembly
        │
        ▼
  export_site.py       runs core/ once, writes columnar JSON + backbone PDBs
        │
        ▼
  site/                static console — no server, no framework, no CDN
```

**Key call: the core is a library, not a Streamlit script.** `locus.core` takes dataframes and returns dataframes with zero I/O and zero UI imports. This is what makes the triage logic unit-testable against hand-checked variants, and what makes the React port a re-skin instead of a rewrite. Every previous prototype that fused logic into the UI had to be rebuilt to change surface.

**Store: DuckDB.** Already on this machine. Reads the 442 MB ClinVar `.gz` directly via `read_csv` without decompressing, joins it to the AlphaMissense hg38 CSV on `(CHROM, POS, REF, ALT)`, and holds the whole demo in one portable file.

**HTTP cache on disk** (`hishel`) so the demo runs with the wifi off. ClinicalTrials.gov threw 500s for several minutes during verification — a live demo cannot depend on it being up.

**3D viewer**: `py3Dmol` + `stmol` for Streamlit; Mol* for the React port. Residues colored from a per-residue value array — same data contract in both, so the viewer swap is isolated.

**No GPU required anywhere.** Nothing is folded, inferred, or trained locally — every model output consumed is precomputed and published. The compute here is joins and rules.

---

## 4. UI surface map

| Surface | Composition |
|---|---|
| **Target list** | Dense linked-record table: symbol, disease associations, global pLDDT, variant count, reclassification-candidate count, status. Every cell traverses. |
| **Target dossier** (entity page) | Header properties → association evidence table → structural coverage strip (AlphaFold span vs PDB spans) → chemical matter → trials → scorecard. Action drawer: shortlist / deprioritize / export. |
| **Variant console** (split view) | Left: filterable variant table (triage class, AM score, ClinVar status, pLDDT band). Right: detail pane — one variant, its residue, its evidence trail, its decision history. |
| **Structure view** | Dual-panel 3D: pathogenicity-colored and pLDDT-colored. Click a residue → filters the variant table to that position. Click a variant → focuses the residue. |
| **Residue profile** | Per-position strip across the sequence: mean AM score, pLDDT, domain annotations, ClinVar density. The one view that makes the disordered-region problem obvious at a glance. |
| **Evidence trail** | Any number on any screen → its source system, URL, retrieval timestamp, and version. |
| **Audit log** | Every triage decision and target status change, with analyst, rationale, and supersession chain. |

Traversal test: variant → residue → structure → target → disease → trial, and back, each in ≤2 clicks.

---

## 5. Build order

**Phase 0 — the join spike. ✅ COMPLETE (2026-08-07). Verdict: GO.**
`phase0_join_spike.py`, results persisted to `locus.duckdb`. Run against BRCA1, TP53, PTEN.

| Result | Finding |
|---|---|
| **Join coverage** | **8,929 / 9,003 = 99.2%** of ClinVar missense SNVs matched an AlphaMissense score (BRCA1 99.0%, TP53 99.1%, PTEN 100%). The join is sound. |
| **Calibration** | Where ClinVar is ≥1★ confident: **94.8% agreement on Pathogenic** (n=698), **86.9% on Benign** (n=612). Contradiction 3.6% / 5.6%. See the circularity caveat below. |
| **Headline** | **6,488** uncertain variants (VUS + conflicting) → **1,308 upgrade candidates** and **4,807 downgrade candidates**. |
| **Confidence gate** | 76.7% / 93.8% / 97.9% of upgrade candidates (BRCA1 / PTEN / TP53) sit at pLDDT ≥ 70. The 3D layer earns its place even on the worst-case gene. |
| **Unmatched** | 74 records (0.8%) — 61 BRCA1, 13 TP53, 0 PTEN. Most likely ClinVar `MC=missense_variant` asserted against a non-MANE transcript where the AlphaMissense MANE-Select transcript has no coding change. Diagnose in Phase 2; not a blocker. |
| **Discordant (≥2★)** | 29 cases where the model contradicts curated clinical calls, including **8 TP53 variants classified Benign by a 3★ expert panel that AlphaMissense calls likely pathogenic**. |

> **Circularity caveat — do not skip this when presenting.** The AlphaMissense class thresholds (0.34 / 0.564) were themselves chosen to achieve ~90% precision *on ClinVar*, calibrated via logistic regression on a balanced ~2,526-variant ClinVar validation set. The 94.8% / 86.9% agreement above is therefore **consistent with published calibration, not independent validation of it.** Locus must never claim to have validated AlphaMissense. The defensible claim is narrower and still useful: *AlphaMissense behaves on these three genes as its published calibration predicts, so its scores can be used to prioritise uncertain variants for expert review.*

> **Input-quality caveat.** 5,193 of the 6,488 uncertain variants carry only 1★ (single submitter). Only 1,229 are ≥2★. Headline counts in the UI should default to the ≥2★ subset, with the 1★ pool available but visibly separated.

**Phase 0b — resolve the 74 unmatched.** Half a day. Confirm the non-MANE-transcript hypothesis by checking ClinVar's `MC` transcript against the AlphaMissense `transcript_id`. Determines whether Locus needs multi-transcript handling or can stand on MANE Select alone.

**Phase 1 — ontology + DuckDB schema. ✅ COMPLETE (2026-08-07). Verdict: PASS.**
`schema.sql`, `locus_db.py`, `connectors/{alphafold,clinvar}.py`, `ingest.py`, `validate_phase1.py` → `locus.duckdb`.

10 tables loaded for BRCA1, TP53, PTEN: 3 targets, 3 structures, 2,659 residues, 30,858 variants, 17,746 variant effects, 22,041 clinical assertions, 30,858 target↔variant edges, 12 evidence rows, 12 events.

Exit criteria, all passing:

| Check | Result |
|---|---|
| Referential + provenance integrity (10 checks) | 0 orphans; every object row resolves to an evidence row |
| **Cross-source validation** — AlphaMissense `ref_aa` vs residue derived from AlphaFold's UniProt sequence | **0 disagreements across 17,746 predictions.** Neither source establishes this alone; agreement confirms the isoform and coordinate mapping. |
| Reproduce Phase 0 from the normalised schema | **exact** — 6,488 uncertain / 1,308 upgrades / 4,807 downgrades, and calibration 94.8% / 86.9% |
| Traversal test | variant → effect → residue → structure → target → assertion → evidence, one query, with retrieval timestamp |

Notes carried forward:
- Ingest is **idempotent per target**: `clear_target()` deletes child-first and only removes a variant once no target still claims it, so reloading one gene cannot delete another's data. Evidence and events are never deleted — they are append-only history, not current state.
- Phase 0's `MC`-based filter happened to produce identical counts here (the delta-diagnostic returned 0), because `MC` **over**-counts missense but never under-counts on this data. The rule still stands: the schema does not depend on that coincidence holding for other genes.
- `CREATE TYPE` has no `IF NOT EXISTS` in DuckDB, and prose in SQL comments contains semicolons — `init_schema` strips comments before splitting on `;` and tolerates "already exists".

**Phase 2 — dossier connectors. ✅ COMPLETE (2026-08-07). Verdict: PASS.**
`connectors/{opentargets,rcsb}.py`, 13 new tables, `validate_phase2.py`. **EGFR (P00533) added** to the gene set — it was already in the planned eight, and without a druggable target the dossier half has nothing to show.

**Four connectors collapsed to two.** Probing before writing DDL showed Open Targets already aggregates the other two layers: `Drug.mechanismsOfAction` carries ChEMBL mechanisms, and `clinicalReports` carries trials (3,987 of 4,205 reports sourced from AACT). ChEMBL and ClinicalTrials.gov were both dropped — the second of which removes the flaky dependency from the critical path.

Open Targets schema facts, all found by introspection after a 400:
- `knownDrugs` **no longer exists** — it is `drugAndClinicalCandidates`, which takes **no arguments**
- `Tractability` has `label`, not `id`
- `max_clinical_stage` uses `APPROVAL`, not `PHASE_4`

| Check | Result |
|---|---|
| Referential + provenance integrity (17 checks) | 0 orphans across the dossier layer |
| Dossier discriminates between targets | BRCA1/PTEN: **0 drugs**, `hasPocket=0`, `hasLigand=0`. EGFR: **82 drugs**, 4,097 trial reports, `SM Approved Drug=True` |
| Traversal target → drug → mechanism → trial | Cetuximab (950 trials), erlotinib (605), gefitinib, osimertinib — all correct for EGFR |
| EGFR oncogenic-variant spot-check | L858R, T790M, G719S/A, L861Q, C797S all present and all `LPath` |

**C797S resolves to two distinct genomic variants** (`7-55181398-T-A` and `7-55181399-G-C`) — the same protein change reached through different codon positions. Direct vindication of keying `variant` genomically rather than by protein change.

**Phase 3 — core logic + tests. ✅ COMPLETE (2026-08-07). Verdict: PASS — 101 tests.**
`core/{triage,confidence,dossier,_util}.py`, `queries.py`, `tests/`, `report.py`.

**Architecture held:** `core/` has zero I/O, zero database and zero UI imports; `queries.py` is the only module that reads DuckDB. The rules are therefore testable without a database — 91 of the 101 tests need no store at all.

**A bug found by looking at the data before writing the rules.** The `OTHER` significance bucket was hiding three defects: 1,114 rows had an *empty* `CLNSIG` (no pathogenicity claim, not "other"), `not_provided` was treated the same way, and `Pathogenic/Likely_pathogenic/Pathogenic,_low_penetrance` — genuinely pathogenic — fell through because the collapse function split on `|` and compared whole strings against a fixed set. Fixed with prefix matching plus a new `NO_ASSERTION` value, distinct from `OTHER`, which triages alongside variants ClinVar has never seen.

**The matrix is total by construction** — materialised over the full cross product of significance × prediction, so no pair can fall through to a default, and a test asserts it. Cases the original spec did not cover: `PATH/BENIGN + Amb` → `model_uninformative`, and `OTHER + *` → `not_triaged` (drug_response and risk_factor are assertions on a different axis, not pathogenicity calls).

**Priority is rule-based, not a weighted score.** A composite number cannot be argued with and hides which evidence carried it. `high` requires an experimentally solved residue **and** ≥2★ ClinVar **and** a model score near a rail; every variant carries the reasons that produced its priority.

| Check | Result |
|---|---|
| Unit tests (matrix totality, boundaries, priority rules, archetypes) | 91 pass, no database |
| Known pathogenic BRCA1 (C61G, T1691K, R1699Q/W, V1736A, M1775R) | all **concordant** |
| Known benign BRCA1 polymorphisms (P871L, E1038G, K1183R, S1613G) | all **concordant** |
| Structural separation | pathogenic all solved (pLDDT 87–98), benign none solved (24–35), `min(path) > max(benign)` |
| EGFR L858R rescued by experimental coverage | pLDDT 51.2, `is_solved`, tier `experimental` |
| Phase 0 regression, scoped to the original 3 genes | **exact** — 6,488 / 1,308 / 4,807 |

**Output:** 7,834 actionable of 25,787 predictions (30.4%), of which **462 at high priority** — a queue a curator could actually work through. 30 discordant cases at ≥2★.

**An honest asymmetry the evidence tiers expose:** upgrade candidates concentrate on solved residues (1,787 experimental vs 86 weak), but *downgrades* concentrate in unsolved, low-confidence regions (2,014 experimental vs 3,827 weak). Benign variation lives in disordered linkers. Downgrade recommendations are therefore structurally far less supportable than upgrades, and the UI must not present them with equal weight.

**Archetypes, all four correct:** EGFR *validated druggable*; TP53 *genetically validated, clinically emerging* (9 drugs, none approved); BRCA1 and PTEN *genetically validated, chemically unexplored*.

**Phase 4 — variant console + 3D. ✅ COMPLETE (2026-08-08). Verdict: PASS, verified in-browser.**
`app.py`, `ui/{theme,structure}.py`, `.streamlit/config.toml`. Run with `uv run streamlit run app.py`.

**Layer discipline holds and is now enforced:** an AST check asserts `core/` imports nothing from `streamlit`, `duckdb`, `httpx`, `py3Dmol`, `ui`, `queries`, `locus_db` or `connectors`. The console computes nothing — every rule comes from `core/`, every fact from `queries.py`.

**Colour, assigned by job and validated with the script rather than by eye:**

| Encoding | Job | Validation |
|---|---|---|
| AlphaMissense class | **diverging** (benign ← ambiguous → pathogenic), blue/red poles + neutral gray | all-pairs on the dark surface: CVD ΔE 19.2, normal-vision 29.0, contrast ≥3:1 |
| Structural evidence tier | **ordinal** (a confidence ordering), one blue hue, monotone lightness | `--ordinal`: monotone L, adjacent ΔL ≥ 0.06, light end 3.23:1 |
| Triage class distribution | **nominal** → every bar takes the same slot-1 hue | colouring nominal bars by value would spend the identity channel on what bar length already shows |

Green/red status tokens were **rejected** for the 3D panels: a structure cannot carry per-residue labels, so red/green confusion would have nothing to fall back on. Two measures on different scales get two charts sharing an x-axis — never a second y-axis.

**Two features that came out of building it:**
- **Deep-linkable state.** `?target=EGFR&find=L858R` restores an exact view. A traversal surface whose state cannot be handed to someone else is only half a surface.
- **Search escapes the filters.** L858R is `not_triaged` (ClinVar calls it `drug_response`), so the default worklist correctly hides it — and a name lookup could not reach it. A search is a lookup, not a filter refinement.

**Five bugs that only rendering could have found:**
1. `.properties()` on a *member* of a layer is silently dropped — the burden chart rendered flat.
2. A `rect` with no `y` encoding has no height; the evidence strip rendered blank three times. Fixed by layering the solved spans *behind* the pLDDT curve, where they inherit a working y-scale and are better placed anyway — coverage is exactly what says how much of that curve to trust.
3. `axis=None` on one layer removes the shared x-axis from the entire layer.
4. A per-residue area over 1,863 points is over-plotting; binned to 180 windows the regional signal appears — burden high at RING (1–100) and BRCT (1650+), low through the disordered middle, mirroring the pLDDT profile exactly.
5. **WebGL context exhaustion.** Two full 1,210-residue viewers rebuilt on every rerun froze the renderer. Fixed by stripping side-chain atoms before embedding (cartoons need only the backbone: EGFR 754 KB → 483 KB, 9,392 → 5,965 atoms) and caching viewer HTML on `(accession, mode, focus)`.

**Verified in-browser:** KPI strip, filters, worklist, detail pane with evidence trail, both structure panels, residue focus, and the sequence profile — for BRCA1 and EGFR. The L858R detail pane reads: ClinVar `OTHER` 3★, raw label `drug_response`, AlphaMissense 0.9968 `LPath`, **residue L858 · pLDDT 51.2, solved: yes, tier `experimentally solved`** — the whole three-tier argument in one panel.

**Phase 5 — target dossier + frontend replacement. ✅ COMPLETE (2026-08-08). Verdict: PASS, verified in-browser.**
`export_site.py`, `site/{index.html,app.js,styles.css}`. Build and serve:

```
uv run python export_site.py
uv run python -m http.server 8080 -d site
```

### Streamlit was replaced, on evidence

Streamlit earned Phase 4 — it proved the analytics and the interactions. It is the wrong host for the finished product, and Phase 4 produced the reasons:

1. **The rerun model fights the core gesture.** Selecting a variant re-executes the whole script and re-mounts both WebGL viewers. The product's most common action ran its slowest path, and it froze the renderer twice during verification.
2. **The data is static.** Nothing changes at runtime, so a Python server earns nothing and a backend would be unjustified complexity.
3. **It cannot be shared.** A portfolio surface should be a link; Streamlit needs a running interpreter.
4. **Layout fought the register** — a vertical stack of blocks behind a lot of corrective CSS, never a dense control surface.

The pure-`core/` architecture was built for exactly this swap: `export_site.py` runs the same rules the 101 tests cover and the browser only renders, so the UI cannot disagree with the tests. `app.py` is kept as a Python-native analytical view but is **superseded** by `site/`.

### Payload design

Columnar with dictionary-encoded categoricals. An array of 12,463 row objects is mostly repeated key names; columns plus small integer codes cut it roughly 5x before the server gzips anything. Whole site **4.3 MB**, largest payload 1.1 MB (BRCA1 variants), 3Dmol vendored so there is **no third-party request at view time** — it runs from a directory, from GitHub Pages, or offline.

### What the dossier shows

Archetype, five banded dimensions with the raw evidence beside each band, disease associations split by datatype (the genetic-vs-literature distinction that decides whether an association is real), drugs with mechanism and trial counts, and experimental structures linked to RCSB.

The contrast is the point: **EGFR** — *validated druggable target*, 82 drugs, 4,097 trial reports, pocket and ligand present. **BRCA1** — *genetically validated, chemically unexplored*, genetic 0.98, **0 drugs, no pocket, no ligand**, and every one of its 20 experimental structures covering only residues 1646–1859. That last fact **is** the 17.6%-solved figure, visible directly.

### Fixed during the build

- **`NaN` is not valid JSON**, and `df.where(df.notna(), None)` does not remove it — assigning `None` into a float column stores `NaN`. The browser rejected an entire dossier file. Fixed by casting to object first, and `allow_nan=False` now makes it a build-time failure instead of a shipped one.
- Reference-line labels printed inside the plot sat on the area fill; the profile now has a left gutter.
- `viewer.resize()` after layout, or the model renders against stale canvas dimensions.
- Table cells are `nowrap` for scannability, which made an empty-state sentence force a horizontal scrollbar.

### Frontend caveats — resolved (2026-08-08)

Six caveats raised after Phase 6 were addressed. Measurements first, then fixes:

| Caveat | Resolution |
|---|---|
| **3D residue-click unverified** | Wiring proven: **403 clickable CA atoms for a 403-residue protein**, each carrying the callback. Handler proven: `selectResidue(130)` selects **R130G** — the commonest PTEN mutation in cancer — and refocuses the structure. Only 3Dmol's own raycast is untested, because synthetic `MouseEvent`s do not drive it. |
| **No rendering tests** | `site/view-model.js` extracted — the browser's counterpart to `core/`. **17 tests** cover decoding, filtering, sorting, colour bucketing, binning, chart geometry and the virtualisation window. Better still, the duplicated logic was *deleted*: per-residue tier and AlphaMissense class now ship from Python, so the pLDDT 70 cut and the 0.34/0.564 cuts exist in exactly one place each. A threshold living in two languages is a threshold that drifts. |
| **500-row cap** | Windowed rendering. Verified at full scale: **End reaches row 12,462 of 12,463** with 25 rows in the DOM, measured row height exactly 26px so the scrollbar does not lie. |
| **Uncompressed payloads** | `reasons` was measured at **33.8% of the payload (374 KB) across 45 distinct values** — dictionary-encoded. Plus `serve.py`, a gzip dev server, because `python -m http.server` sends no `Content-Encoding` and made local testing pessimistic in a way that hid real performance. **1,134,053 → 156,178 bytes over the wire, 7.3×.** |
| **Two WebGL contexts** | Viewers are created lazily on `IntersectionObserver` — verified **0 canvases until the section scrolls into view**. Hit-test targets cut from 2,000 atoms to 403 by making only alpha carbons clickable. |
| **No accessibility audit** | Measured contrast against the surface: `critical` 3.62 and `tier-weak` 3.23, both **below 4.5 as small text**. Fixed at the root — table text now wears ink and a swatch carries the colour, which is what the palette method required anyway. Added roving tabindex with arrow/Home/End/PageUp/PageDown/Enter, `role="grid"` with `aria-rowcount`/`aria-rowindex` (correct for a virtualised grid), `:focus-visible` rings, a status live region, labelled control groups, a skip link, and `sr-only` labels on the action forms. |

Two bugs the work exposed in itself: the rAF scroll throttle latched permanently if `drawWindow` ever threw (the table would silently stop following the scrollbar with no error) — now released in a `finally`; and `serve.py` sent no validator, so a rebuilt payload could go on being served from cache, which looks exactly like a code bug.

**A testing-environment finding worth recording:** the automation tab reports `visibilityState: "hidden"` and `hasFocus: false`, so `requestAnimationFrame` never runs, scroll events are not delivered, and `IntersectionObserver` does not fire. Several "failures" during this work were that, not the code. Keyboard navigation calls `drawWindow()` synchronously and was the route that could verify the virtualiser at all.

`window.__locus` exposes viewers and state for introspection, **gated to localhost** so it is absent from a deployed site.

### Final sweep (2026-08-08)

Rebuilt from scratch and re-verified end to end. `ingest.py` reproduced identical
counts, the Phase 0 baseline held **exact** (6,488 / 1,308 / 4,807), both phase
validators passed, **131 tests** passed, and a new `audit_site.py` re-derived
every displayed figure from the store — **87 checks, all pass**. All 12
target × view combinations loaded with **zero console errors**.

Also added: `README.md`, `.gitignore`, and data-source attribution in the footer
(AlphaFold DB, ClinVar, Open Targets, RCSB PDB, 3Dmol.js) — the demo is built
almost entirely on other people's public work and should say so.

**A real bug the sweep caught.** After a full re-ingest, PTEN came back
`screened` rather than `shortlisted`: `ingest.py` rebuilds target rows and resets
status, but the `event` log is append-only and survives. `import_decisions.py`
judged idempotency from the event alone, saw its own old entry, concluded
"already applied", and **silently dropped a governed decision with no error** —
precisely the failure this system exists to prevent. Idempotency is now judged
from the *effect* (does the target actually hold the intended status?) as well
as the record of it, so a re-import both stays idempotent and self-heals after a
rebuild. The variant branch was already effect-checked and was correct.

Smaller fixes: an unused import; an f-string with no placeholders; and the
"delta explained by" query in `validate_phase1.py` was unscoped while the delta
it explained was scoped to three genes.

### Known limitation

On a HiDPI display the two viewers are 2011×900 canvases, and Chrome DevTools' `Page.captureScreenshot` times out compositing them. **The page itself stays responsive** — verified by evaluating against the live DOM while capture was failing. This affects automated screenshots, not use.

**Phase 6 — governed actions + audit log.** Flag, decide, shortlist, export, with supersession. Note this is the one phase that needs a writable store, so it is also the point at which the static site needs a decision: keep decisions client-side (localStorage, exportable) or reintroduce a small write API.

**Phase 6 — governed actions + audit log. ✅ COMPLETE (2026-08-08). Verdict: PASS — 13 governance tests + verified in-browser.**
`site/decisions.js`, `tests/decisions.test.mjs`, `import_decisions.py`, audit view in `site/`.

**Decisions are held client-side**, which keeps the zero-backend property. The exported log is the durable artifact, and `import_decisions.py` loads it into `triage_decision` — the table declared in Phase 1 and empty until now — so analyst judgements end up in the same store, under the same provenance discipline, as the facts they were made about.

| Rule | How it is enforced |
|---|---|
| Append and supersede, never overwrite | `superseded_by` is the **only** field ever written on an existing record; a test asserts rationale, outcome, analyst and timestamp are immutable after the fact |
| No decision without a rationale | empty or whitespace-only is refused, and the refused decision is not recorded |
| Preconditions state themselves | a blocked action names the actual triage class or the missing evidence sources — a disabled control that does not say what it wants is not a governed action, just a broken one |
| Targets need a real dossier | shortlisting requires ≥4 of 5 evidence sources; BRCA1 with zero drugs still qualifies, which is the point |
| Outcomes are review verdicts, not clinical calls | `endorsed` / `rejected` / `needs_evidence` / `deferred`. A test greps the vocabulary for `pathogenic|benign|vus|likely` and fails if any outcome reads as a clinical assertion. |

**Decisions carry the data build they were made against.** If ClinVar or AlphaMissense is rebuilt, a decision taken against the older evidence is marked *stale build* rather than silently presenting as current — the lineage discipline the system applies to facts, applied to judgements.

**Honesty, stated in the UI and not only here:** browser storage records who *claims* to have decided what. That is attribution, not authentication, and it is not tamper-proof.

Verified end to end: two decisions on PTEN I33N (the second superseding the first) plus a target shortlist taken against an older build → exported → `import_decisions.py` → 2 `triage_decision` rows with supersession intact, `PTEN.status = shortlisted`, and events attributed to the analyst rather than `ingest`.

**A bug this caught in itself:** re-import was idempotent for variants but not for targets, so a second run appended a phantom `target.shortlisted` event — a false attribution in an audit log. Target decisions have no table of their own, so the event log is now what makes them idempotent. The phantom event was removed and the state rebuilt rather than left in place; a corrupted audit record is worse than a targeted correction.

**Test suites (both must pass):**
```
uv run --with pytest pytest -q          # 101 — core rules, schema, integration
node --test tests/decisions.test.mjs    #  13 — governance rules
```

**Phase 7 — React cockpit. ⛔ SUPERSEDED.** The frontend was already replaced in Phase 5, and the replacement is a static site with no framework, no bundler and no dependencies. React would add a build step and a toolchain to a surface that does not need either. All seven phases are closed.

Remaining work is not a phase but a choice: **more targets.** The gene set was eight in the original plan and stands at four. Adding CFTR, KRAS, MLH1 and SCN1A is one edit to `GENES` in `ingest.py` plus a re-export — with the caveat that any protein over ~2,700 residues will trip the deliberate `variant_effect → residue` foreign key, because AlphaFold DB splits those into F1/F2 fragments and this schema assumes one.

Curated demo gene set (mix of well-folded and disordered, well-annotated and sparse): `BRCA1`, `TP53`, `CFTR`, `PTEN`, `KRAS`, `EGFR`, `MLH1`, `SCN1A`.

---

## 6. Risks

| Risk | Mitigation |
|---|---|
| ~~ClinVar↔AlphaMissense join rate is poor~~ | **Retired by Phase 0 — 99.2%** |
| AlphaFold file version drift (`v4`→`v6` already broke) | Never construct file URLs; always read `pdbUrl` from the API |
| ClinicalTrials.gov outages (observed) | Degradable connector; dossier renders and states the gap |
| NCBI throttles to ~101 KB/s | Tabix range fetch; no bulk download in the product at all |
| Class vocabulary differs between the two AlphaMissense files | Normalise on ingest, raise on unrecognised values |
| Overstating the calibration result as validation | Circularity caveat is load-bearing — thresholds were fitted on ClinVar |
| 80% of the uncertain pool is 1★ single-submitter | Default headline counts to ≥2★; separate the 1★ pool visibly |
| Multi-isoform targets (BRCA1 returns 8+ AlphaFold entries) | Filter `isUniProt: true`; pin canonical accession per target |
| AlphaMissense is CC BY-NC-SA | Non-commercial footer; no commercial deployment |
| Over-claiming clinical utility | The tool surfaces reclassification *candidates* for expert review. It never asserts a clinical classification. This framing is load-bearing and belongs in the UI, not just the README. |
