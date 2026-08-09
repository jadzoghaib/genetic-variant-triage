# Locus — target &amp; variant intelligence console

**[Open the console →](https://jadzoghaib.github.io/genetic-variant-triage/)**
 · **[Read the manual →](https://jadzoghaib.github.io/genetic-variant-triage/manual.html)**

> The deployed site is a **snapshot, not a live feed**. It never calls an upstream
> API at view time — that is why it cannot break when one goes down, and why it
> does not update itself. The build date is printed in the console header.
> [How to refresh it →](#refreshing-the-data)


Two linked questions in early drug discovery and clinical genetics, over one
ontology:

1. **Which variants in this gene are misclassified?** AlphaMissense predictions
   reconciled against ClinVar's clinical assertions, mapped onto AlphaFold
   structure, and gated by how much structural evidence actually exists.
2. **Is this target worth pursuing?** Genetic evidence, structural readiness,
   binding-site tractability, chemical matter and clinical activity, assembled
   into one evidenced dossier.

A variant traverses to its residue, its residue to its structure, its structure
to its target, and its target to the diseases and trials that justify working
on it — and every fact on screen resolves to the retrieval that produced it.

> **Non-commercial.** AlphaMissense data is CC BY-NC-SA 4.0. Research and
> portfolio use only.
>
> **Locus surfaces candidates for expert review. It never asserts a clinical
> classification.** That constraint is enforced in the code, not just stated:
> a test fails if any review outcome reads as a clinical call.

---

## Run it

```bash
uv run python ingest.py          # build the ontology from public sources  (~3 min)
uv run python export_site.py     # build the static site
uv run python serve.py           # http://localhost:8080
```

No GPU is required anywhere. Nothing is folded, inferred or trained locally —
every model output consumed is precomputed and published.

### Tests

```bash
uv run --with pytest pytest -q                                   # 101  core rules, schema, integration
node --test tests/decisions.test.mjs tests/view-model.test.mjs   #  30  governance + view-model rules
uv run python validate_phase1.py && uv run python validate_phase2.py
uv run python audit_site.py                                      #  87  exported site vs the store
```

---

## Refreshing the data

Nothing here updates on its own. The site is generated once and served as plain
files, so **new ClinVar assertions or new structures will not appear until you
rebuild and push.** That is a deliberate trade: the deployed console cannot be
broken by an upstream outage, an API change, or a rate limit — but it is only
ever as current as its last build.

The console header prints the build timestamp, so you can always see how old
what you are looking at is.

```bash
uv run python ingest.py          # re-fetch every source          (~6 min for 8 genes)
uv run python export_site.py     # regenerate site/data
uv run python audit_site.py      # confirm the site matches the store
git add -A && git commit -m "Refresh data" && git push
```

The push is the deploy: `.github/workflows/pages.yml` fires on any change under
`site/`, and the new build is live in about twenty seconds.

### What actually changes upstream, and how fast

| Source | Cadence | Worth knowing |
|---|---|---|
| **ClinVar** | weekly | The only source that moves meaningfully week to week — assertions get added, and variants get reclassified out of "uncertain". This is the reason to refresh. |
| **RCSB PDB** | weekly | New depositions can turn a *predicted* residue into a *solved* one, which changes an evidence tier. |
| **Open Targets** | periodic releases | Currently `26.06`. Association scores and drug lists shift between releases. |
| **AlphaFold DB** | infrequent, versioned | Model versions do change — this project lived through `v4` → `v6`. File URLs are always read from the API, never constructed, precisely so that a bump does not 404. |
| **AlphaMissense** | effectively static | A published dataset rather than a moving service. |

### Decisions know when they are stale

If you have recorded review decisions, refreshing does not silently invalidate
them. Every decision stores the data build it was made against, so after a
rebuild any earlier judgement is flagged **stale build** in the audit log rather
than presented as current — the evidence behind it may have moved, and you are
told so instead of having to remember.

### Could this refresh itself?

A scheduled workflow could run the rebuild on a cron. It would work, but weigh
the cost first: each refresh commits several megabytes of regenerated payloads,
so the repository grows with every run, and a silent automated rebuild removes
the chance to look at what changed before it goes live. For a portfolio piece a
manual refresh when you want one is the better trade.

---

## How it is put together

```
connectors/          one module per source; returns records + the evidence to stamp them
  alphafold.py         structure metadata, per-residue pLDDT, AlphaMissense
  clinvar.py           clinical assertions, fetched by genomic region
  opentargets.py       associations, tractability, drugs, mechanisms, trials
  rcsb.py              experimental structures and their UniProt coverage spans
clinvar_regions.py   tabix index + HTTP range reads (see "no bulk downloads")
ingest.py            orchestration -> locus.duckdb
schema.sql           the ontology: 23 tables, provenance and events first-class

core/                PURE rules. No I/O, no database, no UI imports.
  triage.py            the matrix — the product
  confidence.py        structural evidence tiers
  dossier.py           target scorecard and archetypes
queries.py           the ONLY module that reads the database

export_site.py       runs core/ once, writes columnar JSON + backbone PDBs
site/                static console — no framework, no bundler, no build step
  view-model.js        pure browser-side logic (the counterpart of core/)
import_decisions.py  loads an exported decision log back into the ontology
audit_site.py        re-derives every displayed number from the store
serve.py             gzip dev server

app.py               superseded Streamlit console, kept as a Python-native view
report.py            the analytics as a text report, no browser needed
phase0_*.py          the spikes that retired each risk before it was built on;
phase2_probe.py        kept because they record how the design was arrived at
```

The layering is enforced, not merely intended: a check asserts `core/` imports
nothing from `streamlit`, `duckdb`, `httpx`, `py3Dmol`, `ui`, `queries`,
`locus_db` or `connectors`. That is what let the frontend be replaced without
touching a single rule.

---

## Things worth knowing

**No bulk downloads.** NCBI serves the 193 MB ClinVar VCF at ~101 KB/s, and
parallel connections only reach ~177 KB/s — an 18–32 minute wait. The VCF is
bgzip-compressed with a tabix index and NCBI honours byte ranges, so
`clinvar_regions.py` reads the index and fetches only the bytes covering each
gene: **~1.5 MB and seconds per gene.** ClinVar is an on-demand
per-target fetch, exactly like AlphaFold.

**Three tiers of structural evidence, not two.** *experimentally solved >
confidently predicted > predicted but unreliable.* This is a correction, not a
refinement: **EGFR L858R — the most important activating mutation in
non-small-cell lung cancer — sits at pLDDT 51.2**, because AlphaFold models the
activation loop poorly when it is genuinely flexible. A pLDDT-only gate would
have quietly down-weighted the most clinically consequential variant in the
dataset. Crystal structures cover it, so the experimental tier rescues it.

**One place per threshold.** The pLDDT 70 cut and AlphaMissense's 0.34/0.564
cuts live in Python only; the browser receives decided classes as integer
codes. A threshold duplicated across a language boundary is a threshold that
eventually disagrees with itself.

**The calibration result is not independent validation.** AlphaMissense's class
thresholds were themselves fitted to ~90% precision *on ClinVar*. Agreement
with ClinVar is therefore consistent with published calibration, not proof of
it. The defensible claim is narrower and still useful: AlphaMissense behaves on
these genes as its calibration predicts, so its scores can prioritise uncertain
variants for expert review.

**Decisions are attribution, not authentication.** The audit log is append-only
and supersedes rather than overwrites, and every decision records the data
build it was made against — but it lives in browser storage and is not
tamper-proof. Export it; `import_decisions.py` makes it durable.

`SPEC.md` is the living design document and records every phase, measurement
and reversal in full.

---

## Data sources

All public, none requiring an API key.

| Source | Used for | Terms |
|---|---|---|
| [AlphaFold Protein Structure Database](https://alphafold.ebi.ac.uk) (EMBL-EBI &amp; Google DeepMind) | Predicted structures, per-residue pLDDT | CC BY 4.0 |
| [AlphaMissense](https://alphafold.ebi.ac.uk) (Google DeepMind) | Missense pathogenicity predictions | **CC BY-NC-SA 4.0 — non-commercial** |
| [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/) (NCBI) | Clinical significance assertions | Public domain |
| [Open Targets Platform](https://platform.opentargets.org) | Associations, tractability, drugs, mechanisms, trials | CC0 1.0 |
| [RCSB PDB](https://www.rcsb.org) | Experimental structures and coverage | CC0 1.0 |
| [Ensembl VEP](https://rest.ensembl.org) | Transcript consequences (diagnostics only) | Apache 2.0 |
| [3Dmol.js](https://3dmol.csb.pitt.edu) | Structure viewer (vendored) | BSD-3-Clause |

### Licensing note before publishing

The code here is yours to license as you like. **`site/data/` is not** — it
contains values derived from AlphaMissense, which is **CC BY-NC-SA 4.0**. That
carries three obligations that follow the data wherever it goes:

- **Attribution** — satisfied by the footer and the table above.
- **NonCommercial** — a portfolio or research deployment is fine; a commercial
  product is not.
- **ShareAlike** — a derivative of NC-SA data must be offered under the same
  terms. That applies to the exported payloads, not to the pipeline code.

So a public repository and a GitHub Pages deployment are both fine, provided the
non-commercial framing stays visible. If you ever want this under a permissive
licence with no strings, drop `site/data/` from the repository and have viewers
build it themselves with `ingest.py` + `export_site.py` — the pipeline carries
no such restriction.

Cite the underlying work if you build on this:
Cheng *et al.*, *Science* (2023) for AlphaMissense; Jumper *et al.*, *Nature*
(2021) and Varadi *et al.*, *NAR* (2024) for AlphaFold and its database.
