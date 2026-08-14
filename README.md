# Locus — target &amp; variant intelligence console

**[Open the console →](https://jadzoghaib.github.io/genetic-variant-triage/)**
 · **[Read the manual →](https://jadzoghaib.github.io/genetic-variant-triage/manual.html)**

> **Rebuilt from source every Thursday.** The console header shows when the data
> was last refreshed and which ClinVar and Open Targets releases it stands on.
>
> It is still a snapshot, not a live feed — it never calls an upstream API at
> view time, which is why it cannot break when one goes down. **Reloading the
> page does not fetch new data**; it re-downloads the same files the last build
> produced. [How the refresh works →](#refreshing-the-data)


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

### Engineering workflow (implemented tips)

1. **Filter workflow runs first.** Start CI triage with the relevant branch, event
   and status filters so you inspect only the failing run set.
   ```bash
   gh run list --workflow "Build and deploy site" --branch master --event push --status failure
   ```
2. **Pull failed logs before anything else.** Fetch logs for failed jobs first,
   then drill into a specific job only if needed.
   ```bash
   gh run view <run-id> --log-failed
   gh run view <run-id> --job <job-id> --log
   ```
3. **Prefer targeted code discovery tools.** Use file/content search (`glob`,
   `rg`, `view`) for repository exploration instead of broad shell probing.
   ```bash
   rg "symbol_or_rule" .
   ```
4. **Report progress by milestone.** Group updates into clear milestones:
   plan, implementation chunk(s), and validation.
5. **Use one validation gate order.** Run checks in this order before final
   merge: project tests/validators, automated code review, fixes for accepted
   findings, CodeQL scan, then secret scan. If the automated code-review tool is
   unavailable in the current environment, record that clearly and run a manual
   review before merge.

---

## Refreshing the data

`.github/workflows/pages.yml` rebuilds the whole thing **every Thursday at
06:00 UTC**, and on every push that touches code. It runs `ingest.py` against
the live sources, regenerates the site, runs the full test suite plus both phase
validators and `audit_site.py`, and only then deploys. A build that breaks an
invariant leaves the previous deployment live rather than shipping something
wrong. The first end-to-end run took **64 seconds**.

**The payloads are not committed.** Pages deploys an uploaded artifact rather
than a branch, so the workflow generates `site/data/` on each run and throws it
away. Committing ~7 MB of regenerated data weekly would add well over 100 MB a
year and turn the history into data churn; this way the repository stays at
about 0.6 MB and the deployed site is always the output of the pipeline the
tests cover, never a snapshot that has quietly drifted from it.

Nothing is cached on purpose. The connectors key their download cache on
filename with no version in it, so a warm cache would happily serve a superseded
AlphaFold model release forever — exactly the staleness the schedule exists to
prevent.

To refresh by hand, either trigger the workflow (`gh workflow run pages.yml`) or
build locally:

```bash
uv run python ingest.py          # the only step that contacts an upstream source
uv run python export_site.py     # regenerate site/data
uv run python serve.py           # http://localhost:8080
```

One thing that stays true regardless: **reloading the deployed page fetches
nothing new.** Scores, structures and dossiers are baked into `site/data/` at
build time and the browser only ever reads those files. A hard reload helps in
exactly one case — a rebuild has already been published and you suspect a cached
copy.

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

If you have recorded review decisions, a refresh does not silently invalidate
them. Every decision stores the data build it was made against, so after a
rebuild an earlier judgement is flagged **stale build** in the audit log rather
than presented as current — the evidence behind it may have moved, and you are
told rather than left to remember.

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

### Licensing

**The code is MIT** ([LICENSE](LICENSE)) — use it, fork it, lift whatever is
useful.

That is clean because this repository contains only original code. The generated
payloads are not committed (the Pages workflow builds them), so nothing here is
a derivative of a restricted dataset and the pipeline carries no upstream
restriction.

The obligation attaches to the **output**, not the source. Anything
`export_site.py` produces contains values derived from AlphaMissense, which is
**CC BY-NC-SA 4.0**:

- **Attribution** — satisfied by the console footer and the source table above.
- **NonCommercial** — a portfolio or research deployment is fine; a commercial
  product is not.
- **ShareAlike** — a derivative of the *data* must be offered under the same
  terms. This does not reach the code that generated it.

So the deployed site is fine as published, and anyone who clones this and runs
the pipeline inherits the same terms on what they build.

Cite the underlying work if you build on this:
Cheng *et al.*, *Science* (2023) for AlphaMissense; Jumper *et al.*, *Nature*
(2021) and Varadi *et al.*, *NAR* (2024) for AlphaFold and its database.
