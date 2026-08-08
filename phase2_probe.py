# /// script
# requires-python = ">=3.11"
# dependencies = ["httpx"]
# ///
"""
Phase 2 probe — inspect real payloads before designing the dossier tables.

Hypothesis worth testing: Open Targets `knownDrugs` may already aggregate
ChEMBL mechanisms AND clinical trial ids, which would collapse two of the four
planned connectors into enrichment rather than primary sources.

Probe targets deliberately mix a classically druggable kinase (EGFR) with a
tumour suppressor (BRCA1) — the dossier must behave sensibly when there is no
chemical matter at all, which is the honest answer for BRCA1.
"""

import json

import httpx

OT = "https://api.platform.opentargets.org/api/v4/graphql"
RCSB_SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"
RCSB_GQL = "https://data.rcsb.org/graphql"
CHEMBL = "https://www.ebi.ac.uk/chembl/api/data"
CTGOV = "https://clinicaltrials.gov/api/v2/studies"

TARGETS = {"EGFR": ("ENSG00000146648", "P00533"), "BRCA1": ("ENSG00000012048", "P38398")}


def hr(t):
    print("\n" + "=" * 78 + f"\n{t}\n" + "=" * 78)


def gql(url, query, variables=None):
    r = httpx.post(url, json={"query": query, "variables": variables or {}}, timeout=120)
    r.raise_for_status()
    d = r.json()
    if "errors" in d:
        print("  GraphQL errors:", json.dumps(d["errors"])[:400])
    return d.get("data")


# ───────────────────────────────────────────────────────────── Open Targets ──
hr("1. OPEN TARGETS — target core + tractability")
Q_CORE = """
query T($id: String!) {
  target(ensemblId: $id) {
    id approvedSymbol approvedName biotype
    proteinIds { id source }
    tractability { id modality value }
    associatedDiseases(page: {index: 0, size: 5}) {
      count
      rows { score datatypeScores { id score }
             disease { id name therapeuticAreas { id name } } }
    }
  }
}"""
for sym, (ens, acc) in TARGETS.items():
    d = gql(OT, Q_CORE, {"id": ens})
    t = (d or {}).get("target")
    if not t:
        print(f"  {sym}: no target returned")
        continue
    print(f"\n  {sym} ({t['id']}) {t['approvedName']!r} biotype={t['biotype']}")
    accs = [p["id"] for p in t["proteinIds"] if p["source"] == "uniprot_swissprot"]
    print(f"    swissprot ids: {accs}  (expected {acc} -> {'OK' if acc in accs else 'MISMATCH'})")
    print(f"    tractability buckets: {len(t['tractability'])}")
    for tr in t["tractability"][:6]:
        print(f"      {tr['modality']:14s} {tr['id']:38s} {tr['value']}")
    ad = t["associatedDiseases"]
    print(f"    associated diseases: {ad['count']:,} total, top 5:")
    for row in ad["rows"]:
        dts = {x["id"]: round(x["score"], 3) for x in row["datatypeScores"]}
        ta = row["disease"]["therapeuticAreas"]
        print(f"      {row['score']:.3f}  {row['disease']['name'][:44]:46s} "
              f"TA={ta[0]['name'][:20] if ta else '-':22s} {dts}")

hr("2. OPEN TARGETS — knownDrugs (does it carry mechanism AND trial ids?)")
Q_DRUGS = """
query D($id: String!) {
  target(ensemblId: $id) {
    approvedSymbol
    knownDrugs(size: 5) {
      count uniqueDrugs uniqueDiseases
      rows { drugId prefName drugType mechanismOfAction phase status ctIds
             disease { id name } }
    }
  }
}"""
for sym, (ens, _acc) in TARGETS.items():
    d = gql(OT, Q_DRUGS, {"id": ens})
    kd = ((d or {}).get("target") or {}).get("knownDrugs")
    if not kd:
        print(f"\n  {sym}: knownDrugs returned nothing")
        continue
    print(f"\n  {sym}: {kd['count']:,} rows, {kd['uniqueDrugs']} unique drugs, "
          f"{kd['uniqueDiseases']} diseases")
    for r in kd["rows"]:
        print(f"    {r['drugId']:14s} {(r['prefName'] or '')[:22]:24s} ph{r['phase']} "
              f"{(r['status'] or '-')[:12]:14s} nct={len(r['ctIds'] or [])} "
              f"moa={(r['mechanismOfAction'] or '')[:40]}")

# ────────────────────────────────────────────────────────────────── RCSB ──
hr("3. RCSB — experimental structures for a UniProt accession")
for sym, (_ens, acc) in TARGETS.items():
    q = {
        "query": {"type": "terminal", "service": "text", "parameters": {
            "attribute": "rcsb_polymer_entity_container_identifiers"
                         ".reference_sequence_identifiers.database_accession",
            "operator": "exact_match", "value": acc}},
        "return_type": "polymer_entity",
        "request_options": {"paginate": {"start": 0, "rows": 5},
                            "results_content_type": ["experimental"]},
    }
    r = httpx.post(RCSB_SEARCH, json=q, timeout=120)
    if r.status_code != 200:
        print(f"  {sym}: HTTP {r.status_code} {r.text[:200]}")
        continue
    res = r.json()
    ids = [h["identifier"] for h in res.get("result_set", [])]
    print(f"\n  {sym} ({acc}): {res.get('total_count', 0):,} experimental entities; first: {ids}")
    if ids:
        gq = """query E($ids: [String!]!) {
          polymer_entities(entity_ids: $ids) {
            rcsb_id
            rcsb_polymer_entity_container_identifiers { entry_id }
            entity_poly { rcsb_sample_sequence_length }
            rcsb_polymer_entity_align { aligned_regions { entity_beg_seq_id ref_beg_seq_id length } }
            rcsb_entry_container_identifiers { entry_id }
          }
        }"""
        d = gql(RCSB_GQL, gq, {"ids": ids})
        for e in (d or {}).get("polymer_entities", []) or []:
            al = e.get("rcsb_polymer_entity_align") or []
            regions = al[0]["aligned_regions"][0] if al and al[0].get("aligned_regions") else None
            print(f"      {e['rcsb_id']:10s} len={(e.get('entity_poly') or {}).get('rcsb_sample_sequence_length')} "
                  f"aligned={regions}")

# ──────────────────────────────────────────────────────────────── ChEMBL ──
hr("4. ChEMBL — target and mechanisms by accession")
for sym, (_ens, acc) in TARGETS.items():
    r = httpx.get(f"{CHEMBL}/target.json", params={"target_components__accession": acc,
                                                   "limit": 5}, timeout=120)
    tg = r.json().get("targets", []) if r.status_code == 200 else []
    print(f"\n  {sym} ({acc}): HTTP {r.status_code}, {len(tg)} targets")
    for t in tg[:3]:
        print(f"    {t['target_chembl_id']:14s} {t['target_type']:22s} {t['pref_name'][:40]}")
    if tg:
        cid = tg[0]["target_chembl_id"]
        m = httpx.get(f"{CHEMBL}/mechanism.json",
                      params={"target_chembl_id": cid, "limit": 5}, timeout=120)
        mech = m.json().get("mechanisms", []) if m.status_code == 200 else []
        print(f"    mechanisms for {cid}: {len(mech)}")
        for x in mech[:4]:
            print(f"      {x.get('molecule_chembl_id'):14s} {(x.get('mechanism_of_action') or '')[:50]}")

# ──────────────────────────────────────────────────────── ClinicalTrials ──
hr("5. CLINICALTRIALS.GOV v2")
r = httpx.get(CTGOV, params={"query.intr": "osimertinib", "pageSize": 3, "format": "json",
                             "fields": "NCTId,BriefTitle,Phase,OverallStatus,LeadSponsorName"},
              timeout=120)
print(f"  HTTP {r.status_code}")
if r.status_code == 200:
    d = r.json()
    print(f"  totalCount key present: {'totalCount' in d} | studies: {len(d.get('studies', []))}")
    print("  sample:", json.dumps(d.get("studies", [])[:1])[:600])
