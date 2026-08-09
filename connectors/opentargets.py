"""Open Targets Platform connector — the dossier's primary source.

Probing (phase2_probe.py) established that Open Targets already aggregates the
layers originally planned as separate connectors:

  * `Drug.mechanismsOfAction` carries the ChEMBL mechanism data
  * `clinicalReports` carries trial ids, phase, status and URLs, sourced mostly
    from AACT (Aggregate Analysis of ClinicalTrials.gov)

so neither ChEMBL nor ClinicalTrials.gov needs its own connector.

Schema notes learned the hard way, all verified by introspection:
  * `knownDrugs` no longer exists — it is `drugAndClinicalCandidates`
  * `drugAndClinicalCandidates` takes NO arguments (returns everything)
  * `Tractability` has `label`, not `id`
  * `associatedDiseases` takes `page: {index, size}`
"""

from __future__ import annotations

import json

import httpx

API = "https://api.platform.opentargets.org/api/v4/graphql"
SOURCE = "open_targets"

# The dossier shows a ranked shortlist; the true total is stored on the target
# so the UI can say "50 of 6,459".
TOP_DISEASES = 50

QUERY = """
query Dossier($id: String!, $size: Int!) {
  target(ensemblId: $id) {
    id approvedSymbol approvedName biotype
    proteinIds { id source }
    tractability { label modality value }
    prioritisation { items { key value } }
    associatedDiseases(page: {index: 0, size: $size}) {
      count
      rows {
        score novelty
        datatypeScores { id score }
        disease { id name therapeuticAreas { id name } }
      }
    }
    drugAndClinicalCandidates {
      count
      rows {
        id maxClinicalStage
        drug {
          id name drugType
          mechanismsOfAction { rows { mechanismOfAction targetName actionType } }
        }
        clinicalReports { id url source trialPhase trialOverallStatus }
      }
    }
  }
}
"""


def release() -> str | None:
    """The Open Targets data release this build drew on, e.g. "26.06".

    Association scores and drug lists move between releases, so "which release"
    is a more useful statement of currency than the moment we happened to ask.
    """
    try:
        r = httpx.post(API, json={"query": "{ meta { dataVersion { year month } } }"},
                       timeout=60)
        v = r.json()["data"]["meta"]["dataVersion"]
        return f"{v['year']}.{v['month']}"
    except Exception:
        return None          # a missing version label must not fail an ingest


def fetch(ensembl_id: str, uniprot_acc: str) -> tuple[dict, bytes, str]:
    """Fetch a target dossier and verify it is the target we asked for.

    The Ensembl id is an external assumption; checking that the returned
    SwissProt accessions contain ours turns a silent wrong-gene bug into a
    loud failure.
    """
    r = httpx.post(API, json={"query": QUERY,
                              "variables": {"id": ensembl_id, "size": TOP_DISEASES}},
                   timeout=300)
    r.raise_for_status()
    payload = r.json()
    if "errors" in payload:
        raise RuntimeError(f"Open Targets errors: {json.dumps(payload['errors'])[:400]}")

    t = payload["data"]["target"]
    if t is None:
        raise RuntimeError(f"Open Targets returned no target for {ensembl_id}")

    swissprot = {p["id"] for p in t["proteinIds"] if p["source"] == "uniprot_swissprot"}
    if uniprot_acc not in swissprot:
        raise RuntimeError(
            f"{ensembl_id} maps to SwissProt {sorted(swissprot)}, expected {uniprot_acc}"
        )
    return t, r.content, f"{API}#{ensembl_id}"


def _as_float(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def flatten(t: dict, uniprot_acc: str) -> dict[str, list[dict]]:
    """Turn the nested response into flat records, one list per table."""
    out: dict[str, list[dict]] = {
        "disease": [], "therapeutic_area": [], "disease_therapeutic_area": [],
        "association": [], "association_datatype_score": [],
        "tractability": [], "target_prioritisation": [],
        "drug": [], "drug_mechanism": [], "target_drug": [],
        "clinical_report": [], "target_drug_report": [],
    }

    for tr in t.get("tractability") or []:
        out["tractability"].append({
            "uniprot_acc": uniprot_acc, "modality": tr["modality"],
            "label": tr["label"], "value": bool(tr["value"]),
        })

    prio = (t.get("prioritisation") or {}).get("items") or []
    for item in prio:
        out["target_prioritisation"].append({
            "uniprot_acc": uniprot_acc, "metric_key": item["key"],
            "metric_value": _as_float(item["value"]),
        })

    for row in (t.get("associatedDiseases") or {}).get("rows") or []:
        d = row["disease"]
        out["disease"].append({"efo_id": d["id"], "name": d["name"]})
        for ta in d.get("therapeuticAreas") or []:
            out["therapeutic_area"].append({"ta_id": ta["id"], "name": ta["name"]})
            out["disease_therapeutic_area"].append({"efo_id": d["id"], "ta_id": ta["id"]})
        out["association"].append({
            "uniprot_acc": uniprot_acc, "efo_id": d["id"],
            "overall_score": row["score"], "novelty": row.get("novelty"),
        })
        for s in row.get("datatypeScores") or []:
            out["association_datatype_score"].append({
                "uniprot_acc": uniprot_acc, "efo_id": d["id"],
                "datatype_id": s["id"], "score": s["score"],
            })

    for row in (t.get("drugAndClinicalCandidates") or {}).get("rows") or []:
        drug = row["drug"]
        out["drug"].append({"chembl_id": drug["id"], "name": drug.get("name"),
                            "drug_type": drug.get("drugType")})
        out["target_drug"].append({
            "uniprot_acc": uniprot_acc, "chembl_id": drug["id"],
            "max_clinical_stage": row.get("maxClinicalStage"),
        })
        for m in ((drug.get("mechanismsOfAction") or {}).get("rows") or []):
            if not m.get("mechanismOfAction"):
                continue
            out["drug_mechanism"].append({
                "chembl_id": drug["id"],
                "mechanism_of_action": m["mechanismOfAction"],
                "action_type": m.get("actionType"), "target_name": m.get("targetName"),
            })
        for rep in row.get("clinicalReports") or []:
            out["clinical_report"].append({
                "report_id": rep["id"], "source": rep.get("source"),
                "url": rep.get("url"), "trial_phase": rep.get("trialPhase"),
                "trial_status": rep.get("trialOverallStatus"),
            })
            out["target_drug_report"].append({
                "uniprot_acc": uniprot_acc, "chembl_id": drug["id"],
                "report_id": rep["id"],
            })

    return out
