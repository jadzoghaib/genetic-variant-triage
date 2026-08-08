/* Governed actions — the decision store.
 *
 * Append and supersede, never update in place. The only field ever written on
 * an existing record is `superseded_by`; a record's class, rationale, author
 * and timestamp are immutable once appended. That mirrors the triage_decision
 * table exactly, so an exported log re-ingests into the ontology unchanged.
 *
 * Every decision records the data build it was made against. If the underlying
 * ClinVar or AlphaMissense data is rebuilt, a decision taken against the older
 * evidence must not silently present itself as current — it is marked stale.
 * That is the same lineage discipline the rest of the system applies to facts,
 * applied to judgements.
 *
 * HONESTY: this is browser storage. It records who claims to have decided what,
 * which is attribution, not authentication, and it is not tamper-proof. The
 * exported file is the durable artifact.
 */

const KEY = 'locus.decisions.v1';

/** Review outcomes, deliberately NOT clinical classifications. Locus surfaces
 *  candidates for expert review; it never asserts pathogenicity itself. */
export const VARIANT_OUTCOMES = {
  endorsed: 'Endorse — model call looks right, escalate to formal curation',
  rejected: 'Reject — model call looks wrong',
  needs_evidence: 'Needs more evidence',
  deferred: 'Defer',
};

export const TARGET_OUTCOMES = {
  shortlisted: 'Shortlist for pursuit',
  deprioritized: 'Deprioritize',
};

let log = [];

/* ── persistence ───────────────────────────────────────────────────── */

export function load() {
  try {
    log = JSON.parse(localStorage.getItem(KEY) || '[]');
  } catch {
    log = [];
  }
  return log;
}

function persist() {
  localStorage.setItem(KEY, JSON.stringify(log));
}

export function all() { return log; }

/** The live decision for an object, if one has not been superseded. */
export function current(objectKey) {
  return log.find((d) => d.object_key === objectKey && d.superseded_by === null) || null;
}

export function history(objectKey) {
  return log.filter((d) => d.object_key === objectKey)
            .sort((a, b) => b.decided_at.localeCompare(a.decided_at));
}

/* ── preconditions ─────────────────────────────────────────────────── */

const ACTIONABLE = new Set(['reclass_upgrade', 'reclass_downgrade', 'discordant']);

/** Whether a variant may be decided on, and — when it may not — why.
 *  The reason is surfaced in the UI: a disabled control that does not say what
 *  it wants is not a governed action, just a broken one. */
export function variantPreconditions(row) {
  const reasons = [];
  if (!row) reasons.push('no variant selected');
  else if (!ACTIONABLE.has(row.triage)) {
    reasons.push(`triage class is "${row.triage}" — only reclassification `
      + 'candidates and model/curator conflicts are reviewable');
  }
  return { ok: reasons.length === 0, reasons };
}

/** A target may be shortlisted once its dossier actually rests on something:
 *  at least four of the five evidence sources returned data. */
export function targetPreconditions(dossier) {
  const present = [
    ['disease associations', dossier.associations?.length],
    ['tractability', dossier.tractability?.length],
    ['prioritisation', Object.keys(dossier.prioritisation || {}).length],
    ['experimental structures', dossier.structures?.length],
    ['chemical matter', dossier.drugs?.length],
  ];
  const missing = present.filter(([, n]) => !n).map(([name]) => name);
  const ok = present.length - missing.length >= 4;
  return {
    ok,
    reasons: ok ? [] : [`dossier incomplete — no ${missing.join(', no ')}`],
    missing,
  };
}

/* ── recording ─────────────────────────────────────────────────────── */

function id() {
  return 'd_' + Math.random().toString(36).slice(2, 10) + Date.now().toString(36);
}

/**
 * Append a decision, superseding any live decision on the same object.
 * Returns the new record.
 */
export function record({ objectType, objectKey, symbol, outcome, rationale,
                         analyst, priorClass, dataBuild }) {
  const text = (rationale || '').trim();
  if (!text) throw new Error('a rationale is required — a decision without one is not reviewable');
  if (!analyst) throw new Error('an analyst name is required');

  const prev = current(objectKey);
  const rec = {
    id: id(),
    object_type: objectType,
    object_key: objectKey,
    symbol,
    outcome,
    rationale: text,
    analyst,
    prior_class: priorClass ?? null,
    supersedes: prev ? prev.id : null,
    superseded_by: null,
    decided_at: new Date().toISOString(),
    data_build: dataBuild,
  };
  // The one permitted mutation of an existing record.
  if (prev) prev.superseded_by = rec.id;
  log.unshift(rec);
  persist();
  return rec;
}

/* ── lineage ───────────────────────────────────────────────────────── */

export function isStale(rec, dataBuild) {
  return !!rec.data_build && !!dataBuild && rec.data_build !== dataBuild;
}

export function staleCount(dataBuild) {
  return log.filter((d) => d.superseded_by === null && isStale(d, dataBuild)).length;
}

/* ── portability ───────────────────────────────────────────────────── */

export function exportPayload(dataBuild) {
  return {
    kind: 'locus.decisions',
    version: 1,
    exported_at: new Date().toISOString(),
    data_build: dataBuild,
    decisions: log,
  };
}

/** Merge an exported log. Records are matched by id, so re-importing the same
 *  file is a no-op rather than a duplication. */
export function importPayload(payload) {
  if (!payload || payload.kind !== 'locus.decisions') {
    throw new Error('not a Locus decision log');
  }
  const seen = new Set(log.map((d) => d.id));
  const added = (payload.decisions || []).filter((d) => !seen.has(d.id));
  log = [...added, ...log].sort((a, b) => b.decided_at.localeCompare(a.decided_at));
  persist();
  return added.length;
}

export function clear() { log = []; persist(); }
