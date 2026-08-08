/* Governance rules, tested.
 *
 *   node --test tests/
 *
 * The decision store is the one place where a user's judgement enters the
 * system, so its invariants — append-only, supersede-never-overwrite, no
 * decision without a rationale, preconditions that state themselves — deserve
 * the same standard as the triage rules in core/.
 */

import test from 'node:test';
import assert from 'node:assert/strict';

// decisions.js persists to localStorage, which Node does not have. Stub it
// before the module is imported; a dynamic import runs after this executes,
// whereas a static import would be hoisted above it.
globalThis.localStorage = {
  _s: new Map(),
  getItem(k) { return this._s.has(k) ? this._s.get(k) : null; },
  setItem(k, v) { this._s.set(k, String(v)); },
};

const DEC = await import('../site/decisions.js');

const BUILD = '2026-08-07T23:16:16+00:00';
const base = {
  objectType: 'variant', objectKey: 'v1', symbol: 'PTEN',
  outcome: 'endorsed', analyst: 'tester', priorClass: 'reclass_upgrade',
  dataBuild: BUILD,
};

function reset() { DEC.clear(); DEC.load(); }

/* ── append and supersede ──────────────────────────────────────────── */

test('a second decision supersedes the first, and both remain', () => {
  reset();
  const a = DEC.record({ ...base, rationale: 'first call' });
  const b = DEC.record({ ...base, outcome: 'rejected', rationale: 'revised' });

  assert.equal(DEC.all().length, 2, 'nothing is removed');
  assert.equal(a.superseded_by, b.id);
  assert.equal(b.supersedes, a.id);
  assert.equal(b.superseded_by, null);
  assert.equal(DEC.current('v1').id, b.id);
});

test('superseding never rewrites the superseded record', () => {
  reset();
  const a = DEC.record({ ...base, rationale: 'original wording' });
  const snapshot = { ...a };
  DEC.record({ ...base, outcome: 'rejected', rationale: 'different wording' });

  for (const field of ['rationale', 'outcome', 'analyst', 'decided_at', 'id']) {
    assert.equal(a[field], snapshot[field], `${field} must be immutable`);
  }
});

test('decisions on different objects do not supersede each other', () => {
  reset();
  DEC.record({ ...base, objectKey: 'v1', rationale: 'one' });
  DEC.record({ ...base, objectKey: 'v2', rationale: 'two' });
  assert.equal(DEC.current('v1').superseded_by, null);
  assert.equal(DEC.current('v2').superseded_by, null);
});

/* ── a decision must be justifiable ────────────────────────────────── */

test('a decision without a rationale is refused', () => {
  reset();
  for (const rationale of ['', '   ', '\n', undefined]) {
    assert.throws(() => DEC.record({ ...base, rationale }), /rationale is required/);
  }
  assert.equal(DEC.all().length, 0, 'a refused decision is not recorded');
});

test('a decision without an analyst is refused', () => {
  reset();
  assert.throws(() => DEC.record({ ...base, rationale: 'ok', analyst: '' }),
    /analyst name is required/);
});

/* ── preconditions state themselves ────────────────────────────────── */

test('only reviewable triage classes may be decided, and blocking says why', () => {
  for (const triage of ['reclass_upgrade', 'reclass_downgrade', 'discordant']) {
    assert.equal(DEC.variantPreconditions({ triage }).ok, true, triage);
  }
  for (const triage of ['concordant', 'not_triaged', 'unasserted']) {
    const pre = DEC.variantPreconditions({ triage });
    assert.equal(pre.ok, false, triage);
    assert.ok(pre.reasons[0].includes(triage), 'the reason names the actual class');
  }
  assert.equal(DEC.variantPreconditions(null).ok, false);
});

test('a target needs four of five evidence sources before it can be shortlisted', () => {
  const full = {
    associations: [1], tractability: [1], prioritisation: { a: 1 },
    structures: [1], drugs: [1],
  };
  assert.equal(DEC.targetPreconditions(full).ok, true);

  // BRCA1's real shape: everything but chemical matter.
  const noDrugs = { ...full, drugs: [] };
  assert.equal(DEC.targetPreconditions(noDrugs).ok, true,
    'four of five is enough — a tumour suppressor with no drugs is still decidable');

  const thin = { ...full, drugs: [], structures: [] };
  const pre = DEC.targetPreconditions(thin);
  assert.equal(pre.ok, false);
  assert.ok(pre.reasons[0].includes('chemical matter'));
  assert.ok(pre.reasons[0].includes('experimental structures'));
});

/* ── lineage ───────────────────────────────────────────────────────── */

test('a decision taken against an older data build is marked stale', () => {
  reset();
  DEC.record({ ...base, rationale: 'current', dataBuild: BUILD });
  DEC.record({ ...base, objectKey: 'v2', rationale: 'old',
               dataBuild: '2026-01-01T00:00:00+00:00' });

  assert.equal(DEC.staleCount(BUILD), 1);
  assert.equal(DEC.isStale(DEC.current('v2'), BUILD), true);
  assert.equal(DEC.isStale(DEC.current('v1'), BUILD), false);
});

test('a superseded stale decision is not counted — only live ones matter', () => {
  reset();
  DEC.record({ ...base, rationale: 'old', dataBuild: '2020-01-01T00:00:00Z' });
  DEC.record({ ...base, rationale: 'redone against current data', dataBuild: BUILD });
  assert.equal(DEC.staleCount(BUILD), 0);
});

/* ── portability ───────────────────────────────────────────────────── */

test('re-importing the same log is a no-op', () => {
  reset();
  DEC.record({ ...base, rationale: 'one' });
  DEC.record({ ...base, objectKey: 'v2', rationale: 'two' });
  const payload = DEC.exportPayload(BUILD);

  assert.equal(DEC.importPayload(payload), 0, 'nothing new on re-import');
  assert.equal(DEC.all().length, 2);

  reset();
  assert.equal(DEC.importPayload(payload), 2, 'restores into an empty store');
  assert.equal(DEC.all().length, 2);
});

test('a foreign file is rejected rather than silently merged', () => {
  reset();
  assert.throws(() => DEC.importPayload({ kind: 'something.else', decisions: [] }),
    /not a Locus decision log/);
  assert.throws(() => DEC.importPayload(null), /not a Locus decision log/);
});

test('the export carries the build its decisions were made against', () => {
  reset();
  DEC.record({ ...base, rationale: 'x' });
  const p = DEC.exportPayload(BUILD);
  assert.equal(p.kind, 'locus.decisions');
  assert.equal(p.data_build, BUILD);
  assert.equal(p.decisions.length, 1);
  assert.ok(p.exported_at);
});

/* ── vocabulary ────────────────────────────────────────────────────── */

test('outcomes are review verdicts, not clinical classifications', () => {
  // Locus surfaces candidates for expert review and never asserts
  // pathogenicity, so no outcome may read as a clinical call.
  const forbidden = /pathogenic|benign|vus|likely/i;
  for (const key of Object.keys(DEC.VARIANT_OUTCOMES)) {
    assert.ok(!forbidden.test(key), `outcome "${key}" reads as a clinical assertion`);
  }
});
