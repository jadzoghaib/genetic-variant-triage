/* View-model tests — the browser-side logic, without a browser.
 *
 *   node --test tests/view-model.test.mjs
 *
 * These cover the parts of the console that can be wrong in an interesting way:
 * payload decoding, filtering, binning, chart geometry, and the virtualisation
 * window. Rendering itself is left to the eye; arithmetic is not.
 */

import test from 'node:test';
import assert from 'node:assert/strict';

import * as VM from '../site/view-model.js';

/* ── payload decoding ──────────────────────────────────────────────── */

const payload = {
  n: 3,
  levels: {
    amc: ['LPath', 'LBen'], sig: ['VUS', 'PATH'], raw: ['Uncertain_significance'],
    triage: ['reclass_upgrade', 'concordant'], priority: ['high', 'none'],
    tier: ['experimental', 'predicted_weak'], reasons: ['class=reclass_upgrade'],
  },
  cols: {
    pv: ['I33N', 'C61G', 'P871L'], vid: ['v1', 'v2', 'v3'], pos: [33, 61, 871],
    am: [0.9999, 0.99, 0.042], plddt: [98.6, 87.6, 34.6],
    stars: [2, 3, -1], solved: [1, 1, 0],
    amc: [0, 0, 1], sig: [0, 1, -1], raw: [0, -1, -1],
    triage: [0, 1, 1], priority: [0, 1, 1], tier: [0, 0, 1],
    reasons: [0, -1, -1],
  },
};

test('columnar decode restores rows, mapping -1 to null', () => {
  const r = VM.rows(payload);
  assert.equal(r.length, 3);
  assert.deepEqual(
    { pv: r[0].pv, amc: r[0].amc, sig: r[0].sig, tier: r[0].tier, stars: r[0].stars },
    { pv: 'I33N', amc: 'LPath', sig: 'VUS', tier: 'experimental', stars: 2 });
  assert.equal(r[2].sig, null, 'absent ClinVar assertion decodes to null');
  assert.equal(r[2].stars, null, 'sentinel -1 for stars is null, not 0');
  assert.equal(r[1].reasons, null);
  assert.equal(r[0].reasons, 'class=reclass_upgrade');
  assert.equal(r[0].i, 0, 'row keeps its index into the unfiltered array');
});

/* ── filtering ─────────────────────────────────────────────────────── */

const F = { actionable: false, priority: '', tier: '', minStars: 0, search: '' };
const all = VM.rows(payload);

test('actionable-only keeps exactly the reviewable classes', () => {
  const out = VM.applyFilters(all, { ...F, actionable: true });
  assert.deepEqual(out.map((r) => r.pv), ['I33N']);
});

test('priority filter is a threshold, not an equality match', () => {
  assert.equal(VM.applyFilters(all, { ...F, priority: 'high' }).length, 1);
  assert.equal(VM.applyFilters(all, { ...F, priority: 'low' }).length, 1,
    'only the high row clears a low threshold; the others are "none"');
  assert.equal(VM.applyFilters(all, { ...F, priority: '' }).length, 3);
});

test('minimum stars treats an absent rating as zero, not as passing', () => {
  assert.deepEqual(VM.applyFilters(all, { ...F, minStars: 3 }).map((r) => r.pv),
    ['C61G']);
  assert.equal(VM.applyFilters(all, { ...F, minStars: 1 }).length, 2,
    'the unrated row is excluded');
});

test('search bypasses every other filter', () => {
  // P871L is concordant and unrated: actionable-only plus a star floor would
  // hide it. A name lookup must still find it.
  const out = VM.applyFilters(all,
    { actionable: true, priority: 'high', tier: 'experimental', minStars: 4,
      search: 'p871l' });
  assert.deepEqual(out.map((r) => r.pv), ['P871L'], 'and it is case-insensitive');
});

test('search matches substrings, so a residue number finds its variants', () => {
  assert.deepEqual(VM.applyFilters(all, { ...F, search: '33' }).map((r) => r.pv),
    ['I33N']);
});

test('worklist order is priority first, then model confidence', () => {
  const sorted = VM.sortWorklist(all);
  assert.equal(sorted[0].pv, 'I33N', 'high priority leads');
  assert.equal(sorted[1].pv, 'C61G', 'then the stronger model call');
  assert.deepEqual(all.map((r) => r.pv), ['I33N', 'C61G', 'P871L'],
    'sorting does not mutate the input');
});

/* ── 3D colouring ──────────────────────────────────────────────────── */

const profile = {
  n: 4,
  levels: { amc: ['LPath', 'LBen'], tier: ['experimental', 'predicted_weak'] },
  amc: [0, 1, 1, -1],
  tier: [0, 0, 1, 1],
};

test('colour buckets group residues and use 1-based numbering', () => {
  const b = VM.colorBuckets(profile, 'am',
    { LPath: '#e66767', LBen: '#3987e5', _none: '#2c2c2a' });
  assert.deepEqual(b['#e66767'], [1]);
  assert.deepEqual(b['#3987e5'], [2, 3]);
  assert.deepEqual(b['#2c2c2a'], [4], 'a residue with no prediction is not coloured');
});

test('tier colouring reads the tier the pipeline decided, not a pLDDT cut', () => {
  const b = VM.colorBuckets(profile, 'tier',
    { experimental: '#cde2fb', predicted_weak: '#256abf', _none: '#000' });
  assert.deepEqual(b['#cde2fb'], [1, 2]);
  assert.deepEqual(b['#256abf'], [3, 4]);
});

test('a focused view frames a window, never a bare residue', () => {
  const w = VM.focusWindow(858, 30);
  assert.equal(w.length, 61);
  assert.equal(w[0], 828);
  assert.equal(w.at(-1), 888);
  assert.equal(VM.focusWindow(5, 30)[0], 1, 'clamped at the N-terminus');
});

/* ── profile geometry ──────────────────────────────────────────────── */

test('binning averages and ignores gaps', () => {
  const b = VM.binSeries([1, 3, null, 5], 4, 2);
  assert.equal(b.length, 2);
  assert.equal(b[0].mean, 2, '(1+3)/2');
  assert.equal(b[1].mean, 5, 'the null is skipped, not counted as zero');
});

test('an all-empty bin is zero rather than NaN', () => {
  assert.equal(VM.binSeries([null, null], 2, 1)[0].mean, 0);
});

test('the area path closes on the baseline and stays inside the plot', () => {
  const binned = VM.binSeries([0, 1], 2, 2);
  const d = VM.areaPath(binned, { pad: 38, plot: 1000, y0: 12, height: 70, max: 1 });
  assert.ok(d.startsWith('M38,82'), 'opens on the baseline at the left gutter');
  assert.ok(d.endsWith('Z'), 'and closes');
  const xs = [...d.matchAll(/([\d.]+),[\d.]+/g)].map((m) => Number(m[1]));
  assert.ok(Math.min(...xs) >= 38 && Math.max(...xs) <= 1038,
    'no point escapes the gutter or the right edge');
});

test('ticks are round numbers inside the sequence', () => {
  const t = VM.ticks(1863, 38, 1000);
  assert.ok(t.every((x) => x.value % 50 === 0), 'every tick is a multiple of 50');
  assert.ok(t.every((x) => x.value < 1863 && x.x > 38));
});

/* ── spatial neighbourhood ─────────────────────────────────────────── */

// A deliberately contrived fold: residues 1-3 sit together, and residue 200 is
// folded back against them despite being 197 positions away in sequence.
// Residue 500 is off on its own.
const atoms = [
  { resi: 1, x: 0, y: 0, z: 0 },
  { resi: 2, x: 3, y: 0, z: 0 },
  { resi: 3, x: 0, y: 3, z: 0 },
  { resi: 200, x: 0, y: 0, z: 5 },
  { resi: 500, x: 50, y: 50, z: 50 },
];

test('neighbours are found by distance, not by sequence position', () => {
  assert.deepEqual(VM.spatialNeighbours(atoms, 1, 8), [2, 3, 200]);
  assert.deepEqual(VM.spatialNeighbours(atoms, 500, 8), [],
    'an isolated residue has no contacts');
});

test('the cutoff is respected exactly at the boundary', () => {
  assert.deepEqual(VM.spatialNeighbours(atoms, 1, 3), [2, 3],
    'a residue exactly at the cutoff counts');
  assert.deepEqual(VM.spatialNeighbours(atoms, 1, 2.9), []);
});

test('a residue with no coordinates returns null, not an empty answer', () => {
  // Missing is not the same as isolated, and the caller must be able to tell.
  assert.equal(VM.spatialNeighbours(atoms, 999, 8), null);
  assert.equal(VM.neighbourhood(atoms, 999, () => 'LPath'), null);
});

test('the neighbourhood separates damaging contacts the sequence hides', () => {
  // Everything nearby is damaging; only residue 200 is far away in sequence.
  const n = VM.neighbourhood(atoms, 1, () => 'LPath');
  assert.deepEqual(n.near, [2, 3, 200]);
  assert.deepEqual(n.damaging, [2, 3, 200]);
  assert.deepEqual(n.distant, [200],
    'the whole point: a contact 199 positions away is invisible in a profile');
  assert.equal(n.span, 200);
});

test('benign neighbours are not counted as damaging', () => {
  const classOf = (r) => (r === 200 ? 'LPath' : 'LBen');
  const n = VM.neighbourhood(atoms, 1, classOf);
  assert.deepEqual(n.damaging, [200]);
  assert.deepEqual(n.distant, [200]);
});

test('an adjacent-only cluster reports nothing distant', () => {
  // Residues 1-3 touching each other is what any profile already shows, so it
  // must not be presented as a 3D discovery.
  const local = atoms.filter((a) => a.resi < 10);
  const n = VM.neighbourhood(local, 1, () => 'LPath');
  assert.deepEqual(n.distant, []);
});

/* ── virtualisation ────────────────────────────────────────────────── */

test('the window draws a screenful and pads the rest', () => {
  const w = VM.windowRows(12463, 0, 430, 26, 8);
  assert.equal(w.start, 0);
  assert.ok(w.end < 60, 'a screenful plus overscan, not 12,463 rows');
  assert.equal(w.padTop, 0);
  assert.equal(w.padBottom, (12463 - w.end) * 26);
});

test('scrolled window keeps total height constant, so the scrollbar is honest', () => {
  const H = 26, total = 12463;
  for (const scrollTop of [0, 500, 9000, total * H - 430]) {
    const w = VM.windowRows(total, scrollTop, 430, H, 8);
    const drawn = (w.end - w.start) * H;
    assert.equal(w.padTop + drawn + w.padBottom, total * H,
      `heights must sum at scrollTop=${scrollTop}`);
    assert.ok(w.start >= 0 && w.end <= total, 'window stays in range');
  }
});

test('an empty worklist virtualises to nothing', () => {
  assert.deepEqual(VM.windowRows(0, 0, 430, 26),
    { start: 0, end: 0, padTop: 0, padBottom: 0 });
});
