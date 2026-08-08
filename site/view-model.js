/* View model — pure functions over the exported payloads.
 *
 * No DOM, no fetch, no 3Dmol. This is the browser-side counterpart of core/ in
 * Python: the parts of the console that can be wrong in an interesting way,
 * separated from the parts that only move pixels, so they can be tested
 * without a browser.
 *
 * Nothing here re-derives a scientific threshold. Tiers and AlphaMissense
 * classes arrive already decided from the Python pipeline; this module maps
 * codes to labels and colours and does layout arithmetic. That boundary is the
 * point — a threshold living in two languages is a threshold that drifts.
 */

export const TIER_LABEL = {
  experimental: 'experimentally solved',
  predicted_confident: 'confidently predicted',
  predicted_weak: 'predicted, unreliable',
};
export const TIER_ORDER = ['experimental', 'predicted_confident', 'predicted_weak'];
export const PRIORITY_RANK = { high: 3, medium: 2, low: 1, none: 0 };
export const ACTIONABLE = new Set(
  ['reclass_upgrade', 'reclass_downgrade', 'discordant']);

/** Expand a columnar payload into row objects, decoding categorical codes. */
export function rows(v) {
  const c = v.cols, L = v.levels, out = new Array(v.n);
  const dec = (levels, code) => (code < 0 ? null : levels[code]);
  for (let i = 0; i < v.n; i++) {
    out[i] = {
      i,
      pv: c.pv[i], vid: c.vid[i], pos: c.pos[i],
      am: c.am[i], plddt: c.plddt[i],
      stars: c.stars[i] < 0 ? null : c.stars[i],
      solved: !!c.solved[i],
      amc: dec(L.amc, c.amc[i]),
      sig: dec(L.sig, c.sig[i]),
      raw: dec(L.raw, c.raw[i]),
      triage: dec(L.triage, c.triage[i]),
      priority: dec(L.priority, c.priority[i]),
      tier: dec(L.tier, c.tier[i]),
      reasons: dec(L.reasons, c.reasons[i]),
    };
  }
  return out;
}

/**
 * Filter the worklist.
 *
 * A name lookup is not a filter refinement. Someone typing "L858R" wants that
 * variant even when the worklist filters would hide it — ClinVar classifies
 * L858R as drug_response, so it is not_triaged and carries no priority. Search
 * therefore bypasses every other filter rather than intersecting with them.
 */
export function applyFilters(all, f) {
  if (f.search) {
    const q = f.search.trim().toUpperCase();
    return all.filter((r) => r.pv.toUpperCase().includes(q));
  }
  return all.filter((r) => {
    if (f.actionable && !ACTIONABLE.has(r.triage)) return false;
    if (f.priority && PRIORITY_RANK[r.priority] < PRIORITY_RANK[f.priority]) return false;
    if (f.tier && r.tier !== f.tier) return false;
    if (f.minStars && (r.stars ?? 0) < f.minStars) return false;
    return true;
  });
}

/** Worklist order: priority descending, then model confidence descending. */
export function sortWorklist(view) {
  return [...view].sort((a, b) =>
    (PRIORITY_RANK[b.priority] - PRIORITY_RANK[a.priority])
    || ((b.am ?? 0) - (a.am ?? 0)));
}

/* ── 3D colouring ──────────────────────────────────────────────────── */

/** Group residue indices by colour, so the viewer takes one styling call per
 *  bucket instead of one per residue. */
export function colorBuckets(profile, mode, colors) {
  const codes = mode === 'am' ? profile.amc : profile.tier;
  const levels = mode === 'am' ? profile.levels.amc : profile.levels.tier;
  const buckets = {};
  for (let i = 0; i < profile.n; i++) {
    const key = codes[i] < 0 ? null : levels[codes[i]];
    const color = colors[key] ?? colors._none;
    (buckets[color] ||= []).push(i + 1);   // residue numbering is 1-based
  }
  return buckets;
}

/** The residue window a focused view should frame. Fitting to a single residue
 *  loses every bit of structural context that makes the view worth showing. */
export function focusWindow(pos, span = 30) {
  const lo = Math.max(1, pos - span);
  return Array.from({ length: span * 2 + 1 }, (_, k) => lo + k);
}

/* ── sequence profile geometry ─────────────────────────────────────── */

/** Average a per-residue series into `bins` windows.
 *  A protein is longer than the plot is wide; drawing 1,863 points into ~1,000
 *  pixels is over-plotting, and the regional signal the profile exists to show
 *  disappears into a hairball. */
export function binSeries(values, n, bins) {
  const size = Math.ceil(n / bins);
  const out = [];
  for (let s = 0; s < n; s += size) {
    let sum = 0, count = 0;
    for (let i = s; i < Math.min(s + size, n); i++) {
      if (values[i] != null) { sum += values[i]; count++; }
    }
    out.push({ x: s / n, mean: count ? sum / count : 0 });
  }
  return out;
}

/** SVG path for a filled area over binned values. */
export function areaPath(binned, { pad, plot, y0, height, max }) {
  const xat = (frac) => pad + frac * plot;
  const pts = binned.map((b, i) =>
    `${xat(i / (binned.length - 1)).toFixed(1)},`
    + `${(y0 + height - (b.mean / max) * height).toFixed(1)}`);
  return `M${pad},${y0 + height} L${pts.join(' L')} L${pad + plot},${y0 + height} Z`;
}

/** Axis tick positions, rounded to a readable step. */
export function ticks(n, pad, plot, target = 10) {
  const step = Math.max(50, Math.round(n / target / 50) * 50);
  const out = [];
  for (let v = step; v < n; v += step) {
    out.push({ value: v, x: pad + (v / n) * plot });
  }
  return out;
}

/* ── virtualisation ────────────────────────────────────────────────── */

/**
 * Which rows to draw for a scroll position, plus the spacer heights that keep
 * the scrollbar honest. Removes the row cap: the table can hold every variant
 * while only ever putting a screenful of <tr> in the document.
 */
export function windowRows(total, scrollTop, viewportHeight, rowHeight, overscan = 8) {
  if (total === 0) return { start: 0, end: 0, padTop: 0, padBottom: 0 };
  const visible = Math.ceil(viewportHeight / rowHeight);
  const start = Math.max(0, Math.floor(scrollTop / rowHeight) - overscan);
  const end = Math.min(total, start + visible + overscan * 2);
  return {
    start,
    end,
    padTop: start * rowHeight,
    padBottom: Math.max(0, (total - end) * rowHeight),
  };
}
