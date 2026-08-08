/* Locus console — rendering only.
 *
 * Every rule shown here was decided in Python and covered by the test suite;
 * the arithmetic is in view-model.js and covered too. This file moves pixels.
 *
 * Two things it deliberately does NOT do: re-derive a scientific threshold
 * (tiers and AlphaMissense classes arrive already decided), and rebuild the 3D
 * viewers per interaction (they are created once, lazily, and restyled).
 */

import * as DEC from './decisions.js';
import * as VM from './view-model.js';

const PALETTE = {
  am: { LBen: '#3987e5', Amb: '#383835', LPath: '#e66767', _none: '#2c2c2a' },
  tier: {
    experimental: '#cde2fb', predicted_confident: '#6da7ec',
    predicted_weak: '#256abf', _none: '#2c2c2a',
  },
  priority: { high: '#d03b3b', medium: '#ec835a', low: '#fab219', none: '#898781' },
  series1: '#3987e5', muted: '#898781', ink2: '#c3c2b7', surface: '#1a1a19',
};

const ROW_H = 26;        // must match tbody td height in styles.css
const state = {
  view: 'console', target: null, selected: null, cursor: 0,
  analyst: localStorage.getItem('locus.analyst') || '',
  filters: { actionable: true, priority: 'medium', tier: '', minStars: 0, search: '' },
};

const cache = new Map();
let manifest = null;
let worklistView = [];
let worklistAll = [];
const viewers = { am: null, tier: null };
let currentStructureId = null;
let structureVisible = false;
let pendingStructure = null;

const $ = (s) => document.querySelector(s);
const fmt = (n) => n.toLocaleString('en-US');
const esc = (s) => String(s ?? '').replace(/[&<>"]/g,
  (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));

/* ── data ──────────────────────────────────────────────────────────── */

async function load(name) {
  if (!cache.has(name)) {
    const r = await fetch(`data/${name}`);
    if (!r.ok) throw new Error(`${name}: HTTP ${r.status}`);
    cache.set(name, name.endsWith('.pdb') ? await r.text() : await r.json());
  }
  return cache.get(name);
}

/* ── worklist table (virtualised) ──────────────────────────────────── */

/** A coloured mark carries identity; the text beside it stays in ink.
 *  Series colours as small text fail contrast — measured 3.23:1 for the
 *  weakest tier against this surface, against a 4.5 requirement. */
const swatch = (color) => `<i class="dot" style="background:${color}"></i>`;

function rowHTML(r, viewIndex) {
  const isCursor = viewIndex === state.cursor;
  return `<tr data-i="${r.i}" data-v="${viewIndex}" role="row"
      aria-rowindex="${viewIndex + 2}" aria-selected="${state.selected === r.i}"
      tabindex="${isCursor ? 0 : -1}">
    <td>${r.pv}</td>
    <td>${r.triage}</td>
    <td>${swatch(PALETTE.priority[r.priority])}${r.priority}</td>
    <td>${r.sig ?? '<span class="dim">no assertion</span>'}</td>
    <td class="num">${r.stars ?? ''}</td>
    <td class="num"><i class="bar" style="width:${Math.round((r.am ?? 0) * 42)}px"></i>${(r.am ?? 0).toFixed(3)}</td>
    <td class="num">${r.plddt?.toFixed(1) ?? ''}</td>
    <td>${swatch(PALETTE.tier[r.tier])}${VM.TIER_LABEL[r.tier]}</td>
  </tr>`;
}

function drawWindow() {
  const wrap = $('#wl-wrap');
  const w = VM.windowRows(worklistView.length, wrap.scrollTop,
                          wrap.clientHeight || 430, ROW_H);
  const spacer = (h) => (h ? `<tr class="spacer" style="height:${h}px" aria-hidden="true"><td colspan="8"></td></tr>` : '');
  $('#worklist').setAttribute('aria-rowcount', String(worklistView.length + 1));
  $('#worklist tbody').innerHTML = worklistView.length
    ? spacer(w.padTop)
      + worklistView.slice(w.start, w.end).map((r, k) => rowHTML(r, w.start + k)).join('')
      + spacer(w.padBottom)
    : `<tr><td colspan="8" class="empty">No variants match these filters.</td></tr>`;
}

function moveCursor(i) {
  const n = worklistView.length;
  if (!n) return;
  state.cursor = Math.max(0, Math.min(n - 1, i));
  const wrap = $('#wl-wrap');
  const top = state.cursor * ROW_H;
  if (top < wrap.scrollTop) wrap.scrollTop = top;
  else if (top + ROW_H > wrap.scrollTop + wrap.clientHeight) {
    wrap.scrollTop = top + ROW_H - wrap.clientHeight;
  }
  drawWindow();
  $(`#worklist tbody tr[data-v="${state.cursor}"]`)?.focus();
}

function selectRow(allIndex) {
  state.selected = state.selected === allIndex ? null : allIndex;
  render();
}

/* ── detail ────────────────────────────────────────────────────────── */

function chip(text, color) {
  return `<span class="chip" style="color:${color};border-color:${color}55">${text}</span>`;
}

function renderDetail(r) {
  if (!r) {
    $('#detail').className = 'empty';
    $('#detail').innerHTML = 'Select a variant to see its residue, its structural '
      + 'evidence, and the retrieval that produced every field.';
    return;
  }
  $('#detail').className = 'detail';
  $('#detail').innerHTML = `
    <div class="title">${r.pv}</div>
    <div class="ident">${r.vid}</div>
    ${chip(r.triage, PALETTE.priority[r.priority] || PALETTE.muted)}
    ${chip(VM.TIER_LABEL[r.tier], PALETTE.tier[r.tier])}
    <div class="kv">
      <div class="k">ClinVar</div><div class="v">${r.sig ?? 'no assertion'} (${r.stars ?? 0}★)</div>
      <div class="k">raw label</div><div class="v">${esc(r.raw) || '—'}</div>
      <div class="k">AlphaMissense</div><div class="v">${r.am?.toFixed(4)} · ${r.amc}</div>
      <div class="k">residue</div><div class="v">${r.pos} · pLDDT ${r.plddt?.toFixed(1)}</div>
      <div class="k">solved</div><div class="v">${r.solved ? 'yes' : 'no'}</div>
      <div class="k">priority</div><div class="v">${r.priority}</div>
    </div>
    ${r.reasons ? `<div class="ident" style="margin-top:.6rem">${esc(r.reasons)}</div>` : ''}
    <div id="v-action"></div>`;
  renderVariantActions(r);
}

/* ── governed actions ──────────────────────────────────────────────── */

function decisionBlock(rec, build) {
  const stale = DEC.isStale(rec, build);
  return `<div class="decision${rec.superseded_by ? ' superseded' : ''}">
    <b>${esc(DEC.VARIANT_OUTCOMES[rec.outcome] || DEC.TARGET_OUTCOMES[rec.outcome]
             || rec.outcome)}</b>
    ${rec.superseded_by ? '<span class="tag">superseded</span>' : ''}
    ${stale ? '<span class="tag stale">earlier data build</span>' : ''}
    <div class="why">${esc(rec.rationale)}</div>
    <div class="meta">${esc(rec.analyst)} · ${rec.decided_at.replace('T', ' ').slice(0, 19)}Z</div>
  </div>`;
}

const outcomeOptions = (map) => Object.entries(map)
  .map(([k, label]) => `<option value="${k}">${esc(label)}</option>`).join('');

function renderVariantActions(r) {
  const el = $('#v-action');
  if (!el) return;
  const pre = DEC.variantPreconditions(r);
  const past = DEC.history(r.vid);
  const build = manifest.generated_at;

  el.innerHTML = `
    <div class="eyebrow" style="margin-top:1rem">Governed action</div>
    ${pre.ok ? `
      <label class="sr-only" for="v-outcome">Review outcome</label>
      <select id="v-outcome">${outcomeOptions(DEC.VARIANT_OUTCOMES)}</select>
      <div class="row">
        <label class="sr-only" for="v-rationale">Rationale</label>
        <textarea id="v-rationale" placeholder="Rationale — required. What did you check, and what did it show?"></textarea>
      </div>
      <div class="row">
        <button class="btn primary" id="v-record">Record decision</button>
        <span class="stamp" id="v-msg" role="status"></span>
      </div>`
    : `<div class="blocked">Not reviewable: ${pre.reasons.map(esc).join('; ')}.</div>`}
    ${past.length ? `<div class="eyebrow" style="margin-top:0.9rem">Decision history</div>
      ${past.map((d) => decisionBlock(d, build)).join('')}` : ''}`;

  if (!pre.ok) return;
  $('#v-record').addEventListener('click', () => {
    try {
      DEC.record({
        objectType: 'variant', objectKey: r.vid, symbol: state.target,
        outcome: $('#v-outcome').value, rationale: $('#v-rationale').value,
        analyst: requireAnalyst(), priorClass: r.triage,
        dataBuild: manifest.generated_at,
      });
      render();
    } catch (err) { $('#v-msg').textContent = err.message; }
  });
}

function renderTargetActions(t, dossier) {
  const pre = DEC.targetPreconditions(dossier);
  const past = DEC.history(t.acc);
  const build = manifest.generated_at;

  $('#d-action').innerHTML = `
    ${pre.ok ? `
      <label class="sr-only" for="t-outcome">Decision</label>
      <select id="t-outcome">${outcomeOptions(DEC.TARGET_OUTCOMES)}</select>
      <div class="row">
        <label class="sr-only" for="t-rationale">Rationale</label>
        <textarea id="t-rationale" placeholder="Rationale — required."></textarea>
      </div>
      <div class="row">
        <button class="btn primary" id="t-record">Record decision</button>
        <span class="stamp" id="t-msg" role="status"></span>
      </div>`
    : `<div class="blocked">Blocked: ${pre.reasons.map(esc).join('; ')}. A dossier
        resting on fewer than four of five sources is not a basis for a
        portfolio decision.</div>`}
    ${past.length ? `<div class="eyebrow" style="margin-top:0.9rem">Decision history</div>
      ${past.map((d) => decisionBlock(d, build)).join('')}` : ''}`;

  if (!pre.ok) return;
  $('#t-record').addEventListener('click', () => {
    try {
      DEC.record({
        objectType: 'target', objectKey: t.acc, symbol: t.symbol,
        outcome: $('#t-outcome').value, rationale: $('#t-rationale').value,
        analyst: requireAnalyst(), priorClass: null,
        dataBuild: manifest.generated_at,
      });
      render();
    } catch (err) { $('#t-msg').textContent = err.message; }
  });
}

function requireAnalyst() {
  if (state.analyst) return state.analyst;
  const name = (prompt('Analyst name — recorded against every decision:') || '').trim();
  if (name) { state.analyst = name; localStorage.setItem('locus.analyst', name); }
  return name;
}

/* ── 3D ────────────────────────────────────────────────────────────── */

function paint(viewer, buckets, focus) {
  viewer.setStyle({}, { cartoon: { color: PALETTE.am._none } });
  for (const [color, resi] of Object.entries(buckets)) {
    viewer.setStyle({ resi }, { cartoon: { color } });
  }
  if (focus != null) {
    viewer.addStyle({ resi: [focus] }, { stick: { radius: 0.32, color: '#ffffff' } });
    viewer.addStyle({ resi: [focus] },
      { sphere: { radius: 1.0, color: '#ffffff', opacity: 0.55 } });
    viewer.zoomTo({ resi: VM.focusWindow(focus) });
  } else {
    viewer.zoomTo();
  }
  viewer.render();
}

async function renderStructure(t, profile, focus) {
  $('#struct-h').textContent = `Structure · ${t.structure_id}`
    + (focus != null ? ` · focused on residue ${focus}` : '');

  // Two WebGL contexts are the most expensive thing on the page, so they are
  // not created until the section is actually scrolled into view.
  if (!structureVisible) { pendingStructure = [t, profile, focus]; return; }

  if (!viewers.am) {
    viewers.am = $3Dmol.createViewer($('#mol-am'), { backgroundColor: PALETTE.surface });
    viewers.tier = $3Dmol.createViewer($('#mol-tier'), { backgroundColor: PALETTE.surface });
  }
  if (currentStructureId !== t.structure_id) {
    const pdb = await load(`${t.structure_id}.pdb`);
    for (const v of [viewers.am, viewers.tier]) {
      v.clear();
      v.addModel(pdb, 'pdb');
      // Only alpha carbons are clickable: one hit-test target per residue
      // rather than five, which both cuts the cost of every mouse move and
      // makes a click land on the residue the user aimed at.
      v.setClickable({ atom: 'CA' }, true, (atom) => selectResidue(atom.resi));
    }
    currentStructureId = t.structure_id;
  }
  viewers.am.resize();
  viewers.tier.resize();
  paint(viewers.am, VM.colorBuckets(profile, 'am', PALETTE.am), focus);
  paint(viewers.tier, VM.colorBuckets(profile, 'tier', PALETTE.tier), focus);
}

/** Structure → variant traversal: clicking a residue selects its most damaging
 *  variant in the worklist. */
function selectResidue(pos) {
  const hit = worklistAll.filter((r) => r.pos === pos)
    .sort((a, b) => (b.am ?? 0) - (a.am ?? 0))[0];
  if (!hit) return;
  state.filters.search = hit.pv;
  $('#f-search').value = hit.pv;
  state.selected = hit.i;
  render();
}

function observeStructure() {
  if (!('IntersectionObserver' in window)) { structureVisible = true; return; }
  const io = new IntersectionObserver((entries) => {
    if (!entries.some((e) => e.isIntersecting)) return;
    structureVisible = true;
    io.disconnect();
    if (pendingStructure) renderStructure(...pendingStructure);
  }, { rootMargin: '250px' });
  io.observe($('#mols'));
}

/* ── sequence profile ──────────────────────────────────────────────── */

function profileSVG(p) {
  const PAD = 38, W = 1038, PLOT = W - PAD, BINS = 200;
  const am = VM.binSeries(p.maxam, p.n, BINS);
  const pl = VM.binSeries(p.plddt, p.n, BINS);

  const shade = p.runs.filter((r) => r.tier === 'experimental').map((r) => {
    const x = PAD + (r.start / p.n) * PLOT, w = ((r.end - r.start) / p.n) * PLOT;
    return `<rect x="${x.toFixed(1)}" y="100" width="${w.toFixed(1)}" height="70"
             fill="${PALETTE.tier.experimental}" opacity="0.16"/>`;
  }).join('');

  const dash = (y, label) => `
    <line x1="${PAD}" y1="${y}" x2="${W}" y2="${y}" stroke="${PALETTE.muted}"
          stroke-width="0.7" stroke-dasharray="3 3" opacity="0.8"/>
    <text x="${PAD - 4}" y="${y + 3}" fill="${PALETTE.muted}" font-size="8.5"
          text-anchor="end">${label}</text>`;

  const ticks = VM.ticks(p.n, PAD, PLOT).map((t) =>
    `<text x="${t.x.toFixed(1)}" y="188" fill="${PALETTE.muted}" font-size="9"
           text-anchor="middle">${t.value}</text>`).join('');

  return `<svg viewBox="0 0 ${W} 198" width="100%" style="display:block"
       role="img" aria-label="Variant burden and structural confidence along the sequence">
    <text x="${PAD}" y="8" fill="${PALETTE.muted}" font-size="9">variant burden (mean max AlphaMissense)</text>
    <path d="${VM.areaPath(am, { pad: PAD, plot: PLOT, y0: 12, height: 70, max: 1 })}"
          fill="${PALETTE.series1}" opacity="0.9"/>
    ${dash(12 + 70 - 0.564 * 70, '0.564')}${dash(12 + 70 - 0.34 * 70, '0.34')}
    <text x="${PAD}" y="96" fill="${PALETTE.muted}" font-size="9">structural confidence (mean pLDDT)</text>
    ${shade}
    <path d="${VM.areaPath(pl, { pad: PAD, plot: PLOT, y0: 100, height: 70, max: 100 })}"
          fill="${PALETTE.series1}" opacity="0.9"/>
    ${dash(100 + 70 - 0.7 * 70, '70')}
    ${ticks}
    <text x="${PAD + PLOT / 2}" y="198" fill="${PALETTE.muted}" font-size="9"
          text-anchor="middle">residue position</text>
  </svg>`;
}

/* ── dossier ───────────────────────────────────────────────────────── */

const BAND_COLOR = {
  strong: '#0ca30c', moderate: '#fab219', weak: '#ec835a', absent: '#898781',
  approved: '#0ca30c', clinical: '#fab219', preclinical: '#ec835a', none: '#898781',
};

function renderDossier(t, d) {
  const c = d.card, dim = c.dimensions;
  $('#d-title').textContent = `${t.symbol} · ${t.acc}`;
  $('#d-sub').textContent = c.approved_name || '';
  const bands = [['genetic', dim.genetic_evidence.band],
                 ['structure', dim.structural_readiness.band],
                 ['pocket', dim.binding_site.band],
                 ['chemistry', dim.chemical_matter.band]]
    .map(([k, v]) => `${k}=<b style="color:${BAND_COLOR[v]}">${v}</b>`).join(' · ');
  $('#d-archetype').innerHTML =
    `<span class="a">${c.archetype}</span><span class="sub">${bands}</span>`;

  const meter = (v) => `<div class="meter"><i style="width:${Math.round((v ?? 0) * 100)}%"></i></div>`;
  const g = dim.genetic_evidence, s = dim.structural_readiness,
        b = dim.binding_site, ch = dim.chemical_matter, vb = dim.variant_burden;

  $('#d-cards').innerHTML = `
    <div class="card"><div class="k">Genetic evidence</div>
      <div class="band" style="color:${BAND_COLOR[g.band]}">${g.band}</div>
      <div class="d">max genetic score <b>${g.max_genetic_score ?? '—'}</b>${meter(g.max_genetic_score)}
        <div style="margin-top:.35rem">${esc(g.top_disease) || ''}</div>
        <div>${fmt(g.n_associated_diseases ?? 0)} associated diseases</div></div></div>
    <div class="card"><div class="k">Structural readiness</div>
      <div class="band" style="color:${BAND_COLOR[s.band]}">${s.band}</div>
      <div class="d"><b>${s.pct_residues_solved?.toFixed(1) ?? '—'}%</b> of residues solved
        ${meter((s.pct_residues_solved ?? 0) / 100)}
        <div style="margin-top:.35rem">${fmt(s.n_pdb_entities ?? 0)} PDB entities</div>
        <div>global pLDDT ${s.global_plddt?.toFixed(1) ?? '—'}</div></div></div>
    <div class="card"><div class="k">Binding site</div>
      <div class="band" style="color:${BAND_COLOR[b.band]}">${b.band}</div>
      <div class="d">pocket <b>${b.has_pocket ? 'yes' : 'no'}</b><br>
        ligand <b>${b.has_ligand ? 'yes' : 'no'}</b></div></div>
    <div class="card"><div class="k">Chemical matter</div>
      <div class="band" style="color:${BAND_COLOR[ch.band]}">${ch.band}</div>
      <div class="d"><b>${fmt(ch.n_drugs ?? 0)}</b> drugs<br>
        <b>${fmt(ch.n_trials ?? 0)}</b> trial reports</div></div>
    <div class="card"><div class="k">Variant burden</div>
      <div class="band" style="color:${PALETTE.series1}">${fmt(vb.upgrade_candidates ?? 0)}</div>
      <div class="d">reclassification-upgrade candidates<br>
        <b>${fmt(vb.discordant ?? 0)}</b> discordant with a curator</div></div>`;

  $('#d-assoc-n').textContent = `${d.associations.length} shown`;
  $('#d-assoc tbody').innerHTML = d.associations.map((a) => `
    <tr><td>${esc(a.disease)}</td><td class="num">${a.overall ?? ''}</td>
      <td class="num">${a.genetic ?? '—'}</td>
      <td class="num">${a.literature ?? '—'}</td>
      <td class="num">${a.known_drug ?? '—'}</td></tr>`).join('');

  $('#d-drugs-n').textContent = d.drugs.length ? `${d.drugs.length} drugs` : '';
  $('#d-drugs tbody').innerHTML = d.drugs.length
    ? d.drugs.map((x) => `
        <tr><td>${esc(x.drug ?? x.chembl_id)}</td><td>${esc(x.drug_type) || '—'}</td>
          <td>${esc(x.stage) || '—'}</td><td class="num">${fmt(x.trials ?? 0)}</td>
          <td class="dim">${esc((x.mechanism ?? '').slice(0, 52))}</td></tr>`).join('')
    : `<tr><td colspan="5" class="empty">No chemical matter. Open Targets reports
        no pocket and no ligand for this target — the honest answer for a tumour
        suppressor, and the reason the archetype is not "druggable".</td></tr>`;

  $('#d-struct tbody').innerHTML = d.structures.map((x) => `
    <tr><td><a href="https://www.rcsb.org/structure/${x.structure_id.split('_')[0]}"
              target="_blank" rel="noopener">${x.structure_id}</a></td>
      <td>${esc(x.method) || '—'}</td><td class="num">${x.resolution ?? '—'}</td>
      <td class="num">${x.coverage_start ?? '?'}–${x.coverage_end ?? '?'}</td>
      <td class="dim">${esc((x.title ?? '').slice(0, 70))}</td></tr>`).join('');
}

/* ── audit ─────────────────────────────────────────────────────────── */

function renderAudit() {
  const build = manifest.generated_at;
  const log = DEC.all();
  const stale = DEC.staleCount(build);

  $('#analyst').value = state.analyst;
  $('#a-count').textContent = `${log.length} decision${log.length === 1 ? '' : 's'}`;
  $('#audit-warn').innerHTML = stale
    ? `${stale} live decision${stale === 1 ? ' was' : 's were'} taken against an
       earlier data build than the one loaded (${build}). The evidence behind
       ${stale === 1 ? 'it' : 'them'} may have moved — re-review before relying
       on ${stale === 1 ? 'it' : 'them'}.`
    : '';

  $('#audit tbody').innerHTML = log.length ? log.map((d) => `
    <tr class="${d.superseded_by ? 'superseded' : ''}">
      <td>${d.decided_at.replace('T', ' ').slice(0, 19)}Z</td>
      <td>${esc(d.analyst)}</td>
      <td>${d.object_type}<br><span class="stamp">${esc(d.object_key)}</span></td>
      <td>${esc(d.symbol)}</td>
      <td>${esc(DEC.VARIANT_OUTCOMES[d.outcome] || DEC.TARGET_OUTCOMES[d.outcome]
                || d.outcome).split(' — ')[0]}</td>
      <td class="dim">${esc(d.prior_class ?? '—')}</td>
      <td class="wrap">${esc(d.rationale)}</td>
      <td>${d.superseded_by ? '<span class="tag">superseded</span>'
            : DEC.isStale(d, build) ? '<span class="tag stale">stale build</span>'
            : '<span class="tag">live</span>'}</td>
    </tr>`).join('')
    : `<tr><td colspan="8" class="empty">No decisions recorded. Select a
        reclassification candidate in the console and record one — the log is
        append-only, so nothing you record here is ever overwritten.</td></tr>`;
}

/* ── routing + render ──────────────────────────────────────────────── */

function writeHash() {
  const p = new URLSearchParams({ target: state.target });
  if (state.filters.search) p.set('find', state.filters.search);
  const next = `#/${state.view}?${p}`;
  if (location.hash !== next) history.replaceState(null, '', next);
}

function readHash() {
  const m = location.hash.match(/^#\/(console|dossier|audit)\??(.*)$/);
  if (!m) return;
  state.view = m[1];
  const p = new URLSearchParams(m[2]);
  const t = p.get('target');
  if (t && manifest.targets.some((x) => x.symbol === t)) state.target = t;
  const find = p.get('find');
  if (find != null) { state.filters.search = find; $('#f-search').value = find; }
}

async function render() {
  const t = manifest.targets.find((x) => x.symbol === state.target);

  document.querySelectorAll('#targets button').forEach((b) =>
    b.setAttribute('aria-pressed', String(b.dataset.target === state.target)));
  document.querySelectorAll('#views button').forEach((b) =>
    b.setAttribute('aria-pressed', String(b.dataset.view === state.view)));
  $('#view-console').hidden = state.view !== 'console';
  $('#view-dossier').hidden = state.view !== 'dossier';
  $('#view-audit').hidden = state.view !== 'audit';
  writeHash();

  const n = DEC.all().length;
  $('#audit-n').textContent = n ? `(${n})` : '';

  if (state.view === 'audit') { renderAudit(); return; }
  if (state.view === 'dossier') {
    const dossier = await load(`dossier_${t.symbol}.json`);
    renderDossier(t, dossier);
    renderTargetActions(t, dossier);
    return;
  }

  $('#c-title').textContent = `${t.symbol} · ${t.acc}`;
  $('#c-sub').textContent = t.name || '';
  const k = t.kpi;
  $('#kpi').innerHTML = [
    [fmt(k.predictions), 'predictions'], [fmt(k.actionable), 'actionable'],
    [fmt(k.high), 'high priority'], [fmt(k.upgrades), 'upgrade candidates'],
    [`${k.pct_solved}%`, 'residues solved'], [k.global_plddt, 'global pLDDT'],
  ].map(([v, kk]) => `<div><div class="v">${v}</div><div class="k">${kk}</div></div>`).join('');

  const payload = await load(`variants_${t.symbol}.json`);
  worklistAll = VM.rows(payload);
  worklistView = VM.sortWorklist(VM.applyFilters(worklistAll, state.filters));

  $('#wl-scope').textContent = state.filters.search
    ? `Worklist · search "${state.filters.search}" · filters bypassed`
    : `Worklist · ${fmt(worklistView.length)} of ${fmt(worklistAll.length)}`;
  $('#wl-cap').textContent = '';
  $('#wl-live').textContent = `${worklistView.length} variants listed`;
  $('#f-note').textContent = `transcript ${payload.transcript ?? ''}`;

  if (state.selected == null && state.filters.search && worklistView.length === 1) {
    state.selected = worklistView[0].i;
  }
  state.cursor = Math.max(0, worklistView.findIndex((r) => r.i === state.selected));
  drawWindow();
  renderDetail(state.selected != null ? worklistAll[state.selected] : null);

  const profile = await load(`profile_${t.symbol}.json`);
  $('#profile').innerHTML = profileSVG(profile);
  await renderStructure(t, profile,
    state.selected != null ? worklistAll[state.selected].pos : null);
}

/* ── wiring ────────────────────────────────────────────────────────── */

function legend(el, items) {
  $(el).innerHTML = items.map(([label, color]) =>
    `<span><i style="background:${color}"></i>${label}</span>`).join('');
}

function wireAudit() {
  $('#analyst').addEventListener('change', (e) => {
    state.analyst = e.target.value.trim();
    localStorage.setItem('locus.analyst', state.analyst);
  });
  $('#a-export').addEventListener('click', () => {
    const blob = new Blob(
      [JSON.stringify(DEC.exportPayload(manifest.generated_at), null, 1)],
      { type: 'application/json' });
    const a = document.createElement('a');
    a.href = URL.createObjectURL(blob);
    a.download = `locus-decisions-${new Date().toISOString().slice(0, 10)}.json`;
    a.click();
    URL.revokeObjectURL(a.href);
  });
  $('#a-import-file').addEventListener('change', async (e) => {
    const file = e.target.files[0];
    if (!file) return;
    try {
      const added = DEC.importPayload(JSON.parse(await file.text()));
      $('#audit-warn').textContent = `Imported ${added} new decision(s).`;
      render();
    } catch (err) { $('#audit-warn').textContent = `Import failed: ${err.message}`; }
    e.target.value = '';
  });
  $('#a-clear').addEventListener('click', () => {
    if (confirm('Clear the local decision log?\n\nThis cannot be undone. Export '
              + 'first if you want to keep it.')) { DEC.clear(); render(); }
  });
}

async function init() {
  DEC.load();
  manifest = await load('manifest.json');
  state.target = manifest.targets[0].symbol;
  $('#stamp').textContent =
    `built ${manifest.generated_at.replace('T', ' ').replace('+00:00', 'Z')}`;

  $('#targets').innerHTML = manifest.targets.map((t) =>
    `<button data-target="${t.symbol}">${t.symbol}</button>`).join('');
  $('#f-tier').innerHTML = '<option value="">any</option>'
    + VM.TIER_ORDER.map((t) => `<option value="${t}">${VM.TIER_LABEL[t]}</option>`).join('');

  legend('#legend-am', [['likely benign', PALETTE.am.LBen], ['ambiguous', PALETTE.am.Amb],
    ['likely pathogenic', PALETTE.am.LPath], ['no prediction', PALETTE.am._none]]);
  legend('#legend-tier', VM.TIER_ORDER.map((t) => [VM.TIER_LABEL[t], PALETTE.tier[t]]));
  legend('#legend-profile', [['shaded = covered by an experimental structure',
    PALETTE.tier.experimental]]);

  readHash();

  $('#targets').addEventListener('click', (e) => {
    const b = e.target.closest('button'); if (!b) return;
    state.target = b.dataset.target; state.selected = null;
    state.filters.search = ''; $('#f-search').value = '';
    render();
  });
  $('#views').addEventListener('click', (e) => {
    const b = e.target.closest('button'); if (!b) return;
    state.view = b.dataset.view; render();
  });

  const body = $('#worklist tbody');
  body.addEventListener('click', (e) => {
    const tr = e.target.closest('tr[data-i]'); if (!tr) return;
    state.cursor = Number(tr.dataset.v);
    selectRow(Number(tr.dataset.i));
  });
  body.addEventListener('keydown', (e) => {
    const keys = ['ArrowDown', 'ArrowUp', 'Home', 'End', 'PageDown', 'PageUp',
                  'Enter', ' '];
    if (!keys.includes(e.key)) return;
    e.preventDefault();
    if (e.key === 'Enter' || e.key === ' ') {
      const r = worklistView[state.cursor];
      if (r) selectRow(r.i);
      return;
    }
    const page = Math.floor(($('#wl-wrap').clientHeight || 430) / ROW_H);
    const delta = { ArrowDown: 1, ArrowUp: -1, PageDown: page, PageUp: -page }[e.key];
    moveCursor(e.key === 'Home' ? 0
      : e.key === 'End' ? worklistView.length - 1
      : state.cursor + delta);
  });

  // The rAF latch must be released in a finally: if drawWindow ever throws,
  // `ticking` would stay true and the table would silently stop following the
  // scrollbar for the rest of the session — a dead surface with no error.
  let ticking = false;
  $('#wl-wrap').addEventListener('scroll', () => {
    if (ticking) return;
    ticking = true;
    requestAnimationFrame(() => {
      try { drawWindow(); } finally { ticking = false; }
    });
  });

  const bind = (sel, ev, fn) => $(sel).addEventListener(ev, (e) => {
    fn(e); state.selected = null; render();
  });
  bind('#f-actionable', 'change', (e) => { state.filters.actionable = e.target.checked; });
  bind('#f-priority', 'change', (e) => { state.filters.priority = e.target.value; });
  bind('#f-tier', 'change', (e) => { state.filters.tier = e.target.value; });
  bind('#f-stars', 'input', (e) => {
    state.filters.minStars = Number(e.target.value);
    $('#f-stars-v').textContent = e.target.value;
  });
  let debounce;
  $('#f-search').addEventListener('input', (e) => {
    clearTimeout(debounce);
    debounce = setTimeout(() => {
      state.filters.search = e.target.value; state.selected = null; render();
    }, 200);
  });

  wireAudit();
  observeStructure();

  // Read-only introspection hook, local development only. 3Dmol does not
  // register viewers made with createViewer() anywhere reachable, so without
  // this there is no way to check from outside whether the structure→variant
  // click is actually wired — and an unverifiable feature is an unproven one.
  if (['localhost', '127.0.0.1'].includes(location.hostname)) {
    window.__locus = { viewers, state, VM, DEC, selectResidue };
  }

  await render();
}

init().catch((err) => {
  document.querySelector('main').innerHTML =
    `<div class="empty">Failed to load: ${err.message}. Serve this directory over
     HTTP — <code>uv run python -m http.server 8080 -d site</code> — since
     fetch() is blocked on file:// URLs.</div>`;
  console.error(err);
});
