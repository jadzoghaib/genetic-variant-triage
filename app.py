"""Locus — variant triage console.

A rendering surface only. Every rule it displays comes from `core/`, every fact
from `queries.py`; this module computes nothing of its own.

Layout follows the traversal rule: from any variant it is at most two clicks to
its residue, its structure, its target, and the retrieval that produced it.
"""

from __future__ import annotations

import altair as alt
import pandas as pd
import streamlit as st

import queries
from core import confidence as C
from core import triage as T
from ui import structure as S
from ui import theme

st.set_page_config(page_title="Locus — variant triage console",
                   layout="wide", initial_sidebar_state="expanded")
st.markdown(theme.CSS, unsafe_allow_html=True)


# ── data access (cached) ──────────────────────────────────────────────────

@st.cache_resource
def _con():
    return queries.connect(read_only=True)


@st.cache_data(show_spinner=False)
def load_targets() -> pd.DataFrame:
    return queries.targets(_con())


@st.cache_data(show_spinner=False)
def load_worklist(symbol: str) -> pd.DataFrame:
    return T.assign(C.assign(queries.worklist(_con(), symbol)))


@st.cache_data(show_spinner=False)
def load_profile(acc: str) -> pd.DataFrame:
    return queries.residue_profile(_con(), acc)


@st.cache_data(show_spinner=False)
def load_structure(acc: str):
    meta = queries.predicted_structure(_con(), acc)
    if not meta:
        return None, None
    return meta, S.backbone_only(
        S.fetch_structure(meta["file_url"], meta["structure_id"]))


@st.cache_data(show_spinner=False)
def structure_html(acc: str, mode: str, focus: int | None) -> str:
    """Viewer HTML, cached on the only three things that change it.

    Without this every rerun rebuilds both viewers and hands the browser two
    fresh WebGL contexts, which it does not reliably release.
    """
    meta, pdb = load_structure(acc)
    profile = load_profile(acc)
    colors = (S.am_class_colors(profile) if mode == "am"
              else S.evidence_tier_colors(profile))
    return S.render(pdb, colors, focus)


@st.cache_data(show_spinner=False)
def load_evidence(variant_id: str) -> pd.DataFrame:
    return queries.variant_evidence(_con(), variant_id)


# ── sidebar ───────────────────────────────────────────────────────────────

targets = load_targets()
symbols = list(targets["symbol"])

# Console state lives in the URL, so any view can be linked to. A traversal
# surface whose state cannot be handed to someone else is only half a surface.
qp = st.query_params
qp_target, qp_find = qp.get("target"), qp.get("find", "")

with st.sidebar:
    st.markdown('<div class="locus-eyebrow">Target</div>', unsafe_allow_html=True)
    symbol = st.selectbox("Target", symbols, label_visibility="collapsed",
                          index=symbols.index(qp_target) if qp_target in symbols else 0)
    acc = targets.set_index("symbol").loc[symbol, "uniprot_acc"]
    st.caption(targets.set_index("symbol").loc[symbol, "approved_name"])

    st.divider()
    st.markdown('<div class="locus-eyebrow">Filters</div>', unsafe_allow_html=True)
    only_actionable = st.toggle("Actionable only", value=True,
                                help="Reclassification candidates and cases where "
                                     "the model contradicts a curator.")
    priorities = st.multiselect("Review priority", ["high", "medium", "low", "none"],
                                default=["high", "medium"])
    tiers = st.multiselect("Structural evidence", list(C.TIER_ORDER),
                           default=list(C.TIER_ORDER),
                           format_func=lambda t: theme.TIER_LABELS[t])
    min_stars = st.slider("Minimum ClinVar review status (stars)", 0, 4, 0)
    search = st.text_input("Find protein variant", value=qp_find,
                           placeholder="e.g. L858R")
    st.divider()
    show_3d = st.toggle("Render structure", value=True,
                        help="Disable for faster filtering on large proteins.")

# Write state back only when it actually changed, or each write triggers a
# rerun that triggers another write.
if qp.get("target") != symbol:
    qp["target"] = symbol
if (qp.get("find") or "") != search:
    if search:
        qp["find"] = search
    else:
        qp.pop("find", None)


# ── load + filter ─────────────────────────────────────────────────────────

df = load_worklist(symbol)
profile = load_profile(acc)

if search:
    # A name lookup is not a filter refinement. Someone typing "L858R" wants
    # that variant, and the filters that define a triage worklist would hide it:
    # ClinVar classifies L858R as drug_response, so it is `not_triaged` and
    # carries no priority. Searching therefore escapes the filters entirely.
    view = df[df["protein_variant"].str.contains(search.strip(), case=False,
                                                 na=False)]
else:
    view = df
    if only_actionable:
        view = view[view["triage_class"].isin(T.ACTIONABLE)]
    if priorities:
        view = view[view["review_priority"].isin(priorities)]
    if tiers:
        view = view[view["evidence_tier"].isin(tiers)]
    if min_stars:
        view = view[view["stars"].fillna(0) >= min_stars]

view = view.sort_values(["review_priority", "am_pathogenicity"],
                        ascending=[True, False])


# ── header + KPI strip ────────────────────────────────────────────────────

actionable = df[df["triage_class"].isin(T.ACTIONABLE)]
pct_solved = 100.0 * profile["is_solved"].mean()

st.markdown(f'<div class="locus-eyebrow">Variant triage console</div>'
            f'<h1>{symbol} · {acc}</h1>', unsafe_allow_html=True)

st.markdown(f"""
<div class="locus-kpi">
  <div><div class="v">{len(df):,}</div><div class="k">predictions</div></div>
  <div><div class="v">{len(actionable):,}</div><div class="k">actionable</div></div>
  <div><div class="v">{(actionable.review_priority == 'high').sum():,}</div>
       <div class="k">high priority</div></div>
  <div><div class="v">{(df.triage_class == T.UPGRADE).sum():,}</div>
       <div class="k">upgrade candidates</div></div>
  <div><div class="v">{pct_solved:.0f}%</div><div class="k">residues solved</div></div>
</div>
""", unsafe_allow_html=True)

st.write("")

# ── split view: worklist | detail ─────────────────────────────────────────

left, right = st.columns([0.63, 0.37], gap="medium")

with left:
    scope = (f'search "{search}" · filters bypassed' if search
             else f'{len(view):,} of {len(df):,}')
    st.markdown(f'<div class="locus-eyebrow">Worklist · {scope}</div>',
                unsafe_allow_html=True)
    table = view[["protein_variant", "triage_class", "review_priority",
                  "significance", "stars", "am_pathogenicity", "am_class",
                  "plddt", "evidence_tier", "variant_id"]]
    event = st.dataframe(
        table, height=430, hide_index=True, on_select="rerun",
        selection_mode="single-row",
        column_config={
            "protein_variant": st.column_config.TextColumn("Variant", width="small"),
            "triage_class": st.column_config.TextColumn("Triage", width="medium"),
            "review_priority": st.column_config.TextColumn("Priority", width="small"),
            "significance": st.column_config.TextColumn("ClinVar", width="small"),
            "stars": st.column_config.NumberColumn("★", width="small", format="%d"),
            "am_pathogenicity": st.column_config.ProgressColumn(
                "AlphaMissense", min_value=0.0, max_value=1.0, format="%.3f"),
            "am_class": st.column_config.TextColumn("Class", width="small"),
            "plddt": st.column_config.NumberColumn("pLDDT", format="%.1f",
                                                   width="small"),
            "evidence_tier": st.column_config.TextColumn("Structure", width="medium"),
            "variant_id": st.column_config.TextColumn("Genomic", width="medium"),
        },
    )

rows = event.selection.rows if event and event.selection else []
if rows:
    selected = view.iloc[rows[0]]
elif search and len(view) == 1:
    # A search that resolves to exactly one variant is already a selection;
    # making the reader click the row they just named is friction.
    selected = view.iloc[0]
else:
    selected = None

with right:
    st.markdown('<div class="locus-eyebrow">Detail</div>', unsafe_allow_html=True)
    if selected is None:
        st.markdown('<div class="locus-panel locus-sub">Select a variant to see its '
                    'residue, its structural evidence, and the retrieval that '
                    'produced every field.</div>', unsafe_allow_html=True)
    else:
        tier = selected["evidence_tier"]
        st.markdown(
            f'<div class="locus-panel">'
            f'<div style="font-size:1.1rem;font-weight:600;color:{theme.INK}">'
            f'{selected["protein_variant"]}</div>'
            f'<div class="locus-sub" style="margin-bottom:.6rem">'
            f'{selected["variant_id"]} · {selected["transcript_id"]}</div>'
            f'{theme.chip(selected["triage_class"], theme.STATUS[theme.PRIORITY_STATUS.get(selected["review_priority"], "warning")] if selected["review_priority"] != "none" else theme.INK_MUTED)}'
            f' {theme.chip(theme.TIER_LABELS[tier], theme.TIER_COLORS[tier])}'
            f'<div class="locus-kv" style="margin-top:.7rem">'
            f'<div class="k">ClinVar</div><div class="v">{selected["significance"] or "no assertion"}'
            f' ({0 if pd.isna(selected["stars"]) else int(selected["stars"])}★)</div>'
            f'<div class="k">raw label</div><div class="v">{selected["significance_raw"] or "—"}</div>'
            f'<div class="k">AlphaMissense</div><div class="v">{selected["am_pathogenicity"]:.4f} · {selected["am_class"]}</div>'
            f'<div class="k">model strength</div><div class="v">{selected["model_strength"]}</div>'
            f'<div class="k">residue</div><div class="v">{selected["ref_aa"]}{int(selected["aa_pos"])} · pLDDT {selected["plddt"]:.1f}</div>'
            f'<div class="k">solved</div><div class="v">{"yes" if selected["is_solved"] else "no"}</div>'
            f'<div class="k">priority</div><div class="v">{selected["review_priority"]}</div>'
            f'</div>'
            f'<div class="locus-sub" style="margin-top:.6rem;font-size:.75rem">'
            f'{selected["priority_reasons"] or "—"}</div>'
            f'</div>', unsafe_allow_html=True)

        with st.expander("Evidence trail", expanded=False):
            ev = load_evidence(selected["variant_id"])
            st.dataframe(ev, hide_index=True, height=140)

st.write("")

# ── structure ─────────────────────────────────────────────────────────────

meta, pdb_text = load_structure(acc)
if show_3d and pdb_text:
    focus = int(selected["aa_pos"]) if selected is not None else None
    st.markdown(f'<div class="locus-eyebrow">Structure · {meta["structure_id"]}'
                f'{" · focused on residue " + str(focus) if focus else ""}</div>',
                unsafe_allow_html=True)
    sc1, sc2 = st.columns(2, gap="medium")
    with sc1:
        st.markdown('<div class="locus-sub">Predicted effect — most damaging '
                    'substitution per residue</div>', unsafe_allow_html=True)
        st.markdown(theme.legend(S.AM_LEGEND), unsafe_allow_html=True)
        st.components.v1.html(structure_html(acc, "am", focus), height=440)
    with sc2:
        st.markdown('<div class="locus-sub">Structural evidence — how much weight '
                    'that reading can bear</div>', unsafe_allow_html=True)
        st.markdown(theme.legend(S.TIER_LEGEND), unsafe_allow_html=True)
        st.components.v1.html(structure_html(acc, "tier", focus), height=440)

st.write("")

# ── residue profile ───────────────────────────────────────────────────────
# Two measures on different scales, so two charts sharing an x-axis — never a
# second y-axis on one chart.

st.markdown('<div class="locus-eyebrow">Sequence profile</div>',
            unsafe_allow_html=True)

prof = profile.copy()
prof["tier"] = [C.evidence_tier(bool(s), float(p))
                for s, p in zip(prof["is_solved"], prof["plddt"])]

AXIS = dict(labelColor=theme.INK_MUTED, titleColor=theme.INK_MUTED,
            gridColor=theme.GRID, domainColor=theme.BASELINE,
            tickColor=theme.BASELINE, labelFontSize=10, titleFontSize=10)


def _style(ch):
    return (ch.configure_view(strokeWidth=0)
              .configure_axis(**AXIS)
              .configure_legend(labelColor=theme.INK_SECONDARY,
                                titleColor=theme.INK_MUTED, labelFontSize=10,
                                titleFontSize=10))


XSCALE = alt.Scale(nice=False, domainMin=1, domainMax=int(prof["position"].max()))
x_bare = alt.X("position:Q", title=None, scale=XSCALE,
               axis=alt.Axis(labels=False, ticks=False, domain=False))
x_axis = alt.X("position:Q", title="residue position", scale=XSCALE)

tips = [alt.Tooltip("position:Q", title="residue"),
        alt.Tooltip("wt_aa:N", title="wild type"),
        alt.Tooltip("max_am:Q", title="max AlphaMissense", format=".3f"),
        alt.Tooltip("plddt:Q", title="pLDDT", format=".1f"),
        alt.Tooltip("tier:N", title="evidence")]


# A protein is longer than the plot is wide — 1,863 residues over ~1,000px is
# over-plotting, and a per-residue area renders as a hairball. Binning to a
# window near one bin per pixel keeps the regional signal (which domains carry
# variant burden) that the profile is actually for.
BINS = 180


def _area(x_enc, field: str, agg: str, title: str, domain: list, rules: list[float]):
    """One measure over sequence position, with its own y-scale.

    Two measures on different scales get two charts sharing an x-axis — never a
    second y-axis on one chart. Properties belong on the layer, not on a layer
    member, or the sub-chart's sizing is silently dropped.
    """
    area = (alt.Chart(prof)
            .mark_area(color=theme.SERIES_1, opacity=0.9,
                       line={"color": theme.SERIES_1, "size": 1})
            .encode(x=x_enc.copy().bin(maxbins=BINS),
                    y=alt.Y(f"{field}:Q", title=title, aggregate=agg,
                            scale=alt.Scale(domain=domain, nice=False)),
                    tooltip=[alt.Tooltip("position:Q", title="residue", bin=True),
                             alt.Tooltip(f"{field}:Q", title=title, aggregate=agg,
                                         format=".2f")]))
    rule = (alt.Chart(pd.DataFrame({"y": rules}))
            .mark_rule(color=theme.INK_MUTED, strokeDash=[3, 3], size=1, opacity=0.8)
            .encode(y=alt.Y("y:Q", scale=alt.Scale(domain=domain, nice=False))))
    return alt.layer(area, rule).properties(height=88, width="container")


# Reference lines are AlphaMissense's own published cuts and AlphaFold's
# confident band — lines the reader can name, not decorative gridlines.
am_chart = _area(x_bare, "max_am", "mean", "variant burden", [0, 1], [0.34, 0.564])

# Experimentally solved spans are shaded behind the confidence curve rather
# than shown as a separate band: coverage is precisely what says how much of
# that curve you need to rely on, so the two belong in one frame. Explicit
# y/y2 give the rect its height — a rect with no y encoding has none, which is
# why a standalone strip renders blank.
solved = pd.DataFrame([r for r in S.tier_runs(prof) if r["tier"] == "experimental"])
plddt_layers = []
if len(solved):
    solved["y0"], solved["y1"] = 0, 100
    plddt_layers.append(
        alt.Chart(solved)
        .mark_rect(color=theme.TIER_COLORS["experimental"], opacity=0.16)
        # Both layers must declare the same x axis: Vega-Lite merges axes across
        # a layer, and an `axis=None` here silently removes the shared x-axis
        # labels from the chart underneath.
        .encode(x=alt.X("start:Q", scale=XSCALE, title="residue position"),
                x2="end:Q",
                y=alt.Y("y0:Q", scale=alt.Scale(domain=[0, 100], nice=False),
                        title=None),
                y2="y1:Q",
                tooltip=[alt.Tooltip("start:Q", title="solved from"),
                         alt.Tooltip("end:Q", title="to")]))

# Taller than the burden chart because the x-axis labels and title come out of
# the same box; equal `height` would squash the plotting area.
plddt_chart = alt.layer(
    *plddt_layers, _area(x_axis, "plddt", "mean", "pLDDT", [0, 100], [70])
).properties(height=120, width="container")

st.altair_chart(_style(am_chart), use_container_width=True)
st.altair_chart(_style(plddt_chart), use_container_width=True)
# Shading is never the sole carrier of meaning — the band is named in text.
st.markdown(
    f'<div class="locus-legend"><span>'
    f'<span class="sw" style="background:{theme.TIER_COLORS["experimental"]};'
    f'opacity:.45"></span>shaded = covered by an experimental structure; '
    f'dashed lines mark AlphaMissense\'s published thresholds (0.34 / 0.564) '
    f'and AlphaFold\'s confident band (70)</span></div>',
    unsafe_allow_html=True)

st.caption("AlphaMissense data is CC BY-NC-SA 4.0 — non-commercial use only. "
           "Locus surfaces candidates for expert review; it does not assert "
           "clinical classifications.")
