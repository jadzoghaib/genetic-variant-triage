"""Visual tokens and page chrome.

The console commits deliberately to a single dark look — mission control, not a
dashboard — so only the dark column of the reference palette is used. Every
value here is a documented palette hex; none were eyeballed.

Colour assignments, each doing exactly one job:

  * AlphaMissense class  -> DIVERGING (polarity: benign <- ambiguous -> pathogenic).
    blue/red poles with a neutral gray midpoint. Validated all-pairs on the dark
    surface: CVD dE 19.2, normal-vision 29.0, both clear. Green/red status
    tokens were rejected here — a 3D structure cannot carry per-residue labels,
    so the red/green confusion would have nothing to fall back on.
  * Evidence tier -> ORDINAL (a confidence ordering, so the reader must see the
    order in the colour). One blue hue, monotone lightness. Validated with
    --ordinal on the dark surface: monotone L, adjacent dL >= 0.06, light end
    3.23:1 vs surface.
  * Triage class distribution -> NOMINAL, so every bar takes the same slot-1
    hue. Colouring nominal bars by their value would spend the identity channel
    re-encoding what bar length already shows.
"""

from __future__ import annotations

# ── surfaces and ink (dark column) ────────────────────────────────────────
PAGE = "#0d0d0d"
SURFACE = "#1a1a19"
INK = "#ffffff"
INK_SECONDARY = "#c3c2b7"
INK_MUTED = "#898781"
GRID = "#2c2c2a"
BASELINE = "#383835"
BORDER = "rgba(255,255,255,0.10)"

SERIES_1 = "#3987e5"          # categorical slot 1, dark

# ── diverging: AlphaMissense class ────────────────────────────────────────
AM_COLORS = {
    "LBen": "#3987e5",        # cool pole
    "Amb": "#383835",         # neutral midpoint
    "LPath": "#e66767",       # warm pole
}
NO_DATA = "#2c2c2a"

# ── ordinal: structural evidence tier (brightest = strongest evidence) ────
TIER_COLORS = {
    "experimental": "#cde2fb",
    "predicted_confident": "#6da7ec",
    "predicted_weak": "#256abf",
}

TIER_LABELS = {
    "experimental": "experimentally solved",
    "predicted_confident": "confidently predicted",
    "predicted_weak": "predicted, unreliable",
}

# ── status (fixed scale, never themed) ────────────────────────────────────
STATUS = {"good": "#0ca30c", "warning": "#fab219",
          "serious": "#ec835a", "critical": "#d03b3b"}

PRIORITY_STATUS = {"high": "critical", "medium": "serious", "low": "warning"}

FONT = 'system-ui, -apple-system, "Segoe UI", sans-serif'


CSS = f"""
<style>
  .stApp {{ background: {PAGE}; }}
  html, body, [class*="css"] {{ font-family: {FONT}; }}

  /* Dense but calm: tighten Streamlit's generous default rhythm. */
  .block-container {{ padding-top: 2.2rem; padding-bottom: 2rem; max-width: 100%; }}
  section[data-testid="stSidebar"] {{ background: {SURFACE};
                                      border-right: 1px solid {BORDER}; }}

  h1, h2, h3 {{ font-weight: 600; letter-spacing: -0.01em; color: {INK}; }}
  h1 {{ font-size: 1.35rem !important; margin-bottom: 0.1rem; }}

  .locus-eyebrow {{ font-size: 0.7rem; text-transform: uppercase;
                    letter-spacing: 0.09em; color: {INK_MUTED}; }}
  .locus-sub {{ color: {INK_SECONDARY}; font-size: 0.85rem; }}

  .locus-panel {{ background: {SURFACE}; border: 1px solid {BORDER};
                  border-radius: 6px; padding: 0.85rem 1rem; }}

  /* KPI strip — figures use tabular numerals only where they stack. */
  .locus-kpi {{ display: flex; gap: 0; border: 1px solid {BORDER};
                border-radius: 6px; overflow: hidden; background: {SURFACE}; }}
  .locus-kpi > div {{ flex: 1; padding: 0.7rem 1rem;
                      border-right: 1px solid {BORDER}; }}
  .locus-kpi > div:last-child {{ border-right: none; }}
  .locus-kpi .v {{ font-size: 1.5rem; font-weight: 600; color: {INK};
                   line-height: 1.15; }}
  .locus-kpi .k {{ font-size: 0.68rem; text-transform: uppercase;
                   letter-spacing: 0.08em; color: {INK_MUTED}; }}

  .locus-chip {{ display: inline-block; padding: 0.1rem 0.5rem;
                 border-radius: 3px; font-size: 0.72rem; font-weight: 600;
                 border: 1px solid {BORDER}; }}

  .locus-kv {{ display: grid; grid-template-columns: 8.5rem 1fr; gap: 0.25rem 0.6rem;
               font-size: 0.82rem; }}
  .locus-kv .k {{ color: {INK_MUTED}; }}
  .locus-kv .v {{ color: {INK}; font-variant-numeric: tabular-nums; }}

  .locus-legend {{ display: flex; gap: 0.9rem; flex-wrap: wrap;
                   font-size: 0.72rem; color: {INK_SECONDARY};
                   margin: 0.15rem 0 0.35rem 0; }}
  .locus-legend span.sw {{ display: inline-block; width: 9px; height: 9px;
                           border-radius: 2px; margin-right: 0.35rem; }}

  div[data-testid="stDataFrame"] {{ border: 1px solid {BORDER}; border-radius: 6px; }}
  hr {{ border-color: {BORDER}; }}
</style>
"""


def chip(label: str, color: str) -> str:
    """A colour swatch is never the only carrier of identity — the label rides
    with it."""
    return (f'<span class="locus-chip" style="color:{color};'
            f'border-color:{color}55">{label}</span>')


def legend(items: list[tuple[str, str]]) -> str:
    parts = "".join(
        f'<span><span class="sw" style="background:{c}"></span>{lab}</span>'
        for lab, c in items)
    return f'<div class="locus-legend">{parts}</div>'
