"""3D structure panels.

Residues are grouped into colour buckets and issued as one `setStyle` call per
bucket rather than one per residue — a 1,863-residue protein would otherwise
mean 1,863 calls per render.

The viewer is fed a plain per-residue colour mapping, which is the same data
contract Mol* would take. Swapping viewers later touches only this module.
"""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path

import httpx
import py3Dmol

from ui import theme

CACHE = Path(__file__).parent.parent / "data"


#: Cartoon rendering needs only the backbone trace. Side-chain atoms are ~55%
#: of an AlphaFold PDB and contribute nothing to a cartoon, so they are dropped
#: before the text is embedded in the page. EGFR's file falls from ~800 KB to
#: ~360 KB, and two viewers re-embedding it on every rerun is the difference
#: between a responsive page and an unresponsive renderer.
BACKBONE = {"N", "CA", "C", "O", "CB"}


def fetch_structure(url: str, structure_id: str) -> str:
    """AlphaFold PDB text, cached on disk so the console runs offline."""
    CACHE.mkdir(exist_ok=True)
    path = CACHE / f"{structure_id}.pdb"
    if not path.exists():
        r = httpx.get(url, timeout=180)
        r.raise_for_status()
        path.write_bytes(r.content)
    return path.read_text(encoding="utf-8", errors="replace")


def backbone_only(pdb_text: str) -> str:
    """Strip side-chain atoms, keeping everything a cartoon needs."""
    keep = []
    for line in pdb_text.splitlines():
        if line.startswith(("ATOM", "HETATM")):
            if line[12:16].strip() in BACKBONE:
                keep.append(line)
        elif not line.startswith("ANISOU"):
            keep.append(line)
    return "\n".join(keep)


def am_class_colors(profile) -> dict[int, str]:
    """Colour each residue by the AlphaMissense class of its most damaging
    substitution.

    Buckets use AlphaMissense's own published thresholds (0.34 / 0.564) rather
    than arbitrary cuts, so the colour boundaries mean what the model means.
    """
    colors: dict[int, str] = {}
    for pos, max_am in zip(profile["position"], profile["max_am"]):
        if max_am is None or max_am != max_am:
            colors[int(pos)] = theme.NO_DATA
        elif max_am >= 0.564:
            colors[int(pos)] = theme.AM_COLORS["LPath"]
        elif max_am <= 0.34:
            colors[int(pos)] = theme.AM_COLORS["LBen"]
        else:
            colors[int(pos)] = theme.AM_COLORS["Amb"]
    return colors


def evidence_tier_colors(profile) -> dict[int, str]:
    from core.confidence import evidence_tier

    return {
        int(pos): theme.TIER_COLORS[evidence_tier(bool(solved), float(plddt))]
        for pos, solved, plddt in zip(
            profile["position"], profile["is_solved"], profile["plddt"])
    }


def render(pdb_text: str, colors: dict[int, str], selected: int | None = None,
           height: int = 430) -> str:
    """Standalone HTML for one structure panel."""
    view = py3Dmol.view(width="100%", height=height)
    view.addModel(pdb_text, "pdb")
    view.setStyle({}, {"cartoon": {"color": theme.NO_DATA}})

    buckets: dict[str, list[int]] = defaultdict(list)
    for pos, color in colors.items():
        buckets[color].append(pos)
    for color, positions in buckets.items():
        view.setStyle({"resi": positions}, {"cartoon": {"color": color}})

    if selected is not None:
        # A 2px surface ring equivalent: the selection reads as a distinct mark,
        # not merely a different fill.
        view.addStyle({"resi": [int(selected)]},
                      {"stick": {"radius": 0.32, "color": theme.INK}})
        view.addStyle({"resi": [int(selected)]},
                      {"sphere": {"radius": 0.9, "color": theme.INK, "opacity": 0.55}})
        # Frame a window around the residue rather than the residue alone.
        # zoomTo() on a single position fits the viewport to one side chain,
        # which loses every bit of structural context that makes the view
        # worth showing.
        lo, hi = max(1, int(selected) - 30), int(selected) + 30
        view.zoomTo({"resi": list(range(lo, hi + 1))})
    else:
        view.zoomTo()

    view.setBackgroundColor(theme.SURFACE)
    return view._make_html()


def tier_runs(profile) -> list[dict]:
    """Collapse a per-residue tier column into contiguous spans.

    A rect mark needs a width. With a continuous x and no x2 every mark
    collapses to zero pixels and the band silently disappears, so the strip is
    drawn from explicit start/end spans instead of 1,863 point marks — a
    handful of rects rather than one per residue, and guaranteed to render.
    """
    from core.confidence import evidence_tier

    tiers = [evidence_tier(bool(s), float(p))
             for s, p in zip(profile["is_solved"], profile["plddt"])]
    positions = [int(p) for p in profile["position"]]
    if not tiers:
        return []

    runs, start, current = [], positions[0], tiers[0]
    for pos, tier in zip(positions[1:], tiers[1:]):
        if tier != current:
            runs.append({"start": start, "end": pos, "tier": current})
            start, current = pos, tier
    runs.append({"start": start, "end": positions[-1] + 1, "tier": current})
    return runs


AM_LEGEND = [
    ("likely benign", theme.AM_COLORS["LBen"]),
    ("ambiguous", theme.AM_COLORS["Amb"]),
    ("likely pathogenic", theme.AM_COLORS["LPath"]),
    ("no prediction", theme.NO_DATA),
]

TIER_LEGEND = [(theme.TIER_LABELS[k], v) for k, v in theme.TIER_COLORS.items()]
