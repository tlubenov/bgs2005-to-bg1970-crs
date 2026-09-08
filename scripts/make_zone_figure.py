#!/usr/bin/env python3
"""
Render the four KS1970 zone extents from bg_ext/bg_zones1970.shp as an SVG figure
for paper/recovering_ks1970.html, and print the geometry facts the paper quotes.

The zone outlines are the authoritative thing bg_ext gives us: the paper's map is
drawn from them rather than from the bounding boxes of the sample points, so what
the reader sees is the real zone layout (including its staircase map-sheet
boundaries) and not an artefact of where points happened to be sampled.

Everything is drawn in BGS2005 / CCS2005 (EPSG:7801) metres, reprojected from the
shapefile's WGS 84 / UTM zone 35N. A graticule is included because the figure has
to make one specific fact legible: the discontinuity that build_tin_grid.py finds
in AGKK's own transformation sits exactly on the 24 deg E meridian, in the two
western zones only.

The output SVG uses currentColor and CSS custom properties, so it inherits the
paper's light/dark palette instead of carrying baked-in colours.

It also renders the curvature profile that makes the second finding visible: the
second difference of AGKK's output along each lattice row, per easting column,
which is flat for K5/K7 and spikes by two orders of magnitude at 24 deg E for
K3/K9.

Usage:
    python3 make_zone_figure.py    # writes paper/figures/zones.svg and seam.svg
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
from pyproj import Transformer

from bg_points import QUADRANT, ZONES, load_zone_polygon, load_matched_points, FIT_STEP

ROOT = Path(__file__).resolve().parent.parent
FIG_DIR = ROOT / "paper" / "figures"

W, H = 900, 520
PAD = 46
SEAM_LON = 24.0
MERIDIANS = [22, 23, 24, 25, 26, 27, 28]
PARALLELS = [41.5, 42.0, 42.5, 43.0, 43.5, 44.0]

# Map fills: the paper's own muted family, one per zone. Colour is redundant here --
# every polygon carries its zone name as a label -- so these stay in the document's
# register rather than being pushed to a chart-grade chroma.
ZONE_COLOR = {"k3": "#1E6E82", "k5": "#9C6B3E", "k7": "#4B7B4A", "k9": "#7A5A8C"}

# Chart series colours. Only two zones carry a signal, so only two hues are spent;
# the flat zones are drawn as neutral reference lines. Both hues were checked with
# the palette validator against this paper's light (#EEF1EA) and dark (#111E1C)
# surfaces: lightness band, chroma floor, protan/deutan separation and contrast all
# pass in both themes.
SERIES_COLOR = {"k3": ("#0086A4", "#00A2C2"), "k9": ("#C06E00", "#D07802")}

# The paper is bilingual and the figures sit inside its articles, so they are
# generated per language -- including the zone names, which the Bulgarian text
# writes with a Cyrillic К.
STRINGS = {
    "en": {
        "zone": "K{n}",
        "quadrant": {"k3": "northwest", "k5": "southeast", "k7": "northeast", "k9": "southwest"},
        "seam_legend": "24°E: discontinuity in AGKK\u2019s transform",
        "seam_mark": "24°E",
        "x_axis": "longitude of the lattice column",
        "leg_rough": "K3, K9 — discontinuous at 24°E",
        "leg_flat": "K5, K7 — smooth throughout",
        "map_title": "The four KS1970 zones over Bulgaria",
        "map_desc": "Map of the four Coordinate System 1970 Lambert zones K3, K5, K7 and K9, "
                    "drawn from their official extents, with the 24 degrees east meridian "
                    "highlighted as the line along which the transformation is discontinuous "
                    "in the two western zones.",
        "seam_title": "Curvature of AGKK\u2019s transformation across each zone",
        "seam_desc": "Second difference of the transformed coordinates along each lattice row, "
                     "per easting column, on a logarithmic scale. Zones K5 and K7 are flat at "
                     "about half a millimetre and two and a half millimetres respectively. Zones "
                     "K3 and K9 are equally flat except for a single spike reaching about 100 and "
                     "130 millimetres exactly at the 24 degrees east meridian.",
    },
    "bg": {
        "zone": "К{n}",
        "quadrant": {"k3": "северозапад", "k5": "югоизток", "k7": "североизток", "k9": "югозапад"},
        "seam_legend": "24° и.д.: прекъсване в трансформацията на АГКК",
        "seam_mark": "24° и.д.",
        "x_axis": "дължина на колоната от решетката",
        "leg_rough": "К3, К9 — прекъснати на 24° и.д.",
        "leg_flat": "К5, К7 — гладки навсякъде",
        "map_title": "Четирите зони на КС1970 върху територията на България",
        "map_desc": "Карта на четирите конични зони К3, К5, К7 и К9 на Координатна система "
                    "1970 г., начертани по официалните им граници, с открояване на меридиана "
                    "24° и.д. — линията, по която трансформацията е прекъсната в двете западни "
                    "зони.",
        "seam_title": "Кривина на трансформацията на АГКК по зони",
        "seam_desc": "Втора разлика на преобразуваните координати по редовете на решетката, по "
                     "колони, в логаритмична скала. Зоните К5 и К7 са плоски съответно при около "
                     "половин милиметър и два и половина милиметра. К3 и К9 са също толкова "
                     "плоски, с изключение на един връх, достигащ около 100 и 130 милиметра "
                     "точно на меридиана 24° и.д.",
    },
}


def projector():
    return Transformer.from_crs("EPSG:7798", "EPSG:7801", always_xy=True)


def make_scaler(bounds):
    """Map CCS2005 (E, N) metres onto SVG pixels, preserving aspect ratio and
    flipping the northing axis (SVG y grows downward)."""
    e0, e1, n0, n1 = bounds
    sx = (W - 2 * PAD) / (e1 - e0)
    sy = (H - 2 * PAD) / (n1 - n0)
    s = min(sx, sy)
    ox = PAD + ((W - 2 * PAD) - s * (e1 - e0)) / 2
    oy = PAD + ((H - 2 * PAD) - s * (n1 - n0)) / 2

    def to_px(E, N):
        E, N = np.asarray(E, dtype=float), np.asarray(N, dtype=float)
        return ox + (E - e0) * s, H - (oy + (N - n0) * s)

    return to_px, s


def path_d(x, y):
    pts = "".join(f"{'M' if i == 0 else 'L'}{px:.1f},{py:.1f}"
                  for i, (px, py) in enumerate(zip(x, y)))
    return pts + "Z"


def build_svg(sfx: str = "", lang: str = "en") -> tuple[str, dict]:
    """`sfx` suffixes every element id. The paper shows the same figure inside both
    its Bulgarian and its English article, so each copy needs its own ids even
    though only one is visible at a time; `lang` picks the label language."""
    T = STRINGS[lang]
    rings = {z: load_zone_polygon(z) for z in ZONES}
    all_pts = np.vstack([r for z in ZONES for r in rings[z]])
    bounds = (all_pts[:, 0].min(), all_pts[:, 0].max(), all_pts[:, 1].min(), all_pts[:, 1].max())
    to_px, scale = make_scaler(bounds)
    fwd = projector()

    out = [
        f'<svg viewBox="0 0 {W} {H}" width="100%" role="img" xmlns="http://www.w3.org/2000/svg" '
        f'aria-labelledby="zonefig-title{sfx} zonefig-desc{sfx}" style="max-width:100%;height:auto;">',
        f'<title id="zonefig-title{sfx}">{T["map_title"]}</title>',
        f'<desc id="zonefig-desc{sfx}">{T["map_desc"]}</desc>',
        '<style>'
        '.grat{stroke:currentColor;stroke-opacity:.18;stroke-width:1;fill:none}'
        '.gratlbl{font:10px var(--sans,system-ui,sans-serif);fill:currentColor;'
        'fill-opacity:.45}'
        '.zlbl{font:600 15px var(--sans,system-ui,sans-serif);fill:currentColor}'
        '.zsub{font:10px var(--sans,system-ui,sans-serif);fill:currentColor;fill-opacity:.55}'
        '.seam{stroke:#C0392B;stroke-width:2;stroke-dasharray:7 4;fill:none}'
        '.seamlbl{font:600 11px var(--sans,system-ui,sans-serif);fill:#C0392B}'
        '.zone{stroke-width:1.6;fill-opacity:.14}'
        '.bar{stroke:currentColor;stroke-width:1.4;fill:none}'
        '.barlbl{font:10px var(--sans,system-ui,sans-serif);fill:currentColor;fill-opacity:.6}'
        '</style>',
    ]

    # --- graticule, clipped to the drawing area ---
    lat_span = np.linspace(41.0, 44.5, 60)
    lon_span = np.linspace(21.8, 29.0, 60)
    for lon in MERIDIANS:
        E, N = fwd.transform(np.full_like(lat_span, lon), lat_span)
        x, y = to_px(E, N)
        out.append(f'<path class="grat" d="{path_d(x, y)[:-1]}"/>')
        out.append(f'<text class="gratlbl" x="{x[-1]:.1f}" y="{PAD-24:.1f}" '
                   f'text-anchor="middle">{lon}°E</text>')
    for lat in PARALLELS:
        E, N = fwd.transform(lon_span, np.full_like(lon_span, lat))
        x, y = to_px(E, N)
        out.append(f'<path class="grat" d="{path_d(x, y)[:-1]}"/>')
        out.append(f'<text class="gratlbl" x="{PAD-38:.1f}" y="{y[0]+3:.1f}">{lat}°N</text>')

    # --- zone polygons ---
    facts = {}
    for z in ZONES:
        col = ZONE_COLOR[z]
        d = " ".join(path_d(*to_px(r[:, 0], r[:, 1])) for r in rings[z])
        out.append(f'<path class="zone" d="{d}" fill="{col}" stroke="{col}"/>')
        pts = np.vstack(rings[z])
        cx, cy = to_px(pts[:, 0].mean(), pts[:, 1].mean())
        _, src, _ = load_matched_points(z, FIT_STEP)
        facts[z] = {"n_1km": len(src)}
        out.append(f'<text class="zlbl" x="{cx:.1f}" y="{cy:.1f}" text-anchor="middle">'
                   f'{T["zone"].format(n=z[1:])}</text>')
        out.append(f'<text class="zsub" x="{cx:.1f}" y="{cy+14:.1f}" text-anchor="middle">'
                   f'{T["quadrant"][z]}</text>')

    # --- the 24 deg E seam, drawn only across the two western zones ---
    lat_seam = np.linspace(41.15, 44.3, 40)
    E, N = fwd.transform(np.full_like(lat_seam, SEAM_LON), lat_seam)
    x, y = to_px(E, N)
    out.append(f'<path class="seam" d="{path_d(x, y)[:-1]}"/>')

    # --- scale bar (bottom left) ---
    bar_m = 100_000
    bar_px = bar_m * scale
    bx, by = PAD, H - 16
    out.append(f'<path class="bar" d="M{bx},{by}h{bar_px:.1f}m0,-5v10m{-bar_px:.1f},-10v10"/>')
    out.append(f'<text class="barlbl" x="{bx + bar_px/2:.1f}" y="{by-8:.1f}" '
               f'text-anchor="middle">100 km</text>')

    # --- legend (bottom right), kept clear of the map body and the scale bar ---
    # The label is end-anchored at the right margin so it can never overflow the
    # viewBox. The dashed swatch is then placed from the label's estimated width
    # rather than at a fixed offset: the Bulgarian string is a third longer than
    # the English one, and a fixed offset puts the swatch underneath it.
    lx, ly = W - PAD, H - 20
    label_w = len(T["seam_legend"]) * 6.0
    swatch_end = lx - label_w - 10
    out.append(f'<path class="seam" d="M{swatch_end-26:.0f},{ly}h26"/>')
    out.append(f'<text class="seamlbl" x="{lx}" y="{ly+4}" text-anchor="end">'
               f'{T["seam_legend"]}</text>')

    out.append("</svg>")
    return "\n".join(out), facts



# --------------------------------------------------------------------------------
# Figure 2: curvature profile, the evidence for the 24 deg E discontinuity
# --------------------------------------------------------------------------------

CW, CH = 900, 380
CL, CR, CT, CB = 62, 26, 34, 52          # plot margins
Y_TICKS = [0.1, 1, 10, 100]              # mm, log scale
LON_LO, LON_HI = 22.2, 28.8


def curvature_profile(zone: str):
    """Median over lattice rows of |second difference| of AGKK's output, per easting
    column, in millimetres, against the column's longitude.

    A smooth projection sampled on a 1 km lattice has a second difference of a
    couple of millimetres -- that is just the conic's own curvature. Anything far
    above that is a kink in the transformation. Taking the median across every row
    of the zone (rather than reading one row) keeps a single noisy row from
    inventing or hiding a feature."""
    from collections import defaultdict

    _, src, dst = load_matched_points(zone, FIT_STEP)
    rows = defaultdict(list)
    for i in range(len(src)):
        rows[src[i, 1]].append(i)

    per_column = defaultdict(list)
    for northing, idx in rows.items():
        idx = np.array(idx)
        order = idx[np.argsort(src[idx, 0])]
        if len(order) < 5:
            continue
        X = src[order, 0]
        curv = np.hypot(np.diff(dst[order, 0], 2), np.diff(dst[order, 1], 2))
        # Only trust a triple whose two spacings are both exactly one lattice step:
        # across a gap in the row the "second difference" is meaningless.
        contiguous = (np.diff(X)[:-1] == 1000) & (np.diff(X)[1:] == 1000)
        for e, v, ok in zip(X[1:-1], curv, contiguous):
            if ok:
                per_column[e].append(v)

    eastings = np.array(sorted(per_column))
    values = np.array([np.median(per_column[e]) for e in eastings]) * 1000.0
    mid_n = np.median(src[:, 1])
    fwd = Transformer.from_crs("EPSG:7801", "EPSG:7798", always_xy=True)
    lons = np.asarray(fwd.transform(eastings, np.full_like(eastings, mid_n))[0])
    return lons, values


def build_seam_svg(sfx: str = "", lang: str = "en") -> tuple[str, dict]:
    """`sfx` suffixes every element id, `lang` picks the label language -- see build_svg()."""
    T = STRINGS[lang]
    series = {z: curvature_profile(z) for z in ZONES}

    def sx(lon):
        return CL + (np.asarray(lon) - LON_LO) / (LON_HI - LON_LO) * (CW - CL - CR)

    def sy(mm):
        t = (np.log10(np.clip(mm, Y_TICKS[0], None)) - np.log10(Y_TICKS[0])) / (
            np.log10(Y_TICKS[-1]) - np.log10(Y_TICKS[0]))
        return CH - CB - t * (CH - CT - CB)

    out = [
        f'<svg viewBox="0 0 {CW} {CH}" width="100%" role="img" '
        f'xmlns="http://www.w3.org/2000/svg" aria-labelledby="seamfig-title{sfx} seamfig-desc{sfx}" '
        f'class="seamchart" style="max-width:100%;height:auto;">',
        f'<title id="seamfig-title{sfx}">{T["seam_title"]}</title>',
        f'<desc id="seamfig-desc{sfx}">{T["seam_desc"]}</desc>',
        '<style>'
        # Series tokens are defined on the chart itself, in all three theme states,
        # so the file is correct standalone and cannot overwrite the page's tokens.
        '.seamchart{--series-k3:#0086A4;--series-k9:#C06E00}'
        '@media (prefers-color-scheme:dark){:root:not([data-theme="light"]) .seamchart'
        '{--series-k3:#00A2C2;--series-k9:#D07802}}'
        ':root[data-theme="dark"] .seamchart{--series-k3:#00A2C2;--series-k9:#D07802}'
        '.cgrid{stroke:currentColor;stroke-opacity:.14;stroke-width:1;fill:none}'
        '.caxis{stroke:currentColor;stroke-opacity:.35;stroke-width:1;fill:none}'
        '.ctick{font:11px var(--sans,system-ui,sans-serif);fill:currentColor;fill-opacity:.55;'
        'font-variant-numeric:tabular-nums}'
        '.cttl{font:600 11px var(--sans,system-ui,sans-serif);fill:currentColor;'
        'fill-opacity:.55;letter-spacing:.06em;text-transform:uppercase}'
        '.cline{fill:none;stroke-width:2;stroke-linejoin:round;stroke-linecap:round}'
        '.cflat{fill:none;stroke:currentColor;stroke-opacity:.42;stroke-width:2;'
        'stroke-dasharray:2 3}'
        '.clbl{font:600 12px var(--sans,system-ui,sans-serif)}'
        '.cflatlbl{font:11px var(--sans,system-ui,sans-serif);fill:currentColor;fill-opacity:.55}'
        '.cmark{stroke:var(--warn,#C0392B);stroke-width:1.5;stroke-dasharray:6 4;fill:none}'
        '.cmarklbl{font:600 11px var(--sans,system-ui,sans-serif);fill:var(--warn,#C0392B)}'
        '</style>',
    ]

    # gridlines + y ticks (log decades)
    for t in Y_TICKS:
        y = float(sy(t))
        out.append(f'<path class="cgrid" d="M{CL},{y:.1f}H{CW-CR}"/>')
        lab = f"{t:g}"
        out.append(f'<text class="ctick" x="{CL-10}" y="{y+4:.1f}" text-anchor="end">{lab}</text>')
    out.append(f'<text class="cttl" x="{CL-10}" y="{CT-8}" text-anchor="end">mm</text>')

    # x ticks (whole degrees of longitude)
    for lon in range(23, 29):
        x = float(sx(lon))
        out.append(f'<path class="cgrid" d="M{x:.1f},{CT}V{CH-CB}"/>')
        out.append(f'<text class="ctick" x="{x:.1f}" y="{CH-CB+18}" '
                   f'text-anchor="middle">{lon}°E</text>')
    out.append(f'<path class="caxis" d="M{CL},{CH-CB}H{CW-CR}"/>')
    out.append(f'<text class="cttl" x="{(CL+CW-CR)/2:.0f}" y="{CH-CB+40}" '
               f'text-anchor="middle">{T["x_axis"]}</text>')

    # the 24 deg E marker, behind the data
    xm = float(sx(SEAM_LON))
    out.append(f'<path class="cmark" d="M{xm:.1f},{CT}V{CH-CB}"/>')
    out.append(f'<text class="cmarklbl" x="{xm-7:.1f}" y="{CH-CB-8}" '
               f'text-anchor="end">{T["seam_mark"]}</text>')

    facts = {}
    # flat zones first, as neutral reference lines: they carry no signal, so they
    # spend no hue.
    for z in ("k5", "k7"):
        lon, mm = series[z]
        d = "".join(f"{'M' if i == 0 else 'L'}{float(sx(a)):.1f},{float(sy(b)):.1f}"
                    for i, (a, b) in enumerate(zip(lon, mm)))
        out.append(f'<path class="cflat" d="{d}"/>')
        out.append(f'<text class="cflatlbl" x="{float(sx(lon[-1]))+6:.1f}" '
                   f'y="{float(sy(mm[-1]))+4:.1f}">{T["zone"].format(n=z[1:])}</text>')
        facts[z] = {"flat_mm": float(np.median(mm)), "peak_mm": float(mm.max())}

    # K9 peaks higher than K3, so its label sits above and K3's below; without the
    # offsets the two would overlap, the peaks being within a hair of each other in
    # both longitude and log height.
    LABEL_DY = {"k9": -6, "k3": 15}
    for z in ("k3", "k9"):
        lon, mm = series[z]
        col = SERIES_COLOR[z][0]
        d = "".join(f"{'M' if i == 0 else 'L'}{float(sx(a)):.1f},{float(sy(b)):.1f}"
                    for i, (a, b) in enumerate(zip(lon, mm)))
        out.append(f'<path class="cline" style="stroke:var(--series-{z},{col})" d="{d}"/>')
        k = int(np.argmax(mm))
        out.append(f'<circle cx="{float(sx(lon[k])):.1f}" cy="{float(sy(mm[k])):.1f}" r="3.5" '
                   f'fill="var(--series-{z},{col})"/>')
        ly_lbl = max(float(sy(mm[k])) + 4 + LABEL_DY[z], CT + 10)
        out.append(f'<text class="clbl" style="fill:var(--series-{z},{col})" '
                   f'x="{float(sx(lon[k]))+9:.1f}" y="{ly_lbl:.1f}">'
                   f'{T["zone"].format(n=z[1:])} — {mm[k]:.0f} mm</text>')
        facts[z] = {"flat_mm": float(np.median(mm)), "peak_mm": float(mm.max()),
                    "peak_lon": float(lon[k])}

    # Legend. Identity is already carried by the direct labels on every line; the
    # legend states the grouping those labels do not: which zones carry the seam.
    lx, ly = CW - CR - 262, float(sy(0.25))
    out.append(f'<path class="cline" style="stroke:var(--series-k3,#0086A4)" '
               f'd="M{lx},{ly}h16"/>')
    out.append(f'<path class="cline" style="stroke:var(--series-k9,#C06E00)" '
               f'd="M{lx+20},{ly}h16"/>')
    out.append(f'<text class="cflatlbl" x="{lx+42}" y="{ly+4}">'
               f'{T["leg_rough"]}</text>')
    out.append(f'<path class="cflat" d="M{lx},{ly+18}h36"/>')
    out.append(f'<text class="cflatlbl" x="{lx+42}" y="{ly+22}">{T["leg_flat"]}</text>')

    out.append("</svg>")
    return "\n".join(x for x in out if x), facts


def inject_into_paper(paper: Path, figures: dict[str, str]) -> bool:
    """Replace the marked figure blocks in the paper with the freshly generated SVGs.

    The paper is hand-written prose but its figures are generated, so they are kept
    in sync by substitution between markers rather than by remembering to paste.
    Inlining (rather than <img src>) matters because the paper has to work as a
    single self-contained file."""
    if not paper.exists():
        return False
    html = paper.read_text()
    changed = False
    for name, svg in figures.items():
        start, end = f"<!--FIG:{name}-->", f"<!--/FIG:{name}-->"
        i, j = html.find(start), html.find(end)
        if i == -1 or j == -1:
            print(f"  note: no {start} ... {end} markers in {paper.name}, skipped")
            continue
        html = html[:i + len(start)] + "\n" + svg + "\n" + html[j:]
        changed = True
    if changed:
        paper.write_text(html)
    return changed


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, default=FIG_DIR / "zones.svg")
    ap.add_argument("--no-inject", action="store_true",
                    help="only write the .svg files, do not update paper/recovering_ks1970.html")
    args = ap.parse_args()

    svg, facts = build_svg()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(svg)
    print(f"wrote {args.out} ({len(svg)/1024:.1f} KB)")

    seam_svg, seam_facts = build_seam_svg()
    seam_path = args.out.parent / "seam.svg"
    seam_path.write_text(seam_svg)
    print(f"wrote {seam_path} ({len(seam_svg)/1024:.1f} KB)")
    print("\nCurvature of AGKK's transformation (median over lattice rows, per column):")
    for z in ZONES:
        f = seam_facts[z]
        peak = (f" -- peaks at {f['peak_mm']:.0f} mm on longitude {f['peak_lon']:.4f}"
                if "peak_lon" in f else " -- no spike anywhere in the zone")
        print(f"  K{z[1:]}: flat at {f['flat_mm']:.2f} mm{peak}")

    if not args.no_inject:
        paper = ROOT / "paper" / "recovering_ks1970.html"
        blocks = {}
        for lang in ("bg", "en"):
            blocks[f"zones-{lang}"] = build_svg(f"-{lang}", lang)[0]
            blocks[f"seam-{lang}"] = build_seam_svg(f"-{lang}", lang)[0]
        if inject_into_paper(paper, blocks):
            print(f"\ninjected both figures into {paper}, for both languages")

    print("\nZone facts drawn from bg_ext/bg_zones1970.shp:")
    from osgeo import ogr
    ds = ogr.Open(str(ROOT / "bg_ext" / "bg_zones1970.shp"))
    layer = ds.GetLayer()
    total = None
    areas = {}
    for feat in layer:
        z = f"k{feat.GetField('CLIST')}"
        g = feat.GetGeometryRef().Clone()
        areas[z] = g.GetArea() / 1e6
        total = g.Clone() if total is None else total.Union(g)
    for z in ZONES:
        print(f"  K{z[1:]} ({QUADRANT[z]:<9s}): {areas[z]:8.0f} km2, "
              f"{facts[z]['n_1km']:6d} control points on the 1 km lattice")
    print(f"  sum of zone areas {sum(areas.values()):.0f} km2 vs union "
          f"{total.GetArea()/1e6:.0f} km2 -- the zones tile the country, they do not overlap")


if __name__ == "__main__":
    main()
