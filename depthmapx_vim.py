#!/usr/bin/env python3
"""depthmapx-vim: Visibility Integration Mismatch (VIM) for depthmapX VGA exports.

VIM localises intelligibility. Intelligibility is one correlation for a whole plan
(connectivity vs. visual integration); VIM asks, point by point, whether a location is
more or less integrated than other locations with a similar visual field.

    1. Valid points (connectivity and integration both >= 0) are ranked by connectivity
       and split into equal-sized quantile bands (40 by default). Points with identical
       connectivity always share a band: the tie group goes to the band of its middle rank.
       On small plans the number of bands drops so each band has at least ~10 points.
    2. Within each band b: median(Int)_b and sample standard deviation SD(Int)_b.
    3. VIM_i = (Int_i - median(Int)_b) / SD(Int)_b      (VIM = 0 when SD = 0)

    VIM >= 2       strong hidden centre (more integrated than its view suggests)
    1 <= VIM < 2   moderate hidden centre
    VIM <= -1      false centre (looks open but is configurationally peripheral)

    Strong hidden-centre points are grouped into clusters of touching grid cells
    (8-neighbour) and clusters are ranked by summed VIM (combines strength and size).

Input: a point-map export from any depthmapX version (File > Export, .csv or .txt) after
visibility analysis with global measures. Needs Connectivity, Visual Integration [HH]
and either Ref or x/y. Pure standard library; Python 3.8+.

This is a line-by-line port of the C++ implementation (salalib/vgamodules/vgavim.cpp);
on the same data it gives the same values to depthmapX's export precision.

Usage:
    python depthmapx_vim.py plan_vga.csv
    python depthmapx_vim.py plan_vga.txt -o results.csv --bands 40
    python depthmapx_vim.py --help

Copyright (C) 2026 Jennat Jennifer Jogy and Govind Dev.
Licensed under the GNU General Public License v3.0 or later; see LICENSE.
"""

import argparse
import csv
import math
import os
import struct
import sys

__version__ = "1.0.0"

DEFAULT_BANDS = 40
DEFAULT_MIN_PER_BAND = 10
CONNECTIVITY = "Connectivity"
INTEGRATION_PREFIX = "Visual Integration [HH]"
MISSING = -1.0  # depthmapX's "no value"

CLASS_NAMES = {
    0: "false centre",
    1: "as expected",
    2: "moderate hidden centre",
    3: "strong hidden centre",
}


# --------------------------------------------------------------------------- helpers

def to_float32(value):
    """Round a Python float to the nearest float32, as depthmapX stores attribute values."""
    return struct.unpack("f", struct.pack("f", value))[0]


def fmt(value):
    """Format a number the way depthmapX writes exports (8 significant digits)."""
    text = "%.8g" % value
    return "0" if text == "-0" else text


def median(values):
    v = sorted(values)
    n = len(v)
    if n == 0:
        return 0.0
    return v[n // 2] if n % 2 else 0.5 * (v[n // 2 - 1] + v[n // 2])


def sample_sd(values):
    # same summation order as the C++ version, so results agree to the last bit
    if len(values) < 2:
        return 0.0
    mean = 0.0
    for x in values:
        mean += x
    mean /= float(len(values))
    ss = 0.0
    for x in values:
        ss += (x - mean) * (x - mean)
    return math.sqrt(ss / float(len(values) - 1))


def pearson(a, b):
    n = len(a)
    if n < 2:
        return float("nan")
    ma, mb = sum(a) / n, sum(b) / n
    cov = sum((x - ma) * (y - mb) for x, y in zip(a, b))
    va = sum((x - ma) ** 2 for x in a)
    vb = sum((y - mb) ** 2 for y in b)
    return cov / math.sqrt(va * vb) if va > 0 and vb > 0 else float("nan")


# --------------------------------------------------------------------------- core

class VimPoint:
    __slots__ = ("index", "gx", "gy", "x", "y", "connectivity", "integration",
                 "band", "vim", "cls", "cluster", "cluster_size")

    def __init__(self, index, gx, gy, connectivity, integration, x=None, y=None):
        self.index = index            # row number in the input file
        self.gx, self.gy = gx, gy     # integer grid cell
        self.x = float(gx) if x is None else x   # plan coordinates, for reports
        self.y = float(gy) if y is None else y
        self.connectivity = connectivity
        self.integration = integration
        self.band = 0
        self.vim = 0.0
        self.cls = 1
        self.cluster = 0
        self.cluster_size = 0


def compute_vim(points, num_bands=DEFAULT_BANDS, min_per_band=DEFAULT_MIN_PER_BAND,
                strong=2.0, moderate=1.0, false_centre=-1.0):
    """Compute VIM in place on a list of VimPoint (valid points only).

    Points must be ordered by grid column then row (gx, gy), which is the order
    depthmapX visits them; this only affects floating-point summation order.
    Returns the list of clusters as (rank, member_points, vim_sum), rank 1 first.
    """
    n = len(points)
    if n == 0:
        return []
    num_bands = max(1, num_bands)
    bands = max(1, min(num_bands, n // max(1, min_per_band)))

    # 1. rank by connectivity; ties share the band of the tie group's middle rank
    order = sorted(range(n), key=lambda i: points[i].connectivity)  # stable, like std::stable_sort
    k = 0
    while k < n:
        e = k
        while e + 1 < n and points[order[e + 1]].connectivity == points[order[k]].connectivity:
            e += 1
        mid_rank = 0.5 * float(k + e)
        band = int(math.floor(mid_rank * bands / float(n)))
        band = min(band, bands - 1)
        for t in range(k, e + 1):
            points[order[t]].band = band
        k = e + 1

    # 2. median and sample SD of integration per band
    band_values = [[] for _ in range(bands)]
    for p in points:
        band_values[p.band].append(p.integration)
    band_median = [median(v) for v in band_values]
    band_sd = [sample_sd(v) for v in band_values]

    # 3. VIM and class
    index = {}
    for q, p in enumerate(points):
        sd = band_sd[p.band]
        p.vim = (p.integration - band_median[p.band]) / sd if sd > 0.0 else 0.0
        if p.vim >= strong:
            p.cls = 3
        elif p.vim >= moderate:
            p.cls = 2
        elif p.vim <= false_centre:
            p.cls = 0
        else:
            p.cls = 1
        index[(p.gx, p.gy)] = q

    # 4. cluster strong hidden centres (8-neighbour flood fill)
    clusters = []
    visited = [False] * n
    for q in range(n):
        if points[q].cls != 3 or visited[q]:
            continue
        members, vim_sum = [], 0.0
        stack = [q]
        visited[q] = True
        while stack:
            cur = stack.pop()
            members.append(cur)
            vim_sum += points[cur].vim
            cx, cy = points[cur].gx, points[cur].gy
            for dx in (-1, 0, 1):
                for dy in (-1, 0, 1):
                    if dx == 0 and dy == 0:
                        continue
                    nb = index.get((cx + dx, cy + dy))
                    if nb is not None and not visited[nb] and points[nb].cls == 3:
                        visited[nb] = True
                        stack.append(nb)
        clusters.append((members, vim_sum))
    clusters.sort(key=lambda c: c[1], reverse=True)  # rank 1 = largest summed VIM

    ranked = []
    for r, (members, vim_sum) in enumerate(clusters, start=1):
        for m in members:
            points[m].cluster = r
            points[m].cluster_size = len(members)
        ranked.append((r, [points[m] for m in members], vim_sum))
    return ranked


# --------------------------------------------------------------------------- I/O

class InputError(Exception):
    pass


def read_table(path):
    """Read a depthmapX export (comma or tab separated). Returns (header, rows, delimiter)."""
    for encoding in ("utf-8-sig", "latin-1"):
        try:
            with open(path, newline="", encoding=encoding) as fh:
                text = fh.read()
            break
        except UnicodeDecodeError:
            continue
    first = text.split("\n", 1)[0]
    delimiter = "\t" if first.count("\t") > first.count(",") else ","
    reader = csv.reader(text.splitlines(), delimiter=delimiter)
    rows = [r for r in reader if r]
    if not rows:
        raise InputError("the file is empty")
    header = [h.strip() for h in rows[0]]
    return header, rows[1:], delimiter


def parse_number(text):
    try:
        v = float(text)
    except (TypeError, ValueError):
        return None
    return v if math.isfinite(v) else None


def grid_cells(header, rows):
    """Integer grid cell per row: decoded from depthmapX's Ref, else derived from x/y."""
    if "Ref" in header:
        col = header.index("Ref")
        cells = []
        for r in rows:
            ref = parse_number(r[col])
            if ref is None:
                break
            ref = int(ref)
            cells.append((ref >> 16, ref & 0xFFFF))  # PixelRef: (x << 16) + y
        else:
            return cells, "Ref"
    if "x" not in header or "y" not in header:
        raise InputError("needs a 'Ref' column or 'x' and 'y' columns to locate points on the grid")
    xi, yi = header.index("x"), header.index("y")
    xs = [parse_number(r[xi]) for r in rows]
    ys = [parse_number(r[yi]) for r in rows]
    if any(v is None for v in xs + ys):
        raise InputError("some x/y values are not numbers")
    ux, uy = sorted(set(xs)), sorted(set(ys))
    gaps = [b - a for a, b in zip(ux, ux[1:]) if b - a > 1e-9] + \
           [b - a for a, b in zip(uy, uy[1:]) if b - a > 1e-9]
    spacing = min(gaps) if gaps else 1.0
    x0, y0 = ux[0], uy[0]
    return [(int(round((x - x0) / spacing)), int(round((y - y0) / spacing))) for x, y in zip(xs, ys)], \
        "x/y (grid spacing %g)" % spacing


def integration_columns(header):
    return [h for h in header if h.startswith(INTEGRATION_PREFIX)
            and (h == INTEGRATION_PREFIX or h[len(INTEGRATION_PREFIX):].startswith(" R"))]


def vim_column_names(suffix):
    return ["VIM" + suffix, "VIM Map (VIM+5)" + suffix, "VIM Class" + suffix,
            "VIM Band" + suffix, "VIM Hidden Centre Cluster" + suffix, "VIM Cluster Size" + suffix]


# --------------------------------------------------------------------------- report

def summary_text(source, integ_col, cells_from, n_rows, points, ranked, bands_used, args):
    n = len(points)
    lines = []
    w = lines.append
    w("Visibility Integration Mismatch (VIM) - depthmapx-vim %s" % __version__)
    w("Input:            %s" % source)
    w("Integration:      %s   (connectivity: %s)" % (integ_col, args.connectivity))
    w("Grid from:        %s" % cells_from)
    w("Parameters:       bands=%d (used %d), min points/band=%d, strong>=%g, moderate>=%g, false<=%g"
      % (args.bands, bands_used, args.min_per_band, args.strong, args.moderate, args.false_centre))
    w("Points:           %d valid of %d rows" % (n, n_rows))
    if n == 0:
        return "\n".join(lines) + "\n"
    sizes = {}
    for p in points:
        sizes[p.band] = sizes.get(p.band, 0) + 1
    w("Bands:            %d, %d-%d points each" % (len(sizes), min(sizes.values()), max(sizes.values())))
    vims = [p.vim for p in points]
    w("VIM range:        %.3f to %.3f" % (min(vims), max(vims)))
    w("")
    w("Class counts")
    for c in (0, 1, 2, 3):
        k = sum(1 for p in points if p.cls == c)
        w("  %d %-24s %8d  %5.1f%%" % (c, CLASS_NAMES[c], k, 100.0 * k / n))
    w("")
    conn = [p.connectivity for p in points]
    integ = [p.integration for p in points]
    w("Correlations (Pearson r)")
    w("  intelligibility r(Connectivity, Integration) = %.3f" % pearson(conn, integ))
    w("  r(VIM, Integration)                          = %.3f" % pearson(vims, integ))
    w("  r(VIM, Connectivity)                         = %.3f" % pearson(vims, conn))
    w("")
    single = sum(1 for _, m, _ in ranked if len(m) == 1)
    w("Hidden-centre clusters: %d (single-cell: %d, 5+ cells: %d)"
      % (len(ranked), single, sum(1 for _, m, _ in ranked if len(m) >= 5)))
    if ranked:
        w("  rank  cells  summed VIM  max VIM   centre x     centre y")
        for r, members, vim_sum in ranked[:args.top]:
            w("  %4d  %5d  %10.2f  %7.2f   %-11s  %-11s"
              % (r, len(members), vim_sum, max(m.vim for m in members),
                 fmt(sum(m.x for m in members) / len(members)),
                 fmt(sum(m.y for m in members) / len(members))))
        if len(ranked) > args.top:
            w("  ... %d more in the clusters file" % (len(ranked) - args.top))
    return "\n".join(lines) + "\n"


# --------------------------------------------------------------------------- map figures (SVG)
#
# Same layout and colours as the web page's PNG export. Diverging blue <-> red with a neutral
# gray midpoint; the four class steps were validated all-pairs for colour-vision deficiency.

CLASS_COLOURS = ["#2a78d6", "#d3d2cb", "#e86e6d", "#8f1f1e"]
NO_VALUE_COLOUR = "#efeee9"
RAMP_STOPS = [(-5.0, "#104281"), (-2.5, "#3987e5"), (0.0, "#d3d2cb"), (2.5, "#e86e6d"), (5.0, "#8f1f1e")]
INK, INK_2, MUTED, SURFACE = "#0b0b0b", "#52514e", "#898781", "#fcfcfb"
FONT = "system-ui, -apple-system, 'Segoe UI', Helvetica, Arial, sans-serif"


def _to_lab(hex_colour):
    def lin(c):
        return c / 12.92 if c <= 0.04045 else ((c + 0.055) / 1.055) ** 2.4
    r, g, b = (lin(int(hex_colour[i:i + 2], 16) / 255.0) for i in (1, 3, 5))
    l_ = (0.4122214708 * r + 0.5363325363 * g + 0.0514459929 * b) ** (1 / 3)
    m_ = (0.2119034982 * r + 0.6806995451 * g + 0.1073969566 * b) ** (1 / 3)
    s_ = (0.0883024619 * r + 0.2817188376 * g + 0.6299787005 * b) ** (1 / 3)
    return (0.2104542553 * l_ + 0.7936177850 * m_ - 0.0040720468 * s_,
            1.9779984951 * l_ - 2.4285922050 * m_ + 0.4505937099 * s_,
            0.0259040371 * l_ + 0.7827717662 * m_ - 0.8086757660 * s_)


def _from_lab(lab):
    L, a, b = lab
    l_ = (L + 0.3963377774 * a + 0.2158037573 * b) ** 3
    m_ = (L - 0.1055613458 * a - 0.0638541728 * b) ** 3
    s_ = (L - 0.0894841775 * a - 1.2914855480 * b) ** 3
    rgb = (4.0767416621 * l_ - 3.3077115913 * m_ + 0.2309699292 * s_,
           -1.2684380046 * l_ + 2.6097574011 * m_ - 0.3413193965 * s_,
           -0.0041960863 * l_ - 0.7034186147 * m_ + 1.7076147010 * s_)

    def delin(c):
        return 12.92 * c if c <= 0.0031308 else 1.055 * c ** (1 / 2.4) - 0.055
    return "#%02x%02x%02x" % tuple(int(round(min(1.0, max(0.0, delin(c))) * 255)) for c in rgb)


def _ramp_lut():
    labs = [(v, _to_lab(h)) for v, h in RAMP_STOPS]
    lut = []
    for i in range(201):
        v = -5.0 + i * 0.05
        j = 0
        while j < len(labs) - 2 and v > labs[j + 1][0]:
            j += 1
        (v0, a), (v1, b) = labs[j], labs[j + 1]
        t = min(1.0, max(0.0, (v - v0) / (v1 - v0)))
        lut.append(_from_lab(tuple(x + (y - x) * t for x, y in zip(a, b))))
    return lut


_RAMP = None


def ramp_colour(vim):
    global _RAMP
    if _RAMP is None:
        _RAMP = _ramp_lut()
    return _RAMP[int(round((max(-5.0, min(5.0, vim)) + 5.0) / 0.05))]


def _esc(text):
    return str(text).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;").replace('"', "&quot;")


def map_svg(view, cells, points, ranked, title, subtitle, args, units="m"):
    """A print-ready figure (title, map, scale bar, legend) as an SVG string.

    view: "class" (VIM Class) or "vim" (continuous VIM, fixed -5..+5 scale).
    cells: grid cell of every input row (rows without valid values are drawn as 'no value').
    """
    gx_all = [c[0] for c in cells]
    gy_all = [c[1] for c in cells]
    min_x, max_x, min_y, max_y = min(gx_all), max(gx_all), min(gy_all), max(gy_all)
    cols, rows = max_x - min_x + 1, max_y - min_y + 1
    s = float(max(2, min(12, 1100 // cols)))  # whole pixels, so neighbouring cells never show seams
    w0, h0 = cols * s, rows * s
    m, fs = 40.0, 15.0
    width = max(w0, 760.0) + 2 * m
    by_cell = {(p.gx, p.gy): p for p in points}

    # legend layout (wrapping), so the height fits
    if view == "class":
        items = [(CLASS_COLOURS[0], "False centre (VIM ≤ %g)" % args.false_centre),
                 (CLASS_COLOURS[1], "As expected"),
                 (CLASS_COLOURS[2], "Moderate hidden centre (%g ≤ VIM < %g)" % (args.moderate, args.strong)),
                 (CLASS_COLOURS[3], "Strong hidden centre (VIM ≥ %g)" % args.strong)]
        lines, used = [[]], 0.0
        for colour, text in items:
            w = fs * 1.3 + len(text) * fs * 0.55 + fs * 1.6
            if used + w > width - 2 * m and lines[-1]:
                lines.append([])
                used = 0.0
            lines[-1].append((colour, text, used))
            used += w
    else:
        lines = [[]]
    spacing = None
    for p in points[1:]:
        if p.gx != points[0].gx:
            spacing = abs((p.x - points[0].x) / float(p.gx - points[0].gx))
            break
    title_h, scale_h, line_h = fs * 4.2, (fs * 2.6 if spacing else fs), fs * 1.8
    height = m + title_h + h0 + scale_h + len(lines) * line_h + fs * 2.4 + m * 0.4
    mx, my = float(round(m + (width - 2 * m - w0) / 2.0)), float(round(m + title_h))

    out = ['<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d" viewBox="0 0 %.1f %.1f" '
           'font-family="%s">' % (round(width), round(height), width, height, FONT),
           '<rect width="100%%" height="100%%" fill="%s"/>' % SURFACE,
           '<text x="%.1f" y="%.1f" font-size="%.1f" font-weight="650" fill="%s">%s</text>'
           % (m, m + fs * 1.3, fs * 1.6, INK, _esc(title)),
           '<text x="%.1f" y="%.1f" font-size="%.1f" fill="%s">%s</text>'
           % (m, m + fs * 2.9, fs, INK_2, _esc(subtitle)),
           '<g shape-rendering="crispEdges">']
    # map: merge horizontal runs of one colour into a single rect
    grid = {}
    for gx, gy in cells:
        p = by_cell.get((gx, gy))
        if p is None:
            colour = NO_VALUE_COLOUR
        elif view == "class":
            colour = CLASS_COLOURS[p.cls]
        else:
            colour = ramp_colour(p.vim)
        grid[(gx - min_x, max_y - gy)] = colour
    for row in range(rows):
        col = 0
        while col < cols:
            colour = grid.get((col, row))
            if colour is None:
                col += 1
                continue
            end = col
            while end + 1 < cols and grid.get((end + 1, row)) == colour:
                end += 1
            out.append('<rect x="%.2f" y="%.2f" width="%.2f" height="%.2f" fill="%s"/>'
                       % (mx + col * s, my + row * s, (end - col + 1) * s, s, colour))
            col = end + 1
    out.append("</g>")

    # numbered top clusters
    r = fs * 0.78
    for rank, members, _ in ranked[:args.top]:
        cx = mx + (sum(p.gx for p in members) / len(members) - min_x + 0.5) * s
        cy = my + (max_y - sum(p.gy for p in members) / len(members) + 0.5) * s
        out.append('<circle cx="%.1f" cy="%.1f" r="%.1f" fill="%s" stroke="%s" stroke-width="1.5"/>'
                   % (cx, cy, r, SURFACE, INK))
        out.append('<text x="%.1f" y="%.1f" font-size="%.1f" font-weight="600" fill="%s" '
                   'text-anchor="middle" dominant-baseline="central">%d</text>' % (cx, cy, fs * 0.85, INK, rank))

    # scale bar, right-aligned under the map
    if spacing:
        target = cols * spacing / 5.0
        power = 10.0 ** math.floor(math.log10(target))
        nice = min((k * power for k in (1, 2, 5, 10)), key=lambda v: abs(v - target))
        bar = nice / spacing * s
        x1, y = mx + w0, my + h0 + fs * 0.9
        x0 = x1 - bar
        out.append('<path d="M%.1f %.1fH%.1fM%.1f %.1fV%.1fM%.1f %.1fV%.1f" stroke="%s" stroke-width="2" fill="none"/>'
                   % (x0, y, x1, x0, y - 5, y + 5, x1, y - 5, y + 5, INK))
        out.append('<text x="%.1f" y="%.1f" font-size="%.1f" fill="%s" text-anchor="end" dominant-baseline="central">'
                   '%s</text>' % (x0 - 8, y, fs * 0.9, INK_2, _esc(("%g" % float("%.3g" % nice)) + (" " + units if units else ""))))

    # legend
    ly0 = my + h0 + scale_h + line_h / 2
    if view == "class":
        for li, line in enumerate(lines):
            ly = ly0 + li * line_h
            for colour, text, off in line:
                out.append('<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" rx="2" fill="%s"/>'
                           % (m + off, ly - fs * 0.45, fs * 0.9, fs * 0.9, colour))
                out.append('<text x="%.1f" y="%.1f" font-size="%.1f" fill="%s" dominant-baseline="central">%s</text>'
                           % (m + off + fs * 1.3, ly, fs, INK_2, _esc(text)))
    else:
        label_w = fs * 7.2
        out.append('<text x="%.1f" y="%.1f" font-size="%.1f" fill="%s" dominant-baseline="central">'
                   '−5 false centre</text>' % (m, ly0, fs, INK_2))
        out.append('<defs><linearGradient id="vimramp">%s</linearGradient></defs>'
                   % "".join('<stop offset="%.3f" stop-color="%s"/>' % (i / 20.0, ramp_colour(-5 + i * 0.5))
                             for i in range(21)))
        out.append('<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" rx="2" fill="url(#vimramp)"/>'
                   % (m + label_w, ly0 - fs * 0.45, fs * 14, fs * 0.9))
        out.append('<text x="%.1f" y="%.1f" font-size="%.1f" fill="%s" dominant-baseline="central">'
                   '+5 hidden centre · gray = as expected (VIM 0)</text>'
                   % (m + label_w + fs * 14.6, ly0, fs, INK_2))
    out.append('<text x="%.1f" y="%.1f" font-size="%.1f" fill="%s">depthmapx-vim %s · VIM = (integration '
               '− band median) ÷ band SD</text>' % (m, height - m * 0.5, fs * 0.8, MUTED, __version__))
    out.append("</svg>")
    return "\n".join(out) + "\n"


# --------------------------------------------------------------------------- main

def run(args):
    header, rows, delimiter = read_table(args.input)
    if args.connectivity not in header:
        raise InputError("column '%s' not found. Run visibility analysis (Calculate visibility "
                         "relationships) before exporting." % args.connectivity)
    integ_cols = args.integration or integration_columns(header)
    if not integ_cols:
        raise InputError("no '%s' column found. Run visibility analysis with 'Include global "
                         "measures' before exporting, or name the column with --integration."
                         % INTEGRATION_PREFIX)
    for c in integ_cols:
        if c not in header:
            raise InputError("column '%s' not found" % c)

    cells, cells_from = grid_cells(header, rows)
    xi = header.index("x") if "x" in header else None
    yi = header.index("y") if "y" in header else None
    conn_i = header.index(args.connectivity)

    # drop VIM columns the input may already have (e.g. from the modified depthmapX build)
    out_header = list(header)
    out_rows = [list(r) + [""] * (len(header) - len(r)) for r in rows]
    for integ in integ_cols:
        suffix = integ[len(INTEGRATION_PREFIX):] if integ.startswith(INTEGRATION_PREFIX) else ""
        for name in vim_column_names(suffix):
            if name in out_header:
                i = out_header.index(name)
                out_header.pop(i)
                for r in out_rows:
                    r.pop(i)

    first_vim_col = len(out_header)  # VIM columns are appended from here
    base, _ = os.path.splitext(args.output or args.input)
    out_path = args.output or base + "_vim" + (".txt" if delimiter == "\t" else ".csv")
    stem = os.path.splitext(out_path)[0]
    reports = []

    for integ in integ_cols:
        suffix = integ[len(INTEGRATION_PREFIX):] if integ.startswith(INTEGRATION_PREFIX) else ""
        integ_i = header.index(integ)
        points = []
        for idx, r in enumerate(rows):
            c = parse_number(r[conn_i]) if conn_i < len(r) else None
            v = parse_number(r[integ_i]) if integ_i < len(r) else None
            if c is None or v is None:
                continue
            c, v = to_float32(c), to_float32(v)
            if c < 0.0 or v < 0.0:
                continue
            x = parse_number(r[xi]) if xi is not None and xi < len(r) else None
            y = parse_number(r[yi]) if yi is not None and yi < len(r) else None
            points.append(VimPoint(idx, cells[idx][0], cells[idx][1], c, v, x, y))
        points.sort(key=lambda p: (p.gx, p.gy))  # depthmapX's visiting order

        ranked = compute_vim(points, args.bands, args.min_per_band,
                             args.strong, args.moderate, args.false_centre)
        bands_used = max(1, min(max(1, args.bands), len(points) // max(1, args.min_per_band))) if points else 0

        new_vals = [[fmt(MISSING)] * 6 for _ in rows]
        for p in points:
            mapped = max(-5.0, min(5.0, p.vim)) + 5.0
            new_vals[p.index] = [fmt(to_float32(p.vim)), fmt(to_float32(mapped)), str(p.cls),
                                 str(p.band + 1), str(p.cluster), str(p.cluster_size)]
        out_header += vim_column_names(suffix)
        for r, vals in zip(out_rows, new_vals):
            r.extend(vals)
        reports.append((suffix, integ, points, ranked, bands_used))

    with open(out_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, delimiter=delimiter, lineterminator="\n")
        writer.writerow(out_header)
        writer.writerows(out_rows)
    written = [out_path]

    # Ref + VIM columns only, for attaching back onto the VGA map in depthmapX:
    #   depthmapXcli -m IMPORT -f plan.graph -o plan_vim.graph -if <this file> -iaa
    if not args.no_join and cells_from == "Ref":
        ref_i = header.index("Ref")
        path = stem + "_join.csv"
        with open(path, "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh, lineterminator="\n")
            w.writerow(["Ref"] + out_header[first_vim_col:])
            for src, r in zip(rows, out_rows):
                w.writerow([src[ref_i]] + r[first_vim_col:])
        written.append(path)

    for suffix, integ, points, ranked, bands_used in reports:
        tag = suffix.strip().replace(" ", "_")
        tag = "_" + tag if tag else ""
        if not args.no_summary:
            path = "%s_summary%s.txt" % (stem, tag)
            text = summary_text(os.path.basename(args.input), integ, cells_from, len(rows),
                                points, ranked, bands_used, args)
            with open(path, "w", encoding="utf-8") as fh:
                fh.write(text)
            written.append(path)
            if not args.quiet:
                print(text)
        if not args.no_clusters:
            path = "%s_clusters%s.csv" % (stem, tag)
            with open(path, "w", newline="", encoding="utf-8") as fh:
                w = csv.writer(fh, lineterminator="\n")
                w.writerow(["Cluster Rank", "Cells", "Summed VIM", "Mean VIM", "Max VIM",
                            "Centre x", "Centre y", "Min x", "Max x", "Min y", "Max y"])
                for r, members, vim_sum in ranked:
                    xs = [m.x for m in members]
                    ys = [m.y for m in members]
                    w.writerow([r, len(members), fmt(vim_sum), fmt(vim_sum / len(members)),
                                fmt(max(m.vim for m in members)),
                                fmt(sum(xs) / len(xs)), fmt(sum(ys) / len(ys)),
                                fmt(min(xs)), fmt(max(xs)), fmt(min(ys)), fmt(max(ys))])
            written.append(path)
        if not args.no_map and points:
            subtitle = "%s · %s · %s points · %d connectivity bands" % (
                os.path.basename(args.input), integ, "{:,}".format(len(points)), bands_used)
            for view, title in (("class", "VIM Class"), ("vim", "Visibility Integration Mismatch (VIM)")):
                path = "%s_map_%s%s.svg" % (stem, view, tag)
                with open(path, "w", encoding="utf-8") as fh:
                    fh.write(map_svg(view, cells, points, ranked, title, subtitle, args, args.units))
                written.append(path)

    if not args.quiet:
        print("Wrote:")
        for p in written:
            print("  " + p)
    return written


def build_parser():
    ap = argparse.ArgumentParser(
        prog="depthmapx_vim",
        description="Compute Visibility Integration Mismatch (VIM) from a depthmapX VGA export.",
        epilog="Export from depthmapX with File > Export (.csv or .txt) after visibility "
               "analysis with global measures. Writes <input>_vim.csv (all original columns "
               "plus VIM columns), a summary and a clusters table.")
    ap.add_argument("input", help="depthmapX point-map export (.csv or tab-delimited .txt)")
    ap.add_argument("-o", "--output", help="output table (default: <input>_vim.csv)")
    ap.add_argument("--integration", action="append", metavar="COLUMN",
                    help="integration column to use (repeatable). Default: every "
                         "'Visual Integration [HH]' column, including radius versions like ' R3'")
    ap.add_argument("--connectivity", default=CONNECTIVITY, metavar="COLUMN",
                    help="connectivity column (default: %(default)s)")
    ap.add_argument("--bands", type=int, default=DEFAULT_BANDS,
                    help="number of connectivity bands (default: %(default)s)")
    ap.add_argument("--min-per-band", type=int, default=DEFAULT_MIN_PER_BAND,
                    help="reduce bands on small plans to keep at least this many points per band "
                         "(default: %(default)s)")
    ap.add_argument("--strong", type=float, default=2.0, help="strong hidden centre: VIM >= this (default: 2)")
    ap.add_argument("--moderate", type=float, default=1.0, help="moderate hidden centre: VIM >= this (default: 1)")
    ap.add_argument("--false-centre", type=float, default=-1.0, help="false centre: VIM <= this (default: -1)")
    ap.add_argument("--top", type=int, default=10, help="clusters listed in the summary (default: 10)")
    ap.add_argument("--no-summary", action="store_true", help="do not write the summary file")
    ap.add_argument("--no-clusters", action="store_true", help="do not write the clusters file")
    ap.add_argument("--no-join", action="store_true",
                    help="do not write the Ref + VIM join file used to attach VIM back onto the map")
    ap.add_argument("--no-map", action="store_true", help="do not write the SVG map figures")
    ap.add_argument("--units", default="m", help="plan units for the scale bar (default: m; '' for none)")
    ap.add_argument("-q", "--quiet", action="store_true", help="print nothing")
    ap.add_argument("--version", action="version", version="%(prog)s " + __version__)
    return ap


def main(argv=None):
    args = build_parser().parse_args(argv)
    try:
        run(args)
    except (InputError, OSError) as e:
        print("depthmapx_vim: error: %s" % e, file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
