# depthmapx-vim

**Visibility Integration Mismatch (VIM)** for [depthmapX](https://github.com/SpaceGroupUCL/depthmapX)
visibility graph analysis (VGA).

Intelligibility is a single correlation for a whole plan: do places that *look* open (high
connectivity) also turn out to be *configurationally central* (high visual integration)? VIM
localises that question to every point on the grid. Each point is compared only with points that
have a similar-sized visual field, so VIM shows **where** a plan is unintelligible:

- **hidden centres:** more integrated than their view suggests (central, but you can't see it)
- **false centres:** look open but sit at the edge of the layout

It works with the depthmapX you already have: export your VGA results from depthmapX, then get VIM
values, a summary and publication-ready maps back. You don't need a special build of depthmapX.
There are two ways to run it, and both give identical results:

- **In your browser, nothing to install:** **[open the VIM web page](https://jejenjo.github.io/depthmapx-vim/)**
  and drop in your export. It runs entirely on your computer; the file is never uploaded. (You can
  also download `docs/index.html` and open it offline.)
- **As a Python script** (below), for batch processing, scripting or reproducible pipelines.

---

## Requirements

- **Python 3.8 or newer** ([python.org](https://www.python.org/downloads/)). Nothing else to
  install; the script only uses Python's standard library.
- **depthmapX**: any version whose VGA export contains `Connectivity` and
  `Visual Integration [HH]` columns. Tested with depthmapX 0.9.1.

## Quick start

### 1. Run visibility analysis in depthmapX (as usual)

1. Import your plan, **Set Grid** and fill the spaces.
2. **Tools → Visibility → Make Visibility Graph**.
3. **Tools → Visibility → Run Visibility Graph Analysis**: tick **Calculate visibility
   relationships** and **Include global measures** (radius `n`).
4. Save your graph (**File → Save**), e.g. `plan.graph`.
5. With the visibility graph displayed, **File → Export** and save as `.csv` (or the default
   tab-delimited `.txt`), e.g. `plan_vga.csv`.

### 2. Compute VIM

Open a terminal (Windows: *Command Prompt* or *PowerShell*; macOS: *Terminal*) in the folder that
contains `depthmapx_vim.py` and run:

```
python depthmapx_vim.py path/to/plan_vga.csv
```

(On macOS/Linux use `python3`.) The script prints a summary and writes, next to your export:

| File | Contents |
|---|---|
| `plan_vga_vim.csv` | Your full export with the six VIM columns added: for Excel, QGIS, R, SPSS… |
| `plan_vga_vim_summary.txt` | Class counts, band sizes, correlations and the top hidden-centre clusters |
| `plan_vga_vim_clusters.csv` | Every hidden-centre cluster: rank, size, summed/mean/max VIM, centre and extent |
| `plan_vga_vim_join.csv` | `Ref` + VIM columns only, for attaching VIM back onto the map (step 3) |
| `plan_vga_vim_map_class.svg` | Map figure of VIM Class, with legend, numbered top clusters and scale bar |
| `plan_vga_vim_map_vim.svg` | Map figure of continuous VIM on a fixed −5…+5 scale |

The map figures are SVG, a vector format that stays sharp at any print size. They open in any web
browser and in Inkscape or Illustrator, and Word and PowerPoint (Microsoft 365) can insert them
directly. For a PNG, use the web page's **Map PNG** button, which gives the same figure at print
resolution.

### 3. See VIM on the map in depthmapX (optional)

depthmapX's command-line companion, `depthmapXcli`, can attach the join file to the visibility
graph in your saved `.graph` file:

```
depthmapXcli -m IMPORT -f plan.graph -o plan_vim.graph -if plan_vga_vim_join.csv -iaa
```

Open `plan_vim.graph` in depthmapX. The VIM columns are now in the attribute list of your
visibility graph, and you can colour the map with them like any built-in measure. (Check that your
`depthmapXcli -h` lists `-iaa` under the IMPORT options.)

Without `depthmapXcli`, you can map `plan_vga_vim.csv` in QGIS (*Layer → Add Delimited Text
Layer*, X = `x`, Y = `y`) or in any GIS or plotting tool.

---

## Output columns

| Column | Values | Use |
|---|---|---|
| `VIM` | z-score, negative and positive | The value for analysis and statistics |
| `VIM Map (VIM+5)` | 0–10: VIM clipped to ±5, plus 5 (5 = as expected) | Colouring the map inside depthmapX |
| `VIM Class` | 0 false centre (VIM ≤ −1) · 1 as expected · 2 moderate hidden (1 ≤ VIM < 2) · 3 strong hidden (VIM ≥ 2) | Categorical map |
| `VIM Band` | 1–40 | The connectivity band the point was compared within |
| `VIM Hidden Centre Cluster` | 0 = not in a cluster; 1 = top-ranked cluster, 2 = next, … | Locate and rank hidden centres |
| `VIM Cluster Size` | number of grid cells in that point's cluster | Filter out single-cell clusters |

Points without a valid connectivity or integration value get `-1` in every VIM column, depthmapX's
own "no value".

If you analysed with a radius (e.g. 3), the columns get the same suffix (`VIM R3`, …). When an
export has several `Visual Integration [HH]` columns, VIM is calculated for each.

**Map colours.** The figures use a diverging blue ↔ red scale with a neutral gray midpoint: blue =
false centres, gray = as expected, red = hidden centres, darker = stronger. The scale is fixed
(gray is always VIM = 0, the ends are ±5), so the same colour means the same VIM on every plan. The
four class colours were checked for colour-vision deficiency: every pair stays distinguishable under
simulated protanopia and deuteranopia.

**Displaying VIM in depthmapX.** depthmapX treats every negative number as "no value", so false
centres show grey or blank when you colour by `VIM`. Colour by `VIM Map (VIM+5)` or `VIM Class`
instead. depthmapX also stretches colours between each column's minimum and maximum, so for colours
that mean the same thing on every plan, set **View → Colour Range** to 0–10 for `VIM Map (VIM+5)`.

## Method

1. **Valid points** are those with connectivity ≥ 0 and integration ≥ 0. Visual integration is
   `Visual Integration [HH]`.
2. **Bands.** Points are ranked by connectivity and split into *B* equal-sized quantile bands
   (*B* = 40 by default). Points with identical connectivity always share a band: the whole tie
   group goes to the band of its middle rank, so bands are only approximately equal in size. On
   small plans *B* is reduced to ⌊*n*/10⌋ so each band keeps at least ~10 points. Ranking on a log
   scale gives the same order, so log-transforming connectivity changes nothing.
3. **Expected integration.** Within each band *b*: the median and the sample standard deviation
   (*n* − 1) of integration.
4. **VIM** for point *i* in band *b*:

   VIM<sub>i</sub> = (Int<sub>i</sub> − median(Int)<sub>b</sub>) / SD(Int)<sub>b</sub>
   (VIM = 0 when SD = 0)

5. **Classes:** strong hidden centre VIM ≥ 2; moderate hidden centre 1 ≤ VIM < 2; false centre
   VIM ≤ −1; otherwise as expected.
6. **Clusters.** Strong hidden-centre points are grouped into clusters of touching grid cells,
   diagonals included (8-neighbour). Clusters are ranked by summed VIM, so a cluster ranks high by
   being strong, large or both.

Because comparison happens *within* connectivity bands, VIM is independent of visual-field size by
construction. On the depthmapX sample gallery plan, r(VIM, connectivity) = 0.07 while
r(VIM, integration) = 0.51.

### Options

```
python depthmapx_vim.py plan_vga.csv [options]

  -o, --output FILE       output table (default: <input>_vim.csv)
  --integration COLUMN    integration column (repeatable; default: every Visual Integration [HH] column)
  --connectivity COLUMN   connectivity column (default: Connectivity)
  --bands N               number of connectivity bands (default: 40)
  --min-per-band N        minimum points per band on small plans (default: 10)
  --strong X              strong hidden centre threshold (default: 2)
  --moderate X            moderate hidden centre threshold (default: 1)
  --false-centre X        false centre threshold (default: -1)
  --top N                 clusters listed in the summary and numbered on the maps (default: 10)
  --units TEXT            plan units for the map scale bar (default: m; "" for none)
  --no-summary, --no-clusters, --no-join, --no-map   skip those output files
  -q, --quiet             print nothing
```

The defaults are the method as published; the band and threshold options are there for
sensitivity analysis.

## Validation

The script is a line-by-line port of a C++ implementation of VIM built into a modified depthmapX
0.9.1. It follows the same summation order and the same 32-bit storage of values, so on the same
data the two give the **same values to depthmapX's export precision**:

| Test | Points | VIM values identical to C++ |
|---|---|---|
| depthmapX sample gallery plan, radius n (CSV) | 4,332 | 25,992 / 25,992 |
| same, tab-delimited `.txt` export | 4,332 | 25,992 / 25,992 |
| same, without the `Ref` column (grid rebuilt from x/y) | 4,332 | 25,992 / 25,992 |
| sample gallery plan, radius 3 (`VIM R3`) | 4,332 | 25,992 / 25,992 |
| 116 × 33 m hospital floor plan, 0.5 m grid | 13,672 | 82,032 / 82,032 |

Attaching the join file to a standard VGA graph with `depthmapXcli -iaa` and exporting again gives
a file identical, cell for cell, to the C++ build's own export (4,332 rows × 19 columns).

The web page uses a JavaScript port of the same code. On all five tests above, its results table,
join file and clusters table are byte-for-byte identical to the Python script's.

Reference figures for the depthmapX sample gallery plan (`testdata/gallery_connected.graph` in the
depthmapX repository, radius n): 4,332 points in 40 bands of 90–127 points; classes 645 false
centres, 2,880 as expected, 645 moderate and 162 strong hidden centres; 39 clusters.

### Running the tests

```
python -m unittest discover tests
```

## Citation

If you use VIM in published work, please cite:

> Jogy, J. J., & Dev, G. (2026). *Visibility Integration Mismatch: A Space Syntax Measure in
> depthmapX for Locating Hospital Circulation Conflicts, a Case of Thrissur.* Unpublished
> manuscript.

See `CITATION.cff` for a machine-readable citation. GitHub's **Cite this repository** button uses
it.

## License

Copyright © 2026 Jennat Jennifer Jogy and Govind Dev.

This program is free software: you can redistribute it and/or modify it under the terms of the
GNU General Public License as published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version. It is distributed in the hope that it will be
useful, but WITHOUT ANY WARRANTY. See [`LICENSE`](LICENSE) for the full text.
