# Recovering Bulgaria's unpublished "1970 system" zones (K3/K5/K7/K9) from BGS2005

## The problem

Bulgaria's cadastre agency (AGKK) uses a historical "Coordinate system 1970"
(KS1970) split into four zones — **K3, K5, K7, K9** — each a Lambert Conformal
Conic projection on the **Pulkovo 1942(58) datum (Krasovsky 1940 ellipsoid)**.
Its defining parameters were never published: the `.prj` files AGKK's own tool
writes name the datum and the projection *type* and then carry **no numeric
parameter at all**.

That is not a documentation gap you can look up. Bulgaria's 1970-system zone
parameters and its datum-shift parameters from BGS2005 were never officially
released, and existing open-source work (e.g. `bojko108/transformations`) had to
reverse-engineer a correction grid by sampling AGKK's tool rather than working
from a formula.

This repository recovers the transformation from data — twice, by two
independent methods — and reports honestly what the data does and does not
determine.

## The data

`source_data/` holds, for each zone `n` ∈ {3,5,7,9} and each sampling step
`r` ∈ {1km, 2km}:

| file | contents |
|---|---|
| `bg_k{n}_{r}_sel.shp` | points in **EPSG:7801** (`BGS2005 / CCS2005`), a fully documented Lambert Conformal Conic (2SP) on GRS80: central meridian 25.5°E, standard parallels 42°/43°20′, latitude of false origin 42.6678756833333°, false easting 500 000 m, false northing 4 725 824.3591 m — identical across all four zones |
| `bg_k{n}_{r}_sel_tr.shp` | the same points run through AGKK's official tool into KS1970 zone K{n}. Its `.prj` gives the datum and projection type only |
| `bg_k{n}_{r}_sel_tr.log` | AGKK's own report: *"Target zone width: 3 degrees … Plane transform accuracy - 0.14 meters"* |

`bg_ext/bg_zones1970.shp` holds the four official zone extents (one polygon per
zone, field `CLIST` = 3/5/7/9) in WGS 84 / UTM zone 35N. It is what the paper's
map is drawn from.

Each `_sel` / `_sel_tr` pair shares an integer `id`, so every point is an exact
correspondence — no manually surveyed control points are involved.

| zone | quadrant | 1 km points | 2 km points | zone area (`bg_ext`) |
|---|---|---|---|---|
| K3 | northwest | 23 680 | 5 907 | 23 695 km² |
| K5 | southeast | 31 912 | 7 987 | 31 891 km² |
| K7 | northeast | 30 401 | 7 599 | 30 378 km² |
| K9 | southwest | 37 294 | 9 319 | 37 333 km² |

The four zone polygons **tile** the country rather than overlapping: their areas
sum to 123 297 km² against a union of 123 281 km². The point counts track the
areas almost exactly, as a 1 km lattice should.

### Why the two sampling steps matter

The 1 km and 2 km sets are **disjoint** — they share no ids, and their lattices
are offset from each other by 120 m, so no 2 km point coincides with any 1 km
point (verified: zero coordinate collisions in all four zones).

So every model here is built on the 1 km set and scored on the 2 km set, which
is a genuinely independent sample rather than a random split of one lattice.
For the interpolating model that distinction matters a great deal: a random
hold-out of a single lattice leaves test points sitting on the *edges* of
triangles built from their own immediate neighbours, which flatters the result.
The 2 km points instead land in the **interior** of the 1 km cells — the worst
case for piecewise-linear interpolation. Every number below is therefore
conservative.

## Two solutions per zone

| | `fit_analytic_lcc.py` | `build_tin_grid.py` |
|---|---|---|
| Output | closed-form PROJ pipeline (LCC + 2D similarity) | `tinshift` triangulated grid (JSON) |
| Coverage | anywhere (extrapolates) | only within the zone's sampled territory |
| Portability | one short PROJ string, embeddable anywhere | one 2–4 MB JSON file PROJ must load |
| Accuracy | 0.09–0.14 m RMSE | 0.06–0.24 **mm** RMSE |
| Use when | you need a compact formula, or points outside the sampled area | you want survey-grade accuracy (recommended) |

### Closed-form model — out-of-sample accuracy (2 km set)

| zone | cone constant *n* | tangent parallel | rotation | scale | RMSE | p95 | max |
|---|---|---|---|---|---|---|---|
| K3 | 0.687702042 | 43.448480° | +5704.8″ | −6.47 ppm | 0.100 m | 0.189 m | 0.298 m |
| K5 | 0.675158040 | 42.466424° | −2161.3″ | −3.08 ppm | 0.143 m | 0.267 m | 0.466 m |
| K7 | 0.689025227 | 43.552997° | −1812.4″ | +1.34 ppm | 0.120 m | 0.237 m | 0.330 m |
| K9 | 0.672864894 | 42.288566° | +5041.6″ | −8.63 ppm | 0.085 m | 0.169 m | 0.270 m |

This **matches AGKK's own declared 0.14 m plane accuracy** for its tool. The
fit-set and out-of-sample RMSE agree to three decimals in every zone, so this
residual is model bias, not overfitting — a single global conic simply cannot
express the local structure of the historical 1970 network adjustment.

### Grid model — out-of-sample accuracy (2 km set)

| zone | vertices | triangles | size | covered | median | p95 | max | residuals > 1 mm |
|---|---|---|---|---|---|---|---|---|
| K3 | 23 680 | 46 491 | 2.21 MB | 5 817/5 907 | 0.211 mm | 0.213 mm | 121.1 mm | 28 |
| K5 | 31 912 | 62 938 | 3.01 MB | 7 868/7 987 | 0.055 mm | 0.055 mm | 0.055 mm | 0 |
| K7 | 30 401 | 59 976 | 2.86 MB | 7 530/7 599 | 0.240 mm | 0.242 mm | 0.242 mm | 0 |
| K9 | 37 294 | 73 671 | 3.53 MB | 9 219/9 319 | 0.104 mm | 0.105 mm | 111.3 mm | 40 |

Uncovered points are 2 km points falling just outside the 1 km lattice's outline
(the two sets do not span exactly the same rectangle); they are excluded rather
than extrapolated. Round-trip error, forward then inverse, is under 3 nm.

## Two findings worth stating plainly

### 1. The conic parameters are not recoverable — only the cone constant is

It is tempting to fit the textbook parameter set (λ₀, φ₁, φ₂) plus a similarity
and publish the fitted λ₀/φ₁/φ₂ as "the recovered zone parameters". **That would
be wrong**, and an earlier version of this project did it. Those seven
parameters are not identifiable from planar correspondences:

- changing λ₀ rotates the projection plane about the cone apex, which the
  similarity's own rotation θ absorbs exactly — only θ − n·λ₀ is determined;
- changing (φ₁, φ₂) at constant cone constant *n* changes only the projection's
  overall scale factor, which the similarity's scale *s* absorbs exactly.

Measured, on K3 (`fit_analytic_lcc.py --identifiability`):

| start λ₀ | fitted λ₀ | φ₁ | φ₂ | scale−1 | cone *n* | RMSE |
|---|---|---|---|---|---|---|
| 23.0° | 23.119184° | 43.509° | 43.388° | −5.9 ppm | 0.687702042 | 0.1000 m |
| 25.5° | 24.354318° | 39.381° | 47.425° | +2455.2 ppm | 0.687702042 | 0.1000 m |
| 27.5° | 25.345188° | 53.044° | 33.302° | +15010.6 ppm | 0.687702042 | 0.1000 m |

Parameters differing by **degrees** and by 15 000 ppm of scale produce
coordinates agreeing to **under one micrometre**, and agree on *n* to nine
significant figures. So this repository fits the identifiable 5-parameter form
directly — *n*, θ, *s*, tₓ, t_y, with λ₀ pinned to 25.5° as a stated convention —
and reports *n* as the recovered quantity. Same accuracy, no flat directions,
and nothing claimed that the data does not determine.

### 2. AGKK's own transformation is discontinuous along 24°E

Along any straight row of the lattice, a smooth projection produces a second
difference of a couple of millimetres at 1 km spacing. In K3 and K9 there is a
localized spike two samples wide, sitting exactly on the 24°E meridian:

| zone | seam longitude | source easting | kink | typical curvature |
|---|---|---|---|---|
| K3 | 24.0025° | 378 880 m | 88 mm | 2.2 mm |
| K9 | 23.9971° | 375 880 m | 127 mm | 1.0 mm |

This is a property of AGKK's transformation, not of either dataset or of the
triangulation here — building the grid from the 2 km set and scoring the 1 km
set reproduces it at the same eastings. It is almost certainly a panel boundary
in the tool's own internal correction grid. K5 and K7 show nothing comparable.

The consequence: **all 68 out-of-sample residuals above 1 mm in the entire
project lie on this line.** The grid reproduces the kink exactly at sampled
points but must interpolate straight across it in the one cell that straddles
it. Excluding those points, every zone's out-of-sample RMSE is its median:
0.211, 0.055, 0.240 and 0.104 mm. No closed-form model can represent the kink at
all. `build_tin_grid.py` detects and reports this automatically.

## Directory layout

```
source_data/     the 16 input shapefiles (bg_k{3,5,7,9}_{1km,2km}_sel[_tr]) + AGKK logs
bg_ext/          bg_zones1970.shp -- the four official zone extents (UTM 35N)
scripts/
  bg_points.py            shared loader: matched points by id, per zone and step;
                           also reads the bg_ext zone polygons
  lcc_math.py             ellipsoidal LCC (Snyder 1987), both the 2SP and the
                           cone-constant form, + 2D similarity. Self-testing.
  optim.py                Levenberg-Marquardt least squares in numpy. Self-testing.
  lattice.py              triangulation of a regular ragged-edged lattice. Self-testing.
  fit_analytic_lcc.py     fits the identifiable closed-form model per zone
  build_tin_grid.py       builds the tinshift grid per zone, detects AGKK seams
  transform_points.py     CLI: transform points forward/reverse, either model
  reproject_shapefile.py  CLI: reproject a whole vector file, ready for QGIS
  build_qgis_crs_wkt.py   generates a paste-into-QGIS Custom CRS per zone
  make_zone_figure.py     renders paper/figures/*.svg from bg_ext, injects them
                           into the paper, and reports the seam measurements
  make_paper_pdf.py       renders the paper to PDF, one file per language
output/
  lcc_affine_fit_k{n}.json      fitted parameters + accuracy report, per zone
  bg_k{n}_tinshift.json         the PROJ tinshift grid (both directions), per zone
  bg_k{n}_ks1970_fitted.wkt     QGIS Custom CRS wrapping the analytic pipeline
  tinshift_accuracy.json        grid accuracy + seam diagnostics, all zones
paper/
  recovering_ks1970.html        the write-up (Bulgarian + English)
  recovering_ks1970_bg.pdf      Bulgarian edition, A4
  recovering_ks1970_en.pdf      English edition, A4
  figures/zones.svg             zone map, generated from bg_ext
  figures/seam.svg              curvature profile showing the 24°E discontinuity
```

## 1. Installing PROJ and its Python bindings

`tinshift` requires **PROJ ≥ 7.2.0**. Check what you have:

```bash
proj    # first line prints the release, e.g. "Rel. 9.7.1"
```

Dependencies are deliberately minimal: **numpy and pyproj**. There is no SciPy
dependency — the Levenberg-Marquardt fit and the triangulation are implemented
directly in `scripts/optim.py` and `scripts/lattice.py`, both of which are
self-testing (run them directly). GDAL is needed only to read the shapefiles.

### Option A — system packages (Debian/Ubuntu, what this repo was built against)

```bash
sudo apt-get install proj-bin libproj-dev gdal-bin libgdal-dev python3-gdal python3-numpy
pip install pyproj
```

`proj-bin` gives you the `proj`/`cs2cs`/`cct` CLI tools; `libproj-dev` is needed
if anything compiles against PROJ; `gdal-bin`/`libgdal-dev`/`python3-gdal` are
needed only because `bg_points.py` and `reproject_shapefile.py` use GDAL's OGR
Python bindings. If you only ever use the already-built `output/*.json` files,
you can skip GDAL entirely.

### Option B — conda/mamba (easiest way to get matching GDAL+PROJ+pyproj)

```bash
conda create -n bgcrs -c conda-forge python=3.12 gdal proj pyproj numpy
conda activate bgcrs
```

This sidesteps the most common PROJ headache: GDAL's Python bindings and
`pyproj` each embed their own copy of PROJ, and if those versions disagree about
where the PROJ resource data lives you get obscure `proj_create` failures.
conda-forge keeps them in lockstep.

### Option C — macOS (Homebrew)

```bash
brew install proj gdal
pip install pyproj numpy
```

### Option D — build PROJ from source (only if your OS ships something too old)

```bash
git clone https://github.com/OSGeo/PROJ.git && cd PROJ
git checkout 9.7.1   # or newer
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build . -j$(nproc)
sudo cmake --install .
```

## 2. Reproducing everything

```bash
# the numeric building blocks check themselves
python3 scripts/lcc_math.py
python3 scripts/optim.py
python3 scripts/lattice.py

# inspect the inputs
python3 scripts/bg_points.py

# fit / build (all four zones; add --zone k5 for one)
python3 scripts/fit_analytic_lcc.py      # -> output/lcc_affine_fit_k{n}.json
python3 scripts/build_tin_grid.py        # -> output/bg_k{n}_tinshift.json

# the identifiability argument from section 1 above, reproduced from the data
python3 scripts/fit_analytic_lcc.py --identifiability --zone k3

# QGIS-ready outputs (run the fit first, they read its JSON)
python3 scripts/build_qgis_crs_wkt.py    # -> output/bg_k{n}_ks1970_fitted.wkt

# the paper's figures, from bg_ext; also re-injects them into the paper
python3 scripts/make_zone_figure.py      # -> paper/figures/{zones,seam}.svg

# the paper as PDF, one file per language
python3 scripts/make_paper_pdf.py        # -> paper/recovering_ks1970_{bg,en}.pdf
python3 scripts/make_paper_pdf.py --lang en    # just one
```

The PDF build needs a Chrome/Chromium binary (pass `--browser` if it is not on
your `PATH`). It is the same rendering you get from Ctrl+P in a browser: page
setup, printing on white, keeping the figure fills and not splitting a figure or
table across a page all live in the paper's own `@media print` block. The figures
stay vector, so text in them remains selectable and searchable in the PDF. The
page picks its language from `?lang=bg` / `?lang=en` in the URL, which is how the
two editions are produced — and is also a convenient way to link someone straight
to one language.

`fit_analytic_lcc.py` and `build_tin_grid.py` both print a full out-of-sample
accuracy report, so every number in this README regenerates from scratch.

## 3. Using the results

### Via the provided CLI (easiest)

```bash
# CCS2005 -> KS1970 K3, single point, uses the grid (falls back to the analytic
# model automatically for any point outside the grid's coverage)
python3 scripts/transform_points.py --zone k3 --xy 300000 4800000

# KS1970 K5 -> CCS2005 (reverse)
python3 scripts/transform_points.py --zone k5 --direction reverse --xy 9500000 4600000

# batch, CSV with columns x,y (any extra columns are preserved)
python3 scripts/transform_points.py --zone k9 --csv in.csv --out out.csv

# force the closed-form model (e.g. for a point far outside that zone)
python3 scripts/transform_points.py --zone k7 --method analytic --xy 500000 4800000
```

### Via pyproj directly

```python
from pyproj import Transformer

# grid (recommended, sub-millimetre inside the zone's coverage)
tr = Transformer.from_pipeline(
    "+proj=pipeline +step +proj=tinshift +file=output/bg_k3_tinshift.json"
)
E, N = tr.transform(300000, 4800000)            # forward: CCS2005 -> KS1970 K3
x, y = tr.transform(E, N, direction="INVERSE")  # reverse: KS1970 K3 -> CCS2005
```

The analytic pipeline string is printed by `fit_analytic_lcc.py` and saved as
`proj_pipeline` in `output/lcc_affine_fit_k{n}.json`; it maps **geodetic lon/lat
in BGS2005 (EPSG:7798)** to KS1970 E/N. `transform_points.py` prepends the fully
known inverse of EPSG:7801 so it can be driven directly from CCS2005 E/N — see
`build_analytic_pipeline()` there for the composed string.

### Via the `cct` command-line tool

`cct` (part of `proj-bin`) applies an arbitrary PROJ pipeline to coordinates from
stdin. It expects 4 columns (X Y Z T) even for a purely 2D shift:

```bash
echo "300000 4800000 0 0" | cct +proj=pipeline +step +proj=tinshift +file=output/bg_k3_tinshift.json
echo "9500000 4600000 0 0" | cct +proj=pipeline +step +proj=tinshift +file=output/bg_k5_tinshift.json +inv
```

## 4. Using this in QGIS

Two ways, trading convenience for accuracy.

### Recommended for cadastral-grade work: pre-reproject the file

Run the grid-accurate transform once, up front, and load the *result* — no
custom CRS or coordinate-operation setup needed in QGIS at all, since the output
coordinates are already correct:

```bash
python3 scripts/reproject_shapefile.py --zone k3 --direction forward \
    --in source_data/bg_k3_1km_sel.shp --out bg_k3_as_ks1970.gpkg
```

Works for points/lines/polygons, `.shp`/`.gpkg`/`.geojson` in or out. The output
carries no CRS metadata (the true destination CRS has no official definition to
attach), so QGIS will prompt for one or default to the project CRS. That is fine
for display — the coordinate *values* are already correct to the grid's
sub-millimetre accuracy.

Note: **`ogr2ogr -ct "<pipeline>"` does not work reliably for this.** It
internally disables the whole coordinate-transform object after the first few
points near the grid's domain edge fail, which then silently drops every feature
in the file. `gdaltransform` and `reproject_shapefile.py` (which drives pyproj
directly) do not have this problem.

### Convenient for general GIS use: import a Custom CRS

`output/bg_k{n}_ks1970_fitted.wkt` is a self-contained CRS definition embedding
the entire fitted analytic pipeline as one PROJ-based WKT2 `DERIVEDPROJCRS`, so
QGIS can transform to and from it on the fly like any ordinary CRS:

1. **Settings ▸ Custom Projections…**
2. Click **+**, give it a name (e.g. `BG 1970 K3 (fitted)`)
3. Paste the contents of `output/bg_k3_ks1970_fitted.wkt` into the definition
   field (QGIS accepts WKT here, not just Proj strings)
4. It is now selectable anywhere QGIS asks for a CRS

This is the **analytic model**, not the grid — expect ~0.1 m, not
sub-millimetre. Good for visualization and general GIS work; use the
pre-reprojected file above when the numbers themselves must be trustworthy.

### Command-line tools

```bash
# cct (PROJ) -- needs all 4 columns X Y Z T even for a 2D shift
echo "300000 4800000 0 0" | cct +proj=pipeline +step +proj=tinshift +file=output/bg_k3_tinshift.json

# gdaltransform (GDAL) -- equivalent, only needs X Y
echo "300000 4800000" | gdaltransform -ct "+proj=pipeline +step +proj=tinshift +file=output/bg_k3_tinshift.json"
```

## 5. Caveats

- Neither artifact is the official AGKK transformation. Each is an empirical
  reconstruction from that zone's point sample; no authoritative parameters or
  BGS2005 ↔ Pulkovo 1942(58) datum-shift constants are publicly published.
- The conic parameters in `output/lcc_affine_fit_k{n}.json` are a **convention,
  not a recovered definition** — see finding 1 above. The cone constant *n* is
  the part that is actually determined by the data. Do not quote the tangent
  parallel or λ₀ as though they were AGKK's.
- The grid only interpolates inside the sampled territory of its zone. Outside
  it PROJ returns `inf`/`inf`; `transform_points.py` falls back to the analytic
  model. Because the zones tile rather than overlap, a point just outside one
  zone's coverage will usually not be covered by a neighbour's file either.
- The PROJ triangulation-file schema has an optional `fallback_strategy` field
  meant to handle that out-of-domain case gracefully, but some PROJ versions
  reject the whole file if it is present at all (any value, `none` included).
  It is deliberately omitted; re-test with your own PROJ if you want it.
- The destination datum is assumed to be Krasovsky 1940 / Pulkovo 1942(58)
  throughout, per every `_sel_tr.prj`. Nothing here needs to model the datum
  shift explicitly: both methods work from point correspondences, not from
  geodetic datum parameters.
- The `.wkt` files parse and transform correctly with pyproj and `projinfo`, but
  `projinfo` prints a non-fatal `unexpected CS, expecting ID` warning: the
  derived CRS has no EPSG authority code (correctly — it is not registered) and
  strict WKT2 grammar wants one. Harmless.
- The 24°E discontinuity (finding 2) is inherited from AGKK's tool. Survey work
  crossing that meridian in K3 or K9 should expect ~0.1 m of disagreement there
  regardless of which method is used.
