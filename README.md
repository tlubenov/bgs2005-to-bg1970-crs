# bgks1970 — recovering Bulgaria's unpublished "1970 system" zones from BGS2005

[![status](https://img.shields.io/badge/status-stable-brightgreen)](#)
An installable Python library and CLI that recovers Bulgaria's KS1970 Lambert
zones (K3/K5/K7/K9) from BGS2005/CCS2005, and rebuilds every artifact in this
repository with **one command**:

```bash
pip install -e .
bgks1970 reconstruct
```


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

> **No unpacking needed.** The sampled data ships as zip archives — the raw trees
> are ~87 MB against ~14 MB zipped — and GDAL reads the shapefiles straight out of
> them through its `/vsizip/` virtual filesystem. If you do unpack an archive, the
> unpacked directory is used in preference; either way the code addresses it the
> same. Extracted directories are gitignored.

`source_data_buff20km/` is the primary dataset. For each zone `n` ∈ {3,5,7,9} and
each sampling step `r` ∈ {1km, 2km}:

| file | contents |
|---|---|
| `bg_k{n}_{r}_sel.shp` | points in **EPSG:7801** (`BGS2005 / CCS2005`), a fully documented Lambert Conformal Conic (2SP) on GRS80: central meridian 25.5°E, standard parallels 42°/43°20′, latitude of false origin 42.6678756833333°, false easting 500 000 m, false northing 4 725 824.3591 m — identical across all four zones |
| `bg_k{n}_{r}_sel_tr.shp` | the same points run through AGKK's official tool into KS1970 zone K{n}. Its `.prj` gives the datum and projection type only |
| `bg_k{n}_{r}_sel_tr.log` | AGKK's own report: *"Target zone width: 3 degrees … Plane transform accuracy - 0.14 meters"* |

Each `_sel` / `_sel_tr` pair shares an integer `id`, so every point is an exact
correspondence — no manually surveyed control points are involved. The `id` is a
per-file row key only; the join is always made within one pair, so the fact that a
zone's 1 km and 2 km files reuse some id values is harmless.

Sampling covers each zone's extent **buffered 20 km outward**. That matters in
practice: a grid sampled exactly to the zone boundary stops working precisely where
it is most often needed — along the boundaries between zones, and for features that
cross them. Measured on K9's own official boundary vertices, the un-buffered grid
covered **9 of 126** of them (7%); the buffered grid covers **all 126**.

| zone | quadrant | 1 km points | 2 km points | zone area | 
|---|---|---|---|---|
| K3 | northwest | 47 804 | 11 939 | 23 695 km² |
| K5 | southeast | 56 594 | 14 159 | 31 891 km² |
| K7 | northeast | 53 912 | 13 478 | 30 378 km² |
| K9 | southwest | 62 497 | 15 620 | 37 333 km² |

Supporting data:

- `source_data/` — the earlier, un-buffered sampling. Kept because it is a
  complete, independent set of positions, which makes it a usable fall-back
  validation set.
- `source_csv/points_{1,2}km_bgs7801.*` — the master national lattice in EPSG:7801
  (225 216 points at 1 km, 56 448 at 2 km) that the per-zone `_sel` selections were
  cut from.
- `bg_ext/bg_zones1970.shp` — the four official zone extents (field `CLIST` =
  3/5/7/9) in WGS 84 / UTM zone 35N.
- `bg_ext/bg_zones1970_buff20km.geojson` — the same, buffered by 20 km, already in
  EPSG:7801. This is the region the primary dataset samples; its union runs
  E 206 222…780 705, N 4 539 859…4 928 413.

The four official zone polygons **tile** the country rather than overlapping: their
areas sum to 123 297 km² against a union of 123 281 km².

### A note on K5 — resolved

K5's first 1 km export was interrupted. BGSTrans writes the transformed shapefile
last, so it left a `.dbf` and `.shx` behind with no `.shp`, a `.prj` that was never
written, and a 32-byte `.log` against 510 for every other zone — a state that looks
complete to a casual `ls`.

The orphaned `.dbf` was still readable and showed the run had covered 56 633 of
56 636 records, with the three missing ids being the *last three in file order* — a
clean truncation at the end of the run rather than three points the transform
rejected. (A tempting theory, that AGKK's domain stops at 41°N, was disproved by the
data: points further south in the same row transformed fine.) The export has since
been re-run and all four zones now come from complete data.

`bgks1970 check` exists to make that failure mode impossible to miss — and
`bgks1970 reconstruct` runs it first, refusing to fit truncated data.

### Why the two sampling steps matter

The 1 km and 2 km lattices of a zone are offset by exactly **500 m**, so a 2 km
point sits at the *centre* of a 1 km cell and no 2 km point coincides with any 1 km
point (verified: zero coordinate collisions in all four zones).

So every model is built on the 1 km set and scored on the 2 km set, which is a
genuinely independent sample rather than a random split of one lattice. For the
interpolating model that distinction matters a great deal: a random hold-out of a
single lattice leaves test points sitting on the *edges* of triangles built from
their own immediate neighbours, which flatters the result. Cell centres are the
worst possible position for piecewise-linear interpolation, so every number below is
an upper bound rather than a favourable sample.

## Two solutions per zone

| | `bgks1970 fit` | `bgks1970 grid` |
|---|---|---|
| Output | closed-form PROJ pipeline (LCC + 2D similarity) | `tinshift` triangulated grid (JSON) |
| Coverage | anywhere (extrapolates) | only within the zone's sampled territory |
| Portability | one short PROJ string, embeddable anywhere | one 2–4 MB JSON file PROJ must load |
| Accuracy | 0.14–0.21 m RMSE | 0.14–0.61 **mm** median |
| Use when | you need a compact formula, or points outside the sampled area | you want survey-grade accuracy (recommended) |

### Closed-form model — out-of-sample accuracy

| zone | cone constant *n* | tangent parallel | rotation | scale | RMSE | p95 | max |
|---|---|---|---|---|---|---|---|
| K3 | 0.687723149 | 43.450145° | +5704.9″ | −6.55 ppm | 0.172 m | 0.329 m | 0.445 m |
| K5 | 0.675170794 | 42.467415° | −2161.3″ | −3.05 ppm | 0.213 m | 0.401 m | 0.641 m |
| K7 | 0.689029005 | 43.553295° | −1812.4″ | +1.41 ppm | 0.190 m | 0.352 m | 0.494 m |
| K9 | 0.672869214 | 42.288900° | +5041.6″ | −8.81 ppm | 0.139 m | 0.282 m | 0.497 m |

These are of the order of AGKK's own declared 0.14 m plane accuracy, ranging from
0.14 to 0.21 m by zone. Fit-set and out-of-sample RMSE agree to three decimals in
every zone, so this residual is model bias, not overfitting.

Against the earlier un-buffered sampling the same **cone constants reproduce to five
decimals** (K3 0.687723 vs 0.687702, K9 0.672869 vs 0.672865) while the RMSE rises —
as expected, since one global conic now has to cover 20 km more territory in every
direction. That cross-dataset agreement is the strongest evidence that *n* is a real
recovered property rather than a fitting artefact.

### Grid model — out-of-sample accuracy

| zone | vertices | triangles | size | covered | median | p95 | max | residuals > 1 mm |
|---|---|---|---|---|---|---|---|---|
| K3 | 47 804 | 94 569 | 4.54 MB | 11 810/11 939 | 0.535 mm | 0.544 mm | 79.3 mm | 36 |
| K5 | 56 594 | 112 120 | 5.39 MB | 14 025/14 159 | 0.137 mm | 0.139 mm | 0.139 mm | 0 |
| K7 | 53 912 | 106 798 | 5.13 MB | 13 343/13 478 | 0.607 mm | 0.615 mm | 0.616 mm | 0 |
| K9 | 62 497 | 123 915 | 5.96 MB | 15 497/15 620 | 0.262 mm | 0.266 mm | 182.3 mm | 61 |

Uncovered points fall just outside the built lattice's outline; they are excluded
rather than extrapolated. Round-trip error, forward then inverse, is under 3 nm (one
K9 point on the domain edge does not invert and is reported separately rather than
being allowed to poison the statistic).

Medians are higher than the previous un-buffered run (0.055–0.240 mm) for a good
reason: the 2 km points now sit at exact cell centres instead of 120 m off-node, so
this is the true worst case rather than a lucky sample.

Before trusting a fresh export, run `bgks1970 check` — it exits non-zero on an
incomplete one.

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

Measured, on K3 (`bgks1970 fit --identifiability --zone k3`):

| start λ₀ | fitted λ₀ | φ₁ | φ₂ | scale−1 | cone *n* | RMSE |
|---|---|---|---|---|---|---|
| 23.0° | 23.119220° | 43.347° | 43.554° | −4.9 ppm | 0.687723149 | 0.1720 m |
| 25.5° | 24.353738° | 50.062° | 36.583° | +6936.3 ppm | 0.687723149 | 0.1720 m |
| 27.5° | 25.341094° | 47.275° | 39.541° | +2268.2 ppm | 0.687723149 | 0.1720 m |

Parameters differing by **degrees** and by nearly 7000 ppm of scale produce
coordinates agreeing to **under one micrometre**, and agree on *n* to nine
significant figures. So this repository fits the identifiable 5-parameter form
directly — *n*, θ, *s*, tₓ, t_y, with λ₀ pinned to 25.5° as a stated convention —
and reports *n* as the recovered quantity. Same accuracy, no flat directions,
and nothing claimed that the data does not determine.

### 2. AGKK's own transformation is discontinuous along 24°E

Along any straight row of the lattice, a smooth projection produces a second
difference of a couple of millimetres per km² — that is just the conic's own
curvature. In K3 and K9 there is a localized spike two orders of magnitude higher,
sitting exactly on the 24°E meridian:

| zone | seam longitude | peak | typical curvature | residuals > 1 mm |
|---|---|---|---|---|
| K3 | 24.0000° | 130 mm/km² | 2.2 mm/km² | 36 |
| K9 | 24.0000° | 150 mm/km² | 1.0 mm/km² | 61 |

**How this is aggregated is critical.** The kink sits on a *meridian*, and a meridian
is not a column of this projection's grid: over a zone's ~250 km of latitude the
easting of the 24°E meridian drifts by about 4 km. Bin the curvature by easting
column and the kink is smeared across eight columns, present in under half the rows
of each — and a median over rows hides it completely. Bin by longitude and every
row's kink lands in the same bin, pinning the seam to 24.0000°.

This is a property of AGKK's transformation, not of either dataset or of the
triangulation here — building the grid from one set and scoring the other reproduces
it at the same longitude, and the earlier independent sampling shows it in the same
place. It is almost certainly a panel boundary in the tool's own internal correction
grid. K5 and K7 show nothing comparable.

The consequence: **all 97 out-of-sample residuals above 1 mm in the entire project
lie on this line.** The grid reproduces the kink exactly at sampled points but must
interpolate straight across it in the one cell that straddles it. Excluding those
points, every zone's out-of-sample RMSE is its median: 0.535, 0.137, 0.607 and
0.262 mm. No closed-form model can represent the kink at all. `bgks1970 grid`
detects and reports this automatically.

## Layout

```
bgks1970/                  the library
  paths.py            locating the project's data, output and paper directories
  datasets.py         reading layers out of a zip (/vsizip/) or an unpacked dir
  data.py             zones, the fit/test split, matched points, zone extents
  lcc.py              ellipsoidal LCC (Snyder 1987): 2SP and cone-constant forms
  optim.py            Levenberg-Marquardt least squares, in numpy
  lattice.py          triangulation of a regular ragged-edged lattice
  check.py            audits the input archives
  fit.py              fits the identifiable closed-form model per zone
  grid.py             builds the tinshift grids, detects AGKK seams
  crs.py              QGIS-importable Custom CRS per zone
  transform.py        point transforms, forward and reverse, either model
  reproject.py        reprojects a whole vector file
  figures.py          the paper's figures, generated and injected
  paper.py            renders the paper to PDF, one file per language
  pipeline.py         `reconstruct()` -- every stage, in dependency order
  cli.py              the `bgks1970` command
tests/                unit + integration tests (pytest, or `python tests/run.py`)
source_data_buff20km.zip   PRIMARY input data -- read in place, no unpacking
source_data.zip            the earlier un-buffered sampling
source_csv.zip             the master national 1 km / 2 km lattice
bg_ext/                    zone extents: official, and buffered by 20 km
output/
  lcc_affine_fit_k{n}.json      fitted parameters + accuracy report, per zone
  bg_k{n}_tinshift.json         the PROJ tinshift grid (both directions)
  bg_k{n}_ks1970_fitted.wkt     QGIS Custom CRS wrapping the analytic pipeline
  tinshift_accuracy.json        grid accuracy + seam diagnostics, all zones
paper/
  recovering_ks1970.html        the write-up (Bulgarian + English)
  recovering_ks1970_{bg,en}.pdf A4 editions
  figures/{zones,seam}.svg      generated from bg_ext and the data
pyproject.toml
```

## 1. Installing PROJ and its Python bindings

`tinshift` requires **PROJ ≥ 7.2.0**. Check what you have:

```bash
proj    # first line prints the release, e.g. "Rel. 9.7.1"
```

Dependencies are deliberately minimal: **numpy and pyproj**. There is no SciPy
dependency — the Levenberg-Marquardt fit and the triangulation are implemented
directly in `bgks1970/optim.py` and `bgks1970/lattice.py`, both covered by the
test suite. GDAL is needed only to read the shapefiles.

### Option A — system packages (Debian/Ubuntu, what this repo was built against)

```bash
sudo apt-get install proj-bin libproj-dev gdal-bin libgdal-dev python3-gdal
pip install -e .
```

`proj-bin` gives you the `proj`/`cs2cs`/`cct` CLI tools; `libproj-dev` is needed
if anything compiles against PROJ; `gdal-bin`/`libgdal-dev`/`python3-gdal` are
needed only because the library reads shapefiles through GDAL's OGR Python
bindings. If you only ever use the already-built `output/*.json` files,
you can skip GDAL entirely.

### Option B — conda/mamba (easiest way to get matching GDAL+PROJ+pyproj)

```bash
conda create -n bgcrs -c conda-forge python=3.12 gdal proj pyproj numpy
conda activate bgcrs
pip install -e .
```

This sidesteps the most common PROJ headache: GDAL's Python bindings and
`pyproj` each embed their own copy of PROJ, and if those versions disagree about
where the PROJ resource data lives you get obscure `proj_create` failures.
conda-forge keeps them in lockstep.

### Option C — macOS (Homebrew)

```bash
brew install proj gdal
pip install -e .
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

## 2. Installing and reproducing everything

```bash
pip install -e .          # or: pip install -e ".[dev]" for the test suite
```

GDAL is deliberately not a dependency — pip-installing it demands an exact match
with your system libgdal, which usually fails. Get it from your OS package manager
(`apt install python3-gdal`) or conda-forge, as above.

### One command

```bash
bgks1970 reconstruct
```

That audits the input archives, fits the closed-form model, builds the grids,
writes the QGIS CRS definitions, regenerates the paper's figures and renders its
PDF editions — about a minute in total, reading straight from the zips. It stops
at the audit if an export is incomplete, because every number downstream would
otherwise be quietly built on truncated data.

```
======================================================================
== reconstruction summary
======================================================================
  ok    audit        6.9s  all layers complete and consistent
  ok    fit          7.5s  4 zones fitted
  ok    grid         9.5s  4 grids built
  ok    crs          1.8s  WKT written
  ok    figures     33.8s  figures written and injected
  ok    paper        1.5s  2 PDF(s)
        total       60.9s
```

Useful variants:

```bash
bgks1970 reconstruct --zone k5              # one zone
bgks1970 reconstruct --stages fit,grid      # a subset of stages
bgks1970 reconstruct --keep-going           # don't stop at a failing stage
                                            # (e.g. no browser for the PDF stage)
```

### Individual stages

```bash
bgks1970 check                    # audit the input archives; non-zero if unhappy
bgks1970 fit                      # -> output/lcc_affine_fit_k{n}.json
bgks1970 fit --identifiability --zone k3    # reproduce finding 1 from the data
bgks1970 grid                     # -> output/bg_k{n}_tinshift.json
bgks1970 crs                      # -> output/bg_k{n}_ks1970_fitted.wkt
bgks1970 figures                  # -> paper/figures/*.svg, injected into the paper
bgks1970 paper --lang en          # -> paper/recovering_ks1970_en.pdf
```

### Tests

```bash
python -m pytest            # with pytest installed
python tests/run.py         # same tests, no dependencies
```

### Pointing it at other data

The project directory is found from the current directory, then from the installed
repository. Override it, or the output and paper locations, with `BGKS1970_ROOT`,
`BGKS1970_OUTPUT` and `BGKS1970_PAPER`.

## 3. Using the results

### As a library

```python
import bgks1970 as bg

E, N = bg.transform([300000], [4800000], "k3")          # CCS2005 -> KS1970 K3
x, y = bg.transform(E, N, "k3", direction="reverse")    # and back

bg.reconstruct()                                        # rebuild every artifact
```

`transform()` uses the grid by default and retries any point outside its coverage
with the analytic model; pass `method="analytic"` to force the formula, or
`fallback=False` to get `inf` instead of a fallback.

### Via the CLI

```bash
# CCS2005 -> KS1970 K3, single point
bgks1970 transform --zone k3 --xy 300000 4800000

# KS1970 K5 -> CCS2005 (reverse)
bgks1970 transform --zone k5 --direction reverse --xy 9500000 4600000

# batch, CSV with columns x,y (any extra columns are preserved)
bgks1970 transform --zone k9 --csv in.csv --out out.csv

# force the closed-form model (e.g. for a point far outside that zone)
bgks1970 transform --zone k7 --method analytic --xy 500000 4800000
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

The analytic pipeline string is printed by `bgks1970 fit` and saved as
`proj_pipeline` in `output/lcc_affine_fit_k{n}.json`; it maps **geodetic lon/lat in
BGS2005 (EPSG:7798)** to KS1970 E/N. `bgks1970.transform.build_analytic_pipeline()`
prepends the fully known inverse of EPSG:7801 so it can be driven directly from
CCS2005 E/N.

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
bgks1970 reproject --zone k3 --direction forward \
    --in my_layer.shp --out my_layer_ks1970.gpkg
```

It reads zipped input too, e.g.
`--in /vsizip/source_data_buff20km.zip/source_data_buff20km/bg_k3_1km_sel.shp`.

Works for points/lines/polygons, `.shp`/`.gpkg`/`.geojson` in or out. The output
carries no CRS metadata (the true destination CRS has no official definition to
attach), so QGIS will prompt for one or default to the project CRS. That is fine
for display — the coordinate *values* are already correct to the grid's
sub-millimetre accuracy.

Note: **`ogr2ogr -ct "<pipeline>"` does not work reliably for this.** It
internally disables the whole coordinate-transform object after the first few
points near the grid's domain edge fail, which then silently drops every feature
in the file. `gdaltransform` and `bgks1970 reproject` (which drives pyproj directly) do not
have this problem.

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
  it PROJ returns `inf`/`inf`; `bgks1970 transform` falls back to the analytic
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
- The pipeline reads `source_data_buff20km.zip` by default; set `BGKS1970_ROOT` to
  point it at another project directory. Unpacked directories take precedence over
  archives when both are present.
