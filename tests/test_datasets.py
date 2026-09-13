"""Reading the sampled data out of a zip archive, or an unpacked directory."""
import zipfile

import numpy as np

from _util import raises
from bgks1970.data import PRIMARY_DATASET, ZONES, load_matched_points, resolve_split
from bgks1970.datasets import Dataset


def test_primary_dataset_is_addressable():
    d = Dataset.find(PRIMARY_DATASET)
    assert d.exists
    assert d.has_layer("bg_k3_1km_sel")
    assert d.layer_path("bg_k3_1km_sel").endswith("bg_k3_1km_sel.shp")


def test_zip_addressing_uses_vsizip_and_dir_addressing_does_not():
    d = Dataset.find(PRIMARY_DATASET)
    path = d.layer_path("bg_k3_1km_sel")
    assert path.startswith("/vsizip/") == (not d.unpacked)


def test_has_layer_requires_every_shapefile_part(tmp_path=None):
    """An interrupted BGSTrans export leaves .dbf/.shx without .shp -- that must
    not count as an available layer."""
    import tempfile, pathlib
    root = pathlib.Path(tempfile.mkdtemp())
    name = "fake_ds"
    zpath = root / f"{name}.zip"
    with zipfile.ZipFile(zpath, "w") as zf:
        for ext in ("dbf", "shx"):          # deliberately no .shp / .prj
            zf.writestr(f"{name}/bg_k5_1km_sel_tr.{ext}", b"x")
        for ext in ("shp", "shx", "dbf", "prj"):
            zf.writestr(f"{name}/bg_k5_1km_sel.{ext}", b"x")
    d = Dataset(name=name, root=root)
    assert d.has_layer("bg_k5_1km_sel")
    assert not d.has_layer("bg_k5_1km_sel_tr")
    assert sorted(d.missing_parts("bg_k5_1km_sel_tr")) == ["prj", "shp"]


def test_missing_dataset_reports_clearly():
    import tempfile, pathlib
    d = Dataset(name="nope", root=pathlib.Path(tempfile.mkdtemp()))
    with raises(FileNotFoundError):
        d.require()


def test_every_zone_resolves_to_a_complete_split():
    for zone in ZONES:
        split = resolve_split(zone)
        assert split.fit.step in ("1km", "2km")
        assert split.test, f"{zone} has no scoring set"


def test_matched_points_align_and_are_disjoint_between_steps():
    """The two lattices must not share positions -- that independence is what
    makes the 2 km set a real out-of-sample set."""
    ids, src, dst = load_matched_points("k5", "1km")
    assert len(ids) == len(src) == len(dst)
    assert np.isfinite(src).all() and np.isfinite(dst).all()

    _, src2, _ = load_matched_points("k5", "2km")
    a = {tuple(np.round(p, 3)) for p in src}
    b = {tuple(np.round(p, 3)) for p in src2}
    assert not (a & b), "1 km and 2 km lattices share positions"
