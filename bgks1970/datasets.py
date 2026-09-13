"""
Reading the sampled shapefiles, whether they are zipped or unpacked.

The datasets are committed as zip archives (~14 MB against ~87 MB unpacked), and
GDAL reads shapefiles straight out of a zip through its `/vsizip/` virtual
filesystem, so nothing has to be extracted to run the pipeline. An unpacked
directory is still accepted and preferred when present, which keeps a working
copy usable and makes the archive an implementation detail rather than a
requirement.

One `Dataset` wraps that choice: ask it for a layer path and it hands back
whatever OGR needs -- a plain path, or a `/vsizip/` one.
"""
from __future__ import annotations

import zipfile
from dataclasses import dataclass
from functools import cached_property
from pathlib import Path

from .paths import project_root

SHAPEFILE_PARTS = ("shp", "shx", "dbf", "prj")


@dataclass(frozen=True)
class Dataset:
    """One sampled dataset, e.g. `source_data_buff20km`, zipped or unpacked."""

    name: str
    root: Path

    @classmethod
    def find(cls, name: str, root: str | Path | None = None) -> "Dataset":
        return cls(name=name, root=project_root(root))

    # -- where it physically is ---------------------------------------------
    @property
    def directory(self) -> Path:
        return self.root / self.name

    @property
    def archive(self) -> Path:
        return self.root / f"{self.name}.zip"

    @property
    def unpacked(self) -> bool:
        return self.directory.is_dir()

    @property
    def exists(self) -> bool:
        return self.unpacked or self.archive.exists()

    @property
    def label(self) -> str:
        return self.name + ("" if self.unpacked else ".zip")

    # -- what is inside it ---------------------------------------------------
    @cached_property
    def _members(self) -> frozenset[str]:
        """Member names relative to the dataset directory.

        The archives store their own directory as a top-level prefix, so
        `source_data_buff20km.zip` contains `source_data_buff20km/bg_k3_...`.
        Stripping that prefix makes zipped and unpacked datasets addressable the
        same way."""
        if self.unpacked:
            return frozenset(p.name for p in self.directory.iterdir() if p.is_file())
        if not self.archive.exists():
            return frozenset()
        with zipfile.ZipFile(self.archive) as zf:
            names = set()
            for raw in zf.namelist():
                if raw.endswith("/"):
                    continue
                prefix = f"{self.name}/"
                names.add(raw[len(prefix):] if raw.startswith(prefix) else raw)
            return frozenset(names)

    def has(self, filename: str) -> bool:
        return filename in self._members

    def has_layer(self, stem: str) -> bool:
        """True only when every part of the shapefile is present. An interrupted
        BGSTrans export leaves the .dbf and .shx behind without the .shp."""
        return all(self.has(f"{stem}.{ext}") for ext in SHAPEFILE_PARTS)

    def missing_parts(self, stem: str) -> list[str]:
        return [ext for ext in SHAPEFILE_PARTS if not self.has(f"{stem}.{ext}")]

    # -- how to read it ------------------------------------------------------
    def path(self, filename: str) -> str:
        """An OGR-openable path for a member of this dataset."""
        if self.unpacked:
            return str(self.directory / filename)
        return f"/vsizip/{self.archive}/{self.name}/{filename}"

    def layer_path(self, stem: str) -> str:
        return self.path(f"{stem}.shp")

    def read_bytes(self, filename: str) -> bytes:
        if self.unpacked:
            return (self.directory / filename).read_bytes()
        with zipfile.ZipFile(self.archive) as zf:
            try:
                return zf.read(f"{self.name}/{filename}")
            except KeyError:
                return zf.read(filename)

    def read_text(self, filename: str, errors: str = "replace") -> str:
        return self.read_bytes(filename).decode("latin1", errors=errors)

    def require(self) -> "Dataset":
        if self.exists:
            return self
        raise FileNotFoundError(
            f"dataset {self.name!r} not found in {self.root}: expected "
            f"{self.archive.name} or an unpacked {self.name}/"
        )
