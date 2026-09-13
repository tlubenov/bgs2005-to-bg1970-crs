"""
Where the project's data, outputs and paper live.

The package is installable, so it cannot assume it is being run from inside the
repository. Every location is resolved at call time from, in order: an explicit
argument, an environment variable, the current directory, and finally the
repository the package was installed from (which is the common case for an
editable install). Nothing is captured at import time, so a caller can chdir or
set the environment after importing and still get the right answer.
"""
from __future__ import annotations

import os
from pathlib import Path

#: Datasets the project ships, in the order the pipeline prefers them.
DATASET_NAMES = ("source_data_buff20km", "source_data", "source_csv")

ENV_ROOT = "BGKS1970_ROOT"
ENV_OUTPUT = "BGKS1970_OUTPUT"
ENV_PAPER = "BGKS1970_PAPER"

#: The repo this package was installed from, if the layout still looks like it.
_PACKAGE_PARENT = Path(__file__).resolve().parent.parent


def _looks_like_root(path: Path) -> bool:
    """A project root is any directory holding the primary dataset, in either
    form -- the committed `.zip` or an unpacked directory."""
    if not path or not path.is_dir():
        return False
    primary = DATASET_NAMES[0]
    return (path / f"{primary}.zip").exists() or (path / primary).is_dir()


def project_root(explicit: str | Path | None = None) -> Path:
    """Directory holding the dataset archives, bg_ext/ and the paper."""
    candidates: list[Path] = []
    if explicit:
        candidates.append(Path(explicit))
    if os.environ.get(ENV_ROOT):
        candidates.append(Path(os.environ[ENV_ROOT]))
    candidates += [Path.cwd(), _PACKAGE_PARENT]

    for cand in candidates:
        cand = cand.expanduser()
        if _looks_like_root(cand):
            return cand.resolve()

    tried = ", ".join(str(c) for c in candidates)
    raise FileNotFoundError(
        f"Could not find the project data. Looked for {DATASET_NAMES[0]}.zip (or an "
        f"unpacked {DATASET_NAMES[0]}/) in: {tried}. Run from the project directory, "
        f"or set {ENV_ROOT}."
    )


def output_dir(root: str | Path | None = None) -> Path:
    """Where fitted parameters and grids are written."""
    if os.environ.get(ENV_OUTPUT):
        path = Path(os.environ[ENV_OUTPUT]).expanduser()
    else:
        path = project_root(root) / "output"
    path.mkdir(parents=True, exist_ok=True)
    return path


def paper_dir(root: str | Path | None = None) -> Path:
    """Where the paper, its figures and its PDF editions live."""
    if os.environ.get(ENV_PAPER):
        path = Path(os.environ[ENV_PAPER]).expanduser()
    else:
        path = project_root(root) / "paper"
    path.mkdir(parents=True, exist_ok=True)
    return path


def ext_dir(root: str | Path | None = None) -> Path:
    """bg_ext/ -- the zone extents. Small enough that it is not archived."""
    return project_root(root) / "bg_ext"
