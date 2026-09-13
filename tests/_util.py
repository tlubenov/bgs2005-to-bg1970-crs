"""
Tiny helpers so the suite is valid pytest but needs nothing installed to run.

`pip install -e .[dev]` gets you pytest and `python -m pytest` works as usual;
without it, `python tests/run.py` runs the same functions. That matters here
because the environment this project is developed in has GDAL and PROJ but no
pip, and a test suite you cannot run is not a test suite.
"""
from __future__ import annotations

from contextlib import contextmanager


@contextmanager
def raises(exc_type):
    """Minimal stand-in for pytest.raises."""
    try:
        yield
    except exc_type:
        return
    except Exception as other:  # noqa: BLE001
        raise AssertionError(f"expected {exc_type.__name__}, got {type(other).__name__}: {other}")
    raise AssertionError(f"expected {exc_type.__name__}, nothing raised")
