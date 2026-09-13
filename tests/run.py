#!/usr/bin/env python3
"""Run the test suite without pytest. `python tests/run.py [pattern]`."""
from __future__ import annotations

import importlib.util
import sys
import traceback
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE))


def load(path: Path):
    spec = importlib.util.spec_from_file_location(path.stem, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def main(argv: list[str]) -> int:
    pattern = argv[0] if argv else ""
    passed = failed = 0
    failures = []
    for path in sorted(HERE.glob("test_*.py")):
        mod = load(path)
        for name in sorted(dir(mod)):
            if not name.startswith("test_") or pattern not in f"{path.stem}.{name}":
                continue
            fn = getattr(mod, name)
            if not callable(fn):
                continue
            try:
                fn()
                passed += 1
                print(f"  ok   {path.stem}.{name}")
            except Exception:
                failed += 1
                failures.append((f"{path.stem}.{name}", traceback.format_exc()))
                print(f"  FAIL {path.stem}.{name}")
    for name, tb in failures:
        print(f"\n===== {name} =====\n{tb}")
    print(f"\n{passed} passed, {failed} failed")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
