"""Make the package and the `_util` helpers importable under pytest, matching
what tests/run.py does for the no-dependency path."""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
for p in (HERE.parent, HERE):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))
