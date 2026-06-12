"""Test bootstrap: make scripts/ and the repo root importable."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
for p in (ROOT, ROOT / "scripts", ROOT / "analysis"):
    sp = str(p)
    if sp not in sys.path:
        sys.path.insert(0, sp)
