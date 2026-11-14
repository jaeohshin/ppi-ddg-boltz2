"""Make `mmcsm_bench` importable when pytest is run from the repository root.

This half of the repository is self-contained and has its own pyproject, so it is
not on the path by virtue of the root package being installed.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
