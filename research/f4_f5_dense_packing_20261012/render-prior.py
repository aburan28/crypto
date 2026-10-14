#!/usr/bin/env python3
"""Re-render the frozen prior report, adding PNGs from the same figures."""
import argparse
from pathlib import Path
import runpy
import sys
import tempfile
import os

parser = argparse.ArgumentParser()
parser.add_argument("--out", required=True, type=Path)
args = parser.parse_args()
with tempfile.TemporaryDirectory(prefix="f4-f5-prior-render-") as cache:
    os.environ.setdefault("MPLCONFIGDIR", cache)
    os.environ.setdefault("XDG_CACHE_HOME", cache)
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib.figure import Figure

    original = Figure.savefig

    def with_png(self, filename, *positional, **kwargs):
        original(self, filename, *positional, **kwargs)
        if Path(filename).suffix == ".svg":
            original(self, Path(filename).with_suffix(".png"), dpi=150,
                     metadata={"Software": "matplotlib"})

    Figure.savefig = with_png
    source = Path(__file__).resolve().parent.parent / "f4_f5_fused_readback_20261012/render.py"
    sys.argv = [str(source), "--out", str(args.out)]
    runpy.run_path(str(source), run_name="__main__")
