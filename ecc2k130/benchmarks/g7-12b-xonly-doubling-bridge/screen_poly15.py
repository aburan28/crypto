#!/usr/bin/env python3
from runbench import ROOT, bench

arms = [
    ("poly12", ROOT / "build/g7-xonly-15b-poly12", 524288, 16),
    ("arith", ROOT / "build/g7-xonly-15b-arith", 524288, 16),
    ("control", ROOT / "build/g7-xonly-15b-control", 524288, 16),
]
for rep in range(3):
    order = arms if rep % 2 == 0 else list(reversed(arms))
    for label, binary, workers, batch in order:
        bench(binary, f"poly15-r{rep}-{label}", workers=workers, batch=batch, run_id=62590)
