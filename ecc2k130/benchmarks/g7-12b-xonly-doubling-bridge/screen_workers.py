#!/usr/bin/env python3
from runbench import ROOT, bench

binary = ROOT / "build/g7-xonly-sparse-bridge-block"
arms = [
    ("w131", 131072, 16),
    ("w262", 262144, 8),
    ("w524", 524288, 4),
    ("w1048", 1048576, 2),
]
for rep in range(3):
    order = arms if rep % 2 == 0 else list(reversed(arms))
    for label, workers, launches in order:
        bench(binary, f"workers-r{rep}-{label}", workers=workers,
              steps=1024, launches=launches, batch=16, run_id=62573)
