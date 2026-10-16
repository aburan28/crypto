#!/usr/bin/env python3
from runbench import ROOT, bench

arms = [
    ("skipempty", ROOT / "build/g7-xonly-mod72-skipempty", 524288, 16),
    ("control", ROOT / "build/g7-xonly-mod72-control", 524288, 16),
]
for rep in range(3):
    order = arms if rep % 2 == 0 else list(reversed(arms))
    for label, binary, workers, batch in order:
        bench(binary, f"skipempty-r{rep}-{label}", workers=workers, batch=batch, run_id=62580)
