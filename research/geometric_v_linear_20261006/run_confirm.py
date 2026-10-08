"""Confirmatory run for PROTOCOL.md (G1-G3).  Usage: python3 run_confirm.py OUTDIR SEED"""
import json
import random
import sys

import glin

PLANTED = 30


def trials_for(n, l, kind):
    want = 24 * 2 ** max(0, n - 2 * l + 1)          # ~24 expected successes
    return int(min(want, 2 ** 15 if kind == "geometric" else 2 ** 13))


def main(out_dir, seed):
    rng = random.Random(seed)
    cells = []
    for n in (19, 23, 29, 31):
        for l in range(4, (n + 1) // 3 + 3):
            for kind in ("geometric", "random"):
                r = glin.cell(n, l, kind, trials_for(n, l, kind), PLANTED, rng)
                cells.append(r)
                print(json.dumps(r), flush=True)
    with open(f"{out_dir}/confirm.json", "w") as fh:
        json.dump(dict(seed=seed, cells=cells), fh, indent=1)


if __name__ == "__main__":
    main(sys.argv[1], int(sys.argv[2]))
