"""Confirmatory run for PROTOCOL.md (X1, X2).  Usage: python3 run_confirm.py results SEED"""
import json
import random
import sys
from math import comb

import xl

GRID_N = (19, 23, 29, 31)
BIDEGREES = ((2, 1), (3, 1), (1, 2), (2, 2))
TRIALS = 20
PASS = 18          # a b is resolved when >= 18 of 20 planted instances resolve


def Ms(k, D):
    return sum(comb(k, i) for i in range(D + 1)) if D >= 0 else 0


def tight(n, l, D1, D2):
    ok = [b for b in range(0, l + 1)
          if n * Ms(l, D1 - 1) * Ms(b, D2 - 1) >= Ms(l, D1) * Ms(b, D2) - 1]
    return max(ok) if ok else 0


def main(out_dir, seed):
    rng = random.Random(seed)
    cells, raw = [], []
    for n in GRID_N:
        for l in range(6, n // 2 + 2, 2):
            for D1, D2 in BIDEGREES:
                pred = tight(n, l, D1, D2)
                bs = range(max(1, pred - 3), min(l, pred + 3) + 1)
                ok = {}
                for b in bs:
                    r = xl.xl_cell(n, l, b, D1, D2, TRIALS, rng)
                    raw.append(r)
                    ok[b] = r["resolved"] >= PASS
                    if r["inconsistent"]:
                        print("INCONSISTENT on planted:", r, flush=True)
                resolved = [b for b in bs if ok[b]]
                bmax = max(resolved) if resolved else 0
                cell = dict(n=n, l=l, D1=D1, D2=D2, tight=pred, bmax=bmax,
                            tested=list(bs), at_cap=(bmax == l))
                cells.append(cell)
                print(json.dumps(cell), flush=True)
    with open(f"{out_dir}/confirm.json", "w") as fh:
        json.dump(dict(seed=seed, cells=cells, raw=raw), fh, indent=1)


if __name__ == "__main__":
    main(sys.argv[1], int(sys.argv[2]))
