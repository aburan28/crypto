"""Confirmatory run for PROTOCOL.md (T1, T2).  Usage: python3 run_confirm.py results SEED"""
import json
import random
import sys

import count
import tri

TRIALS, PASS, MAXCOLS = 5, 4, 12000


def cols_of(n, a1, l, b, caps):
    return count.counts(n, a1, l, b, caps)[1]


def main(out_dir, seed):
    rng = random.Random(seed)
    part_a, part_b, raw = [], [], []
    # Part A: smallest t-degree that resolves a1 bits of X1, b = 1
    for n in (11, 13, 17):
        for l in (3, 4):
            for dA in (1, 2):
                for a1 in range(1, l + 1):
                    capsf = lambda dT, a1=a1, dA=dA: (min(dA, a1), 1, 1, dT)
                    pred = next((dT for dT in range(1, 8) if count.predicted(n, a1, l, 1, capsf(dT))), None)
                    got, censored_at = None, None
                    for dT in range(1, 4):
                        caps = capsf(dT)
                        if cols_of(n, a1, l, 1, caps) > MAXCOLS:
                            censored_at = dT
                            break
                        r = tri.cell(n, a1, l, 1, caps, TRIALS, rng)
                        raw.append(r)
                        if r["resolved"] >= PASS:
                            got = dT
                            break
                    row = dict(n=n, l=l, dA=dA, a1=a1, count_min_dT=pred, measured_min_dT=got,
                               censored_at_dT=censored_at)
                    part_a.append(row)
                    print("A", json.dumps(row), flush=True)
    # Part B: largest b of X3 resolved with t unknown, a1 = 1, caps (1,1,1,2)
    for n in (11, 13, 17):
        for l in (3, 4, 5, 6):
            best, tested = 0, []
            for b in range(1, l + 1):
                caps = (1, 1, 1, 2)
                if cols_of(n, 1, l, b, caps) > MAXCOLS:
                    break
                r = tri.cell(n, 1, l, b, caps, TRIALS, rng)
                raw.append(r)
                tested.append(b)
                if r["resolved"] >= PASS:
                    best = b
            pred = max([b for b in range(0, l + 1) if count.predicted(n, 1, l, b, (1, 1, 1, 2))], default=0)
            row = dict(n=n, l=l, count_bmax=pred, measured_bmax=best, tested=tested)
            part_b.append(row)
            print("B", json.dumps(row), flush=True)
    with open(f"{out_dir}/confirm.json", "w") as fh:
        json.dump(dict(seed=seed, part_a=part_a, part_b=part_b, raw=raw), fh, indent=1)


if __name__ == "__main__":
    main(sys.argv[1], int(sys.argv[2]))
