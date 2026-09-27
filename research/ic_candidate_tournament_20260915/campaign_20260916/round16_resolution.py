#!/usr/bin/env python3
"""Reproduce the two tables in ROUND16-single-target.md from frozen receipts.

    python3 campaign_20260916/round16_resolution.py [--round runs/round-0015]

Reads nothing but `runs/<round>/runs/confirmation/**/receipt.json` and re-runs
the same nested bootstrap the evaluator uses in `comparison()`: resample the C
curve cells with replacement, then the fixtures inside each chosen cell. The
first line of output is the reproduction check against the interval the round
actually recorded; it must match to the printed precision or the rest is void.

The point of the first table is that only the inner level responds to buying
more fixtures, so the upper limit converges to a value the panel fixes. The
point of the second is that subtracting a constant measurement floor moves a
small effect's upper limit the WRONG way: shrinking both denominators inflates
the per-case log-ratio spread faster than it separates the arms.
"""
import argparse
import collections
import json
import math
import os
import random
import statistics
from pathlib import Path

DRAWS = 2000
SEED = 751203  # the evaluator's own bootstrap seed


def walls(root):
    """case -> arm -> (median native wall over repetitions, cell)."""
    out = collections.defaultdict(dict)
    base = Path(root) / 'runs/confirmation'
    for case in sorted(os.listdir(base)):
        for arm in sorted(os.listdir(base / case)):
            seconds, cell = [], None
            for rep in sorted(os.listdir(base / case / arm)):
                r = json.loads((base / case / arm / rep / 'receipt.json').read_text())
                assert r['status'] == 'VERIFIED', (case, arm, rep, r['status'])
                seconds.append(r['native_process']['process_wall_seconds'])
                cell = r['cell']
            out[case][arm] = (statistics.median(seconds), cell)
    return out


def logs(data, candidate, baseline, floor_us=0.0):
    per_cell = collections.defaultdict(list)
    for case in sorted(data):
        b, cell = data[case][candidate]
        a, _ = data[case][baseline]
        a, b = a - floor_us * 1e-6, b - floor_us * 1e-6
        if a <= 0 or b <= 0:
            return None
        per_cell[cell].append(math.log(b / a))
    return per_cell


def interval(per_cell, fixtures=None, inner=True):
    rng = random.Random(SEED)
    cells = sorted(per_cell)
    draws = []
    for _ in range(DRAWS):
        means = []
        for _ in cells:
            xs = per_cell[rng.choice(cells)]
            if inner:
                k = len(xs) if fixtures is None else fixtures
                idx = rng.choices(range(len(xs)), k=k)
                means.append(statistics.mean(xs[i] for i in idx))
            else:
                means.append(statistics.mean(xs))
        draws.append(math.exp(statistics.mean(means)))
    draws.sort()
    point = math.exp(statistics.mean(statistics.mean(v) for v in per_cell.values()))
    return draws[int(.025 * DRAWS)], draws[min(DRAWS - 1, int(.975 * DRAWS))], point


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--round', default='runs/round-0015')
    ap.add_argument('--candidate', default='scan_io')
    ap.add_argument('--baseline', default='incumbent')
    args = ap.parse_args()

    data = walls(args.round)
    recorded = json.loads((Path(args.round) / 'decision.json').read_text())
    native = recorded['confirmation']['native_wall_ci95']

    per_cell = logs(data, args.candidate, args.baseline)
    lo, hi, _ = interval(per_cell)
    print(f'reproduction: [{lo:.8f}, {hi:.8f}] against recorded '
          f'[{native[0]:.8f}, {native[1]:.8f}] -- '
          f'{"MATCH" if abs(lo - native[0]) < 5e-8 and abs(hi - native[1]) < 5e-8 else "MISMATCH"}')

    print(f'\nresolution budget, {args.candidate}/{args.baseline} native, '
          f'{len(per_cell)} cells')
    print(f'{"fixtures per cell":>18} | {"half-width":>10} | {"upper limit":>11}')
    for k in (3, 6, 12, 30, 100, 400):
        lo, hi, _ = interval(per_cell, fixtures=k)
        print(f'{k:>18} | {(hi - lo) / 2:10.5f} | {hi:11.5f}')
    lo, hi, _ = interval(per_cell, inner=False)
    print(f'{"infinite":>18} | {(hi - lo) / 2:10.5f} | {hi:11.5f}')

    print('\nfloor removal, at the round\'s own replication')
    print(f'{"floor removed":>13} | {"cand point":>10} {"upper":>8} | '
          f'{"winner/rho point":>16} {"upper":>8}')
    for floor in (0, 108, 250, 404):
        row = f'{floor:>10} us |'
        for cand, base in ((args.candidate, args.baseline), (args.baseline, 'rho')):
            pc = logs(data, cand, base, floor_us=floor)
            if pc is None:
                row += f' {"n/a":>10} {"":>8} |'
                continue
            lo, hi, point = interval(pc)
            row += f' {point:10.4f} {hi:8.4f} |'
        print(row)


if __name__ == '__main__':
    main()
