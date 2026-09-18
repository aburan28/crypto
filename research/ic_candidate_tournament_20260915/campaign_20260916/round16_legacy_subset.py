#!/usr/bin/env python3
"""Recompute a round's verdict over a SUBSET of its cells.

    python3 campaign_20260916/round16_legacy_subset.py --round runs/round-0016 \
        --cells n13a0,n17a1,n19a0,n19a1,n23a0

Round 0016 widened the panel from five cells to eight. A ratio measured over
eight cells is not comparable to one measured over five, so the widened round
would otherwise break every cross-round comparison the campaign has published.
It does not, because the legacy five are a subset of the eight: this restricts
the frozen receipts to those cells and re-runs the evaluator's own nested
bootstrap over them, giving the number round 0015 would have reported.

This writes nothing. It is an analysis of committed receipts, and its output
belongs in a report next to the round's own verdict, never in place of it.
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


def comparison(root, stage, candidate, baseline, keep):
    """The evaluator's comparison(), restricted to `keep` cells."""
    base = Path(root) / 'runs' / stage
    per_case = collections.defaultdict(dict)
    for case in sorted(os.listdir(base)):
        for arm in (candidate, baseline):
            path = base / case / arm
            if not path.is_dir():
                continue
            ops, wall, cell = [], [], None
            for rep in sorted(os.listdir(path)):
                r = json.loads((path / rep / 'receipt.json').read_text())
                if r['status'] != 'VERIFIED' or r['total_operations'] is None:
                    raise SystemExit(f'unverified or unpriced: {case}/{arm}/{rep}')
                ops.append(r['total_operations'])
                wall.append(r['native_process']['process_wall_seconds'])
                cell = r['cell']
            per_case[case][arm] = (statistics.median(ops), statistics.median(wall), cell)

    logs, native = collections.defaultdict(list), collections.defaultdict(list)
    cases = 0
    for case in sorted(per_case):
        if candidate not in per_case[case] or baseline not in per_case[case]:
            continue
        ob, wb, cell = per_case[case][candidate]
        oa, wa, _ = per_case[case][baseline]
        if cell not in keep:
            continue
        logs[cell].append(math.log(ob / oa))
        native[cell].append(math.log(wb / wa))
        cases += 1
    if not logs:
        raise SystemExit('no cases survived the cell filter')

    rng = random.Random(SEED)
    cells = sorted(logs)
    boot, native_boot = [], []
    for _ in range(DRAWS):
        means, wall_means = [], []
        for _ in cells:
            cell = rng.choice(cells)
            idx = rng.choices(range(len(logs[cell])), k=len(logs[cell]))
            means.append(statistics.mean(logs[cell][i] for i in idx))
            wall_means.append(statistics.mean(native[cell][i] for i in idx))
        boot.append(math.exp(statistics.mean(means)))
        native_boot.append(math.exp(statistics.mean(wall_means)))
    boot.sort()
    native_boot.sort()
    lo, hi = int(.025 * DRAWS), min(DRAWS - 1, int(.975 * DRAWS))
    return {
        'candidate': candidate, 'baseline': baseline, 'stage': stage,
        'cells': cells, 'paired_cases': cases,
        'candidate_over_baseline': math.exp(statistics.mean(
            statistics.mean(v) for v in logs.values())),
        'ci95': [boot[lo], boot[hi]],
        'per_cell': {c: math.exp(statistics.mean(v)) for c, v in logs.items()},
        'native_wall_candidate_over_baseline': math.exp(statistics.mean(
            statistics.mean(v) for v in native.values())),
        'native_wall_ci95': [native_boot[lo], native_boot[hi]],
        'native_wall_per_cell': {c: math.exp(statistics.mean(v)) for c, v in native.items()},
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--round', required=True)
    ap.add_argument('--cells', default='n13a0,n17a1,n19a0,n19a1,n23a0')
    ap.add_argument('--stages', default='confirmation,replay')
    ap.add_argument('--pairs', default='incumbent:rho,scan_io:incumbent,scan_io:rho')
    args = ap.parse_args()
    keep = {c.strip() for c in args.cells.split(',') if c.strip()}
    for stage in [s.strip() for s in args.stages.split(',')]:
        for pair in [p.strip() for p in args.pairs.split(',')]:
            candidate, baseline = pair.split(':')
            try:
                result = comparison(args.round, stage, candidate, baseline, keep)
            except SystemExit as reason:
                print(f'{stage} {pair}: {reason}')
                continue
            print(json.dumps(result, indent=1, sort_keys=True))


if __name__ == '__main__':
    main()
