#!/usr/bin/env python3
"""Registered readout for the sized-cache timing (PREREGISTRATION.md section 4).  python3 analyze.py runs/<label>"""

import json
import math
import random
import statistics as st
import sys
from collections import defaultdict
from pathlib import Path

CELLS = ['n13a0', 'n17a1', 'n19a0', 'n19a1', 'n23a0']
SMALL = CELLS[:4]
METRICS = {'wall': 'wall_seconds', 'worker': 'worker_elapsed_seconds'}
PAIRS = [('sized', 'rho'), ('incumbent', 'rho'), ('incumbent', 'sized'), ('sized_aa', 'sized'), ('inc_new', 'incumbent')]


def value(r, metric):
    if metric == 'cpu':
        return r['user_seconds'] + r['system_seconds']
    return r[METRICS[metric]]


def ratios(rows, num, den, metric):
    by = defaultdict(lambda: defaultdict(list))
    for r in rows:
        by[r['case']][r['arm']].append(value(r, metric))
    return {case: st.median(a[num]) / st.median(a[den]) for case, a in by.items() if a.get(num) and a.get(den)}


def geo(values):
    return math.exp(st.mean(math.log(v) for v in values))


def boot(values, draws=10000, seed=20260930):
    rng = random.Random(seed)
    vals = list(values)
    samples = sorted(geo(rng.choices(vals, k=len(vals))) for _ in range(draws))
    return samples[int(.025 * draws)], samples[int(.975 * draws) - 1]


def main(directory):
    rows_all = [json.loads(line) for line in (Path(directory) / 'runs.jsonl').read_text().splitlines()]
    timed = [r for r in rows_all if r['timed']]
    contended = [r for r in timed if r['contended']]
    rows = [r for r in timed if not r['contended']]
    print(f'timed runs {len(timed)}, contended {len(contended)} (excluded), '
          f'answer mismatches {sum(not r["answer_matches"] for r in rows_all)}')
    print()
    verdict = {}
    print(f'{"cell":7} {"metric":7} ' + ' '.join(f'{a + "/" + b:>22}' for a, b in PAIRS))
    for cell in CELLS:
        cell_rows = [r for r in rows if r['cell'] == cell]
        for metric in ('wall', 'cpu', 'worker'):
            cols = []
            for num, den in PAIRS:
                per_case = ratios(cell_rows, num, den, metric)
                g, (lo, hi) = geo(per_case.values()), boot(per_case.values())
                verdict[cell, metric, num, den] = (g, lo, hi)
                cols.append(f'{g:.3f} [{lo:.3f},{hi:.3f}]')
            print(f'{cell:7} {metric:7} ' + ' '.join(f'{c:>22}' for c in cols))
    print()
    print('median minor faults per run')
    for cell in CELLS:
        parts = []
        for arm in ('incumbent', 'rho', 'sized'):
            parts.append(f'{arm} {st.median(r["minor_faults"] for r in rows if r["cell"] == cell and r["arm"] == arm):.0f}')
        print(f'  {cell}: ' + '  '.join(parts))
    print()
    w = lambda c, n, d, i=None: verdict[c, 'wall', n, d] if i is None else verdict[c, 'wall', n, d][i]
    p1 = all(verdict[c, 'worker', 'sized', 'rho'][2] < 1 for c in SMALL)
    lead_gone = [c for c in SMALL if w(c, 'incumbent', 'sized', 2) >= 1]
    p2 = len(lead_gone) >= 3
    strict = [c for c in SMALL if w(c, 'incumbent', 'sized', 1) > 1]
    p3 = all(0.92 <= w(c, 'sized_aa', 'sized', 1) and w(c, 'sized_aa', 'sized', 2) <= 1.09 for c in CELLS)
    p4 = len(contended) <= 0.05 * len(timed)
    p5 = all(0.92 <= w(c, 'inc_new', 'incumbent', 1) and w(c, 'inc_new', 'incumbent', 2) <= 1.09 for c in CELLS)
    print('REGISTERED READOUT')
    print(f'  P1 sized/rho worker-timed upper limit below 1 at all four small cells: {p1}')
    print(f'  P2 IC lead gone (incumbent/sized wall upper limit >= 1) at >= 3 of 4 small cells: {p2}  {lead_gone}')
    print(f'  P3 A/A sized_aa/sized wall inside [0.92, 1.09] at every cell: {p3}')
    print(f'  P4 at most 5% of timed runs contended: {p4}')
    print(f'  P5 layout control inc_new/incumbent wall inside [0.92, 1.09] at every cell: {p5}')
    print(f'  reported, not predicted: cells where sized is strictly faster than IC in wall time (lower limit > 1): {strict}')


if __name__ == '__main__':
    main(sys.argv[1])
