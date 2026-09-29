#!/usr/bin/env python3
"""Registered readout for the isolated re-timing (PREREGISTRATION.md section 4).

    python3 analyze.py runs/<label>
"""

import json
import math
import random
import statistics as st
import sys
from collections import defaultdict
from pathlib import Path

CELLS = ['n13a0', 'n17a1', 'n19a0', 'n19a1', 'n23a0']
METRICS = {'wall': 'wall_seconds', 'cpu': None, 'worker': 'worker_elapsed_seconds'}


def value(r, metric):
    if metric == 'cpu':
        return r['user_seconds'] + r['system_seconds']
    return r[METRICS[metric]]


def ratios(rows, num, den, metric):
    """Per case: median over rounds of num divided by median over rounds of den."""
    by = defaultdict(lambda: defaultdict(list))
    for r in rows:
        by[r['case']][r['arm']].append(value(r, metric))
    out = {}
    for case, arms in by.items():
        if arms.get(num) and arms.get(den):
            out[case] = st.median(arms[num]) / st.median(arms[den])
    return out


def geo(values):
    return math.exp(st.mean(math.log(v) for v in values))


def boot(values, draws=10000, seed=20260929):
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
    print(f'median involuntary switches per run: {st.median(r["involuntary_switches"] for r in rows)}; '
          f'max other-process CPU during a kept run: {max(r["other_cpu_seconds"] for r in rows):.3f} s')
    print()
    header = f'{"cell":7} {"metric":7} {"inc/rho":>18} {"switch/rho":>18} {"rho_aa/rho (A/A)":>18}'
    print(header)
    verdict = {}
    for cell in CELLS:
        cell_rows = [r for r in rows if r['cell'] == cell]
        for metric in ('wall', 'cpu', 'worker'):
            cols = []
            for num in ('incumbent', 'switch', 'rho_aa'):
                per_case = ratios(cell_rows, num, 'rho', metric)
                g = geo(per_case.values())
                lo, hi = boot(per_case.values())
                cols.append(f'{g:.3f} [{lo:.3f},{hi:.3f}]')
                verdict[cell, metric, num] = (g, lo, hi)
            print(f'{cell:7} {metric:7} {cols[0]:>18} {cols[1]:>18} {cols[2]:>18}')
    print()
    faults = defaultdict(list)
    for r in rows:
        faults[r['cell'], r['arm']].append((r['minor_faults'], r['system_seconds']))
    print('median minor faults / system ms per run')
    for cell in CELLS:
        parts = []
        for arm in ('incumbent', 'switch', 'rho'):
            f = faults[cell, arm]
            parts.append(f'{arm} {st.median(x for x, _ in f):.0f}/{1e3 * st.median(y for _, y in f):.3f}')
        print(f'  {cell}: ' + '  '.join(parts))
    print()
    small = CELLS[:4]
    p1 = all(verdict[c, 'wall', 'incumbent'][2] < 1 for c in small)
    aa = max(max(abs(math.log(verdict[c, 'wall', 'rho_aa'][i])) for i in (1, 2)) for c in CELLS)
    p2 = all(verdict[c, 'worker', 'incumbent'][0] < verdict[c, 'wall', 'incumbent'][0] for c in small)
    p3 = verdict['n23a0', 'wall', 'incumbent'][2] >= 1
    p4 = all(0.92 <= verdict[c, 'wall', 'rho_aa'][1] and verdict[c, 'wall', 'rho_aa'][2] <= 1.09 for c in CELLS)
    p5 = len(contended) <= 0.05 * len(timed)
    print('REGISTERED READOUT')
    print(f'  P1 small-cell inc/rho wall upper limit below 1 at all four cells: {p1}')
    print(f'  P2 worker-timed inc/rho below whole-process inc/rho at all four cells: {p2}')
    print(f'  P3 n23a0 inc/rho wall upper limit at or above 1 (not below rho): {p3}')
    print(f'  P4 A/A rho_aa/rho wall interval inside [0.92, 1.09] at every cell: {p4}')
    print(f'  P5 at most 5% of timed runs contended: {p5}')
    print(f'  A/A: largest |log| bound of rho_aa/rho interval over cells: {aa:.3f} ({math.exp(aa):.3f}x)')


if __name__ == '__main__':
    main(sys.argv[1])
