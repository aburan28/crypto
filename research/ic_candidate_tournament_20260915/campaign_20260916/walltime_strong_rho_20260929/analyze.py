#!/usr/bin/env python3
"""Registered readout (PREREGISTRATION.md section 4).  python3 analyze.py runs/<label>"""

import json
import math
import random
import statistics as st
import sys
from collections import defaultdict
from pathlib import Path

CELLS = ['n13a0', 'n17a1', 'n19a0', 'n19a1', 'n23a0']
RHOS = ['lean_rho', 'strong32', 'strong1']
KEY = {'inproc': 'in_process_seconds', 'wall': 'wall_seconds'}


def geo(v):
    return math.exp(st.mean(math.log(x) for x in v))


def boot(v, draws=10000, seed=20260929):
    rng = random.Random(seed)
    v = list(v)
    s = sorted(geo(rng.choices(v, k=len(v))) for _ in range(draws))
    return s[int(.025 * draws)], s[int(.975 * draws) - 1]


def main(directory):
    rows = [json.loads(line) for line in (Path(directory) / 'runs.jsonl').read_text().splitlines()]
    timed = [r for r in rows if r['timed']]
    kept = [r for r in timed if not r['contended']]
    print(f'timed {len(timed)}, contended {len(timed) - len(kept)} (excluded), '
          f'answer failures {sum(not r["answer_ok"] for r in rows)}')
    med = defaultdict(dict)  # (cell, metric) -> {(case index, arm): median}
    by = defaultdict(list)
    for r in kept:
        by[r['cell'], r['case'].split('-')[1], r['arm']].append(r)
    for (cell, idx, arm), rs in by.items():
        for m, k in KEY.items():
            med[cell, m][idx, arm] = st.median(x[k] for x in rs)
    verdict = {}
    for m in ('inproc', 'wall'):
        print(f'\n## {m}: absolute geometric means (ms) and ratios, 95% bootstrap over cases')
        for cell in CELLS:
            d = med[cell, m]
            arms = RHOS + ['incumbent', 'strong32_aa']
            every = sorted({i for i, _ in d})
            idxs = [i for i in every if all((i, a) in d for a in arms)]
            if len(idxs) < len(every):
                print(f'  {cell}: {len(every) - len(idxs)} case(s) dropped, an arm had only contended runs')
            absolute = {a: geo(d[i, a] for i in idxs) * 1e3 for a in RHOS + ['incumbent', 'strong32_aa']}
            strongest = min(RHOS, key=lambda a: absolute[a])
            ratio = [d[i, 'incumbent'] / d[i, strongest] for i in idxs]
            aa = [d[i, 'strong32_aa'] / d[i, 'strong32'] for i in idxs]
            g, (lo, hi) = geo(ratio), boot(ratio)
            ga, (alo, ahi) = geo(aa), boot(aa)
            verdict[cell, m] = (strongest, g, lo, hi, alo, ahi)
            parts = ' '.join(f'{a} {absolute[a]:.3f}' for a in ['incumbent'] + RHOS)
            print(f'  {cell}: {parts} | strongest rho: {strongest} | IC/strongest {g:.3f} [{lo:.3f},{hi:.3f}]'
                  f' | A/A {ga:.3f} [{alo:.3f},{ahi:.3f}]')
    small = CELLS[:4]
    p1 = {c: verdict[c, 'inproc'][0] for c in CELLS}
    p2 = all(verdict[c, 'inproc'][3] < 1 for c in small)
    p3 = all(0.9 <= verdict[c, 'inproc'][4] and verdict[c, 'inproc'][5] <= 1.1 for c in CELLS)
    p4 = len(timed) - len(kept) <= 0.05 * len(timed)
    print('\nREGISTERED READOUT')
    print(f'  P1 (reported) strongest rho in-process per cell: {p1}')
    print(f'  P2 IC / strongest rho in-process upper limit below 1 at all four small cells: {p2}')
    print(f'  P3 strong32 A/A in-process interval inside [0.9, 1.1] at every cell: {p3}')
    print(f'  P4 at most 5% of timed runs contended: {p4}')


if __name__ == '__main__':
    main(sys.argv[1])
