#!/usr/bin/env python3
"""Registered readout for the RR norm-form ladder (PREREGISTRATION.md section 5).

    python3 analyze.py RUNS_DIR
"""

import json
import statistics as st
import sys
from pathlib import Path

LOWER = {'rr': 8, 'ctrl': 8, 'x4': 10}  # d_max + 1 per arm at l <= 5; l = 6 scans rr/ctrl to 6


def full_scan_bound(arm, ell):
    return 7 if arm in ('rr', 'ctrl') and ell >= 6 else LOWER[arm]
ARMS = ('rr', 'x4', 'ctrl')


def fit(points):
    xs, ys = zip(*points)
    mx, my = st.mean(xs), st.mean(ys)
    return sum((x - mx) * (y - my) for x, y in points) / sum((x - mx) ** 2 for x in xs)


def degree(o):
    """(value, exact): a resolved degree, or a lower bound (at_least, exact=False); None for caps."""
    if o is None:
        return None
    if o['kind'] == 'resolved':
        return (o['degree'], True)
    if o['kind'] == 'at_least':
        return (o['degree'], False)
    return None  # caps or satisfiable


def show(d):
    return 'caps' if d is None else (str(d[0]) if d[1] else f'>={d[0]}')


def cell_median(degs):
    """Over >= 3 measured draws: '>=b' (b the smallest bound) if lower bounds are the
    majority, else the low median with a minority of bounds counted at their bound
    (the sym-lever convention); None if fewer than 3 draws or any caps."""
    if len(degs) < 3 or None in degs:
        return None
    lowers = [v for v, e in degs if not e]
    if len(lowers) * 2 > len(degs):
        return f'>={min(lowers)}'
    return st.median_low([v for v, _ in degs])


def main(runs):
    cells = {}
    for path in sorted(Path(runs).glob('*.jsonl')):
        rows = [json.loads(line) for line in path.read_text().splitlines() if line.strip()]
        if rows:
            cells[path.stem] = {'rows': rows, 'n': rows[0]['n'], 'ell': rows[0]['ell'], 'a': rows[0]['a']}
        else:
            cells[path.stem] = {'rows': [], 'n': None, 'ell': None, 'a': None}
    print(f"{'cell':10} {'draws':>5} {'rr-sat':>6} {'x4-sat':>6}  {'rr degrees':16} {'x4 degrees':16} {'ctrl degrees':16}  med rr  med x4  med ctrl")
    by_curve = {arm: {} for arm in ARMS}
    lower_below = {arm: {} for arm in ARMS}
    paired = []
    for name, c in sorted(cells.items(), key=lambda kv: (kv[1]['n'] or 0, kv[1]['ell'] or 0)):
        rows = c['rows']
        # one line per (draw, arm), provisional lower bounds included; the last wins
        per = {}
        for r in rows:
            per.setdefault(r['draw'], {})[r['arm']] = r
        degs = {arm: [] for arm in ARMS}
        for d, arms in sorted(per.items()):
            for arm in ARMS:
                r = arms.get(arm)
                if r is not None and r['outcome']['kind'] != 'satisfiable':
                    degs[arm].append(degree(r['outcome']))
            if 'rr' in arms and 'x4' in arms:
                d_rr, d_x4 = degree(arms['rr']['outcome']), degree(arms['x4']['outcome'])
                if arms['rr']['outcome']['kind'] != 'satisfiable' and arms['x4']['outcome']['kind'] != 'satisfiable' \
                        and d_rr is not None and d_x4 is not None and d_rr[1] and d_x4[1]:
                    paired.append((name, d_rr[0], d_x4[0]))
        sat = {arm: sum(1 for arms in per.values() if arm in arms and arms[arm]['outcome']['kind'] == 'satisfiable')
               for arm in ('rr', 'x4')}
        med = {arm: cell_median(degs[arm]) for arm in ARMS}
        print(f"{name:10} {len(per):5} {sat['rr']:6} {sat['x4']:6}  "
              + ' '.join(f"{' '.join(show(d) for d in degs[arm]):16}" for arm in ARMS)
              + f"  {str(med['rr']):6}  {str(med['x4']):6}  {med['ctrl']}")
        if c['n'] is None:
            continue
        key = (c['a'], c['n'])
        for arm in ARMS:
            if isinstance(med[arm], str) and int(med[arm][2:]) >= full_scan_bound(arm, c['ell']):
                lower_below[arm].setdefault(key, []).append(c['ell'])
            elif isinstance(med[arm], int):
                by_curve[arm].setdefault(key, []).append((c['ell'], med[arm]))
    print()
    slopes = {arm: {} for arm in ARMS}
    for arm in ARMS:
        for key, pts in sorted(by_curve[arm].items()):
            if len({x for x, _ in pts}) >= 3:
                slopes[arm][key] = fit(pts)
                print(f'{arm:4} K_{key[0]}/2^{key[1]}: s_n = {slopes[arm][key]:.3f} over l = {sorted(x for x, _ in pts)}')
            else:
                print(f'{arm:4} K_{key[0]}/2^{key[1]}: not fitted (fewer than 3 retained l)')
    print()
    if paired:
        diffs = [a - b for _, a, b in paired]
        print(f'paired rr - x4 over {len(paired)} draws resolved on both: '
              f'mean {st.mean(diffs):+.2f}, min {min(diffs):+d}, max {max(diffs):+d}; '
              f'rr below x4 on {sum(d < 0 for d in diffs)}, equal {sum(d == 0 for d in diffs)}, above {sum(d > 0 for d in diffs)}')
    curves = {(c['a'], c['n']) for c in cells.values() if c['n'] is not None}
    print()
    for arm in ('rr', 'x4'):
        s = slopes[arm]
        if len(s) < 2:
            print(f'{arm}: fewer than two curves fitted')
            continue
        mean = st.mean(s.values())
        early_lower_everywhere = all(any(e <= 6 for e in lower_below[arm].get(k, [])) for k in curves)
        early_lower_any = any(e < 6 for v in lower_below[arm].values() for e in v)
        print(f'{arm}: s_bar = {mean:.3f}; lower-bound cells: {dict(lower_below[arm])}')
        if arm == 'rr':
            if mean >= 0.35 or early_lower_everywhere:
                verdict = 'CONSTANT LEVER (the norm form does not flatten the degree in l)'
            elif mean <= 0.15 and not early_lower_any:
                verdict = 'SLOPE LEVER (rr degree flat in l)'
            else:
                verdict = 'INCONCLUSIVE'
            print(f'REGISTERED VERDICT: {verdict}')


if __name__ == '__main__':
    main(sys.argv[1])
