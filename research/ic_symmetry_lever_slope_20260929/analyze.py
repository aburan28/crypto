#!/usr/bin/env python3
"""Registered readout for the symmetry-lever slope test (PREREGISTRATION.md section 5).

    python3 analyze.py RUNS_DIR
"""

import json
import statistics as st
import sys
from pathlib import Path

LOWER = 8  # d_max + 1


def fit(points):
    xs, ys = zip(*points)
    mx, my = st.mean(xs), st.mean(ys)
    return sum((x - mx) * (y - my) for x, y in points) / sum((x - mx) ** 2 for x in xs)


def main(runs):
    cells = {}
    for path in sorted(Path(runs).glob('*.jsonl')):
        rows = [json.loads(line) for line in path.read_text().splitlines() if line.strip()]
        if not rows:
            cells[path.stem] = {'rows': [], 'n': None, 'ell': None}
            continue
        cells[path.stem] = {'rows': rows, 'n': rows[0]['n'], 'ell': rows[0]['ell'], 'a': rows[0]['a']}
    print(f"{'cell':10} {'draws':>5} {'sat':>4} {'unsat':>5}  degrees                 median")
    by_curve = {}
    lower_below = {}
    for name, c in sorted(cells.items(), key=lambda kv: (kv[1]['n'] or 0, kv[1]['ell'] or 0)):
        rows = c['rows']
        unsat = [r for r in rows if r['outcome']['kind'] != 'satisfiable']
        sat = len(rows) - len(unsat)
        degs = []
        for r in unsat:
            o = r['outcome']
            if o['kind'] == 'resolved':
                degs.append(o['degree'])
            elif o['kind'] == 'at_least':
                degs.append(LOWER)
            else:
                degs.append(None)
        shown = ' '.join('caps' if d is None else (f'>={d}' if d == LOWER else str(d)) for d in degs)
        median = None
        if len(unsat) >= 3 and None not in degs:
            lowers = sum(d == LOWER for d in degs)
            median = f'>={LOWER}' if lowers * 2 > len(degs) else st.median_low(degs)
        print(f'{name:10} {len(rows):5} {sat:4} {len(unsat):5}  {shown:22}  {median}')
        if c['n'] is None:
            continue
        key = (c['a'], c['n'])
        if median == f'>={LOWER}':
            lower_below.setdefault(key, []).append(c['ell'])
        elif isinstance(median, int):
            by_curve.setdefault(key, []).append((c['ell'], median))
    print()
    slopes = {}
    for key, pts in sorted(by_curve.items()):
        if len({x for x, _ in pts}) >= 3:
            slopes[key] = fit(pts)
            print(f'K_{key[0]}/2^{key[1]}: s_n = {slopes[key]:.3f} over l = {sorted(x for x, _ in pts)}')
        else:
            print(f'K_{key[0]}/2^{key[1]}: not fitted (fewer than 3 retained l)')
    curves = {(c['a'], c['n']) for c in cells.values() if c['n'] is not None}
    print()
    if len(slopes) < 2:
        verdict = 'INCONCLUSIVE (fewer than two curves fitted)'
    else:
        mean = st.mean(slopes.values())
        early_lower_everywhere = all(any(e <= 6 for e in lower_below.get(k, [])) for k in curves)
        early_lower_any = any(e < 6 for v in lower_below.values() for e in v)
        print(f's_bar = {mean:.3f}; lower-bound cells: {dict(lower_below)}')
        if mean >= 0.35 or early_lower_everywhere:
            verdict = 'CONSTANT LEVER (symmetrisation does not flatten the degree)'
        elif mean <= 0.15 and not early_lower_any:
            verdict = 'SLOPE LEVER (degree flat in l)'
        else:
            verdict = 'INCONCLUSIVE'
    print(f'REGISTERED VERDICT: {verdict}')


if __name__ == '__main__':
    main(sys.argv[1])
