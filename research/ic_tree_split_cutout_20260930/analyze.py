#!/usr/bin/env python3
"""Readout and decision for the cut-out degree test (PREREGISTRATION.md §4).  python3 analyze.py RUNS_DIR"""

import json
import statistics
import sys
from pathlib import Path


def ds(arm):
    """(value, censored): D_s, or the first degree not reached when the arm stopped first."""
    if arm['D_s'] is not None:
        return arm['D_s'], False
    done = [r['D'] for r in arm['rows'] if 'rank' in r]
    return (max(done) + 1 if done else 1), True


def slope(xs, ys):
    mx, my = statistics.mean(xs), statistics.mean(ys)
    return sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sum((x - mx) ** 2 for x in xs)


def main():
    runs = Path(sys.argv[1])
    cells, failed = {}, []
    for f in sorted(runs.glob('*.json')):
        code = (runs / (f.stem + '.exit')).read_text().strip()
        if code != '0' or not f.read_text().strip():
            failed.append((f.stem, code))
            continue
        d = json.loads(f.read_text())
        cells[f.stem] = d
    print('censored cells (exit != 0):', failed or 'none')
    series = {}
    for name, d in cells.items():
        key = (d['object'], d['n'], d['a'], d['seed'])
        series.setdefault(key, []).append(d)
    verdicts = {}
    for key in sorted(series):
        obj = key[0]
        arms = ['curve', 'ZR', 'R'] if obj == 'Z' else ['curve', 'R']
        rows = sorted(series[key], key=lambda d: d['ell'])
        print(f'\n{obj} n={key[1]} a={key[2]} seed={key[3]}')
        print(f"  {'ell':>3} {'|F|':>4} {'|T2|':>6} {'size':>7}  " + '  '.join(f'{a:>9}' for a in arms)
              + '   excess false-zero fraction of curve at D_s(R)')
        pts = []
        for d in rows:
            vals = {a: ds(d[a]) for a in arms}
            dr = vals['R'][0]
            fz = next((r['false_zeros'] / r['samples'] for r in d['curve']['rows']
                       if r.get('D') == dr and 'false_zeros' in r), None)
            print(f"  {d['ell']:>3} {d['factor_base_points']:>4} {d['T2_size']:>6} {d['curve']['size']:>7}  "
                  + '  '.join(f"{('>' if c else '') + str(v - c):>9}" if c else f'{v:>9}' for v, c in vals.values())
                  + f"   {'n/a' if fz is None else round(fz, 4)}")
            if any(c for _, c in vals.values()):
                print(f"  ell={d['ell']} has a censored arm: excluded from the verdict (§4)")
            else:
                pts.append((d['ell'], vals))
        if not pts:
            continue
        ells = [e for e, _ in pts]
        gap = {e: min(v['curve'][0] - v[a][0] for a in arms if a != 'curve') for e, v in pts}
        s = slope(ells, [v['curve'][0] for _, v in pts]) if len(ells) > 1 else None
        s_r = slope(ells, [v['R'][0] for _, v in pts]) if len(ells) > 1 else None
        top = max(ells)
        verdicts[key] = {'gap_top': gap[top], 'min_gap': min(gap.values()), 'slope': s}
        print(f"  slope curve {s if s is None else round(s, 3)}  slope R {s_r if s_r is None else round(s_r, 3)}  "
              f"gap at top ell {gap[top]}  min gap {min(gap.values())}  ells {ells}")
    for obj in ('T2', 'Z'):
        vs = [v for k, v in verdicts.items() if k[0] == obj]
        if not vs:
            continue
        slopes = [v['slope'] for v in vs if v['slope'] is not None]
        med = statistics.median(slopes) if slopes else None
        alive = all(v['gap_top'] <= -2 for v in vs) and (obj == 'Z' or (med is not None and med <= 0.15))
        closed = all(v['min_gap'] >= -1 for v in vs)
        verdict = 'ALIVE' if alive else 'CLOSED' if closed else 'INCONCLUSIVE'
        print(f'\n{obj}: {len(vs)} series, median curve slope {med if med is None else round(med, 3)}  ->  {verdict}')


if __name__ == '__main__':
    main()
