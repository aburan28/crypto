#!/usr/bin/env python3
"""Readout for PROTOCOL.md.  python3 analyze.py runs/grid ../walltime_strong_rho_20260929/runs/registered/runs.jsonl"""
import json, math, statistics as st, sys
from collections import defaultdict

CELLS = ['n13a0', 'n17a1', 'n19a0', 'n19a1', 'n23a0']
PEN = {'mispredict': 15, 'd1': 12, 'll': 150}


def penalties(t, scale=1.0):
    mis = t['Bcm'] + t['Bim']
    ll = t['ILmr'] + t['DLmr'] + t['DLmw']
    d1 = t['I1mr'] + t['D1mr'] + t['D1mw'] - ll
    return scale * (PEN['mispredict'] * mis + PEN['d1'] * d1 + PEN['ll'] * ll)


def main(grid, timing):
    rows = json.load(open(f'{grid}/counts.json'))
    agg = defaultdict(lambda: defaultdict(int))
    for r in rows:
        for k, v in r['totals'].items():
            agg[r['cell'], r['arm']][k] += v
    med = defaultdict(list)
    for line in open(timing):
        r = json.loads(line)
        if r['timed'] and not r['contended'] and r['arm'] in ('incumbent', 'lean_rho'):
            med[r['cell'], r['arm'], r['case']].append(r['in_process_seconds'])
    print(f"{'cell':6} {'arm':4} {'mispred/kI':>10} {'D1miss/kI':>9} {'LLmiss/kI':>9} {'pen cyc/I':>9}")
    for c in CELLS:
        for a in ('ic', 'rho'):
            t = agg[c, a]
            k = t['Ir'] / 1000
            ll = t['ILmr'] + t['DLmr'] + t['DLmw']
            print(f"{c:6} {a:4} {(t['Bcm'] + t['Bim']) / k:10.2f} {(t['I1mr'] + t['D1mr'] + t['D1mw']) / k:9.2f} "
                  f"{ll / k:9.2f} {penalties(t) / t['Ir']:9.3f}")
    print()
    print('Measured throughput T = (Ir_ic/Ir_rho) / (time_ic/time_rho); time is in-process, lean rho.')
    print('c0* = common base cycles/instruction at which the penalty model alone reproduces T.')
    print(f"{'cell':6} {'Ir ratio':>8} {'time ratio':>10} {'T':>6} {'c0* (x0.5, x1, x2 penalties)':>32}")
    for c in CELLS:
        ic, rho = agg[c, 'ic'], agg[c, 'rho']
        ir_ratio = ic['Ir'] / rho['Ir']
        cases = {k[2] for k in med if k[0] == c}
        ratios = [st.median(med[c, 'incumbent', x]) / st.median(med[c, 'lean_rho', x]) for x in cases
                  if med[c, 'incumbent', x] and med[c, 'lean_rho', x]]
        time_ratio = math.exp(st.mean(map(math.log, ratios)))
        T = ir_ratio / time_ratio
        out = []
        for s in (0.5, 1.0, 2.0):
            p_ic, p_rho = penalties(ic, s) / ic['Ir'], penalties(rho, s) / rho['Ir']
            out.append((p_rho - T * p_ic) / (T - 1))
        print(f"{c:6} {ir_ratio:8.3f} {time_ratio:10.3f} {T:6.2f} " + '  '.join(f'{x:8.3f}' for x in out))
    print()
    print('Top functions by mispredicts (rho, n19a0-000):')
    for r in rows:
        if r['case'] == 'n19a0-000' and r['arm'] == 'rho':
            for f in r['top_mispredict_fns'][:3]:
                print(f"  {f['Bcm'] + f['Bim']:6d} mispredicts / {f['Bc'] + f['Bi']:7d} branches  {f['fn'][:110]}")
        if r['case'] == 'n19a0-000' and r['arm'] == 'ic':
            for f in r['top_mispredict_fns'][:2]:
                print(f"  IC: {f['Bcm'] + f['Bim']:6d} / {f['Bc'] + f['Bi']:7d}  {f['fn'][:110]}")


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
