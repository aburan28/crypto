#!/usr/bin/env python3
"""Round 0019's confirmation allocation and the tables it rests on.

    python3 campaign_20260916/round19_allocate.py [--target 0.01] [--floor 12]

Reads FROZEN ROUNDS 0017 AND 0018b ONLY -- 24 independent fixtures a cell under
two independent seeds, for the two executables round 0019 runs -- and emits

  A. the per-cell allocation, as the smallest multiple of four whose one-sided
     failure probability against the strict rho gate is at most `--target` in
     BOTH metrics, floored at the profile's flat count so no cell is ever
     measured less than before;
  B. where rho's variance lives, per cell;
  C. the phase ledger, and the same-degree contrast that is the only clean
     subgroup-order-only comparison the panel offers.

The strict rho gate asks every cell's winner/rho ratio to be below one in both
metrics.  That per-cell test is a point estimate of a geometric mean over the
cell's fixtures, so its resolution is the cell's own spread over its own sample
count -- and the spread differs ninefold across this panel, because rho's
collision search is where the variance is and it dominates rho's cost only at
the large-subgroup cells.  A flat count therefore buys 1e-8 at one end of the
panel and 9% at the other.  Round 0018b spent that 9%.

Nothing here reads round 0019.  The allocation raises only: more fixtures move
a cell's estimate toward its true value in whichever direction that lies, so
this removes a coin flip without choosing how the coin was weighted.
"""
import argparse
import collections
import glob
import json
import math
import statistics as st
from math import erf, sqrt
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
# The same executable under two names: round 0017 called it `orbits`, the
# challenger it promoted; round 0018b called it `incumbent`, the winner it
# retained.  Identical worker sha256 in both rounds.
WINNER = {'round-0017': 'orbits', 'round-0018b': 'incumbent'}
CELLS = [('n13a0', 13, 2003), ('n17a1', 17, 65587), ('n19a0', 19, 130873),
         ('n19a1', 19, 262543), ('n23a0', 23, 2095853), ('n23a1', 23, 4196903),
         ('n29a1', 29, 42457), ('n31a0', 31, 1439393)]
BY_R = sorted(CELLS, key=lambda c: c[2])
SHARED = {'curve_and_targets', 'final_verification', 'startup_and_input', 'reporting_and_cleanup'}

P = lambda z: 0.5 * (1 - erf(z / sqrt(2)))
gm = lambda v: math.exp(st.mean([math.log(x) for x in v]))


def frozen():
    """(round, stage, cell) -> (arm, case) -> {'ir': [...], 'w': [...], 'phase': [...]}."""
    d = collections.defaultdict(lambda: collections.defaultdict(dict))
    for rnd, winner_arm in WINNER.items():
        for stage in ('confirmation', 'replay'):
            for p in glob.glob(str(ROOT / f'runs/{rnd}/runs/{stage}/*/*/rep-*/receipt.json')):
                r = json.load(open(p))
                if r['status'] != 'VERIFIED':
                    continue
                arm = {'rho': 'rho', winner_arm: 'winner', 'both': 'both'}.get(r['arm'])
                if arm is None:
                    continue
                e = d[(rnd, stage, r['cell'])][(arm, r['case'])]
                e.setdefault('ir', []).append(r['total_operations'])
                e.setdefault('w', []).append(r['native_process']['process_wall_seconds'])
                e.setdefault('phase', []).append(r['phase_costs'])
    if not d:
        raise SystemExit('no frozen receipts; restore rounds 0017 and 0018b from evidence/ first')
    return d


def log_ratios(D, num, den, metric):
    """Per-case log ratios, exactly as the tournament forms them: median over
    repetitions within a case, then one log ratio per case."""
    out = collections.defaultdict(list)
    for (_rnd, _stage, cell), e in D.items():
        for c in sorted({c for (a, c) in e if a == num} & {c for (a, c) in e if a == den}):
            out[cell].append(math.log(st.median(e[(num, c)][metric]) / st.median(e[(den, c)][metric])))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--target', type=float, default=0.01,
                    help='per-cell, per-metric one-sided failure probability to design for')
    ap.add_argument('--floor', type=int, default=12, help="the pilot profile's flat count")
    ap.add_argument('--json', type=Path, help='also write the allocation here')
    args = ap.parse_args()

    D = frozen()
    WR = {m: log_ratios(D, 'winner', 'rho', m) for m in ('ir', 'w')}
    BW = {m: log_ratios(D, 'both', 'winner', m) for m in ('ir', 'w')}

    alloc, rows = {}, {}
    print('## A. Per-cell resolution and the allocation\n')
    print('| cell | r | Ir `both`/rho | per-case sd | native `both`/rho | per-case sd '
          '| P(fail) at 12 | cases | P(fail) at that count |')
    print('|:--|--:|--:|--:|--:|--:|--:|--:|--:|')
    for cell, _deg, r in CELLS:
        est = {}
        for metric, key in (('ir', 'ir'), ('nat', 'w')):
            # `both` is the incumbent plus two patches, so its expected ratio to
            # rho is the pooled winner/rho ratio times the measured both/winner
            # factor.  The spread is the winner/rho spread: it dominates, since
            # the two IC arms differ by a near-constant factor within a cell.
            est[metric] = (st.mean(WR[key][cell]) + st.mean(BW[key][cell]), st.stdev(WR[key][cell]))
        n = args.floor
        while n < 256 and not all(P(-mu / (sd / math.sqrt(n))) <= args.target for mu, sd in est.values()):
            n += 4
        alloc[cell] = n
        rows[cell] = est
        p12 = max(P(-mu / (sd / math.sqrt(args.floor))) for mu, sd in est.values())
        pn = max(P(-mu / (sd / math.sqrt(n))) for mu, sd in est.values())
        (mi, si), (mn, sn) = est['ir'], est['nat']
        print(f'| `{cell}` | {r:,} | {math.exp(mi):.4f} | {si:.3f} | {math.exp(mn):.4f} | {sn:.3f} '
              f'| {p12:.4f} | **{n}** | {pn:.4f} |')
    total = sum(alloc.values())
    print(f'\nConfirmation cases: **{total}** against the flat panel\'s {args.floor * len(CELLS)}. '
          f'No cell is measured less than before.\n')
    flag = ','.join(f'{k}={v}' for k, v in alloc.items() if v != args.floor)
    print(f'`--confirmation-cases {flag}`\n')

    print('## B. Where rho\'s variance lives\n')
    print('| cell | r | rho solve gm (Ir) | sd(log solve) | solve share of rho | sd(log winner/rho) |')
    print('|:--|--:|--:|--:|--:|--:|')
    for cell, _deg, r in BY_R:
        solve, tot = [], []
        for (_rnd, _stage, c2), e in D.items():
            if c2 != cell:
                continue
            for (arm, _case), v in e.items():
                if arm == 'rho':
                    solve.append(v['phase'][0]['rho_solve'])
                    tot.append(sum(v['phase'][0].values()))
        print(f'| `{cell}` | {r:,} | {gm(solve):,.0f} | {st.stdev([math.log(x) for x in solve]):.3f} '
              f'| {st.mean(solve) / st.mean(tot):.3f} | {st.stdev(WR["ir"][cell]):.3f} |')

    print('\n## C. Phase ledger, and the same-degree contrast\n')
    print('| cell | r | fb+tables | collection | individual log | shared fixed | rho solve | winner/rho |')
    print('|:--|--:|--:|--:|--:|--:|--:|--:|')
    phase = {}
    for cell, _deg, r in BY_R:
        acc = collections.defaultdict(list)
        for (_rnd, _stage, c2), e in D.items():
            if c2 != cell:
                continue
            for (arm, _case), v in e.items():
                pc = v['phase'][0]
                if arm == 'winner':
                    acc['fb'].append(pc['factor_base_and_tables'])
                    acc['col'].append(pc['collection_and_decomposition'])
                    acc['il'].append(pc['individual_log'])
                    acc['fx'].append(sum(x for k, x in pc.items() if k in SHARED))
                elif arm == 'rho':
                    acc['rho'].append(pc['rho_solve'])
        ratio = math.exp(st.mean(WR['ir'][cell]))
        phase[cell] = {k: gm(v) for k, v in acc.items()} | {'ratio': ratio}
        print(f'| `{cell}` | {r:,} | {gm(acc["fb"]):,.0f} | {gm(acc["col"]):,.0f} | {gm(acc["il"]):,.0f} '
              f'| {gm(acc["fx"]):,.0f} | {gm(acc["rho"]):,.0f} | {ratio:.4f} |')

    print('\nThe two same-degree pairs are the only clean subgroup-order-only contrast the')
    print('panel offers: `size` and the field arithmetic are identical within a pair.\n')
    print('| pair | r ratio | collection ratio | rho solve ratio | winner/rho ratio |')
    print('|:--|--:|--:|--:|--:|')
    for lo, hi in (('n19a0', 'n19a1'), ('n23a0', 'n23a1')):
        rl = dict((c, r) for c, _d, r in CELLS)
        a, b = phase[lo], phase[hi]
        print(f'| `{lo}` → `{hi}` | {rl[hi] / rl[lo]:.3f} | {b["col"] / a["col"]:.3f} '
              f'| {b["rho"] / a["rho"]:.3f} | {b["ratio"] / a["ratio"]:.3f} |')
    b23 = phase['n23a1']['ratio'] * (phase['n23a1']['ratio'] / phase['n23a0']['ratio'])
    print(f'\nOne further doubling of `r` at degree 23, at the measured rate, reads '
          f'{b23:.3f}: above rho.')

    if args.json:
        args.json.write_text(json.dumps({'allocation': alloc, 'target': args.target,
                                         'floor': args.floor, 'confirmation_cases': total},
                                        indent=1) + '\n')
        print('\nwrote', args.json)


if __name__ == '__main__':
    main()
