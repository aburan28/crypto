#!/usr/bin/env python3
"""Is the winner/rho ratio one number per cell, or does it move with the seed?

    python3 campaign_20260916/round19_seed_variance.py

Round 0019 pre-registered `both`/rho at n23a1 in [0.78, 1.00] with a point
prediction of 0.881, taken from rounds 0017 and 0018b pooled.  It measured
0.7498.  That is a large miss, and the interesting question is which of two
things it means:

  (a) the pooled two-seed estimate was simply noisy, and all three rounds are
      draws from one per-cell value, or
  (b) the ratio has a component that moves with the seed, over and above the
      per-case spread within a round -- in which case a per-cell number is not
      a property of the executable and panel at all, and round 0019's
      allocation bought precision about the wrong quantity.

The test is the standard one for combining estimates: if each round's per-cell
log ratio is an independent draw around one true value with its own standard
error, then the weighted scatter of the three is chi-square on two degrees of
freedom.  Excess scatter is (b).

This reads the SAME EXECUTABLE in all three rounds -- round 0017 promoted it
as `orbits`, round 0018b retained it as `incumbent`, round 0019 ran it as
`incumbent`, and its worker sha256 is identical in all three.

**The denominator needed a control.**  `rho` is not a fixed binary: the
tournament synthesises it from each round's baseline arm, and round 0017's
baseline was the PRE-orbits incumbent (`7ca9953d…`) where rounds 0018b and
0019 both use the orbits winner (`e9f263b8…`).  So round 0017's rho column is
produced by a different executable, and a cross-round ratio could be reading
that rather than the fixture draw.  Measured directly, on identical fixtures,
three cases at each of the eight cells: the two binaries return the same rho
instruction count to within 2 parts in 10,000 (geometric mean 0.99999, worst
cell 1.0002).  The orbits change never reached rho's code path, so the columns
are comparable and round 0018b's reading of its own n23a1 flip as the fixture
draw stands -- now with a control behind it instead of an assumption.

**The test is anti-conservative here, and that matters more than the p-value.**
rho's cost is a collision time, so its distribution is right-skewed with a long
tail.  A twelve-case sample standard deviation therefore underestimates the
spread, and the sampling distribution of a twelve-case mean is not normal.
Both push chi-square toward rejecting.  With three seeds, eight cells and
n=12 in two of the three rounds, treat a flagged cell as a question to put to
a fourth seed, not as an established fact.  Round 0019's n=40 column is the
only one whose mean is well resolved.
"""
import collections
import glob
import json
import math
import statistics as st
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
ROUNDS = [('round-0017', 'orbits', 2026091717),
          ('round-0018b', 'incumbent', 2026091818),
          ('round-0019', 'incumbent', 2026092119)]
CELLS = [('n13a0', 2003), ('n29a1', 42457), ('n17a1', 65587), ('n19a0', 130873),
         ('n19a1', 262543), ('n31a0', 1439393), ('n23a0', 2095853), ('n23a1', 4196903)]


def chi2_sf(x, df):
    """Upper tail of chi-square. Exact for the even df this script uses."""
    if df % 2:
        raise ValueError('even df only')
    half, term, total = df // 2, 1.0, 1.0
    for k in range(1, half):
        term *= (x / 2) / k
        total += term
    return math.exp(-x / 2) * total


def per_case_log_ratios(rnd, winner_arm, stage='confirmation'):
    """cell -> [log(winner/rho) per case], the tournament's own estimator:
    median over repetitions within a case, then one log ratio per case."""
    got = collections.defaultdict(lambda: collections.defaultdict(dict))
    for p in glob.glob(str(ROOT / f'runs/{rnd}/runs/{stage}/*/*/rep-*/receipt.json')):
        r = json.load(open(p))
        if r['status'] != 'VERIFIED' or r['total_operations'] is None:
            continue
        arm = {'rho': 'rho', winner_arm: 'winner'}.get(r['arm'])
        if arm is None:
            continue
        got[r['cell']][(arm, r['case'])].setdefault('ir', []).append(r['total_operations'])
    out = {}
    for cell, e in got.items():
        cases = sorted({c for (a, c) in e if a == 'winner'} & {c for (a, c) in e if a == 'rho'})
        out[cell] = [math.log(st.median(e[('winner', c)]['ir']) / st.median(e[('rho', c)]['ir']))
                     for c in cases]
    return out


def main():
    data = {}
    for rnd, arm, seed in ROUNDS:
        d = per_case_log_ratios(rnd, arm)
        if not d:
            raise SystemExit(f'no confirmation receipts for {rnd}; restore it from evidence/ first')
        data[rnd] = d

    print('winner/rho instruction ratio, the same executable under three seeds\n')
    head = '  '.join(f'{r.replace("round-", ""):>18}' for r, _a, _s in ROUNDS)
    print(f'{"cell":8} {"r":>9}  {head}   {"combined":>8} {"chi2(2)":>8} {"p":>7}  excess?')
    flagged = []
    for cell, r in CELLS:
        cols, weights, logs = [], [], []
        for rnd, _arm, _seed in ROUNDS:
            xs = data[rnd].get(cell) or []
            if len(xs) < 2:
                cols.append(f'{"-":>18}')
                continue
            n, mu, sd = len(xs), st.mean(xs), st.stdev(xs)
            se = sd / math.sqrt(n)
            cols.append(f'{math.exp(mu):8.4f}(n={n:2d})'.rjust(18))
            weights.append(1 / se ** 2)
            logs.append(mu)
        if len(logs) < 3:
            print(f'{cell:8} {r:9d}  ' + '  '.join(cols))
            continue
        total = sum(weights)
        combined = sum(w * x for w, x in zip(weights, logs)) / total
        chi2 = sum(w * (x - combined) ** 2 for w, x in zip(weights, logs))
        p = chi2_sf(chi2, 2)
        mark = 'YES' if p < 0.05 else 'no'
        if p < 0.05:
            flagged.append((cell, chi2, p))
        print(f'{cell:8} {r:9d}  ' + '  '.join(cols) +
              f'   {math.exp(combined):8.4f} {chi2:8.2f} {p:7.4f}  {mark}')

    print('\nEach column is one round: the geometric mean of that round\'s per-case ratios,')
    print('with its case count.  `combined` is the inverse-variance weighted value across')
    print('all three.  chi2(2) is the weighted scatter of the three about it; under the')
    print('null that a cell has ONE ratio and each round samples it, that scatter is')
    print('chi-square on two degrees of freedom.')
    print(f'\ncells with excess scatter at p < 0.05: '
          f'{[c for c, _x, _p in flagged] if flagged else "none"}')
    print(f'expected by chance at p < 0.05 over {len(CELLS)} cells: {0.05 * len(CELLS):.1f}')
    print('\nThree seeds is a small sample for this test and one round can carry the whole')
    print('statistic, so the per-round column matters as much as the p-value: read which')
    print('round is the outlier before concluding a cell is heterogeneous.  And see the')
    print('module docstring: rho\'s cost is right-skewed, so a twelve-case sample sd')
    print('underestimates the spread and this test over-rejects.  A flagged cell is a')
    print('question for a fourth seed, not a finding.')


if __name__ == '__main__':
    main()
