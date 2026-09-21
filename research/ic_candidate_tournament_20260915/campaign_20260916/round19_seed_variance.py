#!/usr/bin/env python3
"""Is the winner/rho ratio one number per cell, or does it move with the seed?

    python3 campaign_20260916/round19_seed_variance.py

Round 0019 pre-registered `both`/rho at n23a1 in [0.78, 1.00] with a point
prediction of 0.881, taken from rounds 0017 and 0018b pooled.  It measured
0.7498.  That is a large miss, and the interesting question is which of two
things it means:

  (a) the pooled two-seed estimate was simply noisy, and every round is a
      draw from one per-cell value, or
  (b) the ratio has a component that moves with the seed, over and above the
      per-case spread within a round -- in which case a per-cell number is not
      a property of the executable and panel at all, and round 0019's
      allocation bought precision about the wrong quantity.

The test is the standard one for combining estimates: if each round's per-cell
log ratio is an independent draw around one true value with its own standard
error, then their weighted scatter is chi-square on one fewer degree of freedom
than there are rounds.  Excess scatter is (b).

This reads the SAME EXECUTABLE in every round -- round 0017 promoted it as
`orbits` and later rounds retained it as `incumbent` -- and does not take the
name map on trust: `check_same_executable` hashes each round's named arm and
refuses to run if any of them is not the reference binary, because the
`arm_sha256` recorded in a receipt is round-scoped and cannot recognise it.
A round whose binaries are not restored, or whose confirmation stage has not
finished, is skipped with its reason printed rather than contributing a
partial column.

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
Both push chi-square toward rejecting.  With few seeds, eight cells and n=12 in
most rounds, treat a flagged cell as a question to put to another seed, not as
an established fact.  Only the n=40 columns -- round 0019 onward, at n23a1 --
carry a well-resolved mean.
"""
import collections
import glob
import hashlib
import json
import math
import statistics as st
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
# (round, the arm running the reference executable, seed).  Adding a round is
# one line here -- and `check_same_executable` below refuses to proceed if the
# named arm is not in fact the reference, because `arm_sha256` in a receipt is
# round-scoped and cannot be used to recognise it.
ROUNDS = [('round-0017', 'orbits', 2026091717),
          ('round-0018b', 'incumbent', 2026091818),
          ('round-0019', 'incumbent', 2026092119),
          ('round-0020', 'incumbent', 2026092120)]
# The worker every round above must be running for its columns to be comparable:
# round 0017 promoted it as `orbits`, rounds 0018b onward retained it.
REFERENCE_WORKER = 'e9f263b839d439424f0d96fb16092cfd9bd21bbf6772fd94c1cb706156c5539c'


def worker_path(rnd, arm):
    """Where a round keeps the binary for one of its arms."""
    baseline = ROOT / f'runs/{rnd}/worker'
    candidate = ROOT / f'runs/{rnd}/source_candidates/{arm}/worker'
    return candidate if candidate.exists() else baseline


def check_same_executable():
    """Every round's named arm must be the SAME BINARY, or the comparison is
    between executables rather than between fixture draws.

    This is not hypothetical.  `rho` is synthesised from each round's baseline
    arm, and round 0017's baseline was the pre-orbits incumbent where later
    rounds use the orbits winner -- measured equal to within 2 parts in 10,000
    on rho's code path, but only because that change never reached it.  A silent
    name-map error on the numerator would not be so kind.

    Binaries are gitignored, so a round whose evidence has not been restored is
    reported and skipped rather than assumed.
    """
    usable, skipped = [], []
    for rnd, arm, seed in ROUNDS:
        p = worker_path(rnd, arm)
        if not p.exists():
            skipped.append((rnd, 'binary not restored'))
            continue
        # A confirmation stage still being written would contribute a column
        # built from however many cases happened to exist when this ran, which
        # is not a seed's value. The stage summary appears only once the stage
        # is complete, so it is the gate.
        if not (ROOT / f'runs/{rnd}/summaries/confirmation.json').exists():
            skipped.append((rnd, 'confirmation stage not finished'))
            continue
        got = hashlib.sha256(p.read_bytes()).hexdigest()
        if got != REFERENCE_WORKER:
            raise SystemExit(
                f'{rnd}: arm {arm!r} hashes to {got[:16]}…, not the reference '
                f'{REFERENCE_WORKER[:16]}…. Either the arm name is wrong for this round or '
                f'the round did not run the reference executable; in both cases its column '
                f'would not be comparable.')
        usable.append((rnd, arm, seed))
    for rnd, why in skipped:
        print(f'  note: {rnd} skipped -- {why}')
    return usable
CELLS = [('n13a0', 2003), ('n29a1', 42457), ('n17a1', 65587), ('n19a0', 130873),
         ('n19a1', 262543), ('n31a0', 1439393), ('n23a0', 2095853), ('n23a1', 4196903)]


def chi2_sf(x, df):
    """Upper tail of chi-square, for any positive integer df.

    The first version of this took a shortcut that was exact for even df and
    raised on odd, which was fine for three rounds and wrong for four.  This is
    the regularised upper incomplete gamma Q(df/2, x/2): series below the
    transition, Lentz's continued fraction above, which is the standard split.
    """
    if df < 1:
        raise ValueError('df must be at least 1')
    a, z = df / 2.0, x / 2.0
    if z <= 0:
        return 1.0
    log_pref = a * math.log(z) - z - math.lgamma(a)
    if z < a + 1.0:                       # series for the LOWER tail P(a, z)
        term = total = 1.0 / a
        for n in range(1, 1000):
            term *= z / (a + n)
            total += term
            if abs(term) < abs(total) * 1e-15:
                break
        return 1.0 - math.exp(log_pref + math.log(total))
    tiny = 1e-300                          # continued fraction for Q(a, z)
    b, c, d = z + 1.0 - a, 1.0 / tiny, 1.0 / (z + 1.0 - a)
    h = d
    for i in range(1, 1000):
        an = -i * (i - a)
        b += 2.0
        d = an * d + b
        if abs(d) < tiny:
            d = tiny
        c = b + an / c
        if abs(c) < tiny:
            c = tiny
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < 1e-15:
            break
    return math.exp(log_pref + math.log(h))


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
    rounds = check_same_executable()
    if len(rounds) < 3:
        raise SystemExit(f'need at least three rounds with restored binaries, have {len(rounds)}')
    data = {}
    for rnd, arm, _seed in rounds:
        d = per_case_log_ratios(rnd, arm)
        if not d:
            raise SystemExit(f'no confirmation receipts for {rnd}; restore it from evidence/ first')
        data[rnd] = d
    print(f'  verified the same executable ({REFERENCE_WORKER[:16]}...) in '
          f'{len(rounds)} rounds: {", ".join(r for r, _a, _s in rounds)}\n')

    print(f'winner/rho instruction ratio, the same executable under {len(rounds)} seeds\n')
    head = '  '.join(f'{r.replace("round-", ""):>18}' for r, _a, _s in rounds)
    print(f'{"cell":8} {"r":>9}  {head}   {"combined":>8} {"chi2":>8} {"p":>7}  excess?')
    flagged = []
    for cell, r in CELLS:
        cols, weights, logs = [], [], []
        for rnd, _arm, _seed in rounds:
            xs = data[rnd].get(cell) or []
            if len(xs) < 2:
                cols.append(f'{"-":>18}')
                continue
            n, mu, sd = len(xs), st.mean(xs), st.stdev(xs)
            se = sd / math.sqrt(n)
            cols.append(f'{math.exp(mu):8.4f}(n={n:2d})'.rjust(18))
            weights.append(1 / se ** 2)
            logs.append(mu)
        if len(logs) < len(rounds):
            print(f'{cell:8} {r:9d}  ' + '  '.join(cols))
            continue
        total = sum(weights)
        combined = sum(w * x for w, x in zip(weights, logs)) / total
        chi2 = sum(w * (x - combined) ** 2 for w, x in zip(weights, logs))
        p = chi2_sf(chi2, len(rounds) - 1)
        mark = 'YES' if p < 0.05 else 'no'
        if p < 0.05:
            flagged.append((cell, chi2, p))
        print(f'{cell:8} {r:9d}  ' + '  '.join(cols) +
              f'   {math.exp(combined):8.4f} {chi2:8.2f} {p:7.4f}  {mark}')

    print('\nEach column is one round: the geometric mean of that round\'s per-case ratios,')
    print('with its case count.  `combined` is the inverse-variance weighted value across')
    print(f'all {len(rounds)}.  The chi2 column is the weighted scatter of the {len(rounds)} about it; under')
    print('the null that a cell has ONE ratio and each round samples it, that scatter is')
    print(f'chi-square on {len(rounds) - 1} degrees of freedom.')
    print(f'\ncells with excess scatter at p < 0.05: '
          f'{[c for c, _x, _p in flagged] if flagged else "none"}')
    print(f'expected by chance at p < 0.05 over {len(CELLS)} cells: {0.05 * len(CELLS):.1f}')
    print(f'\n{len(rounds)} seeds is a small sample for this test and one round can carry the whole')
    print('statistic, so the per-round column matters as much as the p-value: read which')
    print('round is the outlier before concluding a cell is heterogeneous.  And see the')
    print('module docstring: rho\'s cost is right-skewed, so a twelve-case sample sd')
    print('underestimates the spread and this test over-rejects.  A flagged cell is a')
    print('question for a fourth seed, not a finding.')


if __name__ == '__main__':
    main()
