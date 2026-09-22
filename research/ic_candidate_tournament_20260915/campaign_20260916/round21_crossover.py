#!/usr/bin/env python3
"""Measure where this collector crosses rho, on the cell that was out of reach.

    *** CORRECTED BY round22_crossover.py -- see ROUND22-budget-and-curvature.md ***

    This script checks that the IC arm completed and then measures BOTH arms.
    `instructions()` reads callgrind's `Collected:` line and never sees the
    worker's JSON status, so a rho run cut off at `max_trials` contributes a
    TRUNCATED instruction count as the denominator.  On this script's own
    draw, 4 of the 64 rho runs at `n37a0` did not complete (64 of 64 did at
    `n23a1`).  An under-charged denominator inflates IC/rho, so the 1.533
    [1.238, 1.899] below is an upper-biased reading and the correction runs
    downward.  `n23a1`'s 0.831 stands.

    The code below is LEFT EXACTLY AS IT RAN.  Fixing it here would make the
    published numbers unreproducible, which is the one thing a correction
    must not do; the fix, and the measurement that supersedes these numbers,
    are in round22_crossover.py.

    python3 campaign_20260916/round21_crossover.py WORKER [fixtures]

WORKER is a worker built from the round-0020 winner plus
`round21-wide-pair-table.patch`; `round21_build.sh` produces one. The patched
worker is required: the promoted worker refuses degree 37 outright.

Rounds 0019 and 0020 put this collector at `Theta(r^(2/3))` against rho's
`Theta(r^(1/2))` from the solver's own cost model, and measured the ratio
rising 16.6% per doubling of `r` at degree 23 -- an extrapolation from two
points, because `n23a1` at r = 4,196,903 was the largest cell the panel could
reach. PROBE-degree-ceiling.md then found why it was the largest: the pair
table packed its coordinates into `u32`, so `koblitz_tiny_ic::MAX_DEGREE` was
31 and everything above it fell back to a general solver that emits no descent
relation. The field has always been one `u64`; the packing was the whole
ceiling.

With the packing widened, `n37a0` (r = 230,603,167, 55x the panel's largest)
runs the real collector and the crossover stops being an extrapolation.

Every IC report here is verified by `oracle.py`, which is what makes these
numbers admissible where PROBE-degree-ceiling.md's were not: that probe
measured the general fallback path, whose output the checker rejects as "not
certified as index calculus".

Two things this script does NOT do, and both matter:

  * It does not compare at `n41a0` (r = 5.5e11). IC completes there with a
    large enough base, and rho does not complete at all -- returning in 0.06s,
    far too fast for the ~72,600 steps it would need, so it is the frozen
    trial budget cutting rho off rather than rho being out-run. An arm that
    does not complete has no cost, and by the campaign's own rule a budget
    exhaustion is never evidence about an algorithm. There is no comparison to
    report at that cell.

  * It does not claim the factor base was tuned for `n37a0`. It sweeps the
    base and reports every size, because the interesting part is that bigger
    is monotonically WORSE -- the same thing round 0019's sweep found at all
    eight panel cells, now reproduced two decades of `r` higher.
"""
import json
import math
import random
import re
import statistics as st
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
import tournament as T  # noqa: E402
from oracle import verify  # noqa: E402

CONFIG = json.loads((ROOT / 'runs/round-0020/candidates.json').read_text())[0]['config']
# (label, degree, curve_a, r, requested point counts). The point count only
# selects how many orbit batches the sampler builds; the achieved orbit count
# is read back from the report and is what the table reports.
CELLS = [('n23a1', 23, 1, 4_196_903, (138,)),
         ('n37a0', 37, 0, 230_603_167, (222, 1777, 2400, 3600))]
Z = 1.959963985  # two-sided 95%
SEEDS = (20260922, 4242)
gm = lambda v: math.exp(st.mean([math.log(x) for x in v]))


def run(worker, job, env, timeout=1200):
    done = subprocess.run([worker], input=json.dumps(job), text=True,
                          capture_output=True, env=env, timeout=timeout)
    return json.loads(done.stdout) if done.stdout else {}


def instructions(worker, job, env, timeout=2400):
    done = subprocess.run(
        ['valgrind', '--tool=callgrind', '--callgrind-out-file=/dev/null', worker],
        input=json.dumps(job), text=True, capture_output=True, env=env, timeout=timeout)
    found = re.search(r'Collected\s*:\s*(\d+)', done.stderr)
    if not found:
        raise RuntimeError(f'no Ir collected: {done.stderr[-300:]}')
    return int(found.group(1))


def main():
    if len(sys.argv) < 2:
        raise SystemExit('usage: round21_crossover.py WORKER [fixtures]')
    worker = str(Path(sys.argv[1]).resolve())
    # Sixty-four, pooled over two independent seed streams, and the number was
    # arrived at the hard way. Three fixtures read n23a1 at 0.882 and then
    # 1.097 on a second draw. Sixteen read n37a0 at 1.782 [1.394, 2.278] on one
    # stream and 1.116 [0.578, 2.154] on another -- the second band contains
    # one, so sixteen does not answer the question the script exists to ask.
    # n37a0's per-case log spread is 0.87, nearly triple n23a1's 0.30, which is
    # round 0019's finding that the spread grows with the cell, holding two
    # decades of r further out than round 0019 could reach. Sixty-four brings
    # both standard errors under a tenth and both bands clear of one.
    fixtures = int(sys.argv[2]) if len(sys.argv) > 2 else 64
    env = T.child_env()

    probe = run(worker, {'mode': 'fixture', 'degree': 37, 'curve_a': 0, 'target_seeds': [1],
                         'algorithm_seed': 1, 'config': CONFIG,
                         'factor_base': {'kind': 'subgroup_orbits', 'seed': 43, 'points': 222}}, env)
    if probe.get('status') != 'fixture':
        raise SystemExit(f'this worker refuses degree 37 ({probe.get("reason")!r}); '
                         'apply round21-wide-pair-table.patch and rebuild')

    print(f'instructions over {fixtures} fixtures a configuration, every IC report '
          f'checked by oracle.py\n')
    print(f'{"cell":8} {"r":>14} {"orbits":>7} {"IC/rho":>8} {"sd(log)":>8} {"95% band":>18} {"ok":>6}')
    table = {}
    for label, degree, curve_a, r, point_counts in CELLS:
        for points in point_counts:
            ic_ir, rho_ir, verified, orbits = [], [], 0, 0
            # Two streams rather than one: a single stream's draw is what made
            # the sixteen-fixture readings disagree with each other.
            draws = [(s, fixtures // len(SEEDS)) for s in SEEDS]
            for stream, count in draws:
              rng = random.Random(stream)
              for _ in range(count):
                job = {'degree': degree, 'curve_a': curve_a,
                       'target_seeds': [rng.getrandbits(64)],
                       'algorithm_seed': rng.getrandbits(64), 'config': CONFIG,
                       'factor_base': {'kind': 'subgroup_orbits', 'seed': 43, 'points': points}}
                fixture = run(worker, dict(job, mode='fixture'), env)
                if fixture.get('status') != 'fixture':
                    continue
                report = run(worker, dict(job, mode='ic'), env)
                if report.get('status') != 'complete':
                    continue
                orbits = len(report.get('factor_base_orbits') or [])
                verify(report, fixture['fixture'], expected_mode='ic', summands=3)
                verified += 1
                ic_ir.append(instructions(worker, dict(job, mode='ic'), env))
                rho_ir.append(instructions(worker, dict(job, mode='rho'), env))
            if not ic_ir:
                print(f'{label:8} {r:>14,} {orbits:>7} {"did not complete":>39}')
                continue
            logs = [math.log(a / b) for a, b in zip(ic_ir, rho_ir)]
            ratio, spread = math.exp(st.mean(logs)), st.stdev(logs) if len(logs) > 1 else 0.0
            err = spread / math.sqrt(len(logs))
            band = f'[{ratio * math.exp(-Z * err):.3f}, {ratio * math.exp(Z * err):.3f}]'
            table.setdefault(label, []).append((orbits, ratio, band))
            print(f'{label:8} {r:>14,} {orbits:>7} {ratio:>8.3f} {spread:>8.3f} {band:>18} '
                  f'{str(verified) + "/" + str(fixtures):>6}')

    print('\nbest achievable IC/rho at each cell, and what it says about the boundary')
    best = {c: min(v, key=lambda x: x[1]) for c, v in table.items()}
    for label, (orbits, ratio, band) in best.items():
        rr = dict((c, r) for c, _d, _a, r, _p in CELLS)[label]
        verdict = 'IC ahead' if ratio < 1 else 'RHO AHEAD'
        print(f'  {label:8} r = {rr:>14,}  {ratio:6.3f} {band} at {orbits} orbits   {verdict}')
    if len(best) == 2:
        lo, hi = CELLS[0][0], CELLS[1][0]
        if best[lo][1] < 1 < best[hi][1]:
            grew = CELLS[1][3] / CELLS[0][3]
            exponent = math.log(best[hi][1] / best[lo][1]) / math.log(grew)
            print(f'\nThe collector crosses rho between r = {CELLS[0][3]:,} and '
                  f'r = {CELLS[1][3]:,}.')
            print('Measured, on cells where both arms complete and every IC answer is')
            print('certified -- not extrapolated from the same-degree contrast.')
            print(f'\nThe ratio grows as r^{exponent:.3f} between the two ({grew:.0f}x in r,'
                  f' {best[hi][1] / best[lo][1]:.2f}x in the ratio).')
            print('round19_model.py derives r^(1/6) = r^0.167 for the balanced optimum.')
            print('Two cells fix one exponent and no curvature, so this is a rate between')
            print('two points, not a fit -- the agreement is worth no more than that.')


if __name__ == '__main__':
    main()
