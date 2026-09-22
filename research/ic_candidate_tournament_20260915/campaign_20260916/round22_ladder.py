#!/usr/bin/env python3
"""The ladder with each cell at its own measured best base.

    python3 campaign_20260916/round22_ladder.py WORKER [fixtures]

`round22_crossover.py` measures all three cells at the sampler's default single
batch, because that is what round 0021 did and the corrected `n37a0` number has
to be comparable with the published one.  `round22_base_sweep.py` then showed
the default is NOT the best base at two of the three cells, and that reporting
it overstates the candidate's loss:

    cell     best orbits    default (8)    best
    n23a1              8          0.766   0.766     the default IS best here
    n37a0             32          1.782   1.491
    n43a1             48          5.260   1.661

So a ladder built on the default is measuring a handicap that widens with `r`,
and its exponent is not the method's.  This script measures each cell at the
base the sweep chose for it, which is the comparison the cost model is actually
about: `F = (c*#E*t/k)^(1/3)` is where the model says the base SHOULD be, and
`r^(1/6)` is the rate it predicts FOR A BASE THAT IS THERE.

Two things this does not do.  It does not re-choose the base per fixture -- the
size is fixed per cell from an independent sweep, so the draw measured here
cannot select its own configuration.  And it does not claim the swept size is
the true optimum: the sampler builds in batches of 8, so the grid is coarse and
the curve is flat near its minimum.  A size that is merely the best of six
measured is reported as that and not as an optimum.

Both arms must complete on a fixture for it to count, which is the check round
0021 omitted for rho.
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

FROZEN = json.loads((ROOT / 'runs/round-0020/candidates.json').read_text())[0]['config']
CONFIG = dict(FROZEN, max_trials=65536)
# (label, degree, curve_a, r, points, the orbit count those points build).
# The point counts come from round22_base_sweep.py and nowhere else.
CELLS = [('n23a1', 23, 1, 4_196_903, 222, 8),
         ('n37a0', 37, 0, 230_603_167, 1777, 32),
         ('n43a1', 43, 1, 4_644_189_029, 3600, 48)]
SEEDS = (20260922, 4242)
Z = 1.959963985


def run(worker, job, env, timeout=1800):
    done = subprocess.run([worker], input=json.dumps(job), text=True,
                          capture_output=True, env=env, timeout=timeout)
    return json.loads(done.stdout) if done.stdout else {}


def instructions(worker, job, env, timeout=7200):
    done = subprocess.run(
        ['valgrind', '--tool=callgrind', '--callgrind-out-file=/dev/null', worker],
        input=json.dumps(job), text=True, capture_output=True, env=env, timeout=timeout)
    found = re.search(r'Collected\s*:\s*(\d+)', done.stderr)
    if not found:
        raise RuntimeError(f'no Ir collected: {done.stderr[-300:]}')
    return int(found.group(1))


def main():
    if len(sys.argv) < 2:
        raise SystemExit('usage: round22_ladder.py WORKER [fixtures]')
    worker = str(Path(sys.argv[1]).resolve())
    fixtures = int(sys.argv[2]) if len(sys.argv) > 2 else 64
    env = T.child_env()

    print(f'{fixtures} fixtures a cell over seeds {SEEDS}, max_trials='
          f'{CONFIG["max_trials"]}, base fixed per cell by round22_base_sweep.py\n')
    print(f'{"cell":8} {"r":>16} {"orbits":>7} {"IC/rho":>8} {"sd(log)":>8} '
          f'{"95% band":>18} {"complete":>10}')
    table = []
    for label, degree, curve_a, r, points, want in CELLS:
        ic_ir, rho_ir, dropped, built = [], [], [], 0
        for stream in SEEDS:
            rng = random.Random(stream)
            for _ in range(fixtures // len(SEEDS)):
                job = {'degree': degree, 'curve_a': curve_a,
                       'target_seeds': [rng.getrandbits(64)],
                       'algorithm_seed': rng.getrandbits(64), 'config': CONFIG,
                       'factor_base': {'kind': 'subgroup_orbits', 'seed': 43,
                                       'points': points}}
                fixture = run(worker, dict(job, mode='fixture'), env)
                if fixture.get('status') != 'fixture':
                    dropped.append('no fixture')
                    continue
                report = run(worker, dict(job, mode='ic'), env)
                rho = run(worker, dict(job, mode='rho'), env)
                if report.get('status') != 'complete':
                    dropped.append('IC incomplete')
                    continue
                if rho.get('status') != 'complete':
                    dropped.append('rho incomplete')
                    continue
                built = len(report.get('factor_base_orbits') or [])
                verify(report, fixture['fixture'], expected_mode='ic', summands=3)
                verify(rho, fixture['fixture'], expected_mode='rho')
                ic_ir.append(instructions(worker, dict(job, mode='ic'), env))
                rho_ir.append(instructions(worker, dict(job, mode='rho'), env))
        n = len(ic_ir)
        if dropped or n < 2:
            reasons = {x: dropped.count(x) for x in sorted(set(dropped))}
            print(f'{label:8} {r:>16,} {built:>7}   NO RATIO -- {reasons} ({n} usable)')
            continue
        if built != want:
            # The sweep picked a size by the orbit count it achieved, so a
            # different count here means this is not the configuration that was
            # chosen and the number below would be mislabelled.
            print(f'{label:8} {r:>16,} {built:>7}   REFUSED -- sweep chose {want} orbits')
            continue
        logs = [math.log(a / b) for a, b in zip(ic_ir, rho_ir)]
        ratio, spread = math.exp(st.mean(logs)), st.stdev(logs)
        err = spread / math.sqrt(n)
        lo, hi = ratio * math.exp(-Z * err), ratio * math.exp(Z * err)
        table.append((label, r, built, ratio, lo, hi))
        print(f'{label:8} {r:>16,} {built:>7} {ratio:>8.3f} {spread:>8.3f} '
              f'{f"[{lo:.3f}, {hi:.3f}]":>18} {f"{n}/{fixtures}":>10}')

    if len(table) >= 2:
        print('\nrate between consecutive cells, each at its own best measured base')
        rates = []
        for (c0, r0, _o0, v0, *_), (c1, r1, _o1, v1, *_) in zip(table, table[1:]):
            e = math.log(v1 / v0) / math.log(r1 / r0)
            rates.append(e)
            print(f'  {c0} -> {c1}: r^{e:.3f} ({r1 / r0:.0f}x in r, {v1 / v0:.2f}x in the ratio)')
        print('round19_model.py derives r^(1/6) = r^0.167 for the balanced optimum,')
        print('which is the rate for a base AT that optimum -- which is what this')
        print('ladder holds and the default-base ladder does not.')
        if len(rates) >= 2:
            print(f'\nThe two rates are {rates[0]:.3f} and {rates[1]:.3f}. They are free to')
            print('disagree; whether they do is the first curvature this campaign can see.')
    print('\nEvery ratio is between two arms that both completed, at a base size fixed')
    print('from an independent sweep rather than chosen by this draw.')


if __name__ == '__main__':
    main()
