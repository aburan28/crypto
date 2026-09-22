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

**The argmin of a sweep is not a base size, it is a selection.**  At 16
fixtures `n43a1` reads 1.894, 1.886 and 1.877 at 24, 32 and 48 orbits -- three
sizes that are statistically indistinguishable, whose bands overlap almost
entirely.  Taking the lowest of them and quoting one exponent would be six
chances to find a favourable number at the cell that drives the rate, and it
would bias the ladder toward the candidate exactly where the ladder is most
sensitive.  So every size inside the flat region is measured and the exponent
is reported as a RANGE across them.  If that range is wide, the exponent is not
resolved and saying so is the result.

It also does not re-choose the base per fixture: the size is fixed per cell
before the draw, so a fixture cannot select its own configuration.  And no
swept size is called an optimum -- the sampler builds in batches of 8, so the
grid is coarse.  Worth recording against the model: `F = (c*#E*t/k)^(1/3)`
anchored at `n23a1`'s 8 orbits predicts 30 at `n37a0`, which measures best at
32, and 83 at `n43a1`, where 88 orbits measures 2.870 against 1.877 at 48.  The
model's growth holds for the first step and breaks by the second.

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
# (label, degree, curve_a, r, [(points, orbits built), ...]).  The sizes come
# from round22_base_sweep.py, with the sizes that sweep SKIPPED filled in: its
# `points` grid jumped straight from 8 orbits to 24 and 32, so 16 was never
# measured at any cell, and at 64 fixtures `n43a1`'s minimum sat at the lowest
# size tested -- which is how a minimum outside the grid announces itself.
# The flat region is then selected from the measurements by a stated rule
# (every size whose band overlaps the minimum's), because picking the single
# argmin of a flat region is a selection rather than a measurement.
CELLS = [('n23a1', 23, 1, 4_196_903, [(222, 8), (400, 16), (800, 24)]),
         ('n37a0', 37, 0, 230_603_167,
          [(222, 8), (600, 16), (1400, 24), (1777, 32), (2400, 40)]),
         ('n43a1', 43, 1, 4_644_189_029,
          [(222, 8), (800, 16), (1400, 24), (2400, 32), (3600, 48)])]
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

    print(f'{fixtures} fixtures a configuration over seeds {SEEDS}, max_trials='
          f'{CONFIG["max_trials"]}; base sizes fixed by round22_base_sweep.py\n')
    print(f'{"cell":8} {"r":>16} {"orbits":>7} {"IC/rho":>8} {"sd(log)":>8} '
          f'{"95% band":>18} {"complete":>10}')
    table = {}
    for label, degree, curve_a, r, sizes in CELLS:
        for points, want in sizes:
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
                # The sweep chose a size by the orbit count it achieved, so a
                # different count here is a different configuration and the
                # number would be mislabelled.
                print(f'{label:8} {r:>16,} {built:>7}   REFUSED -- sweep chose {want}')
                continue
            logs = [math.log(a / b) for a, b in zip(ic_ir, rho_ir)]
            ratio, spread = math.exp(st.mean(logs)), st.stdev(logs)
            err = spread / math.sqrt(n)
            lo, hi = ratio * math.exp(-Z * err), ratio * math.exp(Z * err)
            table.setdefault(label, []).append((built, ratio, lo, hi))
            print(f'{label:8} {r:>16,} {built:>7} {ratio:>8.3f} {spread:>8.3f} '
                  f'{f"[{lo:.3f}, {hi:.3f}]":>18} {f"{n}/{fixtures}":>10}')

    # The flat region, by a rule fixed before the numbers: the minimum, plus
    # every size whose 95% band overlaps the minimum's. A size that is
    # resolvably worse is excluded; one that is not is carried, because the
    # ladder must not depend on picking between sizes the data cannot separate.
    flat = {}
    for label, rows in table.items():
        best = min(rows, key=lambda x: x[1])
        flat[label] = [x for x in rows if x[2] <= best[3] and best[2] <= x[3]]
        kept = ', '.join(f'{o} orbits {v:.3f}' for o, v, *_ in flat[label])
        edge = (rows[0][0] == best[0] or rows[-1][0] == best[0])
        print(f'\n{label}: flat region = {kept}'
              + ('   MINIMUM AT THE EDGE OF THE GRID -- it may lie outside' if edge else ''))
    table = flat
    order = [(c, r) for c, _d, _a, r, _s in CELLS if table.get(c)]
    if len(order) >= 2:
        print('\nrate between consecutive cells. Every size inside each cell\'s flat')
        print('region is carried, so each step is a RANGE and not a point.')
        rates = []
        for (c0, r0), (c1, r1) in zip(order, order[1:]):
            es = [math.log(v1 / v0) / math.log(r1 / r0)
                  for _o0, v0, *_ in table[c0] for _o1, v1, *_ in table[c1]]
            rates.append((min(es), max(es)))
            print(f'  {c0} -> {c1}: r^[{min(es):.3f}, {max(es):.3f}] '
                  f'over {len(table[c0])}x{len(table[c1])} base choices '
                  f'({r1 / r0:.0f}x in r)')
        print('\nround19_model.py derives r^(1/6) = r^0.167 for the balanced optimum,')
        print('which is the rate for a base AT that optimum. This ladder holds one;')
        print('a ladder pinned at the sampler default does not, which is why round')
        print('0021\'s exponent was measuring a widening handicap rather than a method.')
        if len(rates) >= 2:
            (a0, a1), (b0, b1) = rates[0], rates[1]
            overlap = not (a1 < b0 or b1 < a0)
            print(f'\nThe two steps are r^[{a0:.3f}, {a1:.3f}] and r^[{b0:.3f}, {b1:.3f}].')
            print('They OVERLAP, so this ladder sees no curvature: one rate fits both'
                  if overlap else
                  'They are DISJOINT, so the rate genuinely changes between the steps')
            print('steps, and the base choice inside each cell does not decide that.'
                  if overlap else 'and the base choice inside each cell cannot explain it.')
    print('\nEvery ratio is between two arms that both completed, at a base size fixed')
    print('from an independent sweep rather than chosen by this draw.')


if __name__ == '__main__':
    main()
