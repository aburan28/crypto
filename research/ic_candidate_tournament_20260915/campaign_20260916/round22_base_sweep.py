#!/usr/bin/env python3
"""Is the candidate handicapped by its default factor base at the new cells?

    python3 campaign_20260916/round22_base_sweep.py WORKER [fixtures]

`round22_crossover.py` measures every cell at the base the sampler builds by
default -- one batch, 8 orbits -- on round 0021's finding that bigger is
monotonically WORSE, which round 0019 found at all eight panel cells and round
0021 reproduced at `n37a0` up to 56 orbits.

**That finding cannot simply be carried to `n43a1`.**  The pair-table cost
model puts the balanced optimum at `F = (c*#E*t/k)^(1/3)`, which GROWS with
`r`: `n43a1` is 20x `n37a0` in `r`, so the optimum sits about 2.7x higher.  A
cell whose optimum has moved past the default would be measured with the
candidate handicapped, and "index calculus loses by 15x here" would then be a
statement about a badly chosen base rather than about the method.

So the default is a choice that has to be re-earned at every new cell, and this
sweep is where it is earned or lost.  Note which way the risk runs: reporting
the default when a larger base is better OVERSTATES the candidate's loss, so
this check protects the honesty of a result that already goes against the
candidate -- which is exactly when it is easiest to skip.

The sampler builds orbits in batches of 8 and `points` only requests a target,
so the achieved orbit count is read back from the report and is what the table
reports.  Every IC report is checked by `oracle.py`; an uncertified report is
not a cheaper configuration, it is a different algorithm.
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
# Same point counts round 0021 swept at `n37a0`, extended upward: the optimum
# is expected to move up with `r`, so the sweep has to be able to find it there.
POINTS = (222, 1777, 2400, 3600, 7200, 14400)
CELLS = [('n37a0', 37, 0, 230_603_167), ('n43a1', 43, 1, 4_644_189_029)]


def run(worker, job, env, timeout=3600):
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
        raise SystemExit('usage: round22_base_sweep.py WORKER [fixtures]')
    worker = str(Path(sys.argv[1]).resolve())
    fixtures = int(sys.argv[2]) if len(sys.argv) > 2 else 8
    env = T.child_env()

    print(f'{fixtures} fixtures a configuration, max_trials={CONFIG["max_trials"]}, '
          f'every IC report checked by oracle.py\n')
    for label, degree, curve_a, r in CELLS:
        print(f'{label}  r = {r:,}')
        print(f'  {"points":>7} {"orbits":>7} {"IC/rho":>8} {"95% band":>18} {"complete":>10}')
        best = None
        for points in POINTS:
            ic_ir, rho_ir, orbits, dropped = [], [], 0, 0
            rng = random.Random(20260922)
            for _ in range(fixtures):
                job = {'degree': degree, 'curve_a': curve_a,
                       'target_seeds': [rng.getrandbits(64)],
                       'algorithm_seed': rng.getrandbits(64), 'config': CONFIG,
                       'factor_base': {'kind': 'subgroup_orbits', 'seed': 43,
                                       'points': points}}
                fixture = run(worker, dict(job, mode='fixture'), env)
                if fixture.get('status') != 'fixture':
                    dropped += 1
                    continue
                report = run(worker, dict(job, mode='ic'), env)
                rho = run(worker, dict(job, mode='rho'), env)
                if report.get('status') != 'complete' or rho.get('status') != 'complete':
                    dropped += 1
                    continue
                orbits = len(report.get('factor_base_orbits') or [])
                verify(report, fixture['fixture'], expected_mode='ic', summands=3)
                verify(rho, fixture['fixture'], expected_mode='rho')
                ic_ir.append(instructions(worker, dict(job, mode='ic'), env))
                rho_ir.append(instructions(worker, dict(job, mode='rho'), env))
            n = len(ic_ir)
            if n < 2:
                print(f'  {points:>7} {orbits:>7} {"no usable pair":>8} '
                      f'{"":>18} {f"{n}/{fixtures}":>10}')
                continue
            logs = [math.log(a / b) for a, b in zip(ic_ir, rho_ir)]
            ratio = math.exp(st.mean(logs))
            err = st.stdev(logs) / math.sqrt(n)
            lo, hi = ratio * math.exp(-1.959963985 * err), ratio * math.exp(1.959963985 * err)
            print(f'  {points:>7} {orbits:>7} {ratio:>8.3f} '
                  f'{f"[{lo:.3f}, {hi:.3f}]":>18} {f"{n}/{fixtures}":>10}')
            if best is None or ratio < best[1]:
                best = (orbits, ratio)
        if best:
            note = ('the default is the cheapest measured size, so reporting it is not a handicap'
                    if best[0] <= 8 else
                    f'THE DEFAULT IS NOT BEST HERE: {best[0]} orbits reads {best[1]:.3f}')
            print(f'  best: {best[0]} orbits at {best[1]:.3f} -- {note}\n')


if __name__ == '__main__':
    main()
