#!/usr/bin/env python3
"""What does `max_trials` actually do to a measurement?

    python3 campaign_20260916/round22_budget_effect.py WORKER [reps]

**This script was written as an inertness gate, and the gate failed.**  It is
kept as the record of that, because the failure is the useful result.

Round 0021's `n37a0` ratio is upper-biased: 4 of its 64 rho runs were cut off
at `max_trials = 4096` and charged their truncated instruction counts (see the
correction banner on `round21_crossover.py`).  The obvious repair is to raise
the cap -- the worker accepts up to 65,536 -- and the argument for doing it
freely was:

    a run that COMPLETED under the old cap takes the same path, returns the
    same answer and costs the same instructions under the new one, because a
    cap that is never reached cannot be observed,

which would have made the raise an additive amendment, leaving every published
measurement readable unchanged.  **That prediction is false, and this script is
what falsified it.**  Answers and factor-base hashes are identical at every
cell, and the instruction count is not -- at every cell, on BOTH arms.  The
source says why:

    max_iterations_per_restart: cfg.max_trials,      // examples/ic_tournament_worker.rs

`max_trials` is rho's iterations-per-RESTART, so raising it does not merely
stop truncating: it changes how often rho abandons a walk and starts over,
which is a parameter of the algorithm being measured.  On the candidate's side
it is passed to `TinyIc::new` at construction.  It is a protocol parameter, not
a safety valve, and no amount of it being "never reached" makes it invisible.

So the cap cannot be quietly raised, and this script's job is the other one:
say how much it is worth, per cell and per arm, so that a ladder measured at
one cap can be honestly related to a panel measured at another.  The quantity
that matters is not the instruction count but the RATIO, since a cost that
lands on both arms largely cancels there -- and that, too, is measured rather
than asserted.

A fixture that did not complete at the old cap is counted and skipped: there is
nothing to compare, and it is the reason the question arose.
"""
import hashlib
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

FROZEN = json.loads((ROOT / 'runs/round-0020/candidates.json').read_text())[0]['config']
RAISED = dict(FROZEN, max_trials=65536)
CELLS = [('n13a0', 13, 0, 48), ('n17a1', 17, 1, 96), ('n19a0', 19, 0, 108),
         ('n19a1', 19, 1, 108), ('n23a0', 23, 0, 138), ('n23a1', 23, 1, 138),
         ('n29a1', 29, 1, 174), ('n31a0', 31, 0, 186),
         ('n37a0', 37, 0, 222), ('n43a1', 43, 1, 222)]


def run(worker, job, env, timeout=1800):
    done = subprocess.run([worker], input=json.dumps(job), text=True,
                          capture_output=True, env=env, timeout=timeout)
    return json.loads(done.stdout) if done.stdout else {}


def instructions(worker, job, env, timeout=3600):
    done = subprocess.run(
        ['valgrind', '--tool=callgrind', '--callgrind-out-file=/dev/null', worker],
        input=json.dumps(job), text=True, capture_output=True, env=env, timeout=timeout)
    found = re.search(r'Collected\s*:\s*(\d+)', done.stderr)
    if not found:
        raise RuntimeError(f'no Ir collected: {done.stderr[-300:]}')
    return int(found.group(1))


def answer(report):
    solutions = sorted(report.get('solutions', []), key=lambda s: s['index'])
    base = report.get('factor_base_orbits') or report.get('factor_base') or []
    return (tuple(s['recovered'] for s in solutions),
            hashlib.sha256(json.dumps(base, separators=(',', ':')).encode()).hexdigest())


def main():
    if len(sys.argv) < 2:
        raise SystemExit('usage: round22_budget_effect.py WORKER [reps]')
    worker = str(Path(sys.argv[1]).resolve())
    reps = int(sys.argv[2]) if len(sys.argv) > 2 else 3
    env = T.child_env()
    rng = random.Random(20260922)

    print(f'max_trials {FROZEN["max_trials"]} vs {RAISED["max_trials"]}, '
          f'{reps} fixtures a cell, identical inputs\n')
    print(f'{"cell":8} {"same answer":>12} {"Ir raised/frozen":>18} {"Ir raised/frozen":>18} '
          f'{"IC/rho frozen":>14} {"IC/rho raised":>14} {"skipped":>8}')
    print(f'{"":8} {"":12} {"IC":>18} {"rho":>18} {"":14} {"":14} {"":8}')
    moved_answer = []
    for label, degree, curve_a, points in CELLS:
        shifts = {'ic': [], 'rho': []}
        ratios = {'frozen': [], 'raised': []}
        same_answer = compared = skipped = 0
        for _ in range(reps):
            targets, algo = rng.getrandbits(64), rng.getrandbits(64)
            base = {'degree': degree, 'curve_a': curve_a, 'target_seeds': [targets],
                    'algorithm_seed': algo,
                    'factor_base': {'kind': 'subgroup_orbits', 'seed': 43, 'points': points}}
            # Both arms must complete at BOTH caps, or there is nothing to compare.
            reports = {(m, tag): run(worker, dict(base, mode=m, config=cfg), env)
                       for m in ('ic', 'rho') for tag, cfg in (('frozen', FROZEN), ('raised', RAISED))}
            if any(r.get('status') != 'complete' for r in reports.values()):
                skipped += 1
                continue
            compared += 1
            if all(answer(reports[(m, 'frozen')]) == answer(reports[(m, 'raised')])
                   for m in ('ic', 'rho')):
                same_answer += 1
            else:
                moved_answer.append(label)
            ir = {(m, tag): instructions(worker, dict(base, mode=m, config=cfg), env)
                  for m in ('ic', 'rho') for tag, cfg in (('frozen', FROZEN), ('raised', RAISED))}
            for m in ('ic', 'rho'):
                shifts[m].append(ir[(m, 'raised')] / ir[(m, 'frozen')])
            for tag in ('frozen', 'raised'):
                ratios[tag].append(ir[('ic', tag)] / ir[('rho', tag)])
        if not compared:
            print(f'{label:8} {"-":>12} {"-":>18} {"-":>18} {"-":>14} {"-":>14} {skipped:>8}')
            continue
        gm = lambda v: math.exp(st.mean([math.log(x) for x in v]))
        print(f'{label:8} {f"{same_answer}/{compared}":>12} {gm(shifts["ic"]):>18.4f} '
              f'{gm(shifts["rho"]):>18.4f} {gm(ratios["frozen"]):>14.4f} '
              f'{gm(ratios["raised"]):>14.4f} {skipped:>8}')

    print('\nThe answer columns are the part that held: raising the cap never changed a')
    print('recovered logarithm or a factor base. The instruction columns are the part')
    print('that did not, and they are why the raise is a change of protocol rather')
    print('than an extension of one.')
    print('\nThe last two columns are what a cross-cap comparison actually rests on: a')
    print('cost that lands on both arms cancels in the ratio, and how nearly it')
    print('cancels is a measurement, not an argument. Read them before relating any')
    print('ladder measured at one cap to a panel measured at another.')
    if moved_answer:
        print(f'\nCELLS WHERE THE CAP MOVED AN ANSWER: {sorted(set(moved_answer))}')
        print('That is not a budget effect and not admissible as one. Stop here.')
        raise SystemExit(1)


if __name__ == '__main__':
    main()
