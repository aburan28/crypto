#!/usr/bin/env python3
"""Is there a third certified cell between n37a0 and n41a0?

    python3 campaign_20260916/round22_cell_census.py WORKER

Round 0021 measured the crossover between two cells: `n23a1` at
r = 4,196,903 (IC/rho 0.831) and `n37a0` at r = 230,603,167 (1.533).  Two
points fix a rate and no curvature, and the two cells differ in DEGREE as well
as in `r`, so the r^0.153 that came out of them cannot separate the two
variables.  A third certified cell would.

`next-proposal-single-target.json` records that "degrees 37 and 41 are the only
primes the sampler admits above 31", and concludes a third cell needs either a
protocol amendment or a new curve family.  **That sentence was written against a
degree bound of 41** -- the probe build PROBE-degree-ceiling.md censused with,
which is not the bound that exists now.  `round21-wide-pair-table.patch` lifted
the ceiling to 61, and 43, 47, 53, 59 and 61 are prime and have never been asked
for a fixture.  This script asks.

A cell is usable for the crossover measurement only if all four hold:

  * the sampler returns a fixture at all (most (degree, a) pairs do not -- the
    curve simply has no subgroup this collector can work in, which is
    mathematics and not a guard);
  * the degree is PRIME, so `koblitz_tiny_ic::supports()` takes it and the
    report carries a descent `relation`.  A composite degree falls back to the
    general sparse path, whose output `oracle.py` rejects as "not certified as
    index calculus" -- that is exactly what made PROBE-degree-ceiling.md's
    5.46 and 5.81 inadmissible;
  * IC completes, and the report VERIFIES against the independent checker;
  * **rho completes too.**  This is the one that killed `n41a0`: at
    r = 5.5e11 rho returns in 0.06s against the ~72,600 steps it needs, so the
    frozen `max_trials` is cutting it off rather than IC out-running it, and an
    arm that does not complete has no cost.  A cell where rho is cut off yields
    no comparison at any sample size.

Nothing here is a measurement of anything.  It is a census: which cells exist,
and which of them could carry a measurement.  The ratios it prints are one
fixture each and are diagnostics for deciding where to spend, never evidence.
"""
import argparse
import json
import math
import random
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import tournament as T  # noqa: E402
from oracle import verify, InvalidEvidence  # noqa: E402

CONFIG = json.loads((ROOT / 'runs/round-0020/candidates.json').read_text())[0]['config']
# The window a third cell has to land in to add curvature: above `n37a0`, and
# below the order at which rho stops completing.  `n41a0` sits above it.
WINDOW = (230_603_167, 549_756_390_943)


def prime(n):
    return n > 1 and all(n % d for d in range(2, math.isqrt(n) + 1))


def run(worker, job, timeout=600):
    done = subprocess.run([worker], input=json.dumps(job), text=True,
                          capture_output=True, env=T.child_env(), timeout=timeout)
    return json.loads(done.stdout) if done.stdout else {}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('worker')
    ap.add_argument('--lo', type=int, default=33)
    ap.add_argument('--hi', type=int, default=61)
    ap.add_argument('--points', type=int, default=222)
    ap.add_argument('--seed', type=int, default=20260922)
    args = ap.parse_args()
    worker = str(Path(args.worker).resolve())
    rng = random.Random(args.seed)

    print(f'{"cell":8} {"r":>16} {"cofactor":>10} {"deg prime":>10} '
          f'{"IC":>12} {"rho":>12} {"in window":>10}')
    rows = []
    for degree in range(args.lo | 1, args.hi + 1, 2):
        for a in (0, 1):
            job = {'degree': degree, 'curve_a': a,
                   'target_seeds': [rng.getrandbits(64)],
                   'algorithm_seed': rng.getrandbits(64), 'config': CONFIG,
                   'factor_base': {'kind': 'subgroup_orbits', 'seed': 43,
                                   'points': args.points}}
            fixture = run(worker, dict(job, mode='fixture'))
            if fixture.get('status') != 'fixture':
                continue
            f = fixture['fixture']
            r, h = int(f['subgroup_order']), int(f['cofactor'])
            ic = run(worker, dict(job, mode='ic'))
            rho = run(worker, dict(job, mode='rho'))
            # "complete" is the worker's own word for it; the checker then says
            # whether what completed is admissible as index calculus.
            ic_state = ic.get('status', 'no output')
            if ic_state == 'complete':
                try:
                    verify(ic, f, expected_mode='ic', summands=3)
                    ic_state = 'certified'
                except InvalidEvidence as e:
                    ic_state = f'rejected'
                    if 'not certified as index calculus' in str(e):
                        ic_state = 'uncertified'
            rho_state = rho.get('status', 'no output')
            if rho_state == 'complete':
                try:
                    verify(rho, f, expected_mode='rho')
                except InvalidEvidence:
                    rho_state = 'rejected'
            window = WINDOW[0] < r < WINDOW[1]
            rows.append((f'n{degree}a{a}', r, h, prime(degree), ic_state, rho_state, window))
            print(f'n{degree}a{a:<6} {r:>16,} {h:>10,} {str(prime(degree)):>10} '
                  f'{ic_state:>12} {rho_state:>12} {str(window):>10}')

    certified = [x for x in rows if x[3] and x[4] == 'certified' and x[6]]
    usable = [x for x in certified if x[5] == 'complete']
    print(f'\n{len(rows)} cells exist in degrees {args.lo}..{args.hi}; '
          f'{len(certified)} carry a certified IC path inside the window; '
          f'{len(usable)} also complete rho at max_trials={CONFIG["max_trials"]}.')
    for cell, r, *_ in certified:
        blocked = '' if (cell, r) in [(c, rr) for c, rr, *_ in usable] else '   rho cut off'
        print(f'  {cell:8} r = {r:>16,}{blocked}')
    if certified and not usable:
        # The census's own answer to the question it was written to ask, and it
        # is not the one expected.  New prime cells above 41 DO exist and IC
        # certifies at them; what blocks every one is the same frozen
        # max_trials that `n41a0` hit.  rho needs roughly sqrt(pi*r/2A) steps:
        # about 2,200 at `n37a0`, which is under 4096 and is exactly why that
        # cell looked like the boundary, and about 9,200 at `n43a1`.
        #
        # That cap is a config value, not a ceiling.  The worker accepts
        # max_trials up to 65,536 ("invalid collection limits" above that), so
        # the amendment the next-proposal calls for is available without
        # touching the collector -- and raising a cap can only ever help rho,
        # so it cannot flatter the candidate.
        print('\n  Every one is blocked by the same thing: rho is cut off by')
        print(f'  max_trials={CONFIG["max_trials"]}, which the worker itself would allow up to')
        print('  65,536. The block is a budget, not mathematics -- see')
        print('  round22_budget_inertness.py for whether raising it is admissible.')
    print('\nOne fixture a cell. Nothing here is evidence about either algorithm;')
    print('rho not completing is the frozen trial budget, never a cost.')


if __name__ == '__main__':
    main()
