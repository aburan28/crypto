#!/usr/bin/env python3
"""Census every curve cell the worker admits, with its subgroup order.

    python3 campaign_20260916/round16_cell_census.py [--worker runs/round-0015/worker]

The panel has been five cells since round 0006 and was never chosen; it was
inherited. This asks the frozen worker for a fixture at every odd degree it
accepts (5..31) and both `curve_a` values, and keeps the pairs that yield a
usable subgroup. A pair that does not is not a failure to record as evidence --
the curve simply has no subgroup this collector can work in.
"""
import argparse
import json
import random
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import tournament as T  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--worker', default='runs/round-0015/worker')
    ap.add_argument('--seed', type=int, default=20260916)
    args = ap.parse_args()
    rng = random.Random(args.seed)
    usable = []
    for degree in range(5, 32, 2):
        for a in (0, 1):
            job = {'mode': 'fixture', 'degree': degree, 'curve_a': a,
                   'target_seeds': [rng.getrandbits(64)],
                   'algorithm_seed': rng.getrandbits(64),
                   'factor_base': {'kind': 'subgroup_orbits', 'seed': 43,
                                   'points': 6 * degree},
                   'config': T.BASE_CONFIG}
            done = subprocess.run([args.worker], input=json.dumps(job), text=True,
                                  capture_output=True, timeout=60)
            report = json.loads(done.stdout)
            if report.get('status') != 'fixture':
                continue
            usable.append((degree, a, int(report['fixture']['subgroup_order'])))
    print(f'{"cell":8} {"subgroup order":>16}')
    for degree, a, order in usable:
        print(f'n{degree}a{a:<6} {order:>16,}')
    print(f'\n{len(usable)} usable of {len(range(5, 32, 2)) * 2} (degree, curve_a) pairs')


if __name__ == '__main__':
    main()
