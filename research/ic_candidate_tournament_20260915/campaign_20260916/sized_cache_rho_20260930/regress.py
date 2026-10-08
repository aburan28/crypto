#!/usr/bin/env python3
"""Regression of the sized-cache rho against round 0025's lean rho (PREREGISTRATION.md section 4).

Runs the candidate worker in rho mode on every fixture job of round 0025 that has a recorded lean-rho
answer, and compares the whole `solutions` array with the recorded one.  Not timed.

    python3 regress.py --round-dir ../../runs/round-0025 --worker WORKER --out regression.json
"""

import argparse
import json
import subprocess
from pathlib import Path


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--round-dir', required=True)
    p.add_argument('--worker', required=True)
    p.add_argument('--out', required=True)
    args = p.parse_args()
    rd = Path(args.round_dir).resolve()
    rows = []
    for stage in ('development', 'confirmation', 'replay', 'aa', 'smoke', 'selection'):
        for job in sorted((rd / 'runs' / stage).glob('*/rho/rep-*/job.json')):
            recorded = json.loads((job.parent / 'native/stdout.json').read_text())
            with job.open('rb') as stdin:
                run = subprocess.run([args.worker], stdin=stdin, capture_output=True, check=True)
            got = json.loads(run.stdout)
            rows.append({'stage': stage, 'case': job.parents[2].name, 'rep': job.parent.name,
                         'degree': recorded['fixture']['degree'],
                         'same_solutions': got['solutions'] == recorded['solutions'],
                         'same_recovered': [s['recovered'] for s in got['solutions']]
                         == [s['recovered'] for s in recorded['solutions']],
                         'all_verified': all(s['verified'] for s in got['solutions']),
                         'iterations': [s['iterations'] for s in got['solutions']],
                         'lean_iterations': [s['iterations'] for s in recorded['solutions']],
                         'restarts': [s['restarts'] for s in got['solutions']],
                         'lean_restarts': [s['restarts'] for s in recorded['solutions']]})
    Path(args.out).write_text(json.dumps(rows, indent=1))
    n = len(rows)
    print(f'{n} rho runs; same solutions {sum(r["same_solutions"] for r in rows)}; '
          f'same recovered log {sum(r["same_recovered"] for r in rows)}; all verified {sum(r["all_verified"] for r in rows)}')
    for r in rows:
        if not r['same_solutions']:
            print('DIFF', r['stage'], r['case'], r['rep'], 'n', r['degree'], r['iterations'], 'vs', r['lean_iterations'],
                  'restarts', r['restarts'], 'vs', r['lean_restarts'])


if __name__ == '__main__':
    main()
