#!/usr/bin/env python3
"""Cache/branch simulation of both arms (PROTOCOL.md).  python3 run.py --round-dir ../../runs/round-0025 --out runs/grid"""

import argparse
import json
import subprocess
from pathlib import Path

CELLS = ['n13a0', 'n17a1', 'n19a0', 'n19a1', 'n23a0']
ARMS = {'ic': 'incumbent', 'rho': 'rho'}
EVENTS = ['Ir', 'I1mr', 'ILmr', 'Dr', 'D1mr', 'DLmr', 'Dw', 'D1mw', 'DLmw', 'Bc', 'Bcm', 'Bi', 'Bim']


def parse(path: Path):
    events, totals, per_fn, fn = None, None, {}, None
    for line in path.read_text().splitlines():
        if line.startswith('events:'):
            events = line.split()[1:]
        elif line.startswith('fn='):
            fn = line[3:]
        elif line.startswith('summary:'):
            totals = dict(zip(events, map(int, line.split()[1:])))
        elif fn and line[:1].isdigit():
            vals = list(map(int, line.split()[1:]))
            acc = per_fn.setdefault(fn, [0] * len(events))
            for i, v in enumerate(vals):
                acc[i] += v
    top = sorted(per_fn.items(), key=lambda kv: -kv[1][events.index('Bcm')])[:5]
    return totals, [{'fn': f, **dict(zip(events, v))} for f, v in top]


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--round-dir', required=True)
    p.add_argument('--out', required=True)
    args = p.parse_args()
    rd = Path(args.round_dir).resolve()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=False)
    rows = []
    for cell in CELLS:
        for idx in range(4):
            case = f'{cell}-{idx:03d}'
            for arm, src in ARMS.items():
                job = rd / 'runs/confirmation' / case / src / 'rep-0' / 'job.json'
                cg = out / f'{case}-{arm}.cg'
                with job.open('rb') as stdin:
                    subprocess.run(['valgrind', '--tool=cachegrind', '--cache-sim=yes', '--branch-sim=yes',
                                    f'--cachegrind-out-file={cg}', str(rd / 'worker')],
                                   stdin=stdin, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
                totals, top = parse(cg)
                rows.append({'cell': cell, 'case': case, 'arm': arm, 'totals': totals, 'top_mispredict_fns': top})
                cg.unlink()
    (out / 'counts.json').write_text(json.dumps(rows, indent=1))


if __name__ == '__main__':
    main()
