#!/usr/bin/env python3
"""Isolated native timing of the sized-cache rho against the IC arm (PREREGISTRATION.md).

Every timed run goes through tools/isolated_bench.py's pinned runner inside one reservation.

    python3 run.py --round-dir ../../runs/round-0025 --sized ../../runs/round-0025/sized/worker \
        --out runs/<label> [--rounds 20] [--cells n13a0,...]
"""

import argparse
import hashlib
import json
import os
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO / 'tools'))
import isolated_bench as ib  # noqa: E402

CELLS = ['n13a0', 'n17a1', 'n19a0', 'n19a1', 'n23a0']
ARMS = ['incumbent', 'rho', 'sized', 'sized_aa', 'inc_new']
# (binary, job of the round's recorded arm, how the answer is checked).  `sized` arms run the changed
# worker in rho mode; their walk may differ from the lean rho's by a few steps, so only the recovered
# logarithm and its verification are compared.  `inc_new` runs the changed worker's IC mode, whose
# answer must match exactly.
SPEC = {'incumbent': ('round', 'incumbent', 'exact'), 'rho': ('round', 'rho', 'exact'),
        'sized': ('sized', 'rho', 'recovered'), 'sized_aa': ('sized', 'rho', 'recovered'),
        'inc_new': ('sized', 'incumbent', 'exact')}
# A Williams square for five arms needs ten rows: the standard construction and its reversal.
BASE = [0, 1, 4, 2, 3]
WILLIAMS = [[(b + i) % 5 for b in BASE] for i in range(5)]
WILLIAMS += [row[::-1] for row in WILLIAMS]


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument('--round-dir', required=True)
    p.add_argument('--sized', required=True)
    p.add_argument('--out', required=True)
    p.add_argument('--cpu', type=int, default=3)
    p.add_argument('--rounds', type=int, default=20)
    p.add_argument('--cells', default=','.join(CELLS))
    p.add_argument('--max-other-cpu', type=float, default=0.10)
    args = p.parse_args()
    args.settle, args.max_psi = 3.0, 5.0

    cells = args.cells.split(',')
    round_dir = Path(args.round_dir).resolve()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=False)
    binaries = {'round': round_dir / 'worker', 'sized': Path(args.sized).resolve()}
    case_ids = [c for cell in cells for c in sorted(p.name for p in (round_dir / 'runs/confirmation').glob(f'{cell}-*'))]
    jobs, expected = {}, {}
    for case in case_ids:
        for arm, (_, source, _) in SPEC.items():
            trial = round_dir / 'runs/confirmation' / case / source / 'rep-0'
            jobs[case, arm] = trial / 'job.json'
            expected[case, arm] = json.loads((trial / 'native/stdout.json').read_text())['solutions']
    (out / 'plan.json').write_text(json.dumps({
        'cells': cells, 'arms': ARMS, 'spec': SPEC, 'cases': case_ids, 'rounds': args.rounds, 'cpu': args.cpu,
        'binaries': {k: {'path': str(v), 'sha256': sha(v)} for k, v in binaries.items()},
        'order': 'per round, per case, Williams row (round + case index) mod 10',
        'warmup': 'one untimed run of every arm on the first case of each cell'}, indent=1, sort_keys=True))

    def matches(case, arm, report):
        want, got = expected[case, arm], report.get('solutions')
        if SPEC[arm][2] == 'exact':
            return got == want
        return (got is not None and [s['recovered'] for s in got] == [s['recovered'] for s in want]
                and all(s['verified'] for s in got))

    cpus = {args.cpu}
    ib.check_cpus(cpus)
    records = (out / 'runs.jsonl').open('w')
    scratch = out / 'stdout.json'
    with ib.locked(ib.DEFAULT_LOCK, wait=False):
        me = {os.getpid()}
        pre = ib.preflight(args, me)
        moved = ib.evict(cpus, me)
        os.sched_setaffinity(0, os.sched_getaffinity(0) - cpus)
        header = {'host': ib.host(), 'preflight': pre, 'threads_moved': len(moved),
                  'left_on_reserved': ib.unmovable_on(cpus, me)}
        (out / 'conditions.json').write_text(json.dumps(header, indent=1, sort_keys=True))
        try:
            def once(case, arm, rnd, timed=True):
                with jobs[case, arm].open('rb') as stdin, scratch.open('wb') as stdout:
                    r = ib.run_pinned(args, cpus, [str(binaries[SPEC[arm][0]])], stdin, stdout)
                report = json.loads(scratch.read_text())
                r.update(case=case, cell=case.split('-')[0], arm=arm, round=rnd, timed=timed,
                         worker_elapsed_seconds=report.get('elapsed_seconds'),
                         answer_matches=matches(case, arm, report),
                         iterations=[s.get('iterations') for s in report.get('solutions', [])])
                records.write(json.dumps(r, sort_keys=True) + '\n')
                if r['exit_status'] != 0 or not r['answer_matches']:
                    raise SystemExit(f'{case} {arm}: exit {r["exit_status"]}, answer match {r["answer_matches"]}')

            for cell in cells:
                first = next(c for c in case_ids if c.startswith(cell))
                for arm in ARMS:
                    once(first, arm, -1, timed=False)
            for rnd in range(args.rounds):
                for index, case in enumerate(case_ids):
                    for position in WILLIAMS[(rnd + index) % len(WILLIAMS)]:
                        once(case, ARMS[position], rnd)
        finally:
            ib.restore(moved)
            records.close()
            scratch.unlink(missing_ok=True)
    return 0


if __name__ == '__main__':
    sys.exit(main())
