#!/usr/bin/env python3
"""Isolated native re-timing of round 0025's small cells (see PREREGISTRATION.md).

Every timed run goes through tools/isolated_bench.py's pinned runner inside one
reservation: one lock, one preflight, every other movable thread moved off the
benchmark CPU for the whole batch, and per-run contention counters.

    python3 run.py --round-dir ../../runs/round-0025 --out runs/<label> [--rounds 15]
"""

import argparse
import json
import os
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO / 'tools'))
import isolated_bench as ib  # noqa: E402

CELLS = ['n13a0', 'n17a1', 'n19a0', 'n19a1', 'n23a0']
ARMS = ['incumbent', 'switch', 'rho', 'rho_aa']
# Williams balanced Latin square for four arms: over the four rows every arm
# runs in every position once and directly follows every other arm exactly
# once.  A plain rotation would put rho_aa straight after rho three times in
# four, handing the A/A copy a warm binary, job file and predictors (the
# one-round smoke measured rho_aa/rho at 0.84-0.90 for exactly that reason).
WILLIAMS = [[0, 1, 3, 2], [1, 2, 0, 3], [2, 3, 1, 0], [3, 0, 2, 1]]
SOURCE_ARM = {'incumbent': 'incumbent', 'switch': 'switch', 'rho': 'rho', 'rho_aa': 'rho'}


def cases(round_dir: Path) -> list[str]:
    out = []
    for cell in CELLS:
        found = sorted(p.name for p in (round_dir / 'runs/confirmation').glob(f'{cell}-*'))
        out += found
    return out


def binary(round_dir: Path, arm: str) -> Path:
    candidates = {a['id']: a for a in json.loads((round_dir / 'candidates.json').read_text())}
    source = SOURCE_ARM[arm]
    relative = candidates[source]['binary_relative'] if source in candidates else 'worker'
    return round_dir / relative


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument('--round-dir', required=True)
    p.add_argument('--out', required=True)
    p.add_argument('--cpu', type=int, default=3)
    p.add_argument('--rounds', type=int, default=15)
    p.add_argument('--max-other-cpu', type=float, default=0.10)
    args = p.parse_args()
    args.settle, args.max_psi = 3.0, 5.0

    round_dir = Path(args.round_dir).resolve()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=False)
    case_ids = cases(round_dir)
    jobs, expected = {}, {}
    for case in case_ids:
        for arm in ARMS:
            trial = round_dir / 'runs/confirmation' / case / SOURCE_ARM[arm] / 'rep-0'
            jobs[case, arm] = trial / 'job.json'
            expected[case, arm] = json.loads((trial / 'native/stdout.json').read_text())['solutions']
    binaries = {arm: binary(round_dir, arm) for arm in ARMS}
    (out / 'plan.json').write_text(json.dumps({
        'cells': CELLS, 'arms': ARMS, 'cases': case_ids, 'rounds': args.rounds, 'cpu': args.cpu,
        'binaries': {a: {'path': str(b), 'sha256': ib_sha(b)} for a, b in binaries.items()},
        'order': 'per round, per case, Williams square row (round + case index) mod 4',
        'warmup': 'one untimed run of every arm on the first case of each cell'}, indent=1, sort_keys=True))

    cpus = {args.cpu}
    ib.check_cpus(cpus)
    records = (out / 'runs.jsonl').open('w')
    scratch = out / 'stdout.json'
    with ib.locked(ib.DEFAULT_LOCK, wait=False):
        me = {os.getpid()}
        pre = ib.preflight(args, me)
        moved = ib.evict(cpus, me)
        # Keep this driver off the benchmark CPU between runs; run_pinned
        # pins each child there for its spawn only.
        os.sched_setaffinity(0, os.sched_getaffinity(0) - cpus)
        header = {'host': ib.host(), 'preflight': pre, 'threads_moved': len(moved),
                  'left_on_reserved': ib.unmovable_on(cpus, me)}
        (out / 'conditions.json').write_text(json.dumps(header, indent=1, sort_keys=True))
        try:
            def once(case, arm, rnd, timed=True):
                with jobs[case, arm].open('rb') as stdin, scratch.open('wb') as stdout:
                    r = ib.run_pinned(args, cpus, [str(binaries[arm])], stdin, stdout)
                report = json.loads(scratch.read_text())
                r.update(case=case, cell=case.split('-')[0], arm=arm, round=rnd, timed=timed,
                         worker_elapsed_seconds=report.get('elapsed_seconds'),
                         answer_matches=report.get('solutions') == expected[case, arm])
                records.write(json.dumps(r, sort_keys=True) + '\n')
                if r['exit_status'] != 0 or not r['answer_matches']:
                    raise SystemExit(f'{case} {arm}: exit {r["exit_status"]}, answer match {r["answer_matches"]}')

            for cell in CELLS:
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


def ib_sha(path: Path) -> str:
    import hashlib
    return hashlib.sha256(path.read_bytes()).hexdigest()


if __name__ == '__main__':
    sys.exit(main())
