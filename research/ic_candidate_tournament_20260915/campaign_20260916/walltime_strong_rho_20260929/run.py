#!/usr/bin/env python3
"""Small-cell re-timing against main's strongest rho (see PREREGISTRATION.md).

Arms per case:
  incumbent, lean_rho  -- round 0025's worker on the case's recorded job.json
  strong32, strong1    -- examples/koblitz_rho_batch_ks_strong.rs rung 3 with
                          KIC_RHO_LANES=32 / 1, KIC_RHO_DP_BITS=2, one derived
                          target per process, batch seed = case number + 1
  strong32_aa          -- strong32 again (A/A)
Every run goes through tools/isolated_bench.py's pinned runner inside one
reservation of CPU 3; each answer is checked (worker: recorded solutions;
strong rho: its own all_verified flag).

    python3 run.py --round-dir ../../runs/round-0025 --strong BIN --out runs/<label> [--rounds 7]
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
ARMS = ['incumbent', 'lean_rho', 'strong32', 'strong1', 'strong32_aa']
# A Williams square for five arms needs ten rows; this is the standard
# construction (rows i and its reversal), so each arm follows every other
# equally often over the ten rows.
BASE = [0, 1, 4, 2, 3]
WILLIAMS = [[(b + i) % 5 for b in BASE] for i in range(5)]
WILLIAMS += [row[::-1] for row in WILLIAMS]


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument('--round-dir', required=True)
    p.add_argument('--strong', required=True)
    p.add_argument('--out', required=True)
    p.add_argument('--cpu', type=int, default=3)
    p.add_argument('--rounds', type=int, default=7)
    p.add_argument('--max-other-cpu', type=float, default=0.10)
    args = p.parse_args()
    args.settle, args.max_psi = 3.0, 5.0

    round_dir = Path(args.round_dir).resolve()
    strong = Path(args.strong).resolve()
    worker = round_dir / 'worker'
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=False)
    cases = []
    for cell in CELLS:
        cases += sorted(q.name for q in (round_dir / 'runs/confirmation').glob(f'{cell}-*'))
    expected = {}
    for case in cases:
        for arm, src in (('incumbent', 'incumbent'), ('lean_rho', 'rho')):
            trial = round_dir / 'runs/confirmation' / case / src / 'rep-0'
            expected[case, arm] = (trial / 'job.json',
                                   json.loads((trial / 'native/stdout.json').read_text())['solutions'])
    (out / 'plan.json').write_text(json.dumps({
        'cells': CELLS, 'arms': ARMS, 'cases': cases, 'rounds': args.rounds, 'cpu': args.cpu,
        'worker': {'path': str(worker), 'sha256': sha(worker)},
        'strong': {'path': str(strong), 'sha256': sha(strong), 'rung': 3, 'dp_bits': 2,
                   'lanes': {'strong32': 32, 'strong32_aa': 32, 'strong1': 1},
                   'targets': 'one derived target per process; batch seed = index of the case within its cell + 1'},
        'order': 'Williams square for five arms, row (round + case index) mod 10',
        'warmup': 'one untimed run of every arm on the first case of each cell'}, indent=1, sort_keys=True))

    cpus = {args.cpu}
    ib.check_cpus(cpus)
    records = (out / 'runs.jsonl').open('w')
    scratch = out / 'stdout.json'
    with ib.locked(ib.DEFAULT_LOCK, wait=False):
        me = {os.getpid()}
        pre = ib.preflight(args, me)
        moved = ib.evict(cpus, me)
        os.sched_setaffinity(0, os.sched_getaffinity(0) - cpus)
        (out / 'conditions.json').write_text(json.dumps(
            {'host': ib.host(), 'preflight': pre, 'threads_moved': len(moved),
             'left_on_reserved': ib.unmovable_on(cpus, me)}, indent=1, sort_keys=True))
        try:
            def once(case, arm, rnd, timed=True):
                cell, index = case.split('-')
                if arm in ('incumbent', 'lean_rho'):
                    job, answer = expected[case, arm]
                    with job.open('rb') as stdin, scratch.open('wb') as stdout:
                        r = ib.run_pinned(args, cpus, [str(worker)], stdin, stdout)
                    report = json.loads(scratch.read_text())
                    ok = report.get('solutions') == answer
                    inproc = report.get('elapsed_seconds')
                else:
                    n, a = cell[1:].split('a')
                    lanes = 1 if arm == 'strong1' else 32
                    cmd = ['env', '-i', f'KIC_RHO_LANES={lanes}', 'KIC_RHO_DP_BITS=2', 'KIC_RHO_RUNG=3',
                           str(strong), n, a, 'signed_frobenius', '1', str(int(index) + 1)]
                    with scratch.open('wb') as stdout:
                        r = ib.run_pinned(args, cpus, cmd, None, stdout)
                    summary = json.loads(scratch.read_text().strip().splitlines()[-1])
                    ok = bool(summary.get('all_verified'))
                    inproc = summary['in_process_ms'] / 1e3
                r.update(case=case, cell=cell, arm=arm, round=rnd, timed=timed,
                         in_process_seconds=inproc, answer_ok=ok)
                records.write(json.dumps(r, sort_keys=True) + '\n')
                if r['exit_status'] != 0 or not ok:
                    raise SystemExit(f'{case} {arm}: exit {r["exit_status"]}, answer ok {ok}')

            for cell in CELLS:
                first = next(c for c in cases if c.startswith(cell))
                for arm in ARMS:
                    once(first, arm, -1, timed=False)
            for rnd in range(args.rounds):
                for index, case in enumerate(cases):
                    for position in WILLIAMS[(rnd + index) % len(WILLIAMS)]:
                        once(case, ARMS[position], rnd)
        finally:
            ib.restore(moved)
            records.close()
            scratch.unlink(missing_ok=True)
    return 0


if __name__ == '__main__':
    sys.exit(main())
