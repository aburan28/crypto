#!/usr/bin/env python3
"""Decompose a measured native wall into the segments that make it up.

    # Build the instrumented worker first, from this directory
    # (research/ic_candidate_tournament_20260915). `patch -d` keeps the patch
    # path relative to here while applying it over there, and cargo must run
    # inside the tree so the snapshot's .cargo/config.toml selects the musl
    # static target -- a --manifest-path build would read neither.
    cp -a runs/round-0015/source /tmp/probe
    patch -d /tmp/probe -p1 -i "$PWD/campaign_20260916/round16-floor-probe.patch"
    ( cd /tmp/probe && cargo build --release --example ic_tournament_worker )
    python3 campaign_20260916/round16_floor_probe.py --worker \
      /tmp/probe/target/x86_64-unknown-linux-musl/release/examples/ic_tournament_worker

The probe patch adds `CLOCK_REALTIME` marks inside the worker at the entry to
`main`, after the job is parsed, after the algorithm returns, after the report
bytes exist and after stdout is flushed, and prints them on stderr. This harness
takes the same clock on the evaluator side around `Popen` and around the reap,
so every segment of the wall the tournament charges is attributable.

The measurement this exists to settle: how much of the wall is the workload, how
much is genuine cold-process cost both arms pay, and how much is the evaluator's
own Python. Only the last is an artifact, and ROUND16-single-target.md records
why removing it would make the challenger gate harder rather than easier.
"""
import argparse
import collections
import json
import os
import resource
import statistics
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import tournament as T  # noqa: E402

# The child execs and can reach `main` while the evaluator is still inside
# Python, so parent-side and child-side spans OVERLAP and do not partition the
# wall. Only these two views are serial, and they are what gets reported:
#   in-main   t1..t5, split by the child's own clock
#   outside   wall - in-main, reported whole: Python, the spawn syscall, exec
#             and page-in, exit and reap, with their overlaps unresolved
IN_MAIN = [('input', 't1', 't2'), ('algorithm', 't2', 't3'),
           ('serialise', 't3', 't4'), ('write', 't4', 't5')]


def once(binary, job, scratch, cpu, memory):
    scratch.mkdir(parents=True, exist_ok=True)
    payload = json.dumps(job, sort_keys=True).encode()
    clock = lambda: time.clock_gettime(time.CLOCK_REALTIME)
    with (scratch / 'out.json').open('w') as out, (scratch / 'err.txt').open('w') as err:
        inherited = os.sched_getaffinity(0)
        os.sched_setaffinity(0, {cpu})
        marks = {'pre': clock()}
        try:
            process = subprocess.Popen([binary], stdin=subprocess.PIPE, stdout=out,
                                       stderr=err, env=T.child_env(), start_new_session=True)
        finally:
            os.sched_setaffinity(0, inherited)
        marks['spawned'] = clock()
        resource.prlimit(process.pid, resource.RLIMIT_AS, (memory, memory))
        resource.prlimit(process.pid, resource.RLIMIT_CORE, (0, 0))
        marks['capped'] = clock()
        process.communicate(payload)
        marks['end'] = clock()
    assert process.returncode == 0, process.returncode
    line = (scratch / 'err.txt').read_text().split('\n')[0].split()
    assert line and line[0] == 'MARK', 'worker is not the instrumented build'
    for name, value in zip(('t1', 't2', 't3', 't4', 't5'), line[1:6]):
        marks[name] = float(value)
    return marks


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--worker', required=True)
    ap.add_argument('--round', default='runs/round-0015')
    ap.add_argument('--cases', type=int, default=12)
    ap.add_argument('--repetitions', type=int, default=3)
    ap.add_argument('--cpu', type=int, default=3)
    args = ap.parse_args()

    root = Path(args.round)
    fixtures = json.loads((root / 'fixtures.json').read_text())['confirmation']
    config = json.loads((root / 'candidates.json').read_text())[0]['config']
    scratch = Path(os.environ.get('TMPDIR', '/tmp')) / 'round16-floor-probe'
    by_cell = collections.defaultdict(list)
    for case in fixtures:
        by_cell[case['cell']].append(case)

    rows = collections.defaultdict(list)
    for cell in sorted(by_cell):
        for mode in ('ic', 'rho'):
            for case in by_cell[cell][:args.cases]:
                job = dict(case['job'], mode=mode, config=config)
                for _ in range(args.repetitions):
                    rows[(cell, mode)].append(
                        once(args.worker, job, scratch, args.cpu, 8 * 1024 ** 3))

    print(f'{"cell":8} {"mode":5} {"wall":>9} {"in-main":>9} {"outside":>9} '
          f'{"algo share":>11}' + ''.join(f'{n:>11}' for n, _, _ in IN_MAIN))
    for key in sorted(rows):
        marks = rows[key]
        median = lambda a, b: statistics.median((m[b] - m[a]) * 1e6 for m in marks)
        wall, in_main = median('pre', 'end'), median('t1', 't5')
        line = (f'{key[0]:8} {key[1]:5} {wall:9.1f} {in_main:9.1f} '
                f'{wall - in_main:9.1f} {median("t2", "t3") / wall:11.3f}')
        line += ''.join(f'{median(a, b):11.1f}' for _, a, b in IN_MAIN)
        print(line)
    print('\nin-main and outside are serial and sum to the wall; the parts of '
          '"outside"\n(Python, spawn syscall, exec and page-in, exit, reap) '
          'overlap the child and are\nnot separable from this harness.')


if __name__ == '__main__':
    main()
