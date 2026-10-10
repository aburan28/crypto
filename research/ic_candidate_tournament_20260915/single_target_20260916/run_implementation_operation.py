"""Bounded preflight, fixed-support tournament, audit and publication controller."""
import json
import os
from pathlib import Path
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
ROUND = ROOT/'runs/round-0008-single-implementation'


def status(phase, **extra):
    (HERE/'implementation-operation-status.json').write_text(json.dumps(
        dict(phase=phase, pid=os.getpid(), updated_unix=time.time(), **extra), indent=2)+'\n')


def run(phase, command, timeout, cwd=None, env=None):
    status(phase)
    with (HERE/('implementation-'+phase+'.log')).open('w') as log:
        subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, timeout=timeout,
                       check=True, cwd=cwd, env=env)


def main():
    os.sched_setaffinity(0, {0, 1, 7})
    parent = json.loads((HERE/'policy-operation-status.json').read_text())
    assert parent['phase'] == 'audited' and parent['decision']['winner'] == 'orbits2'
    run('parent-report', ['taskset', '-c', '0,1', sys.executable, str(ROOT/'report.py'),
        '--round', str(ROOT/'runs/round-0007-single-policy'), '--scoreboard',
        str(ROOT.parents[1]/'docs/index-calculus-scoreboard.html')], 1800)
    env = dict(os.environ, CARGO_TARGET_DIR=str(HERE/'preflight-build'))
    run('preflight', ['taskset', '-c', '0,1', 'cargo', 'test', '--offline', '--lib',
                     'single_target_', '--', '--nocapture'], 1800,
        cwd=HERE/'implementation-sources/combined', env=env)
    run('prepare', [sys.executable, str(HERE/'prepare_policy.py'), 'prepare',
        '--out', str(ROUND), '--source-root', str(HERE/'policy-sources/policy'),
        '--candidates', str(HERE/'implementation-candidates.json'), '--seed', '2026091608',
        '--cpu', '7', '--targets', '1', '--require-native-progress',
        '--comparison-kind', 'fixed-support'], 1800)
    run('tournament', [sys.executable, str(ROUND/'evaluator/tournament.py'),
                       'run', '--round', str(ROUND)], 7200)
    run('report', ['taskset', '-c', '0,1', sys.executable, str(ROOT/'report.py'),
        '--round', str(ROUND), '--scoreboard',
        str(ROOT.parents[1]/'docs/index-calculus-scoreboard.html')], 1800)
    status('audited', decision=json.loads((ROUND/'decision.json').read_text()))


if __name__ == '__main__':
    try:
        main()
    except Exception as error:
        status('failed', error=repr(error))
        raise
