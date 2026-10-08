"""Bounded single-target tournament; preserve status and logs across chat turns."""
import json
import os
from pathlib import Path
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
ROUND = ROOT/'runs/round-0006-single'


def status(phase, **extra):
    (HERE/'operation-status.json').write_text(json.dumps(
        {'phase': phase, 'pid': os.getpid(), 'updated_unix': time.time(), **extra}, indent=2)+'\n')


def run(phase, command, timeout):
    status(phase)
    with (HERE/(phase+'.log')).open('w') as log:
        subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, timeout=timeout, check=True)


def main():
    os.sched_setaffinity(0, {0, 1})
    status('waiting_for_preflight')
    deadline = time.monotonic()+1200
    while time.monotonic() < deadline:
        log = (HERE/'preflight.log').read_text()
        if 'test result: ok. 1 passed' in log:
            break
        if 'test result: FAILED' in log or 'error: could not compile' in log:
            raise RuntimeError('preflight failed; preserve its log')
        time.sleep(3)
    else:
        raise TimeoutError('preflight did not complete within the bounded wait')
    run('prepare', [sys.executable, str(ROOT/'tournament.py'), 'prepare',
                   '--out', str(ROUND), '--source-root', str(ROOT/'runs/round-0004/source_candidates/folded_lift/source'),
                   '--candidates', str(HERE/'candidates.json'), '--seed', '2026091606',
                   '--cpu', '7', '--targets', '1', '--require-native-progress'], 1800)
    run('tournament', [sys.executable, str(ROUND/'evaluator/tournament.py'), 'run', '--round', str(ROUND)], 7200)
    run('audit', [sys.executable, str(ROUND/'evaluator/tournament.py'), 'verify', '--round', str(ROUND)], 1800)
    status('audited', decision=json.loads((ROUND/'decision.json').read_text()))


if __name__ == '__main__':
    try:
        main()
    except Exception as error:
        status('failed', error=repr(error))
        raise
