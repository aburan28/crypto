"""Run the separately declared single-target policy round after its parent ends."""
import json
import os
from pathlib import Path
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
PARENT = ROOT/'runs/round-0006b-single'
ROUND = ROOT/'runs/round-0007-single-policy'


def status(phase, **extra):
    (HERE/'policy-operation-status.json').write_text(json.dumps(
        {'phase': phase, 'pid': os.getpid(), 'updated_unix': time.time(), **extra}, indent=2)+'\n')


def run(phase, command, timeout):
    status(phase)
    with (HERE/('policy-'+phase+'.log')).open('w') as log:
        subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, timeout=timeout, check=True)


def main():
    os.sched_setaffinity(0, {0, 1, 7})
    status('waiting_for_parent')
    deadline = time.monotonic()+7200
    while time.monotonic() < deadline:
        parent = json.loads((HERE/'operation-status.json').read_text())
        if parent['phase'] == 'audited':
            break
        if parent['phase'] == 'failed':
            raise RuntimeError('parent failed; preserve this proposal without measuring it')
        time.sleep(5)
    else:
        raise TimeoutError('parent did not complete in the bounded wait')
    decision = json.loads((PARENT/'decision.json').read_text())
    if decision['winner'] != 'incumbent':
        raise RuntimeError('parent promoted a different baseline; rebase proposal before measuring')
    run('prepare', [sys.executable, str(HERE/'prepare_policy.py'), 'prepare',
                   '--out', str(ROUND), '--source-root', str(HERE/'policy-sources/incumbent'),
                   '--candidates', str(HERE/'policy-candidates.json'), '--seed', '2026091607',
                   '--cpu', '7', '--targets', '1', '--require-native-progress',
                   '--comparison-kind', 'factor-base-policy'], 1800)
    run('tournament', [sys.executable, str(ROUND/'evaluator/tournament.py'), 'run', '--round', str(ROUND)], 7200)
    run('audit', [sys.executable, str(ROUND/'evaluator/tournament.py'), 'verify', '--round', str(ROUND)], 1800)
    status('audited', decision=json.loads((ROUND/'decision.json').read_text()))


if __name__ == '__main__':
    try:
        main()
    except Exception as error:
        status('failed', error=repr(error))
        raise
