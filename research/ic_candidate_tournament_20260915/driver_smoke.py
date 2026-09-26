#!/usr/bin/env python3
"""Exercise the existing drivers with public points; no qualification or winner."""
import argparse
import json
from pathlib import Path
import subprocess
import sys

from oracle import require
from tournament import read, write

HERE = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prepared', required=True, type=Path)
    parser.add_argument('--out', required=True, type=Path)
    args = parser.parse_args()
    out, prepared = args.out.resolve(), args.prepared.resolve()
    out.mkdir(parents=True, exist_ok=False)
    config = dict(solver='pair_table', linear_algebra='tiny_gauss', summands=3,
                  batch_trials=1, max_trials=65536)
    registry = [dict(id='incumbent', config=config),
                dict(id='bounded', config=dict(config, max_trials=32768))]
    write(out/'candidates.json', registry, exclusive=True)
    commands = []

    def run(command):
        index = len(commands)
        commands.append([str(v) for v in command])
        write(out/'commands.json', commands)
        with (out/f'command-{index}.log').open('w') as log:
            subprocess.run(commands[-1], stdout=log, stderr=subprocess.STDOUT, check=True)

    native = out/'native'
    run([sys.executable, HERE/'autolab.py', 'screen', '--source-root', prepared/'source',
         '--out', native, '--candidates', out/'candidates.json', '--cells', '13a0,17a1',
         '--cases', '1', '--repetitions', '1', '--seed', '2026092530', '--timeout', '180'])
    run([sys.executable, native/'evaluator/autolab.py', 'verify', '--round', native])
    require(read(native/'summary.json')['all_jobs_verified'], 'native driver retained failures')

    failure_registry = [registry[0], dict(id='budget', config=dict(config, max_trials=1))]
    write(out/'failure-candidates.json', failure_registry, exclusive=True)
    failure = out/'native-failure'
    run([sys.executable, HERE/'autolab.py', 'screen', '--source-screen', native,
         '--out', failure, '--candidates', out/'failure-candidates.json', '--cells', '13a0',
         '--cases', '1', '--repetitions', '1', '--seed', '2026092530', '--timeout', '180'])
    run([sys.executable, failure/'evaluator/autolab.py', 'verify', '--round', failure])
    failed = read(failure/'trials/n13a0-000/budget/rep-0/run.json')
    require(failed['status'] == 'error' and failed['total_operations'] is None
            and failed['native_timing'] is None and failed['certificate'] is None,
            'expected collection exhaustion acquired a verified measurement')
    ids = [{read(p)['run_id'] for p in root.glob('trials/**/run.json')} for root in (native, failure)]
    require(len(ids[0]) == 12 and len(ids[1]) == 6 and not ids[0] & ids[1],
            'repeated workload has colliding run identities')

    campaign = out/'tournament'
    run([sys.executable, HERE/'tournament.py', 'prepare', '--source-root', prepared/'source',
         '--out', campaign, '--candidates', out/'candidates.json', '--cells', '13a0',
         '--holdout-cells', '17a1', '--profile', 'pilot', '--seed', '2026092530',
         '--timeout', '180', '--max-processes', '200', '--selection-width', '2',
         '--exploration-slots', '1'])
    for stage in ('aa', 'smoke'):
        run([sys.executable, campaign/'evaluator/tournament.py', 'run', '--round', campaign, '--stage', stage])
        result = read(campaign/'summaries'/f'{stage}.json')
        require(result['runs'] == result['verified_runs'], 'profiled driver retained failures')
    run([sys.executable, campaign/'evaluator/tournament.py', 'verify', '--round', campaign])
    require(not (campaign/'decision.json').exists(), 'integration control must not select a winner')
    qualification = out/'qualification-control'
    run([sys.executable, HERE/'tournament.py', 'prepare', '--source-root', prepared/'source',
         '--out', qualification, '--candidates', out/'candidates.json', '--cells', '13a0',
         '--holdout-cells', '17a1', '--profile', 'pilot', '--qualification', '--seed', '2026092540',
         '--timeout', '180', '--max-processes', '100', '--selection-width', '2', '--exploration-slots', '1'])
    run([sys.executable, qualification/'evaluator/tournament.py', 'run', '--round', qualification])
    run([sys.executable, qualification/'evaluator/tournament.py', 'verify', '--round', qualification])
    selected=read(qualification/'qualification.json')
    require(selected['status']=='DEVELOPMENT_REFERENCES_SELECTED' and len(selected['table'])==5,
            'qualification control did not retain both IC arms and three rho widths')
    require(set(read(qualification/'fixtures.json'))=={'aa','smoke','development'}
            and not (qualification/'decision.json').exists(), 'qualification accessed confirmation or promotion')
    summary = dict(qualification_control_pairs=66,
        scope='driver wiring controls; not a reference qualification or improvement round',
        native_verified=17, native_expected_failures=1, profile_native_pairs_verified=15, promotion_eligible=False,
        executed_stages=['aa', 'smoke'], prepared_but_unexecuted=['development', 'selection', 'confirmation', 'replay'])
    write(out/'summary.json', summary, exclusive=True)
    print(json.dumps(summary))


if __name__ == '__main__':
    main()
