"""Run the registered mixed-adapter control through the existing tournament."""
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
    parser.add_argument('--prepared', type=Path, required=True)
    parser.add_argument('--generic-build', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    config = dict(solver='pair_table', linear_algebra='tiny_gauss', summands=3,
                  batch_trials=1, max_trials=65536)
    registry = [dict(id='incumbent', config=config)]
    for alias, la, limit in (('generic_dense','dense',65536), ('generic_sparse','sparse',65536),
                             ('generic_budget','dense',1)):
        registry.append(dict(id=alias, adapter='generic-v1', generic_build=str(args.generic_build.resolve()),
            config=dict(config, linear_algebra=la, max_trials=limit)))
    write(out/'candidates.json', registry, exclusive=True)
    commands = []

    def run(command):
        number = len(commands)
        commands.append(list(map(str, command)))
        write(out/'commands.json', commands)
        with (out/f'command-{number}.log').open('x') as log:
            subprocess.run(commands[-1], stdout=log, stderr=subprocess.STDOUT, check=True)

    campaign = out/'tournament'
    run([sys.executable, HERE/'tournament.py', 'prepare', '--qualification',
         '--qualification-protocol', HERE/'goal_20260924/generic-qualification-adapter/PROTOCOL.md',
         '--source-root', args.prepared.resolve()/'source', '--out', campaign,
         '--candidates', out/'candidates.json', '--cells', '13a0', '--holdout-cells', '17a1',
         '--profile', 'pilot', '--seed', '2026092662', '--timeout', '180', '--max-processes', '120',
         '--selection-width', '2', '--exploration-slots', '1', '--qualification-widths', '1', '4',
         '--comparison-kind', 'factor-base-policy'])
    evaluator = campaign/'evaluator/tournament.py'
    run([sys.executable, evaluator, 'run', '--round', campaign])
    run([sys.executable, evaluator, 'verify', '--round', campaign])
    fixtures, report = read(campaign/'fixtures.json'), read(campaign/'qualification.json')
    require(set(fixtures) == {'aa','smoke','development'} and not (campaign/'decision.json').exists(),
            'mixed control entered confirmation or promotion')
    require(len(report['table']) == 8 and not report['promotion_eligible'], 'changed mixed control panel')
    require(report['eligible_for_improvement'] is False and report['observer_qualification'] is None,
            'wiring control acquired research or observer qualification')
    rows = [read(p) for p in sorted((campaign/'runs').glob('**/receipt.json'))]
    complete = [r for r in rows if r['status'] == 'VERIFIED']
    failed = [r for r in rows if r['status'] != 'VERIFIED']
    require(len(rows) == 102 and len(complete) == 90 and len(failed) == 12, 'unexpected mixed outcomes')
    require(all(r['arm']=='generic_budget' and r['status']=='INVALID_OR_INCOMPLETE'
                and r['measurement']['total_operations'] is None
                and r['measurement']['native_timing'] is None
                and r['measurement']['native_process_status']=='NOT_RUN' for r in failed),
            'intentional incomplete control acquired a result')
    require(all(r['qualified'] == (r['alias'] != 'generic_budget') for r in report['table']),
            'mixed eligibility does not match retained outcomes')
    keys = []
    for r in rows:
        run = r['measurement']
        keys.append(run['run_id'])
        if run.get('adapter') == 'generic-v1':
            keys.append(run['profile_execution']['run_id'])
    require(len(keys) == len(set(keys)), 'mixed execution IDs collide')
    for stage in ('smoke', 'development'):
        for case in fixtures[stage]:
            selected = [r for r in rows if r['stage']==stage and r['case']==case['id']]
            require(len({r['measurement']['workload_id'] for r in selected}) == 1,
                    'mixed adapters silently changed workload identity')
    summary = dict(status='PASS', trial_slots=102, verified_native_profile_pairs=90,
        retained_expected_incomplete_profiles=12, distinct_execution_keys=len(keys),
        frozen_replay=True, promotion_eligible=False, performance_qualified=False,
        scope='mixed-adapter wiring control; no research reference or observer qualification')
    write(out/'summary.json', summary, exclusive=True)
    print(json.dumps(summary))


if __name__ == '__main__':
    main()
