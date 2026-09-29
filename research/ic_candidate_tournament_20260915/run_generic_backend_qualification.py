#!/usr/bin/env python3
"""Execute the registered F4/F5/SAT complete-solve qualification once.

This is a development comparison, not a fourth improvement attempt. Restore
and verify all sealed prior points before generating any new public point.
"""
import argparse
import json
from pathlib import Path
import shutil
import subprocess
import sys

import campaign_rules_v2 as rules
from generic_build import verify_build_record
from oracle import require
from run_improvement_v3 import EXPOSED, PANEL as ROUND3_PANEL
from tournament import digest, read, write

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
REGISTRATION = HERE / 'goal_20260924/generic-backend-qualification'
PANEL = REGISTRATION / 'panel.json'
PANEL_SHA256 = '83c640a03b4239b918851f6f1b8450e2fe27dc710fb99f306e3481a99d8875cf'
HISTORY = HERE / 'goal_20260924/improvement/target-history.json'
ARCHIVES = ('ic-improvement-round1-20260925', 'ic-improvement-round2-20260928',
            'ic-improvement-round3-20260928')
IC_ALIASES = ('incumbent', 'prepared_both', 'generic_pair_dense',
              'generic_pair_sparse', 'generic_f4_dense', 'generic_f4_sparse',
              'generic_f5_dense', 'generic_sat_xor_dense',
              'generic_sat_cnf_dense', 'generic_inherited_f4_dense')
GENERIC_SOURCE_OBJECTS = {
    'src': '6caf5dd2de704de77fd3cae5affab8b0ec3b6911',
    'examples/ic_tournament_worker.rs': '745247f3d89d46cdffb6d1f28778971caae7aa21',
    'Cargo.toml': '4177c1daa3b7f779abcf5fdeecb3b284d93b0f19',
    'research/ic_candidate_tournament_20260915/ci/Cargo.lock':
        'a3181ddde6d0a7460f0d9a1c6e87e3bead801980',
}


def check_generic_source(root=ROOT):
    """Reject main-branch worker drift before any fixture is generated."""
    for path, expected in GENERIC_SOURCE_OBJECTS.items():
        actual = subprocess.check_output(['git', 'rev-parse', 'HEAD:'+path],
                                         cwd=root, text=True).strip()
        require(actual == expected, 'registered generic worker source changed: '+path)


def registry(panel, generic_build):
    """Map the public registration to executable arms without alias duplication.

    `prepared_both` is the accepted `ic_online` reference, not a second
    candidate with identical method identity. Factor-base construction belongs
    to each job; the panel's descriptive 6n recipe is checked here.
    """
    require([row['id'] for row in panel['candidates']] == list(IC_ALIASES),
            'changed registered IC aliases')
    require(panel['candidates'][0]['config'] == rules.CONFIG
            and panel['candidates'][1]['config'] == rules.CONFIG,
            'changed prepared IC references')
    rows = [dict(id='incumbent', config=rules.CONFIG)]
    for row in panel['candidates'][2:]:
        require(row['adapter'] == 'generic-v1' and row['source'] == 'controlled-generic',
                'unknown generic adapter or source')
        config = dict(row['config'])
        require(config.pop('factor_base') == dict(recipe='subgroup_orbits', seed=43,
                                                  requested_points='6n'),
                'changed registered factor-base policy')
        require(config['summands'] == 3 and config['batch_trials'] == 1
                and config['max_trials'] == 65536,
                'changed registered generic query policy')
        rows.append(dict(id=row['id'], adapter='generic-v1',
                         generic_build=str(generic_build), config=config))
    require(len(rows) == 9 and len({row['id'] for row in rows}) == 9,
            'registered candidate panel is incomplete')
    return rows


def registration_check(panel):
    require(digest(PANEL) == PANEL_SHA256, 'changed premeasurement panel bytes')
    require(panel['schema_version'] == 1 and panel['seed'] == 2026092901
            and panel['stages'] == ['aa', 'smoke', 'development']
            and panel['cells'] == ['n17a1', 'n19a0', 'n23a0', 'n23a1', 'n31a0']
            and panel['holdout_cells_excluded'] == ['n29a1']
            and panel['repetitions'] == 3 and panel['timeout_seconds'] == 300
            and panel['memory_bytes'] == 8 * 1024**3
            and panel['scheduled_pair_bound'] == 750 and panel['pair_cap'] == 900,
            'changed registered schedule or resources')
    require(panel['rho_arms'] == [
        dict(id='rho', role='cold_rho_reference', config=dict(rho_parallel_walks=8)),
        dict(id='rho_online', role='online_rho_reference',
             config=dict(rho_parallel_walks=16, linear_algebra='dense'))],
        'changed registered matched rho roles')


def panel_result(campaign):
    contract, fixtures = read(campaign/'contract.json'), read(campaign/'fixtures.json')
    require(contract['purpose'] == 'reference-qualification'
            and contract['qualification_reference_schema'] == 1
            and contract['reference_qualification'] == rules.expected_binding()
            and contract['seed'] == 2026092901 and contract['repetitions'] == 3
            and contract['cells'] == ['n17a1','n19a0','n23a0','n23a1','n31a0']
            and [a['id'] for a in contract['reference_arms']] == list(rules.REFERENCE_ROLES),
            'changed source-bound qualification contract')
    require({stage:len(cases) for stage,cases in fixtures.items()}
            == dict(aa=5, smoke=5, development=15)
            and not (campaign/'decision.json').exists(),
            'qualification generated a final-stage workload')
    rows = [read(path) for path in sorted(campaign.glob('runs/**/receipt.json'))]
    require(len(rows) == 750, 'incomplete schedule; preserve raw work without claiming qualification')
    keys = []
    for row in rows:
        measured = row['measurement']
        keys.append(measured['run_id'])
        if measured.get('adapter') == 'generic-v1':
            keys.append(measured['profile_execution']['run_id'])
        if row['status'] != 'VERIFIED':
            require(measured['total_operations'] is None
                    and measured['native_timing'] is None,
                    'failed job acquired complete competitive costs')
    require(len(keys) == len(set(keys)), 'colliding execution identities')
    for stage, cases in fixtures.items():
        for case in cases:
            paired = [row for row in rows if row['stage'] == stage and row['case'] == case['id']]
            require(len({row['measurement']['workload_id'] for row in paired}) == 1,
                    'candidate and references used different public workloads')
    qualified = read(campaign/'qualification.json')
    return dict(status='AUDITED_WITH_FAILURES' if any(row['status'] != 'VERIFIED' for row in rows)
                else 'AUDITED_COMPLETE', trial_slots=750,
                verified_native_profile_pairs=sum(row['status'] == 'VERIFIED' for row in rows),
                retained_failures=sum(row['status'] != 'VERIFIED' for row in rows),
                exposed_points=25, distinct_record_ids=len(keys),
                selected_development_aliases=dict(ic_cold=qualified['selected_ic_cold'],
                    ic_online=qualified['selected_ic_online'],
                    rho_cold=qualified['selected_rho_cold'],
                    rho_online=qualified['selected_rho_online']),
                promotion_eligible=False, improvement_rounds_used=0,
                scope='five-cell fully charged development comparison; no held-out or global claim')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    panel = read(PANEL)
    registration_check(panel)
    check_generic_source()
    require(digest(HISTORY) == rules.HISTORY_SHA256, 'changed original point census')
    expected_exposures = read(ROUND3_PANEL)['exposed_fixture_sha256']
    for path in EXPOSED:
        relative = str(path.relative_to(ROOT))
        require(path.is_file() and digest(path) == expected_exposures[relative],
                'changed exposed fixture corpus: '+relative)
    write(out/'registered-panel.json', panel, exclusive=True)
    shutil.copy2(REGISTRATION/'PROTOCOL.md', out/'PROTOCOL.md')
    shutil.copy2(Path(__file__), out/'registered-runner.py')
    commands = []

    def run(command):
        number = len(commands)
        commands.append(list(map(str, command)))
        write(out/'commands.json', commands)
        with (out/f'command-{number}.log').open('x') as log:
            completed = subprocess.run(commands[-1], stdout=log,
                                       stderr=subprocess.STDOUT, check=False)
        write(out/f'command-{number}.status.json',
              dict(exit_code=completed.returncode), exclusive=True)
        require(completed.returncode == 0,
                f'command {number} failed; preserved log and partial evidence')

    run([sys.executable, HERE/'target_history.py', '--history', HISTORY,
         '--repository', ROOT])
    references_root = out/'reference-evidence'
    run([sys.executable, HERE/'evidence/restore.py', '--archive',
         'ic-generic-reference-qualification-20260926', '--out', references_root])
    bundle = references_root/'ic-generic-reference-qualification'
    registry_path = out/'reference-registry.json'
    run([sys.executable, HERE/'reference_registry_v2.py', '--bundle', bundle,
         '--out', registry_path])
    declarations = read(registry_path)
    rules.validate_reference_declarations(declarations)
    from generic_reference_qualification import prepared_sources
    prepared = prepared_sources(bundle/'controls')
    require(read(prepared['pairinv']/'preparation.json')['source_manifest_sha256']
            == rules.IC_SOURCE, 'weakened accepted incumbent')
    for archive in ARCHIVES:
        run([sys.executable, HERE/'evidence/restore.py', '--archive', archive,
             '--out', out/archive])
    previous = [out/archive/archive.replace('-20260925','').replace('-20260928','')
                /'round/tournament' for archive in ARCHIVES]
    for path in previous:
        require((path/'contract.json').is_file() and (path/'decision.json').is_file(),
                'sealed prior round not restored')
    generic_build = out/'new-generic-build'
    run([sys.executable, HERE/'generic_build.py', '--root', ROOT, '--out', generic_build])
    build, source = read(generic_build/'build-record.json'), read(generic_build/'source-manifest.json')
    verify_build_record(build, source)
    require(digest(generic_build/'worker') == build['worker_sha256'], 'changed generic executable')
    for source_root in (prepared['pairinv']/'source', prepared['both']/'source'):
        run(['cargo', 'fetch', '--locked', '--manifest-path', source_root/'Cargo.toml'])
    candidates = out/'candidates.json'
    write(candidates, registry(panel, generic_build), exclusive=True)
    campaign = out/'tournament'
    command = [sys.executable, HERE/'tournament.py', 'prepare', '--qualification',
               '--qualification-protocol', REGISTRATION/'PROTOCOL.md',
               '--source-root', prepared['pairinv']/'source',
               '--reference-registry', registry_path,
               '--qualified-report', bundle/'tournament/qualification.json',
               '--qualified-observer', bundle/'observer/summary.json',
               '--target-history', HISTORY]
    for path in previous:
        command.extend(['--prior-round', path])
    for path in EXPOSED:
        command.extend(['--exposed-fixtures', path])
    command.extend(['--out', campaign, '--candidates', candidates,
                    '--cells', '17a1,19a0,23a0,23a1,31a0', '--holdout-cells', '29a1',
                    '--profile', 'pilot', '--seed', str(panel['seed']), '--timeout', '300',
                    '--max-processes', str(panel['pair_cap']),
                    '--selection-width', '6', '--exploration-slots', '1',
                    '--comparison-kind', 'factor-base-policy'])
    run(command)
    evaluator = campaign/'evaluator/tournament.py'
    run([sys.executable, evaluator, 'run', '--round', campaign])
    run([sys.executable, evaluator, 'verify', '--round', campaign])
    run([sys.executable, campaign/'evaluator/generic_backend_yield.py',
         '--round', campaign, '--out', out/'natural-yield.json'])
    result = panel_result(campaign)
    result['natural_yield_sha256'] = digest(out/'natural-yield.json')
    write(out/'summary.json', result, exclusive=True)
    print(json.dumps(result), flush=True)


if __name__ == '__main__':
    main()
