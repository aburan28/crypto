#!/usr/bin/env python3
"""Execute the registered F4/F5 subspace smoke qualification once.

Seed 2026093001 is reserved. Lost v1/v2 exposures are excluded. SAT is out of
scope. Without --out this module only locks the registration. Seeds 2026092901
and 2026092902 are not retried.
"""
import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

import campaign_rules_v2 as rules
from generic_build import verify_build_record
from generic_solver_feasibility import assess
from oracle import require
from run_improvement_v3 import EXPOSED, PANEL as ROUND3_PANEL, PANEL_SHA256 as ROUND3_PANEL_SHA256
from tournament import digest, read, write

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
REGISTRATION = HERE / 'goal_20260924/generic-backend-qualification-v3'
PANEL = REGISTRATION / 'panel.json'
PANEL_SHA256 = 'df92d5507785446a2a5b333bd7a04776781de4f54c63c5ea3ffa921bc99ad2e4'
V2 = HERE / 'goal_20260924/generic-backend-qualification-v2'
LOST_V1 = V2 / 'lost-campaign-exposures.json'
LOST_V1_SHA256 = 'a728677b199eac02800d8338204d5306f391ec5da757c910bec1e51955fe7b41'
LOST_V2 = V2 / 'lost-v2-campaign-exposures.json'
LOST_V2_SHA256 = '0cc792cceb8c7190a533e6f4e665e8486911153ad57609ac8d1243af282b8e54'
GENERIC_ROOT = Path(os.environ.get('IC_GENERIC_SOURCE_ROOT', ROOT)).resolve()
HISTORY = HERE / 'goal_20260924/improvement/target-history.json'
ARCHIVES = ('ic-improvement-round1-20260925', 'ic-improvement-round2-20260928',
            'ic-improvement-round3-20260928')
IC_ALIASES = ('incumbent', 'prepared_both', 'generic_pair_subspace_dense',
              'generic_f4_subspace_dense', 'generic_f5_subspace_dense')
GENERIC_SOURCE_OBJECTS = {
    'src': 'f87695b9adbc21b2e1af15a3923675dcae1c3dd7',
    'examples/ic_tournament_worker.rs': '745247f3d89d46cdffb6d1f28778971caae7aa21',
    'Cargo.toml': '4177c1daa3b7f779abcf5fdeecb3b284d93b0f19',
    'research/ic_candidate_tournament_20260915/ci/Cargo.lock':
        'a3181ddde6d0a7460f0d9a1c6e87e3bead801980',
}
# This commit adds the single campaign path. Do not redispatch after one run.
DISPATCH_AUTHORIZED = True


def check_generic_source(root=GENERIC_ROOT):
    """Require the registered worker checkout before any fixture is generated."""
    require(subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=root, text=True).strip()
            == '765c3c5f19032bd852163805f257c56babef2040',
            'generic source is not the registered frozen checkout')
    for path, expected in GENERIC_SOURCE_OBJECTS.items():
        actual = subprocess.check_output(['git', 'rev-parse', 'HEAD:'+path],
                                         cwd=root, text=True).strip()
        require(actual == expected, 'registered generic worker source changed: '+path)


def registry(panel, generic_build):
    """Map the public registration to executable arms without alias duplication."""
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
        require(config.get('factor_base') == dict(kind='standard_subspace', dimension=6),
                'changed registered subspace factor-base policy')
        require(config['summands'] == 3, 'changed registered generic summand count')
        rows.append(dict(id=row['id'], adapter='generic-v1',
                         generic_build=str(generic_build), config=config))
    require(len(rows) == 4 and len({row['id'] for row in rows}) == 4,
            'registered candidate panel is incomplete')
    return rows


def registration_check(panel):
    """Reject a panel that is not the frozen planning registration."""
    require(digest(PANEL) == PANEL_SHA256, 'changed premeasurement v3 panel bytes')
    require(panel['schema_version'] == 1 and panel['status'] == 'REGISTERED_PLANNING'
            and panel['seed'] == 2026093001
            and panel['comparison_kind'] == 'factor-base-policy'
            and panel['stages'] == ['aa', 'smoke']
            and panel['cells'] == ['n17a1', 'n19a0', 'n23a0', 'n23a1', 'n31a0']
            and panel['holdout_cells_excluded'] == ['n29a1']
            and panel['repetitions'] == 1 and panel['timeout_seconds'] == 300
            and panel['memory_bytes'] == 8 * 1024**3
            and panel['measure_timeout_minutes'] == 180
            and panel['pack_timeout_minutes'] == 35
            and panel['job_timeout_minutes'] == 240
            and panel['scheduled_pair_bound'] == 45 and panel['pair_cap'] == 60
            and panel['sat_in_scope'] is False and panel['promotion_eligible'] is False
            and panel['worker_commit'] == '765c3c5f19032bd852163805f257c56babef2040',
            'changed registered v3 schedule or resources')
    ids = [row['id'] for row in panel['candidates']]
    require(ids == list(IC_ALIASES), 'changed registered v3 candidate ids')
    for row in panel['candidates'][2:]:
        base = row['config']['factor_base']
        require(row['adapter'] == 'generic-v1' and row['source'] == 'controlled-generic'
                and base == {'kind': 'standard_subspace', 'dimension': 6},
                'changed subspace arm '+row['id'])
    for row in panel['candidates'][3:]:
        config = row['config']
        require(config['summands'] == 3 and config['groebner_degree'] == 3
                and config['node_budget'] == 4096 and config['max_trials'] == 256
                and config['batch_trials'] == 8 and config['linear_algebra'] == 'dense',
                'changed algebraic budget '+row['id'])
    require(panel['rho_arms'] == [
        {'id': 'rho', 'role': 'cold_rho_reference', 'config': {'rho_parallel_walks': 8}},
        {'id': 'rho_online', 'role': 'online_rho_reference',
         'config': {'rho_parallel_walks': 16, 'linear_algebra': 'dense'}}],
        'changed registered v3 rho roles')
    require(digest(LOST_V1) == LOST_V1_SHA256, 'changed frozen first-run exposures')
    layout = assess(panel)
    require(layout['status'] == 'PASS_STATIC_LAYOUT_ONLY'
            and layout['impossible_algebraic_cells'] == 0,
            'v3 algebraic layout is outside the encoder cap')
    return layout


def exposure_block():
    """Why seed 2026093001 must not sample points on this checkout."""
    if not LOST_V2.is_file():
        return 'v2 exposure census is not on this checkout'
    if digest(LOST_V2) != LOST_V2_SHA256:
        return 'v2 exposure census hash does not match the sealed registration'
    lost = read(LOST_V2)
    if lost.get('seed') != 2026092902:
        return 'v2 exposure census is not seed 2026092902'
    return None


def status_report(panel):
    census_block = exposure_block()
    dispatch_block = census_block
    if census_block is None and not DISPATCH_AUTHORIZED:
        dispatch_block = 'dispatch path is not in this checker'
    return {
        'seed': panel['seed'],
        'panel_sha256': PANEL_SHA256,
        'static_layout': 'PASS_STATIC_LAYOUT_ONLY',
        'v2_exposure_sha256': LOST_V2_SHA256,
        'v2_exposure_present': census_block is None,
        'dispatch_block': dispatch_block,
        'dispatch_authorized': DISPATCH_AUTHORIZED,
        'measurement': 'not_run',
        'qualification_schedule': 'smoke',
    }


def panel_result(campaign):
    contract, fixtures = read(campaign/'contract.json'), read(campaign/'fixtures.json')
    require(contract['purpose'] == 'reference-qualification'
            and contract['qualification_reference_schema'] == 1
            and contract.get('qualification_schedule') == 'smoke'
            and contract['reference_qualification'] == rules.expected_binding()
            and contract['seed'] == 2026093001 and contract['repetitions'] == 1
            and contract['cells'] == ['n17a1', 'n19a0', 'n23a0', 'n23a1', 'n31a0']
            and contract['stages'] == ['aa', 'smoke']
            and [a['id'] for a in contract['reference_arms']] == list(rules.REFERENCE_ROLES),
            'changed source-bound smoke qualification contract')
    require({stage: len(cases) for stage, cases in fixtures.items()}
            == dict(aa=5, smoke=5)
            and 'development' not in fixtures
            and not (campaign/'decision.json').exists()
            and not (campaign/'qualification.json').exists(),
            'smoke schedule generated development or a promotion decision')
    rows = [read(path) for path in sorted(campaign.glob('runs/**/receipt.json'))]
    require(len(rows) == 45, 'incomplete schedule; preserve raw work without claiming qualification')
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
    return dict(status='AUDITED_WITH_FAILURES' if any(row['status'] != 'VERIFIED' for row in rows)
                else 'AUDITED_COMPLETE', trial_slots=45, repetitions=1,
                verified_native_profile_pairs=sum(row['status'] == 'VERIFIED' for row in rows),
                retained_failures=sum(row['status'] != 'VERIFIED' for row in rows),
                exposed_points=10, distinct_record_ids=len(keys),
                prior_censored_exposures_sha256=LOST_V1_SHA256,
                prior_v2_censored_exposures_sha256=LOST_V2_SHA256,
                registration_panel_sha256=PANEL_SHA256,
                qualification_schedule='smoke', sat_in_scope=False,
                promotion_eligible=False, improvement_rounds_used=0,
                scope='five-cell F4/F5 subspace smoke qualification; SAT out of scope; '
                      'no held-out confirmation or global claim')


def run_campaign(out):
    out = out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    panel = read(PANEL)
    registration_check(panel)
    require(DISPATCH_AUTHORIZED, 'v3 dispatch is not authorized in this checker')
    require(exposure_block() is None, exposure_block() or 'v2 exposure census missing')
    check_generic_source()
    require(digest(HISTORY) == rules.HISTORY_SHA256, 'changed original point census')
    require(digest(ROUND3_PANEL) == ROUND3_PANEL_SHA256,
            'changed sealed third-round supplemental exposure registry')
    expected_exposures = read(ROUND3_PANEL)['exposed_fixture_sha256']
    for path in EXPOSED:
        relative = str(path.relative_to(ROOT))
        require(path.is_file() and digest(path) == expected_exposures[relative],
                'changed exposed fixture corpus: '+relative)
    write(out/'registered-panel.json', panel, exclusive=True)
    shutil.copy2(LOST_V1, out/'prior-censored-exposures.json')
    shutil.copy2(LOST_V2, out/'prior-v2-censored-exposures.json')
    shutil.copy2(REGISTRATION/'PROTOCOL.md', out/'PROTOCOL.md')
    shutil.copy2(Path(__file__), out/'registered-runner.py')
    archive_manifest = read(HERE/'evidence/manifest.json')
    needed = {'ic-generic-reference-qualification-20260926', *ARCHIVES}
    inventory = {entry['file'].removesuffix('.tar.zst'): entry
                 for entry in archive_manifest['archives']
                 if entry['file'].removesuffix('.tar.zst') in needed}
    require(set(inventory) == needed, 'missing sealed archive in evidence manifest')
    write(out/'restored-archives.json', inventory, exclusive=True)
    commands = []

    def run(command):
        number = len(commands)
        commands.append(list(map(str, command)))
        write(out/'commands.json', commands)
        print(json.dumps({'command_number': number, 'argv': commands[-1]}), flush=True)
        with (out/f'command-{number}.log').open('x') as log:
            process = subprocess.Popen(commands[-1], stdout=subprocess.PIPE,
                                       stderr=subprocess.STDOUT, text=True, bufsize=1)
            for line in process.stdout:
                log.write(line)
                log.flush()
                print(line.rstrip(), flush=True)
            completed = process.wait()
        write(out/f'command-{number}.status.json',
              dict(exit_code=completed), exclusive=True)
        require(completed == 0,
                f'command {number} failed; preserved log and partial evidence')

    run([sys.executable, HERE/'target_history.py', '--history', HISTORY,
         '--repository', ROOT])
    references_root = out.parent/(out.name+'-reference-evidence')
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
             '--out', out.parent/(out.name+'-'+archive)])
    previous = [out.parent/(out.name+'-'+archive)/archive.replace('-20260925','').replace('-20260928','')
                /'round/tournament' for archive in ARCHIVES]
    for path in previous:
        require((path/'contract.json').is_file() and (path/'decision.json').is_file(),
                'sealed prior round not restored')
    generic_build = out/'new-generic-build'
    run([sys.executable, HERE/'generic_build.py', '--root', GENERIC_ROOT, '--out', generic_build])
    build, source = read(generic_build/'build-record.json'), read(generic_build/'source-manifest.json')
    verify_build_record(build, source)
    require(digest(generic_build/'worker') == build['worker_sha256'], 'changed generic executable')
    for source_root in (prepared['pairinv']/'source', prepared['both']/'source'):
        run(['cargo', 'fetch', '--locked', '--manifest-path', source_root/'Cargo.toml'])
    candidates = out/'candidates.json'
    write(candidates, registry(panel, generic_build), exclusive=True)
    campaign = out/'tournament'
    command = [sys.executable, HERE/'tournament.py', 'prepare', '--qualification',
               '--qualification-schedule', 'smoke',
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
    command.extend(['--exposed-fixtures', LOST_V1, '--exposed-fixtures', LOST_V2])
    command.extend(['--out', campaign, '--candidates', candidates,
                    '--cells', '17a1,19a0,23a0,23a1,31a0', '--holdout-cells', '29a1',
                    '--profile', 'pilot', '--seed', str(panel['seed']),
                    '--repetitions', '1', '--timeout', '300',
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
    shutil.copy2(HERE/'generic_backend_gate_v3.py', out/'generic_backend_gate_v3.py')
    run([sys.executable, HERE/'generic_backend_gate_v3.py',
         '--bundle', out, '--out', out/'family-gate.json'])
    print(json.dumps(result), flush=True)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path,
                        help='Measure once into this empty directory. Omit for the lock check.')
    args = parser.parse_args()
    panel = read(PANEL)
    registration_check(panel)
    report = status_report(panel)
    require(LOST_V2.is_file() and report['v2_exposure_present'] is True,
            report['dispatch_block'] or 'v2 exposure census is not on this checkout')
    require(report['dispatch_authorized'] is True, 'v3 dispatch flag is not enabled')
    require(report['dispatch_block'] is None, 'v3 checker still blocks sampling')
    if args.out is None:
        print(json.dumps(report, sort_keys=True), flush=True)
        return 0
    run_campaign(args.out)
    return 0


if __name__ == '__main__':
    sys.exit(main())
