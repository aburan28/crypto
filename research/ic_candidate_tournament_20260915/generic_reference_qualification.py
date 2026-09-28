"""Execute the registered generic/reference development panel, never a final round."""
import argparse
import copy
import json
from pathlib import Path
import shutil
import subprocess
import sys

from campaign_rules import IC_SOURCE, COLD_RHO_SOURCE, CELLS, CONFIG
from generic_build import verify_build_record
from identity import sha256
from oracle import require
from tournament import digest, read, write

HERE = Path(__file__).resolve().parent
PROTOCOL = HERE / 'goal_20260924/generic-reference-qualification/PROTOCOL.md'
WIDTHS = (1, 2, 4, 8, 16, 32)
SEED = 2026092663
SLOTS = 1350


def prepared_sources(artifacts):
    prepared = {}
    for reference, expected in (('pairinv', IC_SOURCE), ('both', COLD_RHO_SOURCE)):
        matches = list(artifacts.glob(f'ic-producer-{reference}-*/ic-producer'))
        require(len(matches) == 1, 'expected one admitted artifact for ' + reference)
        directory = matches[0]
        receipt, manifest = read(directory / 'preparation.json'), read(directory / 'source-manifest.json')
        require(receipt['reference'] == reference and receipt['instrumented'] is True
                and receipt['source_manifest_sha256'] == sha256(manifest) == expected,
                'changed qualified source binding: ' + reference)
        actual = {p.relative_to(directory / 'source').as_posix(): digest(p)
                  for p in (directory / 'source').rglob('*') if p.is_file()}
        require(actual == manifest, 'qualified source bytes differ: ' + reference)
        prepared[reference] = directory
    return prepared


def registry(prepared, generic_build):
    return [dict(id='incumbent', config=copy.deepcopy(CONFIG)),
        dict(id='prepared_both', config=copy.deepcopy(CONFIG), source_root=str(prepared['both'] / 'source'))] + [
        dict(id='generic_' + mode, adapter='generic-v1', generic_build=str(generic_build),
             config=dict(CONFIG, linear_algebra=mode)) for mode in ('dense', 'sparse')]


def panel_result(campaign):
    contract, fixtures = read(campaign / 'contract.json'), read(campaign / 'fixtures.json')
    require(contract['purpose'] == 'reference-qualification' and contract['seed'] == SEED
            and contract['cells'] == CELLS and contract['repetitions'] == 3,
            'changed registered development panel')
    require(set(fixtures) == {'aa', 'smoke', 'development'}
            and [len(fixtures[s]) for s in ('aa', 'smoke', 'development')] == [5, 5, 15]
            and not (campaign / 'decision.json').exists(), 'qualification entered a final stage')
    arms = read(campaign / 'candidates.json')
    require([a['id'] for a in arms] == ['incumbent', 'prepared_both', 'generic_dense', 'generic_sparse']
            and arms[0]['source_manifest_sha256'] == IC_SOURCE
            and arms[1]['source_manifest_sha256'] == COLD_RHO_SOURCE,
            'qualified incumbents were weakened or substituted')
    references = contract['reference_arms']
    require(len(references) == 18, 'incomplete rho source/width panel')
    for source in {a['source_manifest_sha256'] for a in arms}:
        require(sorted(a['config']['rho_parallel_walks'] for a in references
                       if a['source_manifest_sha256'] == source) == list(WIDTHS),
                'changed rho width screen')
    report = read(campaign / 'qualification.json')
    require(len(report['table']) == 22 and report['promotion_eligible'] is False
            and report['eligible_for_improvement'] is False
            and report['observer_qualification'] is None, 'development acquired unearned qualification')
    rows = [read(p) for p in sorted(campaign.glob('runs/**/receipt.json'))]
    require(len(rows) == SLOTS, 'partial schedule; retained data cannot masquerade as a complete panel')
    keys = []
    for row in rows:
        measurement = row['measurement']
        keys.append(measurement['run_id'])
        if measurement.get('adapter') == 'generic-v1':
            keys.append(measurement['profile_execution']['run_id'])
        if row['status'] != 'VERIFIED':
            require(measurement['total_operations'] is None and measurement['native_timing'] is None,
                    'incomplete execution acquired a measured result')
    require(len(keys) == len(set(keys)), 'execution identifiers collide')
    for stage, cases in fixtures.items():
        for case in cases:
            paired = [r for r in rows if r['stage'] == stage and r['case'] == case['id']]
            require(len({r['measurement']['workload_id'] for r in paired}) == 1,
                    'references and candidates silently changed workload identity')
    verified = sum(r['status'] == 'VERIFIED' for r in rows)
    return dict(status='COMPLETE_PANEL' if verified == SLOTS else 'PANEL_WITH_RETAINED_FAILURES',
        trial_slots=SLOTS, verified_native_profile_pairs=verified, retained_failures=SLOTS-verified,
        qualification_sha256=digest(campaign / 'qualification.json'),
        fixture_sha256=digest(campaign / 'fixtures.json'), exposed_points=25,
        distinct_record_ids=len(keys), promotion_eligible=False, improvement_rounds_used=0,
        scope='fully charged instrumented development comparison; observer study and new reference binding still required')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--artifacts', type=Path, required=True)
    parser.add_argument('--generic-build', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    prepared = prepared_sources(args.artifacts.resolve())
    build_dir = args.generic_build.resolve()
    build, source = read(build_dir / 'build-record.json'), read(build_dir / 'source-manifest.json')
    verify_build_record(build, source)
    require(digest(build_dir / 'worker') == build['worker_sha256'], 'changed controlled generic worker')
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    shutil.copy2(PROTOCOL, out / 'PROTOCOL.md')
    shutil.copy2(Path(__file__), out / 'registered-runner.py')
    commands = []

    def run(command):
        number = len(commands)
        commands.append(list(map(str, command)))
        write(out / 'commands.json', commands)
        with (out / f'command-{number}.log').open('x') as log:
            subprocess.run(commands[-1], stdout=log, stderr=subprocess.STDOUT, check=True)

    for reference in ('pairinv', 'both'):
        run(['cargo', 'fetch', '--locked', '--manifest-path', prepared[reference] / 'source/Cargo.toml'])
    write(out / 'candidates.json', registry(prepared, build_dir), exclusive=True)
    campaign = out / 'tournament'
    run([sys.executable, HERE / 'tournament.py', 'prepare', '--qualification',
        '--qualification-protocol', PROTOCOL, '--source-root', prepared['pairinv'] / 'source',
        '--out', campaign, '--candidates', out / 'candidates.json',
        '--cells', ','.join(c.removeprefix('n') for c in CELLS), '--holdout-cells', '29a1',
        '--profile', 'pilot', '--seed', str(SEED), '--timeout', '180', '--max-processes', '1400',
        '--qualification-widths', *map(str, WIDTHS), '--selection-width', '4', '--exploration-slots', '1',
        '--comparison-kind', 'factor-base-policy'])
    observer = out / 'observer'
    run([sys.executable, HERE / 'generic_observer.py', 'prepare', '--round', campaign, '--out', observer])
    run([sys.executable, campaign / 'evaluator/tournament.py', 'run', '--round', campaign])
    run([sys.executable, campaign / 'evaluator/tournament.py', 'verify', '--round', campaign])
    result = panel_result(campaign)
    write(out / 'comparison-summary.json', result, exclusive=True)
    run([sys.executable, observer / 'evaluator/generic_observer.py', 'run', '--out', observer])
    run([sys.executable, observer / 'evaluator/generic_observer.py', 'verify', '--out', observer])
    observed = read(observer / 'summary.json')
    require(observed['pairs'] == 360 and observed['failed_pairs'] == 0, 'incomplete observer study')
    result.update(observer_summary_sha256=digest(observer / 'summary.json'),
                  observer_pairs=360, observer_complete_pairs=observed['complete_pairs'],
                  observer_incomplete_pairs=observed['incomplete_pairs'],
                  accepted_reference_binding_changed=False)
    write(out / 'summary.json', result, exclusive=True)
    print(json.dumps(result))


if __name__ == '__main__':
    main()
