"""Replay the retained mixed-adapter control without executing a worker."""
import argparse
import collections
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tarfile
import tempfile

HERE = Path(__file__).resolve().parent


def require(condition, message):
    if not condition:
        raise ValueError(message)


def read(path):
    return json.loads(path.read_text())


def restore(destination):
    manifest = read(HERE / 'EVIDENCE.json')
    archive = HERE / manifest['archive']
    data = archive.read_bytes()
    require(len(data) == manifest['archive_bytes'] and
            hashlib.sha256(data).hexdigest() == manifest['archive_sha256'], 'changed archive')
    with tarfile.open(archive) as source:
        members = source.getmembers()
        expected = {manifest['archive_root'] + '/' + p for p in manifest['files']}
        require(len(members) == len(expected) and {m.name for m in members} == expected
                and all(m.isfile() for m in members), 'changed archive member set')
        source.extractall(destination, filter='data')
    root = destination / manifest['archive_root']
    for name, expected in manifest['files'].items():
        data = (root / name).read_bytes()
        require(len(data) == expected['bytes'] and hashlib.sha256(data).hexdigest() == expected['sha256'],
                'changed evidence member: ' + name)
    require(len(manifest['files']) == manifest['file_count'] and
            sum(v['bytes'] for v in manifest['files'].values()) == manifest['expanded_bytes'],
            'changed archive census')
    require((root / 'tournament/fixtures.json').read_bytes() == (HERE / 'fixtures.json').read_bytes(),
            'changed exposed fixture export')
    return root


def control_result(root):
    campaign = root / 'tournament'
    fixtures, qualification = read(campaign / 'fixtures.json'), read(campaign / 'qualification.json')
    require(set(fixtures) == {'aa', 'smoke', 'development'} and not (campaign / 'decision.json').exists(),
            'control entered a research final stage')
    require(qualification['eligible_for_improvement'] is False and
            qualification['observer_qualification'] is None and
            qualification['promotion_eligible'] is False, 'control acquired research qualification')
    rows = [read(p) for p in sorted(campaign.glob('runs/**/receipt.json'))]
    require(len(rows) == 102, 'changed trial census')
    counts = collections.Counter((r['arm'], r['status']) for r in rows)
    aliases = ['incumbent', 'aa_control', 'generic_dense', 'generic_sparse', 'generic_budget',
               'rho_incumbent_1', 'rho_incumbent_4', 'rho_generic_dense_1', 'rho_generic_dense_4']
    expected = {(a, 'INVALID_OR_INCOMPLETE' if a == 'generic_budget' else 'VERIFIED'):
                15 if a == 'incumbent' else 3 if a == 'aa_control' else 12 for a in aliases}
    require(counts == expected, 'changed complete/incomplete control outcomes')
    keys, unstarted_native_ids = [], []
    table = []
    for alias in aliases:
        selected = [r for r in rows if r['arm'] == alias]
        ids, inventory = set(), None
        for row in selected:
            run = row['measurement']
            keys.append(run['run_id'])
            ids.add(run.get('candidate_id', run.get('reference_id')))
            if run.get('adapter') == 'generic-v1':
                keys.append(run['profile_execution']['run_id'])
            if alias == 'generic_budget':
                require(run['total_operations'] is None and run['native_timing'] is None
                        and run['native_process_status'] == 'NOT_RUN' and run['certificate'] is None,
                        'incomplete profile acquired a verified native result')
                unstarted_native_ids.append(run['run_id'])
            admission = read(campaign / 'admissions' / row['stage'] / row['case'] / alias / 'admission.json')
            if 'candidate' in admission:
                base = admission['candidate']['record']['factor_base']['inventory']
                inventory = dict(usable_points=base['usable_point_count'], folded_columns=base['effective_columns'])
                require(inventory == dict(usable_points=182, folded_columns=7), 'changed actual base census')
        require(len(ids) == 1 and None not in ids, 'unstable control pipeline identity')
        table.append(dict(alias=alias, identity=next(iter(ids)), base=inventory,
            planned_slots=len(selected), verified_pairs=0 if alias == 'generic_budget' else len(selected),
            retained_incomplete_profiles=len(selected) if alias == 'generic_budget' else 0))
    require(len(keys) == 162 and len(set(keys)) == len(keys), 'execution identifiers collide')
    for stage, cases in fixtures.items():
        for case in cases:
            paired = [r for r in rows if r['stage'] == stage and r['case'] == case['id']]
            require(len({r['measurement']['workload_id'] for r in paired}) == 1,
                    'adapters changed the paired workload identity')
    require(all(r['qualified'] == (r['alias'] != 'generic_budget') for r in qualification['table'])
            and len(qualification['table']) == 8, 'failed control became eligible')
    raw = read(root / 'summary.json')
    require(raw == dict(status='PASS', trial_slots=102, verified_native_profile_pairs=90,
        retained_expected_incomplete_profiles=12, distinct_execution_keys=162, frozen_replay=True,
        promotion_eligible=False, performance_qualified=False,
        scope='mixed-adapter wiring control; no research reference or observer qualification'),
        'changed original control summary')
    return dict(schema_version=1, status='PASS', scope=raw['scope'],
        workflow_run=36280983990, source_commit='717879287fa389c1f2eb3817316446dc86664f97',
        trial_slots=102, verified_native_profile_pairs=90, retained_expected_incomplete_profiles=12,
        distinct_record_ids=162, unstarted_native_record_ids=len(unstarted_native_ids),
        recorded_id_note='Prepared records describe pairs; generic records distinguish native/profile, including reserved NOT_RUN native IDs.',
        exposed_points=sum(len(c['fixture']['targets']) for cases in fixtures.values() for c in cases),
        cell='n13a0', subgroup_order=2003, table=table,
        promotion_eligible=False, performance_qualified=False, improvement_rounds_used=0,
        comparative_claims=dict(online_speedup=None, cold_speedup=None, normalized_S=None),
        next_gate='Registered five-cell generic/reference and observer qualification; no improvement round consumed.')


def replay(root):
    campaign = root / 'tournament'
    process = subprocess.run([sys.executable, str(campaign / 'evaluator/tournament.py'),
        'verify', '--round', str(campaign)], capture_output=True, text=True, timeout=300)
    require(process.returncode == 0, process.stdout + process.stderr)
    report = json.loads(process.stdout)
    require(report == dict(status='VERIFIED', trial_receipts=102, source_files=969),
            'incomplete frozen replay')
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    require(sys.version_info[:2] == (3, 12), 'use Python 3.12 for frozen summaries')
    with tempfile.TemporaryDirectory(prefix='ic-adapter-archive-') as temporary:
        root = restore(Path(temporary))
        audit = replay(root)
        result = control_result(root)
        require(result == read(HERE / 'RESULTS.json'), 'changed result export')
    with args.out.open('x') as stream:
        json.dump(dict(status='VERIFIED', frozen_replay=audit, result=result), stream, indent=2)
        stream.write('\n')
    print(json.dumps(dict(status='VERIFIED', trial_receipts=102, verified_pairs=90,
                         retained_incomplete_profiles=12, promotion_eligible=False)))


if __name__ == '__main__':
    main()
