#!/usr/bin/env python3
"""Replay the retained diagnostic panel without executing any producer."""
import argparse
from collections import Counter
import hashlib
import io
import json
import os
from pathlib import Path, PurePosixPath
import subprocess
import sys
import tarfile
import tempfile
import time

from frozen_sat_runtime import SOLVER_BUILD_BUNDLE, verified_sources
from campaign_rules import IC_SOURCE
from identity import curve_record, sha256, write_immutable
from oracle import require
from producer.evidence import check_build_identity
from register_paired_generic import identities
from run_local_pairinv import audit_ic, audit_rho
from run_paired_generic import audit_f5, audit_rho as audit_generic_rho
from tournament import read

HERE = Path(__file__).resolve().parent
EVIDENCE = HERE/'goal_20260924/paired-fresh-n17a1/results-20260929'


def retained_files(evidence):
    evidence = Path(evidence)
    receipt = read(evidence/'receipt.json')
    data = (evidence/receipt['archive_file']).read_bytes()
    require(len(data) == receipt['archive_bytes']
            and hashlib.sha256(data).hexdigest() == receipt['archive_sha256'],
            'paired evidence archive changed')
    files = {}
    with tarfile.open(fileobj=io.BytesIO(data), mode='r:gz') as tar:
        for item in tar:
            name = PurePosixPath(item.name)
            require(item.isfile() and not name.is_absolute()
                    and '..' not in name.parts and item.name not in files,
                    'unsafe or duplicate paired evidence member')
            files[item.name] = tar.extractfile(item).read()
    require({name: dict(bytes=len(value),
                        sha256=hashlib.sha256(value).hexdigest())
             for name, value in files.items()} == receipt['inventory'],
            'paired evidence inventory changed')
    require(receipt['promotion_eligible'] is False
            and receipt['online_speedup'] is None,
            'diagnostic evidence claims unsupported promotion')
    return files


def replay_sat(version, run, evidence):
    sources, _ = verified_sources()
    auditor_receipt = read(evidence/'auditor-receipt.json')
    for name, expected in auditor_receipt['files'].items():
        data = (evidence/'audit-source'/name).read_bytes()
        require(hashlib.sha256(data).hexdigest() == expected,
                'recovered independent SAT auditor changed')
        sources[(Path('research/ic_candidate_tournament_20260915')/name)
                .as_posix()] = data
    suffix = '' if version == 'v1' else '_v2'
    script = '''import importlib,json,sys
from pathlib import Path
module=importlib.import_module('audit_static_sat_full'+sys.argv[1])
print(json.dumps(module.audit(Path(sys.argv[2])),sort_keys=True))
'''
    with tempfile.TemporaryDirectory(prefix='paired-sat-audit-') as temporary:
        root = Path(temporary)
        for name, data in sources.items():
            path = root/name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(data)
        directory = root/'research/ic_candidate_tournament_20260915'
        for name in ('goal_20260924', 'evidence'):
            (directory/name).symlink_to(HERE/name, target_is_directory=True)
        bundle = root/SOLVER_BUILD_BUNDLE
        bundle.parent.mkdir(parents=True, exist_ok=True)
        bundle.symlink_to(HERE.parents[1]/SOLVER_BUILD_BUNDLE,
                          target_is_directory=True)
        environment = {key: value for key, value in os.environ.items()
                       if key not in ('PYTHONPATH', 'PYTHONSTARTUP')}
        process = subprocess.run(
            [sys.executable, '-c',
             'import sys;sys.path.insert(0,sys.argv.pop(1));'+script,
             str(directory), suffix, str(run)],
            cwd=root, env=environment, capture_output=True, text=True,
            timeout=300)
        require(process.returncode == 0,
                'historical independent SAT audit failed: '+process.stderr[-6000:])
        return json.loads(process.stdout)


def rank_gap(f5, sat):
    """Independent field elimination identifies the missing relation direction."""
    matrix = f5['relation_matrix']
    count, modulus = f5['columns'], int(matrix['modulus'])
    require([list(map(int, point)) for point in sat['matrix']['column_points']]
            == [list(map(int, point)) for point in matrix['column_points']],
            'paired PDP relation columns differ')
    pivots = {}
    for item in matrix['rows']:
        row = [0]*count
        for column, value in item['entries']:
            row[column] = int(value) % modulus
        for column, basis in sorted(pivots.items()):
            if row[column]:
                coefficient = row[column]
                row = [(x-coefficient*y) % modulus for x, y in zip(row, basis)]
        column = next((index for index, value in enumerate(row) if value), None)
        if column is not None:
            inverse = pow(row[column], -1, modulus)
            pivots[column] = [value*inverse % modulus for value in row]
    free = [column for column in range(count) if column not in pivots]
    require(len(free) == 1, 'diagnostic F5 rank shortfall changed')
    vector = [0]*count
    vector[free[0]] = 1
    for column, basis in sorted(pivots.items(), reverse=True):
        vector[column] = -sum(x*y for x, y in zip(basis, vector)) % modulus
    require(all(sum(int(value)*vector[column] for column, value in item['entries'])
                % modulus == 0 for item in matrix['rows']),
            'independent F5 nullspace does not annihilate every row')
    by_trial = {item['trial']: item for batch in f5['collection_reports']
                for item in batch['attempts']}
    paired = Counter()
    by_scalar = {}
    for item in sat['collection']:
        trial = item['trial']
        require(by_trial[trial]['a'] == item['probe_scalar'],
                'paired ordinary query differs')
        by_scalar[item['probe_scalar']] = trial
        paired[(by_trial[trial]['pdp']['outcome'], item['status'])] += 1
    resolving = []
    for row in sat['matrix']['rows']:
        dot = sum(int(value)*vector[column] for column, value in row['entries']) % modulus
        if dot:
            trial = by_scalar[row['scalar']]
            resolving.append(dict(trial=trial, nullspace_dot=dot,
                                  f5_outcome=by_trial[trial]['pdp']['outcome']))
    return dict(f5_rank=len(pivots), columns=count, free_column=free[0],
                nullspace_vector=vector, sat_rows_resolving_gap=resolving,
                paired_query_count=len(sat['collection']),
                paired_statuses=[dict(f5=first, sat=second, count=value)
                                 for (first, second), value in sorted(paired.items())])


def replay(evidence=EVIDENCE):
    evidence = Path(evidence)
    files = retained_files(evidence)
    fixtures = [json.loads(files[alias+'/stdout.json'])['fixture']
                for alias in ('f5', 'incumbent', 'rho-selected', 'rho-generic')]
    require(all(curve_record(fixture) == curve_record(fixtures[0])
                and fixture['targets'] == fixtures[0]['targets']
                and fixture['target_seeds'] == fixtures[0]['target_seeds']
                and fixture['target_scalar_constructed'] is False
                for fixture in fixtures),
            'paired arm changes curve, supplied point or fixture provenance')
    with tempfile.TemporaryDirectory(prefix='paired-n17-evidence-') as temporary:
        root = Path(temporary)
        for name, data in files.items():
            path = root/name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(data)
        sat_audits = {}
        for version in ('v1', 'v2'):
            sat_audits[version] = replay_sat(version, root/('sat-'+version), evidence)
            require(sat_audits[version] == read(root/'reviews'/('sat-'+version+'-audit.json')),
                    'retained SAT audit does not reproduce')
        require(read(root/'sat-v2/registered-panel.json')['target_input']['point']
                    == list(map(int, fixtures[0]['targets'][0]))
                and read(root/'sat-v2/candidate.json')['record']['curve']['curve_id']
                    == curve_record(fixtures[0])['curve']['curve_id'],
                'SAT does not solve the same paired public point and curve')
        data = identities()
        f5 = read(root/'f5/stdout.json')
        generic_audits = {}
        original_external_audit_times = {}
        for arm, alias in (('f5', 'f5'), ('rho', 'rho-generic')):
            directory = root/alias
            report, job, process = map(read, (directory/'stdout.json',
                directory/'executed-job.json', directory/'process.json'))
            audited = (audit_f5(report, job, process, data, directory/'worker')
                       if arm == 'f5' else
                       audit_generic_rho(report, job, process, data, directory/'worker'))
            original = read(directory/'admission.json')
            if arm == 'f5':
                # Repeating an external audit has its own elapsed time. Only
                # that observation is excluded from exact receipt equality.
                original_time = original['run']['independent_audit_wall_ns']
                replay_time = audited['run']['independent_audit_wall_ns']
                require(type(original_time) is int and original_time >= 0
                        and type(replay_time) is int and replay_time >= 0,
                        'generic independent audit elapsed time invalid')
                original_external_audit_times[arm] = original_time
                original = original | {'run': original['run'] | {
                    'independent_audit_wall_ns': replay_time}}
            else:
                original_external_audit_times['rho-generic'] = None
            require(audited == original,
                    'retained generic admission does not reproduce')
            generic_audits[arm] = audited
        native_audits = {}
        for arm, alias in (('ic', 'incumbent'), ('rho', 'rho-selected')):
            directory = root/alias
            build, panel, seal = map(read, (directory/'build-record.json',
                directory/'panel.json', directory/'seal.json'))
            require(sha256(build) == panel['build_record_sha256']
                    and sha256(build['source_manifest']) == IC_SOURCE
                    and build['accepted_source_manifest_sha256'] == IC_SOURCE
                    and hashlib.sha256((directory/'worker').read_bytes()).hexdigest()
                        == build['worker_sha256']
                    and hashlib.sha256((directory/'registered-runner.py').read_bytes()).hexdigest()
                        == seal['controller_sha256']
                    and hashlib.sha256((directory/'registered-registration.py').read_bytes()).hexdigest()
                        == seal['registration_sha256']
                    and all(hashlib.sha256((directory/(key.replace('_','-')+'.json'))
                        .read_bytes()).hexdigest() == value
                            for key, value in seal['files'].items()),
                    'native source/binary/registration binding changed')
            require(all(hashlib.sha256(files['pairinv-source/'+name]).hexdigest()
                        == expected for name, expected in build['source_manifest'].items()),
                    'retained native source differs from accepted build manifest')
            bound = dict(inventory=read(directory/'inventory/stdout.json'),
                         job=read(directory/'executed-job.json'),
                         method=read(directory/'method.json'),
                         candidate=read(directory/'candidate.json'))
            report, process = read(directory/'stdout.json'), read(directory/'process.json')
            check_build_identity(report, IC_SOURCE, panel['field_kernel'])
            audited = (audit_ic(report, process, bound) if arm == 'ic' else
                       audit_rho(report, process, bound))
            require(audited == read(directory/'admission.json'),
                    'retained native admission does not reproduce')
            native_audits[arm] = audited
        gap = rank_gap(f5, read(root/'sat-v2/summary.json'))
    return dict(schema_version=1, status='AUDITED_LOCAL_DIAGNOSTIC_PANEL',
                evidence_archive_sha256=read(evidence/'receipt.json')['archive_sha256'],
                sat_audits=sat_audits,
                f5_audit_status=generic_audits['f5']['status'],
                f5_report_status=generic_audits['f5']['run']['status'],
                native_verified_targets={arm: audit['certificate']['verified_targets']
                                         for arm, audit in native_audits.items()},
                generic_rho_verified_targets=generic_audits['rho']['certificate']['verified_targets'],
                same_point_rho_audited=True,
                rank_gap=gap, complete_sat_preexecution_python_coverage=False,
                calibrated_host_conditions=False, arm_order_followed=False,
                original_external_audit_timing_ns=original_external_audit_times | {
                    'sat-v1': None, 'sat-v2': None,
                    'incumbent': None, 'rho-selected': None},
                promotion_eligible=False, headline_online_admissible=False,
                online_speedup=None)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--evidence', type=Path, default=EVIDENCE)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    started = time.monotonic_ns()
    result = replay(args.evidence)
    result['posthoc_replay_wall_ns'] = time.monotonic_ns()-started
    write_immutable(args.out, result)
    print(json.dumps({key: result[key] for key in
                      ('status', 'online_speedup', 'promotion_eligible')}, sort_keys=True))
