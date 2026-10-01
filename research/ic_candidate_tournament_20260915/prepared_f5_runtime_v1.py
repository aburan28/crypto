"""Source-bound registration for one prepared n17 F5 target control.

The retained macOS worker from the consumed F5 v2 solve does not contain the
prepared-target entrypoint. This registrar admits only a different build whose
worker source contains that entrypoint, and only the already disclosed point.
It does not sample, dispatch, or qualify a fresh target. Build logs and asset
seals stay in the invocation; method identity binds mathematical preparation
and executed code.
"""
import argparse
import copy
import hashlib
import io
import json
from pathlib import Path
import platform
import tarfile

from f5_runtime_inputs_v1 import (SOURCE_MANIFEST as CONSUMED_SOURCE_MANIFEST,
                                 WORKER_SHA256 as CONSUMED_WORKER_SHA256)
from generic_admission import method_record
from generic_bases import verify_base
from generic_build import verify_build_record
from generic_solver_feasibility import assess as assess_layout
from generic_stages import DEFAULTS, effective_config
from identity import candidate_manifest, run_id, sha256, workload_manifest
from oracle import require
from prepared_ic_state_v1 import STATE_SHA256
from prepared_target_v1 import audit_native_target, load_preparation, native_job
from sat_runtime_execution_v3 import audit_execution, binding, execute, read, register as register_runtime
from static_sat_assets_v3 import check_extracted_assets, freeze_assets, verified_assets
from static_sat_native_v3 import audit_meter, meter
from tournament import write

QUESTION = 'prepared-development-source-control'
INPUT_LAW = 'one-disclosed-public-point; fixture-construction-excluded'
POINT = [52411, 72106]
MARKER = b'fn run_prepared_target'
WORKER_ROLE = 'examples/ic_tournament_worker.rs'
SETTINGS = {'question', 'resources', 'target_input', 'algorithm_seed',
            'max_descent_queries', 'run_number'}
ROLES = {'bin/worker', 'build/build-record.json', 'build/build-policy.json',
         'build/build-exit.json', 'build/build.log', 'rust/source-manifest.json',
         'rust/root-source.tar.gz', 'rust/dependency-source.tar.gz',
         'fixture.json', 'inventory.json'}
PLATFORM = {('linux', 'x86_64'): ('Linux', 'x86_64'),
            ('linux', 'aarch64'): ('Linux', 'aarch64'),
            ('macos', 'aarch64'): ('Darwin', 'arm64'),
            ('macos', 'x86_64'): ('Darwin', 'x86_64')}


def digest(data):
    return hashlib.sha256(data).hexdigest()


def source_archive(data, expected):
    observed = {}
    with tarfile.open(fileobj=io.BytesIO(data), mode='r:gz') as archive:
        for item in archive:
            require(item.isfile() and item.name not in observed,
                    'non-file or duplicate prepared F5 source member')
            observed[item.name] = digest(archive.extractfile(item).read())
    require(observed == expected, 'prepared F5 source bytes differ from the build manifest')


def admit_native(files):
    """Check a new prepared-capable build. Reject the consumed v2 worker."""
    require(set(files) == ROLES, 'prepared F5 native input roles missing or extra')
    record = json.loads(files['build/build-record.json'])
    source = json.loads(files['rust/source-manifest.json'])
    identity = verify_build_record(record, source)
    worker_sha = digest(files['bin/worker'])
    source_sha = sha256(source)
    require(worker_sha == record['worker_sha256'] != CONSUMED_WORKER_SHA256
            and source_sha == record['source_manifest_sha256'] != CONSUMED_SOURCE_MANIFEST,
            'prepared F5 refuses the consumed v2 worker and source manifest')
    require(json.loads(files['build/build-policy.json']) == record['build']
            and json.loads(files['build/build-exit.json']) == {'exit_code': 0},
            'prepared F5 build policy or exit receipt differs')
    require(WORKER_ROLE in source['root_files'], 'prepared F5 source omits the worker entrypoint')
    source_archive(files['rust/root-source.tar.gz'], source['root_files'])
    with tarfile.open(fileobj=io.BytesIO(files['rust/root-source.tar.gz']), mode='r:gz') as archive:
        worker_source = archive.extractfile(WORKER_ROLE).read()
    require(digest(worker_source) == source['root_files'][WORKER_ROLE] and MARKER in worker_source,
            'prepared F5 worker source lacks the prepared-target entrypoint')
    dependencies = {item['package']+'-'+item['version']+'/'+name: expected
                    for item in source['dependencies'] for name, expected in item['files'].items()}
    source_archive(files['rust/dependency-source.tar.gz'], dependencies)
    target = (record['build']['target_os'], record['build']['target_arch'])
    require(target in PLATFORM, 'prepared F5 build target has no platform binding')
    fixture = json.loads(files['fixture.json'])
    inventory = json.loads(files['inventory.json'])
    require(fixture['degree'] == 17 and fixture['curve_a'] == 1
            and fixture['targets'] == [] and fixture['target_seeds'] == []
            and fixture['target_scalar_constructed'] is False
            and inventory['fixture'] == fixture
            and inventory['columns'] == 29
            and len(inventory['factor_base']) == 63,
            'prepared F5 reusable fixture contains a target or changed geometry')
    probe = dict(mode='ic', degree=17, curve_a=1,
                 factor_base=dict(kind='standard_subspace', dimension=6), config={})
    checked = verify_base(inventory, fixture, probe)
    native = dict(build_identity=identity, source_manifest_sha256=source_sha,
                  worker_sha256=worker_sha, build_sha256=record['build_sha256'],
                  build_record_sha256=sha256(record),
                  worker_source_sha256=source['root_files'][WORKER_ROLE],
                  target_os=target[0], target_arch=target[1],
                  host_class=target[0]+'-'+target[1])
    return fixture, inventory, checked, record, source, native


def validate_panel(panel, native):
    require(set(panel) == SETTINGS and panel['question'] == QUESTION,
            'prepared F5 fresh qualification requires the reviewed paired protocol')
    require(type(panel['algorithm_seed']) is int and 0 <= panel['algorithm_seed'] < 2**64
            and type(panel['run_number']) is int and panel['run_number'] >= 0
            and type(panel['max_descent_queries']) is int
            and 1 <= panel['max_descent_queries'] <= 8,
            'prepared F5 seed, run or attempt cap is outside the frozen control')
    resources = panel['resources']
    require(resources == dict(host_class=native['host_class'], cpu_workers=1, target_count=1,
                memory_limit_bytes=None,
                total_wall_limit_seconds=resources.get('total_wall_limit_seconds'))
            and type(resources['cpu_workers']) is int and type(resources['target_count']) is int
            and type(resources['total_wall_limit_seconds']) is int
            and 30 < resources['total_wall_limit_seconds'] <= 900,
            'prepared F5 requires the one-worker build host and a bounded watchdog')
    target_input = panel['target_input']
    require(target_input == dict(point=POINT, seed=None, input_law=INPUT_LAW,
                point_was_previously_supplied=True, known_scalar_supplied=False),
            'prepared F5 control requires the already disclosed public point without a scalar')
    return resources


def mathematical_registration(panel, spec, files, document, certificate_sha256):
    fixture, inventory, checked, record, _source, native = admit_native(files)
    resources = validate_panel(panel, native)
    _curve, imported_base, _matrix, _logs, preparation = load_preparation(
        document, certificate_sha256, family='f5')
    require([[str(x), str(y)] for x, y in imported_base] == inventory['factor_base']
            and document['record_sha256'] == STATE_SHA256,
            'prepared F5 certificate geometry differs from the admitted base')
    require(spec['entrypoint'] == dict(module='prepared_f5_runtime_v1', callable='run')
            and spec['runtime_watchdog_seconds'] == resources['total_wall_limit_seconds']
            and 'asset_manifest_sha256' in spec['binding']
            and 'asset_archive_sha256' in spec['binding'],
            'prepared F5 entrypoint, watchdog or asset seal differs')
    job = native_job(document, certificate_sha256, point=POINT,
                     algorithm_seed=panel['algorithm_seed'],
                     max_attempts=panel['max_descent_queries'])
    layout = assess_layout(dict(cells=['n17a1'], candidates=[dict(
        id='prepared-development-arm',
        config=dict(job['config'], factor_base=job['factor_base']))]))
    require(layout['status'] == 'PASS_STATIC_LAYOUT_ONLY'
            and layout['rows'][0]['boolean_variables'] <= layout['max_boolean_variables'],
            'prepared F5 encoder layout exceeds the static cap')
    declared = copy.deepcopy(inventory)
    declared['effective_config'] = dict(copy.deepcopy(DEFAULTS), **job['config'])
    effective_config(job, declared)
    method = method_record(job, declared, {'base': checked}, record,
                           collector_plan=inventory['collector_dispatch'])
    code_binding = {key: value for key, value in spec['binding'].items()
                    if not key.startswith('asset_')}
    code_binding['native_components'] = dict(
        worker_sha256=native['worker_sha256'],
        source_manifest_sha256=native['source_manifest_sha256'],
        build_sha256=native['build_sha256'],
        worker_source_sha256=native['worker_source_sha256'])
    method['relation_collection']['query_distribution'] = 'reusable-target-independent-preparation'
    method['relation_collection']['query_rule'] = dict(
        policy='import-certified-ordinary-rows; no new ordinary queries', window='none')
    method['relation_collection']['stop_rule'] = dict(
        batch_trials=job['config']['batch_trials'], max_trials=job['config']['max_trials'],
        success='independent-full-rank-and-all-column-log-replay-before-online')
    method['implementation']['components'].append(dict(
        role='complete-frozen-Python-interpreter-native-binding', sha256=sha256(code_binding)))
    method['implementation']['flags'].update(
        algorithm_execution_binding=code_binding,
        mathematical_preparation_sha256=document['record_sha256'],
        preparation_policy='replay-all-ordinary-rows-and-logs-before-online-v1',
        online_interval='first-target-query-through-independent-scalar-replay',
        native_thread_environment='all-listed-pools-one-thread',
        native_watchdog_group='inherit-controller-group-no-native-fork-or-setsid',
        observer_cost='target-native-wrapper-and-source-checks-charged-online',
        stdin_policy='canonical-registered-prepared-job-UTF8')
    candidate = candidate_manifest(fixture, declared, method)
    target_fixture = dict(fixture, targets=[POINT], target_seeds=[None])
    workload = workload_manifest(target_fixture, input_law=INPUT_LAW,
        algorithm_seed=panel['algorithm_seed'], resource_envelope=resources, cache_policy='warm')
    workload['record'].update(question=panel['question'], point_was_previously_supplied=True)
    workload['record_sha256'] = sha256(workload['record'])
    workload['workload_id'] = workload['record_sha256'][:12]
    static_layout = dict(status=layout['status'],
        boolean_variables=layout['rows'][0]['boolean_variables'],
        max_boolean_variables=layout['max_boolean_variables'],
        scope='encoder preflight only; not a claim that this build is an older reviewed worker')
    source_row = dict(schema_version=1, execution_binding=spec['binding'], native=native)
    seal = dict(registration_stage='before-execution', panel_sha256=sha256(panel),
        candidate_id=candidate['candidate_id'], workload_id=workload['workload_id'],
        run_id=run_id(candidate['candidate_id'], workload['workload_id'], panel['run_number']),
        candidate_sha256=sha256(candidate), workload_sha256=sha256(workload),
        source_sha256=sha256(source_row), method_sha256=sha256(method),
        preparation_certificate_sha256=certificate_sha256)
    return dict(panel=copy.deepcopy(panel), job=job, source=source_row, method=method,
        candidate=candidate, workload=workload, seal=seal, static_layout=static_layout,
        preparation=copy.deepcopy(document), preparation_certificate_sha256=certificate_sha256,
        preparation_receipt=preparation)


def package_assets(files, output):
    """Seal an already admitted build directory. Does not compile or run it."""
    admit_native(files)
    manifest, seal = freeze_assets(files, {'bin/worker'}, output)
    return dict(manifest_sha256=sha256(manifest), archive_sha256=seal['archive_sha256'])


def register(repository, assets, panel, document, certificate_sha256, output):
    assets = Path(assets)
    files = verified_assets(assets, read(assets/'manifest.json'), read(assets/'seal.json'))
    return register_runtime(repository, output, module='prepared_f5_runtime_v1', action='run',
        arguments=None, timeout_seconds=panel['resources']['total_wall_limit_seconds'],
        asset_snapshot=assets,
        arguments_factory=lambda spec: mathematical_registration(
            panel, spec, files, document, certificate_sha256))


def expected_platform(native):
    return PLATFORM[(native['target_os'], native['target_arch'])]


def run(arguments, out):
    out = Path(out).resolve()
    execution = out.parent
    spec = read(execution/'execution.json')
    require(out == execution/'entry-output' and not (out/'summary.json').exists(),
            'prepared F5 target output reused or outside its sealed execution')
    files = check_extracted_assets(execution/'asset-files', spec['asset_manifest'])
    panel, document, certificate_sha = (arguments[key] for key in
        ('panel', 'preparation', 'preparation_certificate_sha256'))
    require(mathematical_registration(panel, spec, files, document, certificate_sha) == arguments
            and spec['arguments'] == arguments, 'prepared F5 mathematical invocation changed')
    native = arguments['source']['native']
    require((platform.system(), platform.machine()) == expected_platform(native),
            'prepared F5 worker must run on the platform it was built for')
    for name in ('panel', 'job', 'source', 'method', 'candidate', 'workload', 'seal'):
        write(out/(name+'.json'), arguments[name], exclusive=True)
    preflight = meter(execution, 'bin/worker', ['--build-identity'], out, 'build_identity', 10)
    require(not preflight['timed_out'] and preflight['returncode'] == 0
            and read(out/'build_identity.stdout') == native['build_identity'],
            'prepared F5 worker embeds a different build identity')
    seconds = panel['resources']['total_wall_limit_seconds']-30
    result = meter(execution, 'bin/worker', [], out, 'pipeline', seconds, stdin_argument='job')
    if result['timed_out']:
        write(out/'summary.json', dict(status='NATIVE_TIMEOUT',
            candidate_id=arguments['seal']['candidate_id'], workload_id=arguments['seal']['workload_id'],
            run_id=arguments['seal']['run_id'], scalar_verified=False, recovered_scalar=None,
            ordinary_queries_executed=0, online_wall_ns=None, online_speedup=None,
            source_bound_execution_admitted=False, fresh_paired_qualification=False,
            promotion_eligible=False, build_identity=preflight, pipeline=result), exclusive=True)
        return dict(status='NATIVE_TIMEOUT', run_id=arguments['seal']['run_id'])
    report = read(out/'pipeline.stdout')
    require((result['returncode'] == 0 and report['status'] == 'complete')
            or (result['returncode'] == 2 and report['status'] == 'incomplete'),
            'prepared F5 native process result is inconsistent or failed')
    require(report.get('generic_build') == native['build_identity'],
            'prepared F5 report names a different build')
    audited = audit_native_target(report, arguments['job'], document, certificate_sha)
    summary = dict(status=audited['status'], ordinary_queries_executed=0,
        candidate_id=arguments['seal']['candidate_id'], workload_id=arguments['seal']['workload_id'],
        run_id=arguments['seal']['run_id'], recovered_scalar=audited['recovered_scalar'],
        scalar_verified=audited['scalar_verified'], online_wall_ns=audited['online_wall_ns'],
        online_attempt_wall_ns=audited['online_attempt_wall_ns'],
        online_phases_ns=audited['online_phases_ns'], preparation=audited['preparation'],
        target_input=POINT, promotion_eligible=False, online_speedup=None,
        fresh_paired_qualification=False, source_bound_execution_admitted=False,
        headline_online_admissible=False, build_identity=preflight, pipeline=result)
    write(out/'summary.json', summary, exclusive=True)
    return dict(status=summary['status'], run_id=summary['run_id'],
                scalar_verified=summary['scalar_verified'])


def audit(execution, expected_spec):
    """Replay source gates and the mathematical target audit. No worker runs."""
    execution = Path(execution)
    source_audit = audit_execution(execution, expected_spec)
    require(source_audit['entrypoint_succeeded'], 'partial prepared F5 execution remains a failure')
    files = check_extracted_assets(execution/'asset-files', expected_spec['asset_manifest'])
    arguments = expected_spec['arguments']
    panel, document, certificate_sha = (arguments[key] for key in
        ('panel', 'preparation', 'preparation_certificate_sha256'))
    require(mathematical_registration(panel, expected_spec, files, document, certificate_sha) == arguments,
            'audited prepared F5 registration differs')
    root = execution/'entry-output'
    result = read(root/'summary.json')
    report = read(root/'pipeline.stdout')
    native = arguments['source']['native']
    require(report.get('generic_build') == native['build_identity'],
            'audited prepared F5 report names a different build')
    audited = audit_native_target(report, arguments['job'], document, certificate_sha)
    seconds = panel['resources']['total_wall_limit_seconds']-30
    preflight = audit_meter(execution, root, 'build_identity', asset_role='bin/worker',
                            arguments=['--build-identity'], seconds=10)
    pipeline = audit_meter(execution, root, 'pipeline', asset_role='bin/worker',
                           arguments=[], seconds=seconds, stdin_argument='job')
    require(result == dict(status=audited['status'], ordinary_queries_executed=0,
            candidate_id=arguments['seal']['candidate_id'], workload_id=arguments['seal']['workload_id'],
            run_id=arguments['seal']['run_id'], recovered_scalar=audited['recovered_scalar'],
            scalar_verified=audited['scalar_verified'], online_wall_ns=audited['online_wall_ns'],
            online_attempt_wall_ns=audited['online_attempt_wall_ns'],
            online_phases_ns=audited['online_phases_ns'], preparation=audited['preparation'],
            target_input=POINT, promotion_eligible=False, online_speedup=None,
            fresh_paired_qualification=False, source_bound_execution_admitted=False,
            headline_online_admissible=False, build_identity=preflight, pipeline=pipeline),
            'prepared F5 summary differs from the independent mathematical audit')
    return dict(status='ADMITTED_'+audited['status']+'_PREPARED_F5_CONTROL',
        candidate_id=result['candidate_id'], workload_id=result['workload_id'], run_id=result['run_id'],
        source_audit=source_audit, source_bound_execution_admitted=True,
        scalar_verified=audited['scalar_verified'], recovered_scalar=audited['recovered_scalar'],
        online_wall_ns=result['online_wall_ns'], online_attempt_wall_ns=result['online_attempt_wall_ns'],
        online_phases_ns=result['online_phases_ns'], fresh_paired_qualification=False,
        headline_online_admissible=False, promotion_eligible=False, online_speedup=None)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    freeze = sub.add_parser('register')
    for flag in ('repository', 'assets', 'panel', 'certificate', 'out'):
        freeze.add_argument('--'+flag, type=Path, required=True)
    freeze.add_argument('--expected-certificate-sha256', required=True)
    launch = sub.add_parser('execute')
    launch.add_argument('--registration', type=Path, required=True)
    launch.add_argument('--out', type=Path, required=True)
    launch.add_argument('--expected-execution-sha256', required=True)
    args = parser.parse_args()
    if args.command == 'register':
        spec = register(args.repository, args.assets, read(args.panel), read(args.certificate),
                        args.expected_certificate_sha256, args.out)
        print(json.dumps(dict(execution_sha256=sha256(spec), **spec['arguments']['seal']), sort_keys=True))
    else:
        spec = read(args.registration/'execution.json')
        require(sha256(spec) == args.expected_execution_sha256,
                'prepared F5 externally frozen invocation differs')
        receipt = execute(args.registration, args.out, expected_spec=spec,
                          timeout_seconds=spec['runtime_watchdog_seconds'])
        print(json.dumps(receipt, sort_keys=True))


if __name__ == '__main__':
    main()
