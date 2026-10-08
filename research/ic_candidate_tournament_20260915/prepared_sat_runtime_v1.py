"""New one-use source-bound SAT target adapter on certified n17 preparation.

Only a disclosed development control is admitted by this version. It cannot
sample or qualify fresh targets, extend the old SAT invocation or assert a win.
Preparation evidence enters invocation data; method identity binds math state
and executed sources, never a hash of old preparation seeds or measurements.
"""
import argparse
from collections import Counter
import copy
import hashlib
import json
from pathlib import Path
import platform

from audit_static_sat_full_v3 import verify_query
from generic_query_law import descent_coefficients
from identity import candidate_manifest, run_id, sha256, workload_manifest
from oracle import require
from prepared_target_v1 import (ONLINE_PHASES, load_preparation,
                               solve_sat_target, supplied_target)
from run_generic_exact_yield_audit import pair_index
from sat_runtime_execution_v3 import (audit_execution, execute, read,
                                     register as register_runtime)
from static_sat_assets_v3 import check_extracted_assets, verified_assets
from static_sat_inputs_v3 import native_admission
from static_sat_native_v3 import audit_meter, meter
from static_sat_pipeline_v3 import save_progress
from static_sat_query_v3 import one_query
from static_sat_registration_v2 import method_record
from tournament import write

SETTINGS = {'cms_conflict_budget', 'cms_timeout_seconds', 'export_timeout_seconds',
            'max_descent_queries', 'descent_query_seed', 'export_nonce', 'target_input',
            'resources', 'question', 'run_number'}


def mathematical_registration(panel, spec, files, document, certificate_sha256):
    fixture, inventory, curve, base, native = native_admission(files)
    _, imported_base, _, _, preparation = load_preparation(document, certificate_sha256, family='sat')
    require(tuple(base) == imported_base and set(panel) == SETTINGS,
            'prepared SAT native geometry or panel fields differ')
    require(panel['question'] == 'prepared-development-source-control',
            'prepared SAT fresh qualification requires the new reviewed paired protocol')
    for key in SETTINGS-{'target_input', 'resources', 'question'}:
        require(type(panel[key]) is int and 0 <= panel[key] < 2**64,
                'prepared SAT setting is not an exact bounded integer: '+key)
    require(1 <= panel['max_descent_queries'] <= 8 and panel['cms_conflict_budget'] == 100000
            and 1 <= panel['cms_timeout_seconds'] <= 60
            and 1 <= panel['export_timeout_seconds'] <= 60,
            'prepared SAT control exceeds its frozen native limits')
    target_input = panel['target_input']
    require(target_input == dict(point=[52411,72106], seed=None,
                input_law='one-disclosed-public-point; fixture-construction-excluded',
                point_was_previously_supplied=True, known_scalar_supplied=False)
            and target_input['point_was_previously_supplied'] is True
            and target_input['known_scalar_supplied'] is False,
            'prepared SAT control requires the already disclosed public point without a scalar')
    target = supplied_target(curve, target_input['point'])
    resources = panel['resources']
    require(resources == dict(host_class='physical-macos-arm64', cpu_workers=1,
                target_count=1, memory_limit_bytes=None,
                total_wall_limit_seconds=resources.get('total_wall_limit_seconds'))
            and type(resources['cpu_workers']) is int and type(resources['target_count']) is int
            and type(resources['total_wall_limit_seconds']) is int
            and 30 <= resources['total_wall_limit_seconds'] <= 900,
            'prepared SAT requires the bounded one-worker development envelope')
    require(spec['entrypoint'] == dict(module='prepared_sat_runtime_v1', callable='run')
            and spec['runtime_watchdog_seconds'] == resources['total_wall_limit_seconds'],
            'prepared SAT source entrypoint or watchdog differs')
    source = dict(execution_binding=spec['binding'], native=native)
    # Build logs/receipts and asset transport are evidence, not method. Bind
    # every executed component without hashing their incidental build history
    # into candidate identity. The full evidence remains sealed in arguments.
    code_binding = {key:value for key,value in spec['binding'].items()
                    if not key.startswith('asset_')}
    code_binding['native_components'] = dict(
        exporter_binary_sha256=native['exporter_binary_sha256'],
        exporter_source_sha256=native['exporter_source_sha256'],
        rust_manifest_sha256=sha256(json.loads(files['rust/source-manifest.json'])),
        cms_binary_sha256=native['cms_binary_sha256'],
        cms_source_sha256=hashlib.sha256(files['cms/source.tar']).hexdigest(),
        cadical_source_sha256=hashlib.sha256(files['cms/cadical.tar']).hexdigest(),
        cadiback_source_sha256=hashlib.sha256(files['cms/cadiback.tar']).hexdigest())
    parameters = dict(panel, max_relation_queries=0, cms_max_models_per_query=1,
        cms_executable_sha256=native['cms_binary_sha256'],
        cms_build_receipt_sha256=native['cms_build_receipt_sha256'],
        exporter_source_sha256=native['exporter_source_sha256'],
        source_encoding='wide-symmetrised-S4-circuit-XOR-DIMACS')
    method = method_record(parameters, sha256(code_binding), fixture)
    del method['implementation']['flags']['cms_build_receipt_sha256']
    method['relation_collection'].update(
        query_distribution='reusable-target-independent-preparation',
        query_rule='import-certified-ordinary-rows; no new ordinary queries',
        stop_rule='independent-full-rank-and-all-column-log-replay-before-online')
    method['implementation']['flags'].update(
        algorithm_execution_binding=code_binding, mathematical_preparation_sha256=document['record_sha256'],
        preparation_policy='replay-all-ordinary-rows-and-logs-before-online-v1',
        online_interval='first-target-query-through-independent-scalar-replay',
        native_thread_environment='all-listed-pools-one-thread; CMS-one-thread',
        native_watchdog_group='inherit-controller-group-no-native-fork-or-setsid',
        observer_cost='all-target-native-wrapper-source-checks-charged-to-PDP')
    candidate = candidate_manifest(fixture, inventory, method)
    target_fixture = dict(fixture, targets=[list(target)], target_seeds=[None])
    workload = workload_manifest(target_fixture, input_law=target_input['input_law'],
        algorithm_seed=panel['descent_query_seed'], resource_envelope=resources, cache_policy='warm')
    workload['record'].update(algorithm_seeds=dict(descent=panel['descent_query_seed'], exporter=panel['export_nonce']),
                              question=panel['question'], point_was_previously_supplied=True)
    workload['record_sha256'] = sha256(workload['record'])
    workload['workload_id'] = workload['record_sha256'][:12]
    seal = dict(registration_stage='before-execution', panel_sha256=sha256(panel),
        candidate_id=candidate['candidate_id'], workload_id=workload['workload_id'],
        run_id=run_id(candidate['candidate_id'], workload['workload_id'], panel['run_number']),
        candidate_sha256=sha256(candidate), workload_sha256=sha256(workload),
        source_sha256=sha256(source), method_sha256=sha256(method),
        preparation_certificate_sha256=certificate_sha256)
    return dict(panel=copy.deepcopy(panel), source=source, method=method, candidate=candidate,
        workload=workload, seal=seal, preparation=copy.deepcopy(document),
        preparation_certificate_sha256=certificate_sha256, preparation_receipt=preparation)


def register(repository, assets, panel, document, certificate_sha256, output):
    assets = Path(assets)
    files = verified_assets(assets, read(assets/'manifest.json'), read(assets/'seal.json'))
    return register_runtime(repository, output, module='prepared_sat_runtime_v1', action='run',
        arguments=None, timeout_seconds=panel['resources']['total_wall_limit_seconds'], asset_snapshot=assets,
        arguments_factory=lambda spec:mathematical_registration(panel, spec, files, document, certificate_sha256))


def run(arguments, out):
    out = Path(out).resolve()
    execution = out.parent
    spec = read(execution/'execution.json')
    require(out == execution/'entry-output' and not (out/'summary.json').exists(),
            'prepared SAT target output reused or outside its sealed execution')
    files = check_extracted_assets(execution/'asset-files', spec['asset_manifest'])
    panel, document, certificate_sha = (arguments[key] for key in
        ('panel', 'preparation', 'preparation_certificate_sha256'))
    require(mathematical_registration(panel, spec, files, document, certificate_sha) == arguments
            and spec['arguments'] == arguments, 'prepared SAT mathematical invocation changed')
    native = native_admission(files)[-1]
    require(native['platform'] == dict(system=platform.system(), machine=platform.machine()),
            'prepared SAT native binaries require their validated physical platform')
    for name in ('panel', 'source', 'method', 'candidate', 'workload', 'seal'):
        write(out/(name+'.json'), arguments[name], exclusive=True)
    preflight = meter(execution, 'bin/cms', ['--version'], out, 'cms_preflight', 10)
    require(not preflight['timed_out'] and preflight['returncode'] == 0
            and 'CryptoMiniSat version 5.14.7' in (out/'cms_preflight.stdout').read_text(),
            'prepared SAT pinned native solver failed version preflight')
    (out/'target').mkdir()

    def query(settings, item, curve, base):
        return one_query(settings, item, execution, curve, base, out/'target')

    try:
        result = solve_sat_target(document, certificate_sha, point=panel['target_input']['point'],
            panel=panel, query=query,
            progress=lambda row:save_progress(out, 'target.progress.jsonl', row))
        # Final reporting follows the scalar-replay endpoint. Intermediate
        # failure reports are already charged by the adapter core.
        if result['scalar_verified']:
            save_progress(out, 'target.progress.jsonl', result['target_attempts'][-1])
        result.update(candidate_id=arguments['seal']['candidate_id'],
            workload_id=arguments['seal']['workload_id'], run_id=arguments['seal']['run_id'],
            cms_preflight=preflight)
        write(out/'summary.json', result, exclusive=True)
        return dict(status=result['status'], run_id=result['run_id'], scalar_verified=result['scalar_verified'])
    except Exception as error:
        write(out/'failure.json', dict(status='ERROR', error_type=type(error).__name__, error=str(error),
            run_id=arguments['seal']['run_id'], online_wall_ns=None, online_speedup=None,
            progress='target.progress.jsonl; native trial directories retained'), exclusive=True)
        raise


def audit(execution, expected_spec):
    """No solver executes; replay source gates, source models and target rows."""
    execution = Path(execution)
    source_audit = audit_execution(execution, expected_spec)
    require(source_audit['entrypoint_succeeded'], 'partial prepared SAT execution remains a failure')
    files = check_extracted_assets(execution/'asset-files', expected_spec['asset_manifest'])
    arguments = expected_spec['arguments']
    panel, document, certificate_sha = (arguments[key] for key in
        ('panel', 'preparation', 'preparation_certificate_sha256'))
    require(mathematical_registration(panel, expected_spec, files, document, certificate_sha) == arguments,
            'audited prepared SAT registration differs')
    root = execution/'entry-output'
    result = read(root/'summary.json')
    require(result['status'] in ('COMPLETE', 'INCOMPLETE_TARGET')
            and result['ordinary_queries_executed'] == 0
            and result['candidate_id'] == arguments['seal']['candidate_id']
            and result['workload_id'] == arguments['seal']['workload_id']
            and result['run_id'] == arguments['seal']['run_id']
            and result['promotion_eligible'] is False and result['online_speedup'] is None
            and result['fresh_paired_qualification'] is False
            and result['source_bound_execution_admitted'] is False,
            'prepared SAT summary altered its admitted identity or scientific scope')
    curve, base, matrix, logs, prep = load_preparation(document, certificate_sha, family='sat')
    require(result['preparation'] == prep and result['target_input'] == panel['target_input']['point'],
            'prepared SAT result imported a different preparation or target')
    target = supplied_target(curve, result['target_input'])
    rows = result['target_attempts']
    require([json.loads(line) for line in (root/'target.progress.jsonl').read_text().splitlines()] == rows,
            'prepared SAT progress dropped or changed target attempts')
    require(0 < len(rows) <= panel['max_descent_queries'], 'prepared SAT target attempt budget changed')
    stream = descent_coefficients(panel['descent_query_seed'], curve.r, walked=False)
    pairs = pair_index(curve, base)
    scalar = None
    for trial, row in enumerate(rows):
        a,b = next(stream)
        require(row['trial'] == trial and row['target_query_index'] == trial
                and row['a'] == a and row['b'] == b and row['probe_scalar'] == a,
                'prepared SAT target query chronology or coefficients differ')
        point = curve.add(curve.mul(curve.g,a), curve.mul(target,b))
        require(row['public_point'] == (None if point is None else list(point)),
                'prepared SAT target query point differs')
        _, witness, status = verify_query(execution, panel, root, f'target/trial-{trial:02d}/',
                                          row, point, curve, base, pairs)
        require(type(row['query_wall_ns']) is int and type(row['verification_wall_ns']) is int
                and row['query_wall_ns'] >= row['verification_wall_ns'] >= 0
                and row['pdp_wall_ns'] == row['query_wall_ns']-row['verification_wall_ns'],
                'prepared SAT query split overlaps')
        if status == 'VALID_POINT_WITNESS':
            scalar = matrix.recover_target(logs, target, a,b,witness)
            require(trial+1 == len(rows) and row['candidate_scalar'] == str(scalar)
                    and row['scalar_replay_verified'] is True, 'prepared SAT scalar or stop event differs')
    complete = scalar is not None
    require(type(result['scalar_verified']) is bool
            and complete == (result['status'] == 'COMPLETE') == result['scalar_verified']
            and result['recovered_scalar'] == scalar, 'prepared SAT completion/scalar changed')
    require(complete or len(rows) == panel['max_descent_queries'], 'incomplete SAT target stopped before cap')
    require(result['target_status_mix'] == dict(Counter(row['status'] for row in rows))
            and result['online_stop_event'] == ('independent-scalar-replay' if complete else 'frozen-target-attempt-cap')
            and type(result['online_bookkeeping_assigned_to_target_query_ns']) is int
            and result['online_bookkeeping_assigned_to_target_query_ns'] >= 0,
            'prepared SAT status census or endpoint differs')
    phases = result['online_phases_ns']
    require(set(phases) == set(ONLINE_PHASES)
            and all(type(v) is int and v >= 0 for v in phases.values()), 'prepared SAT target phase slots changed')
    attempted = result['online_attempt_wall_ns']
    require(type(attempted) is int and attempted > 0 and sum(phases.values()) == attempted
            and result['online_stop_monotonic_ns']-result['online_start_monotonic_ns'] == attempted
            and result['online_wall_ns'] == (attempted if complete else None)
            and phases['target_pdp'] == sum(row['pdp_wall_ns'] for row in rows)
            and phases['target_relation_check'] >= sum(row['verification_wall_ns'] for row in rows)
            and phases['target_query'] >= result['online_bookkeeping_assigned_to_target_query_ns'],
            'prepared SAT online interval or charged PDP attempts differ')
    preflight = audit_meter(execution, root, 'cms_preflight', asset_role='bin/cms', arguments=['--version'], seconds=10)
    require(preflight == result['cms_preflight'], 'prepared SAT preflight receipt differs')
    return dict(status='ADMITTED_COMPLETE_PREPARED_SAT_CONTROL' if complete else 'ADMITTED_INCOMPLETE_PREPARED_SAT_CONTROL',
        candidate_id=result['candidate_id'], workload_id=result['workload_id'], run_id=result['run_id'],
        source_audit=source_audit, source_bound_execution_admitted=True,
        scalar_verified=complete, recovered_scalar=scalar, target_attempt_count=len(rows),
        online_wall_ns=result['online_wall_ns'], online_attempt_wall_ns=attempted,
        online_phases_ns=phases, fresh_paired_qualification=False, headline_online_admissible=False,
        promotion_eligible=False, online_speedup=None)


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
        require(sha256(spec) == args.expected_execution_sha256, 'prepared SAT externally frozen invocation differs')
        receipt = execute(args.registration, args.out, expected_spec=spec,
                          timeout_seconds=spec['runtime_watchdog_seconds'])
        print(json.dumps(receipt, sort_keys=True))


if __name__ == '__main__':
    main()
