"""Version-three independent F5 admission with exact bounded negative proofs.

No worker executes here. The exact group oracle is postexecution only, and
timings on the development host cannot substantiate a competitive speedup.
"""
import argparse
from collections import Counter
import json
from pathlib import Path

from audit_static_cms_s4_natural import wilson
from f5_runtime_inputs_v1 import native_admission
from f5_runtime_registration_v2 import mathematical_registration
from generic_admission_exact_v1 import admit
from generic_admission import method_record
from identity import write_immutable
from oracle import require
from run_generic_exact_yield_audit import exact_three_sum, pair_index
from sat_runtime_execution_v3 import audit_execution, read
from static_sat_assets_v3 import verified_assets
from static_sat_native_v3 import audit_meter
from f5_runtime_receipts_v2 import publication_process
from generic_build import digest


def natural_queries(report, curve):
    base = [curve.decode(point) for point in report['factor_base']]
    pairs = pair_index(curve, base)
    rows = []
    for batch in report['collection_reports']:
        for attempt in batch['attempts']:
            require(attempt['b'] == 0, 'ordinary F5 query depends on supplied target')
            point = curve.mul(curve.g, attempt['a'])
            indices = exact_three_sum(curve, base, pairs, point)
            pdp = attempt['pdp']
            witness = pdp['outcome'] == 'witness'
            require(not witness or indices is not None, 'native witness has no exact group decomposition')
            if witness:
                total = None
                for index in pdp['points']:
                    require(type(index) is int and 0 <= index < len(base), 'F5 witness index leaves base')
                    total = curve.add(total, base[index])
                require(total == point, 'native F5 witness fails independent group addition')
            stats = pdp['stats']
            detail = stats['stats'] if stats['family'] == 'groebner' else None
            rows.append(dict(trial=attempt['trial'], scalar=attempt['a'], point=None if point is None else list(point),
                pdp_status=pdp['outcome'], exact_group_feasible=indices is not None,
                exact_group_witness_indices=indices, verified_native_witness=witness,
                feasible_but_no_witness=indices is not None and not witness,
                native_stats=detail))
    require(len(rows) == report['trials'], 'F5 ordinary query census differs from charged attempt count')
    return dict(attempts=len(rows), status_mix=dict(Counter(r['pdp_status'] for r in rows)),
        exact_group_feasible=sum(r['exact_group_feasible'] for r in rows),
        verified_native_witnesses=sum(r['verified_native_witness'] for r in rows),
        feasible_but_no_witness=sum(r['feasible_but_no_witness'] for r in rows),
        witness_rate_wilson95=wilson(sum(r['verified_native_witness'] for r in rows), len(rows))
            if rows else None,
        rows=rows, scope='all actual ordinary queries; exact group feasibility only after native execution')


def audit(execution, expected_spec):
    execution = Path(execution)
    source_audit = audit_execution(execution, expected_spec)
    require(source_audit['entrypoint_succeeded'], 'failed F5 entrypoint requires failure retention')
    files = verified_assets(execution/'assets', expected_spec['asset_manifest'], expected_spec['asset_seal'])
    fixture, inventory, curve, base, build, rust_source, native = native_admission(files, check_host=False)
    arguments = expected_spec['arguments']
    require(arguments == mathematical_registration(arguments['panel'], expected_spec, files, check_host=False),
            'F5 mathematical registration differs from independent reconstruction')
    root = execution/'entry-output'
    require(all(read(root/(key+'.json')) == value for key, value in arguments.items()),
            'executed F5 invocation records differ')
    host = read(root/'host.json')
    require(host['system'] == 'Darwin' and host['machine'] == 'arm64'
            and type(host['cpu_count']) is int and host['cpu_count'] > 0
            and type(host['os_release']) is str and host['os_release']
            and host['cpu_model'] is None
            and host['scope'] == 'uncontrolled local development correctness; no performance claim',
            'F5 execution host differs from admitted development class')
    preflight = audit_meter(execution, root, 'build_identity', asset_role='bin/worker',
                            arguments=['--build-identity'], seconds=10)
    require(preflight['returncode'] == 0 and not preflight['timed_out']
            and read(root/'build_identity.stdout') == native['build_identity'],
            'F5 worker build preflight failed or substituted native identity')
    process = audit_meter(execution, root, 'pipeline', asset_role='bin/worker', arguments=[],
        seconds=arguments['panel']['resources']['total_wall_limit_seconds']-30, stdin_argument='job')
    common = dict(schema_version=3,
        independent_negative_adapter='exact-group-three-sum-n17-v1; postexecution auditor extension', candidate_id=arguments['seal']['candidate_id'],
        workload_id=arguments['seal']['workload_id'], run_id=arguments['seal']['run_id'],
        source_audit=source_audit, native_pipeline=publication_process(process),
        raw_native_process_sha256=digest(root/'pipeline.metrics.json'),
        native_resource_serialization='roundtrip decimal strings; raw diagnostics retained, not calibrated cost',
        complete_ic_admitted=False,
        native_binding_complete=True, online_speedup=None, promotion_eligible=False,
        headline_online_admissible=False, qualification=None)
    if process['timed_out']:
        require(read(root/'summary.json') == dict(status='NATIVE_TIMEOUT',
            candidate_id=arguments['seal']['candidate_id'], workload_id=arguments['seal']['workload_id'],
            run_id=arguments['seal']['run_id'], complete_ic_admitted=False,
            final_rank=None, verified_target_count=None, online_wall_ns=None, online_speedup=None),
            'F5 timeout summary differs')
        return dict(common, status='AUDITED_NATIVE_TIMEOUT', verified_target_count=None,
                    final_rank=None, online_wall_ns=None, natural_query_audit=None)
    report = read(root/'pipeline.stdout')
    require(report['fixture'] == arguments['fixture'] and report['mode'] == 'ic'
            and ((process['returncode'] == 0 and report['status'] == 'complete')
                 or (process['returncode'] == 2 and report['status'] == 'incomplete')),
            'F5 report substituted fixture, mode or terminal outcome')
    require(read(root/'summary.json') == dict(status=report['status'].upper(),
        run_id=arguments['seal']['run_id'], independent_mathematical_audit='pending',
        complete_ic_admitted=False, online_speedup=None), 'F5 producer summary differs')
    admitted = admit(report, arguments['fixture'], arguments['job'], build, rust_source,
        executable=execution/'asset-files/bin/worker', process_wall_ns=process['native_wall_ns'],
        resources=arguments['panel']['resources'], number=arguments['panel']['run_number'])
    require(method_record(arguments['job'], report, admitted['stages'], build) == arguments['native_method'],
            'F5 native observed method differs from declared pipeline')
    require(admitted['workload'] == arguments['workload'], 'F5 observed workload identity differs')
    # This elapsed audit diagnostic changes on replay and is not a measured
    # algorithm phase. Keep it out of the deterministic mathematical receipt.
    admitted['run'].pop('independent_audit_wall_ns')
    # Preserve the native-only admission as a subsidiary diagnostic. The
    # primary row identifies the full preregistered controller/native method.
    run = dict(admitted['run'], candidate_id=arguments['seal']['candidate_id'],
               workload_id=arguments['seal']['workload_id'], run_id=arguments['seal']['run_id'])
    complete = report['status'] == 'complete'
    verified = admitted['stages']['certificate']
    require(not complete or verified is not None, 'complete F5 pipeline lacks a verified IC certificate')
    natural = natural_queries(report, curve)
    return dict(common, status='AUDITED_COMPLETE_IC' if complete else 'AUDITED_INCOMPLETE_IC',
        complete_ic_admitted=complete, verified_target_count=1 if complete else 0,
        run=run, native_only_admission=admitted, natural_query_audit=natural,
        online_wall_ns=admitted['run']['online_wall_ns'], online_phases_ns=admitted['run']['online_phases_ns'],
        final_rank=admitted['stages']['matrix']['rank'],
        primary_timing='native first target-dependent work through native general-group scalar replay',
        external_scalar_and_pipeline_audit='postexecution Python check, excluded from native online interval')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--execution', type=Path, required=True)
    parser.add_argument('--expected-spec', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    write_immutable(args.out, audit(args.execution, read(args.expected_spec)))
