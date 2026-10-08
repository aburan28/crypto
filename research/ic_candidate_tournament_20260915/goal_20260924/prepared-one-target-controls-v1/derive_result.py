#!/usr/bin/env python3
"""Close these sealed, disclosed toy controls; never launch a native solver."""
import argparse
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))
from identity import curve_record, sha256, write_immutable  # noqa: E402
from oracle import Curve, require  # noqa: E402
from sat_runtime_execution_v3 import digest, read  # noqa: E402
from diagnose_controls import SEALS  # noqa: E402


def derive(root, restored, diagnosis, restored_diagnosis, out):
    root, restored, diagnosis, restored_diagnosis, out = map(
        Path, (root, restored, diagnosis, restored_diagnosis, out))
    require(not out.exists(), 'result output already exists')
    analysis = read(diagnosis)
    require(diagnosis.read_bytes() == restored_diagnosis.read_bytes()
            and analysis['native_solvers_executed'] == 0
            and analysis['fresh_targets_generated'] == 0,
            'restored independent diagnosis differs')
    out.mkdir(parents=True)
    files = []

    def retain(source, role):
        dest = out/role
        dest.parent.mkdir(parents=True, exist_ok=True)
        with dest.open('xb') as stream:
            stream.write(source.read_bytes())
        require(digest(source) == digest(dest), 'compact copy changed')
        files.append(dict(role=role, bytes=dest.stat().st_size, sha256=digest(dest)))

    retain(diagnosis, 'diagnosis.json')
    target_points, query_points, rows, replay = set(), set(), {}, {}
    fixtures = []
    for family, seal in SEALS.items():
        execution = root/(family+'-execution')
        spec = read(execution/'execution.json')
        require(sha256(spec) == seal and spec == read(restored/(family+'-execution')/'execution.json'),
                'control is outside the external registration')
        fixture = spec['arguments']['preparation']['certificate']['inputs']['fixture']
        fixtures.append(fixture)
        curve = Curve(fixture)
        require(curve.n == 17, 'diagnosis is restricted to the disclosed toy curve')
        target = curve.decode(spec['arguments']['panel']['target_input']['point'])
        require(target == (52411, 72106), 'disclosed public point changed')
        target_points.add(target)
        arm = analysis['arms'][family]
        require(arm['execution_sha256'] == seal and arm['attempt_count'] == 8
                and arm['false_negative_claims'] == 0, 'attempt diagnosis changed')
        query_points.update(curve.decode(attempt['point']) for attempt in arm['attempts'])
        process = read(execution/'process.json')
        require(process['exit_code'] == 0 and not process['timed_out'], 'controller was not terminal')
        registration = root/(family+'-registration')
        claim = read(registration/'execution-claim.json')
        require(sha256(claim) == process['execution_claim_sha256'], 'one-shot claim changed')
        original = root/(family+'-audit')
        relocated = Path(str(restored)+'-'+family+'-audit')
        transport, transported = read(original/'transport.json'), read(relocated/'transport.json')
        require(all(transport[key] == transported[key] for key in (
            'status', 'exit_code', 'timed_out', 'binding', 'execution_sha256', 'native_solvers_executed')),
            'relocation changed frozen transport outcome')
        require(transport['native_solvers_executed'] == 0 and not transport['timed_out'],
                'audit attempted native work or timed out')
        summary = read(execution/'entry-output/summary.json')
        rows[family] = dict(candidate_id=summary['candidate_id'], workload_id=summary['workload_id'],
                            run_id=summary['run_id'], execution_sha256=seal,
                            controller_exit_code=0, controller_timed_out=False,
                            registration_consumed=True, native_retry_allowed=False,
                            attempt_count=8, exact_geometric_feasible_queries=arm['exact_feasible_count'],
                            false_negative_claims=0, original_transport_status=transport['status'],
                            scalar_verified=False, recovered_scalar=None, verified_online_wall_ns=None,
                            cold_wall_ns=None, rho_online_wall_ns=None, online_speedup=None,
                            fresh_paired_qualification=False, headline_online_admissible=False,
                            promotion_eligible=False)
        replay[family] = dict(original_transport_sha256=digest(original/'transport.json'),
                              restored_transport_sha256=digest(relocated/'transport.json'),
                              original_before_sha256=digest(original/'before.json'),
                              restored_before_sha256=digest(relocated/'before.json'),
                              original_transport_status=transport['status'],
                              restored_transport_status=transported['status'],
                              native_solvers_executed=0)
        require(replay[family]['original_before_sha256'] == replay[family]['restored_before_sha256'],
                'relocation changed initial loaded-module gate')
        if family == 'f5':
            require(transport['status'] == 'REJECTED_PREPARED_TRANSPORT'
                    and transport['exit_code'] == 1
                    and 'oracle.InvalidEvidence: missing query schema' in (original/'stderr.txt').read_text()
                    and 'oracle.InvalidEvidence: missing query schema' in (relocated/'stderr.txt').read_text(),
                    'original F5 header rejection was not reproduced')
            report = read(execution/'entry-output/pipeline.stdout')
            require('query_schema_version' not in report and arm['exact_feasible_count'] == 0
                    and arm['view_replay']['status'] == 'REPLAYED_POSTEXECUTION_VIEW_ONLY'
                    and not arm['view_replay']['original_transport_admitted'],
                    'postexecution view cannot rewrite original admission')
            rows[family].update(status='INCOMPLETE_AND_FROZEN_AUDIT_REJECTED',
                                native_status_mix={'proved_unsat': 8},
                                original_mathematical_admission=False,
                                original_failure='missing query schema',
                                postexecution_view_status='REPLAYED_POSTEXECUTION_VIEW_ONLY')
            replay[family]['stderr_difference'] = 'retained traceback locations change on relocation'
            retain(execution/'entry-output/pipeline.stdout', family+'/pipeline.stdout')
        else:
            require(transport['status'] == 'PASS_FROZEN_PREPARED_TRANSPORT'
                    and original.joinpath('admission.json').read_bytes() == relocated.joinpath('admission.json').read_bytes()
                    and arm['exact_feasible_count'] == 1, 'SAT incomplete admission changed')
            admission = read(original/'admission.json')
            require(admission['status'] == 'ADMITTED_INCOMPLETE_PREPARED_SAT_CONTROL'
                    and not admission['scalar_verified'] and admission['online_wall_ns'] is None,
                    'incomplete SAT evidence cannot become a solved target')
            rows[family].update(status='ADMITTED_INCOMPLETE_PREPARED_SAT_CONTROL',
                                native_status_mix={'CONFLICT_BUDGET_INCONCLUSIVE': 8},
                                original_mathematical_admission=True)
            replay[family].update(admission_byte_identical=True,
                                  admission_sha256=digest(original/'admission.json'))
            retain(original/'admission.json', family+'/admission.json')
            retain(original/'after.json', family+'/audit-after.json')
        for source, role in (
            (execution/'process.json', family+'/process.json'),
            (execution/'entry-output/summary.json', family+'/summary.json'),
            (registration/'execution-claim.json', family+'/execution-claim.json'),
            (original/'transport.json', family+'/transport.json'),
            (original/'before.json', family+'/audit-before.json'),
            (original/'stderr.txt', family+'/audit-stderr.txt'),
            (relocated/'transport.json', 'restored/'+family+'/transport.json'),
            (relocated/'before.json', 'restored/'+family+'/audit-before.json'),
            (relocated/'stderr.txt', 'restored/'+family+'/audit-stderr.txt')):
            retain(source, role)
    require(fixtures[0] == fixtures[1], 'control curves differ')
    curve = Curve(fixtures[0])
    exclusions = set()
    for point in target_points | query_points:
        current = point
        for _ in range(curve.n):
            exclusions.update((current, curve.neg(current)))
            current = curve.frob(current)
        require(current == point, 'Frobenius orbit did not close')
    # No sampling: only the already exposed point and the recorded query points.
    require(None not in exclusions, 'unexpected identity exposure')
    exposure = dict(schema_version=1, status='PARTIAL_CURRENT_CONTROL_EXPOSURE_CENSUS',
                    curve=curve_record(fixtures[0]), execution_sha256=SEALS,
                    diagnosis_sha256=digest(diagnosis), generator_script_sha256=digest(__file__),
                    public_targets=[list(p) for p in sorted(target_points)],
                    recorded_query_points=[list(p) for p in sorted(query_points)],
                    excluded_points=[list(p) for p in sorted(exclusions)],
                    excluded_point_count=len(exclusions),
                    closure='negation and all 17 Frobenius powers of target and recorded queries',
                    fresh_targets_generated=0,
                    complete_historical_exposure_census=False,
                    scope='merge with all historical/current and preparation-only exclusions before fresh sampling')
    write_immutable(out/'new-exposures.json', exposure)
    write_immutable(out/'restoration-replay.json', dict(
        schema_version=1, status='REPRODUCED_ORIGINAL_OUTCOMES_AFTER_RELOCATION', arms=replay,
        diagnosis_byte_identical=True, diagnosis_sha256=digest(diagnosis), native_solvers_executed=0,
        scope='local archive relocation using the same bound macOS interpreter; no cross-host claim'))
    result = dict(schema_version=1, status='CLOSED_INCOMPLETE_ONE_SHOT_CONTROLS', arms=rows,
                  curve_id=exposure['curve']['curve']['curve_id'], target_point=[52411, 72106],
                  generator_script_sha256=digest(__file__), compact_raw_inventory=files,
                  native_prepared_controller_invocations=2, further_native_invocations_by_analysis=0,
                  fresh_targets_generated=0, complete_ic_admitted=False,
                  promotion_eligible=False, online_speedup=None,
                  full_goal_complete=False,
                  uncertainty='no population rate or comparative timing is estimated from these fixed eight-query controls')
    write_immutable(out/'result.json', result)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for flag in ('root', 'restored', 'diagnosis', 'restored-diagnosis', 'out'):
        parser.add_argument('--'+flag, type=Path, required=True)
    args = parser.parse_args()
    result = derive(args.root, args.restored, args.diagnosis, args.restored_diagnosis, args.out)
    print(json.dumps({family: row['status'] for family, row in result['arms'].items()}, sort_keys=True))
