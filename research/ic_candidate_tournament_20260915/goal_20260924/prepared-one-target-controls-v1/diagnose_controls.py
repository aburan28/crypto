#!/usr/bin/env python3
"""Postexecution diagnosis of these two disclosed toy controls; no native call.

The original frozen F5 audit rejection is preserved. A labelled in-memory
header reconstruction locates further interface failures; it never changes
the raw report, reclassifies that transport, or qualifies a solved target.
"""
import argparse
import copy
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))
from identity import sha256, write_immutable  # noqa: E402
from oracle import Curve, require  # noqa: E402
from prepared_target_v1 import audit_native_target  # noqa: E402
from run_generic_exact_yield_audit import exact_three_sum, pair_index  # noqa: E402
from sat_runtime_execution_v3 import audit_execution, digest, read  # noqa: E402

SEALS = {
    'f5': '8c5afe8d3f355010b4b293405f69e76a2362db25ae163490a1685955883b0737',
    'sat': '43539f7d440289dae1ad4867ba1ca951bd4664bf01147b9070eaa68f3bad07cb',
}


def diagnose(root, out):
    root, out = Path(root).resolve(), Path(out).resolve()
    require(not out.exists() and not out.is_relative_to(root), 'diagnosis requires a new separate output')
    result = dict(schema_version=1, status='POSTEXECUTION_DISCLOSED_CONTROL_DIAGNOSIS',
                  script_sha256=digest(__file__), native_solvers_executed=0,
                  fresh_targets_generated=0, complete_ic_admitted=False,
                  promotion_eligible=False, online_speedup=None, arms={})
    for family in ('f5', 'sat'):
        execution = root/(family+'-execution')
        spec = read(execution/'execution.json')
        require(sha256(spec) == SEALS[family], 'diagnosis is limited to the original frozen controls')
        source = audit_execution(execution, spec)
        require(source['entrypoint_succeeded'], 'partial controller cannot supply this diagnosis')
        document = spec['arguments']['preparation']
        inputs = document['certificate']['inputs']
        curve = Curve(inputs['fixture'])
        base = tuple(curve.decode(p) for p in inputs['base'])
        require(curve.n == 17 and len(base) == 63,
                'diagnosis is restricted to the sealed toy geometry')
        target = curve.decode(spec['arguments']['panel']['target_input']['point'])
        pairs = pair_index(curve, base)
        if family == 'f5':
            report_path = execution/'entry-output/pipeline.stdout'
            report = read(report_path)
            attempts = report['solutions'][0]['attempts']
        else:
            report_path = execution/'entry-output/summary.json'
            report = read(report_path)
            attempts = report['target_attempts']
        require(len(attempts) == 8, 'diagnosis lost a frozen failed attempt')
        rows = []
        for trial, attempt in enumerate(attempts):
            require(attempt['trial'] == trial, 'diagnosis chronology changed')
            point = curve.add(curve.mul(curve.g, attempt['a']), curve.mul(target, attempt['b']))
            witness = exact_three_sum(curve, base, pairs, point)
            outcome = attempt['pdp']['outcome'] if family == 'f5' else attempt['status']
            rows.append(dict(trial=trial, point=None if point is None else list(point),
                             native_outcome=outcome, exact_geometric_feasible=witness is not None,
                             independently_readded_witness=witness,
                             false_negative_claim=outcome in ('proved_unsat', 'SOURCE_UNSAT')
                                and witness is not None))
        arm = dict(execution_sha256=SEALS[family], raw_report_sha256=digest(report_path),
                   source_binding_audit=source, attempts=rows, attempt_count=len(rows),
                   exact_feasible_count=sum(row['exact_geometric_feasible'] for row in rows),
                   false_negative_claims=sum(row['false_negative_claim'] for row in rows),
                   scope='exact feasibility of these recorded target queries; no natural ordinary-yield estimate')
        if family == 'f5':
            require('query_schema_version' not in report, 'original missing-header finding changed')
            rejected = read(root/'f5-audit/transport.json')
            require(rejected['status'] == 'REJECTED_PREPARED_TRANSPORT'
                    and rejected['exit_code'] == 1 and not rejected['timed_out'],
                    'original frozen audit rejection must remain retained')
            arm['original_frozen_transport'] = rejected
            arm['original_failure_stderr_sha256'] = digest(root/'f5-audit/stderr.txt')
            view = copy.deepcopy(report)
            view['query_schema_version'] = 1
            arm['postexecution_view_change'] = {'query_schema_version': 1}
            try:
                replay = audit_native_target(view, spec['arguments']['job'], document,
                                             spec['arguments']['preparation_certificate_sha256'])
            except Exception as exception:
                arm['view_replay'] = dict(status='REJECTED_POSTEXECUTION_VIEW',
                                         error=dict(type=type(exception).__name__, message=str(exception)))
            else:
                arm['view_replay'] = dict(status='REPLAYED_POSTEXECUTION_VIEW_ONLY', mathematical=replay,
                                         original_transport_admitted=False, complete_ic_admitted=False)
        else:
            admitted = read(root/'sat-audit/admission.json')
            require(admitted['status'] == 'ADMITTED_INCOMPLETE_PREPARED_SAT_CONTROL'
                    and admitted['scalar_verified'] is False and admitted['online_wall_ns'] is None,
                    'incomplete SAT control cannot become a complete result')
            arm['original_admission_sha256'] = sha256(admitted)
        result['arms'][family] = arm
    out.mkdir(parents=True)
    write_immutable(out/'diagnosis.json', result)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    result = diagnose(args.root, args.out)
    print(json.dumps({name: {'attempts': arm['attempt_count'], 'feasible': arm['exact_feasible_count'],
                             'false_negative_claims': arm['false_negative_claims'],
                             'view_replay': arm.get('view_replay', {}).get('status')}
                      for name, arm in result['arms'].items()}, sort_keys=True))
