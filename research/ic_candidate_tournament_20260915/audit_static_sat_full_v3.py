"""Independent replay of the v3 source-bound SAT development pipeline.

No native solver or measured producer runs during this audit. Independent
field/group arithmetic, exact three-sum enumeration and modular rank replay
validate every ordinary and target-dependent attempt, including failures.
"""
import argparse
from collections import Counter
import json
from pathlib import Path

from audit_static_sat_full import (contents, digest, orbit_columns, raw_json, relation_row,
                                   verify_query as verify_old_query)
from audit_static_cms_s4_natural import wilson
from generic_query_law import descent_coefficients, probe_scalar
from oracle import rank, require
from run_generic_exact_yield_audit import pair_index, exact_three_sum
from sat_runtime_execution_v3 import audit_execution
from static_sat_assets_v3 import check_extracted_assets
from static_sat_inputs_v3 import native_admission
from static_sat_native_v3 import audit_meter
from static_sat_registration_v3 import mathematical_registration
from tournament import read, write


def verify_query(execution,panel,root,prefix,row,point,curve,base,pairs):
    directory=root/prefix
    if point is None:
        require(row['status']=='IDENTITY_QUERY' and row['public_point'] is None
                and row['cms'] is None and row['exporter'] is None
                and row['point_witness'] is None and row['verification_wall_ns']==0,
                'identity query claims native work or a point witness')
        return exact_three_sum(curve,base,pairs,None) is not None,None,'IDENTITY_QUERY'
    # Reconstruct exact native argv independently of the producer's builder.
    exporter_args=['17','6','standard',str(panel['export_nonce']),
                   str(panel['cms_conflict_budget']),'instance','1','0',
                   '--target-x',str(point[0]),'--target-y',str(point[1]),
                   '--blind-instance-id',f"control-{row['trial']:02d}",'--export-only']
    export=audit_meter(execution,directory,'export',asset_role='bin/exporter',
                       arguments=exporter_args,seconds=panel['export_timeout_seconds'])
    require(export==row['exporter'],'SAT query exporter receipt changed')
    native_wall=export['native_wall_ns']
    if row['status']=='EXPORT_FAILURE':
        require(export['timed_out'] or export['returncode']!=0,
                'SAT export failure has a successful native receipt')
    if row['cms'] is not None:
        cms=audit_meter(execution,directory,'cms',asset_role='bin/cms',
                        arguments=['--verb','1','--threads','1','--random','1',
                                   '--maxsol','1','--maxconfl',
                                   str(panel['cms_conflict_budget']),'instance/instance.xor.cnf'],
                        seconds=panel['cms_timeout_seconds'])
        require(cms==row['cms'],'SAT query solver receipt changed')
        native_wall+=cms['native_wall_ns']
    require(row['query_wall_ns']>=native_wall,
            'SAT query interval omits sequential native process time')
    if row['status'] not in ('EXPORT_FAILURE','INVALID_EXPORT'):
        manifest=raw_json(root,prefix+'instance/manifest.json')
        require(manifest['kind']=='binary_koblitz_pdp_cross_solver_instance'
                and manifest['target_mode']=='explicit_affine'
                and manifest['ell']==6 and manifest['m']==3
                and manifest['irreducible_low_terms']==[0,3]
                and manifest['factor_base_basis_bitmasks']==['1','2','4','8','16','32']
                and manifest['seed']==panel['export_nonce']
                and manifest['blind_instance_id']==f"control-{row['trial']:02d}"
                and manifest['native_sat']['status']=='not_run_in_export_process'
                and manifest['direct_meet_in_the_middle']['status']=='not_run_in_export_process',
                'SAT export encoding or instance settings changed')
    return verify_old_query(root,prefix,row,point,curve,base,pairs)


def audit(execution,expected_spec):
    execution=Path(execution).resolve()
    source_audit=audit_execution(execution,expected_spec)
    require(source_audit['entrypoint_succeeded'],
            'SAT failed/partial execution needs a failure receipt, never complete admission')
    assets=check_extracted_assets(execution/'asset-files',expected_spec['asset_manifest'])
    arguments=expected_spec['arguments']
    panel=arguments['panel']
    require(mathematical_registration(panel,expected_spec,assets)==arguments,
            'SAT mathematical registration does not reconstruct')
    fixture,report,curve,base,native=native_admission(assets)
    target=curve.decode(panel['target_input']['point'])
    columns,mapping=orbit_columns(curve,base)
    root=execution/'entry-output'
    source,method,candidate,workload,seal=(arguments[name] for name in (
        'source','method','candidate','workload','seal'))
    for name in ('panel','source','method','candidate','workload','seal'):
        require(raw_json(root,name+'.json')==arguments[name],
                'SAT retained mathematical record differs: '+name)
    preflight=audit_meter(execution,root,'cms_preflight',asset_role='bin/cms',
                          arguments=['--version'],seconds=10)
    require(preflight['returncode']==0 and not preflight['timed_out']
            and b'CryptoMiniSat version 5.14.7' in contents(root,'cms_preflight.stdout'),
            'SAT native version preflight failed')
    result=raw_json(root,'summary.json')
    require(result['status'] in {'COMPLETE','INCOMPLETE_RELATION_RANK','INCOMPLETE_TARGET'}
            and result['candidate_id']==candidate['candidate_id']
            and result['workload_id']==workload['workload_id']
            and result['run_id']==seal['run_id']
            and result['panel_sha256']==seal['panel_sha256']
            and result['source_binding']==source
            and result['exporter_build']==json.loads(assets['exporter/build-record.json'])
            and result['cms_preflight']==preflight
            and result['target_input']==panel['target_input']
            and result['online_speedup'] is None,
            'SAT result identity, source or claim boundary differs')
    collection = result['collection']
    require(0 < len(collection) <= panel['max_relation_queries']
            and [json.loads(line) for line in contents(
                root, 'collection.progress.jsonl').decode().splitlines()]
                == collection,
            'ordinary query chronology differs from retained progress')
    pairs = pair_index(curve, base)
    statuses = Counter()
    exact_feasible = 0
    matrix_rows = []
    trajectory = []
    duplicates = dependencies = 0
    seen = set()
    for trial, row in enumerate(collection):
        scalar = probe_scalar(panel['relation_query_seed'], trial, curve.r)
        point = curve.mul(curve.g, scalar)
        require(row['trial'] == trial and row['probe_scalar'] == scalar,
                'ordinary scalar differs from frozen query law')
        feasible, indices, status = verify_query(
            execution, panel, root, f'collection/trial-{trial:02d}/', row,
            point, curve, base, pairs)
        exact_feasible += feasible
        statuses[status] += 1
        before = rank([item[0] for item in matrix_rows], len(columns), curve.r)
        if indices is not None:
            coeff, rhs = relation_row(
                curve, base, mapping, columns, scalar, indices)
            key = scalar, tuple(sorted(indices))
            if key in seen:
                expected = 'duplicate'
                duplicates += 1
            else:
                seen.add(key)
                matrix_rows.append((coeff, rhs, scalar, indices))
                after = rank([item[0] for item in matrix_rows],
                             len(columns), curve.r)
                expected = 'novel_rank' if after > before else 'dependent'
                dependencies += after == before
            require(row['relation_status'] == expected,
                    'relation novelty differs from independent reconstruction')
        else:
            require(row['relation_status'] == 'no_verified_relation',
                    'failed query contributed a relation')
        actual_rank = rank([item[0] for item in matrix_rows],
                           len(columns), curve.r)
        require(row['rank_after'] == actual_rank
                and row['attempt_wall_ns'] >= row['query_wall_ns'] >= 0
                and row['query_wall_ns'] >= row['verification_wall_ns'] >= 0
                and row['pdp_wall_ns']
                    == row['query_wall_ns']-row['verification_wall_ns']
                and row['attempt_wall_ns'] >= row['matrix_update_wall_ns'] >= 0,
                'rank trajectory or exclusive ordinary-query clock differs')
        trajectory.append(actual_rank)
    require(result['rank_trajectory'] == trajectory,
            'stored rank trajectory changed')
    snapshot = raw_json(root, 'relation-matrix.json')
    rows = [dict(entries=[[j, str(value)] for j, value in enumerate(coeff)
                          if value], rhs=str(rhs), scalar=scalar,
                 indices=list(indices))
            for coeff, rhs, scalar, indices in matrix_rows]
    require(snapshot == result['matrix']
            and snapshot['modulus'] == str(curve.r)
            and snapshot['column_points'] == [list(point) for point in columns]
            and snapshot['columns'] == len(columns)
            and snapshot['rank'] == trajectory[-1]
            and snapshot['accepted_rows'] == len(matrix_rows)
            and snapshot['duplicate_relations'] == duplicates
            and snapshot['dependent_relations'] == dependencies
            and snapshot['rows'] == rows
            and snapshot['rows_sha256']
                == digest(json.dumps(rows, sort_keys=True,
                                     separators=(',', ':')).encode())
            and result['matrix_update_wall_ns']
                == sum(row['matrix_update_wall_ns'] for row in collection)
            and result['collection_wall_ns'] >= result['matrix_update_wall_ns']
            and result['preparation_wall_ns'] >= result['collection_wall_ns'],
            'final matrix or collection costs differ from raw rows')

    common = dict(schema_version=3, complete_source_bound=True,
                  candidate_id=candidate['candidate_id'],
                  workload_id=workload['workload_id'], run_id=seal['run_id'],
                  attempts=len(collection), statuses=dict(statuses),
                  exact_feasible=exact_feasible,
                  verified_relations=statuses['VALID_POINT_WITNESS'],
                  witness_rate_wilson95=wilson(
                      statuses['VALID_POINT_WITNESS'], len(collection)),
                  final_rank=snapshot['rank'], columns=len(columns),
                  same_point_rho_audited=False, scientific_control_only=True,
                  promotion_eligible=False, fresh_target_qualified=False,
                  headline_online_admissible=False, online_speedup=None)
    if snapshot['rank'] < len(columns):
        require(len(collection) == panel['max_relation_queries']
                and result['status'] == 'INCOMPLETE_RELATION_RANK'
                and result['target_attempts'] == []
                and result['online_wall_ns'] is None
                and result['recovered_scalar'] is None,
                'rank-deficient run overclaims one-target recovery')
        return dict(common, status='AUDITED_INCOMPLETE_RANK',
                    solved_targets=0, online_wall_ns=None,
                    online_endpoint_admissible=False)
    require(result['relation_la_wall_ns'] is not None
            and result['relation_la_wall_ns'] >= 0,
            'full-rank matrix lacks final LA cost')
    logs = [int(item['log']) for item in result['column_logs']]
    retained_logs = raw_json(root, 'column-logs.json')
    require(len(logs) == len(columns)
            and all(0 <= log < curve.r for log in logs)
            and [item['point'] for item in result['column_logs']]
                == [list(point) for point in columns]
            and all(curve.mul(curve.g, log) == point
                    for log, point in zip(logs, columns))
            and all(sum(a*b for a, b in zip(coeff, logs)) % curve.r == rhs
                    for coeff, rhs, _, _ in matrix_rows)
            and retained_logs == dict(columns=result['column_logs'],
                                      independently_verified=True),
            'solved logs fail group or matrix replay')
    target_rows = result['target_attempts']
    require(0 < len(target_rows) <= panel['max_descent_queries']
            and [json.loads(line) for line in contents(
                root, 'descent.progress.jsonl').decode().splitlines()]
                == target_rows,
            'target-dependent progress differs from summary')
    stream = descent_coefficients(panel['descent_query_seed'], curve.r,
                                  walked=False)
    target_statuses = Counter()
    recovered = None
    for index, row in enumerate(target_rows):
        a, b = next(stream)
        point = curve.add(curve.mul(curve.g, a), curve.mul(target, b))
        require(row['target_query_index'] == index
                and row['trial'] == panel['max_relation_queries']+index
                and row['probe_scalar'] == a
                and (row['a'], row['b']) == (a, b),
                'target aG+bQ query differs from frozen law')
        _, indices, status = verify_query(
            execution, panel, root, f"descent/trial-{row['trial']:02d}/", row,
            point, curve, base, pairs)
        target_statuses[status] += 1
        require(row['query_wall_ns'] >= row['verification_wall_ns'] >= 0
                and row['pdp_wall_ns']
                    == row['query_wall_ns']-row['verification_wall_ns'],
                'target query clock or phase split changed')
        if indices is not None:
            projected = 0
            for base_index in indices:
                if mapping[base_index] is not None:
                    column, coefficient = mapping[base_index]
                    projected = (projected + coefficient*logs[column]) % curve.r
            recovered = (projected-curve.h*a)*pow(
                curve.h*b % curve.r, -1, curve.r) % curve.r
            require(row['candidate_scalar'] == str(recovered)
                    and row['scalar_replay_verified'] is True
                    and curve.mul(curve.g, recovered) == target
                    and index == len(target_rows)-1,
                    'target descent/scalar certificate fails replay')
    complete = recovered is not None
    phases = result['online_phases_ns']
    require(set(phases) == {'target_query', 'target_pdp',
                            'target_relation_check', 'target_descent',
                            'target_recovery_check'}
            and all(type(value) is int and value >= 0
                    for value in phases.values())
            and phases['target_pdp']
                == sum(row['pdp_wall_ns'] for row in target_rows)
            and phases['target_relation_check'] >= sum(
                row['verification_wall_ns'] for row in target_rows)
            and sum(phases.values()) == result['online_attempt_wall_ns']
            and type(result['online_start_monotonic_ns']) is int
            and type(result['online_stop_monotonic_ns']) is int
            and result['online_stop_monotonic_ns']
                - result['online_start_monotonic_ns'] == result['online_attempt_wall_ns']
            and result['online_bookkeeping_assigned_to_target_query_ns'] >= 0
            and result['online_stop_event'] == (
                'independent-scalar-replay' if complete
                else 'frozen-target-attempt-cap')
            and result['status'] == (
                'COMPLETE' if complete else 'INCOMPLETE_TARGET')
            and result['scalar_verified'] is complete
            and result['recovered_scalar']
                == (str(recovered) if complete else None),
            'one-target endpoint, exclusive phases or scalar status differs')
    require(result['online_wall_ns'] == (result['online_attempt_wall_ns'] if complete else None),
            'failed target claims verified online wall time')
    require(complete or len(target_rows) == panel['max_descent_queries'],
            'incomplete target stopped before frozen attempt cap')
    return dict(common,
                status=('AUDITED_COMPLETE_SINGLE_ARM' if complete
                        else 'AUDITED_INCOMPLETE_TARGET'),
                target_attempts=len(target_rows),
                target_statuses=dict(target_statuses),
                solved_targets=int(complete),
                recovered_scalar=None if recovered is None else str(recovered),
                online_wall_ns=result['online_wall_ns'],
                online_attempt_wall_ns=result['online_attempt_wall_ns'],
                online_phases_ns=phases,
                online_endpoint_admissible=complete)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--execution',type=Path,required=True)
    parser.add_argument('--expected-spec',type=Path,required=True)
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args()
    require(not args.out.exists(),'SAT audit output already exists')
    write(args.out,audit(args.execution,read(args.expected_spec)),exclusive=True)


if __name__=='__main__':
    main()
