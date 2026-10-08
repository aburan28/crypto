"""One-shot source-bound SAT collection, prime LA, and one-target recovery.

The n17 adapter admits a development correctness control. It supplies no
paired/calibrated performance evidence and does not reopen a historical run.
"""
import argparse
import json
import os
from pathlib import Path
import platform
import time

from generic_query_law import descent_coefficients, probe_scalar
from oracle import require
from identity import sha256
from sat_runtime_execution_v3 import execute, read
from static_sat_assets_v3 import check_extracted_assets
from static_sat_inputs_v3 import native_admission
from static_sat_matrix import RelationMatrix
from static_sat_native_v3 import meter
from static_sat_query_v3 import one_query
from static_sat_registration_v3 import mathematical_registration
from tournament import write


def save_progress(out,name,row):
    with (out/name).open('a') as stream:
        stream.write(json.dumps(row,sort_keys=True)+'\n')


def run(arguments,out):
    out=Path(out).resolve()
    try:
        return run_registered(arguments,out)
    except Exception as error:
        if not (out/'summary.json').exists():
            seal=arguments['seal']
            write(out/'summary.json',dict(schema_version=3,status='ERROR',
                  error_type=type(error).__name__,error=str(error),
                  candidate_id=seal['candidate_id'],workload_id=seal['workload_id'],
                  run_id=seal['run_id'],recovered_scalar=None,scalar_verified=False,
                  online_wall_ns=None,online_phases_ns=None,online_speedup=None,
                  progress_files=['collection.progress.jsonl','descent.progress.jsonl']),
                  exclusive=True)
        raise


def run_registered(arguments,out):
    preparation_start=time.monotonic_ns()
    execution=out.parent
    spec=read(execution/'execution.json')
    require(out==execution/'entry-output' and not (out/'summary.json').exists(),
            'SAT pipeline output reused or outside execution')
    assets=check_extracted_assets(execution/'asset-files',spec['asset_manifest'])
    panel=arguments['panel']
    require(mathematical_registration(panel,spec,assets)==arguments
            and spec['arguments']==arguments, 'SAT mathematical registration changed')
    fixture,report,curve,base,native=native_admission(assets)
    require(dict(system=platform.system(),machine=platform.machine())==native['platform'],
            'SAT native binaries require their validated physical platform')
    target=curve.decode(panel['target_input']['point'])
    matrix=RelationMatrix(curve,base)
    seal=arguments['seal']
    binding=arguments['source']
    build_record=json.loads(assets['exporter/build-record.json'])
    for name in ('panel','source','method','candidate','workload','seal'):
        write(out/(name+'.json'),arguments[name],exclusive=True)
    write(out/'host.json',dict(system=platform.system(),machine=platform.machine(),
          cpu_count=os.cpu_count(),scope='uncontrolled local development correctness; no performance claim'),
          exclusive=True)
    version=meter(execution,'bin/cms',['--version'],out,'cms_preflight',10)
    require(version['returncode']==0 and not version['timed_out']
            and 'CryptoMiniSat version 5.14.7' in (out/'cms_preflight.stdout').read_text(),
            'bound static SAT solver failed version preflight')
    (out/'collection').mkdir()
    (out/'descent').mkdir()
    collection=[]
    rank_trajectory=[]
    logs=None
    relation_la_wall_ns=None
    collection_start=time.monotonic_ns()
    matrix_update_total_ns=0
    for trial in range(panel['max_relation_queries']):
        attempt_start=time.monotonic_ns()
        scalar=probe_scalar(panel['relation_query_seed'],trial,curve.r)
        point=curve.mul(curve.g,scalar)
        item=dict(trial=trial,probe_scalar=scalar,point=list(point))
        query_start=time.monotonic_ns()
        row=one_query(panel,item,execution,curve,base,out/'collection')
        query_wall_ns=time.monotonic_ns()-query_start
        require(0<=row['verification_wall_ns']<=query_wall_ns,
                'ordinary query verification interval exceeds query wall')
        row['query_wall_ns']=query_wall_ns
        row['pdp_wall_ns']=query_wall_ns-row['verification_wall_ns']
        if row['status']=='VALID_POINT_WITNESS':
            indices=row['point_witness']['point_indices']
            matrix_start=time.monotonic_ns()
            row['relation_status']=matrix.push(scalar,indices)
            row['matrix_update_wall_ns']=time.monotonic_ns()-matrix_start
            matrix_update_total_ns+=row['matrix_update_wall_ns']
        else:
            row['relation_status']='no_verified_relation'
            row['matrix_update_wall_ns']=0
        row['rank_after']=matrix.rank
        row['attempt_wall_ns']=time.monotonic_ns()-attempt_start
        collection.append(row)
        rank_trajectory.append(matrix.rank)
        save_progress(out,'collection.progress.jsonl',row)
        print(json.dumps(dict(stage='collection',trial=trial,
                              status=row['status'],rank=matrix.rank)),flush=True)
        if matrix.rank==len(matrix.columns):
            la_start=time.monotonic_ns()
            logs=matrix.solve()
            relation_la_wall_ns=time.monotonic_ns()-la_start
            break
    collection_wall_ns=time.monotonic_ns()-collection_start
    snapshot=matrix.snapshot()
    write(out/'relation-matrix.json',snapshot,exclusive=True)
    if logs is None:
        write(out/'summary.json',dict(schema_version=3,
              status='INCOMPLETE_RELATION_RANK',panel_sha256=seal['panel_sha256'],
              source_binding=binding,exporter_build=build_record,
              cms_preflight=version,collection=collection,
              rank_trajectory=rank_trajectory,matrix=snapshot,
              preparation_wall_ns=time.monotonic_ns()-preparation_start,
              collection_wall_ns=collection_wall_ns,
              matrix_update_wall_ns=matrix_update_total_ns,
              relation_la_wall_ns=relation_la_wall_ns,
              target_input=panel['target_input'],target_attempts=[],
              recovered_scalar=None,scalar_verified=False,
              online_wall_ns=None,online_phases_ns=None,
              candidate_id=seal['candidate_id'],
              workload_id=seal['workload_id'],run_id=seal['run_id'],
              online_speedup=None),exclusive=True)
        return dict(status='INCOMPLETE_RELATION_RANK',run_id=seal['run_id'])
    write(out/'column-logs.json',dict(
        columns=[dict(point=list(point),log=str(log))
                 for point,log in zip(matrix.columns,logs)],
        independently_verified=True),exclusive=True)

    # Nothing target-dependent runs before this interval begins.
    phases={name:0 for name in ('target_query','target_pdp',
                                 'target_relation_check','target_descent',
                                 'target_recovery_check')}
    attempts=[]
    recovered=None
    preparation_wall_ns=time.monotonic_ns()-preparation_start
    online_start=time.monotonic_ns()
    online_stop=None
    stream=descent_coefficients(panel['descent_query_seed'],curve.r,walked=False)
    for index in range(panel['max_descent_queries']):
        started=time.monotonic_ns()
        a,b=next(stream)
        query=curve.add(curve.mul(curve.g,a),curve.mul(target,b))
        phases['target_query']+=time.monotonic_ns()-started
        item=dict(trial=panel['max_relation_queries']+index,
                  probe_scalar=a,point=None if query is None else list(query))
        started=time.monotonic_ns()
        row=one_query(panel,item,execution,curve,base,out/'descent')
        query_wall_ns=time.monotonic_ns()-started
        require(0<=row['verification_wall_ns']<=query_wall_ns,
                'target query verification interval exceeds query wall')
        phases['target_pdp']+=query_wall_ns-row['verification_wall_ns']
        phases['target_relation_check']+=row['verification_wall_ns']
        row['query_wall_ns']=query_wall_ns
        row['pdp_wall_ns']=query_wall_ns-row['verification_wall_ns']
        row['a'],row['b'],row['target_query_index']=a,b,index
        if row['status']=='VALID_POINT_WITNESS':
            indices=row['point_witness']['point_indices']
            started=time.monotonic_ns()
            group_sum=None
            for base_index in indices:
                group_sum=curve.add(group_sum,base[base_index])
            require(group_sum==query,'target relation fails independent group re-addition')
            phases['target_relation_check']+=time.monotonic_ns()-started
            started=time.monotonic_ns()
            candidate=matrix.descent_scalar(logs,a,b,indices)
            phases['target_descent']+=time.monotonic_ns()-started
            started=time.monotonic_ns()
            verified=curve.mul(curve.g,candidate)==target
            phases['target_recovery_check']+=time.monotonic_ns()-started
            require(verified,'target descent scalar fails independent replay')
            online_stop=time.monotonic_ns()
            row['candidate_scalar']=str(candidate)
            row['scalar_replay_verified']=verified
            recovered=candidate
        if online_stop is None and index+1==panel['max_descent_queries']:
            online_stop=time.monotonic_ns()
        attempts.append(row)
        save_progress(out,'descent.progress.jsonl',row)
        print(json.dumps(dict(stage='target',trial=index,status=row['status'],
                              verified=recovered is not None)),flush=True)
        if recovered is not None:
            break
    require(online_stop is not None,'one-target online stop event missing')
    online_wall=online_stop-online_start
    residual=online_wall-sum(phases.values())
    require(residual>=0,'exclusive target phase clocks overlap')
    phases['target_query']+=residual
    require(sum(phases.values())==online_wall,'online phase ledger does not close')
    write(out/'summary.json',dict(schema_version=3,
          status='COMPLETE' if recovered is not None else 'INCOMPLETE_TARGET',
          panel_sha256=seal['panel_sha256'],source_binding=binding,
          exporter_build=build_record,cms_preflight=version,
          collection=collection,rank_trajectory=rank_trajectory,
          matrix=snapshot,target_input=panel['target_input'],
          preparation_wall_ns=preparation_wall_ns,
          collection_wall_ns=collection_wall_ns,
          matrix_update_wall_ns=matrix_update_total_ns,
          relation_la_wall_ns=relation_la_wall_ns,
          column_logs=[dict(point=list(point),log=str(log))
                       for point,log in zip(matrix.columns,logs)],
          target_attempts=attempts,
          recovered_scalar=None if recovered is None else str(recovered),
          scalar_verified=recovered is not None,
          online_wall_ns=online_wall if recovered is not None else None,
          online_attempt_wall_ns=online_wall,online_phases_ns=phases,
          online_start_monotonic_ns=online_start,
          online_stop_monotonic_ns=online_stop,
          online_stop_event=('independent-scalar-replay' if recovered is not None
                             else 'frozen-target-attempt-cap'),
          online_bookkeeping_assigned_to_target_query_ns=residual,
          candidate_id=seal['candidate_id'],
          workload_id=seal['workload_id'],run_id=seal['run_id'],
          online_speedup=None),exclusive=True)

    return dict(status='COMPLETE' if recovered is not None else 'INCOMPLETE_TARGET',
                run_id=seal['run_id'],scalar_verified=recovered is not None)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--registration',type=Path,required=True)
    parser.add_argument('--expected-execution-sha256',required=True)
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args()
    spec=read(args.registration/'execution.json')
    require(sha256(spec)==args.expected_execution_sha256,
            'SAT dispatch differs from externally recorded invocation hash')
    receipt=execute(args.registration,args.out,expected_spec=spec,
                    timeout_seconds=spec['runtime_watchdog_seconds'])
    print(json.dumps(receipt,sort_keys=True))
