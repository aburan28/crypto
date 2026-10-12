#!/usr/bin/env python3
"""Independently replay a registered static-SAT one-target IC run directory."""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path

from audit_static_cms_s4_natural import (
    check_formula, complete_model, group_lift, wilson,
)
from generic_query_law import descent_coefficients, probe_scalar
from identity import candidate_manifest, run_id, workload_manifest
from oracle import Curve, rank, require
from run_generic_exact_yield_audit import (PANEL as EXACT_PANEL, exact_three_sum,
                                           load_evidence, pair_index, record)
from static_sat_registration import (PANEL, REGISTRATION, identities,
                                     source_manifest)
from tournament import read, write


def digest(data):
    return hashlib.sha256(data).hexdigest()


def contents(root, name):
    path = root/name
    require(path.is_file() and not path.is_symlink()
            and path.resolve().is_relative_to(root.resolve()),
            'missing or unsafe full-SAT evidence: '+name)
    return path.read_bytes()


def raw_json(root, name):
    return json.loads(contents(root, name))


def hash_to_curve(curve, domain, seed):
    from generic_bases import lifts
    for counter in range(1_000_000):
        data = hashlib.sha256(domain+b'\0'
                              +seed.to_bytes(8,'little')
                              +counter.to_bytes(8,'little')).digest()
        x = int.from_bytes(data[:8],'little') & ((1<<curve.n)-1)
        choices = lifts(curve,x)
        if choices:
            point=curve.mul(choices[data[8] % len(choices)],curve.h)
            if point is not None:
                return counter,point
    raise ValueError('no target from frozen hash-to-curve law')


def target_from_seed(curve, seed):
    return hash_to_curve(curve,b'ic-static-sat-target-v1',seed)


def orbit_columns(curve, base):
    """Build quotient coefficients without importing the producer row builder."""
    def orbit(point):
        seen=set()
        for _ in range(curve.n):
            seen.add(point)
            seen.add(curve.neg(point))
            point=curve.frob(point)
        return seen
    projected=[curve.mul(point,curve.h) for point in base]
    columns=sorted({min(orbit(point)) for point in projected
                    if point is not None})
    mapping={}
    for index,start in enumerate(columns):
        point=start
        coefficient=1
        for _ in range(curve.n):
            for image,value in ((point,coefficient),
                                (curve.neg(point),-coefficient%curve.r)):
                present=mapping.get(image)
                require(present is None or present==(index,value),
                        'ambiguous independently built quotient coefficient')
                mapping[image]=(index,value)
            point=curve.frob(point)
            coefficient=coefficient*curve.lam%curve.r
    require(set(projected)-{None} <= set(mapping),
            'projected factor-base point has no quotient coefficient')
    return columns,[mapping.get(point) for point in projected]


def relation_row(curve, base, mapping, columns, scalar, indices):
    require(len(indices)==3 and all(type(index) is int
            and 0 <= index < len(base) for index in indices),
            'malformed three-summand point witness')
    group=None
    row=[0]*len(columns)
    for index in indices:
        group=curve.add(group,base[index])
        if mapping[index] is not None:
            column,coefficient=mapping[index]
            row[column]=(row[column]+coefficient)%curve.r
    require(group==curve.mul(curve.g,scalar),
            'ordinary relation fails independent group replay')
    return row,curve.h*scalar%curve.r


def audited_status(raw_status, cms, stdout):
    """Keep the producer label, but identify CryptoMiniSat's budget exit."""
    if (raw_status=='SOLVER_ERROR' and cms is not None
            and cms['returncode']==15 and not cms['timed_out']
            and b's INDETERMINATE' in stdout
            and '--maxconfl' in cms['command']):
        return 'CONFLICT_BUDGET_INCONCLUSIVE'
    return raw_status


def verify_query(root, prefix, row, point, curve, base, pairs):
    """Recheck the raw source system, solver receipt and exact group yield."""
    require(row['public_point']==list(point),
            'run row has the wrong public query point')
    export=raw_json(root,prefix+'export.metrics.json')
    require(export==row['exporter'], 'export process receipt differs from row')
    status=row['status']
    exact=exact_three_sum(curve,base,pairs,point)
    cms=row['cms']
    if cms is not None:
        require(cms==raw_json(root,prefix+'cms.metrics.json'),
                'solver process receipt differs from row')
    if status in {'EXPORT_FAILURE','INVALID_EXPORT'}:
        require(cms is None, 'failed export unexpectedly ran SAT')
        return exact is not None,None,status
    manifest_raw=contents(root,prefix+'instance/manifest.json')
    manifest=json.loads(manifest_raw)
    require(row['manifest_sha256']==digest(manifest_raw)
            and manifest['representation']=='symmetrised_s4'
            and manifest['n']==curve.n and manifest['curve_a']==curve.a
            and [int(manifest['target'][key]) for key in ('x','y')]==list(point)
            and manifest['factor_base_geometry']['curve_points']==len(base),
            'source instance differs from the public query or factor base')
    for name,descriptor in manifest['exports'].items():
        require(Path(descriptor['path']).name==descriptor['path'],
                'unsafe source-export descriptor path')
        raw=contents(root,prefix+'instance/'+descriptor['path'])
        require(row['exports'][name]=={'bytes':len(raw),'sha256':digest(raw)},
                'source-export bytes differ from run row')
    require(cms is not None,'source export lacks solver receipt')
    stdout=contents(root,prefix+'cms.stdout')
    interpreted=audited_status(status,cms,stdout)
    if status in {'VALID_POINT_WITNESS','SOURCE_MODEL_NONLIFTING'}:
        require(cms['returncode']==10 and not cms['timed_out']
                and b's SATISFIABLE' in stdout,
                'SAT witness has no solver verdict')
        count=manifest['exports']['cryptominisat_xor_dimacs']['variables']
        model=complete_model(stdout,count)
        check_formula(contents(root,prefix+'instance/instance.xor.cnf'),model)
        xs,lifted=group_lift(curve,base,model,point)
        require(row['source_model_valid'] is True
                and row['source_model_sha256']==digest(bytes(model)),
                'model receipt differs from complete source assignment')
        if status=='VALID_POINT_WITNESS':
            require(lifted is not None and exact is not None
                    and row['point_witness']['x_coordinates']==xs
                    and row['point_witness']['points']==[list(p) for p in lifted]
                    and row['point_witness']['point_indices']
                        ==[base.index(p) for p in lifted],
                    'witness fails source, lift or independent group replay')
            return True,row['point_witness']['point_indices'],interpreted
        require(lifted is None,'reported nonlifting model has a group lift')
    elif status=='SOURCE_UNSAT':
        require(cms['returncode']==20 and not cms['timed_out']
                and b's UNSATISFIABLE' in stdout and exact is None,
                'solver UNSAT contradicts exact three-sum oracle')
    elif status=='TIMEOUT':
        require(cms['timed_out'],'timeout has no process watchdog evidence')
    elif status=='CONFLICT_BUDGET_INCONCLUSIVE':
        require(cms['returncode']==15 and not cms['timed_out']
                and b's INDETERMINATE' in stdout
                and '--maxconfl' in cms['command'],
                'censored conflict-budget row lacks solver evidence')
    else:
        require(status in {'SOLVER_ERROR','UNKNOWN_INCONCLUSIVE',
                           'INVALID_SOURCE_MODEL'},
                'unknown or unaccounted SAT status')
    return exact is not None,None,interpreted


def audit(root):
    root=Path(root)
    require(root.is_dir(),'full-SAT run directory missing')
    panel=read(PANEL)
    seal=read(REGISTRATION/'seal.json')
    source,method,candidate,workload=identities(panel)
    require(source==source_manifest()
            and candidate==read(REGISTRATION/'candidate.json')
            and workload==read(REGISTRATION/'workload.json')
            and method==read(REGISTRATION/'method.json')
            and seal['candidate_id']==candidate['candidate_id']
            and seal['workload_id']==workload['workload_id']
            and seal['run_id']==run_id(candidate['candidate_id'],
                                       workload['workload_id'],0),
            'frozen candidate or workload identity changed')
    for name,registered in (
            ('registered-panel.json','panel.json'),
            ('registered-seal.json','seal.json'),
            ('candidate.json','candidate.json'),
            ('workload.json','workload.json'),
            ('method.json','method.json'),
            ('source-manifest.json','source-manifest.json')):
        require(contents(root,name)==(REGISTRATION/registered).read_bytes(),
                'run used a different registration file: '+name)
    require(contents(root,'registered-runner.py')
            ==(PANEL.parents[2]/'run_static_sat_full.py').read_bytes()
            and contents(root,'static_sat_matrix.py')
                ==(PANEL.parents[2]/'static_sat_matrix.py').read_bytes(),
            'run used a different controller or relation-matrix source')
    require(digest(contents(root,'cms-executable'))
                ==panel['cms_executable_sha256']
            and digest(contents(root,'cms-build-receipt.json'))
                ==panel['cms_build_receipt_sha256']
            and digest(contents(root,'cms-build-bundle-seal.json'))
                ==panel['cms_build_bundle_seal_sha256']
            and '@rpath' not in contents(root,'cms-linkage.txt').decode(),
            'solver binary/build receipt/linkage differs from registration')
    preflight=raw_json(root,'cms-preflight.metrics.json')
    require(preflight['returncode']==0 and not preflight['timed_out']
            and b'CryptoMiniSat version 5.14.7'
                in contents(root,'cms-preflight.stdout'),
            'copied solver did not pass execution preflight')
    binding=raw_json(root,'source-binding.json')
    require(binding['manifest']==source
            and binding['source_commit']
                == 'fa1e9f8fcc68cbc63ce75861f3989bbe16b5862d',
            'complete runner source receipt changed')
    files=load_evidence(read(EXACT_PANEL))
    parent=record(files,'jobs/n17a1/f5/stdout.json')
    parent_source=record(files,'build/source-manifest.json')
    require(raw_json(root,'parent-source-manifest.json')==parent_source,
            'pinned Rust source/dependency manifest changed')
    curve=Curve(parent['fixture'])
    base=[curve.decode(p) for p in parent['factor_base']]
    report=dict(parent,columns=29)
    require(candidate_manifest(parent['fixture'],report,method)==candidate,
            'candidate identity does not replay from actual source factor base')
    counter,target=target_from_seed(curve,panel['target_input']['seed'])
    require(panel['target_input']['counter']==counter
            and panel['target_input']['point']==list(target)
            and workload['record']['targets']==[list(target)],
            'frozen public target does not replay')
    expected_fixture=dict(parent['fixture'],targets=[list(target)],
                          target_seeds=[panel['target_input']['seed']],
                          target_scalar_constructed=False)
    require(workload_manifest(expected_fixture,
            input_law=workload['record']['input_law'],
            algorithm_seed=workload['record']['algorithm_seed'],
            resource_envelope=workload['record']['resource_envelope'],
            cache_policy=workload['record']['cache_policy'])==workload,
            'one-target workload identity does not replay')
    columns,mapping=orbit_columns(curve,base)
    require(len(base)==63 and len(set(curve.mul(p,curve.h) for p in base)-{None})==62
            and len(columns)==29,
            'geometric, usable or folded factor-base count changed')
    result=raw_json(root,'summary.json')
    exporter_build=raw_json(root,'build/build-record.json')
    require(result['exporter_build']==exporter_build
            and exporter_build['source_commit']==panel['source_commit']
            and exporter_build['exporter_source_sha256']
                ==panel['exporter_source_sha256']
            and exporter_build['exporter_executable_sha256']
                ==digest(contents(root,'build/exporter'))
            and exporter_build['build_log_sha256']
                ==digest(contents(root,'build/build.log'))
            and raw_json(root,'build/build-exit.json')['exit_code']==0,
            'source-bound Rust exporter build receipt changed')
    require(result['candidate_id']==candidate['candidate_id']
            and result['workload_id']==workload['workload_id']
            and result['run_id']==seal['run_id']
            and result['panel_sha256']==seal['panel_sha256']
            and result['source_binding']==binding
            and result['cms_preflight']==preflight
            and result['online_speedup'] is None,
            'result identity, source binding or claim boundary differs')
    require(result['status'] in {'COMPLETE','INCOMPLETE_RELATION_RANK',
                                 'INCOMPLETE_TARGET'},
            'failed run requires a separate partial/error audit')
    collection=result['collection']
    require(0<len(collection)<=panel['max_relation_queries']
            and [json.loads(line) for line in
                 contents(root,'collection.progress.jsonl').decode().splitlines()]
                ==collection,
            'ordinary progress differs from summary')
    pairs=pair_index(curve,base)
    statuses=Counter()
    audited_statuses=Counter()
    exact_feasible=0
    rank_trajectory=[]
    matrix_rows=[]
    duplicates=dependencies=0
    keys=set()
    for trial,row in enumerate(collection):
        scalar=probe_scalar(panel['relation_query_seed'],trial,curve.r)
        point=curve.mul(curve.g,scalar)
        require(row['trial']==trial and row['probe_scalar']==scalar,
                'ordinary query scalar differs from frozen law')
        feasible,indices,interpreted=verify_query(
            root,f'collection/trial-{trial:02d}/',row,point,curve,base,pairs)
        exact_feasible+=feasible
        statuses[row['status']]+=1
        audited_statuses[interpreted]+=1
        before=rank([item[0] for item in matrix_rows],len(columns),curve.r)
        if indices is not None:
            coeff,rhs=relation_row(curve,base,mapping,columns,scalar,indices)
            key=scalar,tuple(sorted(indices))
            if key in keys:
                expected='duplicate'
                duplicates+=1
            else:
                keys.add(key)
                matrix_rows.append((coeff,rhs,scalar,indices))
                after=rank([item[0] for item in matrix_rows],len(columns),curve.r)
                expected='novel_rank' if after>before else 'dependent'
                dependencies+=after==before
            require(row['relation_status']==expected,
                    'relation rank or duplicate classification changed')
        else:
            require(row['relation_status']=='no_verified_relation',
                    'failed ordinary query produced a matrix row')
        actual_rank=rank([item[0] for item in matrix_rows],len(columns),curve.r)
        require(row['rank_after']==actual_rank
                and row['attempt_wall_ns']>=row['matrix_update_wall_ns']>=0,
                'ordinary rank or measured attempt interval changed')
        rank_trajectory.append(actual_rank)
    require(result['rank_trajectory']==rank_trajectory,
            'rank trajectory changed')
    snapshot=raw_json(root,'relation-matrix.json')
    rows=[dict(entries=[[index,str(value)]
                        for index,value in enumerate(coeff) if value],
               rhs=str(rhs),scalar=scalar,indices=list(indices))
          for coeff,rhs,scalar,indices in matrix_rows]
    require(snapshot==result['matrix']
            and snapshot['column_points']==[list(p) for p in columns]
            and snapshot['columns']==len(columns)
            and snapshot['rank']==rank_trajectory[-1]
            and snapshot['accepted_rows']==len(matrix_rows)
            and snapshot['duplicate_relations']==duplicates
            and snapshot['dependent_relations']==dependencies
            and snapshot['rows']==rows
            and snapshot['rows_sha256']==digest(json.dumps(
                rows,sort_keys=True,separators=(',',':')).encode()),
            'final relation matrix differs from independent row reconstruction')
    require(result['matrix_update_wall_ns']
                ==sum(row['matrix_update_wall_ns'] for row in collection)
            and result['collection_wall_ns']>=result['matrix_update_wall_ns']
            and result['preparation_wall_ns']>=result['collection_wall_ns'],
            'preparation or matrix diagnostic timing is inconsistent')
    if snapshot['rank']<len(columns):
        require(result['status']=='INCOMPLETE_RELATION_RANK'
                and result['target_attempts']==[]
                and result['online_wall_ns'] is None
                and result['recovered_scalar'] is None,
                'rank-deficient attempt overclaims target recovery')
        return dict(schema_version=1,status='AUDITED_INCOMPLETE_RANK',
                    candidate_id=candidate['candidate_id'],
                    workload_id=workload['workload_id'],run_id=seal['run_id'],
                    attempts=len(collection),statuses=dict(statuses),
                    audited_statuses=dict(audited_statuses),
                    exact_feasible=exact_feasible,
                    verified_relations=statuses['VALID_POINT_WITNESS'],
                    final_rank=snapshot['rank'],columns=len(columns),
                    solved_targets=0,online_wall_ns=None,
                    headline_online_admissible=False,online_speedup=None)
    require(result['relation_la_wall_ns'] is not None
            and result['relation_la_wall_ns']>=0,
            'full-rank matrix has no final LA diagnostic')
    logs=[int(item['log']) for item in result['column_logs']]
    retained_logs=raw_json(root,'column-logs.json')
    require(len(logs)==len(columns)
            and [item['point'] for item in result['column_logs']]
                ==[list(point) for point in columns]
            and all(curve.mul(curve.g,log)==point
                    for log,point in zip(logs,columns))
            and all(sum(a*b for a,b in zip(coeff,logs))%curve.r==rhs
                    for coeff,rhs,_,_ in matrix_rows)
            and retained_logs==dict(columns=result['column_logs'],
                                    independently_verified=True),
            'solved column logs fail independent group or matrix checks')
    target_rows=result['target_attempts']
    require(0<len(target_rows)<=panel['max_descent_queries']
            and [json.loads(line) for line in
                 contents(root,'descent.progress.jsonl').decode().splitlines()]
                ==target_rows,
            'target progress differs from summary')
    stream=descent_coefficients(panel['descent_query_seed'],curve.r,walked=False)
    target_statuses=Counter()
    target_audited_statuses=Counter()
    recovered=None
    for index,row in enumerate(target_rows):
        a,b=next(stream)
        point=curve.add(curve.mul(curve.g,a),curve.mul(target,b))
        require(row['target_query_index']==index
                and row['trial']==panel['max_relation_queries']+index
                and row['probe_scalar']==a
                and (row['a'],row['b'])==(a,b),
                'target query differs from frozen aG+bQ law')
        _,indices,interpreted=verify_query(root,
            f"descent/trial-{row['trial']:02d}/",row,
            point,curve,base,pairs)
        target_statuses[row['status']]+=1
        target_audited_statuses[interpreted]+=1
        if indices is not None:
            projected=0
            for base_index in indices:
                if mapping[base_index] is not None:
                    column,coefficient=mapping[base_index]
                    projected=(projected+coefficient*logs[column])%curve.r
            recovered=(projected-curve.h*a)*pow(curve.h*b%curve.r,-1,curve.r)%curve.r
            require(row['candidate_scalar']==str(recovered)
                    and row['scalar_replay_verified'] is True
                    and curve.mul(curve.g,recovered)==target,
                    'target logarithm fails independent equation or scalar replay')
            require(index==len(target_rows)-1,
                    'target attempts continued after verified recovery')
    phases=result['online_phases_ns']
    require(set(phases)=={'target_query','target_pdp',
                          'target_relation_check','target_descent',
                          'target_recovery_check'}
            and all(type(v) is int and v>=0 for v in phases.values())
            and sum(phases.values())==result['online_wall_ns']
            and result['online_bookkeeping_assigned_to_target_query_ns']>=0,
            'one-target exclusive online phase ledger does not close')
    complete=recovered is not None
    require(result['status']==('COMPLETE' if complete else 'INCOMPLETE_TARGET')
            and result['scalar_verified'] is complete
            and result['recovered_scalar']
                ==(str(recovered) if complete else None),
            'one-target completion status overclaims verified recovery')
    return dict(schema_version=1,
                status=('AUDITED_COMPLETE_CORRECTNESS_ONLY' if complete
                        else 'AUDITED_INCOMPLETE_TARGET'),
                candidate_id=candidate['candidate_id'],
                workload_id=workload['workload_id'],run_id=seal['run_id'],
                attempts=len(collection),statuses=dict(statuses),
                audited_statuses=dict(audited_statuses),
                exact_feasible=exact_feasible,
                verified_relations=statuses['VALID_POINT_WITNESS'],
                witness_rate_wilson95=wilson(
                    statuses['VALID_POINT_WITNESS'],len(collection)),
                final_rank=snapshot['rank'],columns=len(columns),
                target_attempts=len(target_rows),
                target_statuses=dict(target_statuses),
                target_audited_statuses=dict(target_audited_statuses),
                solved_targets=int(complete),
                recovered_scalar=None if recovered is None else str(recovered),
                online_wall_ns=result['online_wall_ns'],
                online_phases_ns=phases,
                online_endpoint_limitation=(
                    'runner stops clock after final progress write, not immediately '
                    'after scalar replay; online wall is overinclusive'),
                headline_online_admissible=False,online_speedup=None)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run-dir',type=Path,required=True)
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args()
    require(not args.out.exists(),'audit output exists')
    write(args.out,audit(args.run_dir),exclusive=True)


if __name__=='__main__':
    main()
