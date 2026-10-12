#!/usr/bin/env python3
"""One-shot complete n17 static-SAT relation/LA/descent IC attempt."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import time

from generic_bases import lifts
from generic_query_law import descent_coefficients, probe_scalar
from generic_solver_feasibility import check_source_checkout
from identity import candidate_manifest, run_id, workload_manifest
from oracle import require
from run_cms_s4_controls import EXPORTER, build_exporter, digest, meter
from run_static_cms_s4_controls import BUNDLE, SCRIPTS
from run_static_cms_s4_natural import (
    PANEL as NATURAL_PANEL, admit as admit_natural, one_query,
)
from static_sat_matrix import RelationMatrix
from static_sat_registration import source_manifest
from tournament import read, write

HERE = Path(__file__).resolve().parent
REGISTRATION = HERE/'goal_20260924/static-sat-full'
PANEL = REGISTRATION/'panel.json'
SEAL = REGISTRATION/'seal.json'


def target_from_seed(curve, seed):
    for counter in range(1_000_000):
        data = hashlib.sha256(
            b'ic-static-sat-target-v1\0'
            +seed.to_bytes(8, 'little')+counter.to_bytes(8, 'little')).digest()
        x = int.from_bytes(data[:8], 'little') & ((1 << curve.n)-1)
        choices = lifts(curve, x)
        if not choices:
            continue
        point = curve.mul(choices[data[8] % len(choices)], curve.h)
        if point is not None:
            return counter, point
    raise ValueError('target hash-to-curve counter limit exhausted')


def admit(panel, *, require_local_binary=True):
    seal = read(SEAL)
    require(digest(PANEL) == seal['panel_sha256'], 'complete SAT panel changed')
    require(digest(REGISTRATION/'candidate.json') == seal['candidate_sha256']
            and digest(REGISTRATION/'workload.json') == seal['workload_sha256']
            and digest(REGISTRATION/'method.json') == seal['method_sha256']
            and digest(REGISTRATION/'source-manifest.json')
                == seal['source_manifest_file_sha256'],
            'complete SAT identity registration changed')
    require(source_manifest() == read(REGISTRATION/'source-manifest.json'),
            'complete SAT executed source changed')
    candidate = read(REGISTRATION/'candidate.json')
    workload = read(REGISTRATION/'workload.json')
    require(panel['candidate_id'] == candidate['candidate_id']
            and panel['workload_id'] == workload['workload_id']
            and panel['run_id'] == seal['run_id']
            and seal['candidate_id'] == candidate['candidate_id']
            and seal['workload_id'] == workload['workload_id']
            and seal['run_id'] == run_id(candidate['candidate_id'],
                                         workload['workload_id'], 0),
            'complete SAT identity label changed')
    parent = read(NATURAL_PANEL)
    files, curve, base, built = admit_natural(
        parent, require_local_binary=require_local_binary)
    require(panel['schema_version'] == 1
            and panel['status'] == 'REGISTERED_BEFORE_EXECUTION'
            and panel['field_degree'] == 17 and panel['curve_a'] == 1
            and panel['irreducible_low_terms'] == [0, 3]
            and panel['factor_base'] == {'kind':'standard_subspace','dimension':6}
            and panel['summands'] == 3
            and panel['source_encoding'] == 'wide-symmetrised-S4-circuit-XOR-DIMACS'
            and panel['relation_query_seed'] == 2026092936
            and panel['max_relation_queries'] == 256
            and panel['descent_query_seed'] == 2026092937
            and panel['max_descent_queries'] == 64
            and panel['export_nonce'] == 2026092939
            and panel['export_timeout_seconds'] == 60
            and panel['cms_timeout_seconds'] == 120
            and panel['cms_conflict_budget'] == 1000000
            and panel['cms_threads'] == panel['cms_random_seed']
                == panel['cms_max_models_per_query'] == 1,
            'complete SAT method or frozen attempt limits changed')
    for name in ('source_commit', 'exporter_source_sha256', 'cms_path',
                 'cms_executable_sha256', 'cms_build_receipt_sha256',
                 'cms_build_bundle_seal_sha256', 'parent_exact_result_sha256'):
        require(panel[name] == parent[name],
                'complete SAT method differs from source-admitted stage: '+name)
    require(digest(NATURAL_PANEL.parent/'evidence.tar.gz')
            == panel['parent_natural_archive_sha256']
            and digest(NATURAL_PANEL.parent/'RESULT.json')
                == panel['parent_natural_result_sha256'],
            'parent natural SAT evidence changed')
    counter, target = target_from_seed(curve, panel['target_input']['seed'])
    require(panel['target_input'] == {
                'law':'sha256-x-lift-cofactor-v1',
                'seed':2026092938, 'counter':counter,
                'point':list(target),
                'point_was_previously_supplied':False,
                'known_scalar_supplied':False},
            'frozen one-target public point does not replay')
    prior = {tuple(item['point']) for item in parent['schedule']}
    prior.update(tuple(item['point']) for item in
                 read(HERE/'goal_20260924/static-cms-s4-controls/panel.json')['schedule'])
    prior.add((853,39791))
    require(target not in prior and curve.mul(target,curve.r) is None,
            'target repeats a previous public point or leaves subgroup')
    matrix = RelationMatrix(curve, base)
    require(len(base) == 63
            and sum(curve.mul(point,curve.h) is not None for point in base) == 62
            and len(matrix.columns) == 29,
            'actual base size or folded matrix width changed')
    from run_generic_exact_yield_audit import record
    report = record(files, 'jobs/n17a1/f5/stdout.json')
    old = read(HERE/'goal_20260924/generic-exact-yield-audit/RESULT.json')
    old_f5 = next(row for row in old['rows']
                  if row['cell']=='n17a1' and row['solver']=='f5')
    prior_scalars={item['a'] for item in old_f5['queries']}
    prior_scalars.update(item['probe_scalar'] for item in parent['schedule'])
    planned=[probe_scalar(panel['relation_query_seed'], trial, curve.r)
             for trial in range(panel['max_relation_queries'])]
    require(len(set(planned))==len(planned)
            and not prior_scalars.intersection(planned),
            'complete SAT collection repeats a prior or internal ordinary query')
    fixture = dict(report['fixture'], targets=[list(target)],
                   target_seeds=[panel['target_input']['seed']],
                   target_scalar_constructed=False)
    require(candidate_manifest(report['fixture'],report,
                               read(REGISTRATION/'method.json')) == candidate
            and workload_manifest(fixture,
                input_law=workload['record']['input_law'],
                algorithm_seed=workload['record']['algorithm_seed'],
                resource_envelope=workload['record']['resource_envelope'],
                cache_policy=workload['record']['cache_policy']) == workload,
            'complete SAT candidate or workload fails independent reconstruction')
    return files, curve, base, target, built, matrix, seal


def source_binding(out):
    manifest = source_manifest()
    require(manifest == read(REGISTRATION/'source-manifest.json'),
            'complete SAT source manifest changed before execution')
    commit = subprocess.check_output(
        ['git','rev-parse','HEAD'], cwd=HERE, text=True).strip()
    dirty = subprocess.check_output(
        ['git','status','--porcelain=v1','--untracked-files=all'],
        cwd=HERE, text=True).strip()
    require(not dirty, 'complete SAT controller requires a clean committed source tree')
    record = dict(source_commit=commit, manifest=manifest)
    write(out/'source-binding.json',record,exclusive=True)
    return record


def save_progress(out, name, row):
    with (out/name).open('a') as stream:
        stream.write(json.dumps(row,sort_keys=True)+'\n')


def run(panel, source_root, out):
    preparation_start = time.monotonic_ns()
    files, curve, base, target, built, matrix, seal = admit(panel)
    check_source_checkout(source_root)
    require(digest(source_root/EXPORTER) == panel['exporter_source_sha256'],
            'pinned Rust source exporter changed')
    require(not out.exists(), 'complete SAT output exists; no job may be retried')
    out.mkdir(parents=True)
    for source, name in ((PANEL,'registered-panel.json'),
                         (SEAL,'registered-seal.json'),
                         (REGISTRATION/'candidate.json','candidate.json'),
                         (REGISTRATION/'workload.json','workload.json'),
                         (REGISTRATION/'method.json','method.json'),
                         (REGISTRATION/'source-manifest.json','source-manifest.json'),
                         (REGISTRATION/'PROTOCOL.md','PROTOCOL.md'),
                         (Path(__file__),'registered-runner.py'),
                         (HERE/'static_sat_matrix.py','static_sat_matrix.py'),
                         (Path(one_query.__code__.co_filename),'sat-stage-runner.py'),
                         (SCRIPTS/'process_meter.py','process_meter.py'),
                         (SCRIPTS/'run_koblitz_pdp_matrix.py','cms_parser.py'),
                         (BUNDLE/'builds/cryptominisat-receipt.json',
                          'cms-build-receipt.json'),
                         (BUNDLE/'bundle-seal.json','cms-build-bundle-seal.json')):
        shutil.copy2(source,out/name)
    binding = source_binding(out)
    write(out/'host.json',dict(system=platform.system(),
                               machine=platform.machine(),
                               cpu_count=os.cpu_count(),
                               scope='one local source-bound full IC attempt'),
          exclusive=True)
    write(out/'cms-build-verification.json',built,exclusive=True)
    from run_generic_exact_yield_audit import record
    parent_source=record(files,'build/source-manifest.json')
    write(out/'parent-source-manifest.json',parent_source,exclusive=True)
    cms=out/'cms-executable'
    shutil.copy2(panel['cms_path'],cms)
    require(digest(cms)==panel['cms_executable_sha256'],
            'copied static SAT solver changed')
    linkage=subprocess.check_output(['otool','-L',str(cms)],text=True)
    (out/'cms-linkage.txt').write_text(linkage)
    require('@rpath' not in linkage and '/opt/homebrew' not in linkage,
            'copied static solver has a nonportable dynamic dependency')
    version=meter([cms,'--version'],out,'cms-preflight',10)
    require(version['returncode']==0 and not version['timed_out']
            and 'CryptoMiniSat version 5.14.7'
                in (out/'cms-preflight.stdout').read_text(),
            'copied solver failed execution preflight')
    build_dir=out/'build'
    build_dir.mkdir()
    exporter,build_record=build_exporter(panel,source_root,build_dir,parent_source)

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
        row=one_query(panel,item,exporter,cms,curve,base,out/'collection')
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
        write(out/'summary.json',dict(schema_version=1,
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
        return
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
    stream=descent_coefficients(panel['descent_query_seed'],curve.r,walked=False)
    for index in range(panel['max_descent_queries']):
        started=time.monotonic_ns()
        a,b=next(stream)
        query=curve.add(curve.mul(curve.g,a),curve.mul(target,b))
        phases['target_query']+=time.monotonic_ns()-started
        item=dict(trial=panel['max_relation_queries']+index,
                  probe_scalar=a,point=list(query))
        started=time.monotonic_ns()
        row=one_query(panel,item,exporter,cms,curve,base,out/'descent')
        phases['target_pdp']+=time.monotonic_ns()-started
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
            row['candidate_scalar']=str(candidate)
            row['scalar_replay_verified']=verified
            require(verified,'target descent scalar fails independent replay')
            recovered=candidate
        attempts.append(row)
        save_progress(out,'descent.progress.jsonl',row)
        print(json.dumps(dict(stage='target',trial=index,status=row['status'],
                              verified=recovered is not None)),flush=True)
        if recovered is not None:
            break
    online_wall=time.monotonic_ns()-online_start
    residual=online_wall-sum(phases.values())
    require(residual>=0,'exclusive target phase clocks overlap')
    phases['target_query']+=residual
    require(sum(phases.values())==online_wall,'online phase ledger does not close')
    write(out/'summary.json',dict(schema_version=1,
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
          online_wall_ns=online_wall,online_phases_ns=phases,
          online_bookkeeping_assigned_to_target_query_ns=residual,
          candidate_id=seal['candidate_id'],
          workload_id=seal['workload_id'],run_id=seal['run_id'],
          online_speedup=None),exclusive=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-root',type=Path,required=True)
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args()
    out=args.out.resolve()
    try:
        run(read(PANEL),args.source_root.resolve(),out)
    except Exception as error:
        if out.is_dir() and not (out/'summary.json').exists():
            seal=read(SEAL)
            write(out/'summary.json',dict(schema_version=1,status='ERROR',
                  error_type=type(error).__name__,error=str(error),
                  candidate_id=seal['candidate_id'],
                  workload_id=seal['workload_id'],run_id=seal['run_id'],
                  progress_files=['collection.progress.jsonl',
                                  'descent.progress.jsonl'],
                  recovered_scalar=None,scalar_verified=False,
                  online_wall_ns=None,online_phases_ns=None,
                  online_speedup=None),exclusive=True)
        raise


if __name__=='__main__':
    main()
