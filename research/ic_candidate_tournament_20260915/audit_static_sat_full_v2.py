#!/usr/bin/env python3
"""Replay the frozen SAT v2 run independently; never infer a paired speedup."""
import argparse
from collections import Counter
import json
from pathlib import Path

from audit_static_sat_full import (contents, digest, hash_to_curve,
                                   orbit_columns, raw_json, relation_row,
                                   verify_query)
from audit_static_cms_s4_natural import wilson
from generic_query_law import descent_coefficients, probe_scalar
from identity import candidate_manifest, run_id, workload_manifest
from oracle import Curve, rank, require
from run_generic_exact_yield_audit import (PANEL as EXACT_PANEL,
                                           load_evidence, pair_index, record)
from static_sat_registration_v2 import (PANEL, REGISTRATION, identities,
                                        source_manifest)
from tournament import read, write

HERE = Path(__file__).resolve().parent
PAIRED_TARGET = HERE/'goal_20260924/paired-fresh-n17a1/target-panel.json'


def audit(root):
    root = Path(root)
    require(root.is_dir(), 'SAT v2 run directory missing')
    panel, seal = read(PANEL), read(REGISTRATION/'seal.json')
    source, method, candidate, workload = identities(panel)
    require(source == source_manifest()
            and source == read(REGISTRATION/'source-manifest.json')
            and method == read(REGISTRATION/'method.json')
            and candidate == read(REGISTRATION/'candidate.json')
            and workload == read(REGISTRATION/'workload.json')
            and panel['candidate_id'] == candidate['candidate_id']
            and panel['workload_id'] == workload['workload_id']
            and panel['run_id'] == seal['run_id']
            and seal['run_id'] == run_id(candidate['candidate_id'],
                                         workload['workload_id'], 0),
            'frozen SAT v2 registration does not reconstruct')
    for copied, registered in (
            ('registered-panel.json', PANEL),
            ('paired-target-panel.json', PAIRED_TARGET),
            ('registered-seal.json', REGISTRATION/'seal.json'),
            ('candidate.json', REGISTRATION/'candidate.json'),
            ('workload.json', REGISTRATION/'workload.json'),
            ('method.json', REGISTRATION/'method.json'),
            ('source-manifest.json', REGISTRATION/'source-manifest.json'),
            ('PROTOCOL.md', REGISTRATION/'PROTOCOL.md'),
            ('registered-runner.py', HERE/'run_static_sat_full_v2.py'),
            ('static_sat_matrix.py', HERE/'static_sat_matrix.py'),
            ('sat-stage-runner.py', HERE/'static_sat_query_v2.py')):
        require(contents(root, copied) == registered.read_bytes(),
                'run differs from registered code or input: '+copied)
    require(digest(contents(root, 'cms-executable'))
                == panel['cms_executable_sha256']
            and digest(contents(root, 'cms-build-receipt.json'))
                == panel['cms_build_receipt_sha256']
            and digest(contents(root, 'cms-build-bundle-seal.json'))
                == panel['cms_build_bundle_seal_sha256']
            and '@rpath' not in contents(root, 'cms-linkage.txt').decode(),
            'retained CryptoMiniSat binary or build receipt differs')
    preflight = raw_json(root, 'cms-preflight.metrics.json')
    require(preflight['returncode'] == 0 and not preflight['timed_out']
            and b'CryptoMiniSat version 5.14.7'
                in contents(root, 'cms-preflight.stdout'),
            'SAT executable has no successful version preflight')
    binding = raw_json(root, 'source-binding.json')
    require(binding['manifest'] == source
            and type(binding['source_commit']) is str
            and len(binding['source_commit']) == 40,
            'source binding differs from frozen implementation')

    files = load_evidence(read(EXACT_PANEL))
    parent = record(files, 'jobs/n17a1/f5/stdout.json')
    parent_source = record(files, 'build/source-manifest.json')
    require(raw_json(root, 'parent-source-manifest.json') == parent_source,
            'retained Rust source/dependency receipt differs')
    curve = Curve(parent['fixture'])
    base = [curve.decode(point) for point in parent['factor_base']]
    require(candidate_manifest(parent['fixture'], parent, method) == candidate,
            'candidate digest does not replay from actual factor base')
    counter, target = hash_to_curve(
        curve, panel['target_input']['domain'].encode(),
        panel['target_input']['seed'])
    paired = read(PAIRED_TARGET)
    require(digest(contents(root, 'paired-target-panel.json'))
                == panel['paired_target_panel_sha256']
            and paired['target_input']['point'] == list(target)
            and paired['target_input']['counter'] == counter
            and panel['target_input']['point'] == list(target)
            and panel['target_input']['counter'] == counter
            and workload['record']['targets'] == [list(target)]
            and curve.mul(target, curve.r) is None,
            'one fresh paired public target does not replay')
    fixture = dict(parent['fixture'], targets=[list(target)],
                   target_seeds=[panel['target_input']['seed']],
                   target_scalar_constructed=False)
    require(workload_manifest(
                fixture, input_law=workload['record']['input_law'],
                algorithm_seed=workload['record']['algorithm_seed'],
                resource_envelope=workload['record']['resource_envelope'],
                cache_policy=workload['record']['cache_policy']) == workload,
            'one-target workload digest does not replay')
    columns, mapping = orbit_columns(curve, base)
    require(len(base) == 63
            and len({curve.mul(point, curve.h) for point in base}-{None}) == 62
            and len(columns) == 29,
            'actual or folded factor-base count changed')
    result = raw_json(root, 'summary.json')
    exporter_build = raw_json(root, 'build/build-record.json')
    require(result['exporter_build'] == exporter_build
            and exporter_build['source_commit'] == panel['source_commit']
            and exporter_build['exporter_source_sha256']
                == panel['exporter_source_sha256']
            and exporter_build['exporter_executable_sha256']
                == digest(contents(root, 'build/exporter'))
            and exporter_build['build_log_sha256']
                == digest(contents(root, 'build/build.log'))
            and raw_json(root, 'build/build-exit.json')['exit_code'] == 0,
            'Rust exporter build/source receipt differs')
    require(result['status'] in {'COMPLETE', 'INCOMPLETE_RELATION_RANK',
                                 'INCOMPLETE_TARGET'}
            and result['candidate_id'] == candidate['candidate_id']
            and result['workload_id'] == workload['workload_id']
            and result['run_id'] == seal['run_id']
            and result['panel_sha256'] == seal['panel_sha256']
            and result['source_binding'] == binding
            and result['cms_preflight'] == preflight
            and result['target_input'] == panel['target_input']
            and result['online_speedup'] is None,
            'result claims a different run or unsupported speedup')

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
            root, f'collection/trial-{trial:02d}/', row,
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

    common = dict(schema_version=1, candidate_id=candidate['candidate_id'],
                  workload_id=workload['workload_id'], run_id=seal['run_id'],
                  attempts=len(collection), statuses=dict(statuses),
                  exact_feasible=exact_feasible,
                  verified_relations=statuses['VALID_POINT_WITNESS'],
                  witness_rate_wilson95=wilson(
                      statuses['VALID_POINT_WITNESS'], len(collection)),
                  final_rank=snapshot['rank'], columns=len(columns),
                  same_point_rho_audited=False,
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
            root, f"descent/trial-{row['trial']:02d}/", row,
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
            and sum(phases.values()) == result['online_wall_ns']
            and type(result['online_start_monotonic_ns']) is int
            and type(result['online_stop_monotonic_ns']) is int
            and result['online_stop_monotonic_ns']
                - result['online_start_monotonic_ns'] == result['online_wall_ns']
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
                online_phases_ns=phases,
                online_endpoint_admissible=complete)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run-dir', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    require(not args.out.exists(), 'audit output already exists')
    write(args.out, audit(args.run_dir), exclusive=True)


if __name__ == '__main__':
    main()
