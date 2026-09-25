#!/usr/bin/env python3
"""Point-only holdout recovery child: receives no k labels or source indices."""
from __future__ import annotations

import argparse
import collections
import importlib.util
import json
import resource
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
OLD_GATE = HERE.parent / 'rotated_subspace_support_20260925/gate.py'
Q = 2003
ARMS = ((3, 5), (3, 6), (7, 5), (7, 6))


def prior():
    spec = importlib.util.spec_from_file_location('rank_recovery_arithmetic', OLD_GATE)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def rss():
    v = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return v if sys.platform == 'darwin' else v * 1024


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--points', type=Path, required=True)
    parser.add_argument('--logs', type=Path, required=True)
    parser.add_argument('--oracle', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    start, cpu = time.monotonic(), time.process_time()
    points = json.loads(args.points.read_text())
    logs = json.loads(args.logs.read_text())
    oracle = json.loads(args.oracle.read_text())
    assert logs['q'] == Q and logs['generator'] == [4793, 2429]
    assert [(a['beta'], a['m']) for a in points['arms']] == list(ARMS)
    assert [(a['beta'], a['m']) for a in logs['arms']] == list(ARMS)
    assert [(a['beta'], a['m']) for a in oracle['arms']] == list(ARMS)
    mod = prior()
    field = mod.Field(13, mod.MODELS[13]['low'])
    curve = mod.Curve(field)
    H = (4793, 2429)
    assert curve.on_curve(H) and curve.scalar(H, Q) is None
    counters = collections.Counter()
    results = []
    for arm_points, arm_logs, arm_oracle in zip(points['arms'], logs['arms'], oracle['arms']):
        assert (arm_points['beta'], arm_points['m']) == (arm_logs['beta'], arm_logs['m'])
        assert (arm_points['beta'], arm_points['m']) == (arm_oracle['beta'], arm_oracle['m'])
        columns = [tuple(c) for c in arm_logs['columns']]
        by_column = {c: i for i, c in enumerate(columns)}
        assert len(by_column) == len(columns)
        base_logs = arm_logs['logs']
        if base_logs is not None:
            assert len(base_logs) == len(columns)
        cases = []
        assert len(arm_points['cases']) == len(arm_oracle['responses']) == 16
        for query, oracle_record in zip(arm_points['cases'], arm_oracle['responses']):
            assert query['case_id'] == oracle_record['case_id']
            assert query['point'] == oracle_record['point']
            Qpoint = tuple(query['point'])
            assert curve.on_curve(Qpoint)
            if base_logs is None:
                result = {'case_id': query['case_id'], 'point': query['point'],
                          'status': 'rank_deficient', 'recovered_k': None}
            elif oracle_record['lookup_status'] != 'hit':
                assert oracle_record['lookup_status'] == 'miss' and oracle_record['response'] is None
                result = {'case_id': query['case_id'], 'point': query['point'],
                          'status': 'oracle_miss', 'recovered_k': None}
            else:
                response = oracle_record['response']
                assert response['point'] == query['point']
                assert 0 <= response['torsion_index'] < 4
                total = 0
                for term in response['row']:
                    index = by_column[tuple(term['column'])]
                    total = (total + term['coefficient'] * base_logs[index]) % Q
                    counters['modular_mul'] += 1
                    counters['modular_add'] += 1
                scalar = (total * pow(4, -1, Q)) % Q
                counters['modular_mul'] += 1
                assert curve.scalar(H, scalar) == Qpoint
                result = {'case_id': query['case_id'], 'point': query['point'],
                          'status': 'recovered_and_group_verified', 'recovered_k': scalar,
                          'oracle_torsion_index': response['torsion_index']}
            cases.append(result)
        results.append({'beta': arm_points['beta'], 'm': arm_points['m'], 'cases': cases})
    report = {'schema': 'rotated_point_only_recovery_v1', 'arms': results,
              'modular_counters': dict(counters),
              'field_operations': dict(field.operations),
              'curve_operations': dict(curve.operations),
              'total_wall_seconds': time.monotonic() - start,
              'total_cpu_seconds': time.process_time() - cpu,
              'peak_rss_bytes': rss()}
    args.out.write_text(json.dumps(report, sort_keys=True, separators=(',', ':')) + '\n')
    assert report['total_wall_seconds'] <= 120 and rss() <= 512 * 1024 * 1024


if __name__ == '__main__':
    main()
