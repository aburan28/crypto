#!/usr/bin/env python3
"""Point-only recovery child: source arguments contain no scalar labels."""
from __future__ import annotations

import argparse
import collections
import json
import resource
import sys
import time
from pathlib import Path

from core import BETAS, H, Q, arithmetic, save, sha


def rss():
    raw = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return raw if sys.platform == 'darwin' else raw * 1024


def forbid_labels(value):
    if isinstance(value, dict):
        assert not {'k', 'source_witness_indices', 'archived_full_sum', 'witness_indices'} & set(value)
        for child in value.values():
            forbid_labels(child)
    elif isinstance(value, list):
        for child in value:
            forbid_labels(child)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--points', type=Path, required=True)
    parser.add_argument('--logs', type=Path, required=True)
    parser.add_argument('--oracle', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    assert not args.out.exists()
    start, cpu = time.monotonic(), time.process_time()
    points = json.loads(args.points.read_text())
    logs = json.loads(args.logs.read_text())
    oracle = json.loads(args.oracle.read_text())
    for value in (points, logs, oracle):
        forbid_labels(value)
    assert points['schema'] == 'rotated_joint_point_only_v1'
    assert logs['schema'] == 'rotated_joint_base_logs_v1'
    assert oracle['schema'] == 'rotated_joint_archive_oracle_v1'
    assert points['q'] == logs['q'] == Q
    assert points['generator'] == logs['generator'] == list(H)
    assert oracle['point_only_sha256'] == sha(args.points)
    assert oracle['beta_order'] == list(BETAS)
    assert [a['name'] for a in logs['arms']] == ['beta3', 'joint4']
    assert len(points['targets']) == len(oracle['cases']) == 64
    mod = arithmetic()
    field = mod.Field(19, [0, 1, 2, 5])
    curve = mod.Curve(field)
    assert curve.on_curve(H) and curve.scalar(H, Q) is None
    cases = []
    counts = collections.Counter()
    modular = collections.Counter()
    for target, oracle_case in zip(points['targets'], oracle['cases']):
        assert target['case_id'] == oracle_case['case_id']
        assert target['point'] == oracle_case['point']
        qpoint = tuple(target['point'])
        assert curve.on_curve(qpoint) and qpoint is not None
        case = {'case_id': target['case_id'], 'point': target['point']}
        for arm in logs['arms']:
            name = arm['name']
            lookup = oracle_case[name]
            assert lookup['case_id'] == target['case_id'] and lookup['point'] == target['point']
            columns = [tuple(c) for c in arm['columns']]
            assert columns == sorted(set(columns))
            by_col = {c: i for i, c in enumerate(columns)}
            base_logs = arm['logs']
            if base_logs is None:
                answer = {'status': 'rank_deficient', 'recovered_k': None}
            elif lookup['status'] == 'miss':
                assert lookup['response'] is None
                answer = {'status': 'archive_oracle_miss', 'recovered_k': None}
            else:
                assert lookup['status'] == 'hit'
                response = lookup['response']
                assert response['point'] == target['point']
                assert tuple(response['projected_point']) == curve.scalar(qpoint, 4)
                assert response['beta'] in (BETAS[:1] if name == 'beta3' else BETAS)
                assert 0 <= response['torsion_index'] < 4
                total = 0
                for term in response['row']:
                    j = by_col[tuple(term['column'])]
                    total = (total + term['coefficient'] * base_logs[j]) % Q
                    modular['mul'] += 1
                    modular['add'] += 1
                scalar = total * pow(4, -1, Q) % Q
                modular['mul'] += 1
                assert curve.scalar(H, scalar) == qpoint
                answer = {'status': 'group_verified', 'recovered_k': scalar,
                          'oracle_beta': response['beta']}
            case[name] = answer
            counts[f'{name}_{answer["status"]}'] += 1
        cases.append(case)
    result = {'schema': 'rotated_joint_point_recovery_v1', 'cases': cases,
              'point_only_sha256': sha(args.points), 'base_logs_sha256': sha(args.logs),
              'oracle_sha256': sha(args.oracle), 'counts': dict(counts),
              'modular_operations': dict(modular),
              'field_operations': dict(field.operations),
              'curve_operations': dict(curve.operations),
              'wall_seconds': time.monotonic() - start,
              'cpu_seconds': time.process_time() - cpu, 'peak_rss_bytes': rss()}
    save(args.out, result)
    assert result['wall_seconds'] <= 120 and result['peak_rss_bytes'] <= 512 * 1024 * 1024


if __name__ == '__main__':
    main()
