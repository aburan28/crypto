#!/usr/bin/env python3
"""Cold first-witness archive lookup on frozen point-only Q targets."""
from __future__ import annotations

import argparse
import collections
import json
import resource
import sys
import time
from pathlib import Path

from core import BETAS, H, Q, Data, save, sha


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


def strip(row):
    if row is None:
        return None
    return {key: row[key] for key in ('beta', 'point', 'projected_point',
                                     'torsion_index', 'row')}


def response(data, beta_order, case):
    target = tuple(case['point'])
    assert data.curve.on_curve(target) and target is not None
    attempts = []
    chosen = None
    for beta in beta_order:
        attempts.append(beta)
        chosen = data.row(beta, target)
        if chosen is not None:
            break
    return {'case_id': case['case_id'], 'point': case['point'],
            'attempted_betas': attempts, 'status': 'hit' if chosen else 'miss',
            'response': strip(chosen)}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--points', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    assert not args.out.exists()
    start, cpu = time.monotonic(), time.process_time()
    points = json.loads(args.points.read_text())
    forbid_labels(points)
    assert points['schema'] == 'rotated_joint_point_only_v1' and points['q'] == Q
    assert points['generator'] == list(H) and len(points['targets']) == 64
    data = Data()
    cases = []
    counts = collections.Counter()
    for index, case in enumerate(points['targets']):
        assert case['case_id'] == f'ho-{index:03d}'
        ctrl = response(data, BETAS[:1], case)
        joint = response(data, BETAS, case)
        assert (ctrl['status'] != 'hit' or joint['status'] == 'hit')
        cases.append({'case_id': case['case_id'], 'point': case['point'],
                      'beta3': ctrl, 'joint4': joint})
        for name, output in (('beta3', ctrl), ('joint4', joint)):
            counts[f'{name}_{output["status"]}'] += 1
            counts[f'{name}_base_probes'] += len(output['attempted_betas'])
    result = {'schema': 'rotated_joint_archive_oracle_v1', 'point_only_sha256': sha(args.points),
              'beta_order': list(BETAS), 'cases': cases, 'counts': dict(counts),
              'archive_uncompressed_bytes_scanned': data.bytes_scanned,
              'field_operations': dict(data.field.operations),
              'curve_operations': dict(data.curve.operations),
              'wall_seconds': time.monotonic() - start, 'cpu_seconds': time.process_time() - cpu,
              'peak_rss_bytes': rss()}
    forbid_labels(result)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    save(args.out, result)
    assert result['wall_seconds'] <= 120 and result['peak_rss_bytes'] <= 512 * 1024 * 1024


if __name__ == '__main__':
    main()
