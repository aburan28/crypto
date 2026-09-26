#!/usr/bin/env python3
"""Exact modular training from frozen, unconditional n19 target stream."""
from __future__ import annotations

import argparse
import collections
import json
import resource
import sys
import time
from pathlib import Path

from core import BETAS, H, Q, Data, save, sha, vector


def rss():
    raw = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return raw if sys.platform == 'darwin' else raw * 1024


def solve(rows, columns, name, curve):
    pivots = {}
    modular = collections.Counter()
    zero = dependent = independent = 0
    first_full = None
    for position, record in enumerate(rows, 1):
        values = vector(record, columns)
        rhs = (4 * record['k']) % Q
        if not any(values):
            zero += 1
        for col in sorted(pivots):
            factor = values[col]
            if factor:
                pv, prhs = pivots[col]
                values = [(a - factor * b) % Q for a, b in zip(values, pv)]
                rhs = (rhs - factor * prhs) % Q
                modular['mul'] += len(values) + 1
                modular['add'] += len(values) + 1
        pivot = next((j for j, value in enumerate(values) if value), None)
        if pivot is None:
            assert rhs == 0, f'{name}: inconsistent row {position}'
            dependent += 1
        else:
            inverse = pow(values[pivot], -1, Q)
            modular['inverse'] += 1
            pivots[pivot] = ([v * inverse % Q for v in values], rhs * inverse % Q)
            modular['mul'] += len(values) + 1
            independent += 1
            if len(pivots) == len(columns) and first_full is None:
                first_full = position
    assert independent + dependent == len(rows) and independent == len(pivots)
    logs = None
    if len(pivots) == len(columns):
        logs = [0] * len(columns)
        for col in reversed(range(len(columns))):
            pv, rhs = pivots[col]
            logs[col] = (rhs - sum(pv[j] * logs[j] for j in range(col + 1, len(columns)))) % Q
            modular['mul'] += len(columns) - col - 1
            modular['add'] += len(columns) - col - 1
        for column, log in zip(columns, logs):
            assert curve.scalar(H, log) == column
        for record in rows:
            assert sum(a * b for a, b in zip(vector(record, columns), logs)) % Q == 4 * record['k'] % Q
            modular['training_equations_checked'] += 1
    return ({'arm': name, 'columns': [list(c) for c in columns], 'column_count': len(columns),
             'rows': len(rows), 'zero_rows': zero, 'dependent_rows_including_zero': dependent,
             'independent_rows': independent, 'rank': len(pivots),
             'nullity': len(columns) - len(pivots),
             'first_full_rank_row': first_full, 'full_rank': logs is not None,
             'modular_operations': dict(modular)}, logs)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--training', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    assert not args.out.exists()
    args.out.mkdir(parents=True)
    start, cpu = time.monotonic(), time.process_time()
    data = Data()
    training = json.loads(args.training.read_text())
    assert training['schema'] == 'rotated_joint_training_v1' and training['q'] == Q
    assert training['generator'] == list(H) and len(training['targets']) == 256
    assert len({row['k'] for row in training['targets']}) == 256
    rows = []
    misses = collections.Counter()
    for pos, target in enumerate(training['targets']):
        assert target['case_id'] == f'tr-{pos:03d}' and 0 < target['k'] < Q
        qpoint = tuple(target['point'])
        assert data.curve.scalar(H, target['k']) == qpoint
        for beta in BETAS:
            row = data.row(beta, qpoint)
            if row is None:
                misses[str(beta)] += 1
                continue
            row['case_id'] = target['case_id']
            row['k'] = target['k']
            rows.append(row)
    base_rows = [r for r in rows if r['beta'] == BETAS[0]]
    base_cols = data.arms[BETAS[0]]['columns']
    base, base_logs = solve(base_rows, base_cols, 'beta3', data.curve)
    joint, joint_logs = solve(rows, data.columns, 'joint4', data.curve)
    assert base['rows'] + misses[str(BETAS[0])] == 256
    assert joint['rows'] + sum(misses.values()) == 256 * len(BETAS)
    overlap = {}
    for i, left in enumerate(BETAS):
        for right in BETAS[i+1:]:
            overlap[f'{left},{right}'] = len(set(data.arms[left]['columns']) & set(data.arms[right]['columns']))
    with (args.out / 'training_rows.jsonl').open('w') as out:
        for row in rows:
            out.write(json.dumps(row, sort_keys=True, separators=(',', ':')) + '\n')
    save(args.out / 'base_logs.json', {'schema': 'rotated_joint_base_logs_v1',
         'q': Q, 'generator': list(H),
         'arms': [{'name': 'beta3', 'columns': base['columns'], 'logs': base_logs},
                  {'name': 'joint4', 'columns': joint['columns'], 'logs': joint_logs}]})
    result = {'schema': 'rotated_joint_training_summary_v1',
              'input_sha256': sha(args.training), 'archive_sha256': [data_sha for data_sha in
              ('39f16990213e2c7525c6dae92c131127b12f914bbfd6182264249da72bc9546c',
               'fe84aef6a2cf7f6f4c950245c9c8e870354fb750666f5482b81f9a997d107140')],
              'beta_order': list(BETAS), 'training_targets': 256,
              'training_misses_by_beta': dict(misses),
              'physical_columns_by_beta': {str(beta): [list(c) for c in data.arms[beta]['columns']]
                                           for beta in BETAS},
              'cross_base_canonical_column_overlaps': overlap,
              'global_columns': [list(c) for c in data.columns],
              'base3': base, 'joint4': joint,
              'archive_uncompressed_bytes_scanned': data.bytes_scanned,
              'field_operations': dict(data.field.operations),
              'curve_operations': dict(data.curve.operations),
              'wall_seconds': time.monotonic() - start,
              'cpu_seconds': time.process_time() - cpu,
              'peak_rss_bytes': rss()}
    save(args.out / 'training_summary.json', result)
    assert result['wall_seconds'] <= 120 and result['peak_rss_bytes'] <= 512 * 1024 * 1024


if __name__ == '__main__':
    main()
