#!/usr/bin/env python3
"""Train modular column logs with all SHA-frozen point holdouts excluded."""
from __future__ import annotations

import argparse
import collections
import hashlib
import importlib.util
import json
import resource
import sys
import tarfile
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
PRIOR = HERE.parent / 'rotated_row_certificate_20260925'
ARCHIVE = PRIOR / 'evidence/raw.tar.gz'
SHA = '863f50ea872ce01b1999d0ef67bc0341fa0c764f2ad4187d16736768b0ab9b3d'
OLD_GATE = HERE.parent / 'rotated_subspace_support_20260925/gate.py'
Q = 2003
ARMS = ((3, 5), (3, 6), (7, 5), (7, 6))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def old():
    spec = importlib.util.spec_from_file_location('rank_prior_gate', OLD_GATE)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def member(tar, name):
    x = tar.getmember('raw/' + name)
    assert x.isfile() and x.size < 50_000_000
    handle = tar.extractfile(x)
    assert handle is not None
    return handle.read().decode()


def rss():
    v = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return v if sys.platform == 'darwin' else v * 1024


def save(path, obj):
    path.write_text(json.dumps(obj, sort_keys=True, separators=(',', ':')) + '\n')


def train_arm(tar, mod, field, curve, beta, m, cases, counters):
    start, cpu = time.monotonic(), time.process_time()
    stem = f'n13-b{beta}-m{m}'
    meta = json.loads(member(tar, stem + '-summary.json'))
    columns = [tuple(c) for c in meta['term_column_coordinates']]
    assert columns == sorted(set(columns)) and len(columns) in (2, 3)
    column_index = {column: i for i, column in enumerate(columns)}
    point_exclude = {tuple(case['point']) for case in cases}
    assert len(point_exclude) == len(cases) == 16
    rows = [json.loads(line) for line in member(tar, stem + '-rows.jsonl').splitlines()]
    assert len(rows) == meta['archive_positive_rows']
    H = (4793, 2429)
    walk = None
    multiples = []
    for _ in range(Q):
        multiples.append(walk)
        walk = curve.add(walk, H)
    assert walk is None
    pivots = {}
    first_full = None
    zero_rows = dependent = independent = excluded = considered = 0
    for row in rows:
        counters['rows_scanned'] += 1
        k = row['k']
        assert 0 <= k < Q
        target = None if row['target_Q'] is None else tuple(row['target_Q'])
        assert target == multiples[k]
        if target in point_exclude:
            excluded += 1
            continue
        considered += 1
        vector = [0] * len(columns)
        for item in row['row']:
            index = column_index[tuple(item['column'])]
            vector[index] = (vector[index] + item['coefficient']) % Q
            counters['modular_add'] += 1
        rhs = (4 * k) % Q
        if not any(vector):
            zero_rows += 1
        for col in sorted(pivots):
            factor = vector[col]
            if factor:
                pv, prhs = pivots[col]
                vector = [(a - factor * b) % Q for a, b in zip(vector, pv)]
                rhs = (rhs - factor * prhs) % Q
                counters['modular_mul'] += len(vector) + 1
                counters['modular_add'] += len(vector) + 1
        pivot = next((j for j, v in enumerate(vector) if v), None)
        if pivot is None:
            assert rhs == 0
            dependent += 1
        else:
            inverse = pow(vector[pivot], -1, Q)
            counters['modular_inverse'] += 1
            pivots[pivot] = ([(v * inverse) % Q for v in vector], (rhs * inverse) % Q)
            counters['modular_mul'] += len(vector) + 1
            independent += 1
            if len(pivots) == len(columns) and first_full is None:
                first_full = considered
    assert independent == len(pivots) and considered == independent + dependent
    logs = None
    if len(pivots) == len(columns):
        logs = [0] * len(columns)
        for col in reversed(range(len(columns))):
            pv, rhs = pivots[col]
            logs[col] = (rhs - sum(pv[j] * logs[j] for j in range(col + 1, len(columns)))) % Q
            counters['modular_mul'] += len(columns) - col - 1
        for row in rows:
            target = None if row['target_Q'] is None else tuple(row['target_Q'])
            if target in point_exclude:
                continue
            lhs = sum(item['coefficient'] * logs[column_index[tuple(item['column'])]]
                      for item in row['row']) % Q
            assert lhs == (4 * row['k']) % Q
            counters['training_equations_checked'] += 1
    total = {'wall_seconds': time.monotonic() - start,
             'cpu_seconds': time.process_time() - cpu,
             'peak_rss_bytes': rss()}
    return ({'beta': beta, 'm': m, 'total_rows': len(rows), 'holdout_points': len(point_exclude),
             'excluded_rows': excluded, 'considered_rows': considered, 'zero_rows': zero_rows,
             'dependent_rows_including_zero': dependent, 'independent_rows': independent,
             'rank': len(pivots), 'nullity': len(columns) - len(pivots),
             'first_full_rank_considered_row': first_full, 'columns': [list(c) for c in columns],
             'full_rank': logs is not None, 'cost': total},
            {'beta': beta, 'm': m, 'columns': [list(c) for c in columns], 'logs': logs})


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--points', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    assert digest(ARCHIVE) == SHA
    points = json.loads(args.points.read_text())
    assert [(a['beta'], a['m']) for a in points['arms']] == list(ARMS)
    start, cpu = time.monotonic(), time.process_time()
    mod = old()
    field = mod.Field(13, mod.MODELS[13]['low'])
    curve = mod.Curve(field)
    assert mod.source_group_order(13) == 4 * Q
    counters = collections.Counter()
    summaries, logs = [], []
    with tarfile.open(ARCHIVE, 'r:gz') as tar:
        for (beta, m), arm in zip(ARMS, points['arms']):
            summary, base_logs = train_arm(tar, mod, field, curve, beta, m,
                                           arm['cases'], counters)
            summaries.append(summary)
            logs.append(base_logs)
    args.out.mkdir(parents=True, exist_ok=False)
    save(args.out / 'training_summary.json',
         {'source_archive_sha256': SHA, 'point_only_sha256': digest(args.points),
          'arms': summaries, 'modular_counters': dict(counters),
          'field_operations': dict(field.operations), 'curve_operations': dict(curve.operations),
          'total_wall_seconds': time.monotonic() - start,
          'total_cpu_seconds': time.process_time() - cpu, 'peak_rss_bytes': rss()})
    save(args.out / 'base_logs.json', {'schema': 'rotated_rank_base_logs_v1',
                                      'q': Q, 'generator': [4793, 2429], 'arms': logs})
    assert time.monotonic() - start <= 120 and rss() <= 512 * 1024 * 1024


if __name__ == '__main__':
    main()
