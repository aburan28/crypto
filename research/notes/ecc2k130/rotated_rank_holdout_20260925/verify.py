#!/usr/bin/env python3
"""Independent bit-serial/Fermat replay of rank, holdout and point-only logs."""
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
OLD_VERIFY = HERE.parent / 'rotated_subspace_support_20260925/verify.py'
DOMAIN = 'ECC2K130-ROTATED-RANK-20260925-v1'
ARMS = ((3, 5), (3, 6), (7, 5), (7, 6))
Q = 2003


def prior():
    spec = importlib.util.spec_from_file_location('rank_independent_verify', OLD_VERIFY)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def member(tar, name):
    x = tar.getmember('raw/' + name)
    assert x.isfile() and x.size < 50_000_000
    handle = tar.extractfile(x)
    assert handle is not None
    return handle.read().decode()


def point(raw):
    return None if raw is None else tuple(raw)


def forbidden_keys(value):
    if isinstance(value, dict):
        assert not {'k', 'source_indices', 'source_points', 'witness_indices'} & set(value)
        for x in value.values():
            forbidden_keys(x)
    elif isinstance(value, list):
        for x in value:
            forbidden_keys(x)


def scan_selection(tar, manifest, point_only, sealed, curve, multiples):
    assert manifest['domain'] == DOMAIN and manifest['q'] == Q
    assert manifest['source_archive_sha256'] == SHA
    assert manifest['generator'] == [4793, 2429]
    assert len(manifest['arms']) == len(point_only['arms']) == len(sealed['arms']) == 4
    for (beta, m), a, p, lab in zip(ARMS, manifest['arms'], point_only['arms'], sealed['arms']):
        assert (a['beta'], a['m']) == (p['beta'], p['m']) == (lab['beta'], lab['m']) == (beta, m)
        roster = [json.loads(line) for line in member(tar, f'n13-b{beta}-m{m}-roster.jsonl').splitlines()]
        assert len(roster) == Q and [r['k'] for r in roster] == list(range(Q))
        expected_scanned, selected = [], []
        for k in sorted(range(1, Q),
                        key=lambda x: (hashlib.sha256(f'{DOMAIN}/{beta}/{m}/{x}'.encode()).digest(), x)):
            record = roster[k]
            assert point(record['point']) == multiples[k]
            cosets = record['coset_witness_indices']
            positive = any(x is not None for x in cosets)
            expected_scanned.append({'k': k,
                'score_sha256': hashlib.sha256(f'{DOMAIN}/{beta}/{m}/{k}'.encode()).hexdigest(),
                'positive': positive, 'positive_cosets': sum(x is not None for x in cosets)})
            if positive:
                selected.append((k, record['point']))
                if len(selected) == 16:
                    break
        assert a['scanned_candidates'] == expected_scanned
        assert a['selection_positive_count'] == 16
        assert a['selection_negative_count'] == sum(not x['positive'] for x in expected_scanned)
        assert a['selected_case_ids'] == [f'b{beta}m{m}-{j:02d}' for j in range(16)]
        assert len(p['cases']) == len(lab['labels']) == 16
        for j, ((k, qpoint), case, label) in enumerate(zip(selected, p['cases'], lab['labels'])):
            assert case == {'case_id': f'b{beta}m{m}-{j:02d}', 'point': qpoint}
            assert label == {'case_id': case['case_id'], 'k': k}
            assert curve.on(point(case['point']))
        assert len({tuple(c['point']) for c in p['cases']}) == 16


def coefficient_vector(row, cols):
    index = {c: j for j, c in enumerate(cols)}
    result = [0] * len(cols)
    for item in row['row']:
        j = index[tuple(item['column'])]
        result[j] = (result[j] + item['coefficient']) % Q
    return result


def replay_arm(tar, beta, m, point_arm, training, logs, oracle, recovery,
               sealed, curve, multiples):
    stem = f'n13-b{beta}-m{m}'
    source_rows = [json.loads(line) for line in member(tar, stem + '-rows.jsonl').splitlines()]
    meta = json.loads(member(tar, stem + '-summary.json'))
    columns = [tuple(c) for c in meta['term_column_coordinates']]
    assert columns == sorted(set(columns))
    assert logs['columns'] == [list(c) for c in columns]
    assert training['columns'] == [list(c) for c in columns]
    excluded = {tuple(c['point']) for c in point_arm['cases']}
    assert len(excluded) == 16
    count_excluded = count_kept = count_zero = dependent = independent = 0
    first_full = None
    basis = {}
    for row in source_rows:
        assert 0 <= row['k'] < Q
        target = point(row['target_Q'])
        assert target == multiples[row['k']]
        vec = coefficient_vector(row, columns)
        rhs = (4 * row['k']) % Q
        if target in excluded:
            count_excluded += 1
            continue
        count_kept += 1
        if not any(vec):
            count_zero += 1
        # Separate unnormalised elimination; the trainer stores normalised pivots.
        for col in sorted(basis):
            if vec[col]:
                pv, prhs = basis[col]
                factor = vec[col] * pow(pv[col], -1, Q) % Q
                vec = [(x - factor * y) % Q for x, y in zip(vec, pv)]
                rhs = (rhs - factor * prhs) % Q
        pivot = next((j for j, x in enumerate(vec) if x), None)
        if pivot is None:
            assert rhs == 0
            dependent += 1
        else:
            basis[pivot] = (vec, rhs)
            independent += 1
            if len(basis) == len(columns) and first_full is None:
                first_full = count_kept
    assert training['total_rows'] == len(source_rows) == meta['archive_positive_rows']
    assert training['holdout_points'] == 16
    assert training['excluded_rows'] == count_excluded
    assert training['considered_rows'] == count_kept
    assert training['zero_rows'] == count_zero
    assert training['dependent_rows_including_zero'] == dependent
    assert training['independent_rows'] == independent
    assert training['rank'] == len(basis)
    assert training['nullity'] == len(columns) - len(basis)
    assert training['first_full_rank_considered_row'] == first_full
    assert training['full_rank'] == (len(basis) == len(columns))
    if len(basis) == len(columns):
        assert logs['logs'] is not None
        # Check the claimed solution against every retained equation.
        base_logs = logs['logs']
        assert len(base_logs) == len(columns)
        for c, scalar in zip(columns, base_logs):
            assert 0 <= scalar < Q and curve.scalar((4793, 2429), scalar) == c
        for row in source_rows:
            if point(row['target_Q']) in excluded:
                continue
            original_vector = coefficient_vector(row, columns)
            assert sum(x * y for x, y in zip(original_vector, base_logs)) % Q == 4 * row['k'] % Q
    else:
        assert logs['logs'] is None
    # Full-archive oracle scan, first row in frozen (k,T) order, no saved k copied.
    first = {}
    for row in source_rows:
        target = point(row['target_Q'])
        if target in excluded and target not in first:
            first[target] = {'point': list(target), 'torsion_index': row['torsion_index'],
                             'row': row['row']}
    assert len(oracle['responses']) == len(recovery['cases']) == len(sealed['labels']) == 16
    hit = recovered = missed = skipped = 0
    for case, oracle_response, result, label in zip(point_arm['cases'], oracle['responses'],
                                                     recovery['cases'], sealed['labels']):
        assert case['case_id'] == oracle_response['case_id'] == result['case_id'] == label['case_id']
        assert case['point'] == oracle_response['point'] == result['point']
        target = tuple(case['point'])
        expected_response = first.get(target)
        assert oracle_response['response'] == expected_response
        assert oracle_response['lookup_status'] == ('hit' if expected_response else 'miss')
        if expected_response is not None:
            hit += 1
        if logs['logs'] is None:
            assert result['status'] == 'rank_deficient' and result['recovered_k'] is None
            skipped += 1
        elif expected_response is None:
            assert result['status'] == 'oracle_miss' and result['recovered_k'] is None
            missed += 1
        else:
            assert result['status'] == 'recovered_and_group_verified'
            assert result['oracle_torsion_index'] == expected_response['torsion_index']
            vector = coefficient_vector(expected_response, columns)
            computed = (sum(x * y for x, y in zip(vector, logs['logs'])) * pow(4, -1, Q)) % Q
            assert result['recovered_k'] == computed
            assert curve.scalar((4793, 2429), computed) == target
            # Sealed scalar comparison only after independent group-law validation.
            assert computed == label['k']
            recovered += 1
    assert oracle['beta'] == recovery['beta'] == training['beta'] == beta
    assert oracle['m'] == recovery['m'] == training['m'] == m
    return {'beta': beta, 'm': m, 'columns': len(columns), 'rank': len(basis),
            'nullity': len(columns) - len(basis), 'rows': len(source_rows),
            'excluded_rows': count_excluded, 'considered_rows': count_kept,
            'zero_rows': count_zero, 'dependent_rows_including_zero': dependent,
            'first_full_rank_considered_row': first_full,
            'holdout_oracle_hits': hit, 'holdout_recovered': recovered,
            'holdout_misses': missed, 'rank_deficient_skips': skipped}


def replay(raw: Path, inputs: Path):
    wall, cpu = time.monotonic(), time.process_time()
    assert digest(ARCHIVE) == SHA
    manifest = json.loads((inputs / 'manifest.json').read_text())
    point_only = json.loads((inputs / 'point_only.json').read_text())
    sealed = json.loads((inputs / 'sealed_labels.json').read_text())
    train = json.loads((raw / 'training' / 'training_summary.json').read_text())
    logs = json.loads((raw / 'training' / 'base_logs.json').read_text())
    oracle = json.loads((raw / 'oracle_responses.json').read_text())
    recovered = json.loads((raw / 'recovery.json').read_text())
    forbidden_keys(point_only)
    forbidden_keys(logs)
    forbidden_keys(oracle)
    assert logs['q'] == Q and logs['generator'] == [4793, 2429]
    assert oracle['source_archive_sha256'] == SHA
    assert oracle['point_only_sha256'] == digest(inputs / 'point_only.json')
    assert train['source_archive_sha256'] == SHA
    assert train['point_only_sha256'] == digest(inputs / 'point_only.json')
    assert len(train['arms']) == len(logs['arms']) == len(oracle['arms']) == len(recovered['arms']) == 4
    mod = prior()
    field = mod.GF(13, mod.P13)
    curve = mod.E(field)
    assert mod.source_order_by_recurrence(13) == 8012
    H = (4793, 2429)
    assert curve.on(H) and curve.scalar(H, Q) is None
    multiples, current = [], None
    for _ in range(Q):
        multiples.append(current)
        current = curve.add(current, H)
    assert current is None
    with tarfile.open(ARCHIVE, 'r:gz') as tar:
        scan_selection(tar, manifest, point_only, sealed, curve, multiples)
        reports = [replay_arm(tar, b, m, point_only['arms'][i], train['arms'][i],
                              logs['arms'][i], oracle['arms'][i], recovered['arms'][i],
                              sealed['arms'][i], curve, multiples)
                   for i, (b, m) in enumerate(ARMS)]
    assert sum(x['archive_rows_scanned'] for x in oracle['scan_cost']['arms']) == sum(r['rows'] for r in reports)
    assert [(x['beta'], x['m']) for x in oracle['scan_cost']['arms']] == list(ARMS)
    report = {'status': 'success', 'source_archive_sha256': SHA,
              'input_hashes': {p.name: digest(p) for p in sorted(inputs.glob('*.json'))},
              'arms': reports, 'total_wall_seconds': time.monotonic() - wall,
              'total_cpu_seconds': time.process_time() - cpu,
              'peak_rss_bytes': mod.peak_rss_bytes(),
              'field_operations': dict(field.ops), 'curve_operations': dict(curve.ops)}
    assert report['total_wall_seconds'] <= 300 and report['peak_rss_bytes'] <= 512 * 1024 * 1024
    return report


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--raw', required=True, type=Path)
    parser.add_argument('--inputs', required=True, type=Path)
    parser.add_argument('--report', required=True, type=Path)
    args = parser.parse_args()
    args.report.write_text(json.dumps(replay(args.raw, args.inputs),
                                      sort_keys=True, separators=(',', ':')) + '\n')


if __name__ == '__main__':
    main()
