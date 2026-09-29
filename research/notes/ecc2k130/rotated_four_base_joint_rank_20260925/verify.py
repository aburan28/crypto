#!/usr/bin/env python3
"""Independent bit-serial/Fermat replay of four-base rank and sealed holdouts."""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import resource
import sys
import tarfile
import time
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
NOTES = HERE.parent
ARITHMETIC = NOTES / 'rotated_subspace_support_20260925/verify.py'
CORPUS = NOTES / 'rotated_pdp_corpus_20260925/evidence/raw.tar.gz'
SWEEP = NOTES / 'rotated_beta_sweep_20260925/evidence/raw.tar.gz'
CORPUS_SHA = '39f16990213e2c7525c6dae92c131127b12f914bbfd6182264249da72bc9546c'
SWEEP_SHA = 'fe84aef6a2cf7f6f4c950245c9c8e870354fb750666f5482b81f9a997d107140'
DOMAIN = 'ECC2K130-ROTATED-JOINT-RANK-20260925-v1'
BETAS = (3, 338435, 303097, 464276)
SUPPORT = (62389, 66179, 66203, 59323)
Q = 130873
H = (385982, 301867)
LAM = 41811
TORSION = (None, (0, 1), (1, 0), (1, 1))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def point(value):
    return None if value is None else tuple(value)


def pjson(value):
    return None if value is None else list(value)


def arithmetic():
    spec = importlib.util.spec_from_file_location('joint_independent_arithmetic', ARITHMETIC)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def peak_rss():
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


def read(tar, name):
    member = tar.getmember(name)
    assert member.isfile() and member.size < 40_000_000
    handle = tar.extractfile(member)
    assert handle is not None
    raw = handle.read()
    assert len(raw) == member.size
    return raw.decode()


def gf2_rank(values):
    pivots = {}
    for value in values:
        while value:
            leading = value.bit_length() - 1
            if leading in pivots:
                value ^= pivots[leading]
            else:
                pivots[leading] = value
                break
    return len(pivots)


def load_archives(mod, curve):
    result = {}
    with tarfile.open(CORPUS, 'r:gz') as corpus, tarfile.open(SWEEP, 'r:gz') as sweep:
        for beta, known in zip(BETAS, SUPPORT):
            archive, prefix = ((corpus, 'raw/n19-m6/') if beta == 3
                               else (sweep, f'raw/beta-{beta}/'))
            # Derive all four factor sets directly, without producer code.
            f = curve.f
            conjugates, value = [], beta
            for _ in range(19):
                conjugates.append(value)
                value = f.square(value)
            assert value == beta and gf2_rank(conjugates) == 19 and f.trace(beta) == 1
            factors = []
            for i in range(6):
                a, b = conjugates[i], conjugates[6 + i]
                xs = (0, a, b, a ^ b)
                factors.append(sorted(p for x in xs for p in independent_lifts(curve, x)))
            archived_factors = [[point(p) for p in slot]
                                for slot in json.loads(read(archive, prefix + 'factors.json'))]
            assert factors == archived_factors and all(len(slot) == 7 for slot in factors)
            hist = {}
            tuple_total = 0
            for line in read(archive, prefix + 'projected_histogram.jsonl').splitlines():
                row = json.loads(line)
                target = point(row['point'])
                assert target not in hist and row['count'] > 0
                assert len(row['witness_indices']) == 6
                hist[target] = row
                tuple_total += row['count']
            assert len(hist) == known and tuple_total == 7 ** 6
            term_maps, universe = [], set()
            for i, slot in enumerate(factors):
                terms = []
                for p in slot:
                    back = p
                    for _ in range((19 - i) % 19):
                        back = curve.tau(back)
                    r = curve.scalar(back, 4)
                    if r is None:
                        assert p in TORSION
                        terms.append((None, 0))
                    else:
                        c = min(r, mod.neg(r) if hasattr(mod, 'neg') else (r[0], r[0] ^ r[1]))
                        coeff = (1 if c == r else -1) * pow(LAM, i, Q) % Q
                        assert curve.scalar(p, 4) == curve.scalar(c, coeff)
                        universe.add(c)
                        terms.append((c, coeff))
                term_maps.append(terms)
            base = {c for c, _ in term_maps[0] if c is not None}
            assert len(base) == 3 and universe == base
            result[beta] = {'factors': factors, 'hist': hist,
                            'terms': term_maps, 'columns': sorted(base)}
    return result


def independent_lifts(curve, x):
    f = curve.f
    if x == 0:
        return [(0, 1)]
    rhs = x ^ f.square(f.inv(x))
    if f.trace(rhs):
        return []
    z, term = 0, rhs
    for _ in range(10):
        z ^= term
        term = f.square(f.square(term))
    assert f.square(z) ^ z == rhs
    return sorted(((x, f.mul(x, z)), (x, f.mul(x, z ^ 1))))


def reconstruct(curve, archive_arm, beta, qpoint):
    r = curve.scalar(qpoint, 4)
    saved = archive_arm['hist'].get(r)
    if saved is None:
        return None
    indices = saved['witness_indices']
    assert len(indices) == 6 and all(0 <= j < 7 for j in indices)
    summation = None
    combined = {}
    for i, j in enumerate(indices):
        summation = curve.add(summation, archive_arm['factors'][i][j])
        col, coeff = archive_arm['terms'][i][j]
        if col is not None:
            combined[col] = (combined.get(col, 0) + coeff) % Q
    assert summation == point(saved['full_sum']) and curve.scalar(summation, 4) == r
    shift = curve.add(summation, None if qpoint is None else (qpoint[0], qpoint[0] ^ qpoint[1]))
    assert shift in TORSION and curve.add(qpoint, shift) == summation
    row = [{'column': list(col), 'coefficient': value}
           for col, value in sorted(combined.items()) if value]
    evaluated = None
    for term in row:
        evaluated = curve.add(evaluated, curve.scalar(tuple(term['column']), term['coefficient']))
    assert evaluated == r
    return {'beta': beta, 'point': pjson(qpoint), 'projected_point': pjson(r),
            'torsion_index': TORSION.index(shift), 'row': row,
            'source_witness_indices': indices, 'archived_full_sum': pjson(summation)}


def eliminate(records, columns):
    idx = {c: i for i, c in enumerate(columns)}
    basis = {}
    zero = dep = independent = 0
    first_full = None
    for position, record in enumerate(records, 1):
        vec = [0] * len(columns)
        for term in record['row']:
            j = idx[tuple(term['column'])]
            vec[j] = (vec[j] + term['coefficient']) % Q
        if not any(vec):
            zero += 1
        rhs = 4 * record['k'] % Q
        for lead in sorted(basis):
            if vec[lead]:
                pv, prhs = basis[lead]
                factor = vec[lead] * pow(pv[lead], -1, Q) % Q
                vec = [(a - factor * b) % Q for a, b in zip(vec, pv)]
                rhs = (rhs - factor * prhs) % Q
        lead = next((j for j, value in enumerate(vec) if value), None)
        if lead is None:
            assert rhs == 0
            dep += 1
        else:
            basis[lead] = (vec, rhs)
            independent += 1
            if len(basis) == len(columns) and first_full is None:
                first_full = position
    assert independent + dep == len(records)
    logs = None
    if len(basis) == len(columns):
        logs = [0] * len(columns)
        for j in reversed(range(len(columns))):
            vec, rhs = basis[j]
            logs[j] = (rhs - sum(vec[k] * logs[k] for k in range(j + 1, len(columns)))) * pow(vec[j], -1, Q) % Q
        for record in records:
            vec = [0] * len(columns)
            for term in record['row']:
                vec[idx[tuple(term['column'])]] = (vec[idx[tuple(term['column'])]] + term['coefficient']) % Q
            assert sum(a * b for a, b in zip(vec, logs)) % Q == 4 * record['k'] % Q
    return {'rows': len(records), 'zero_rows': zero, 'dependent_rows_including_zero': dep,
            'independent_rows': independent, 'rank': len(basis),
            'nullity': len(columns) - len(basis), 'first_full_rank_row': first_full,
            'full_rank': logs is not None}, logs


def replay(raw, inputs):
    started = time.monotonic()
    assert sha(CORPUS) == CORPUS_SHA and sha(SWEEP) == SWEEP_SHA
    train_input = json.loads((inputs / 'training.json').read_text())
    points = json.loads((inputs / 'point_only.json').read_text())
    sealed = json.loads((inputs / 'sealed_labels.json').read_text())
    summary = json.loads((raw / 'training' / 'training_summary.json').read_text())
    logs = json.loads((raw / 'training' / 'base_logs.json').read_text())
    rows = [json.loads(line) for line in (raw / 'training' / 'training_rows.jsonl').read_text().splitlines()]
    oracle = json.loads((raw / 'oracle.json').read_text())
    recovery = json.loads((raw / 'recovery.json').read_text())
    for clean in (points, logs, oracle):
        forbid_labels(clean)
    assert summary['input_sha256'] == sha(inputs / 'training.json')
    assert oracle['point_only_sha256'] == recovery['point_only_sha256'] == sha(inputs / 'point_only.json')
    assert recovery['base_logs_sha256'] == sha(raw / 'training' / 'base_logs.json')
    assert recovery['oracle_sha256'] == sha(raw / 'oracle.json')
    assert train_input['domain'] == points['domain'] == sealed['domain'] == DOMAIN
    assert train_input['q'] == points['q'] == logs['q'] == Q
    assert train_input['generator'] == points['generator'] == logs['generator'] == list(H)
    assert len(train_input['targets']) == 256 and len(points['targets']) == len(sealed['targets']) == 64
    ordered = sorted(range(1, Q), key=lambda k: (hashlib.sha256(f'{DOMAIN}/{k}'.encode()).digest(), k))[:320]
    mod = arithmetic()
    field = mod.GF(19, 0x80027)
    curve = mod.E(field)
    assert curve.on(H) and curve.scalar(H, Q) is None
    assert curve.tau(H) == curve.scalar(H, LAM)
    for i, case in enumerate(train_input['targets']):
        assert case == {'case_id': f'tr-{i:03d}', 'k': ordered[i],
                        'point': list(curve.scalar(H, ordered[i]))}
    for i, (case, label) in enumerate(zip(points['targets'], sealed['targets'])):
        assert case == {'case_id': f'ho-{i:03d}', 'point': list(curve.scalar(H, ordered[256+i]))}
        assert label == {'case_id': case['case_id'], 'k': ordered[256+i]}
    archive = load_archives(mod, curve)
    global_cols = sorted({c for beta in BETAS for c in archive[beta]['columns']})
    assert summary['global_columns'] == [list(c) for c in global_cols]
    assert summary['physical_columns_by_beta'] == {str(b): [list(c) for c in archive[b]['columns']] for b in BETAS}
    expected_rows = []
    misses = Counter()
    for case in train_input['targets']:
        for beta in BETAS:
            row = reconstruct(curve, archive[beta], beta, tuple(case['point']))
            if row is None:
                misses[str(beta)] += 1
            else:
                row['case_id'] = case['case_id']
                row['k'] = case['k']
                expected_rows.append(row)
    assert expected_rows == rows
    assert summary['training_misses_by_beta'] == dict(misses)
    assert summary['joint4']['rows'] + sum(misses.values()) == 1024
    base_rows = [row for row in rows if row['beta'] == 3]
    for name, records, columns, saved in (
        ('beta3', base_rows, archive[3]['columns'], summary['base3']),
        ('joint4', rows, global_cols, summary['joint4'])):
        derived, values = eliminate(records, columns)
        assert all(saved[k] == v for k, v in derived.items())
        claimed = next(arm for arm in logs['arms'] if arm['name'] == name)
        assert claimed['columns'] == [list(c) for c in columns]
        assert claimed['logs'] == values
        if values is not None:
            for col, k in zip(columns, values):
                assert curve.scalar(H, k) == col
    assert len(oracle['cases']) == len(recovery['cases']) == 64
    counts = Counter()
    recovered = Counter()
    for case, saved, answer, label in zip(points['targets'], oracle['cases'], recovery['cases'], sealed['targets']):
        assert case['case_id'] == saved['case_id'] == answer['case_id'] == label['case_id']
        assert case['point'] == saved['point'] == answer['point']
        qpoint = tuple(case['point'])
        for name, choices in (('beta3', BETAS[:1]), ('joint4', BETAS)):
            attempts, first = [], None
            for beta in choices:
                attempts.append(beta)
                first = reconstruct(curve, archive[beta], beta, qpoint)
                if first is not None:
                    break
            clean = None if first is None else {k: first[k] for k in
                ('beta', 'point', 'projected_point', 'torsion_index', 'row')}
            expected = {'case_id': case['case_id'], 'point': case['point'],
                        'attempted_betas': attempts, 'status': 'hit' if first else 'miss',
                        'response': clean}
            assert saved[name] == expected
            counts[f'{name}_{expected["status"]}'] += 1
            counts[f'{name}_base_probes'] += len(attempts)
            arm_logs = next(a for a in logs['arms'] if a['name'] == name)
            result = answer[name]
            if arm_logs['logs'] is None:
                assert result == {'status': 'rank_deficient', 'recovered_k': None}
            elif clean is None:
                assert result == {'status': 'archive_oracle_miss', 'recovered_k': None}
            else:
                by_col = {tuple(c): i for i, c in enumerate(arm_logs['columns'])}
                total = sum(term['coefficient'] * arm_logs['logs'][by_col[tuple(term['column'])]]
                            for term in clean['row']) % Q
                scalar = total * pow(4, -1, Q) % Q
                assert result == {'status': 'group_verified', 'recovered_k': scalar,
                                  'oracle_beta': clean['beta']}
                assert curve.scalar(H, scalar) == qpoint
                # The sealed scalar is opened only after the group-law check.
                assert scalar == label['k']
                recovered[name] += 1
    assert oracle['counts'] == dict(counts)
    assert recovered['joint4'] >= recovered['beta3']
    return {'status': 'PASS', 'source_archive_sha256': [CORPUS_SHA, SWEEP_SHA],
            'training_targets': 256, 'holdout_targets': 64,
            'beta3_rank': summary['base3']['rank'], 'beta3_columns': summary['base3']['column_count'],
            'joint_rank': summary['joint4']['rank'], 'joint_columns': summary['joint4']['column_count'],
            'beta3_recovered': recovered['beta3'], 'joint_recovered': recovered['joint4'],
            'full_rank_gate': summary['joint4']['full_rank'],
            'semantic_followup_gate': (summary['joint4']['full_rank'] and
                                       recovered['joint4'] - recovered['beta3'] >= 12),
            'wall_seconds': time.monotonic() - started,
            'peak_rss_bytes': peak_rss(),
            'field_operations': dict(field.ops), 'curve_operations': dict(curve.ops)}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--raw', type=Path, required=True)
    parser.add_argument('--inputs', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    result = replay(args.raw, args.inputs)
    args.out.write_text(json.dumps(result, sort_keys=True, separators=(',', ':')) + '\n')
    assert result['wall_seconds'] <= 300 and result['peak_rss_bytes'] <= 512 * 1024 * 1024


if __name__ == '__main__':
    main()
