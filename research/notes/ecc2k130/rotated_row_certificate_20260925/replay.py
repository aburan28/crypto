#!/usr/bin/env python3
"""Independent full-point replay of every archived compressed relation row."""
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
PRIOR = HERE.parent / 'rotated_subspace_support_20260925'
OLD_ARCHIVE = PRIOR / 'evidence/raw.tar.gz'
OLD_SHA = 'fad4b4cc48f980349d1ef472d18327368dd0f93f3d1284fc8bbdc20f0d3086d7'
ARMS = ((3, 5, 1799), (3, 6, 6307), (7, 5, 1799), (7, 6, 8012))
DOMAIN = 'ECC2K130-ROTATED-ROW-20260925-v1'


def prior():
    spec = importlib.util.spec_from_file_location('independent_old_verify', PRIOR / 'verify.py')
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def point(value):
    return None if value is None else tuple(value)


def neg(p):
    return None if p is None else (p[0], p[0] ^ p[1])


def many(curve, items):
    answer = None
    for item in items:
        answer = curve.add(answer, item)
    return answer


def json_lines(path):
    return [json.loads(line) for line in path.read_text().splitlines()]


def member(tar, name):
    item = tar.getmember('raw/toy/' + name)
    assert item.isfile() and item.size < 20_000_000
    data = tar.extractfile(item)
    assert data is not None
    return data.read().decode()


def torsion_check(curve):
    torsion = [None, (0, 1), (1, 0), (1, 1)]
    assert all(curve.on(t) for t in torsion)
    assert curve.scalar(torsion[1], 2) is None
    assert curve.scalar(torsion[2], 2) == torsion[1]
    assert curve.scalar(torsion[3], 2) == torsion[1]
    return torsion


def independent_row(curve, n, q, lam, source, target, ti, torsion, indices):
    assert all(curve.on(p) for p in source) and curve.on(target)
    raw = many(curve, source)
    assert raw == curve.add(target, torsion[ti])
    terms = []
    combined = collections.defaultdict(int)
    for i, p in enumerate(source):
        transported_back = p
        for _ in range((n - i) % n):
            transported_back = curve.tau(transported_back)
        projected = curve.scalar(transported_back, 4)
        scalar = pow(lam, i, q)
        assert curve.scalar(p, 4) == curve.scalar(projected, scalar)
        if projected is None:
            assert p in torsion and transported_back in torsion
            terms.append({'slot': i, 'source_index': None if indices is None else indices[i],
                          'source_point': None if p is None else list(p),
                          'back_point': None if transported_back is None else list(transported_back),
                          'projected_point': None, 'canonical_column': None,
                          'sign': 0, 'coefficient': 0,
                          'torsion_index': torsion.index(p)})
        else:
            twin = neg(projected)
            canonical = min(projected, twin)
            sign = +1 if canonical == projected else -1
            coefficient = (sign * scalar) % q
            combined[canonical] = (combined[canonical] + coefficient) % q
            terms.append({'slot': i, 'source_index': None if indices is None else indices[i],
                          'source_point': list(p), 'back_point': list(transported_back),
                          'projected_point': list(projected), 'canonical_column': list(canonical),
                          'sign': sign, 'coefficient': coefficient, 'torsion_index': None})
    counts = collections.Counter(tuple(t['canonical_column']) for t in terms
                                 if t['canonical_column'] is not None)
    row = [{'column': list(column), 'coefficient': coefficient}
           for column, coefficient in sorted(combined.items()) if coefficient]
    result = many(curve, (curve.scalar(tuple(item['column']), item['coefficient']) for item in row))
    assert result == curve.scalar(target, 4)
    assert result == many(curve, (curve.scalar(p, 4) for p in source))
    return {'source_indices': indices,
            'source_points': [None if p is None else list(p) for p in source],
            'target_Q': None if target is None else list(target), 'torsion_index': ti,
            'raw_sum': None if raw is None else list(raw), 'terms': terms,
            'row': row, 'row_point': None if result is None else list(result),
            'weight': len(row), 'zero_terms': sum(t['projected_point'] is None for t in terms),
            'repeated_terms': sum(count - 1 for count in counts.values() if count > 1),
            'cancelled_columns': sum(1 for col, count in counts.items()
                                     if count > 1 and combined[col] == 0)}


def check_arm(mod, tar, raw: Path, curve, by_x, beta: int, m: int, expected: int, H, torsion):
    begin = time.monotonic()
    stem = f'n13-b{beta}-m{m}'
    archive_stem = stem + '-rotated-'
    archive_factors = json.loads(member(tar, archive_stem + 'factors.json'))
    factors = mod.factor_points(curve, by_x, beta, m, False)
    assert archive_factors == [[list(p) for p in factor] for factor in factors]
    archive_summary = json.loads(member(tar, archive_stem + 'summary.json'))
    archive_targets = [json.loads(line) for line in member(tar, archive_stem + 'targets.jsonl').splitlines()]
    roster = json_lines(raw / f'{stem}-roster.jsonl')
    rows = json_lines(raw / f'{stem}-rows.jsonl')
    assert len(roster) == len(archive_targets) == 2003
    assert len(rows) == expected == archive_summary['distinct_full_sums']
    row_iter = iter(rows)
    positive = 0
    coset_hits = [0] * 4
    columns = set()
    term_columns = set()
    weights = collections.Counter()
    zero_terms = repeats = cancellations = 0
    current = None
    for k in range(2003):
        source_record = archive_targets[k]
        frozen_roster = roster[k]
        assert source_record['k'] == frozen_roster['k'] == k
        assert point(source_record['point']) == point(frozen_roster['point']) == current
        assert point(source_record['projected_point']) == curve.scalar(current, 4)
        assert source_record['coset_multiplicities'] == frozen_roster['coset_multiplicities']
        assert source_record['coset_witness_indices'] == frozen_roster['coset_witness_indices']
        for ti, indices in enumerate(source_record['coset_witness_indices']):
            if indices is None:
                assert source_record['coset_multiplicities'][ti] == 0
                continue
            assert source_record['coset_multiplicities'][ti] > 0
            got = next(row_iter)
            source = [factors[i][j] for i, j in enumerate(indices)]
            replayed = independent_row(curve, 13, 2003, 89, source, current, ti, torsion, indices)
            replayed.update({'beta': beta, 'm': m, 'k': k, 'status': 'archive_witness'})
            assert got == replayed
            positive += 1
            coset_hits[ti] += 1
            weights[got['weight']] += 1
            zero_terms += got['zero_terms']
            repeats += got['repeated_terms']
            cancellations += got['cancelled_columns']
            columns.update(tuple(item['column']) for item in got['row'])
            term_columns.update(tuple(item['canonical_column']) for item in got['terms']
                                if item['canonical_column'] is not None)
        current = curve.add(current, H)
    assert current is None and positive == expected
    try:
        next(row_iter)
        raise AssertionError('extra row after frozen roster')
    except StopIteration:
        pass
    assert coset_hits == archive_summary['coset_hits']
    saved_controls = json.loads((raw / f'{stem}-controls.json').read_text())
    p0 = next(p for p in factors[0] if curve.scalar(p, 4) is not None)
    assert len(saved_controls) == 2
    for name, source in (('all_x0', [(0, 1)] * m),
                         ('repeated_signed', [p0, curve.tau(neg(p0))] + [(0, 1)] * (m - 2))):
        raw_sum = many(curve, source)
        target = curve.scalar(curve.scalar(raw_sum, 4), pow(4, -1, 2003))
        tor = curve.add(raw_sum, neg(target))
        assert tor in torsion
        replayed = independent_row(curve, 13, 2003, 89, source,
                                   target, torsion.index(tor), torsion, None)
        replayed.update({'beta': beta, 'm': m, 'control': name})
        got = saved_controls.pop(0)
        assert got == replayed
        if name == 'all_x0':
            assert got['weight'] == 0 and got['zero_terms'] == m
        else:
            assert got['repeated_terms'] == 1 and got['zero_terms'] == m - 2
            assert len(got['terms']) == m
            assert (got['terms'][0]['coefficient'] + got['terms'][1]['coefficient']) % 2003 != 0
    saved = json.loads((raw / f'{stem}-summary.json').read_text())
    assert (saved['beta'], saved['m'], saved['n']) == (beta, m, 13)
    assert saved['archive_positive_rows'] == expected and saved['coset_hits'] == coset_hits
    assert saved['distinct_columns'] == len(columns)
    assert saved['column_coordinates'] == [list(c) for c in sorted(columns)]
    assert saved['distinct_term_columns'] == len(term_columns)
    assert saved['term_column_coordinates'] == [list(c) for c in sorted(term_columns)]
    assert {int(k): v for k, v in saved['row_weights'].items()} == dict(weights)
    assert saved['zero_terms'] == zero_terms and saved['repeated_terms'] == repeats
    assert saved['cancelled_columns'] == cancellations
    assert saved['total']['wall_seconds'] <= 120 and saved['total']['peak_rss_bytes'] <= 512 * 1024 * 1024
    return {'beta': beta, 'm': m, 'positive_rows': positive,
            'coset_hits': coset_hits, 'distinct_columns': len(columns),
            'distinct_term_columns': len(term_columns),
            'row_weights': dict(sorted(weights.items())), 'zero_terms': zero_terms,
            'repeated_terms': repeats, 'cancelled_columns': cancellations,
            'wall_seconds': time.monotonic() - begin}


def rank(vectors):
    pivots = {}
    for item in vectors:
        value = item
        while value:
            bit = value.bit_length() - 1
            if bit in pivots:
                value ^= pivots[bit]
            else:
                pivots[bit] = value
                break
    return len(pivots)


def independent_basis(f, m, d):
    conjugates, value = [], 3
    for _ in range(131):
        conjugates.append(value)
        value = f.square(value)
    assert value == 3 and rank(conjugates) == 131
    assert f.trace(3) == 1
    bases = [[conjugates[m * j + i] for j in range(d)] for i in range(m)]
    assert rank([x for basis in bases for x in basis]) == m * d
    return bases


def masked_x(basis, mask):
    return __import__('functools').reduce(int.__xor__,
                                        (basis[j] for j in range(len(basis)) if mask >> j & 1), 0)


def check_inputs(mod, inputs: dict, curve):
    assert inputs['domain'] == DOMAIN and inputs['kind'] == 'public_synthetic_planted_row_inputs'
    assert inputs['n'] == 131 and inputs['polynomial'] == mod.P131
    assert inputs['q'] == mod.Q131 and inputs['lambda'] == mod.LAMBDA131 and inputs['beta'] == 3
    assert mod.source_order_by_recurrence(131) == 4 * mod.Q131
    torsion = torsion_check(curve)
    assert inputs['torsion'] == [None if t is None else list(t) for t in torsion]
    assert [(c['m'], c['d']) for c in inputs['cells']] == [(5, 25), (6, 21)]
    for cell in inputs['cells']:
        m, d = cell['m'], cell['d']
        assert [r['tuple'] for r in cell['tuples']] == list(range(4))
        bases = independent_basis(curve.f, m, d)
        for witness in cell['tuples']:
            t = witness['tuple']
            slots = witness['slots']
            assert len(slots) == m and [s['slot'] for s in slots] == list(range(m))
            for i, slot in enumerate(slots):
                if t == 1 and i == 0:
                    assert slot == {'slot': 0, 'mask': 0, 'counter': None, 'rejected': [],
                                    'sign_bit': 0, 'base_point': [0, 1], 'source_point': [0, 1]}
                    continue
                if t == 0 and i > 0:
                    assert slot['mask'] == slots[0]['mask']
                    assert slot['counter'] == slots[0]['counter']
                    assert slot['rejected'] == []
                else:
                    assert slot['counter'] is not None and 0 <= slot['counter'] < 128
                    rejected = []
                    for counter in range(slot['counter']):
                        label = f'{DOMAIN}/{m}/{d}/{t}/{i}/{counter}'
                        mask = int.from_bytes(hashlib.sha256(label.encode()).digest(), 'big') % (1 << d)
                        if mask == 0:
                            rejected.append({'counter': counter, 'reason': 'zero'})
                        else:
                            x = masked_x(bases[0], mask)
                            rhs = x ^ curve.f.square(curve.f.inv(x))
                            assert curve.f.trace(rhs) == 1
                            rejected.append({'counter': counter, 'reason': 'no_rational_lift'})
                    assert rejected == slot['rejected']
                    label = f"{DOMAIN}/{m}/{d}/{t}/{i}/{slot['counter']}"
                    expected_mask = int.from_bytes(hashlib.sha256(label.encode()).digest(), 'big') % (1 << d)
                    assert slot['mask'] == expected_mask
                x = masked_x(bases[0], slot['mask'])
                assert x != 0 and masked_x(bases[i], slot['mask']) == point(slot['source_point'])[0]
                rhs = x ^ curve.f.square(curve.f.inv(x))
                assert curve.f.trace(rhs) == 0
                z = mod.independent_half_trace(curve.f, rhs)
                base = (x, curve.f.mul(x, z))
                sign_label = f"{DOMAIN}/{m}/{d}/{t}/{i}/{slot['counter']}/sign"
                bit = hashlib.sha256(sign_label.encode()).digest()[-1] & 1
                assert slot['sign_bit'] == bit
                signed = neg(base) if bit else base
                assert point(slot['base_point']) == signed
                source = signed
                for _ in range(i):
                    source = curve.tau(source)
                assert point(slot['source_point']) == source
            source = [point(slot['source_point']) for slot in slots]
            raw = many(curve, source)
            projected = curve.scalar(raw, 4)
            target = curve.scalar(projected, pow(4, -1, mod.Q131))
            tor = curve.add(raw, neg(target))
            assert tor in torsion and curve.scalar(target, mod.Q131) is None
            assert witness['raw_sum'] == (None if raw is None else list(raw))
            assert witness['synthetic_Q'] == (None if target is None else list(target))
            assert witness['torsion_index'] == torsion.index(tor)
            assert witness['projected_Q'] == (None if projected is None else list(projected))
    return torsion


def check_planted(mod, inputs: dict, raw: Path, curve, torsion):
    reports = []
    for cell in inputs['cells']:
        m, d = cell['m'], cell['d']
        start = time.monotonic()
        stem = f'n131-m{m}-d{d}'
        rows = json_lines(raw / f'{stem}-rows.jsonl')
        summary = json.loads((raw / f'{stem}-summary.json').read_text())
        assert len(rows) == len(cell['tuples']) == summary['rows'] == 4
        columns = set()
        term_columns = set()
        weights = collections.Counter()
        for got, witness in zip(rows, cell['tuples']):
            source = [point(slot['source_point']) for slot in witness['slots']]
            expected = independent_row(curve, 131, mod.Q131, mod.LAMBDA131,
                                       source, point(witness['synthetic_Q']),
                                       witness['torsion_index'], torsion, None)
            expected.update({'n': 131, 'm': m, 'd': d, 'tuple': witness['tuple'],
                             'kind': 'public_synthetic_planted'})
            assert got == expected
            columns.update(tuple(item['column']) for item in got['row'])
            term_columns.update(tuple(item['canonical_column']) for item in got['terms']
                                if item['canonical_column'] is not None)
            weights[got['weight']] += 1
        assert (summary['n'], summary['m'], summary['d']) == (131, m, d)
        assert summary['distinct_columns'] == len(columns)
        assert summary['column_coordinates'] == [list(c) for c in sorted(columns)]
        assert summary['distinct_term_columns'] == len(term_columns)
        assert summary['term_column_coordinates'] == [list(c) for c in sorted(term_columns)]
        assert {int(k): v for k, v in summary['row_weights'].items()} == dict(weights)
        assert summary['total']['wall_seconds'] <= 30 and summary['total']['peak_rss_bytes'] <= 512 * 1024 * 1024
        reports.append({'m': m, 'd': d, 'rows': 4, 'distinct_columns': len(columns),
                        'distinct_term_columns': len(term_columns),
                        'row_weights': dict(sorted(weights.items())),
                        'wall_seconds': time.monotonic() - start})
    return reports


def replay(raw: Path, input_path: Path):
    begin, cpu = time.monotonic(), time.process_time()
    assert digest(OLD_ARCHIVE) == OLD_SHA
    mod = prior()
    field = mod.GF(13, mod.P13)
    curve = mod.E(field)
    by_x, points = mod.roots_and_points(curve)
    assert len(points) + 1 == mod.source_order_by_recurrence(13) == 8012
    torsion = torsion_check(curve)
    H = next(curve.scalar(p, 4) for p in points if curve.scalar(p, 4) is not None)
    assert H is not None and curve.tau(H) == curve.scalar(H, 89)
    with tarfile.open(OLD_ARCHIVE, 'r:gz') as tar:
        toy = [check_arm(mod, tar, raw, curve, by_x, b, m, count, H, torsion)
               for b, m, count in ARMS]
    inputs = json.loads(input_path.read_text())
    field131 = mod.GF(131, mod.P131)
    curve131 = mod.E(field131)
    torsion131 = check_inputs(mod, inputs, curve131)
    planted = check_planted(mod, inputs, raw, curve131, torsion131)
    saved = json.loads((raw / 'summary.json').read_text())
    assert saved['old_archive_sha256'] == OLD_SHA
    assert saved['planted_input_sha256'] == digest(input_path)
    assert [(r['beta'], r['m']) for r in saved['toy']] == [(r['beta'], r['m']) for r in toy]
    assert [(r['m'], r['d']) for r in saved['planted']] == [(r['m'], r['d']) for r in planted]
    assert saved['peak_rss_bytes'] <= 512 * 1024 * 1024
    report = {'status': 'success', 'old_archive_sha256': OLD_SHA,
              'planted_input_sha256': digest(input_path), 'toy': toy, 'planted': planted,
              'total_wall_seconds': time.monotonic() - begin,
              'total_cpu_seconds': time.process_time() - cpu,
              'peak_rss_bytes': mod.peak_rss_bytes(),
              'n13_operations': {'field': dict(field.ops), 'group': dict(curve.ops)},
              'n131_operations': {'field': dict(field131.ops), 'group': dict(curve131.ops)}}
    assert report['total_wall_seconds'] <= 600
    return report


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--raw', type=Path, required=True)
    parser.add_argument('--inputs', type=Path, required=True)
    parser.add_argument('--report', type=Path, required=True)
    args = parser.parse_args()
    args.report.write_text(json.dumps(replay(args.raw, args.inputs),
                                      sort_keys=True, separators=(',', ':')) + '\n')


if __name__ == '__main__':
    main()
