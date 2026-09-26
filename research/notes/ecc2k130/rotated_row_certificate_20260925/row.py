#!/usr/bin/env python3
"""Frozen producer for exact Frobenius-compressed relation-row certificates."""
from __future__ import annotations

import argparse
import collections
import hashlib
import importlib.util
import json
import resource
import signal
import sys
import tarfile
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
PRIOR = HERE.parent / 'rotated_subspace_support_20260925'
OLD_ARCHIVE = PRIOR / 'evidence/raw.tar.gz'
OLD_SHA = 'fad4b4cc48f980349d1ef472d18327368dd0f93f3d1284fc8bbdc20f0d3086d7'
ARMS = ((3, 5), (3, 6), (7, 5), (7, 6))
EXPECTED = (1799, 6307, 1799, 8012)
N131_CELLS = ((5, 25), (6, 21))


def old_module():
    spec = importlib.util.spec_from_file_location('prior_gate', PRIOR / 'gate.py')
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_member(archive: tarfile.TarFile, name: str):
    member = archive.getmember('raw/toy/' + name)
    assert member.isfile() and member.size < 20_000_000
    data = archive.extractfile(member)
    assert data is not None
    return data.read().decode()


def point(value):
    return None if value is None else tuple(value)


def save(path: Path, obj):
    path.write_text(json.dumps(obj, sort_keys=True, separators=(',', ':')) + '\n')


def rss_bytes() -> int:
    raw = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return raw if sys.platform == 'darwin' else raw * 1024


def snapshot(field, curve):
    return (time.monotonic(), time.process_time(), dict(field.operations),
            dict(curve.operations), rss_bytes())


def delta(before, field, curve, name: str):
    wall, cpu, f_ops, c_ops, rss = before
    keys = set(f_ops) | set(c_ops) | set(field.operations) | set(curve.operations)
    ops = {key: field.operations.get(key, 0) - f_ops.get(key, 0) +
           curve.operations.get(key, 0) - c_ops.get(key, 0) for key in sorted(keys)}
    return {'name': name, 'wall_seconds': time.monotonic() - wall,
            'cpu_seconds': time.process_time() - cpu,
            'peak_rss_bytes': rss_bytes(), 'peak_rss_start_bytes': rss,
            'operations': ops}


def deadline(_signum, _frame):
    raise TimeoutError('hard per-cell deadline exceeded')


def start_deadline(seconds: int):
    signal.signal(signal.SIGALRM, deadline)
    signal.setitimer(signal.ITIMER_REAL, seconds)


def stop_deadline():
    signal.setitimer(signal.ITIMER_REAL, 0)


def add_many(curve, values):
    total = None
    for value in values:
        total = curve.add(total, value)
    return total


def split_cofactor(curve, raw, q: int, torsion):
    four = curve.scalar(raw, 4)
    target = curve.scalar(four, pow(4, -1, q))
    tor = curve.add(raw, curve.neg(target))
    assert tor in torsion and curve.scalar(target, q) is None
    assert curve.add(target, tor) == raw
    return target, torsion.index(tor)


def relation(curve, n: int, q: int, lam: int, source, target, ti: int,
             torsion, indices=None):
    m = len(source)
    assert all(curve.on_curve(p) for p in source)
    assert curve.on_curve(target)
    raw = add_many(curve, source)
    assert raw == curve.add(target, torsion[ti])
    terms = []
    coefficients = collections.defaultdict(int)
    for i, p in enumerate(source):
        back = p
        for _ in range((n - i) % n):
            back = curve.tau(back)
        assert curve.on_curve(back)
        transported = curve.scalar(back, 4)
        multiplier = pow(lam, i, q)
        assert curve.scalar(p, 4) == curve.scalar(transported, multiplier)
        if transported is None:
            assert p in torsion and back in torsion
            term = {'slot': i, 'source_index': None if indices is None else indices[i],
                    'source_point': None if p is None else list(p),
                    'back_point': None if back is None else list(back),
                    'projected_point': None, 'canonical_column': None,
                    'sign': 0, 'coefficient': 0,
                    'torsion_index': torsion.index(p)}
        else:
            other = curve.neg(transported)
            column = min(transported, other)
            sign = 1 if column == transported else -1
            coeff = (sign * multiplier) % q
            coefficients[column] = (coefficients[column] + coeff) % q
            term = {'slot': i, 'source_index': None if indices is None else indices[i],
                    'source_point': list(p), 'back_point': list(back),
                    'projected_point': list(transported),
                    'canonical_column': list(column), 'sign': sign,
                    'coefficient': coeff, 'torsion_index': None}
        terms.append(term)
    counts = collections.Counter(tuple(t['canonical_column']) for t in terms
                                 if t['canonical_column'] is not None)
    cancelled = sum(1 for column, count in counts.items()
                    if count > 1 and coefficients[column] == 0)
    row = [{'column': list(column), 'coefficient': coeff}
           for column, coeff in sorted(coefficients.items()) if coeff]
    evaluated = add_many(curve, (curve.scalar(tuple(item['column']), item['coefficient'])
                                 for item in row))
    assert evaluated == curve.scalar(target, 4)
    assert evaluated == add_many(curve, (curve.scalar(p, 4) for p in source))
    return {'source_indices': indices, 'source_points': [None if p is None else list(p) for p in source],
            'target_Q': None if target is None else list(target), 'torsion_index': ti,
            'raw_sum': None if raw is None else list(raw), 'terms': terms,
            'row': row, 'row_point': None if evaluated is None else list(evaluated),
            'weight': len(row), 'zero_terms': sum(t['projected_point'] is None for t in terms),
            'repeated_terms': sum(count - 1 for count in counts.values() if count > 1),
            'cancelled_columns': cancelled}


def arm_data(archive: tarfile.TarFile, beta: int, m: int):
    prefix = f'n13-b{beta}-m{m}-rotated-'
    factors = json.loads(read_member(archive, prefix + 'factors.json'))
    targets = [json.loads(line) for line in read_member(archive, prefix + 'targets.jsonl').splitlines()]
    summary = json.loads(read_member(archive, prefix + 'summary.json'))
    assert len(factors) == m and len(targets) == 2003
    assert (summary['beta'], summary['m'], summary['policy']) == (beta, m, 'rotated')
    return factors, targets, summary


def toy_arm(mod, archive, beta: int, m: int, expected: int, out: Path):
    field = mod.Field(13, mod.MODELS[13]['low'])
    curve = mod.Curve(field)
    torsion = mod.four_torsion(curve)
    start = snapshot(field, curve)
    start_deadline(120)
    try:
        factors, targets, prior = arm_data(archive, beta, m)
        H = point(targets[1]['point'])
        assert mod.source_group_order(13) == 4 * 2003
        assert H is not None and curve.scalar(H, 2003) is None
        assert curve.tau(H) == curve.scalar(H, 89)
        for i, base in enumerate(factors):
            assert base and all(curve.on_curve(point(p)) for p in base)
            for p in base:
                b = point(p)
                for _ in range((13 - i) % 13):
                    b = curve.tau(b)
                assert list(b) in factors[0]
        setup = delta(start, field, curve, 'factor_and_archive_setup')
        run_start = snapshot(field, curve)
        count = 0
        hits = [0] * 4
        columns = set()
        term_columns = set()
        weights = collections.Counter()
        zero_terms = repeated = cancelled = 0
        H_walk = None
        roster_path = out / f'n13-b{beta}-m{m}-roster.jsonl'
        rows_path = out / f'n13-b{beta}-m{m}-rows.jsonl'
        with roster_path.open('w') as roster, rows_path.open('w') as rows:
            for k, record in enumerate(targets):
                assert record['k'] == k and point(record['point']) == H_walk
                project = curve.scalar(H_walk, 4)
                assert point(record['projected_point']) == project
                entries = []
                for ti in range(4):
                    idx = record['coset_witness_indices'][ti]
                    multiplicity = record['coset_multiplicities'][ti]
                    if idx is None:
                        assert multiplicity == 0
                        entries.append(None)
                        continue
                    assert multiplicity > 0 and len(idx) == m
                    source = [point(factors[i][j]) for i, j in enumerate(idx)]
                    row = relation(curve, 13, 2003, 89, source, H_walk, ti, torsion, idx)
                    row.update({'beta': beta, 'm': m, 'k': k, 'status': 'archive_witness'})
                    rows.write(json.dumps(row, sort_keys=True, separators=(',', ':')) + '\n')
                    entries.append(idx)
                    count += 1
                    hits[ti] += 1
                    weights[row['weight']] += 1
                    zero_terms += row['zero_terms']
                    repeated += row['repeated_terms']
                    cancelled += row['cancelled_columns']
                    columns.update(tuple(item['column']) for item in row['row'])
                    term_columns.update(tuple(item['canonical_column']) for item in row['terms']
                                        if item['canonical_column'] is not None)
                roster.write(json.dumps({'k': k, 'point': record['point'],
                                         'coset_multiplicities': record['coset_multiplicities'],
                                         'coset_witness_indices': entries},
                                        sort_keys=True, separators=(',', ':')) + '\n')
                H_walk = curve.add(H_walk, H)
        assert H_walk is None and count == expected == prior['distinct_full_sums']
        assert hits == prior['coset_hits']
        processing = delta(run_start, field, curve, 'all_roster_rows_and_misses')
        control_start = snapshot(field, curve)
        p0 = next(point(p) for p in factors[0] if curve.scalar(point(p), 4) is not None)
        assert p0 is not None
        controls = []
        for name, source in (('all_x0', [(0, 1)] * m),
                             ('repeated_signed', [p0, curve.tau(curve.neg(p0))] + [(0, 1)] * (m - 2))):
            raw = add_many(curve, source)
            qpt, ti = split_cofactor(curve, raw, 2003, torsion)
            row = relation(curve, 13, 2003, 89, source, qpt, ti, torsion)
            row.update({'beta': beta, 'm': m, 'control': name})
            assert (name != 'all_x0' or row['weight'] == 0)
            assert (name != 'repeated_signed' or row['repeated_terms'] >= 1)
            controls.append(row)
        save(out / f'n13-b{beta}-m{m}-controls.json', controls)
        control_cost = delta(control_start, field, curve, 'exception_controls')
        assert rss_bytes() <= 512 * 1024 * 1024
        result = {'beta': beta, 'm': m, 'n': 13, 'archive_positive_rows': count,
                  'coset_hits': hits, 'distinct_columns': len(columns),
                  'column_coordinates': [list(c) for c in sorted(columns)],
                  'distinct_term_columns': len(term_columns),
                  'term_column_coordinates': [list(c) for c in sorted(term_columns)],
                  'row_weights': dict(sorted(weights.items())),
                  'zero_terms': zero_terms, 'repeated_terms': repeated,
                  'cancelled_columns': cancelled,
                  'stages': [setup, processing, control_cost],
                  'total': delta(start, field, curve, 'total')}
        save(out / f'n13-b{beta}-m{m}-summary.json', result)
        return result
    finally:
        stop_deadline()


def n131_cell(mod, cell: dict, out: Path):
    field = mod.Field(131, mod.MODELS[131]['low'])
    curve = mod.Curve(field)
    start = snapshot(field, curve)
    start_deadline(30)
    try:
        torsion = mod.four_torsion(curve)
        assert mod.source_group_order(131) == 4 * mod.Q131
        rows = []
        columns = set()
        term_columns = set()
        weights = collections.Counter()
        for witness in cell['tuples']:
            source = [point(slot['source_point']) for slot in witness['slots']]
            target = point(witness['synthetic_Q'])
            ti = witness['torsion_index']
            row = relation(curve, 131, mod.Q131, mod.LAMBDA131,
                           source, target, ti, torsion)
            assert row['raw_sum'] == witness['raw_sum']
            assert row['row_point'] == witness['projected_Q']
            row.update({'n': 131, 'm': cell['m'], 'd': cell['d'],
                        'tuple': witness['tuple'], 'kind': 'public_synthetic_planted'})
            rows.append(row)
            columns.update(tuple(item['column']) for item in row['row'])
            term_columns.update(tuple(item['canonical_column']) for item in row['terms']
                                if item['canonical_column'] is not None)
            weights[row['weight']] += 1
        assert len(rows) == 4 and rss_bytes() <= 512 * 1024 * 1024
        label = f"n131-m{cell['m']}-d{cell['d']}"
        with (out / f'{label}-rows.jsonl').open('w') as handle:
            for row in rows:
                handle.write(json.dumps(row, sort_keys=True, separators=(',', ':')) + '\n')
        result = {'n': 131, 'm': cell['m'], 'd': cell['d'], 'rows': len(rows),
                  'distinct_columns': len(columns), 'column_coordinates': [list(c) for c in sorted(columns)],
                  'distinct_term_columns': len(term_columns),
                  'term_column_coordinates': [list(c) for c in sorted(term_columns)],
                  'row_weights': dict(sorted(weights.items())),
                  'total': delta(start, field, curve, 'total_cell')}
        save(out / f'{label}-summary.json', result)
        return result
    finally:
        stop_deadline()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--inputs', type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    assert digest(OLD_ARCHIVE) == OLD_SHA
    inputs = json.loads(args.inputs.read_text())
    assert inputs['kind'] == 'public_synthetic_planted_row_inputs'
    assert [(c['m'], c['d']) for c in inputs['cells']] == list(N131_CELLS)
    mod = old_module()
    overall_wall = time.monotonic()
    overall_cpu = time.process_time()
    with tarfile.open(OLD_ARCHIVE, 'r:gz') as archive:
        toy = [toy_arm(mod, archive, b, m, expected, args.out)
               for (b, m), expected in zip(ARMS, EXPECTED)]
    planted = [n131_cell(mod, cell, args.out) for cell in inputs['cells']]
    save(args.out / 'summary.json', {'old_archive_sha256': OLD_SHA,
         'planted_input_sha256': digest(args.inputs), 'toy': toy, 'planted': planted,
         'total_wall_seconds': time.monotonic() - overall_wall,
         'total_cpu_seconds': time.process_time() - overall_cpu,
         'peak_rss_bytes': rss_bytes()})


if __name__ == '__main__':
    main()
