#!/usr/bin/env python3
"""Independent natural-mask and group-law replay of unequal slot censuses."""
from __future__ import annotations

import argparse
import hashlib
import json
import resource
import signal
import sys
import time
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent


def save(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n")


def peak_rss() -> int:
    value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return value if sys.platform == "darwin" else value * 1024


def poly_divmod(a: int, b: int) -> tuple[int, int]:
    if b == 0:
        raise ZeroDivisionError
    quotient = 0
    while a and a.bit_length() >= b.bit_length():
        shift = a.bit_length() - b.bit_length()
        quotient ^= 1 << shift
        a ^= b << shift
    return quotient, a


def poly_product(a: int, b: int) -> int:
    result = 0
    while b:
        if b & 1:
            result ^= a
        a <<= 1
        b >>= 1
    return result


class Field:
    def __init__(self, n: int, poly: int):
        assert poly.bit_length() == n + 1 and poly & 1
        self.n, self.poly = n, poly
        self.ops = Counter()
        self.square_byte = tuple(sum(((byte >> bit) & 1) << (2 * bit)
                                     for bit in range(8)) for byte in range(256))

    def reduce(self, a: int) -> int:
        self.ops["reduce"] += 1
        return poly_divmod(a, self.poly)[1]

    def mul(self, a: int, b: int) -> int:
        self.ops["mul"] += 1
        return self.reduce(poly_product(a, b))

    def square(self, a: int) -> int:
        self.ops["square"] += 1
        expansion, index = 0, 0
        while a:
            expansion ^= self.square_byte[a & 255] << (16 * index)
            a >>= 8
            index += 1
        return self.reduce(expansion)

    def inverse(self, a: int) -> int:
        assert a != 0
        self.ops["inverse"] += 1
        r0, r1 = self.poly, a
        t0, t1 = 0, 1
        while r1:
            quotient, rem = poly_divmod(r0, r1)
            r0, r1 = r1, rem
            t0, t1 = t1, t0 ^ poly_product(quotient, t1)
        assert r0 == 1
        return self.reduce(t0)

    def trace(self, a: int) -> int:
        total, term = 0, a
        for _ in range(self.n):
            total ^= term
            term = self.square(term)
        assert term == a and total in (0, 1)
        return total

    def trace_mask(self) -> int:
        return sum(self.trace(1 << bit) << bit for bit in range(self.n))


def rank(values: list[int]) -> int:
    rows = {}
    for original in values:
        value = original
        while value:
            pivot = value.bit_length() - 1
            if pivot in rows:
                value ^= rows[pivot]
            else:
                rows[pivot] = value
                break
    return len(rows)


def natural_x(row: list[int], mask: int) -> int:
    x, bit = 0, 0
    while mask:
        if mask & 1:
            x ^= row[bit]
        bit += 1
        mask >>= 1
    return x


def half_trace(f: Field, rhs: int) -> int:
    assert f.trace(rhs) == 0
    result, term = 0, rhs
    for _ in range((f.n + 1) // 2):
        result ^= term
        term = f.square(f.square(term))
    assert f.square(result) ^ result == rhs
    return result


class Curve:
    def __init__(self, f: Field):
        self.f = f
        self.ops = Counter()

    def on(self, point) -> bool:
        if point is None:
            return True
        x, y = point
        return self.f.square(y) ^ self.f.mul(x, y) == self.f.mul(self.f.square(x), x) ^ 1

    @staticmethod
    def neg(point):
        return None if point is None else (point[0], point[0] ^ point[1])

    def add(self, a, b):
        self.ops["add"] += 1
        if a is None:
            return b
        if b is None:
            return a
        x, y = a
        u, v = b
        f = self.f
        if x == u:
            if y ^ v == x:
                return None
            assert y == v and x != 0
            slope = x ^ f.mul(y, f.inverse(x))
            new_x = f.square(slope) ^ slope
            new_y = f.square(x) ^ f.mul(slope ^ 1, new_x)
        else:
            slope = f.mul(y ^ v, f.inverse(x ^ u))
            new_x = f.square(slope) ^ slope ^ x ^ u
            new_y = f.mul(slope, x ^ new_x) ^ new_x ^ y
        return new_x, new_y

    def four(self, point):
        return self.add(self.add(point, point), self.add(point, point))


def sample_ordinals(domain: str, m: int, d: int, count: int) -> set[int]:
    result = set()
    counter = 0
    while len(result) < count:
        digest = hashlib.sha256(f"{domain}/{m}/{d}/sample/{counter}".encode()).digest()
        result.add(int.from_bytes(digest, "big") % (1 << d))
        counter += 1
    return result


def sample_point_check(f: Field, curve: Curve, x: int, lift: int, projected: int) -> bool:
    if x == 0:
        assert lift == 1 and projected == -1 and curve.on((0, 1))
        assert curve.four((0, 1)) is None
        return False
    inv = f.inverse(x)
    assert f.mul(x, inv) == 1
    rhs = x ^ f.square(inv)
    assert f.trace(rhs) == (lift == 0)
    if lift == 0:
        assert projected == -2
        return False
    z = half_trace(f, rhs)
    point = (x, f.mul(x, z))
    assert curve.on(point)
    negative = curve.neg(point)
    assert curve.on(negative) and negative != point
    four = curve.four(point)
    assert four is not None and four[0] == projected
    assert curve.four(negative) == curve.neg(four)
    return True


def normalized_bases(f: Field, beta: int, arm: dict):
    m, low_d, high_d, r = arm['m'], arm['d_low'], arm['d_high'], arm['r']
    assert high_d == low_d + 1 and r == 131 - m * low_d
    conjugates, value = [], beta
    for _ in range(131):
        conjugates.append(value)
        value = f.square(value)
    assert value == beta and rank(conjugates) == 131 and f.trace(beta) == 1
    dims = [high_d if i < r else low_d for i in range(m)]
    indices = [[m*j+i for j in range(d)] for i,d in enumerate(dims)]
    assert sorted(j for slot in indices for j in slot) == list(range(131))
    slots = [[conjugates[j] for j in slot] for slot in indices]
    assert all(rank(slot) == len(slot) for slot in slots)
    assert rank([v for slot in slots for v in slot]) == 131
    low = [conjugates[m*j] for j in range(low_d)]
    high = [conjugates[m*j] for j in range(high_d)]
    assert high[:-1] == low
    assert rank(low + [1]) == low_d + 1
    assert rank(high + [1]) == high_d + 1
    for i, slot in enumerate(slots):
        base = high if i < r else low
        expected = []
        for v in base:
            for _ in range(i):
                v = f.square(v)
            expected.append(v)
        assert slot == expected
    return {'low': low, 'high': high}, dims


def verify_space(f, curve, basis, label, arm, data, directory):
    saved = json.loads((directory / 'result.json').read_text())
    chunks = [json.loads(line) for line in (directory / 'chunks.jsonl').read_text().splitlines()]
    trace_mask = f.trace_mask()
    assert all((trace_mask & value).bit_count() & 1 == 1 for value in basis)
    samples = sample_ordinals(data['domain'], arm['m'], len(basis), data['sample_ordinals_per_space'])
    columns = {}
    row_hash, chunk_hash = hashlib.sha256(), hashlib.sha256()
    chunk_liftable = sample_checks = sample_lifts = 0
    zero_x = one_x = liftable = 0
    limit, chunk = 1 << len(basis), data['chunk_rows']
    for ordinal in range(limit):
        mask = ordinal ^ (ordinal >> 1)
        x = natural_x(basis, mask)
        if x == 0:
            assert ordinal == 0
            zero_x += 1
            lift, projected = 1, -1
        else:
            assert x != 1
            inverse_x = f.inverse(x)
            assert f.mul(x, inverse_x) == 1
            criterion = ((trace_mask & x).bit_count() ^
                         (trace_mask & inverse_x).bit_count()) & 1
            if criterion:
                lift, projected = 0, -2
            else:
                liftable += 1
                chunk_liftable += 1
                u = f.square(x) ^ f.square(inverse_x)
                assert u != 0
                projected = f.square(u) ^ f.square(f.inverse(u))
                assert projected != 0
                columns[projected] = columns.get(projected, 0) + 1
                assert columns[projected] <= 4
                lift = 2
        row = f'{mask},{x},{lift},{projected}\n'.encode('ascii')
        row_hash.update(row)
        chunk_hash.update(row)
        if ordinal in samples:
            sample_checks += 1
            sample_lifts += sample_point_check(f, curve, x, lift, projected)
        if (ordinal + 1) % chunk == 0 or ordinal + 1 == limit:
            index = ordinal // chunk
            assert chunks[index] == {'start_ordinal': index * chunk,
                                     'stop_ordinal': ordinal + 1,
                                     'row_sha256': chunk_hash.hexdigest(),
                                     'liftable_nonzero_x': chunk_liftable,
                                     'column_count_so_far': len(columns)}
            assert peak_rss() <= data['caps']['rss_bytes']
            chunk_hash = hashlib.sha256()
            chunk_liftable = 0
    assert sample_checks == len(samples) == data['sample_ordinals_per_space']
    assert len(chunks) == (limit + chunk - 1) // chunk
    multiplicities = Counter(columns.values())
    column_hash = hashlib.sha256()
    for v in sorted(columns):
        column_hash.update(f'{v}\n'.encode('ascii'))
    expected = {'label': label, 'dimension': len(basis), 'total_masks': limit,
                'zero_x_count': zero_x, 'one_x_count': one_x,
                'liftable_nonzero_x': liftable, 'physical_f0_points': 1 + 2 * liftable,
                'nonzero_signed_columns': len(columns),
                'projected_signed_columns_including_O': 1 + len(columns),
                'column_x_multiplicity_histogram': {str(i): multiplicities[i] for i in range(1,5)},
                'row_sha256': row_hash.hexdigest(),
                'column_set_sha256': column_hash.hexdigest(), 'chunk_count': len(chunks)}
    assert all(saved[key] == value for key, value in expected.items())
    assert saved['wall_seconds'] <= data['caps']['producer_wall_seconds']
    assert saved['peak_rss_bytes'] <= data['caps']['rss_bytes']
    return columns, {'label': label, 'dimension': len(basis),
                     'rows_checked': limit, 'sample_ordinals_checked': sample_checks,
                     'sample_rational_lifts_checked': sample_lifts,
                     'row_sha256': row_hash.hexdigest(),
                     'column_set_sha256': column_hash.hexdigest(),
                     'source_result_sha256': hashlib.sha256((directory / 'result.json').read_bytes()).hexdigest()}


def verify(arm: dict, data: dict, source: Path, out: Path) -> None:
    started, cpu_started = time.perf_counter(), time.process_time()
    f = Field(131, data['field_poly'])
    bases, dims = normalized_bases(f, data['beta'], arm)
    curve = Curve(f)
    columns, spaces = {}, {}
    for label in ('low', 'high'):
        columns[label], spaces[label] = verify_space(f, curve, bases[label], label, arm,
                                                      data, source / label)
    assert set(columns['low']) <= set(columns['high'])
    high, low = json.loads((source / 'high/result.json').read_text()), json.loads((source / 'low/result.json').read_text())
    ordered = high['physical_f0_points'] ** arm['r'] * low['physical_f0_points'] ** (arm['m'] - arm['r'])
    threshold = data['threshold']
    followup = (ordered * threshold['minimum_ordered_tuple_count_over_q_den'] >=
                data['q'] * threshold['minimum_ordered_tuple_count_over_q_num'] and
                len(columns['high']) <= threshold['max_nonzero_signed_columns'])
    expected = {'arm': arm, 'domain': data['domain'], 'beta': data['beta'],
                'field_poly': f.poly, 'q': data['q'], 'slot_dimensions': dims,
                'combined_slot_rank': 131, 'normal_basis_rank': 131,
                'normalized_low_columns_subset_high': True,
                'normalized_global_nonzero_signed_columns': len(columns['high']),
                'low_nonzero_signed_columns': len(columns['low']),
                'ordered_physical_tuples': ordered,
                'necessary_support_upper_num': min(ordered, data['q']),
                'necessary_support_upper_den': data['q'],
                'raw_affine_chain_variables': 131 + 131 * (arm['m'] - 2),
                'followup_screen_pass': followup,
                'space_result_sha256': {label: hashlib.sha256((source / label / 'result.json').read_bytes()).hexdigest()
                                        for label in ('low', 'high')}}
    saved = json.loads((source / 'result.json').read_text())
    assert all(saved[key] == value for key, value in expected.items())
    assert saved['total_wall_seconds'] <= data['caps']['producer_wall_seconds']
    assert saved['peak_rss_bytes'] <= data['caps']['rss_bytes']
    result = {'status': 'PASS', 'arm': arm, 'spaces': spaces,
              'normal_basis_rank': 131, 'combined_slot_rank': 131,
              'normalized_low_columns_subset_high': True,
              'ordered_physical_tuples': ordered,
              'global_nonzero_signed_columns': len(columns['high']),
              'followup_screen_pass': followup,
              'field_operations': dict(f.ops), 'curve_operations': dict(curve.ops),
              'source_result_sha256': hashlib.sha256((source / 'result.json').read_bytes()).hexdigest(),
              'wall_seconds': time.perf_counter() - started,
              'cpu_seconds': time.process_time() - cpu_started,
              'peak_rss_bytes': peak_rss()}
    assert result['wall_seconds'] <= data['caps']['verifier_wall_seconds']
    assert result['peak_rss_bytes'] <= data['caps']['rss_bytes']
    save(out, result)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument('--m', type=int, required=True)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    data = json.loads((HERE / 'INPUT.json').read_text())
    arm = next((row for row in data['arms'] if row['m'] == args.m), None)
    assert arm is not None
    started = time.perf_counter()
    signal.signal(signal.SIGALRM, lambda *_: (_ for _ in ()).throw(TimeoutError('verifier wall cap')))
    signal.alarm(data['caps']['verifier_wall_seconds'])
    try:
        verify(arm, data, args.source, args.out)
    except BaseException as error:
        save(args.out.with_suffix('.failure.json'),
             {'arm': arm, 'error': repr(error), 'wall_seconds': time.perf_counter() - started,
              'peak_rss_bytes': peak_rss()})
        raise
    finally:
        signal.alarm(0)


if __name__ == '__main__':
    main()
