#!/usr/bin/env python3
"""Gray-code exact low/high F0 census for full-dimensional rotated slots."""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import resource
import signal
import sys
import time
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
PARENT = HERE.parent / 'rotated_subspace_support_20260925/gate.py'
spec = importlib.util.spec_from_file_location('unequal_parent_gate', PARENT)
assert spec is not None and spec.loader is not None
gate = importlib.util.module_from_spec(spec)
spec.loader.exec_module(gate)


def save(path: Path, value) -> None:
    path.write_text(json.dumps(value, sort_keys=True, separators=(',', ':')) + '\n')


def peak_rss() -> int:
    value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return value if sys.platform == 'darwin' else value * 1024


def counter_delta(after: dict, before: dict) -> dict:
    return {key: after.get(key, 0) - before.get(key, 0) for key in sorted(set(after) | set(before))}


def scan_space(f, basis: list[int], label: str, arm: dict, data: dict, out: Path):
    d = len(basis)
    assert d == arm['d_low'] + (label == 'high')
    out.mkdir()
    started, cpu_started = time.perf_counter(), time.process_time()
    before_ops = dict(f.operations)
    trace_mask = gate.trace_mask(f)
    assert all(gate.trace_fast(trace_mask, b) == 1 for b in basis)
    assert gate.rank(basis) == d and gate.rank(basis + [1]) == d + 1
    columns: dict[int, int] = {}
    row_hash, chunk_hash = hashlib.sha256(), hashlib.sha256()
    x = 0
    zero_x = one_x = liftable = chunk_liftable = chunk_index = 0
    limit, chunk = 1 << d, data['chunk_rows']
    with (out / 'chunks.jsonl').open('w') as stream:
        for ordinal in range(limit):
            if ordinal:
                changed = ordinal & -ordinal
                x ^= basis[changed.bit_length() - 1]
            mask = ordinal ^ (ordinal >> 1)
            if x == 0:
                assert ordinal == 0 and mask == 0
                zero_x += 1
                lift, projected = 1, -1
            else:
                assert x != 1
                inverse_x = f.inverse(x)
                trace = (mask.bit_count() ^ (trace_mask & inverse_x).bit_count()) & 1
                if trace:
                    lift, projected = 0, -2
                else:
                    liftable += 1
                    chunk_liftable += 1
                    u = f.square(x ^ inverse_x)
                    assert u != 0
                    projected = f.square(u ^ f.inverse(u))
                    assert projected != 0
                    columns[projected] = columns.get(projected, 0) + 1
                    assert columns[projected] <= 4
                    lift = 2
            row = f'{mask},{x},{lift},{projected}\n'.encode('ascii')
            row_hash.update(row)
            chunk_hash.update(row)
            if (ordinal + 1) % chunk == 0 or ordinal + 1 == limit:
                assert peak_rss() <= data['caps']['rss_bytes']
                stream.write(json.dumps({'start_ordinal': chunk_index * chunk,
                                         'stop_ordinal': ordinal + 1,
                                         'row_sha256': chunk_hash.hexdigest(),
                                         'liftable_nonzero_x': chunk_liftable,
                                         'column_count_so_far': len(columns)},
                                        sort_keys=True, separators=(',', ':')) + '\n')
                stream.flush()
                chunk_index += 1
                chunk_hash = hashlib.sha256()
                chunk_liftable = 0
    assert zero_x == 1 and one_x == 0 and sum(columns.values()) == liftable
    multiplicities = Counter(columns.values())
    assert not set(multiplicities) - {1, 2, 3, 4}
    column_hash = hashlib.sha256()
    for v in sorted(columns):
        column_hash.update(f'{v}\n'.encode('ascii'))
    result = {'label': label, 'dimension': d, 'total_masks': limit,
              'zero_x_count': zero_x, 'one_x_count': one_x,
              'liftable_nonzero_x': liftable, 'physical_f0_points': 1 + 2 * liftable,
              'nonzero_signed_columns': len(columns),
              'projected_signed_columns_including_O': 1 + len(columns),
              'column_x_multiplicity_histogram': {str(i): multiplicities[i] for i in range(1, 5)},
              'row_sha256': row_hash.hexdigest(),
              'column_set_sha256': column_hash.hexdigest(), 'chunk_count': chunk_index,
              'field_operations': counter_delta(dict(f.operations), before_ops),
              'wall_seconds': time.perf_counter() - started,
              'cpu_seconds': time.process_time() - cpu_started,
              'peak_rss_bytes': peak_rss()}
    assert result['wall_seconds'] <= data['caps']['producer_wall_seconds']
    assert result['peak_rss_bytes'] <= data['caps']['rss_bytes']
    save(out / 'result.json', result)
    return columns, result


def run(arm: dict, data: dict, out: Path):
    started, cpu_started = time.perf_counter(), time.process_time()
    m, low_d, high_d, r = arm['m'], arm['d_low'], arm['d_high'], arm['r']
    assert high_d == low_d + 1 and r == 131 - m * low_d
    f = gate.Field(131, [0, 1, 2, 13])
    assert f.poly == data['field_poly']
    f.rabin_prime_degree()
    assert gate.source_group_order(131) == data['group_order'] == 4 * data['q']
    conjugates = gate.normal_conjugates(f, data['beta'])
    assert gate.rank(conjugates) == 131
    slot_dims = [high_d if i < r else low_d for i in range(m)]
    slot_indices = [[m * j + i for j in range(d)] for i, d in enumerate(slot_dims)]
    flat_indices = [index for slot in slot_indices for index in slot]
    assert len(flat_indices) == 131 and set(flat_indices) == set(range(131))
    slots = [[conjugates[index] for index in indices] for indices in slot_indices]
    assert all(gate.rank(slot) == d for slot, d in zip(slots, slot_dims))
    assert gate.rank([value for slot in slots for value in slot]) == 131
    low_basis = [conjugates[m * j] for j in range(low_d)]
    high_basis = [conjugates[m * j] for j in range(high_d)]
    assert low_basis == high_basis[:-1] and all(gate.rank(b + [1]) == len(b) + 1
                                               for b in (low_basis, high_basis))
    for i, slot in enumerate(slots):
        basis = high_basis if i < r else low_basis
        assert slot == [frobenius(f, value, i) for value in basis]
    setup_wall, setup_cpu = time.perf_counter() - started, time.process_time() - cpu_started
    low_columns, low = scan_space(f, low_basis, 'low', arm, data, out / 'low')
    high_columns, high = scan_space(f, high_basis, 'high', arm, data, out / 'high')
    assert set(low_columns) <= set(high_columns)
    physical_tuple_count = high['physical_f0_points'] ** r * low['physical_f0_points'] ** (m - r)
    threshold = data['threshold']
    followup = (physical_tuple_count * threshold['minimum_ordered_tuple_count_over_q_den'] >=
                data['q'] * threshold['minimum_ordered_tuple_count_over_q_num'] and
                high['nonzero_signed_columns'] <= threshold['max_nonzero_signed_columns'])
    result = {'arm': arm, 'domain': data['domain'], 'beta': data['beta'],
              'field_poly': f.poly, 'q': data['q'],
              'slot_dimensions': slot_dims, 'combined_slot_rank': 131,
              'normal_basis_rank': 131,
              'normalized_low_columns_subset_high': True,
              'normalized_global_nonzero_signed_columns': high['nonzero_signed_columns'],
              'low_nonzero_signed_columns': low['nonzero_signed_columns'],
              'ordered_physical_tuples': physical_tuple_count,
              'necessary_support_upper_num': min(physical_tuple_count, data['q']),
              'necessary_support_upper_den': data['q'],
              'raw_affine_chain_variables': 131 + 131 * (m - 2),
              'followup_screen_pass': followup,
              'space_result_sha256': {label: hashlib.sha256((out / label / 'result.json').read_bytes()).hexdigest()
                                      for label in ('low', 'high')},
              'setup_wall_seconds': setup_wall, 'setup_cpu_seconds': setup_cpu,
              'total_field_operations': dict(f.operations),
              'total_wall_seconds': time.perf_counter() - started,
              'total_cpu_seconds': time.process_time() - cpu_started,
              'peak_rss_bytes': peak_rss()}
    assert result['total_wall_seconds'] <= data['caps']['producer_wall_seconds']
    assert result['peak_rss_bytes'] <= data['caps']['rss_bytes']
    save(out / 'result.json', result)


def frobenius(f, value, i):
    for _ in range(i):
        value = f.square(value)
    return value


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--m', type=int, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    data = json.loads((HERE / 'INPUT.json').read_text())
    arm = next((row for row in data['arms'] if row['m'] == args.m), None)
    assert arm is not None
    args.out.mkdir(parents=True, exist_ok=False)
    started = time.perf_counter()
    signal.signal(signal.SIGALRM, lambda *_: (_ for _ in ()).throw(TimeoutError('producer wall cap')))
    signal.alarm(data['caps']['producer_wall_seconds'])
    try:
        run(arm, data, args.out)
    except BaseException as error:
        save(args.out / 'failure.json', {'arm': arm, 'error': repr(error),
                                        'wall_seconds': time.perf_counter() - started,
                                        'peak_rss_bytes': peak_rss()})
        raise
    finally:
        signal.alarm(0)


if __name__ == '__main__':
    main()
