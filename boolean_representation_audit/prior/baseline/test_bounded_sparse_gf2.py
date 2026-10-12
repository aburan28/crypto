"""Independent correctness checks and bounded memory measurements."""

from array import array
from hashlib import sha256
from itertools import combinations
from pathlib import Path
from random import Random
from time import perf_counter
import json
import platform
import tracemalloc

from bounded_sparse_gf2 import (
    BoundedMonomials, ResourceLimitError, SparseEchelon, xor_sorted,
)


def correctness():
    counts = {}
    for n, d in [(6, 3), (8, 4), (30, 4), (31, 4)]:
        space = BoundedMonomials(n, d)
        if n <= 8:
            expected = [sum(1 << i for i in positions)
                        for k in range(d + 1) for positions in combinations(range(n), k)]
            ranks = [space.rank(mask) for mask in expected]
            assert set(ranks) == set(range(space.count))
            for mask in expected:
                assert space.unrank(space.rank(mask)) == mask
        else:
            rng = Random(131 * n + d)
            indices = {0, 1, space.count - 1} | {
                rng.randrange(space.count) for _ in range(10_000)}
            for index in indices:
                assert space.rank(space.unrank(index)) == index
        counts[f'n{n}_d{d}'] = space.count

    space = BoundedMonomials(8, 3)
    rng = Random(5103)
    matrix = SparseEchelon(space, max_row_terms=100, max_total_terms=5000)
    dense_pivots = {}
    independent = dependent = 0
    for _ in range(180):
        raw = [rng.randrange(space.count) for _ in range(rng.randrange(2, 30))]
        row = space.row_from_indices(raw)
        assert set(row) == {v for v in raw if raw.count(v) % 2}
        dense = 0
        for index in raw: dense ^= 1 << index
        while dense and (dense.bit_length() - 1) in dense_pivots:
            dense ^= dense_pivots[dense.bit_length() - 1]
        expected_independent = bool(dense)
        if dense:
            dense_pivots[dense.bit_length() - 1] = dense
        assert matrix.add(row) == expected_independent
        assert matrix.rank == len(dense_pivots)
        independent += expected_independent
        dependent += not expected_independent
    assert dependent > 0 and independent > 0

    a = space.row_from_masks([1, 3, 3, 2])
    b = space.row_from_masks([1, 5])
    assert [space.unrank(i) for i in xor_sorted(a, b, space.typecode)] == [2, 5]
    quadratic = BoundedMonomials(8, 2)
    low = quadratic.row_from_masks([1, 2, 3])
    product, dropped = quadratic.multiply_row(low, 1)
    assert {quadratic.unrank(i) for i in product} == {1}
    assert dropped == 0
    high = quadratic.row_from_masks([1 << 1 | 1 << 2])
    product, dropped = quadratic.multiply_row(high, 1)
    assert len(product) == 0 and dropped == 1

    limited = SparseEchelon(space, max_row_terms=3, max_total_terms=3)
    try:
        limited.add(space.row_from_indices([0, 1, 2, 3]))
        raise AssertionError('row limit did not reject')
    except ResourceLimitError:
        assert limited.rank == 0
    limited.add(space.row_from_indices([0, 1, 2]))
    try:
        limited.add(space.row_from_indices([3]))
        raise AssertionError('matrix limit did not reject')
    except ResourceLimitError:
        assert limited.rank == 1
    fill_limited = SparseEchelon(space, max_row_terms=3, max_total_terms=10)
    assert fill_limited.add(space.row_from_indices([0, 1, 5]))
    try:
        fill_limited.add(space.row_from_indices([2, 3, 5]))
        raise AssertionError('fill-in limit did not reject')
    except ResourceLimitError:
        assert fill_limited.rank == 1

    # The packed index width upgrades without allocating C(n,d) entries.
    wide = BoundedMonomials(90, 7)
    assert wide.typecode == 'Q'
    boundary = (1 << 89) | (1 << 42) | (1 << 2)
    assert wide.unrank(wide.rank(boundary)) == boundary
    assert wide.row_from_masks([boundary]).itemsize == 8
    return dict(counts=counts, compared_rows=180, independent_rows=independent,
                dependent_rows=dependent, budget_rejections=3,
                wide_index_columns=wide.count)


def measure(n):
    # This measurement includes index construction, deterministic row
    # generation and elimination. A dense 2^n-column row is never allocated.
    tracemalloc.start()
    started = perf_counter()
    space = BoundedMonomials(n, 4)
    rng = Random(20260928 + n)
    matrix = SparseEchelon(space, max_row_terms=512,
                            max_total_terms=256 * 512)
    for _ in range(256):
        row = space.row_from_indices(rng.sample(range(space.count), 64))
        matrix.add(row)
    seconds = perf_counter() - started
    _, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    assert matrix.rank == 256
    return dict(variables=n, max_degree=4, indexable_columns=space.count,
                rows=256, input_terms_per_row=64, rank=matrix.rank,
                stored_pivot_terms=matrix.stored_terms,
                peak_row_terms=matrix.peak_row_terms, xor_steps=matrix.xor_steps,
                index_width_bytes=array(space.typecode).itemsize,
                packed_pivot_payload_bytes=matrix.packed_payload_bytes,
                measured_python_heap_peak_bytes=peak,
                unbounded_dense_row_bytes=(1 << n) // 8,
                bounded_dense_row_bytes=(space.count + 7) // 8,
                elapsed_seconds=round(seconds, 6))


def measure_storage_only():
    """Keep 39,001 packed rows alive without running elimination."""
    tracemalloc.start()
    started = perf_counter()
    space = BoundedMonomials(30, 4)
    rng = Random(39001)
    rows = [space.row_from_indices(rng.sample(range(space.count), 64))
            for _ in range(39_001)]
    seconds = perf_counter() - started
    _, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    assert len(rows) == 39_001 and all(len(row) == 64 for row in rows)
    return dict(variables=30, max_degree=4, rows=len(rows),
                terms_per_row=64, packed_payload_bytes=sum(len(row) * row.itemsize for row in rows),
                measured_python_heap_peak_bytes=peak,
                unbounded_dense_row_bytes=(1 << 30) // 8,
                elapsed_seconds=round(seconds, 6), elimination='not run')


def main():
    checks = correctness()
    measurements = [measure(n) for n in [30, 31]]
    storage = measure_storage_only()
    result = dict(status='PASS', python=platform.python_version(),
                  source_sha256=sha256(Path('bounded_sparse_gf2.py').read_bytes()).hexdigest(),
                  checks=checks, measurements=measurements, storage_only=storage,
                  interpretation='Synthetic matrices; 39,001 rows held without elimination; no full solver or cryptographic workload.')
    Path('bounded_sparse_results.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__': main()
