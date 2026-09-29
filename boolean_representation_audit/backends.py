"""Representation-only comparison for the recovered n<=10 Boolean audit.

The prior engine is loaded twice, so its scheduling code is byte-identical.
The second module receives a packed row/space/reducer backend at import time.
There is no mutable per-call monkeypatch and no application integration.
"""
import importlib.util
from pathlib import Path
import sys

PRIOR = Path(__file__).parent / 'prior'
sys.path.insert(0, str(PRIOR))
from baseline.bounded_sparse_gf2 import BoundedMonomials, ResourceLimitError
from boolean_closure import compute as sparse_compute, normalize, solutions, verify


class PackedRow:
    __slots__ = ('bits',)

    def __init__(self, bits):
        self.bits = bits

    def __len__(self):
        return self.bits.bit_count()

    def __bool__(self):
        return bool(self.bits)

    def __eq__(self, other):
        return isinstance(other, PackedRow) and self.bits == other.bits

    def __iter__(self):
        value = self.bits
        while value:
            low = value & -value
            yield low.bit_length() - 1
            value ^= low

    def __getitem__(self, index):
        if index != -1 or not self.bits:
            raise IndexError('only the leading column is exposed')
        return self.bits.bit_length() - 1


class PackedSpace(BoundedMonomials):
    def row_from_indices(self, indices):
        value = 0
        for index in indices:
            if type(index) is not int or not 0 <= index < self.count:
                raise IndexError('column outside declared space')
            value ^= 1 << index
        return PackedRow(value)

    def multiply_row(self, row, monomial_mask):
        if type(monomial_mask) is not int or not 0 <= monomial_mask < (1 << self.n):
            raise ValueError('multiplier outside variable universe')
        if not isinstance(row, PackedRow):
            raise TypeError('expected packed row')
        value, dropped = 0, 0
        for column in row:
            mask = self.unrank(column) | monomial_mask
            if mask.bit_count() > self.max_degree:
                dropped += 1
            else:
                value ^= 1 << self.rank(mask)
        return PackedRow(value), dropped


class PackedEchelon:
    def __init__(self, space, *, max_row_terms=100_000, max_total_terms=1_000_000):
        if max_row_terms < 1 or max_total_terms < 1:
            raise ValueError('budgets must be positive')
        self.space = space
        self.max_row_terms, self.max_total_terms = max_row_terms, max_total_terms
        self.pivots = {}
        self.stored_terms = self.xor_steps = self.peak_row_terms = 0

    @property
    def rank(self):
        return len(self.pivots)

    @property
    def packed_payload_bytes(self):
        # Minimal per-row bit payload; excludes Python integer/object overhead.
        return sum((r.bits.bit_length() + 7) // 8 for r in self.pivots.values())

    def add(self, row):
        if not isinstance(row, PackedRow):
            raise TypeError('expected packed row')
        if len(row) > self.max_row_terms:
            raise ResourceLimitError('input row exceeds term budget')
        if row.bits.bit_length() > self.space.count:
            raise IndexError('column outside declared space')
        value = row.bits
        while value and value.bit_length() - 1 in self.pivots:
            value ^= self.pivots[value.bit_length() - 1].bits
            self.xor_steps += 1
            terms = value.bit_count()
            self.peak_row_terms = max(self.peak_row_terms, terms)
            if terms > self.max_row_terms:
                raise ResourceLimitError('elimination fill-in exceeded row budget')
        if not value:
            return False
        terms = value.bit_count()
        if self.stored_terms + terms > self.max_total_terms:
            raise ResourceLimitError('matrix term budget exceeded')
        self.peak_row_terms = max(self.peak_row_terms, terms)
        self.pivots[value.bit_length() - 1] = PackedRow(value)
        self.stored_terms += terms
        return True


spec = importlib.util.spec_from_file_location('packed_closure_engine', PRIOR / 'boolean_closure.py')
packed_engine = importlib.util.module_from_spec(spec)
spec.loader.exec_module(packed_engine)
packed_engine.BoundedMonomials = PackedSpace
packed_engine.SparseEchelon = PackedEchelon
packed_compute = packed_engine.compute


def compute(n, generators, *, backend='sparse', **kwargs):
    # The imported engine independently enforces 1 <= n <= 10.
    if backend == 'sparse':
        return sparse_compute(n, generators, **kwargs)
    if backend == 'packed':
        return packed_compute(n, generators, **kwargs)
    raise ValueError('unknown backend')


def semantic_result(result):
    """All outputs/counters that must be identical for a backend-only change."""
    return {**result, 'stats': {k: v for k, v in result['stats'].items()
                               if k != 'packed_payload_bytes'}}


def structure(n, generators):
    """Equation-level primal graph and a min-fill upper bound, not treewidth."""
    generators = normalize(n, generators)
    graph = {i: set() for i in range(n)}
    degrees = {}
    for row in generators:
        support = 0
        for mask in row:
            support |= mask
            degree = mask.bit_count()
            degrees[degree] = degrees.get(degree, 0) + 1
        variables = [i for i in range(n) if support & (1 << i)]
        for i in variables:
            graph[i].update(j for j in variables if j != i)
    edges = [[i, j] for i in graph for j in sorted(graph[i]) if i < j]
    order, fill, width = [], [], 0
    while graph:
        def missing(i):
            ns = sorted(graph[i])
            return [(a, b) for k, a in enumerate(ns) for b in ns[k+1:] if b not in graph[a]]
        vertex = min(graph, key=lambda i: (len(missing(i)), len(graph[i]), i))
        width = max(width, len(graph[vertex]))
        for a, b in missing(vertex):
            graph[a].add(b)
            graph[b].add(a)
            fill.append([a, b])
        for neighbor in graph[vertex]:
            graph[neighbor].remove(vertex)
        del graph[vertex]
        order.append(vertex)
    return {'input_terms': sum(map(len, generators)),
            'max_input_row_terms': max(map(len, generators), default=0),
            'input_degree_histogram': degrees, 'primal_edges': edges,
            'min_fill_order': order, 'min_fill_width_upper_bound': width,
            'fill_edges': fill}
