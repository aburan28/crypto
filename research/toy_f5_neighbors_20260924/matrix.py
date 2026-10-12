"""Bounded reference of the repository's Boolean matrix-F5 criterion.

Hard limit: six Boolean variables / one 64-bit row. Not a general solver.
"""
import hashlib


class Ring:
    def __init__(self, n):
        if not 1 <= n <= 6:
            raise ValueError('toy experiments require 1..6 variables')
        self.n = n
        self.monos = sorted(range(1 << n), key=lambda m: (m.bit_count(), m))
        self.index = {m: i for i, m in enumerate(self.monos)}

    def pack(self, terms):
        return sum(1 << self.index[m] for m in terms)

    def terms(self, row):
        return [m for i, m in enumerate(self.monos) if row >> i & 1]

    def degree(self, row):
        return max((m.bit_count() for m in self.terms(row)), default=0)

    def product(self, row, m):
        out = 0
        for t in self.terms(row):
            out ^= 1 << self.index[t | m]
        return out

    def multipliers(self, d):
        return [m for m in self.monos if m.bit_count() <= d]

    def evaluate(self, row, assignment):
        return sum((t & assignment) == t for t in self.terms(row)) % 2

    def anf(self, values):
        if len(values) != 1 << self.n:
            raise ValueError('wrong truth table size')
        coeff = list(values)
        for j in range(self.n):
            for a in range(1 << self.n):
                if a >> j & 1:
                    coeff[a] ^= coeff[a ^ (1 << j)]
        return self.pack(m for m, bit in enumerate(coeff) if bit)


class Echelon:
    def __init__(self, pivots=None):
        self.pivots = dict(pivots or {})
        self.xors = self.inserted = self.zeros = 0

    def reduce(self, row):
        while row:
            p = row.bit_length() - 1
            if p not in self.pivots:
                break
            row ^= self.pivots[p]
            self.xors += 1
        return row

    def insert(self, row):
        self.inserted += 1
        row = self.reduce(row)
        if row:
            self.pivots[row.bit_length() - 1] = row
        else:
            self.zeros += 1

    def rref(self):
        pivots = dict(self.pivots)
        for low in sorted(pivots):
            for high in sorted(pivots):
                if high > low and pivots[high] >> low & 1:
                    pivots[high] ^= pivots[low]
                    self.xors += 1
        return tuple(pivots[k] for k in sorted(pivots))


def criterion(ring, generators, degree):
    """LM(V_{i-1}(D-d_i) + <s(f_i+1): deg s <= D-2d_i>)."""
    forbidden = [set() for _ in generators]
    cost = dict(xors=0, products=0, rows=0)
    degrees = [ring.degree(f) for f in generators]
    if 0 in degrees:
        return forbidden, cost
    for e in sorted({degree - d for d in degrees if d <= degree}):
        prefix = Echelon()
        for i, (f, d) in enumerate(zip(generators, degrees)):
            if degree - d == e:
                extra = Echelon(prefix.pivots)
                for s in ring.multipliers(degree - 2*d):
                    extra.insert(ring.product(f ^ 1, s))
                    cost['products'] += 1
                forbidden[i] = {ring.monos[p] for p in extra.pivots}
                cost['xors'] += extra.xors
                cost['rows'] += extra.inserted
            for s in ring.multipliers(e - d):
                prefix.insert(ring.product(f, s))
                cost['products'] += 1
        cost['xors'] += prefix.xors
        cost['rows'] += prefix.inserted
    return forbidden, cost


def step(ring, generators, degree, f5):
    forbidden, cost = criterion(ring, generators, degree) if f5 else (
        [set() for _ in generators], dict(xors=0, products=0, rows=0))
    basis = Echelon()
    products = pruned = zero_products = 0
    for i, f in enumerate(generators):
        for s in ring.multipliers(degree - ring.degree(f)):
            if s in forbidden[i]:
                pruned += 1
                continue
            row = ring.product(f, s)
            products += 1
            if row:
                basis.insert(row)
            else:
                zero_products += 1
    rows = basis.rref()
    return rows, dict(degree=degree, rank=len(rows), matrix_rows=basis.inserted,
        pruned_rows=pruned, zero_products=zero_products, zero_reductions=basis.zeros,
        elimination_xors=basis.xors, criterion_xors=cost['xors'],
        total_matrix_xors=basis.xors+cost['xors'], products=products,
        criterion_products=cost['products'], criterion_rows=cost['rows'],
        rowspace_sha256=hashlib.sha256(repr(rows).encode()).hexdigest())


def complete(ring, rows, generators):
    checker = Echelon({r.bit_length()-1: r for r in rows})
    if any(checker.reduce(f) for f in generators):
        return False, checker.xors
    for row in rows:
        for j in range(ring.n):
            if checker.reduce(ring.product(row, 1 << j)):
                return False, checker.xors
    return True, checker.xors


def measure(ring, generators):
    generators = [f for f in generators if f]
    traces = {'f4': [], 'f5': []}
    for degree in range(2*ring.n+1):
        f4, a = step(ring, generators, degree, False)
        f5, b = step(ring, generators, degree, True)
        if f4 != f5:
            raise AssertionError('F5 changed the Macaulay row space')
        done, check_xors = complete(ring, f4, generators)
        for label, report in [('f4', a), ('f5', b)]:
            report['completion_check_xors'] = check_xors
            report['ideal_complete'] = done
            traces[label].append(report)
        if done:
            roots = [a for a in range(1 << ring.n)
                     if all(ring.evaluate(f, a) == 0 for f in f5)]
            direct = [a for a in range(1 << ring.n)
                      if all(ring.evaluate(f, a) == 0 for f in generators)]
            if roots != direct or len(f5) != (1 << ring.n) - len(roots):
                raise AssertionError('Boolean ideal/root certificate failed')
            return dict(completion_degree=degree, roots=roots, traces=traces,
                        status='VERIFIED', max_input_degree=max(map(ring.degree, generators), default=0))
    raise AssertionError('full bounded Boolean ideal should complete by 2*n')
