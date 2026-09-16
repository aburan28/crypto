"""Independent polynomial-basis binary-curve certificate checker (stdlib only)."""
import math
import hashlib
import json


class InvalidEvidence(ValueError):
    pass


def require(condition, message):
    if not condition:
        raise InvalidEvidence(message)


class Curve:
    def __init__(self, fixture):
        self.n = int(fixture['degree'])
        self.a = int(fixture['curve_a'])
        self.r = int(fixture['subgroup_order'])
        require(5 <= self.n <= 31 and self.n % 2 == 1, 'unsupported degree')
        require(self.a in (0, 1), 'unsupported coefficient')
        terms = fixture['irreducible']['low_terms']
        require(fixture['irreducible']['degree'] == self.n, 'field degree mismatch')
        require(len(set(terms)) == len(terms) and all(0 <= x < self.n for x in terms), 'bad modulus')
        self.modulus = (1 << self.n) | sum(1 << t for t in terms)
        require(self.r > 2 and all(self.r % d for d in range(2, math.isqrt(self.r) + 1)), 'nonprime subgroup')
        t = -1 if self.a == 0 else 1
        s0, s1 = 2, t
        for _ in range(1, self.n):
            s0, s1 = s1, t * s1 - 2 * s0
        order = (1 << self.n) + 1 - s1
        require(order == int(fixture['group_order']), 'wrong group order')
        self.h = int(fixture['cofactor'])
        require(order == self.h * self.r, 'wrong cofactor')
        self.g = self.decode(fixture['generator'])
        require(self.g is not None and self.mul(self.g, self.r) is None, 'bad generator')
        self.lam = int(fixture['lambda'])
        require(self.mul(self.g, self.lam) == self.frob(self.g), 'wrong Frobenius eigenvalue')

    def fm(self, a, b):
        z = 0
        while b:
            if b & 1:
                z ^= a
            b >>= 1
            a <<= 1
            if a >> self.n:
                a ^= self.modulus
        return z

    def inv(self, x):
        require(x != 0, 'zero denominator')
        u, v, a, b = x, self.modulus, 1, 0
        while u != 1:
            require(u != 0, 'reducible field modulus')
            shift = u.bit_length() - v.bit_length()
            if shift < 0:
                u, v, a, b = v, u, b, a
                shift = -shift
            u ^= v << shift
            a ^= b << shift
        while a.bit_length() > self.n:
            a ^= self.modulus << (a.bit_length() - self.n - 1)
        return a

    def decode(self, p):
        if p is None:
            return None
        require(isinstance(p, list) and len(p) == 2, 'malformed point')
        x, y = map(int, p)
        require(0 <= x < 1 << self.n and 0 <= y < 1 << self.n, 'point outside field')
        xx = self.fm(x, x)
        require(self.fm(y, y) ^ self.fm(x, y) == self.fm(xx, x) ^ self.fm(self.a, xx) ^ 1,
                'point does not lift to curve')
        return x, y

    def neg(self, p):
        return None if p is None else (p[0], p[0] ^ p[1])

    def frob(self, p):
        return None if p is None else (self.fm(p[0], p[0]), self.fm(p[1], p[1]))

    def add(self, p, q):
        if p is None:
            return q
        if q is None:
            return p
        x, y = p
        u, v = q
        if x == u:
            if y != v or x == 0:
                return None
            lam = x ^ self.fm(y, self.inv(x))
            z = self.fm(lam, lam) ^ lam ^ self.a
            return z, self.fm(x, x) ^ self.fm(lam ^ 1, z)
        lam = self.fm(y ^ v, self.inv(x ^ u))
        z = self.fm(lam, lam) ^ lam ^ x ^ u ^ self.a
        return z, self.fm(lam, x ^ z) ^ z ^ y

    def mul(self, p, k):
        require(k >= 0, 'negative scalar')
        q = None
        while k:
            if k & 1:
                q = self.add(q, p)
            p = self.add(p, p)
            k >>= 1
        return q


def rank(rows, columns, r):
    pivots = {}
    for original in rows:
        row = list(original)
        for j in range(columns):
            if row[j] == 0:
                continue
            if j in pivots:
                c = row[j]
                row = [(a - c*b) % r for a, b in zip(row, pivots[j])]
            else:
                inv = pow(row[j], -1, r)
                pivots[j] = [(a * inv) % r for a in row]
                break
    return len(pivots)


def verify(report, expected_fixture, *, expected_mode='ic', summands=3):
    require(report.get('schema_version') == 1, 'wrong schema')
    require(report.get('status') == 'complete', 'incomplete workload')
    require(report.get('mode') == expected_mode, 'wrong algorithm mode')
    require(report.get('fixture') == expected_fixture, 'changed curve or target')
    require(expected_fixture.get('target_scalar_constructed') is False, 'planted scalar supplied')
    c = Curve(expected_fixture)
    targets = [c.decode(p) for p in expected_fixture['targets']]
    require(all(q is not None and c.mul(q, c.r) is None for q in targets), 'bad subgroup target')
    solutions = report.get('solutions', [])
    require(len(solutions) == len(targets), 'missing or extra solutions')
    require({s.get('index') for s in solutions} == set(range(len(targets))), 'duplicate solution indices')
    for s in solutions:
        require(isinstance(s.get('recovered'), str), 'missing scalar')
        d = int(s['recovered'])
        require(0 <= d < c.r and c.mul(c.g, d) == targets[s['index']], 'incorrect scalar')
    if expected_mode == 'rho':
        require(report.get('automorphism_order') == 2*c.n, 'rho quotient mismatch')
        return {'verified_targets':len(targets), 'verified_relations':0, 'rank':None,
                'solutions':[s['recovered'] for s in sorted(solutions,key=lambda s:s['index'])]}
    require(report.get('rejected_relations') == 0, 'rejected relation in run')
    require(report.get('summands') == summands, 'changed summand count')
    base = [c.decode(p) for p in report['factor_base']]
    require(base and all(p is not None for p in base) and len(set(base)) == len(base), 'invalid factor base')
    logs = report.get('column_logs', [])
    require(len(logs) == report.get('columns') and len(logs) > 1, 'missing/nontrivial log matrix')
    mapping = {}
    for j, entry in enumerate(logs):
        p = c.decode(entry['point'])
        ell = int(entry['log'])
        require(0 <= ell < c.r and p is not None and c.mul(c.g, ell) == p, 'bad column log')
        coeff = 1
        for _ in range(c.n):
            for q, k in ((p, coeff), (c.neg(p), -coeff % c.r)):
                require(q not in mapping or mapping[q] == (j, k), 'overlapping projected columns')
                mapping[q] = j, k
            p = c.frob(p)
            coeff = coeff*c.lam % c.r
    projected = [mapping.get(c.mul(p, c.h)) for p in base]
    # An identity projection contributes zero; every other projection needs a column.
    require(all(x is not None or c.mul(p,c.h) is None for p,x in zip(base,projected)), 'uncovered base column')
    matrix = []
    seen = set()
    rows = report.get('relations', [])
    for rel in rows:
        ids = rel['points']
        require(len(ids) == summands and all(type(i) is int and 0 <= i < len(base) for i in ids), 'bad relation indices')
        a = int(rel['a'])
        require(0 < a < c.r, 'trivial or invalid relation scalar')
        q = None
        row = [0] * len(logs)
        for i in ids:
            q = c.add(q, base[i])
            if projected[i] is not None:
                j, k = projected[i]
                row[j] = (row[j]+k) % c.r
        require(q is not None and q == c.mul(c.g,a), 'incorrect point relation')
        require(sum(x*int(l['log']) for x,l in zip(row,logs)) % c.r == c.h*a % c.r, 'incorrect scalar-field row')
        key = a, tuple(sorted(ids))
        if key not in seen:
            matrix.append(row)
            seen.add(key)
    achieved_rank = rank(matrix, len(logs), c.r)
    require(achieved_rank == len(logs), 'rank-deficient log recovery')
    require(report.get('accepted_relations') == len(seen), 'accepted-row accounting mismatch')
    require(report.get('duplicate_relations') == len(rows)-len(seen), 'duplicate-row accounting mismatch')
    # Index-calculus admission: every target's logarithm must be derived from one
    # relation [a]G + [b]Q = sum of factor-base points, verified in the group, with
    # the scalar the consequence of that relation under the verified column logs.
    degenerate = 0
    for s in solutions:
        rel = s.get('relation')
        require(isinstance(rel, dict), 'missing descent relation: logarithm not certified as index calculus')
        a, b, ids = rel.get('a'), rel.get('b'), rel.get('points')
        require(type(a) is int and type(b) is int and 0 <= a < c.r and 0 < b < c.r, 'invalid descent scalars')
        require(isinstance(ids, list) and len(ids) in (0, summands)
                and all(type(i) is int and 0 <= i < len(base) for i in ids), 'bad descent relation indices')
        q = targets[s['index']]
        probe = c.add(c.mul(c.g, a), c.mul(q, b))
        total = None
        logsum = 0
        for i in ids:
            total = c.add(total, base[i])
            if projected[i] is not None:
                j, k = projected[i]
                logsum = (logsum + k*int(logs[j]['log'])) % c.r
        require(total == probe, 'descent relation does not hold in the group')
        require(c.h*(a + b*int(s['recovered'])) % c.r == logsum, 'logarithm is not the consequence of its descent relation')
        degenerate += not ids
    return {'verified_targets':len(targets), 'verified_relations':len(rows), 'fresh_rows':len(seen),
            'certified_descents':len(solutions), 'degenerate_descents':degenerate,
            'rank':achieved_rank, 'signed_base_size':len(base),
            'factor_base_sha256':hashlib.sha256(json.dumps(sorted(base),separators=(',',':')).encode()).hexdigest(),
            'solutions':[s['recovered'] for s in sorted(solutions,key=lambda s:s['index'])]}
