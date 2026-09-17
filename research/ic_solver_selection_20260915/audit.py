"""Read-only independent polynomial-basis audit of selector query evidence."""
import hashlib
import json
from pathlib import Path

from oracle_pb import Curve, rank, require


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


class Checker:
    def __init__(self, fixture):
        self.fixture = fixture
        self.degree = int(fixture['parameters']['curve']['degree'])
        self.prime = int(fixture['prime'])
        self.summands = int(fixture['parameters']['summands'])
        self.basis = {}
        current = (1 << self.degree) - 1
        modulus = None
        for k in range(self.degree + 1):
            vector, tag = current, 1 << k
            while vector:
                pivot = vector.bit_length() - 1
                if pivot not in self.basis:
                    self.basis[pivot] = vector, tag
                    break
                value, coefficient = self.basis[pivot]
                vector ^= value
                tag ^= coefficient
            if not vector:
                require(k == self.degree, 'normal-basis conversion has low degree')
                modulus = tag
            current = self.normal_mul(current, 1)
        require(modulus is not None and len(self.basis) == self.degree, 'singular normal basis')
        spec = {'degree': self.degree, 'curve_a': 0, 'subgroup_order': str(self.prime),
                'cofactor': 4, 'group_order': str(4 * self.prime),
                'irreducible': {'degree': self.degree,
                                'low_terms': [i for i in range(self.degree) if modulus >> i & 1]},
                'generator': self.convert(fixture['generator_onb']), 'lambda': fixture['eigen']}
        self.curve = Curve(spec)
        self.base = [self.curve.decode(self.convert(p)) for p in fixture['base']]
        reps = [self.curve.decode(self.convert(p)) for p in fixture['representatives']]
        require(len(set(self.base)) == len(self.base), 'duplicate base points')
        self.mapping = {}
        for column, point in enumerate(reps):
            require(self.curve.mul(point, self.prime) is None, 'base point outside subgroup')
            coefficient = 1
            for _ in range(self.degree):
                for p, a in ((point, coefficient), (self.curve.neg(point), -coefficient % self.prime)):
                    require(p not in self.mapping or self.mapping[p] == (column, a), 'overlapping orbits')
                    self.mapping[p] = column, a
                point = self.curve.frob(point)
                coefficient = coefficient * self.curve.lam % self.prime
        require(set(self.base) == set(self.mapping), 'incomplete base orbit expansion')
        self.columns = len(reps)

    def normal_mul(self, a, b):
        # gamma_i gamma_j = gamma_(i+j) + gamma_(i-j), gamma_0 = 0.
        # This direct coordinate implementation shares no arithmetic with the engine.
        result = 0
        for i in range(1, self.degree + 1):
            if not (a >> (i - 1) & 1):
                continue
            for j in range(1, self.degree + 1):
                if b >> (j - 1) & 1:
                    for index in (i + j, i - j):
                        index %= 2 * self.degree + 1
                        index = min(index, 2 * self.degree + 1 - index)
                        if index:
                            result ^= 1 << (index - 1)
        return result

    def coordinate(self, value):
        require(type(value) is int and 0 <= value < 1 << self.degree, 'bad normal coordinate')
        result = 0
        while value:
            pivot = value.bit_length() - 1
            v, tag = self.basis[pivot]
            value ^= v
            result ^= tag
        return result

    def convert(self, point):
        return None if point is None else [self.coordinate(x) for x in point]

    def relation(self, target, witness, expected):
        require(isinstance(witness, list) and len(witness) == self.summands
                and all(type(i) is int and 0 <= i < len(self.base) for i in witness), 'bad witness')
        actual = [self.base[i] for i in witness]
        require(all(p != self.curve.neg(q) for i, p in enumerate(actual) for q in actual[:i]),
                'inverse cancellation')
        total, row = None, [0] * self.columns
        for p in actual:
            total = self.curve.add(total, p)
            column, coefficient = self.mapping[p]
            row[column] = (row[column] + coefficient) % self.prime
        require(total == target and row == expected, 'incorrect decomposition or transported row')
        return row

    def query(self, record):
        target = self.curve.decode(self.convert(record['target']))
        require(target is not None and self.curve.mul(target, self.prime) is None, 'invalid query target')
        # Independently compute the orbit in the normal coordinates for its split ID.
        members = []
        current = record['target']
        for _ in range(self.degree):
            x, y = current
            members.extend(([x, y], [x, x ^ y]))
            current = [self.normal_mul(z, z) for z in current]
        orbit_id = digest([self.degree, str(self.prime), min(members)])
        require(orbit_id == record['orbit'], 'wrong canonical orbit')
        split = ('train', 'train', 'train', 'validation', 'confirmation')[int(orbit_id, 16) % 5]
        require(record['split'] == split, 'orbit leakage')
        direct = record['features']['direct_target']
        require(direct in (0, 1), 'invalid query phase')
        rows = []
        for prior in record['prior_relations']:
            point = self.curve.decode(self.convert(prior['target']))
            row = self.relation(point, prior['witness'], prior['row'])
            rows.append(row + ([0] if direct else []))
        columns = self.columns + direct
        before = rank(rows, columns, self.prime)
        require(before == record['rank_before'], 'incorrect snapshot rank')
        x = record['target'][0]
        require(record['features'] == {'x_weight': x.bit_count(),
                'x_frobenius_distance': (x ^ self.normal_mul(x, x)).bit_count(),
                'rank_fraction': before / columns, 'direct_target': direct}, 'post-query or incorrect features')
        verified = 0
        for result in record['outcomes'].values():
            require(result['status'] != 'error', 'retained solver error prevents admission')
            require(result['elapsed_ns'] > 0, 'missing query cost')
            if result['verified']:
                row = self.relation(target, result['witness'], result['row'])
                gain = int(rank(rows + [row + ([-1 % self.prime] if direct else [])], columns, self.prime) > before)
                require(gain == result['gain'], 'incorrect useful-rank label')
                verified += 1
            else:
                require(result['gain'] == 0 and result['witness'] is None, 'unverified positive label')
                if result['status'] == 'unsat':
                    # Exhaustive small pair/triple-sum oracle with the same cancellation convention.
                    require(not self.reachable(target), 'false no-decomposition certificate')
        return verified

    def reachable(self, target):
        c = self.curve
        if not hasattr(self, 'pair_sums'):
            self.pair_sums = {}
            for i, p in enumerate(self.base):
                for q in self.base[i:]:
                    if p != c.neg(q):
                        self.pair_sums.setdefault(c.add(p, q), []).append((p, q))
        if self.summands == 2:
            return target in self.pair_sums
        for p in self.base:
            for q, r in self.pair_sums.get(c.add(target, c.neg(p)), []):
                if p != c.neg(q) and p != c.neg(r):
                    return True
        return False


def audit(directory):
    path = Path(directory)
    checker = Checker(json.loads((path / 'fixture.json').read_text()))
    queries = verified = 0
    origins = {}
    for line in (path / 'queries.jsonl').read_text().splitlines():
        record = json.loads(line)
        verified += checker.query(record)
        queries += 1
        origins[record['origin']] = origins.get(record['origin'], 0) + 1
    return {'directory': str(path), 'queries': queries, 'verified_policy_relations': verified,
            'origins': origins, 'independent_polynomial_basis_checks': 'passed'}


if __name__ == '__main__':
    import sys
    print(json.dumps([audit(p) for p in sys.argv[1:]], indent=2))
