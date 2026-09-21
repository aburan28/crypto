"""Exact radical unordered-pair systems for finite index-calculus controls.

This module builds public synthetic polynomial systems only.  It neither accepts
an unknown scalar nor implements discrete-log recovery.
"""

from __future__ import annotations

import hashlib
import itertools
import json
import time

from sage.all import EllipticCurve, GF, PolynomialRing, matrix, prod, vector


def _stable_hash(value) -> str:
    blob = json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(blob).hexdigest()


def system_hash(polynomials) -> str:
    ring = polynomials[0].parent()
    return _stable_hash({
        "base": str(ring.base_ring()),
        "order": str(ring.term_order()),
        "variables": list(ring.variable_names()),
        "polynomials": [str(f) for f in polynomials],
    })


def s3_prime(x, y, z, curve_a, curve_b):
    return ((x - y) ** 2 * z**2
            - 2 * (x + y) * (x * y + curve_a) * z
            - 4 * curve_b * z
            + (x * y - curve_a) ** 2
            - 4 * curve_b * (x + y))


def s3_prime_symmetric(s, t, z, curve_a, curve_b):
    return ((s * s - 4 * t) * z**2
            - 2 * s * (t + curve_a) * z
            - 4 * curve_b * z
            + (t - curve_a) ** 2
            - 4 * curve_b * s)


def norm_pair_membership(ring, s, t, roots):
    extension = PolynomialRing(ring, "X")
    x = extension.gen()
    root_poly = prod(x - extension(root) for root in roots)
    remainder = root_poly.mod(x**2 - extension(s) * x + extension(t))
    r0, r1 = ring(remainder[0]), ring(remainder[1])
    return [r1 * s + 2 * r0,
            r1 * r1 * t + r1 * r0 * s + r0 * r0]


def radical_pair_basis(field, roots, *, exhaustive=True):
    """Return the DRL basis of {(x+y,xy): x,y in roots} and phase timings."""
    started = time.perf_counter()
    elimination_ring = PolynomialRing(field, names=("x", "y", "s", "t"), order="lex")
    x, y, s, t = elimination_ring.gens()
    root_poly = lambda value: prod(value - root for root in roots)
    graph = elimination_ring.ideal([
        root_poly(x), root_poly(y), s - x - y, t - x * y,
    ])
    lex_basis = list(graph.groebner_basis(algorithm="libsingular:slimgb"))
    eliminated = [g for g in lex_basis if g.degree(x) == 0 and g.degree(y) == 0]
    lex_seconds = time.perf_counter() - started

    target_ring = PolynomialRing(field, names=("s", "t"), order="degrevlex")
    target_s, target_t = target_ring.gens()
    started = time.perf_counter()
    image_ideal = target_ring.ideal([
        target_ring(g(x=0, y=0, s=target_s, t=target_t)) for g in eliminated
    ])
    basis = list(image_ideal.groebner_basis(algorithm="libsingular:slimgb"))
    drl_seconds = time.perf_counter() - started

    image = {(a + b, a * b) for a in roots for b in roots}
    if exhaustive:
        zeros = {(a, b) for a in field for b in field
                 if all(g(a, b) == 0 for g in basis)}
        assert zeros == image
    assert image_ideal.vector_space_dimension() == len(image)
    assert len(basis) == len(roots) + 1
    assert {int(g.total_degree()) for g in basis} == {len(roots)}

    return target_ring, basis, {
        "lex_elimination_seconds": lex_seconds,
        "drl_conversion_seconds": drl_seconds,
        "total_seconds": lex_seconds + drl_seconds,
        "basis_size": len(basis),
        "basis_degree": len(roots),
        "basis_terms": sum(len(g.dict()) for g in basis),
        "image_size": len(image),
        "exhaustive": exhaustive,
    }


def ordered_quadratic_roots(field, s, t):
    ring = PolynomialRing(field, "X")
    x = ring.gen()
    expanded = []
    for root, multiplicity in (x**2 - s * x + t).roots():
        expanded.extend([root] * multiplicity)
    assert len(expanded) == 2
    if expanded[0] == expanded[1]:
        return [(expanded[0], expanded[1])]
    return [(expanded[0], expanded[1]), (expanded[1], expanded[0])]


class PrimeCell:
    def __init__(self, p, curve_a, curve_b, factor_base_size):
        self.field = GF(p)
        self.a = self.field(curve_a)
        self.b = self.field(curve_b)
        self.curve = EllipticCurve(self.field, [self.a, self.b])
        rational_x = [x for x in self.field
                      if any(y * y == x**3 + self.a * x + self.b for y in self.field)]
        self.roots = rational_x[:factor_base_size]
        assert len(self.roots) == factor_base_size
        assert not any(x**3 + self.a * x + self.b == 0 for x in self.roots)
        self.radical_ring, self.radical_basis, self.preprocessing = radical_pair_basis(
            self.field, self.roots)
        self.pair_image = {(x + y, x * y) for x in self.roots for y in self.roots}
        self.pair_u = {
            (x, y): [u for u in self.field if s3_prime(x, y, u, self.a, self.b) == 0]
            for x in self.roots for y in self.roots
        }
        self.image_u = {
            (s, t): [u for u in self.field
                     if s3_prime_symmetric(s, t, u, self.a, self.b) == 0]
            for s, t in self.pair_image
        }
        assert all(self.pair_u[(x, y)] == self.image_u[(x + y, x * y)]
                   for x in self.roots for y in self.roots)
        self.lifts = {
            x: [self.curve(x, y) for y in self.field
                if y * y == x**3 + self.a * x + self.b]
            for x in self.roots
        }

    @property
    def cell_id(self):
        return f"prime-p{self.field.characteristic()}-b{len(self.roots)}"

    def direct_system(self, target):
        ring = PolynomialRing(
            self.field, names=("x1", "x2", "x3", "x4", "u", "v"),
            order="degrevlex")
        x1, x2, x3, x4, u, v = ring.gens()
        membership = [prod(x - root for root in self.roots)
                      for x in (x1, x2, x3, x4)]
        return membership + [
            s3_prime(x1, x2, u, self.a, self.b),
            s3_prime(x3, x4, v, self.a, self.b),
            s3_prime(u, v, self.field(target), self.a, self.b),
        ]

    def norm_system(self, target):
        names = ("u", "v", "t1", "t2", "s1", "s2")
        ring = PolynomialRing(self.field, names=names, order="degrevlex")
        values = ring.gens_dict()
        u, v, t1, t2, s1, s2 = [values[name] for name in names]
        return (norm_pair_membership(ring, s1, t1, self.roots)
                + norm_pair_membership(ring, s2, t2, self.roots)
                + [s3_prime_symmetric(s1, t1, u, self.a, self.b),
                   s3_prime_symmetric(s2, t2, v, self.a, self.b),
                   s3_prime(u, v, self.field(target), self.a, self.b)])

    def radical_system(self, target):
        names = ("t1", "t2", "s1", "s2", "u", "v")
        ring = PolynomialRing(self.field, names=names, order="degrevlex")
        values = ring.gens_dict()
        t1, t2, s1, s2, u, v = [values[name] for name in names]
        embed = lambda polynomial, s, t: ring(polynomial(s, t))
        return ([embed(g, s1, t1) for g in self.radical_basis]
                + [embed(g, s2, t2) for g in self.radical_basis]
                + [s3_prime_symmetric(s1, t1, u, self.a, self.b),
                   s3_prime_symmetric(s2, t2, v, self.a, self.b),
                   s3_prime(u, v, self.field(target), self.a, self.b)])

    def systems(self, target):
        return {
            "direct": self.direct_system(target),
            "norm_symmetric": self.norm_system(target),
            "radical_symmetric": self.radical_system(target),
        }

    def certificate(self, target):
        target = self.field(target)
        ordered = []
        for x1, x2, x3, x4 in itertools.product(self.roots, repeat=4):
            for u in self.pair_u[(x1, x2)]:
                for v in self.pair_u[(x3, x4)]:
                    if s3_prime(u, v, target, self.a, self.b) == 0:
                        ordered.append((x1, x2, x3, x4, u, v))
        symmetric = {
            (x1 + x2, x1 * x2, x3 + x4, x3 * x4, u, v)
            for x1, x2, x3, x4, u, v in ordered
        }
        independently_symmetric = {
            (s1, t1, s2, t2, u, v)
            for s1, t1 in self.pair_image
            for s2, t2 in self.pair_image
            for u in self.image_u[(s1, t1)]
            for v in self.image_u[(s2, t2)]
            if s3_prime(u, v, target, self.a, self.b) == 0
        }
        assert symmetric == independently_symmetric
        recovered = set()
        for s1, t1, s2, t2, u, v in symmetric:
            for x1, x2 in ordered_quadratic_roots(self.field, s1, t1):
                for x3, x4 in ordered_quadratic_roots(self.field, s2, t2):
                    recovered.add((x1, x2, x3, x4, u, v))
        assert recovered == set(ordered)

        chain_x = {row[:4] for row in ordered}
        finite_group_x = set()
        full_group_x = set()
        for xs in itertools.product(self.roots, repeat=4):
            for points in itertools.product(*(self.lifts[x] for x in xs)):
                left = points[0] + points[1]
                right = points[2] + points[3]
                total = left + right
                if total != self.curve(0) and total[0] == target:
                    full_group_x.add(xs)
                    if left != self.curve(0) and right != self.curve(0):
                        finite_group_x.add(xs)
        assert chain_x == finite_group_x
        normalized = [[int(x) for x in row] for row in sorted(set(ordered))]
        return {
            "ordered_roots": len(ordered),
            "symmetric_roots": len(symmetric),
            "chain_x": len(chain_x),
            "nondegenerate_group_x": len(finite_group_x),
            "full_group_x": len(full_group_x),
            "full_group_coverage": len(chain_x) / len(full_group_x) if full_group_x else 1.0,
            "root_set_sha256": _stable_hash(normalized),
            "presentation_equivalence": True,
            "root_recovery": True,
            "nondegenerate_group_equivalence": True,
        }


def _monomials(ring, degree):
    output = [ring.one()]
    for d in range(1, degree + 1):
        output.extend(prod(ring.gen(i) for i in indices)
                      for indices in itertools.combinations(range(ring.ngens()), d))
    return sorted(output, reverse=True)


def _boolean_image_basis(ring, points):
    fields = [x * x + x for x in ring.gens()]
    for degree in range(1, ring.ngens() + 1):
        monomials = _monomials(ring, degree)
        evaluation = matrix(ring.base_ring(), [[m(*point) for m in monomials]
                                               for point in points])
        relations = [sum((c * m for c, m in zip(row, monomials) if c), ring.zero())
                     for row in evaluation.right_kernel().basis()]
        ideal = ring.ideal(fields + relations)
        if ideal.dimension() == 0 and ideal.vector_space_dimension() == len(points):
            return list(ideal.groebner_basis(algorithm="libsingular:slimgb"))
    raise RuntimeError("Boolean image basis did not close")


def _affine_coordinates(points, base_field):
    origin = vector(base_field, points[0])
    basis = matrix(base_field, [vector(base_field, p) - origin for p in points]) \
        .row_space().basis_matrix()
    coordinates = [tuple(basis.transpose().solve_right(vector(base_field, p) - origin))
                   for p in points]
    return origin, basis, coordinates


class FieldVectors:
    def __init__(self, field):
        self.field = field
        self.base = GF(2)
        self.n = field.degree()
        _, self.from_vector, self.to_vector = field.vector_space(map=True)
        units = [vector(self.base, [int(i == j) for i in range(self.n)])
                 for j in range(self.n)]
        elements = [self.from_vector(unit) for unit in units]
        self.mul_tensor = [[[self.base(0) for _ in range(self.n)]
                            for _ in range(self.n)] for _ in range(self.n)]
        for i in range(self.n):
            for j in range(self.n):
                coordinates = self.to_vector(elements[i] * elements[j])
                for k in range(self.n):
                    self.mul_tensor[i][j][k] = coordinates[k]
        self.square_matrix = matrix(self.base, [self.to_vector(x * x) for x in elements])

    def const(self, ring, value):
        return [ring(c) for c in self.to_vector(value)]

    @staticmethod
    def add(a, b):
        return [x + y for x, y in zip(a, b)]

    def mul(self, a, b):
        ring = a[0].parent()
        output = [ring.zero() for _ in range(self.n)]
        for i, x in enumerate(a):
            if not x:
                continue
            for j, y in enumerate(b):
                if not y:
                    continue
                for k, coefficient in enumerate(self.mul_tensor[i][j]):
                    if coefficient:
                        output[k] += x * y
        return output

    def square(self, value):
        ring = value[0].parent()
        return [sum((self.square_matrix[i, k] * value[i]
                     for i in range(self.n) if self.square_matrix[i, k]), ring.zero())
                for k in range(self.n)]

    def s3(self, x, y, z, a6):
        xx, yy, zz = self.square(x), self.square(y), self.square(z)
        return self.add(self.add(self.add(self.mul(self.add(xx, yy), zz),
                                           self.mul(self.mul(x, y), z)),
                                   self.mul(xx, yy)), a6)

    def s3_symmetric(self, s, t, z, a6):
        return self.add(self.add(self.add(self.mul(self.square(s), self.square(z)),
                                           self.mul(t, z)), self.square(t)), a6)


def s3_binary_numeric(x, y, z, a6):
    return (x * x + y * y) * z * z + x * y * z + x * x * y * y + a6


def s3_binary_symmetric_numeric(s, t, z, a6):
    return s * s * z * z + t * z + t * t + a6


class BinaryCell:
    def __init__(self, k):
        self.k = k
        self.field = GF(2 ** (2 * k), "a")
        self.base = GF(2)
        self.vectors = FieldVectors(self.field)
        self.a6 = self.field(1)
        self.curve = EllipticCurve(self.field, [1, 0, 0, 0, self.a6])
        subfield = [x for x in self.field if x ** (2**k) == x]
        self.roots = [x for x in subfield if x != 0]
        points = [list(self.vectors.to_vector(x)) for x in subfield]
        origin, self.subfield_basis, _ = _affine_coordinates(points, self.base)
        assert not any(origin) and self.subfield_basis.nrows() == k
        self.coordinates = {
            x: tuple(self.subfield_basis.transpose().solve_right(self.vectors.to_vector(x)))
            for x in subfield
        }
        image_points = sorted({self.coordinates[x + y] + self.coordinates[x * y]
                               for x in self.roots for y in self.roots})
        local_ring = PolynomialRing(self.base, names=[f"q{i}" for i in range(2 * k)],
                                    order="degrevlex")
        started = time.perf_counter()
        self.local_basis = _boolean_image_basis(local_ring, image_points)
        self.preprocessing = {
            "total_seconds": time.perf_counter() - started,
            "basis_size": len(self.local_basis),
            "basis_terms": sum(len(g.dict()) for g in self.local_basis),
            "basis_degrees": [int(g.total_degree()) for g in self.local_basis],
            "image_size": len(image_points),
        }
        self.pair_image = {(x + y, x * y) for x in self.roots for y in self.roots}
        self.pair_u = {(x, y): [u for u in self.field
                                if s3_binary_numeric(x, y, u, self.a6) == 0]
                       for x in self.roots for y in self.roots}
        self.image_u = {(s, t): [u for u in self.field
                                 if s3_binary_symmetric_numeric(s, t, u, self.a6) == 0]
                        for s, t in self.pair_image}
        assert all(self.pair_u[(x, y)] == self.image_u[(x + y, x * y)]
                   for x in self.roots for y in self.roots)
        self.lifts = {x: [self.curve(x, y) for y in self.field
                          if y * y + x * y == x**3 + self.a6]
                      for x in self.roots}

    @property
    def cell_id(self):
        return f"binary-k{self.k}-b{len(self.roots)}"

    def _subfield_vector(self, variables, ring):
        return [sum((variables[i] * self.subfield_basis[i, j]
                     for i in range(self.k) if self.subfield_basis[i, j]), ring.zero())
                for j in range(2 * self.k)]

    def direct_system(self, target):
        names = ([f"x{block}_{j}" for block in range(1, 5) for j in range(self.k)]
                 + [f"u{j}" for j in range(2 * self.k)]
                 + [f"v{j}" for j in range(2 * self.k)])
        ring = PolynomialRing(self.base, names=names, order="degrevlex")
        values = ring.gens_dict()
        def xvec(block):
            return self._subfield_vector([values[f"x{block}_{j}"] for j in range(self.k)], ring)
        x1, x2, x3, x4 = [xvec(block) for block in range(1, 5)]
        u = [values[f"u{j}"] for j in range(2 * self.k)]
        v = [values[f"v{j}"] for j in range(2 * self.k)]
        polynomials = [x * x + x for x in ring.gens()]
        for block in range(1, 5):
            polynomials.append(prod(ring.one() + values[f"x{block}_{j}"]
                                    for j in range(self.k)))
        a6 = self.vectors.const(ring, self.a6)
        z = self.vectors.const(ring, self.field(target))
        polynomials += self.vectors.s3(x1, x2, u, a6)
        polynomials += self.vectors.s3(x3, x4, v, a6)
        polynomials += self.vectors.s3(u, v, z, a6)
        return [f for f in dict.fromkeys(polynomials) if f]

    def radical_system(self, target):
        names = ([f"t1_{j}" for j in range(self.k)]
                 + [f"t2_{j}" for j in range(self.k)]
                 + [f"s1_{j}" for j in range(self.k)]
                 + [f"s2_{j}" for j in range(self.k)]
                 + [f"u{j}" for j in range(2 * self.k)]
                 + [f"v{j}" for j in range(2 * self.k)])
        ring = PolynomialRing(self.base, names=names, order="degrevlex")
        values = ring.gens_dict()
        bits = lambda prefix: [values[f"{prefix}_{j}"] for j in range(self.k)]
        s1_bits, t1_bits = bits("s1"), bits("t1")
        s2_bits, t2_bits = bits("s2"), bits("t2")
        s1, t1 = self._subfield_vector(s1_bits, ring), self._subfield_vector(t1_bits, ring)
        s2, t2 = self._subfield_vector(s2_bits, ring), self._subfield_vector(t2_bits, ring)
        u = [values[f"u{j}"] for j in range(2 * self.k)]
        v = [values[f"v{j}"] for j in range(2 * self.k)]
        pair1 = [ring(f(*(s1_bits + t1_bits))) for f in self.local_basis]
        pair2 = [ring(f(*(s2_bits + t2_bits))) for f in self.local_basis]
        polynomials = list(dict.fromkeys([x * x + x for x in ring.gens()] + pair1 + pair2))
        a6 = self.vectors.const(ring, self.a6)
        z = self.vectors.const(ring, self.field(target))
        polynomials += self.vectors.s3_symmetric(s1, t1, u, a6)
        polynomials += self.vectors.s3_symmetric(s2, t2, v, a6)
        polynomials += self.vectors.s3(u, v, z, a6)
        return [f for f in dict.fromkeys(polynomials) if f]

    def systems(self, target):
        return {"direct": self.direct_system(target),
                "radical_symmetric": self.radical_system(target)}

    def certificate(self, target):
        target = self.field(target)
        ordered = []
        for xs in itertools.product(self.roots, repeat=4):
            for u in self.pair_u[(xs[0], xs[1])]:
                for v in self.pair_u[(xs[2], xs[3])]:
                    if s3_binary_numeric(u, v, target, self.a6) == 0:
                        ordered.append(xs + (u, v))
        symmetric = {(x1 + x2, x1 * x2, x3 + x4, x3 * x4, u, v)
                     for x1, x2, x3, x4, u, v in ordered}
        independently_symmetric = {
            (s1, t1, s2, t2, u, v)
            for s1, t1 in self.pair_image
            for s2, t2 in self.pair_image
            for u in self.image_u[(s1, t1)]
            for v in self.image_u[(s2, t2)]
            if s3_binary_numeric(u, v, target, self.a6) == 0
        }
        assert symmetric == independently_symmetric
        recovered = set()
        for s1, t1, s2, t2, u, v in symmetric:
            for x1, x2 in ordered_quadratic_roots(self.field, s1, t1):
                for x3, x4 in ordered_quadratic_roots(self.field, s2, t2):
                    recovered.add((x1, x2, x3, x4, u, v))
        assert recovered == set(ordered)

        chain_x = {row[:4] for row in ordered}
        finite_group_x, full_group_x = set(), set()
        for xs in itertools.product(self.roots, repeat=4):
            for points in itertools.product(*(self.lifts[x] for x in xs)):
                left, right = points[0] + points[1], points[2] + points[3]
                total = left + right
                if total != self.curve(0) and total[0] == target:
                    full_group_x.add(xs)
                    if left != self.curve(0) and right != self.curve(0):
                        finite_group_x.add(xs)
        assert chain_x == finite_group_x
        normalized = [[str(x) for x in row] for row in sorted(ordered, key=str)]
        return {
            "ordered_roots": len(ordered),
            "symmetric_roots": len(symmetric),
            "chain_x": len(chain_x),
            "nondegenerate_group_x": len(finite_group_x),
            "full_group_x": len(full_group_x),
            "full_group_coverage": len(chain_x) / len(full_group_x) if full_group_x else 1.0,
            "root_set_sha256": _stable_hash(normalized),
            "presentation_equivalence": True,
            "root_recovery": True,
            "nondegenerate_group_equivalence": True,
        }
