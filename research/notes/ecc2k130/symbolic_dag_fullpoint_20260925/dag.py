#!/usr/bin/env python3
"""Solver-neutral Boolean DAG for the complete affine K0 group relation.

This is an experimental expression system, not a SAT/ANF exporter.  Every
Boolean input has its own integer identifier.  A model is decoded from an
exact-length sequence of 64-bit limbs so that widths above 64 cannot wrap.
"""
from __future__ import annotations

from dataclasses import dataclass
def _poly_rem(a: int, modulus: int) -> int:
    while a and a.bit_length() >= modulus.bit_length():
        a ^= modulus << (a.bit_length() - modulus.bit_length())
    return a


def _poly_gcd(a: int, b: int) -> int:
    while b:
        a, b = b, _poly_rem(a, b)
    return a


def _square_mod(a: int, modulus: int) -> int:
    out = 0
    while a:
        low = a & -a
        out ^= 1 << (2 * (low.bit_length() - 1))
        a ^= low
    return _poly_rem(out, modulus)


def _prime_divisors(n: int) -> list[int]:
    answer = []
    d = 2
    while d * d <= n:
        if n % d == 0:
            answer.append(d)
            while n % d == 0:
                n //= d
        d += 1
    if n > 1:
        answer.append(n)
    return answer


def check_field(n: int, modulus: int) -> None:
    """Rabin irreducibility check, including canonical polynomial encoding."""
    if not isinstance(n, int) or not 1 <= n <= 4096:
        raise ValueError("field width out of supported range")
    if not isinstance(modulus, int) or modulus.bit_length() != n + 1 or not modulus & 1:
        raise ValueError("modulus must be monic, degree n, and have constant 1")
    x = _poly_rem(2, modulus)
    for divisor in _prime_divisors(n):
        power = x
        for _ in range(n // divisor):
            power = _square_mod(power, modulus)
        if _poly_gcd(power ^ x, modulus) != 1:
            raise ValueError("reducible field modulus")
    power = x
    for _ in range(n):
        power = _square_mod(power, modulus)
    if power != x:
        raise ValueError("reducible field modulus")


@dataclass(frozen=True)
class PackedModel:
    bit_count: int
    limbs: tuple[int, ...]

    @classmethod
    def from_bits(cls, bits: list[int]) -> "PackedModel":
        if any(bit not in (0, 1) for bit in bits):
            raise ValueError("Boolean model contains a non-bit")
        limbs = [0] * ((len(bits) + 63) // 64)
        for i, bit in enumerate(bits):
            limbs[i // 64] |= bit << (i % 64)
        return cls(len(bits), tuple(limbs))

    def validate(self, expected_bits: int) -> None:
        if self.bit_count != expected_bits:
            raise ValueError("Boolean model width mismatch")
        if len(self.limbs) != (expected_bits + 63) // 64:
            raise ValueError("Boolean model limb count mismatch")
        if any(not isinstance(limb, int) or not 0 <= limb < 1 << 64 for limb in self.limbs):
            raise ValueError("Boolean model has an invalid limb")
        if expected_bits % 64 and self.limbs[-1] >> (expected_bits % 64):
            raise ValueError("Boolean model has nonzero unused high bits")

    def bit(self, i: int) -> int:
        return (self.limbs[i // 64] >> (i % 64)) & 1


class Dag:
    """Hash-consed, topologically ordered constants, inputs, XORs and ANDs."""

    def __init__(self) -> None:
        self.nodes: list[tuple[str, int, int]] = [("zero", 0, 0), ("one", 0, 0)]
        self.cache: dict[tuple[str, int, int], int] = {}
        self.names: list[str] = []

    def var(self, name: str) -> int:
        if name in self.names:
            raise ValueError("duplicate Boolean input")
        i = len(self.names)
        self.names.append(name)
        self.nodes.append(("var", i, 0))
        return len(self.nodes) - 1

    def _binary(self, op: str, a: int, b: int) -> int:
        if op == "xor":
            if a == b:
                return 0
            if a == 0:
                return b
            if b == 0:
                return a
        else:
            if a == 0 or b == 0:
                return 0
            if a == 1:
                return b
            if b == 1 or a == b:
                return a
        if a > b:
            a, b = b, a
        key = (op, a, b)
        if key not in self.cache:
            self.cache[key] = len(self.nodes)
            self.nodes.append(key)
        return self.cache[key]

    def xor(self, a: int, b: int) -> int:
        return self._binary("xor", a, b)

    def and_(self, a: int, b: int) -> int:
        return self._binary("and", a, b)

    def not_(self, a: int) -> int:
        return self.xor(a, 1)

    def or_(self, a: int, b: int) -> int:
        return self.not_(self.and_(self.not_(a), self.not_(b)))

    def all_(self, bits: list[int]) -> int:
        out = 1
        for bit in bits:
            out = self.and_(out, bit)
        return out

    def any_(self, bits: list[int]) -> int:
        out = 0
        for bit in bits:
            out = self.or_(out, bit)
        return out

    def eq(self, a: int, b: int) -> int:
        return self.not_(self.xor(a, b))

    def implies(self, guard: int, consequence: int) -> int:
        return self.or_(self.not_(guard), consequence)

    def evaluate(self, model: PackedModel, output: int) -> bool:
        model.validate(len(self.names))
        values = bytearray(len(self.nodes))
        values[1] = 1
        for i, (op, a, b) in enumerate(self.nodes[2:], 2):
            if op == "var":
                values[i] = model.bit(a)
            elif op == "xor":
                values[i] = values[a] ^ values[b]
            elif op == "and":
                values[i] = values[a] & values[b]
            else:
                raise AssertionError("unknown DAG operation")
        return bool(values[output])

    def counts(self) -> dict[str, int]:
        return {
            "variables": len(self.names),
            "xor": sum(op == "xor" for op, _, _ in self.nodes),
            "and": sum(op == "and" for op, _, _ in self.nodes),
            "total_nodes": len(self.nodes),
            "model_limbs": (len(self.names) + 63) // 64,
        }


class Field:
    def __init__(self, dag: Dag, n: int, modulus: int) -> None:
        check_field(n, modulus)
        self.dag, self.n, self.modulus = dag, n, modulus

    def const(self, x: int) -> tuple[int, ...]:
        if not 0 <= x < 1 << self.n:
            raise ValueError("noncanonical field constant")
        return tuple((x >> i) & 1 for i in range(self.n))

    def var(self, prefix: str) -> tuple[int, ...]:
        return tuple(self.dag.var(f"{prefix}_{i}") for i in range(self.n))

    def add(self, a: tuple[int, ...], b: tuple[int, ...]) -> tuple[int, ...]:
        return tuple(self.dag.xor(x, y) for x, y in zip(a, b, strict=True))

    def _reduce(self, raw: list[int]) -> tuple[int, ...]:
        if len(raw) != 2 * self.n - 1:
            raise AssertionError("wrong raw product width")
        for k in range(len(raw) - 1, self.n - 1, -1):
            high = raw[k]
            for low in range(self.n):
                if self.modulus >> low & 1:
                    raw[k - self.n + low] = self.dag.xor(raw[k - self.n + low], high)
        return tuple(raw[: self.n])

    def mul(self, a: tuple[int, ...], b: tuple[int, ...]) -> tuple[int, ...]:
        raw = [0] * (2 * self.n - 1)
        for i, x in enumerate(a):
            for j, y in enumerate(b):
                raw[i + j] = self.dag.xor(raw[i + j], self.dag.and_(x, y))
        return self._reduce(raw)

    def square(self, a: tuple[int, ...]) -> tuple[int, ...]:
        raw = [0] * (2 * self.n - 1)
        for i, x in enumerate(a):
            raw[2 * i] = x
        return self._reduce(raw)

    def zero(self, a: tuple[int, ...]) -> int:
        return self.dag.all_([self.dag.not_(x) for x in a])

    def equal(self, a: tuple[int, ...], b: tuple[int, ...]) -> int:
        return self.zero(self.add(a, b))


@dataclass
class Relation:
    dag: Dag
    field: Field
    output: int
    branches: dict[str, int]

    def model(self, p: tuple[int, int, int], q: tuple[int, int, int],
              r: tuple[int, int, int], lam: int) -> PackedModel:
        n = self.field.n
        def encode_point(point: tuple[int, int, int]) -> list[int]:
            o, x, y = point
            if o not in (0, 1) or not 0 <= x < 1 << n or not 0 <= y < 1 << n:
                raise ValueError("noncanonical point-bit encoding")
            return [o] + [(x >> i) & 1 for i in range(n)] + [(y >> i) & 1 for i in range(n)]
        if not 0 <= lam < 1 << n:
            raise ValueError("noncanonical lambda")
        bits = encode_point(p) + encode_point(q) + encode_point(r)
        bits += [(lam >> i) & 1 for i in range(n)]
        return PackedModel.from_bits(bits)

    def accepts(self, p: tuple[int, int, int], q: tuple[int, int, int],
                r: tuple[int, int, int], lam: int) -> bool:
        return self.dag.evaluate(self.model(p, q, r, lam), self.output)


def build_relation(n: int, modulus: int) -> Relation:
    """Build E: y² + xy = x³ + 1 with existential slope λ.

    Finite same-x points are either inverses, or equal with x != 0.
    A derived branch-cover constraint rejects all other assignments.
    """
    d = Dag()
    f = Field(d, n, modulus)
    points = []
    for label in ("p", "q", "r"):
        points.append((d.var(f"{label}_o"), f.var(f"{label}_x"), f.var(f"{label}_y")))
    lam = f.var("lambda")
    (po, px, py), (qo, qx, qy), (ro, rx, ry) = points

    def point_valid(point: tuple[int, tuple[int, ...], tuple[int, ...]]) -> int:
        o, x, y = point
        infinity_is_canonical = d.implies(o, d.and_(f.zero(x), f.zero(y)))
        lhs = f.add(f.square(y), f.mul(x, y))
        rhs = f.add(f.mul(f.square(x), x), f.const(1))
        affine_is_on_curve = d.implies(d.not_(o), f.equal(lhs, rhs))
        return d.and_(infinity_is_canonical, affine_is_on_curve)

    validity = d.all_([point_valid(point) for point in points])
    finite = d.and_(d.not_(po), d.not_(qo))
    same_x = f.equal(px, qx)
    same_y = f.equal(py, qy)
    inverse = f.equal(f.add(py, qy), px)
    nonzero_x = d.not_(f.zero(px))
    branches = {
        "copy_q": po,
        "copy_p": d.and_(d.not_(po), qo),
        "inverse": d.all_([finite, same_x, inverse]),
        "double": d.all_([finite, same_x, d.not_(inverse), same_y, nonzero_x]),
        "generic": d.and_(finite, d.not_(same_x)),
    }
    cover = d.any_(list(branches.values()))

    def equal_points(a: tuple[int, tuple[int, ...], tuple[int, ...]],
                     b: tuple[int, tuple[int, ...], tuple[int, ...]]) -> int:
        return d.all_([d.eq(a[0], b[0]), f.equal(a[1], b[1]), f.equal(a[2], b[2])])

    copy_q = d.implies(branches["copy_q"], equal_points(points[2], points[1]))
    copy_p = d.implies(branches["copy_p"], equal_points(points[2], points[0]))
    inv_result = d.implies(branches["inverse"], d.and_(ro, d.and_(f.zero(rx), f.zero(ry))))

    double_slope = f.equal(f.mul(lam, px), f.add(f.square(px), py))
    double_x = f.add(f.add(f.square(lam), lam), f.const(0))
    double_y = f.add(f.square(px), f.mul(f.add(lam, f.const(1)), rx))
    double_result = d.implies(branches["double"], d.all_([
        d.not_(ro), double_slope, f.equal(rx, double_x), f.equal(ry, double_y),
    ]))

    generic_slope = f.equal(f.mul(lam, f.add(px, qx)), f.add(py, qy))
    generic_x = f.add(f.add(f.square(lam), lam), f.add(px, qx))
    generic_y = f.add(f.add(f.mul(lam, f.add(px, rx)), rx), py)
    generic_result = d.implies(branches["generic"], d.all_([
        d.not_(ro), generic_slope, f.equal(rx, generic_x), f.equal(ry, generic_y),
    ]))
    output = d.all_([validity, cover, copy_q, copy_p, inv_result,
                     double_result, generic_result])
    return Relation(d, f, output, branches)
