"""Fixed GF(256) ordinary curves and explicit degree-three Velu maps.

Only the eight-bit toy field is supported. No external curve/target inputs.
"""
from itertools import product

O = (-1, -1)
MODULUS = 0x11B


def multiply(a, b):
    result = 0
    while b:
        if b & 1:
            result ^= a
        b >>= 1
        a <<= 1
        if a & 256:
            a ^= MODULUS
    return result


MUL = tuple(tuple(multiply(a, b) for b in range(256)) for a in range(256))
SQUARE = tuple(MUL[a][a] for a in range(256))
INV = (0,) + tuple(next(b for b in range(1, 256) if MUL[a][b] == 1) for a in range(1, 256))


class Curve:
    """y^2+xy=x^3+b, in a fixed polynomial basis."""

    def __init__(self, b=1):
        if not 1 <= b < 256:
            raise ValueError("nonzero GF(256) coefficient required")
        self.b = b

    def contains(self, p):
        if p == O:
            return True
        x, y = p
        return (0 <= x < 256 and 0 <= y < 256
                and SQUARE[y] ^ MUL[x][y] == MUL[SQUARE[x]][x] ^ self.b)

    def points(self):
        return (O,) + tuple((x, y) for x in range(256) for y in range(256)
                            if SQUARE[y] ^ MUL[x][y] == MUL[SQUARE[x]][x] ^ self.b)

    @staticmethod
    def neg(p):
        return O if p == O else (p[0], p[0] ^ p[1])

    @staticmethod
    def add(p, q):
        if p == O:
            return q
        if q == O:
            return p
        x, y = p
        u, v = q
        if x == u:
            if y != v or x == 0:
                return O
            slope = x ^ MUL[y][INV[x]]
            out_x = SQUARE[slope] ^ slope
            return out_x, SQUARE[x] ^ MUL[slope ^ 1][out_x]
        slope = MUL[y ^ v][INV[x ^ u]]
        out_x = SQUARE[slope] ^ slope ^ x ^ u
        return out_x, MUL[slope][x ^ out_x] ^ out_x ^ y

    def scalar(self, k, p):
        result = O
        while k:
            if k & 1:
                result = self.add(result, p)
            p = self.add(p, p)
            k >>= 1
        return result

    def total(self, points):
        result = O
        for p in points:
            result = self.add(result, p)
        return result


class Velu3:
    """Normalized quotient by {O,Q,-Q}, with an explicit degree-three x map."""

    def __init__(self, source, q):
        if q == O or not source.contains(q) or source.scalar(3, q) != O:
            raise ValueError("kernel generator must have exact order three")
        self.source, self.q = source, q
        self.kernel = {O, q, source.neg(q)}
        self.target = Curve(source.b ^ q[0] ^ SQUARE[q[0]])

    def __call__(self, p):
        if p in self.kernel:
            return O
        plus = self.source.add(p, self.q)
        minus = self.source.add(p, self.source.neg(self.q))
        # Velu y followed by y -> y+x(Q), eliminating the codomain a4 term.
        return p[0] ^ plus[0] ^ minus[0], p[1] ^ plus[1] ^ minus[1]

    def rational_x(self, x):
        denominator = x ^ self.q[0]
        if denominator == 0:
            raise ZeroDivisionError("isogeny pole")
        d = MUL[self.q[0]][INV[denominator]]
        return x ^ d ^ SQUARE[d]


def certified_neighbors():
    source = Curve()
    points = source.points()
    h = tuple(p for p in points if source.scalar(32, p) == O)
    three = tuple(p for p in points if source.scalar(3, p) == O)
    if len(points) != 288 or len(h) != 32 or len(three) != 9:
        raise AssertionError("source order/torsion certificate failed")
    kernels = sorted({min(p, source.neg(p)) for p in three if p != O})
    if len(kernels) != 4:
        raise AssertionError("expected four rational cyclic order-three kernels")
    result = []
    for q in kernels:
        phi = Velu3(source, q)
        target = phi.target
        codomain_points = target.points()
        mapped = {p: phi(p) for p in points}
        if len(codomain_points) != 288 or any(not target.contains(v) for v in mapped.values()):
            raise AssertionError("codomain equation or point count failed")
        if {p for p in points if mapped[p] == O} != phi.kernel:
            raise AssertionError("wrong isogeny kernel")
        if any(phi.rational_x(p[0]) != mapped[p][0] for p in points if p not in phi.kernel):
            raise AssertionError("rational x map mismatch")
        # Numerator x^3+(q_x^2+q_x)x has value q_x^2 at its denominator's root.
        if q[0] == 0 or SQUARE[q[0]] == 0:
            raise AssertionError("x-map degree drops through cancellation")
        if len({mapped[p] for p in h}) != 32:
            raise AssertionError("transport is not injective on the chosen subgroup")
        for p, r in product(h, repeat=2):
            if mapped[source.add(p, r)] != target.add(mapped[p], mapped[r]):
                raise AssertionError("homomorphism certificate failed")
        target_three = tuple(p for p in codomain_points if target.scalar(3, p) == O)
        if len(target_three) != 3:
            raise AssertionError("neighbor unexpectedly retains full rational three-torsion")
        dual = Velu3(target, target_three[1])
        if dual.target.b != source.b:
            raise AssertionError("dual codomain is not the original normalized model")
        sign = next((s for s in (1, -1) if all(
            (dual(mapped[p]) if s == 1 else source.neg(dual(mapped[p]))) == source.scalar(3, p)
            for p in points)), None)
        if sign is None:
            raise AssertionError("dual composition is not [3]")
        for p in codomain_points:
            back = dual(p) if sign == 1 else source.neg(dual(p))
            if phi(back) != target.scalar(3, p):
                raise AssertionError("reverse dual composition is not [3]")
        certificate = {"kernel_generator": q, "kernel": sorted(phi.kernel),
            "source_b": source.b, "target_b": target.b, "target_j": INV[target.b],
            "source_order": len(points), "target_order": len(codomain_points),
            "subgroup_order": len(h), "degree": 3, "degree_cancellation_check": SQUARE[q[0]],
            "source_rational_3_torsion": 9, "target_rational_3_torsion": 3,
            "dual_kernel_generator": target_three[1], "dual_sign": sign,
            "dual_forward_checked_points": len(points), "dual_reverse_checked_points": len(codomain_points),
            "homomorphism_checked_pairs": len(h)**2,
            "trace": -31, "frobenius_discriminant": -63, "fundamental_discriminant": -7,
            "frobenius_conductor": 3, "source_endomorphism_conductor": 1,
            "target_endomorphism_conductor": 3,
            "direction": "descending", "status": "VERIFIED_TOY_MAP_AND_TORSION",
            "conductor_reason": "source descends from the maximal order over F2; 3 is inert in Q(sqrt(-7)); -63=-7*3^2; the nonmaximal neighbor is confirmed by rational 3-torsion dropping from 9 to 3"}
        result.append((phi, certificate))
    if len({phi.target.b for phi, _ in result}) != 4:
        raise AssertionError("codomains are not four distinct j-invariants")
    return source, h, result


def semaev3(x, y, z, b):
    pair_sum = MUL[x][y] ^ MUL[x][z] ^ MUL[y][z]
    return SQUARE[pair_sum] ^ MUL[MUL[x][y]][z] ^ b


def semaev4(x, y, z, r, b):
    # Resultant of A*u^2+B*u+C and D*u^2+E*u+F in characteristic two.
    a, bb, c = SQUARE[x ^ y], MUL[x][y], SQUARE[MUL[x][y]] ^ b
    d, e, f = SQUARE[z ^ r], MUL[z][r], SQUARE[MUL[z][r]] ^ b
    return SQUARE[MUL[a][f] ^ MUL[c][d]] ^ MUL[MUL[a][e] ^ MUL[bb][d]][MUL[bb][f] ^ MUL[c][e]]


def decomposition_value(xs, target, b):
    if len(xs) == 2:
        return xs[0] ^ xs[1] if target == O else semaev3(*xs, target[0], b)
    if len(xs) == 3:
        return semaev3(*xs, b) if target == O else semaev4(*xs, target[0], b)
    raise ValueError("only two and three toy summands are supported")
