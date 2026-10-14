# Low Hamming-weight points for a special factor base.
#
# The factor base is the set of affine points whose x-coordinate has
# Hamming weight at most w.
#
#   ECC2K-130   weight in the permuted type-II optimal normal basis.
#               That weight is invariant under Frobenius, so the base is a
#               union of signed Frobenius orbits. This is the weight bound
#               the fixed Koblitz campaign materializes
#               (docs/ic/params/ecc2k130-fixed.json, weight 2:
#               14 orbit columns, 3,668 signed points).
#               --basis polynomial lists polynomial-basis weight instead;
#               squaring does not preserve that weight.
#
#   P-192/224/256   weight of the integer x in [0, p). Negation does not
#               change x, so each stored point stands for {P, -P}.
#
# Points are on the public curve and in the prime-order subgroup. The
# script lists them. It does not search for a discrete logarithm.
#
#   sage low_hw_factor_base.sage ecc2k130 --weight 2 --out ecc2k130-w2.json
#   sage low_hw_factor_base.sage p192  --weight 2 --out p192-w2.json
#   sage low_hw_factor_base.sage p224  --weight 2 --out p224-w2.json
#   sage low_hw_factor_base.sage p256  --weight 2 --out p256-w2.json
#
# This file is preparsed: bitwise xor is written ^^.

import argparse
import itertools
import json
import sys


# ── shared ────────────────────────────────────────────────────────────

def popcount(n):
    return int(n).bit_count()


def hx(n):
    return "0x" + format(int(n), "x")


def masks_of_weight_at_most(width, weight):
    """Integers in [0, 2^width) whose bits lie in positions 0..width-1."""
    width = int(width)
    weight = int(weight)
    yield 0
    for k in range(1, weight + 1):
        for comb in itertools.combinations(range(width), k):
            mask = 0
            for i in comb:
                mask |= 1 << int(i)
            yield mask


def support(mask):
    mask = int(mask)
    bits = []
    i = 0
    while mask:
        if mask & 1:
            bits.append(i)
        mask >>= 1
        i += 1
    return bits


def emit(document, path):
    text = json.dumps(plain(document), indent=2) + "\n"
    if path:
        with open(path, "w") as handle:
            handle.write(text)
        print("wrote %s" % path, file=sys.stderr)
    else:
        sys.stdout.write(text)


def plain(value):
    """Sage Integer literals are not accepted by json.dumps."""
    if isinstance(value, dict):
        return {str(k): plain(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [plain(v) for v in value]
    if isinstance(value, Integer):
        return int(value)
    return value


# ── ECC2K-130, permuted type-II ONB ───────────────────────────────────

M = 131
# x^131 + x^13 + x^2 + x + 1, the polynomial basis of the fixed profile.
POLY = 0x800000000000000000000000000002007
# Image of the polynomial generator in ONB coordinates.
ONB_ROOT = 0x55bced9c666b0fe76e49ce3445ba9eee7
SUBGROUP = 680564733841876926932320129493409985129
GX = 0x051C99BFA6F18DE467C80C23B98C7994AA
GY = 0x042EA2D112ECEC71FCF7E000D7EFC978BD


class Onb:
    """GF(2^m) as symmetric vectors modulo the all-ones vector, n = 2m+1.

    Coordinate bit i-1 is the coefficient of gamma_i = zeta^i + zeta^{-i}.
    Multiplication is cyclic convolution, matching ecc2k130/codegen/field.py.
    """

    def __init__(self, m):
        self.m = int(m)
        self.n = 2 * self.m + 1
        if not is_prime(self.n):
            raise ValueError("2m+1 = %d is not prime; no type-II ONB" % self.n)
        order = Mod(2, self.n).multiplicative_order()
        if order != self.m and order != 2 * self.m:
            raise ValueError("ord_%d(2) = %s, no type-II ONB" % (self.n, order))
        self.all_ones = (1 << self.n) - 1

    def from_coords(self, a):
        a = int(a)
        u = 0
        for i in range(1, self.m + 1):
            if (a >> (i - 1)) & 1:
                u |= (1 << i) | (1 << (self.n - i))
        return u

    def normalize(self, u):
        u = int(u)
        if u & 1:
            u ^^= self.all_ones
        return u

    def to_coords(self, u):
        u = self.normalize(u)
        a = 0
        for i in range(1, self.m + 1):
            if (u >> i) & 1:
                a |= 1 << (i - 1)
        return a

    def rot(self, u, k):
        u = int(u)
        k = int(k) % self.n
        if k == 0:
            return u
        return ((u << k) | (u >> (self.n - k))) & self.all_ones

    def mul(self, a, b):
        r = 0
        bb = int(b)
        i = 0
        while bb:
            if bb & 1:
                r ^^= self.rot(a, i)
            bb >>= 1
            i += 1
        return self.normalize(r)

    def one(self):
        return self.from_coords((1 << self.m) - 1)


def squaring_permutation(m, n):
    """Where ONB coordinate i lands under x |-> x^2. See decomp.sigmaPerm."""
    e = pow(2, 1, int(n))
    perm = []
    for i in range(int(m)):
        t = ((i + 1) * e) % int(n)
        perm.append((t if t <= m else n - t) - 1)
    return perm


def apply_perm(mask, perm):
    mask = int(mask)
    out = 0
    for i, dest in enumerate(perm):
        if (mask >> i) & 1:
            out |= 1 << int(dest)
    return out


def canonical_mask(mask, perm, m):
    best = int(mask)
    cur = best
    for _ in range(int(m) - 1):
        cur = apply_perm(cur, perm)
        if cur < best:
            best = cur
    return best


def invert_columns(columns, m):
    """Rows of M^{-1}. Column j of M is the bitmask columns[j].

    Bit j of M^{-1} * target is the parity of row j with target.
    """
    m = int(m)
    rows = []
    for i in range(m):
        row = 0
        for j in range(m):
            if (int(columns[j]) >> i) & 1:
                row |= 1 << j
        rows.append([row, 1 << i])
    where = [-1] * m
    at = 0
    for col in range(m):
        sel = None
        for i in range(at, m):
            if (rows[i][0] >> col) & 1:
                sel = i
                break
        if sel is None:
            raise RuntimeError("singular ONB/polynomial change of basis")
        rows[at], rows[sel] = rows[sel], rows[at]
        for i in range(m):
            if i != at and ((rows[i][0] >> col) & 1):
                rows[i][0] ^^= rows[at][0]
                rows[i][1] ^^= rows[at][1]
        where[col] = at
        at += 1
    return [rows[where[col]][1] for col in range(m)]


def apply_rows(rows, target):
    target = int(target)
    out = 0
    for j, row in enumerate(rows):
        if popcount(int(row) & target) & 1:
            out |= 1 << j
    return out


def poly_to_onb(columns, poly_mask):
    acc = 0
    poly_mask = int(poly_mask)
    j = 0
    while poly_mask:
        if poly_mask & 1:
            acc ^^= int(columns[j])
        poly_mask >>= 1
        j += 1
    return acc


def koblitz_order(m):
    """#E(F_{2^m}) for y^2 + xy = x^3 + 1."""
    v0, v1 = 2, -1
    for _ in range(int(m) - 1):
        v0, v1 = v1, -v1 - 2 * v0
    return (1 << int(m)) + 1 - v1


def half_trace(c, m):
    acc = c
    t = c
    for _ in range((int(m) - 1) // 2):
        t = t**4
        acc += t
    return acc


def from_mask(z, mask):
    """Polynomial-basis element whose coefficient of z^i is bit i of mask."""
    mask = int(mask)
    acc = z.parent()(0)
    term = z.parent()(1)
    while mask:
        if mask & 1:
            acc += term
        term *= z
        mask >>= 1
    return acc


def to_mask(elt):
    mask = 0
    for i, coeff in enumerate(elt.polynomial().list()):
        if coeff:
            mask |= 1 << i
    return mask


def setup_ecc2k130():
    onb = Onb(M)
    root = onb.from_coords(ONB_ROOT)
    power = onb.one()
    columns = []
    acc = 0
    for i in range(M + 1):
        if (POLY >> i) & 1:
            acc ^^= power
        if i < M:
            columns.append(onb.to_coords(power))
        power = onb.mul(power, root)
    if onb.normalize(acc) != 0:
        raise RuntimeError("onb_root is not a root of the field polynomial")
    inv_rows = invert_columns(columns, M)
    for j in range(M):
        if apply_rows(inv_rows, columns[j]) != (1 << j):
            raise RuntimeError("change of basis failed to invert")

    R = PolynomialRing(GF(2), "X")
    X = R.gen()
    modulus = X**M + X**13 + X**2 + X + 1
    K = GF(2**M, name="z", modulus=modulus)
    E = EllipticCurve(K, [1, 0, 0, 0, 1])
    perm = squaring_permutation(M, onb.n)

    z = K.gen()
    sample = (1 << 0) | (1 << 7) | (1 << 20)
    sample_poly = apply_rows(inv_rows, sample)
    squared = from_mask(z, sample_poly)**2
    got = poly_to_onb(columns, to_mask(squared))
    if got != apply_perm(sample, perm):
        raise RuntimeError("ONB squaring permutation does not match the field")

    G = E(from_mask(z, GX), from_mask(z, GY))
    order = koblitz_order(M)
    if order != 4 * SUBGROUP:
        raise RuntimeError("curve order %d != 4 * subgroup" % order)
    if not is_prime(SUBGROUP):
        raise RuntimeError("subgroup order is not prime")
    if SUBGROUP * G != E(0) or G == E(0):
        raise RuntimeError("generator is not a point of order the subgroup")

    root_mod = Mod(-7, SUBGROUP).sqrt()
    inv2 = inverse_mod(2, SUBGROUP)
    eigen = None
    frob_G = E(G[0]**2, G[1]**2)
    for sign in (1, -1):
        cand = ((-1 + sign * Integer(root_mod)) * inv2) % SUBGROUP
        if Integer(cand) * G == frob_G:
            eigen = int(cand)
            break
    if eigen is None:
        raise RuntimeError("no Frobenius eigenvalue on the prime subgroup")
    return K, E, columns, inv_rows, perm, G, eigen


def lift_binary(E, K, x):
    """One y with y^2 + x y = x^3 + 1, or None. The other y is x + y."""
    if x == 0:
        return None
    c = x + x**(-2)
    if c.trace() != 0:
        return None
    y = x * half_trace(c, M)
    return E(x, y)


def poly_coords(point):
    return to_mask(point[0]), to_mask(point[1])


def signed_orbit_size(E, P, m):
    seen = set()
    Q = P
    for _ in range(int(m)):
        for S in (Q, -Q):
            seen.add(poly_coords(S))
        Q = E(Q[0]**2, Q[1]**2)
        if Q == P:
            break
    return len(seen)


def enumerate_ecc2k130(weight, basis, limit):
    if weight < 0 or weight > 4:
        raise SystemExit("ecc2k130 weight must be between 0 and 4")
    K, E, columns, inv_rows, perm, _G, eigen = setup_ecc2k130()
    reps = []
    seen = set()
    tested = 0
    odd_skipped = 0
    for mask in masks_of_weight_at_most(M, weight):
        if mask == 0:
            continue
        if basis == "onb":
            # Trace equals ONB weight mod 2. Prime-subgroup points have odd
            # order, so they are doubles, and x of a double has trace 0.
            if popcount(mask) & 1:
                odd_skipped += 1
                continue
            canon = canonical_mask(mask, perm, M)
            if canon in seen:
                continue
            seen.add(canon)
            if popcount(canon) > weight:
                continue
            x = from_mask(K.gen(), apply_rows(inv_rows, canon))
            onb_x = canon
        else:
            x = from_mask(K.gen(), mask)
            key = to_mask(x)
            if key in seen:
                continue
            seen.add(key)
            onb_x = poly_to_onb(columns, key)
        P = lift_binary(E, K, x)
        if P is None:
            continue
        tested += 1
        if SUBGROUP * P != E(0):
            continue
        x_poly, y_poly = poly_coords(P)
        neg_y = poly_coords(-P)[1]
        image = to_mask(P[0]**2)
        rec = {
            "x": hx(x_poly),
            "y": hx(y_poly),
            "neg_y": hx(neg_y),
            "weight": popcount(onb_x if basis == "onb" else x_poly),
            "onb_x": hx(onb_x),
            "onb_weight": popcount(onb_x),
            "polynomial_weight": popcount(x_poly),
            "support": support(onb_x if basis == "onb" else x_poly),
            "frobenius_image_polynomial_weight": popcount(image),
        }
        if basis == "onb":
            rec["signed_orbit_size"] = signed_orbit_size(E, P, M)
        reps.append(rec)
        if limit and len(reps) >= limit:
            break
    reps.sort(key=lambda rec: int(rec["x"], 16))
    signed = sum(rec.get("signed_orbit_size", 2) for rec in reps)
    document = {
        "schema": "low-hw-factor-base-v1",
        "curve": "ECC2K-130",
        "equation": "y^2 + x*y = x^3 + 1",
        "field": {
            "kind": "binary",
            "degree": M,
            "polynomial": hx(POLY),
            "basis": "type-ii-onb" if basis == "onb" else "polynomial",
            "onb_root": hx(ONB_ROOT),
        },
        "subgroup_order": str(SUBGROUP),
        "cofactor": 4,
        "frobenius_eigenvalue": str(eigen),
        "weight": weight,
        "identify": "signed-frobenius-orbit" if basis == "onb" else "negation-only",
        "frobenius_stable": basis == "onb",
        "odd_weight_skipped": odd_skipped,
        "abscissae_tested": tested,
        "orbit_columns": len(reps),
        "signed_points": signed,
        "note": (
            "Public curve points in the prime subgroup. "
            "One half-trace y per abscissa; neg_y is the other sign. "
            "ONB orbits are the relation columns of the weight-bounded base."
        ),
        "representatives": [
            {"x": rec["x"], "y": rec["y"]} for rec in reps
        ],
        "points": reps,
    }
    if basis == "onb" and weight == 2 and not limit:
        if len(reps) != 14 or signed != 3668:
            raise SystemExit(
                "weight-2 ONB base has %d orbits and %d signed points; "
                "expected 14 and 3668" % (len(reps), signed)
            )
        document["matches_fixed_campaign"] = True
    return document


# ── NIST prime curves ─────────────────────────────────────────────────

# Parameters copied from src/ecc/curve.rs (P-256) and src/ecc/curve_zoo.rs
# (P-192, P-224). Each has a = -3 and cofactor 1.
NIST = {
    "p192": {
        "name": "P-192",
        "p": 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFFFFFFFFFFFF,
        "a": 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFFFFFFFFFFFC,
        "b": 0x64210519E59C80E70FA7E9AB72243049FEB8DEECC146B9B1,
        "gx": 0x188DA80EB03090F67CBF20EB43A18800F4FF0AFD82FF1012,
        "gy": 0x07192B95FFC8DA78631011ED6B24CDD573F977A11E794811,
        "n": 0xFFFFFFFFFFFFFFFFFFFFFFFF99DEF836146BC9B1B4D22831,
    },
    "p224": {
        "name": "P-224",
        "p": 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF000000000000000000000001,
        "a": 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFFFFFFFFFFFFFFFFFFFE,
        "b": 0xB4050A850C04B3ABF54132565044B0B7D7BFD8BA270B39432355FFB4,
        "gx": 0xB70E0CBD6BB4BF7F321390B94A03C1D356C21122343280D6115C1D21,
        "gy": 0xBD376388B5F723FB4C22DFE6CD4375A05A07476444D5819985007E34,
        "n": 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFF16A2E0B8F03E13DD29455C5C2A3D,
    },
    "p256": {
        "name": "P-256",
        "p": 0xFFFFFFFF00000001000000000000000000000000FFFFFFFFFFFFFFFFFFFFFFFF,
        "a": 0xFFFFFFFF00000001000000000000000000000000FFFFFFFFFFFFFFFFFFFFFFFC,
        "b": 0x5AC635D8AA3A93E7B3EBBD55769886BC651D06B0CC53B0F63BCE3C3E27D2604B,
        "gx": 0x6B17D1F2E12C4247F8BCE6E563A440F277037D812DEB33A0F4A13945D898C296,
        "gy": 0x4FE342E2FE1A7F9B8EE7EB4A7C0F9E162BCE33576B315ECECBB6406837BF51F5,
        "n": 0xFFFFFFFF00000000FFFFFFFFFFFFFFFFBCE6FAADA7179E84F3B9CAC2FC632551,
    },
}


def enumerate_prime(key, weight, limit):
    if weight < 0 or weight > 4:
        raise SystemExit("%s weight must be between 0 and 4" % key)
    spec = NIST[key]
    p = int(spec["p"])
    a = int(spec["a"])
    b = int(spec["b"])
    n = int(spec["n"])
    F = GF(p)
    E = EllipticCurve(F, [a, b])
    G = E(spec["gx"], spec["gy"])
    if not is_prime(n) or n * G != E(0) or G == E(0):
        raise RuntimeError("%s generator does not have prime order n" % spec["name"])
    cardinality = int(E.cardinality())
    if cardinality != n:
        raise RuntimeError("%s has %d points, cofactor is not 1" % (spec["name"], cardinality))

    width = p.bit_length()
    points = []
    by_weight = {}
    scanned = 0
    for mask in masks_of_weight_at_most(width, weight):
        scanned += 1
        if scanned % 200000 == 0:
            print("%s scanned %d" % (spec["name"], scanned), file=sys.stderr)
        x_int = int(mask)
        if x_int >= p:
            continue
        rhs = (x_int * x_int * x_int + a * x_int + b) % p
        if rhs == 0:
            y_int = 0
        elif pow(rhs, (p - 1) // 2, p) != 1:
            continue
        else:
            y_int = int(F(rhs).sqrt())
            if y_int > p - y_int:
                y_int = p - y_int
        P = E(F(x_int), F(y_int))
        wt = popcount(x_int)
        by_weight[wt] = by_weight.get(wt, 0) + 1
        points.append({
            "x": hx(x_int),
            "y": hx(y_int),
            "neg_y": hx(int((-P)[1])),
            "weight": wt,
            "bits": support(x_int),
        })
        if limit and len(points) >= limit:
            break
    points.sort(key=lambda rec: (rec["weight"], int(rec["x"], 16)))
    # Cofactor 1: every affine point has order n. Check three of them.
    for rec in points[:3]:
        Q = E(int(rec["x"], 16), int(rec["y"], 16))
        if n * Q != E(0):
            raise RuntimeError("point failed the order check: %s" % rec["x"])
    return {
        "schema": "low-hw-factor-base-v1",
        "curve": spec["name"],
        "equation": "y^2 = x^3 - 3*x + b",
        "field": {"kind": "prime", "modulus": hx(p), "basis": "integer-bits"},
        "a": hx(a),
        "b": hx(b),
        "subgroup_order": hx(n),
        "cofactor": 1,
        "cardinality": cardinality,
        "weight": weight,
        "identify": "negation",
        "counts_by_weight": {str(k): by_weight[k] for k in sorted(by_weight)},
        "point_count": len(points),
        "note": (
            "Public curve points. y is the square root in [0, (p-1)/2]; "
            "neg_y is the other sign. Each row is one {P, -P} class."
        ),
        "points": points,
    }


def main():
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("--weight", type=int, default=2, help="maximum Hamming weight (default 2)")
    common.add_argument("--out", help="write the JSON document here; default stdout")
    common.add_argument("--limit", type=int, default=0, help="stop after this many points")
    parser = argparse.ArgumentParser(
        description="Enumerate low Hamming-weight points for a special factor base."
    )
    sub = parser.add_subparsers(dest="curve", required=True)
    binary = sub.add_parser(
        "ecc2k130", parents=[common],
        help="Koblitz curve y^2 + xy = x^3 + 1 over F_2^131",
    )
    binary.add_argument("--basis", choices=("onb", "polynomial"), default="onb")
    sub.add_parser("p192", parents=[common], help="NIST P-192 / secp192r1")
    sub.add_parser("p224", parents=[common], help="NIST P-224 / secp224r1")
    sub.add_parser("p256", parents=[common], help="NIST P-256 / secp256r1")
    # Sage preparses the script but leaves sys.argv for the script.
    args = parser.parse_args()
    if args.curve == "ecc2k130":
        document = enumerate_ecc2k130(args.weight, args.basis, args.limit)
        kind = "orbits" if args.basis == "onb" else "points"
        print(
            "ECC2K-130 basis=%s weight<=%d: %d %s, %d signed points"
            % (args.basis, args.weight, document["orbit_columns"], kind, document["signed_points"]),
            file=sys.stderr,
        )
    else:
        document = enumerate_prime(args.curve, args.weight, args.limit)
        print(
            "%s weight<=%d: %d negation classes %s"
            % (document["curve"], args.weight, document["point_count"], document["counts_by_weight"]),
            file=sys.stderr,
        )
    emit(document, args.out)


if __name__ == "__main__":
    main()
