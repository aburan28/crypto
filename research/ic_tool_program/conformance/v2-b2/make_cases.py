#!/usr/bin/env python3
"""Conformance cases for B2: prime and extension fields, the other curve
forms, the importers, estimates and budgets
(`research/ic_tool_program/design/schema-v2.md` §9, C032–C051;
`research/ic_tool_program/rounds/B2-fields-forms-estimates/PROTOCOL.md`).

    python3 make_cases.py           # writes params/, cases.json and SHA256SUMS; refuses to overwrite
    python3 make_cases.py --check   # re-derives every file and compares the bytes

Written at B2's declaration, before any B2 code, as the design requires.
B1's case files (`../v2/`) are frozen, so B2's sit here; `../run.py` runs
every step's cases.

The arithmetic is its own, as B1's is:
- binary fields and binary curves are `../v2/make_cases.py`'s;
- prime fields, `GF(p^3)`, the general binary Weierstrass form and the
  group orders are written here.

None of it is shared with the Rust tool.  Each file's intended property
is asserted when it is built: the field's, the curve's, the order of
the group by an exact method, the subgroup's, the target's, and the
defect a refusal case carries.

Every choice comes from a label, as in B1's generator: SHAKE-256 of the
label gives a field element, and `:1`, `:2`, … are appended until it
works.

**Group orders.** `#E` is found by the orders of points, as follows. For
each of several points `P`, every `m` in the Hasse interval with
`[m]P = O` is found by baby steps and giant steps. The sets are then
intersected until one value is left. That is exact. It is checked
against enumeration on small fields by `selftest`.
"""
from __future__ import annotations

import copy
import hashlib
import importlib.util
import json
import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
V2 = HERE.parent / "v2"
_spec = importlib.util.spec_from_file_location("conformance_v2_cases", V2 / "make_cases.py")
v2 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(v2)

STEP = "B2"
RHO_SEED = 0x230000 + 1
PRIME_TEST = v2.prime_exact  # 13-base Miller-Rabin, exact below psi_13


def shake(label: str, nbytes: int) -> int:
    return int.from_bytes(hashlib.shake_256(label.encode()).digest(nbytes), "big")


def labelled(label: str, i: int) -> str:
    return f"{label}:{i}" if i else label


def next_prime(n: int) -> int:
    while not PRIME_TEST(n):
        n += 1
    return n


def trial_factor(n: int, bound: int) -> tuple[list[int], int]:
    """The prime factors of n up to `bound`, with multiplicity, and what is left."""
    out, d = [], 2
    while d <= bound and d * d <= n:
        while n % d == 0:
            out.append(d)
            n //= d
        d += 1 if d == 2 else 2
    if 1 < n <= bound:
        out.append(n)
        n = 1
    return out, n


# -- the order of a group, from the orders of its points -----------------------------------------

def multiples_in(group, P, lo: int, hi: int) -> set[int]:
    """Every m in [lo, hi] with [m]P = O, by baby steps and giant steps."""
    width = hi - lo
    b = math.isqrt(width) + 1
    table, R = {}, None
    for j in range(b):
        if R is None and j > 0:
            # P's order j is at most b: its multiples in the interval, directly.
            return {m for m in range(lo, hi + 1) if m % j == 0}
        table.setdefault(R, j)
        R = group.add(R, P)
    step = group.neg(group.mul(b, P))
    Q = group.neg(group.mul(lo, P))
    hits = set()
    for g in range(width // b + 2):
        j = table.get(Q)
        if j is not None and g * b + j <= width:
            hits.add(lo + g * b + j)
        Q = group.add(Q, step)
    return hits


def group_order(group, q: int, points) -> int:
    """#E for a curve over a field of q elements: the Hasse interval narrowed by points' orders."""
    s = math.isqrt(q)
    lo, hi = q + 1 - 2 * s - 2, q + 1 + 2 * s + 2
    cands = None
    for P in points:
        hits = multiples_in(group, P, lo, hi)
        cands = hits if cands is None else cands & hits
        if len(cands) == 1:
            return next(iter(cands))
    raise ValueError("the Hasse interval did not narrow to one order")


# -- GF(p) and y^2 = x^3 + ax + b -----------------------------------------------------------------

def sqrt_mod(a: int, p: int):
    """A square root of a modulo an odd prime p, or None (Tonelli-Shanks)."""
    a %= p
    if a == 0:
        return 0
    if pow(a, (p - 1) // 2, p) != 1:
        return None
    q, s = p - 1, 0
    while q % 2 == 0:
        q //= 2
        s += 1
    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1:
        z += 1
    m, c, t, r = s, pow(z, q, p), pow(a, q, p), pow(a, (q + 1) // 2, p)
    while t != 1:
        i, t2 = 0, t
        while t2 != 1:
            t2 = t2 * t2 % p
            i += 1
        b = pow(c, 1 << (m - i - 1), p)
        m, c, t, r = i, b * b % p, t * b * b % p, r * b % p
    return r


class PrimeCurve:
    """y^2 = x^3 + ax + b over GF(p), p > 3; a point is (x, y) or None."""

    def __init__(self, p: int, a: int, b: int):
        self.p, self.a, self.b = p, a % p, b % p
        assert (4 * self.a ** 3 + 27 * self.b ** 2) % p, "non-singular"

    def on_curve(self, P) -> bool:
        if P is None:
            return True
        x, y = P
        return (y * y - x * x * x - self.a * x - self.b) % self.p == 0

    def neg(self, P):
        return None if P is None else (P[0], (-P[1]) % self.p)

    def add(self, P, Q):
        if P is None:
            return Q
        if Q is None:
            return P
        p = self.p
        (x1, y1), (x2, y2) = P, Q
        if x1 == x2:
            if (y1 + y2) % p == 0:
                return None
            lam = (3 * x1 * x1 + self.a) * pow(2 * y1, -1, p) % p
        else:
            lam = (y2 - y1) * pow(x2 - x1, -1, p) % p
        x3 = (lam * lam - x1 - x2) % p
        return x3, (lam * (x1 - x3) - y1) % p

    def mul(self, k: int, P):
        R = None
        for bit in bin(k)[2:] if k > 0 else "":
            R = self.add(R, R)
            if bit == "1":
                R = self.add(R, P)
        return R

    def point(self, label: str):
        nbytes = (self.p.bit_length() + 7) // 8 + 8
        i = 0
        while True:
            x = shake(labelled(label, i), nbytes) % self.p
            y = sqrt_mod(x * x * x + self.a * x + self.b, self.p)
            if y is not None:
                return x, min(y, self.p - y)
            i += 1

    def order(self, label: str) -> int | None:
        try:
            return group_order(self, self.p, (self.point(f"{label}/order/{i}") for i in range(64)))
        except ValueError:
            return None

    def subgroup_point(self, label: str, h: int, r: int):
        i = 0
        while True:
            G = self.mul(h, self.point(labelled(label, i)))
            if G is not None:
                assert self.mul(r, G) is None
                return G
            i += 1


# -- GF(p^3) and y^2 = x^3 + ax + b over it -------------------------------------------------------

class Fp3:
    """GF(p)[t] / (t^3 + c2 t^2 + c1 t + c0); elements are 3-tuples (e0, e1, e2)."""

    def __init__(self, p: int, c: tuple[int, int, int]):
        self.p, self.c = p, tuple(x % p for x in c)
        self.q = p ** 3
        self.zero, self.one = (0, 0, 0), (1, 0, 0)

    def irreducible(self) -> bool:
        """A cubic is irreducible over GF(p) exactly when it has no root there."""
        c0, c1, c2 = self.c
        return all((x * x * x + c2 * x * x + c1 * x + c0) % self.p for x in range(self.p))

    def add(self, u, v):
        return tuple((a + b) % self.p for a, b in zip(u, v))

    def sub(self, u, v):
        return tuple((a - b) % self.p for a, b in zip(u, v))

    def mul(self, u, v):
        p, (c0, c1, c2) = self.p, self.c
        w = [0] * 5
        for i in range(3):
            for j in range(3):
                w[i + j] += u[i] * v[j]
        for k in (4, 3):  # t^3 = -(c2 t^2 + c1 t + c0)
            hi, w[k] = w[k], 0
            w[k - 1] -= hi * c2
            w[k - 2] -= hi * c1
            w[k - 3] -= hi * c0
        return tuple(x % p for x in w[:3])

    def scalar(self, k: int):
        return (k % self.p, 0, 0)

    def pow(self, u, e: int):
        r = self.one
        for bit in bin(e)[2:]:
            r = self.mul(r, r)
            if bit == "1":
                r = self.mul(r, u)
        return r

    def inv(self, u):
        assert u != self.zero
        return self.pow(u, self.q - 2)

    def is_square(self, u) -> bool:
        return u == self.zero or self.pow(u, (self.q - 1) // 2) == self.one

    def sqrt(self, u):
        if u == self.zero:
            return self.zero
        if not self.is_square(u):
            return None
        q, s = self.q - 1, 0
        while q % 2 == 0:
            q //= 2
            s += 1
        z = (2, 1, 0)
        while self.is_square(z):
            z = self.add(z, self.one)
        m, c, t, r = s, self.pow(z, q), self.pow(u, q), self.pow(u, (q + 1) // 2)
        while t != self.one:
            i, t2 = 0, t
            while t2 != self.one:
                t2 = self.mul(t2, t2)
                i += 1
            b = self.pow(c, 1 << (m - i - 1))
            m, c, t, r = i, self.mul(b, b), self.mul(t, self.mul(b, b)), self.mul(r, b)
        return r

    def element(self, label: str):
        return tuple(shake(f"{label}/{i}", 16) % self.p for i in range(3))


class Fp3Curve:
    def __init__(self, k: Fp3, a, b):
        self.k, self.a, self.b = k, a, b
        f = k
        disc = f.add(f.mul(f.scalar(4), f.mul(a, f.mul(a, a))), f.mul(f.scalar(27), f.mul(b, b)))
        assert disc != f.zero, "non-singular"

    def rhs(self, x):
        f = self.k
        return f.add(f.add(f.mul(x, f.mul(x, x)), f.mul(self.a, x)), self.b)

    def on_curve(self, P) -> bool:
        return P is None or self.k.mul(P[1], P[1]) == self.rhs(P[0])

    def neg(self, P):
        return None if P is None else (P[0], self.k.sub(self.k.zero, P[1]))

    def add(self, P, Q):
        f = self.k
        if P is None:
            return Q
        if Q is None:
            return P
        (x1, y1), (x2, y2) = P, Q
        if x1 == x2:
            if f.add(y1, y2) == f.zero:
                return None
            num = f.add(f.mul(f.scalar(3), f.mul(x1, x1)), self.a)
            lam = f.mul(num, f.inv(f.mul(f.scalar(2), y1)))
        else:
            lam = f.mul(f.sub(y2, y1), f.inv(f.sub(x2, x1)))
        x3 = f.sub(f.sub(f.mul(lam, lam), x1), x2)
        return x3, f.sub(f.mul(lam, f.sub(x1, x3)), y1)

    def mul(self, k: int, P):
        R = None
        for bit in bin(k)[2:] if k > 0 else "":
            R = self.add(R, R)
            if bit == "1":
                R = self.add(R, P)
        return R

    def point(self, label: str):
        i = 0
        while True:
            x = self.k.element(labelled(label, i))
            y = self.k.sqrt(self.rhs(x))
            if y is not None:
                return x, y
            i += 1


# -- y^2 + a1 xy + a3 y = x^3 + a2 x^2 + a4 x + a6 over GF(2^n) ------------------------------------

class GeneralBinary:
    """The general Weierstrass form over GF(2^n), Silverman III.2.3's formulas in characteristic 2."""

    def __init__(self, field, a1: int, a2: int, a3: int, a4: int, a6: int):
        self.k = field
        self.a1, self.a2, self.a3, self.a4, self.a6 = a1, a2, a3, a4, a6

    def on_curve(self, P) -> bool:
        if P is None:
            return True
        k, (x, y) = self.k, P
        lhs = k.sqr(y) ^ k.mul(self.a1, k.mul(x, y)) ^ k.mul(self.a3, y)
        rhs = k.mul(k.sqr(x), x) ^ k.mul(self.a2, k.sqr(x)) ^ k.mul(self.a4, x) ^ self.a6
        return lhs == rhs

    def neg(self, P):
        return None if P is None else (P[0], P[1] ^ self.k.mul(self.a1, P[0]) ^ self.a3)

    def add(self, P, Q):
        k = self.k
        if P is None:
            return Q
        if Q is None:
            return P
        (x1, y1), (x2, y2) = P, Q
        if x1 == x2:
            if Q == self.neg(P):
                return None
            den = k.mul(self.a1, x1) ^ self.a3
            lam = k.mul(k.sqr(x1) ^ self.a4 ^ k.mul(self.a1, y1), k.inv(den))
            nu = k.mul(k.mul(k.sqr(x1), x1) ^ k.mul(self.a4, x1) ^ k.mul(self.a3, y1), k.inv(den))
        else:
            inv = k.inv(x1 ^ x2)
            lam = k.mul(y1 ^ y2, inv)
            nu = k.mul(k.mul(y1, x2) ^ k.mul(y2, x1), inv)
        x3 = k.sqr(lam) ^ k.mul(self.a1, lam) ^ self.a2 ^ x1 ^ x2
        y3 = k.mul(lam ^ self.a1, x3) ^ nu ^ self.a3
        return x3, y3

    def mul(self, m: int, P):
        R = None
        for bit in bin(m)[2:] if m > 0 else "":
            R = self.add(R, R)
            if bit == "1":
                R = self.add(R, P)
        return R


class BinaryGroup:
    """B1's y^2 + xy = x^3 + ax^2 + b with the negation the order search needs."""

    def __init__(self, curve):
        self.c = curve

    def add(self, P, Q):
        return self.c.add(P, Q)

    def neg(self, P):
        return None if P is None else (P[0], P[0] ^ P[1])

    def mul(self, m: int, P):
        return self.c.mul(m, P)


def binary_order(curve, n: int, label: str) -> int | None:
    """#E of a binary curve (B1's Curve class) by the orders of its points, or None if they do not
    narrow the Hasse interval to one value."""
    try:
        return group_order(BinaryGroup(curve), 1 << n, (curve.point(f"{label}/order/{i}") for i in range(64)))
    except ValueError:
        return None


def brent_factor(n: int, label: str, steps: int):
    """A nontrivial factor of composite n by Pollard-Brent within about `steps` steps, or None."""
    if n % 2 == 0:
        return 2
    for attempt in range(4):
        c = shake(f"{label}/brent/{attempt}", 16) % (n - 1) + 1
        y, g, r_, q, x, ys = 2, 1, 1, 1, 2, 2
        while g == 1 and r_ <= steps:
            x = y
            for _ in range(r_):
                y = (y * y + c) % n
            k = 0
            while k < r_ and g == 1:
                ys = y
                for _ in range(min(64, r_ - k)):
                    y = (y * y + c) % n
                    q = q * abs(x - y) % n
                g = math.gcd(q, n)
                k += 64
            r_ *= 2
        if g == n:
            g = 1
            while g == 1:
                ys = (ys * ys + c) % n
                g = math.gcd(abs(x - ys), n)
        if 1 < g < n:
            return g
    return None


# -- documents ------------------------------------------------------------------------------------

def dec_point(P) -> dict:
    return {"x": str(P[0]), "y": str(P[1])}


def arr(u) -> list[str]:
    return [str(c) for c in u]


def prime_doc(name: str, p: int, curve: dict, r: int, h: int, generator: dict, target: dict,
              method: dict) -> dict:
    return {"schema_version": 2, "name": name, "field": {"kind": "prime", "p": str(p)}, "curve": curve,
            "subgroup": {"order": str(r), "cofactor": str(h), "generator": generator},
            "target": target, "method": method}


def known_log(label: str, r: int) -> int:
    return shake(label, 32) % (r - 1) + 1


def rho_alone(seed: int = 7) -> dict:
    return {"solve": "rho", "rho": {"pipeline": "auto", "seed": seed}}


def paired_auto(seed: int = 7) -> dict:
    return {"solve": "paired", "index_calculus": {"pipeline": "auto", "recipe": "auto"},
            "rho": {"pipeline": "auto", "seed": seed}}


def prime_curve_with_subgroup(p: int, label: str, h_ok, a: int | None = None):
    """The first labelled y^2 = x^3 + ax + b over GF(p) whose order is h·r with h_ok(h), r prime and
    r > 4·sqrt(p) (design §4.5 method 3 applies)."""
    i = 0
    while True:
        lab = labelled(label, i)
        aa = shake(f"{lab}/a", 16) % p if a is None else a % p
        b = shake(f"{lab}/b", 16) % p
        i += 1
        if (4 * aa ** 3 + 27 * b * b) % p == 0:
            continue
        c = PrimeCurve(p, aa, b)
        order = c.order(lab)
        if order is None:
            continue
        for h in range(1, 65):
            if order % h == 0 and h_ok(h) and PRIME_TEST(order // h) and order // h > 4 * math.isqrt(p) + 4:
                return c, order, order // h, h


def selftest() -> None:
    """The orders against enumeration on small fields, the field arithmetic against its own laws."""
    for p in (101, 1009, 10007):
        for a, b in ((1, 1), (2, 3), (0, 7), (p - 3, 5)):
            if (4 * a ** 3 + 27 * b * b) % p == 0:
                continue
            c = PrimeCurve(p, a, b)
            count = 1 + sum(2 if pow(v, (p - 1) // 2, p) == 1 else (1 if v == 0 else 0)
                            for v in ((x * x * x + a * x + b) % p for x in range(p)))
            assert c.order(f"selftest/{p}/{a}/{b}") == count, (p, a, b)
    k = Fp3(13, (2, 0, 1))  # t^3 + t^2 + 2: no root mod 13
    assert k.irreducible()
    u = k.element("selftest/u")
    assert k.mul(u, k.inv(u)) == k.one and k.pow(u, k.q - 1) == k.one
    s = k.mul(u, u)
    assert k.mul(k.sqrt(s), k.sqrt(s)) == s
    c = Fp3Curve(k, k.element("selftest/a"), k.element("selftest/b"))
    count = 1
    for e0 in range(13):
        for e1 in range(13):
            for e2 in range(13):
                v = c.rhs((e0, e1, e2))
                count += 1 if v == k.zero else (2 if k.is_square(v) else 0)
    assert group_order(c, k.q, (c.point(f"selftest/fp3/{i}") for i in range(64))) == count
    # The general binary form against B1's curve, on the identity map.
    f7 = v2.curve_id.find_irreducible_sparse(7)
    B = v2.Curve(7, f7, 1, 1)
    G = GeneralBinary(B.k, 1, 1, 0, 0, 1)
    P, Q = B.point("selftest/p"), B.point("selftest/q")
    assert G.add(P, Q) == B.add(P, Q) and G.add(P, P) == B.add(P, P)


# -- the cases ------------------------------------------------------------------------------------

def build() -> tuple[dict[str, dict], list[dict]]:
    files: dict[str, dict] = {}
    cases: list[dict] = []

    def case(cid: str, purpose: str, params: dict, argv_tail: list[str] | None, expect: dict,
             command: str = "price", timeout: int = 300, extra_files: dict | None = None,
             until: str | None = None) -> None:
        argv = [command, "--params", "{tmp}/params.json", "--json", "--out", "{tmp}/report.json"]
        if argv_tail:
            argv += argv_tail
        expect = {"json_file": "{tmp}/report.json", **expect}
        out = {"id": cid, "step": STEP, "checks": purpose, "files": {"params.json": params, **(extra_files or {})},
               "argv": argv, "env": {"RAYON_NUM_THREADS": "1"}, "expect": expect, "timeout_s": timeout}
        if until:
            out["until"] = until
        cases.append(out)

    def ok_run(paths: dict, contains: dict | None = None) -> dict:
        e = {"exit": 0, "json_paths": {"status": "complete", "result.verified": True, **paths}}
        if contains:
            e["json_contains"] = contains
        return e

    def refused(code: str, cls: str, status: int, contains: dict | None = None) -> dict:
        e = {"exit": status, "json_paths": {"refusal.code": code, "refusal.class": cls}}
        if contains:
            e["json_contains"] = contains
        return e

    # A prime field near 2^24 and a curve on it with a large prime subgroup.
    p24 = next_prime(1 << 24)
    c32, order32, r32, h32 = prime_curve_with_subgroup(p24, "ic-conformance-b2/C032/curve", lambda h: h <= 8)
    g32 = c32.subgroup_point("ic-conformance-b2/C032/generator", h32, r32)
    k32 = known_log("ic-conformance-b2/C032/known_log", r32)
    sw32 = {"form": "short_weierstrass", "a": str(c32.a), "b": str(c32.b)}
    doc32 = prime_doc("C032: a prime curve near 2^24, rho alone", p24, sw32, r32, h32, dec_point(g32),
                      {"known_log": str(k32)}, rho_alone())
    files["C032-prime-curve-rho.json"] = doc32
    case("C032-prime-curve-rho", "B2's prime fields: rho-negation recovers a known logarithm",
         {"copy": "{here}/C032-prime-curve-rho.json"}, None,
         ok_run({"route.rho.pipeline": "rho-negation", "result.scalar": str(k32)}))

    doc33 = copy.deepcopy(doc32)
    doc33["name"] = "C033: the same prime curve, paired"
    doc33["method"] = paired_auto()
    files["C033-prime-curve-paired.json"] = doc33
    case("C033-prime-curve-paired", "B2's prime pipelines: the study index calculus and rho on one point",
         {"copy": "{here}/C033-prime-curve-paired.json"}, ["--repeats", "1", "--repeats-fast", "1"],
         ok_run({"route.ic.pipeline": "ic-prime-s3", "route.rho.pipeline": "rho-negation",
                 "result.scalar": str(k32), "all_verified": True},
                {"disclosures": [{"code": "study-pipeline"}]}), timeout=900)

    # C034: a Montgomery curve By^2 = x^3 + Ax^2 + x, given in its own coordinates.
    i = 0
    while True:
        lab = labelled("ic-conformance-b2/C034/curve", i)
        i += 1
        A, B = shake(f"{lab}/A", 16) % p24, shake(f"{lab}/B", 16) % p24
        if B == 0 or (A * A - 4) % p24 == 0:
            continue
        inv3, invB = pow(3, -1, p24), pow(B, -1, p24)
        a_w = (3 - A * A) * inv3 * invB * invB % p24
        b_w = (2 * A ** 3 - 9 * A) * pow(27, -1, p24) * pow(B, -3, p24) % p24
        if (4 * a_w ** 3 + 27 * b_w * b_w) % p24 == 0:
            continue
        w = PrimeCurve(p24, a_w, b_w)
        order = w.order(lab)
        if order and order % 4 == 0 and PRIME_TEST(order // 4) and order // 4 > 4 * math.isqrt(p24):
            break
    r34, h34 = order // 4, 4
    gw = w.subgroup_point("ic-conformance-b2/C034/generator", h34, r34)
    # Weierstrass (x, y) = ((3u + A)/(3B), v/B), so u = Bx - A/3 and v = By.
    to_mont = lambda P: ((B * P[0] - A * inv3) % p24, B * P[1] % p24)
    gm = to_mont(gw)
    assert (B * gm[1] ** 2 - (gm[0] ** 3 + A * gm[0] ** 2 + gm[0])) % p24 == 0
    k34 = known_log("ic-conformance-b2/C034/known_log", r34)
    files["C034-montgomery-rho.json"] = prime_doc(
        "C034: a Montgomery curve near 2^24, rho alone", p24,
        {"form": "montgomery", "A": str(A), "B": str(B)}, r34, h34, dec_point(gm), {"known_log": str(k34)},
        rho_alone())
    case("C034-montgomery-rho", "B2's conversions: Montgomery to short Weierstrass, the map recorded",
         {"copy": "{here}/C034-montgomery-rho.json"}, None,
         ok_run({"result.scalar": str(k34), "conversion.from": "montgomery",
                 "conversion.to": "short_weierstrass"}))

    # C035: a twisted Edwards curve ax^2 + y^2 = 1 + dx^2y^2, given in its own coordinates.
    i = 0
    while True:
        lab = labelled("ic-conformance-b2/C035/curve", i)
        i += 1
        ae, de = shake(f"{lab}/a", 16) % p24, shake(f"{lab}/d", 16) % p24
        if ae * de * (ae - de) % p24 == 0:
            continue
        inv = pow(ae - de, -1, p24)
        A, B = 2 * (ae + de) * inv % p24, 4 * inv % p24
        if (A * A - 4) % p24 == 0:
            continue
        inv3, invB = pow(3, -1, p24), pow(B, -1, p24)
        a_w = (3 - A * A) * inv3 * invB * invB % p24
        b_w = (2 * A ** 3 - 9 * A) * pow(27, -1, p24) * pow(B, -3, p24) % p24
        if (4 * a_w ** 3 + 27 * b_w * b_w) % p24 == 0:
            continue
        w = PrimeCurve(p24, a_w, b_w)
        order = w.order(lab)
        if order and order % 4 == 0 and PRIME_TEST(order // 4) and order // 4 > 4 * math.isqrt(p24):
            break
    r35, h35 = order // 4, 4
    gw = w.subgroup_point("ic-conformance-b2/C035/generator", h35, r35)
    u, v = (B * gw[0] - A * inv3) % p24, B * gw[1] % p24  # Montgomery
    xe, ye = u * pow(v, -1, p24) % p24, (u - 1) * pow(u + 1, -1, p24) % p24  # Edwards
    assert (ae * xe * xe + ye * ye - 1 - de * xe * xe * ye * ye) % p24 == 0
    k35 = known_log("ic-conformance-b2/C035/known_log", r35)
    files["C035-twisted-edwards-rho.json"] = prime_doc(
        "C035: a twisted Edwards curve near 2^24, rho alone", p24,
        {"form": "twisted_edwards", "a": str(ae), "d": str(de)}, r35, h35, {"x": str(xe), "y": str(ye)},
        {"known_log": str(k35)}, rho_alone())
    case("C035-twisted-edwards-rho", "B2's conversions: twisted Edwards to short Weierstrass, the map "
         "recorded", {"copy": "{here}/C035-twisted-edwards-rho.json"}, None,
         ok_run({"result.scalar": str(k35), "conversion.from": "twisted_edwards",
                 "conversion.to": "short_weierstrass"}))

    # Curve A's field, for the binary general forms.
    f31 = v2.curve_id.find_irreducible_sparse(31)
    A31 = v2.Curve(31, f31, 0, 1)
    k31 = A31.k

    # C036: y^2 + y = x^3 over GF(2^31), supersingular (a1 = 0); #E = 2^31 + 1 for odd n.
    S = GeneralBinary(k31, 0, 0, 1, 0, 0)
    order36 = (1 << 31) + 1
    factors, rest = trial_factor(order36, 1 << 16)
    r36 = rest if rest > 1 else max(factors)
    h36 = order36 // r36
    assert PRIME_TEST(r36) and h36 % r36 and r36 > h36

    def lift_general(curve: GeneralBinary, label: str):
        i = 0
        while True:
            x = shake(labelled(label, i), 8) & ((1 << 31) - 1)
            i += 1
            c = curve.k.mul(curve.a1, x) ^ curve.a3
            rhs = curve.k.mul(curve.k.sqr(x), x) ^ curve.k.mul(curve.a2, curve.k.sqr(x)) \
                ^ curve.k.mul(curve.a4, x) ^ curve.a6
            if c == 0:
                continue
            z = curve.k.solve_quadratic(curve.k.mul(rhs, curve.k.inv(curve.k.sqr(c))))
            if z is None:
                continue
            P = (x, curve.k.mul(c, z))
            assert curve.on_curve(P)
            return P

    for i in range(4):
        assert S.mul(order36, lift_general(S, f"ic-conformance-b2/C036/check/{i}")) is None
    g36 = None
    i = 0
    while g36 is None:
        g36 = S.mul(h36, lift_general(S, labelled("ic-conformance-b2/C036/generator", i)))
        i += 1
    assert S.mul(r36, g36) is None
    files["C036-supersingular-binary.json"] = {
        "schema_version": 2, "name": "C036: y^2 + y = x^3 over GF(2^31), supersingular",
        "field": {"kind": "binary", "degree": 31, "modulus": v2.hx(f31)},
        "curve": {"form": "general_weierstrass", "a1": "0x0", "a2": "0x0", "a3": "0x1", "a4": "0x0", "a6": "0x0"},
        "subgroup": {"order": str(r36), "cofactor": str(h36), "generator": v2.point_doc(g36)},
        "target": {"known_log": str(known_log("ic-conformance-b2/C036/known_log", r36))},
        "method": rho_alone()}
    case("C036-supersingular-binary", "design §4.3: a binary curve with a1 = 0 is valid but unsupported",
         {"copy": "{here}/C036-supersingular-binary.json"}, None,
         refused("supersingular-binary", "unsupported", 3))

    # C037: curve A (y^2 + xy = x^3 + 1) under a change of variables (u, r, s, t), in the general form.
    c10 = json.loads((V2 / "params" / "C010-curve-a-other-generator.json").read_text())
    g_a = (int(c10["subgroup"]["generator"]["x"], 16), int(c10["subgroup"]["generator"]["y"], 16))
    order_a = v2.koblitz_order(0, 31)
    r_a = v2.largest_prime_factor(order_a)
    mask31 = (1 << 31) - 1
    u_, r_, s_, t_ = (shake(f"ic-conformance-b2/C037/{w_}", 8) & mask31 for w_ in "urst")
    assert u_
    m, sq = k31.mul, k31.sqr
    u2, u3 = sq(u_), m(sq(u_), u_)
    # Silverman Table III.1.2 in characteristic 2, from y'^2 + x'y' = x'^3 + 1 (a1' = 1, a6' = 1).
    a1 = u_
    a2 = m(s_, a1) ^ r_ ^ sq(s_)
    a3 = m(r_, a1)
    a4 = m(s_, a3) ^ m(t_ ^ m(r_, s_), a1) ^ sq(r_)
    a6 = m(m(u3, u3), 1) ^ m(r_, a4) ^ m(sq(r_), a2) ^ m(sq(r_), r_) ^ m(t_, a3) ^ sq(t_) ^ m(m(r_, t_), a1)
    E37 = GeneralBinary(k31, a1, a2, a3, a4, a6)

    def to_general(P):
        x1, y1 = P
        return m(u2, x1) ^ r_, m(u3, y1) ^ m(m(s_, u2), x1) ^ t_

    G37 = to_general(g_a)
    assert E37.on_curve(G37) and E37.mul(r_a, G37) is None
    k37 = known_log("ic-conformance-b2/C037/known_log", r_a)
    assert E37.on_curve(to_general(A31.mul(k37, g_a)))
    recipe = c10["method"]["index_calculus"]["recipe"]
    files["C037-general-binary-is-koblitz.json"] = {
        "schema_version": 2, "name": "C037: curve A under a change of variables, in the general form",
        "field": {"kind": "binary", "degree": 31, "modulus": v2.hx(f31)},
        "curve": {"form": "general_weierstrass", **{k_: v2.hx(v_) for k_, v_ in
                                                     (("a1", a1), ("a2", a2), ("a3", a3), ("a4", a4), ("a6", a6))}},
        "subgroup": {"order": str(r_a), "cofactor": str(order_a // r_a), "generator": v2.point_doc(G37)},
        "target": {"known_log": str(k37)},
        "method": {"solve": "paired", "index_calculus": {"pipeline": "auto", "recipe": recipe},
                   "rho": {"pipeline": "auto", "seed": RHO_SEED}}}
    case("C037-general-binary-is-koblitz", "B2's conversions: a general binary curve with a1 != 0 that is "
         "Koblitz after conversion runs on kic", {"copy": "{here}/C037-general-binary-is-koblitz.json"},
         ["--repeats", "1", "--repeats-fast", "1"],
         ok_run({"route.ic.pipeline": "kic", "result.scalar": str(k37), "all_verified": True,
                 "conversion.from": "general_weierstrass"}))

    # C038-C040: generic binary curves y^2 + xy = x^3 + ax^2 + b, b outside GF(2).
    def generic_binary(n: int, label: str):
        f = v2.curve_id.find_irreducible_sparse(n)
        assert v2.irreducible(f)
        i = 0
        while True:
            lab = labelled(label, i)
            i += 1
            a = shake(f"{lab}/a", 1) & 1
            b = shake(f"{lab}/b", 16) & ((1 << n) - 1)
            if b in (0, 1):
                continue
            c = v2.Curve(n, f, a, b)
            order = binary_order(c, n, lab)
            h = 2 if a == 1 else 4
            if order and order % h == 0 and PRIME_TEST(order // h) and order // h > 4 * math.isqrt(1 << n) + 4:
                return f, c, order, order // h, h

    for cid, n, solve, tail in (("C038", 23, "paired", None), ("C039", 37, "paired", None)):
        f, c, order, r, h = generic_binary(n, f"ic-conformance-b2/{cid}/curve")
        g = c.subgroup_point(f"ic-conformance-b2/{cid}/generator", h, r)
        k_ = known_log(f"ic-conformance-b2/{cid}/known_log", r)
        files[f"{cid}-generic-binary-n{n}.json"] = {
            "schema_version": 2, "name": f"{cid}: a generic binary curve over GF(2^{n}), paired",
            "field": {"kind": "binary", "degree": n, "modulus": v2.hx(f)},
            "curve": {"form": "binary_weierstrass", "a": v2.hx(c.a), "b": v2.hx(c.b)},
            "subgroup": {"order": str(r), "cofactor": str(h), "generator": v2.point_doc(g)},
            "target": {"known_log": str(k_)}, "method": paired_auto()}
        if cid == "C038":
            case("C038-generic-binary-n23-paired", "B2's binary importers: ic-binary-s4 and rho-negation on "
                 "one point", {"copy": "{here}/C038-generic-binary-n23.json"},
                 ["--repeats", "1", "--repeats-fast", "1"],
                 ok_run({"route.ic.pipeline": "ic-binary-s4", "route.rho.pipeline": "rho-negation",
                         "result.scalar": str(k_), "all_verified": True}), timeout=900)
        else:
            k39 = k_
            case("C039-generic-binary-n37-paired", "ic-binary-s4 enumerates GF(2^n) and stops at n = 32",
                 {"copy": "{here}/C039-generic-binary-n37.json"}, ["--repeats", "1", "--repeats-fast", "1"],
                 refused("no-ic-route", "unsupported", 3, {"route.considered": [
                     {"pipeline": "ic-binary-s4", "admitted": False, "gate": "enumeration-bound"}]}))
            case("C040-generic-binary-n37-rho", "the same file, rho alone, on rho-negation",
                 {"copy": "{here}/C039-generic-binary-n37.json", "set": {"method": rho_alone()}}, None,
                 ok_run({"route.rho.pipeline": "rho-negation", "result.scalar": str(k39)}))

    # C041-C043: refusals on C032's document.
    p41 = (1 << 24) + 1
    assert not PRIME_TEST(p41)
    case("C041-composite-p", "design §4.2: p must be prime", {"copy": "{here}/C032-prime-curve-rho.json",
         "set": {"field.p": str(p41)}}, None, refused("p-composite", "invalid", 2))
    a42, b42 = p24 - 3, 2
    assert (4 * a42 ** 3 + 27 * b42 * b42) % p24 == 0
    case("C042-singular-prime-curve", "design §4.3: 4a^3 + 27b^2 = 0 is singular",
         {"copy": "{here}/C032-prime-curve-rho.json",
          "set": {"curve": {"form": "short_weierstrass", "a": str(a42), "b": str(b42)}}}, None,
         refused("curve-singular", "invalid", 2))

    i = 0
    while True:
        lab = labelled("ic-conformance-b2/C043/curve", i)
        i += 1
        a, b = shake(f"{lab}/a", 16) % p24, shake(f"{lab}/b", 16) % p24
        if (4 * a ** 3 + 27 * b * b) % p24 == 0:
            continue
        c = PrimeCurve(p24, a, b)
        order = c.order(lab)
        if order is None:
            continue
        sq_primes = [q for q in (3, 5, 7, 11, 13) if order % (q * q) == 0]
        if sq_primes:
            r43 = sq_primes[-1]
            break
    h43 = order // r43
    assert h43 % r43 == 0
    g43, j = None, 0
    while g43 is None:
        g43 = c.mul(h43, c.point(labelled("ic-conformance-b2/C043/generator", j)))
        j += 1
    assert c.mul(r43, g43) is None
    files["C043-order-squared.json"] = prime_doc(
        "C043: a prime curve whose subgroup order's square divides #E", p24,
        {"form": "short_weierstrass", "a": str(a), "b": str(b)}, r43, h43, dec_point(g43),
        {"known_log": "1"}, rho_alone())
    case("C043-order-squared", "design §4.4: r dividing h is refused", {"copy": "{here}/C043-order-squared.json"},
         None, refused("order-squared", "invalid", 2))

    # C044: p near 2^34, every prime factor of #E below 4·sqrt(p): no exact method applies.
    p34 = next_prime(1 << 34)
    bound = 4 * math.isqrt(p34)
    i = 0
    while True:
        lab = labelled("ic-conformance-b2/C044/curve", i)
        i += 1
        a, b = shake(f"{lab}/a", 16) % p34, shake(f"{lab}/b", 16) % p34
        if (4 * a ** 3 + 27 * b * b) % p34 == 0:
            continue
        c = PrimeCurve(p34, a, b)
        order = c.order(lab)
        if order is None:
            continue
        factors, rest = trial_factor(order, bound)
        if rest == 1 and factors and factors.count(max(factors)) == 1:
            break
    r44 = max(factors)
    h44 = order // r44
    assert r44 < bound and PRIME_TEST(r44) and h44 % r44
    g44 = c.subgroup_point("ic-conformance-b2/C044/generator", h44, r44)
    files["C044-cardinality-unknown.json"] = prime_doc(
        "C044: a curve over p near 2^34 whose subgroup is below 4 sqrt(p)", p34,
        {"form": "short_weierstrass", "a": str(a), "b": str(b)}, r44, h44, dec_point(g44),
        {"known_log": str(known_log("ic-conformance-b2/C044/known_log", r44))}, rho_alone())
    case("C044-cardinality-unknown", "design §4.5: no exact method finds #E, so the instance is unsupported",
         {"copy": "{here}/C044-cardinality-unknown.json"}, None,
         refused("cardinality-unknown", "unsupported", 3))

    # C045: an anomalous curve, #E = p, by complex multiplication with D = -3: 4p - 1 = 3v^2.
    v = math.isqrt((4 << 20) // 3) | 1
    while not PRIME_TEST((3 * v * v + 1) // 4):
        v += 2
    p45 = (3 * v * v + 1) // 4
    b45 = 1
    while True:
        c = PrimeCurve(p45, 0, b45)
        if c.order(f"ic-conformance-b2/C045/{b45}") == p45:
            break
        b45 += 1
    g45 = c.point("ic-conformance-b2/C045/generator")
    assert c.mul(p45, g45) is None
    k45 = known_log("ic-conformance-b2/C045/known_log", p45)
    files["C045-anomalous.json"] = prime_doc(
        "C045: an anomalous curve, #E = p, near 2^20", p45, {"form": "short_weierstrass", "a": "0", "b": str(b45)},
        p45, 1, dec_point(g45), {"known_log": str(k45)}, rho_alone())
    case("C045-anomalous", "design §4.6: an anomalous curve runs, with the disclosure",
         {"copy": "{here}/C045-anomalous.json"}, None,
         ok_run({"route.rho.pipeline": "rho-negation", "result.scalar": str(k45)},
                {"disclosures": [{"code": "anomalous"}]}))

    # C046: a curve over GF(2^4) taken over GF(2^100), a subgroup near 2^20, on rho-bignum.
    n46 = 100
    f100 = v2.curve_id.find_irreducible_sparse(n46)
    K = v2.Field(n46, f100)
    gamma, beta = 2, None
    while beta is None:
        e, base, cand = ((1 << n46) - 1) // 15, gamma, 1
        while e:
            if e & 1:
                cand = K.mul(cand, base)
            base = K.sqr(base)
            e >>= 1
        powers = [1]
        for _ in range(15):
            powers.append(K.mul(powers[-1], cand))
        if powers[15] == 1 and powers[3] != 1 and powers[5] != 1:
            beta = cand
        gamma += 1
    sub = [0] + powers[:15]  # GF(2^4) inside GF(2^100)

    def tr16(cv: int) -> int:
        t, s = cv, cv
        for _ in range(3):
            s = K.sqr(s)
            t ^= s
        assert t in (0, 1)
        return t

    found = None
    for ia in range(16):
        for ib in range(16):
            a, b = sub[ia], sub[ib]
            if b == 0 or K.sqr(K.sqr(b)) == b:  # b in GF(2^2): the least subfield would be smaller
                continue
            count = 2  # O and (0, sqrt b)
            for x in sub[1:]:
                cx = K.mul(K.mul(K.mul(x, x), x) ^ K.mul(a, K.sqr(x)) ^ b, K.inv(K.sqr(x)))
                count += 2 if tr16(cx) == 0 else 0
            t = 17 - count
            s_prev, s_cur = 2, t
            for _ in range(24):
                s_prev, s_cur = s_cur, t * s_cur - 16 * s_prev
            N = (1 << n46) + 1 - s_cur
            # The prime factors near 2^20: split off small ones, then look for one by Pollard-Brent.
            factors, rest = trial_factor(N, 1 << 12)
            pending = [rest] if rest > 1 else []
            while pending and not found:
                m_ = pending.pop()
                if m_ <= (1 << 21) and PRIME_TEST(m_):
                    if (1 << 19) <= m_ and N % (m_ * m_):
                        found = (a, b, N, m_)
                    continue
                if m_ < v2.PSI13 and PRIME_TEST(m_):
                    continue
                g_ = brent_factor(m_, f"ic-conformance-b2/C046/{ia}/{ib}/{m_}", 1 << 14)
                if g_:
                    pending += [g_, m_ // g_]
            if found:
                break
        if found:
            break
    a46, b46, N46, r46 = found
    h46 = N46 // r46
    C46 = v2.Curve(n46, f100, a46, b46)
    for i in range(3):
        assert C46.mul(N46, C46.point(f"ic-conformance-b2/C046/check/{i}")) is None
    g46 = C46.subgroup_point("ic-conformance-b2/C046/generator", h46, r46)
    k46 = known_log("ic-conformance-b2/C046/known_log", r46)
    files["C046-subfield-over-gf2-100-bignum.json"] = {
        "schema_version": 2, "name": "C046: a curve over GF(2^4) taken over GF(2^100), rho alone",
        "field": {"kind": "binary", "degree": n46, "modulus": v2.hx(f100)},
        "curve": {"form": "binary_weierstrass", "a": v2.hx(a46), "b": v2.hx(b46)},
        "subgroup": {"order": str(r46), "cofactor": str(h46), "generator": v2.point_doc(g46)},
        "target": {"known_log": str(k46)}, "method": rho_alone()}
    case("C046-subfield-curve-on-rho-bignum", "B2's rho-bignum takes what no one-word rho admits",
         {"copy": "{here}/C046-subfield-over-gf2-100-bignum.json"}, None,
         ok_run({"route.rho.pipeline": "rho-bignum", "result.scalar": str(k46)}), timeout=900)

    # C047: secp256k1 by name, a target point, rho alone: over budget.
    p256 = (1 << 256) - (1 << 32) - 977
    n256 = int("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141", 16)
    G256 = (int("79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798", 16),
            int("483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8", 16))
    secp = PrimeCurve(p256, 0, 7)
    assert secp.on_curve(G256) and secp.mul(n256, G256) is None
    Q256 = secp.mul(known_log("ic-conformance-b2/C047/scalar", n256), G256)
    files["C047-secp256k1-over-budget.json"] = {
        "schema_version": 2, "name": "C047: secp256k1 by name, a target point, rho alone",
        "named": "secp256k1", "target": {"point": {"x": v2.hx(Q256[0]), "y": v2.hx(Q256[1])}},
        "method": rho_alone()}
    steps_bits = math.ceil(0.5 * math.log2(math.pi * n256 / 4))
    assert steps_bits == 128
    case("C047-secp256k1-over-budget", "design §5.4: the estimate, near 2^128 steps, exceeds the budget",
         {"copy": "{here}/C047-secp256k1-over-budget.json"}, None,
         {"exit": 4, "json_paths": {"refusal.code": "over-budget", "refusal.class": "over_budget",
                                    "estimate.rho.steps_bits": steps_bits}})

    # C048: the suite row at 2^47.2, translated, with a one-second budget.
    row_rel = "params/S/k0n61/M1-T01.json"
    row = v2.v1_row(row_rel)
    f61 = v2.curve_id.find_irreducible_sparse(61)
    curve61 = v2.Curve(61, f61, 0, 1)
    order61 = v2.koblitz_order(0, 61)
    r61 = v2.largest_prime_factor(order61)
    assert v2.curve_id.binary_id(61, f61, 0, 1, order61, end="-7")["slug"] == "icv1-f2m61-t158598901-ab42b6c5"
    doc48 = v2.document("C048: suite v1 row S/k0n61/M1-T01, translated, with a one-second budget", curve61,
                        {"form": "koblitz", "a": 0}, r61, order61 // r61, {"rule": "koblitz_search_v1"},
                        {"public_hash_seed": row["targets"][0]["public_hash_seed"]}, v2.recipe_of(row), RHO_SEED)
    doc48["budget"] = {"wall_seconds": 1}
    files["C048-suite-row-over-budget.json"] = doc48
    case("C048-suite-row-over-budget", "design §3.6 and §5.4: a run whose estimate exceeds its budget is "
         "refused before it starts", {"copy": "{here}/C048-suite-row-over-budget.json"},
         ["--repeats", "1", "--repeats-fast", "1"], {"exit": 4, "json_paths": {
             "refusal.code": "over-budget", "refusal.class": "over_budget"}})

    # C049: the gate file at F2: estimates only.
    case("C049-gate-estimates", "design §5.2 step 5: fidelity F2 runs nothing and estimates both arms",
         {"copy": "{cases}/gate-m83-T001.json", "set": {"method.fidelity": "F2"}},
         ["--repeats", "1", "--repeats-fast", "1"],
         {"exit": 0, "json_paths": {"status": "estimated", "estimate.ic.pipeline": "kic",
                                    "estimate.rho.pipeline": "rho-koblitz"}})

    # C050-C051: GF(p^3).
    p3 = 1009
    i = 0
    while True:
        cs = tuple(shake(f"{labelled('ic-conformance-b2/C050/modulus', i)}/{j}", 8) % p3 for j in range(3))
        i += 1
        K3 = Fp3(p3, cs)
        if K3.irreducible():
            break
    i = 0
    while True:
        lab = labelled("ic-conformance-b2/C050/curve", i)
        i += 1
        a3, b3 = K3.element(f"{lab}/a"), K3.element(f"{lab}/b")
        try:
            E3 = Fp3Curve(K3, a3, b3)
            order = group_order(E3, K3.q, (E3.point(f"{lab}/order/{j}") for j in range(64)))
        except (AssertionError, ValueError):
            continue
        factors, rest = trial_factor(order, 1 << 12)
        if rest > 4 * math.isqrt(K3.q) and PRIME_TEST(rest) and order % (rest * rest):
            break
    r50, h50 = rest, order // rest
    g50, j = None, 0
    while g50 is None:
        g50 = E3.mul(h50, E3.point(labelled("ic-conformance-b2/C050/generator", j)))
        j += 1
    assert E3.mul(r50, g50) is None
    files["C050-cubic-extension.json"] = {
        "schema_version": 2, "name": "C050: a curve over GF(1009^3)",
        "field": {"kind": "prime_extension", "p": str(p3), "degree": 3, "modulus": arr(cs)},
        "curve": {"form": "short_weierstrass", "a": arr(a3), "b": arr(b3)},
        "subgroup": {"order": str(r50), "cofactor": str(h50), "generator": {"x": arr(g50[0]), "y": arr(g50[1])}},
        "target": {"known_log": str(known_log("ic-conformance-b2/C050/known_log", r50))},
        "method": rho_alone()}
    case("C050-cubic-extension", "design §3.1: GF(p^k) is parsed and validated, and routed only from B5",
         {"copy": "{here}/C050-cubic-extension.json"}, None,
         refused("no-pipeline-for-field", "unsupported", 3), until="B5")
    reducible = (p3 - 1, 0, 0)  # t^3 - 1 has the root 1
    assert not Fp3(p3, reducible).irreducible()
    case("C051-reducible-cubic-modulus", "design §4.2: an extension's modulus must be irreducible",
         {"copy": "{here}/C050-cubic-extension.json", "set": {"field.modulus": arr(reducible)}}, None,
         refused("modulus-reducible", "invalid", 2))
    return files, cases


def texts() -> dict[str, str]:
    selftest()
    files, cases = build()
    out = {f"params/{name}": json.dumps(doc, indent=1) + "\n" for name, doc in files.items()}
    out["cases.json"] = json.dumps({
        "suite": "ic tool programme conformance suite v2, B2's cases",
        "design": "research/ic_tool_program/design/schema-v2.md §9 (C032-C051); "
                  "research/ic_tool_program/rounds/B2-fields-forms-estimates/PROTOCOL.md",
        "includes": "../v1/cases.json (B0) and ../v2/cases.json (B1); ../run.py runs every step's cases",
        "rules": [
            "The rules of ../v2/cases.json hold. {here} is this directory's params/; {cases} is "
            "../v2/params/.",
            "B2's expectations are the design's §9 rows, made exact here: the codes, the exit statuses, "
            "the pipelines, and the report keys B2 adds (conversion.from and .to, estimate.<arm>.pipeline "
            "and .steps_bits, the disclosure study-pipeline).",
        ],
        "generator": {"path": "research/ic_tool_program/conformance/v2-b2/make_cases.py",
                      "sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},
        "cases": cases,
    }, indent=1, ensure_ascii=False) + "\n"
    out["SHA256SUMS"] = "".join(f"{hashlib.sha256(t.encode()).hexdigest()}  {rel}\n"
                                for rel, t in sorted(out.items()))
    return out


def main() -> None:
    out = texts()
    if "--check" in sys.argv:
        bad = [rel for rel, t in out.items() if not (HERE / rel).exists() or (HERE / rel).read_text() != t]
        print(json.dumps({"files": len(out), "mismatches": bad}, indent=1))
        raise SystemExit(1 if bad else 0)
    if (HERE / "cases.json").exists():
        raise SystemExit("cases.json exists; B2's cases are frozen (use --check)")
    for rel, t in out.items():
        path = HERE / rel
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, "x") as fh:
            fh.write(t)
    print(f"{len(out)} files written")


if __name__ == "__main__":
    main()
