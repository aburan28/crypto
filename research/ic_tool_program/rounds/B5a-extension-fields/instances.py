#!/usr/bin/env python3
"""B5a's instances (`../../design/extension-fields.md` §4): curves over
`GF(p^k)`, found in Python arithmetic that shares nothing with the Rust
tool.

    python3 instances.py G1 E2 ...     # searches the instances named, one JSON record a line
    python3 instances.py --assemble FILES... > instances.json
    python3 instances.py --check       # repeats every search and compares with instances.json (slow)

Every choice comes from SHAKE-256 of a public label, with `:1`, `:2`, …
appended until it works, as in B2's generator, whose helpers this reuses
(`../../conformance/v2-b2/make_cases.py`):
- **the modulus:** `t³ − c` with `c` the least non-cube for `G` and `H`;
  otherwise the first labelled monic polynomial that Rabin's test finds
  irreducible;
- **the curve:** `y² = x³ + ax + b`, `a` and `b` labelled, the first that
  meets the instance's rule;
- **`#E`:** found exactly, by baby steps and giant steps over the Hasse
  interval, narrowed by the orders of labelled points until one value is
  left (B2's `group_order`);
- **`r`:** for `G`, `#E` itself, prime. Otherwise `r` is `#E`'s largest
  prime factor, which must lie between `4√q` and `2^44` with `r²` not
  dividing `#E`; `H` also needs `h > 1`. Factors come from trial division
  and Brent's method, each certified by B1's exact test.

The search is slow in Python (minutes at `q ≈ 2^70`), so each instance
is searched on its own and the records are assembled. `make_cases.py`
re-checks every record by §4.5's method 3, which is cheap.
"""
from __future__ import annotations

import importlib.util
import json
import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
_spec = importlib.util.spec_from_file_location(
    "conformance_v2_b2_cases", HERE.parents[1] / "conformance" / "v2-b2" / "make_cases.py")
b2 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(b2)
_spec = importlib.util.spec_from_file_location("curve_id", ROOT / "scripts" / "curve_id.py")
curve_id = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(curve_id)

LABEL = "ic-tool-programme/B5a/instances"
R_MAX = 1 << 44


# -- GF(p)[t] and GF(p^k) -------------------------------------------------------------------------

def _trim(a: list[int]) -> list[int]:
    while a and a[-1] == 0:
        a.pop()
    return a


def _pmul(a: list[int], b: list[int], p: int) -> list[int]:
    if not a or not b:
        return []
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                out[i + j] += x * y
    return _trim([v % p for v in out])


def _psub(a: list[int], b: list[int], p: int) -> list[int]:
    n = max(len(a), len(b))
    a, b = a + [0] * (n - len(a)), b + [0] * (n - len(b))
    return _trim([(x - y) % p for x, y in zip(a, b)])


def _pdivmod(a: list[int], b: list[int], p: int) -> tuple[list[int], list[int]]:
    a, q = a[:], [0] * max(len(a) - len(b) + 1, 1)
    lead = pow(b[-1], p - 2, p)
    while len(a) >= len(b) and a:
        c = a[-1] * lead % p
        s = len(a) - len(b)
        q[s] = c
        for i, y in enumerate(b):
            a[s + i] = (a[s + i] - c * y) % p
        _trim(a)
    return _trim(q), a


def _pgcd(a: list[int], b: list[int], p: int) -> list[int]:
    while b:
        a, b = b, _pdivmod(a, b, p)[1]
    return a


class Fpk:
    """GF(p)[t] / (t^k + c_{k-1} t^{k-1} + ... + c_0); an element is a k-tuple (e_0, ..., e_{k-1})."""

    def __init__(self, p: int, mod: list[int]):
        self.p, self.k = p, len(mod)
        self.mod = [c % p for c in mod]
        self.f = self.mod + [1]
        self.q = p ** self.k
        self.zero = (0,) * self.k
        self.one = (1,) + (0,) * (self.k - 1)

    def reduce(self, u: list[int]) -> tuple:
        return tuple(curve_id._fpk_reduce(list(u) or [0], self.p, self.mod))

    def add(self, u, v):
        return tuple((x + y) % self.p for x, y in zip(u, v))

    def sub(self, u, v):
        return tuple((x - y) % self.p for x, y in zip(u, v))

    def mul(self, u, v):
        return self.reduce(_pmul(list(u), list(v), self.p))

    def scalar(self, c: int):
        return (c % self.p,) + (0,) * (self.k - 1)

    def pow(self, u, e: int):
        out = self.one
        for bit in bin(e)[2:] if e > 0 else "":
            out = self.mul(out, out)
            if bit == "1":
                out = self.mul(out, u)
        return out

    def inv(self, u):
        """The extended Euclidean algorithm in GF(p)[t]: s u = 1 mod f."""
        p = self.p
        r0, r1, s0, s1 = self.f[:], _trim(list(u)), [], [1]
        assert r1, "the inverse of zero"
        while r1:
            quo, rem = _pdivmod(r0, r1, p)
            r0, r1 = r1, rem
            s0, s1 = s1, _psub(s0, _pmul(quo, s1, p), p)
        assert len(r0) == 1, "the modulus is reducible"
        c = pow(r0[0], p - 2, p)
        return self.reduce([x * c for x in s0])

    def is_square(self, u) -> bool:
        return u == self.zero or self.pow(u, (self.q - 1) // 2) == self.one

    def non_residue(self):
        """The first labelled non-square, the same for every root this field takes."""
        if getattr(self, "_z", None) is None:
            i = 0
            while self.is_square(z := self.element(b2.labelled(f"{LABEL}/non-residue", i))):
                i += 1
            self._z = z
        return self._z

    def sqrt(self, u):
        """A square root by Tonelli-Shanks, or None."""
        if u == self.zero:
            return self.zero
        if not self.is_square(u):
            return None
        m, s = self.q - 1, 0
        while m % 2 == 0:
            m //= 2
            s += 1
        z = self.non_residue()
        c, t, r = self.pow(z, m), self.pow(u, m), self.pow(u, (m + 1) // 2)
        while t != self.one:
            j, t2 = 0, t
            while t2 != self.one:
                t2 = self.mul(t2, t2)
                j += 1
            bb = self.pow(c, 1 << (s - j - 1))
            s, c, t, r = j, self.mul(bb, bb), self.mul(t, self.mul(bb, bb)), self.mul(r, bb)
        return r

    def element(self, label: str):
        return tuple(b2.shake(f"{label}/{i}", 16) % self.p for i in range(self.k))

    def pack(self, u) -> int:
        """Σ e_i p^i: the integer the tool's keys pack an element into (design §2.2)."""
        return sum(x * self.p ** i for i, x in enumerate(u))


def rabin_irreducible(p: int, mod: list[int]) -> bool:
    """Rabin's test: t^(p^k) = t mod f, and gcd(t^(p^(k/l)) - t, f) = 1 for every prime l | k."""
    f, k = [c % p for c in mod] + [1], len(mod)

    def t_pow(e: int) -> list[int]:
        out, base = [1], [0, 1]
        while e:
            if e & 1:
                out = _pdivmod(_pmul(out, base, p), f, p)[1]
            base = _pdivmod(_pmul(base, base, p), f, p)[1]
            e >>= 1
        return out

    if _psub(t_pow(p ** k), [0, 1], p):
        return False
    for ell in {d for d in range(2, k + 1) if k % d == 0 and all(d % e for e in range(2, d))}:
        if len(_pgcd(f, _psub(t_pow(p ** (k // ell)), [0, 1], p), p)) != 1:
            return False
    return True


class Curve:
    """y^2 = x^3 + ax + b over an Fpk; a point is (x, y) or None, the identity."""

    def __init__(self, F: Fpk, a, b):
        self.F, self.a, self.b = F, a, b

    def singular(self) -> bool:
        F = self.F
        return F.add(F.mul(F.scalar(4), F.mul(self.a, F.mul(self.a, self.a))),
                     F.mul(F.scalar(27), F.mul(self.b, self.b))) == F.zero

    def rhs(self, x):
        F = self.F
        return F.add(F.add(F.mul(x, F.mul(x, x)), F.mul(self.a, x)), self.b)

    def on_curve(self, P) -> bool:
        return P is None or self.F.mul(P[1], P[1]) == self.rhs(P[0])

    def neg(self, P):
        return None if P is None else (P[0], self.F.sub(self.F.zero, P[1]))

    def add(self, P, Q):
        F = self.F
        if P is None:
            return Q
        if Q is None:
            return P
        (x1, y1), (x2, y2) = P, Q
        if x1 == x2:
            if F.add(y1, y2) == F.zero:
                return None
            lam = F.mul(F.add(F.mul(F.scalar(3), F.mul(x1, x1)), self.a), F.inv(F.mul(F.scalar(2), y1)))
        else:
            lam = F.mul(F.sub(y2, y1), F.inv(F.sub(x2, x1)))
        x3 = F.sub(F.sub(F.mul(lam, lam), x1), x2)
        return x3, F.sub(F.mul(lam, F.sub(x1, x3)), y1)

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
            x = self.F.element(b2.labelled(label, i))
            y = self.F.sqrt(self.rhs(x))
            if y is not None:
                return x, y
            i += 1


# -- the instances --------------------------------------------------------------------------------

def least_non_cube(p: int) -> int:
    assert p % 3 == 1
    return next(c for c in range(2, p) if pow(c, (p - 1) // 3, p) != 1)


def labelled_modulus(label: str, p: int, k: int) -> tuple[list[int], int]:
    """The first labelled monic modulus of degree k that Rabin's test finds irreducible, and its index."""
    i = 0
    while True:
        lab = b2.labelled(f"{label}/modulus", i)
        mod = [b2.shake(f"{lab}/{j}", 16) % p for j in range(k)]
        if mod[0] and rabin_irreducible(p, mod):
            return mod, i
        i += 1


def factor(n: int, label: str) -> list[int]:
    """The prime factors of n with multiplicity, each certified by B1's exact test."""
    small, rest = b2.trial_factor(n, 1 << 16)
    out, todo = list(small), [rest] if rest > 1 else []
    while todo:
        m = todo.pop()
        if b2.PRIME_TEST(m):
            out.append(m)
            continue
        d = b2.brent_factor(m, f"{label}/{m}", 1 << 22)
        if d is None:
            raise ValueError(f"could not factor {m}")
        todo += [d, m // d]
    return sorted(out)


def subgroup(order: int, q: int, rule: str, label: str):
    """(r, h) by the instance's rule, or None."""
    floor = 4 * math.isqrt(q) + 4
    if rule == "prime":
        return (order, 1) if b2.PRIME_TEST(order) and order > floor else None
    fs = factor(order, label)
    r = fs[-1]
    if not (floor < r < R_MAX) or order % (r * r) == 0:
        return None
    h = order // r
    if rule == "cofactor" and h == 1:
        return None
    return r, h


# id: (p, k, modulus rule, subgroup rule, purpose)
SPECS = {
    "G1": (271, 3, "binomial", "prime", "ic-gaudry-cubic at n ≈ 2^24 (§11.2's first prime)"),
    "G2": (523, 3, "binomial", "prime", "ic-gaudry-cubic at n ≈ 2^27"),
    "G3": (1039, 3, "binomial", "prime", "ic-gaudry-cubic at n ≈ 2^30"),
    "H1": (271, 3, "binomial", "cofactor", "a cofactor: refused by ic-gaudry-cubic"),
    "E2": ((1 << 31) - 1, 2, "labelled", "largest", "the one-word arithmetic at its widest, q ≈ 2^62"),
    "E5": (2053, 5, "labelled", "largest", "a general degree, q ≈ 2^55"),
    "E11": (37, 11, "labelled", "largest", "past eight coefficients, q ≈ 2^57"),
    "B2": (34359738421, 2, "labelled", "largest", "q ≈ 2^70, past one word: rho-bignum"),
}


def search(iid: str) -> dict:
    p, k, mod_rule, sub_rule, purpose = SPECS[iid]
    assert b2.PRIME_TEST(p)
    label = f"{LABEL}/{iid}"
    if mod_rule == "binomial":
        c = least_non_cube(p)
        mod, mod_index = [(-c) % p] + [0] * (k - 1), None
    else:
        mod, mod_index = labelled_modulus(label, p, k)
    assert rabin_irreducible(p, mod)
    F = Fpk(p, mod)
    i = 0
    while True:
        lab = b2.labelled(f"{label}/curve", i)
        E = Curve(F, F.element(f"{lab}/a"), F.element(f"{lab}/b"))
        i += 1
        if E.singular():
            continue
        try:
            order = b2.group_order(E, F.q, (E.point(f"{lab}/order/{j}") for j in range(64)))
        except ValueError:  # the points' orders left more than one value: not cyclic enough
            continue
        rh = subgroup(order, F.q, sub_rule, lab)
        if rh is None:
            continue
        r, h = rh
        ident = curve_id.extension_id(p, k, mod, list(E.a), list(E.b), order)
        return {"id": iid, "purpose": purpose, "p": str(p), "k": k, "modulus": [str(c) for c in mod],
                "modulus_rule": mod_rule, "modulus_index": mod_index, "curve_index": i - 1,
                "a": [str(x) for x in E.a], "b": [str(x) for x in E.b], "order": str(order),
                "r": str(r), "h": str(h), "log2_q": round(math.log2(F.q), 2), "log2_r": round(math.log2(r), 2),
                "slug": ident["slug"], "icv1": ident["icv1"]}


def main() -> None:
    args = sys.argv[1:]
    if args[:1] == ["--assemble"]:
        recs = {}
        for path in args[1:]:
            for line in Path(path).read_text().splitlines():
                if line.strip():
                    rec = json.loads(line)
                    recs[rec["id"]] = rec
        missing = [i for i in SPECS if i not in recs]
        if missing:
            raise SystemExit(f"missing records: {missing}")
        print(json.dumps({"label": LABEL, "generator": "research/ic_tool_program/rounds/B5a-extension-fields/instances.py",
                          "instances": [recs[i] for i in SPECS]}, indent=1, ensure_ascii=False))
        return
    if args[:1] == ["--check"]:
        want = json.loads((HERE / "instances.json").read_text())["instances"]
        bad = [rec["id"] for rec in want if search(rec["id"]) != rec]
        print(json.dumps({"instances": len(want), "mismatches": bad}))
        raise SystemExit(1 if bad else 0)
    for iid in args or list(SPECS):
        print(json.dumps(search(iid), ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()
