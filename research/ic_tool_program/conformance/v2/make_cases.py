#!/usr/bin/env python3
"""Conformance suite v2: B1's cases and their parameter files
(`research/ic_tool_program/design/schema-v2.md` §9).

    python3 make_cases.py           # writes params/, cases.json and SHA256SUMS; refuses to overwrite
    python3 make_cases.py --check   # re-derives every file and compares the bytes

The arithmetic here is its own. Polynomials over GF(2) are Python
integers (bit i is the coefficient of z^i), and points use the affine
group law of y² + xy = x³ + ax² + b. None of it is shared with the Rust
tool the cases test. Each file's intended property is asserted when the
file is built: whether the point is on the curve, whether it lies in the
subgroup, and whether the cofactor is right or deliberately wrong. A case
therefore cannot pass because the tool and its test share a mistake.
ICV1 strings come from `scripts/curve_id.py`, the reference
implementation.

Points in the files are chosen by one rule, so that anyone can re-derive
them:
- **A field element from a label.** Take SHAKE-256 of the label, as many
  bytes as the degree needs, and keep the low `n` bits.
- **A point.** Lift the abscissa, with the root `z` of `z² + z = c`
  whose bit 0 is 0. If there is none, append `:1`, `:2`, … to the label
  and retry.
- **A point of order r.** Multiply that point by `h`. Retry if the
  result is `O`.
"""
from __future__ import annotations

import copy
import hashlib
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
SUITE = ROOT / "research" / "ic_tool_program" / "suite" / "v1"
sys.path.insert(0, str(ROOT / "scripts"))
import curve_id  # noqa: E402

# Miller-Rabin with the first 13 primes as bases is exact below this
# bound (Sorenson and Webster, Math. Comp. 86 (2017) 985-1003).
PSI13 = 3317044064679887385961981
BASES13 = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41)


# -- GF(2^n) ------------------------------------------------------------------------------------

class Field:
    """GF(2)[z] / (f), elements as integers below 2^n."""

    def __init__(self, n: int, f: int):
        assert f.bit_length() == n + 1 and f & 1, "the modulus has degree n and a constant term"
        self.n, self.f = n, f
        self._kernel = None

    def mul(self, a: int, b: int) -> int:
        r = 0
        while b:
            if b & 1:
                r ^= a
            b >>= 1
            a <<= 1
            if a >> self.n & 1:
                a ^= self.f
        return r

    def sqr(self, a: int) -> int:
        return self.mul(a, a)

    def inv(self, a: int) -> int:
        """Algorithm 2.48 of Hankerson, Menezes and Vanstone."""
        assert a
        u, v, g1, g2 = a, self.f, 1, 0
        while u != 1:
            j = u.bit_length() - v.bit_length()
            if j < 0:
                u, v, g1, g2, j = v, u, g2, g1, -j
            u ^= v << j
            g1 ^= g2 << j
        return g1

    def trace(self, c: int) -> int:
        t, s = 0, c
        for _ in range(self.n):
            t ^= s
            s = self.sqr(s)
        assert t in (0, 1)
        return t

    def sqrt(self, a: int) -> int:
        for _ in range(self.n - 1):
            a = self.sqr(a)
        return a

    def solve_quadratic(self, c: int) -> int | None:
        """The root of z² + z = c whose bit 0 is 0, by linear algebra over GF(2), for either parity of n."""
        if self.trace(c):
            return None
        if self._kernel is None:
            basis = []  # (vector, combination): leading bits distinct, in descending order
            for i in range(self.n):
                v, combo = self.sqr(1 << i) ^ (1 << i), 1 << i
                for bv, bc in basis:
                    if v ^ bv < v:
                        v, combo = v ^ bv, combo ^ bc
                if v:
                    basis.append((v, combo))
                    basis.sort(reverse=True)
            self._kernel = basis
        rest, z = c, 0
        for bv, bc in self._kernel:
            if rest ^ bv < rest:
                rest, z = rest ^ bv, z ^ bc
        assert rest == 0
        z ^= z & 1
        assert self.sqr(z) ^ z == c
        return z


def irreducible(f: int) -> bool:
    """Rabin's test over GF(2)."""
    n = f.bit_length() - 1
    k = Field(n, f)

    def x_pow_2k(e: int) -> int:
        y = 2
        for _ in range(e):
            y = k.sqr(y)
        return y

    def gcd(a: int, b: int) -> int:
        while b:
            while a and a.bit_length() >= b.bit_length():
                a ^= b << (a.bit_length() - b.bit_length())
            a, b = b, a
        return a

    if x_pow_2k(n) != 2:
        return False
    primes = [p for p in range(2, n + 1) if n % p == 0 and all(p % q for q in range(2, p))]
    return all(gcd(x_pow_2k(n // p) ^ 2, f) == 1 for p in primes)


# -- the curve y² + xy = x³ + ax² + b -------------------------------------------------------

class Curve:
    def __init__(self, n: int, f: int, a: int, b: int):
        self.k = Field(n, f)
        self.n, self.f, self.a, self.b = n, f, a, b
        assert b, "non-singular"

    def on_curve(self, p) -> bool:
        if p is None:
            return True
        x, y = p
        k = self.k
        return k.sqr(y) ^ k.mul(x, y) == k.mul(k.sqr(x), x) ^ k.mul(self.a, k.sqr(x)) ^ self.b

    def add(self, p, q):
        if p is None:
            return q
        if q is None:
            return p
        k = self.k
        (x1, y1), (x2, y2) = p, q
        if x1 == x2:
            if y1 ^ y2 == x1:
                return None
            return self.dbl(p)
        lam = k.mul(y1 ^ y2, k.inv(x1 ^ x2))
        x3 = k.sqr(lam) ^ lam ^ x1 ^ x2 ^ self.a
        return x3, k.mul(lam, x1 ^ x3) ^ x3 ^ y1

    def dbl(self, p):
        if p is None or p[0] == 0:
            return None
        k = self.k
        x, y = p
        lam = x ^ k.mul(y, k.inv(x))
        x3 = k.sqr(lam) ^ lam ^ self.a
        return x3, k.sqr(x) ^ k.mul(lam ^ 1, x3)

    def mul(self, e: int, p):
        r = None
        for bit in bin(e)[2:] if e else "":
            r = self.dbl(r)
            if bit == "1":
                r = self.add(r, p)
        return r

    def lift(self, x: int):
        k = self.k
        if x == 0:
            return 0, k.sqrt(self.b)
        c = x ^ self.a ^ k.mul(self.b, k.inv(k.sqr(x)))
        z = k.solve_quadratic(c)
        if z is None:
            return None
        p = (x, k.mul(x, z))
        assert self.on_curve(p)
        return p

    def element(self, label: str) -> int:
        digest = hashlib.shake_256(label.encode()).digest((self.n + 7) // 8)
        return int.from_bytes(digest, "big") & ((1 << self.n) - 1)

    def point(self, label: str):
        for i in range(10_000):
            p = self.lift(self.element(label if i == 0 else f"{label}:{i}"))
            if p is not None:
                return p
        raise AssertionError("no point")

    def subgroup_point(self, label: str, h: int, r: int):
        for i in range(10_000):
            g = self.mul(h, self.point(label if i == 0 else f"{label}#{i}"))
            if g is not None:
                assert self.mul(r, g) is None
                return g
        raise AssertionError("no subgroup point")


def koblitz_order(a: int, n: int) -> int:
    t = -1 if a == 0 else 1
    prev, cur = 2, t
    for _ in range(1, n):
        prev, cur = cur, t * cur - 2 * prev
    return (1 << n) + 1 - cur


def miller_rabin(n: int, bases) -> bool:
    if n < 2:
        return False
    for p in BASES13:
        if n % p == 0:
            return n == p
    d, s = n - 1, 0
    while d % 2 == 0:
        d, s = d // 2, s + 1
    for b in bases:
        x = pow(b, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def prime_exact(n: int) -> bool:
    assert n < PSI13, "beyond the exact range"
    return miller_rabin(n, BASES13)


def largest_prime_factor(m: int) -> int:
    d, last = 2, 1
    while d * d <= m:
        while m % d == 0:
            m, last = m // d, d
        d += 1 if d == 2 else 2
    return max(m, last) if m > 1 else last


# -- documents ------------------------------------------------------------------------------

def hx(v: int) -> str:
    return f"0x{v:x}"


def point_doc(p) -> dict:
    return {"x": hx(p[0]), "y": hx(p[1])}


def v1_row(rel: str) -> dict:
    return json.loads((SUITE / rel).read_text())


RECIPE_KEYS = ("summands", "descent_summands", "collection_window", "collection_aim", "pair_table_bytes",
               "pair_table_tier", "solver", "seed", "max_trials", "linear_algebra", "collection", "factor_base")


def recipe_of(v1: dict) -> dict:
    return {k: v1[k] for k in RECIPE_KEYS if k in v1}


def document(name: str, curve: Curve, form: dict, r: int, h: int, generator, target: dict,
             recipe: dict | str, rho_seed: int) -> dict:
    field = {"kind": "binary", "degree": curve.n, "modulus": hx(curve.f)}
    return {
        "schema_version": 2,
        "name": name,
        "field": field,
        "curve": form,
        "subgroup": {"order": str(r), "cofactor": str(h),
                     "generator": generator if isinstance(generator, dict) else point_doc(generator)},
        "target": target,
        "method": {"solve": "paired",
                   "index_calculus": {"pipeline": "kic", "recipe": recipe},
                   "rho": {"pipeline": "rho-koblitz", "seed": rho_seed}},
    }


def icv1(curve: Curve, order: int) -> str:
    return curve_id.binary_id(curve.n, curve.f, curve.a, curve.b, order, end="-7")["icv1"]


def build() -> tuple[dict[str, dict], list[dict]]:
    """Every parameter file (by name) and every case, in order."""
    files: dict[str, dict] = {}
    cases: list[dict] = []
    smoke_rel = "params/smoke/k0n31/M1-T01.json"
    smoke = v1_row(smoke_rel)
    smoke_recipe = recipe_of(smoke)
    rho_seed = 0x230000 + 1  # the suite's rho seed for target T01

    # Curve A: the smoke curve, under the repository's modulus for degree 31.
    f_a = curve_id.find_irreducible_sparse(31)
    assert f_a == 0x80000009 and irreducible(f_a)
    A = Curve(31, f_a, 0, 1)
    order_a = koblitz_order(0, 31)
    r_a = largest_prime_factor(order_a)
    h_a = order_a // r_a
    assert (order_a, r_a, h_a) == (2147574356, 1439393, 1492) and prime_exact(r_a) and h_a % r_a

    def run_case(cid: str, purpose: str, fname: str, expect_paths: dict, contains: dict | None = None,
                 command: str = "price", exit_code: int = 0, timeout: int = 300) -> dict:
        argv = [command, "--params", "{tmp}/params.json", "--json", "--out", "{tmp}/report.json"]
        if command == "price":
            argv += ["--repeats", "1", "--repeats-fast", "1"]
        expect = {"exit": exit_code, "json_file": "{tmp}/report.json", "json_paths": expect_paths}
        if contains:
            expect["json_contains"] = contains
        return {"id": cid, "step": "B1", "checks": purpose,
                "files": {"params.json": {"copy": f"{{cases}}/{fname}"}},
                "argv": argv, "env": {"RAYON_NUM_THREADS": "1"}, "expect": expect, "timeout_s": timeout}

    def refusal(cid: str, purpose: str, fname: str, code: str) -> dict:
        return run_case(cid, purpose, fname, {"refusal.code": code, "refusal.class": "invalid"},
                        exit_code=2, timeout=60)

    # C009: the v2 translation of the smoke row, by design §8, against the v1 file.
    doc = document("C009: suite v1 smoke row M1-T01, translated to v2", A, {"form": "koblitz", "a": 0},
                   r_a, h_a, {"rule": "koblitz_search_v1"}, {"public_hash_seed": smoke["targets"][0]["public_hash_seed"]},
                   smoke_recipe, rho_seed)
    files["C009-translation-of-smoke-row.json"] = doc
    outputs = ["status", "counts", "certificates.ic.scalar", "certificates.rho.scalar", "rho_counts",
               "all_verified", "ic_and_rho_agree"]
    case = run_case("C009-v2-translation-same-outputs", "design §8: the v2 translation of a v1 row runs to "
                    "identical outputs", "C009-translation-of-smoke-row.json", {"status": "complete"})
    case["files"]["v1.json"] = {"copy": "{suite}/" + smoke_rel}
    case["expect"]["same_outputs_as"] = {
        "argv": ["price", "--params", "{tmp}/v1.json", "--json", "--out", "{tmp}/v1-report.json",
                 "--single-target", "--rho-seed", str(rho_seed), "--repeats", "1", "--repeats-fast", "1"],
        "json_file": "{tmp}/v1-report.json", "paths": outputs}
    cases.append(case)

    # C010: curve A with another generator and a known answer.
    g_a = A.subgroup_point("ic-conformance-v2/curve-a/generator", h_a, r_a)
    k10 = int.from_bytes(hashlib.sha256(b"ic-conformance-v2/C010/known_log").digest(), "big") % (r_a - 1) + 1
    base = document("C010: curve A, another generator, a known answer", A, {"form": "koblitz", "a": 0},
                    r_a, h_a, g_a, {"known_log": str(k10)}, smoke_recipe, rho_seed)
    files["C010-curve-a-other-generator.json"] = base
    cases.append(run_case("C010-curve-a-other-generator", "B1's importer: an explicit generator under the "
                          "repository's modulus", "C010-curve-a-other-generator.json",
                          {"status": "complete", "result.scalar": str(k10), "result.verified": True,
                           "result.known_answer": True, "all_verified": True}))

    # C011: the same equation under another modulus, a different model.
    f_b = 0x80000041  # z^31 + z^6 + 1
    assert irreducible(f_b)
    B = Curve(31, f_b, 0, 1)
    g_b = B.subgroup_point("ic-conformance-v2/curve-a-prime/generator", h_a, r_a)
    k11 = int.from_bytes(hashlib.sha256(b"ic-conformance-v2/C011/known_log").digest(), "big") % (r_a - 1) + 1
    files["C011-curve-a-prime-other-modulus.json"] = document(
        "C011: y^2 + xy = x^3 + 1 under z^31 + z^6 + 1", B, {"form": "koblitz", "a": 0}, r_a, h_a, g_b,
        {"known_log": str(k11)}, smoke_recipe, rho_seed)
    assert icv1(B, order_a) != icv1(A, order_a)
    cases.append(run_case("C011-another-modulus-is-another-model", "B1's importer: a modulus the repository's "
                          "rule would not pick, named as its own model", "C011-curve-a-prime-other-modulus.json",
                          {"status": "complete", "result.scalar": str(k11), "result.verified": True,
                           "curve_id.icv1": icv1(B, order_a)}))

    # C012-C025: one defect each, on C010's valid document.
    def mutant(cid: str, fname: str, purpose: str, code: str, edit) -> None:
        d = copy.deepcopy(base)
        d["name"] = f"{cid.split('-', 1)[0]}: {purpose}"
        edit(d)
        files[fname] = d
        cases.append(refusal(cid, purpose, fname, code))

    off = (g_a[0], g_a[1] ^ 1)
    assert not A.on_curve(off)
    outside = A.point("ic-conformance-v2/curve-a/outside")
    assert A.on_curve(outside) and A.mul(r_a, outside) is not None
    target_outside = A.point("ic-conformance-v2/curve-a/target-outside")
    assert A.on_curve(target_outside) and A.mul(r_a, target_outside) is not None
    target_off = (target_outside[0], target_outside[1] ^ 1)
    assert not A.on_curve(target_off)
    assert not irreducible(0x80000001)
    assert not prime_exact(3 * r_a) and order_a != 1491 * r_a

    mutant("C012-unknown-key", "C012-unknown-key.json", "a key outside the schema, inside subgroup",
           "unknown-key", lambda d: d["subgroup"].update({"comment": "not a schema key"}))
    mutant("C013-schema-version", "C013-schema-version.json", "schema_version 3", "schema-version",
           lambda d: d.update({"schema_version": 3}))
    mutant("C014-two-targets", "C014-two-targets.json", "two target forms in one file", "target-count",
           lambda d: d["target"].update({"public_hash_seed": 23001}))
    mutant("C015-coordinates-without-modulus", "C015-coordinates-without-modulus.json",
           "a generator's coordinates, and no modulus", "modulus-required", lambda d: d["field"].pop("modulus"))
    mutant("C016-integer-syntax", "C016-integer-syntax.json", "an order written with a sign", "integer-syntax",
           lambda d: d["subgroup"].update({"order": "+" + str(r_a)}))
    mutant("C017-reducible-modulus", "C017-reducible-modulus.json", "the modulus z^31 + 1, which is reducible",
           "modulus-reducible", lambda d: d["field"].update({"modulus": hx(0x80000001)}))
    mutant("C018-singular-curve", "C018-singular-curve.json", "b = 0", "curve-singular",
           lambda d: d.update({"curve": {"form": "binary_weierstrass", "a": "0x0", "b": "0x0"}}))
    mutant("C019-generator-off-curve", "C019-generator-off-curve.json", "a generator off the curve",
           "generator-off-curve", lambda d: d["subgroup"].update({"generator": point_doc(off)}))
    mutant("C020-generator-order", "C020-generator-order.json",
           "a generator on the curve, outside the order-r subgroup", "generator-order",
           lambda d: d["subgroup"].update({"generator": point_doc(outside)}))
    mutant("C021-composite-order", "C021-composite-order.json", "order 3r, a composite", "order-composite",
           lambda d: d["subgroup"].update({"order": str(3 * r_a)}))
    mutant("C022-cofactor-mismatch", "C022-cofactor-mismatch.json", "cofactor 1491, where #E = 1492 r",
           "cofactor-mismatch", lambda d: d["subgroup"].update({"cofactor": "1491"}))
    mutant("C023-target-outside-subgroup", "C023-target-outside-subgroup.json",
           "a target on the curve, outside the order-r subgroup", "target-outside-subgroup",
           lambda d: d.update({"target": {"point": point_doc(target_outside)}}))
    mutant("C024-target-off-curve", "C024-target-off-curve.json", "a target off the curve", "target-off-curve",
           lambda d: d.update({"target": {"point": point_doc(target_off)}}))
    mutant("C025-known-log-range", "C025-known-log-range.json", "known_log = r", "known-log-range",
           lambda d: d.update({"target": {"known_log": str(r_a)}}))

    # C026: the identity target.
    d = copy.deepcopy(base)
    d["name"] = "C026: the identity as the target"
    d["target"] = {"point": "identity"}
    files["C026-identity-target.json"] = d
    cases.append(run_case("C026-identity-target", "design §4.4: the identity is a valid target with logarithm "
                          "0, by the trivial route", "C026-identity-target.json",
                          {"status": "complete", "result.scalar": "0", "result.verified": True,
                           "route.ic.pipeline": "trivial"}, timeout=60))

    # C027, C028: the m = 83 gate curve (AGENTS.md §8a) with a frozen generator and a public target.
    f_g = 0x800000000200000000007
    assert irreducible(f_g)
    G83 = Curve(83, f_g, 0, 1)
    order_g = koblitz_order(0, 83)
    r_g, h_g = 2417851639230796216685689, 4
    assert order_g == h_g * r_g and prime_exact(r_g)
    gen83 = G83.subgroup_point("ic-tool-programme/gate/m83/generator", h_g, r_g)
    tgt83 = G83.subgroup_point("ic-tool-programme/gate/m83/T001", h_g, r_g)
    files["gate-m83-T001.json"] = document(
        "the m = 83 confidence gate (AGENTS.md §8a): frozen generator, public target T001", G83,
        {"form": "koblitz", "a": 0}, r_g, h_g, gen83, {"point": point_doc(tgt83)}, "auto", rho_seed)
    c27 = run_case("C027-gate-curve-has-no-route-yet", "design §5: the gate curve is valid and no pipeline "
                   "admits it before B3", "gate-m83-T001.json",
                   {"refusal.code": "no-ic-route", "refusal.class": "unsupported"},
                   {"route.considered": [{"pipeline": "kic", "admitted": False, "gate": "field-wider-than-one-word"},
                                         {"pipeline": "rho-koblitz", "admitted": False,
                                          "gate": "field-wider-than-one-word"}]},
                   exit_code=3, timeout=120)
    c27["until"] = "B3"
    cases.append(c27)
    cases.append(run_case("C028-gate-curve-checks", "design §4: every check passes on the gate file, r's "
                          "primality exactly", "gate-m83-T001.json", {"status": "checks_passed"},
                          {"checks": [{"code": "order-composite", "status": "pass", "exact": True},
                                      {"code": "cofactor-mismatch", "status": "pass", "exact": True},
                                      {"code": "generator-order", "status": "pass"},
                                      {"code": "target-outside-subgroup", "status": "pass"}],
                           "disclosures": [{"code": "frobenius-module", "value": {"ord_n_2": 82, "blocks": 1}}]},
                          command="check", timeout=120))

    # C029: a composite degree, the suite's curve at (a, n) = (1, 45).
    s45 = v1_row("params/S/k1n45/M1-T01.json")
    f_45 = curve_id.find_irreducible_sparse(45)
    assert f_45 == 0x20000000001B and irreducible(f_45)
    K45 = Curve(45, f_45, 1, 1)
    order_45 = koblitz_order(1, 45)
    r_45 = largest_prime_factor(order_45)
    h_45 = order_45 // r_45
    assert r_45 == 29264761 and prime_exact(r_45) and h_45 % r_45
    g_45 = K45.subgroup_point("ic-conformance-v2/k45/generator", h_45, r_45)
    k29 = int.from_bytes(hashlib.sha256(b"ic-conformance-v2/C029/known_log").digest(), "big") % (r_45 - 1) + 1
    files["C029-composite-degree-45.json"] = document(
        "C029: the suite's degree-45 curve, explicit generator, known answer", K45, {"form": "koblitz", "a": 1},
        r_45, h_45, g_45, {"known_log": str(k29)}, recipe_of(s45), rho_seed)
    cases.append(run_case("C029-composite-degree-disclosed", "AGENTS.md §8b: a composite degree runs, with its "
                          "intermediate subfields disclosed", "C029-composite-degree-45.json",
                          {"status": "complete", "result.scalar": str(k29), "result.verified": True},
                          {"disclosures": [{"code": "intermediate-subfields", "value": [3, 5, 9, 15]}]}))

    # C030: an even degree, which kic cannot take.
    f_32 = curve_id.find_irreducible_sparse(32)
    assert irreducible(f_32)
    K32 = Curve(32, f_32, 0, 1)
    order_32 = koblitz_order(0, 32)
    r_32 = largest_prime_factor(order_32)
    h_32 = order_32 // r_32
    # An even degree splits #E through GF(2^16), so r is smaller than h here; v2 allows that.
    assert prime_exact(r_32) and h_32 % r_32 and (r_32, h_32) == (32993, 130176)
    g_32 = K32.subgroup_point("ic-conformance-v2/k32/generator", h_32, r_32)
    k30 = int.from_bytes(hashlib.sha256(b"ic-conformance-v2/C030/known_log").digest(), "big") % (r_32 - 1) + 1
    files["C030-even-degree-32.json"] = document(
        "C030: y^2 + xy = x^3 + 1 over GF(2^32)", K32, {"form": "koblitz", "a": 0}, r_32, h_32, g_32,
        {"known_log": str(k30)}, smoke_recipe, rho_seed)
    cases.append(run_case("C030-even-degree-has-no-ic-route", "design §5.2: kic needs an odd extension degree",
                          "C030-even-degree-32.json", {"refusal.code": "no-ic-route", "refusal.class": "unsupported"},
                          {"route.considered": [{"pipeline": "kic", "admitted": False,
                                                 "gate": "even-extension-degree"}]}, exit_code=3, timeout=120))

    # C031: ECC2K-130, with the challenge's polynomial-basis points from ic fixed's profile.
    fixed = json.loads((ROOT / "docs" / "ic" / "params" / "ecc2k130-fixed.json").read_text())
    f_c = int(fixed["curve"]["polynomial"], 16)
    assert f_c == 0x800000000000000000000000000002007 and irreducible(f_c)
    C130 = Curve(131, f_c, 0, 1)
    r_c, h_c = int(fixed["curve"]["subgroup_order"]), int(fixed["curve"]["cofactor"])
    assert koblitz_order(0, 131) == h_c * r_c and r_c > PSI13 and miller_rabin(r_c, BASES13)
    g_c = (int(fixed["generator"]["x"], 16), int(fixed["generator"]["y"], 16))
    t_c = next(t for t in fixed["targets"] if t["id"] == "ecc2k-130")
    q_c = (int(t_c["x"], 16), int(t_c["y"], 16))
    for p in (g_c, q_c):
        assert C130.on_curve(p) and C130.mul(r_c, p) is None
    files["ecc2k130-challenge.json"] = document(
        "ECC2K-130: the challenge's generator and target, in the polynomial basis "
        "(docs/ic/params/ecc2k130-fixed.json)", C130, {"form": "koblitz", "a": 0}, r_c, h_c, g_c,
        {"point": point_doc(q_c)}, "auto", rho_seed)
    c31 = run_case("C031-challenge-checks", "design §4: the challenge file is valid; r's primality is a screen",
                   "ecc2k130-challenge.json", {"status": "checks_passed"},
                   {"checks": [{"code": "order-composite", "status": "pass", "exact": False}],
                    "disclosures": [{"code": "primality-screen"}],
                    "route.considered": [{"pipeline": "kic", "admitted": False, "gate": "field-wider-than-one-word"}]},
                   command="check", timeout=120)
    c31["until"] = "B4"
    cases.append(c31)
    return files, cases


def selftest() -> None:
    """The arithmetic against counting: on small fields of both parities, the points enumerated equal
    the trace recurrence's #E, [#E]P = O for every point tried, and the group law is associative on them."""
    for n in (7, 8, 9, 10, 11):
        f = curve_id.find_irreducible_sparse(n)
        assert irreducible(f)
        for a in (0, 1):
            c = Curve(n, f, a, 1)
            points = [p for x in range(1 << n) for p in ([c.lift(x)] if c.lift(x) else [])]
            # Each abscissa x != 0 with a root gives two points; x = 0 gives one.
            count = 1 + sum(1 if p[0] == 0 else 2 for p in points)
            assert count == koblitz_order(a, n), (n, a)
            for p in points[:8]:
                assert c.mul(count, p) is None
                q, s = points[-1], points[len(points) // 2]
                assert c.add(c.add(p, q), s) == c.add(p, c.add(q, s))


def texts() -> dict[str, str]:
    selftest()
    files, cases = build()
    out = {f"params/{name}": json.dumps(doc, indent=1) + "\n" for name, doc in files.items()}
    out["cases.json"] = json.dumps({
        "suite": "ic tool programme conformance suite v2",
        "design": "research/ic_tool_program/design/schema-v2.md §9",
        "includes": "../v1/cases.json, run first and unchanged; its cases are step B0's",
        "rules": [
            "A case passes only if the process exits within its timeout, never panics (exit status 101), "
            "and meets every expectation.",
            "{tmp} is the case's own temporary directory; {suite} is research/ic_tool_program/suite/v1; "
            "{cases} is this directory's params/.",
            "exit is a number, or zero / nonzero. json_paths compares values at dotted paths. "
            "json_contains asks, for each listed object, that some element of the list at the path "
            "contains it. same_outputs_as runs a second command and compares the listed paths.",
            "Cases are added, never removed or loosened. A case with `until` expects a refusal that step "
            "lifts; that step changes the expectation with a dated note here, and adds the new run as its "
            "own case.",
        ],
        "generator": {"path": "research/ic_tool_program/conformance/v2/make_cases.py",
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
        raise SystemExit("cases.json exists; the B1 cases are frozen (use --check)")
    for rel, t in out.items():
        path = HERE / rel
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, "x") as fh:
            fh.write(t)
    print(f"{len(out)} files written")


if __name__ == "__main__":
    main()
