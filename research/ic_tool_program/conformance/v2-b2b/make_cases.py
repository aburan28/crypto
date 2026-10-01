#!/usr/bin/env python3
"""Conformance cases for B2b: `kic` on curves over a subfield, `kic`
alone, and order certificates (C059-C070;
`research/ic_tool_program/rounds/B2b-subfield-kic-certificates/PROTOCOL.md`).

    python3 make_cases.py           # writes params/, cases.json and SHA256SUMS; refuses to overwrite
    python3 make_cases.py --check   # re-derives every file and compares the bytes

Written at B2b's declaration, before any B2b code.  `../run.py` runs
the cases with `--steps ...,B2b`.

The arithmetic is B1's for GF(2^n) and its curves (`../v2/make_cases.py`),
B2's for prime curves and for group orders from the orders of points
(`../v2-b2/make_cases.py`), and the following, written here:
- GF(2^k) inside GF(2^n), as the kernel of x -> x^(2^k) + x;
- #E(GF(2^k)) from the trace over GF(2^k), and #E(GF(2^n)) by the trace
  recurrence.  Every curve's order is found a second way, by the orders
  of its points;
- Pocklington certificates: building one from a factorisation, and
  checking one by the protocol's rules.

None of it is shared with the Rust tool.  Every file's intended property
is asserted when it is built.
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


def _load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


v2 = _load("conformance_v2_cases", HERE.parent / "v2" / "make_cases.py")
b2 = _load("conformance_v2_b2_cases", HERE.parent / "v2-b2" / "make_cases.py")

STEP = "B2b"
RHO_SEED = 0x230000 + 1
PSI13 = v2.PSI13
CERT_MAX_DEPTH = 64
CERT_MAX_FACTORS = 1024


# -- GF(2^k) inside GF(2^n) ---------------------------------------------------------------------

def frobenius(F, x: int, k: int) -> int:
    for _ in range(k):
        x = F.sqr(x)
    return x


def subfield(F, k: int) -> list[int]:
    """Every element of GF(2^k) inside F = GF(2^n), ascending: the kernel of x -> x^(2^k) + x."""
    pivots, kernel = [], []  # pivots: (image, combination), leading bits distinct, descending
    for i in range(F.n):
        v, c = frobenius(F, 1 << i, k) ^ (1 << i), 1 << i
        for pv, pc in pivots:
            if v ^ pv < v:
                v, c = v ^ pv, c ^ pc
        if v:
            pivots.append((v, c))
            pivots.sort(reverse=True)
        else:
            kernel.append(c)
    assert len(kernel) == k, "GF(2^k) has dimension k"
    out = set()
    for m in range(1 << k):
        x = 0
        for i, e in enumerate(kernel):
            if m >> i & 1:
                x ^= e
        out.add(x)
    out = sorted(out)
    assert len(out) == 1 << k and all(frobenius(F, x, k) == x for x in out)
    return out


def least_subfield(F, a: int, b: int) -> int:
    """The least d with a, b in GF(2^d)."""
    n = F.n
    return next(d for d in range(1, n + 1)
                if n % d == 0 and frobenius(F, a, d) == a and frobenius(F, b, d) == b)


def trace_over(F, c: int, k: int) -> int:
    """The trace from GF(2^k) to GF(2) of c in GF(2^k)."""
    t, s = 0, c
    for _ in range(k):
        t ^= s
        s = F.sqr(s)
    assert t in (0, 1)
    return t


def count_over_subfield(curve, k: int, elements: list[int]) -> int:
    """#E(GF(2^k)) for y^2 + xy = x^3 + ax^2 + b, a and b in GF(2^k): O, (0, sqrt b), and two points
    for each x != 0 for which z^2 + z = x + a + b/x^2 has a root in GF(2^k), that is, whose trace over
    GF(2^k) is 0."""
    F = curve.k
    count = 2
    for x in elements:
        if x:
            c = x ^ curve.a ^ F.mul(curve.b, F.inv(F.sqr(x)))
            count += 0 if trace_over(F, c, k) else 2
    return count


def order_over_extension(t: int, q: int, e: int) -> int:
    """#E(GF(q^e)) = q^e + 1 - s_e, with s_0 = 2, s_1 = t, s_i = t s_{i-1} - q s_{i-2}."""
    prev, cur = 2, t
    for _ in range(1, e):
        prev, cur = cur, t * cur - q * prev
    return q ** e + 1 - cur


def subfield_curve(n: int, k: int, a: int, b: int, label: str):
    """y^2 + xy = x^3 + ax^2 + b over GF(2^n) under the repository's modulus, with a and b in GF(2^k):
    its #E by the subfield count and again by points' orders, and its largest prime factor r."""
    f = v2.curve_id.find_irreducible_sparse(n)
    assert v2.irreducible(f)
    C = v2.Curve(n, f, a, b)
    assert least_subfield(C.k, a, b) == k
    elements = subfield(C.k, k)
    t = (1 << k) + 1 - count_over_subfield(C, k, elements)
    N = order_over_extension(t, 1 << k, n // k)
    by_points = b2.group_order(b2.BinaryGroup(C), 1 << n, (C.point(f"{label}/order/{i}") for i in range(64)))
    assert by_points == N, (by_points, N)
    r = v2.largest_prime_factor(N)
    h = N // r
    assert v2.prime_exact(r) and h % r, "r is prime and r^2 does not divide #E"
    return C, N, r, h


def first_subfield_curve(n: int, k: int, accept, a_values=None):
    """The first (a, b), ascending over GF(2^k), with least subfield k whose #E's largest prime factor r
    and cofactor h pass `accept(r, h)`."""
    f = v2.curve_id.find_irreducible_sparse(n)
    C0 = v2.Curve(n, f, 0, 1)
    elements = subfield(C0.k, k)
    for a in (a_values if a_values is not None else elements):
        for b in elements:
            if b == 0 or least_subfield(C0.k, a, b) != k:
                continue
            C = v2.Curve(n, f, a, b)
            N = order_over_extension((1 << k) + 1 - count_over_subfield(C, k, elements), 1 << k, n // k)
            r = v2.largest_prime_factor(N)
            if v2.prime_exact(r) and (N // r) % r and accept(r, N // r):
                return a, b
    raise AssertionError("no curve")


def recipe(e: int, r: int, columns: int, m: int = 2, seed: int = 201) -> dict:
    """Ledger §20's rules (research/ic_exponent_20260926/make_params.py) with the q-power Frobenius
    orbit's length e in the place of n: 2e·columns points, window max(1, round(|F|/32)), unit
    clamp(round(2r/(e|F|^2)), 16, 65536), descent cap max(10^5, 64·ceil(2.02r/|F|^2))."""
    points = 2 * e * columns
    return {
        "summands": 3,
        "descent_summands": m,
        "collection_window": max(1, round(points / 32)),
        "collection_aim": True,
        "solver": "pair_table",
        "seed": seed,
        "max_trials": min(max(100_000, 64 * math.ceil(2.02 * r / (points * points))), 2**32 - 1),
        "collection": {"unit_trials": min(65536, max(16, round(2 * r / (e * points * points)))),
                       "units": 1, "max_units": 100_000},
        "factor_base": {"mode": "spec", "spec": {"kind": "subgroup_orbits", "points": points, "seed": seed}},
    }


# -- Pocklington certificates ---------------------------------------------------------------------

def parse_integer(s) -> int | None:
    """The schema's integer syntax: a string, decimal or 0x hexadecimal, unsigned, at most 1,024 bits."""
    if not isinstance(s, str) or not s:
        return None
    if s.startswith("0x"):
        digits, base = s[2:], 16
        ok = digits != "" and all(ch in "0123456789abcdefABCDEF" for ch in digits)
    else:
        digits, base = s, 10
        ok = all(ch in "0123456789" for ch in digits)
    if not ok:
        return None
    v = int(digits, base)
    return v if v.bit_length() <= 1024 else None


def check_certificate(N: int, cert, depth: int = 0, budget: list[int] | None = None) -> str | None:
    """None when `cert` proves N prime by the protocol's rules, else why it does not.

    A certificate is {"kind": "pocklington", "factors": [...]}; each factor is {"prime", "exponent",
    "witness"} and, optionally, "certificate", a certificate for that prime.  With F the product of
    prime^exponent: F divides N - 1, (F + 1)^2 > N, and for each factor q with witness a,
    1 < a < N, a^(N-1) = 1 (mod N) and gcd(a^((N-1)/q) - 1, N) = 1.  Every q is prime: below psi_13
    by Miller-Rabin on the first 13 primes, otherwise by its own certificate.  Then every prime factor
    of N is 1 mod F, so above sqrt(N), and N is prime (Pocklington; Brillhart, Lehmer and Selfridge
    1975)."""
    if budget is None:
        budget = [CERT_MAX_FACTORS]
    if depth >= CERT_MAX_DEPTH:
        return f"deeper than {CERT_MAX_DEPTH} levels"
    if not isinstance(cert, dict) or set(cert) != {"kind", "factors"} or cert["kind"] != "pocklington":
        return "not {\"kind\": \"pocklington\", \"factors\": [...]}"
    factors = cert["factors"]
    if not isinstance(factors, list) or not factors:
        return "factors is not a non-empty list"
    if N < 3:
        return "N is below 3"
    F, seen, parsed = 1, set(), []
    for i, entry in enumerate(factors):
        budget[0] -= 1
        if budget[0] < 0:
            return f"more than {CERT_MAX_FACTORS} factors in all"
        if not isinstance(entry, dict) or not ({"prime", "exponent", "witness"} <= set(entry)
                                               <= {"prime", "exponent", "witness", "certificate"}):
            return f"factor {i}: keys"
        q, a, e = parse_integer(entry["prime"]), parse_integer(entry["witness"]), entry["exponent"]
        if q is None or a is None or isinstance(e, bool) or not isinstance(e, int) or not 1 <= e <= 1024:
            return f"factor {i}: values"
        if q < 2 or q in seen:
            return f"factor {i}: {q} is below 2 or repeated"
        seen.add(q)
        F *= q ** e
        if F > N:
            return "F exceeds N - 1"
        parsed.append((q, a, entry))
    if (N - 1) % F:
        return "F does not divide N - 1"
    if (F + 1) ** 2 <= N:
        return "(F + 1)^2 <= N"
    for q, a, entry in parsed:
        if not 1 < a < N:
            return f"witness {a} for {q} is not in (1, N)"
        if pow(a, N - 1, N) != 1:
            return f"witness {a} for {q}: a^(N-1) != 1"
        if math.gcd(pow(a, (N - 1) // q, N) - 1, N) != 1:
            return f"witness {a} for {q}: gcd(a^((N-1)/q) - 1, N) != 1"
        if "certificate" in entry:
            why = check_certificate(q, entry["certificate"], depth + 1, budget)
            if why:
                return f"factor {q}: {why}"
        elif not (q < PSI13 and v2.prime_exact(q)):
            return f"factor {q} is not below psi_13 and prime, and has no certificate"
    return None


def certificate(N: int, factors: list[tuple[int, int]], nested: dict[int, dict] | None = None) -> dict:
    """A certificate for N from prime powers (q, e): each q with the least witness a >= 2, and each q at
    or above psi_13 with its certificate from `nested`."""
    out = []
    for q, e in factors:
        a = 2
        while not (pow(a, N - 1, N) == 1 and math.gcd(pow(a, (N - 1) // q, N) - 1, N) == 1):
            a += 1
            assert a < 1 << 16, "no small witness"
        entry = {"prime": str(q), "exponent": e, "witness": str(a)}
        if q >= PSI13:
            entry["certificate"] = (nested or {})[q]
        out.append(entry)
    cert = {"kind": "pocklington", "factors": out}
    why = check_certificate(N, cert)
    assert why is None, why
    return cert


def probable_prime(n: int, label: str) -> bool:
    """Miller-Rabin on the first 13 primes and 24 labelled bases: a screen, used only where a
    certificate or the tool's own test decides."""
    extra = [2 + b2.shake(f"{label}/base/{i}", 32) % (n - 3) for i in range(24)] if n > 5 else []
    return v2.miller_rabin(n, v2.BASES13 + tuple(extra))


def factorise(m: int, label: str) -> list[tuple[int, int]]:
    """m's prime factorisation, ascending, by trial division to 2^20 and Pollard-Brent; every factor
    below psi_13."""
    small, rest = b2.trial_factor(m, 1 << 20)
    out = list(small)
    pending = [rest] if rest > 1 else []
    while pending:
        x = pending.pop()
        if x < PSI13 and v2.prime_exact(x):
            out.append(x)
            continue
        g = b2.brent_factor(x, f"{label}/{x}", 1 << 22)
        assert g, f"no factor of {x}"
        pending += [g, x // g]
    return sorted((q, out.count(q)) for q in set(out))


def selftest() -> None:
    """The subfield count against every pair (x, y), and the certificate checker on primes,
    composites and broken certificates."""
    checked = 0
    for k, n in ((2, 10), (2, 14), (3, 9), (3, 15), (4, 12)):
        F = v2.Field(n, v2.curve_id.find_irreducible_sparse(n))
        E = subfield(F, k)
        for a in E[:3]:
            for b in E[1:4]:
                C = v2.Curve(n, F.f, a, b)
                brute = 1 + sum(1 for x in E for y in E if C.on_curve((x, y)))
                assert count_over_subfield(C, k, E) == brute, (k, n, a, b)
                N = order_over_extension((1 << k) + 1 - brute, 1 << k, n // k)
                try:
                    found = b2.group_order(b2.BinaryGroup(C), 1 << n,
                                           (C.point(f"selftest/{n}/{a}/{b}/{i}") for i in range(64)))
                except ValueError:
                    continue
                assert found == N, (k, n, a, b)
                checked += 1
    assert checked >= 20, checked
    # Certificates.
    for N in (3, 5, 7, 97, 1009, 65537, 2**61 - 1, 1439393):
        cert = certificate(N, factorise(N - 1, f"selftest/{N}"))
        assert check_certificate(N, cert) is None
    p61 = 2**61 - 1
    good = certificate(p61, factorise(p61 - 1, "selftest/p61"))
    for edit in (lambda c: c.update({"factors": c["factors"][:1]}),                # F = 2, too small
                 lambda c: c["factors"][0].update({"witness": "1"}),               # witness range
                 lambda c: c["factors"][0].update({"exponent": 0}),                # exponent range
                 lambda c: c["factors"][0].update({"prime": "4"}),                 # not prime, not dividing
                 lambda c: c.update({"kind": "pratt"}),
                 lambda c: c["factors"].append(copy.deepcopy(c["factors"][0]))):   # repeated
        bad = copy.deepcopy(good)
        edit(bad)
        assert check_certificate(p61, bad) is not None
    for composite in (561, 1105, 2**61 + 1, 3215031751):  # Carmichael numbers and a strong pseudoprime
        cert = {"kind": "pocklington", "factors": [{"prime": str(q), "exponent": e, "witness": "2"}
                                                   for q, e in factorise(composite - 1, "selftest/c")]}
        assert check_certificate(composite, cert) is not None


# -- the cases ------------------------------------------------------------------------------------

def build() -> tuple[dict[str, dict], list[dict]]:
    files: dict[str, dict] = {}
    cases: list[dict] = []

    def case(cid: str, purpose: str, params: dict, expect: dict, command: str = "price",
             argv_tail: list[str] | None = None, extra_files: dict | None = None,
             same_outputs_as: dict | None = None, timeout: int = 300) -> None:
        argv = [command, "--params", "{tmp}/params.json", "--json", "--out", "{tmp}/report.json"]
        if argv_tail:
            argv += argv_tail
        expect = {"json_file": "{tmp}/report.json", **expect}
        if same_outputs_as:
            expect["same_outputs_as"] = same_outputs_as
        cases.append({"id": cid, "step": STEP, "checks": purpose,
                      "files": {"params.json": params, **(extra_files or {})}, "argv": argv,
                      "env": {"RAYON_NUM_THREADS": "1"}, "expect": expect, "timeout_s": timeout})

    def refused(code: str, cls: str, status: int) -> dict:
        return {"exit": status, "json_paths": {"refusal.code": code, "refusal.class": cls}}

    one_rep = ["--repeats", "1", "--repeats-fast", "1"]

    def ic_alone_matches(paired_file: str) -> dict:
        return {"argv": ["price", "--params", "{tmp}/paired.json", "--json", "--out", "{tmp}/paired-report.json",
                         *one_rep],
                "json_file": "{tmp}/paired-report.json",
                "paths": ["counts", "certificates.ic.scalar", "certificates.ic.target", "certificates.ic.curve"]}

    def binary_doc(name: str, C, r: int, h: int, G, target: dict, method: dict) -> dict:
        return {"schema_version": 2, "name": name,
                "field": {"kind": "binary", "degree": C.n, "modulus": v2.hx(C.f)},
                "curve": {"form": "binary_weierstrass", "a": v2.hx(C.a), "b": v2.hx(C.b)},
                "subgroup": {"order": str(r), "cofactor": str(h), "generator": v2.point_doc(G)},
                "target": target, "method": method}

    def paired(ic_pipeline: str, ic_recipe, rho_pipeline: str = "auto") -> dict:
        return {"solve": "paired", "index_calculus": {"pipeline": ic_pipeline, "recipe": ic_recipe},
                "rho": {"pipeline": rho_pipeline, "seed": RHO_SEED}}

    # Curve B: y^2 + xy = x^3 + w x^2 + 1 over GF(2^22), w^2 + w + 1 = 0: defined over GF(4), e = 11.
    F22 = v2.Field(22, v2.curve_id.find_irreducible_sparse(22))
    w22 = F22.solve_quadratic(1)
    assert w22 is not None and F22.sqr(w22) ^ w22 == 1
    CB, NB, rB, hB = subfield_curve(22, 2, w22, 1, "ic-conformance-b2b/curve-b")
    assert (rB, hB) == (2097349, 2) and rB > hB
    GB = CB.subgroup_point("ic-conformance-b2b/curve-b/generator", hB, rB)
    kB = b2.known_log("ic-conformance-b2b/curve-b/known_log", rB)
    recipe_b = recipe(11, rB, 24)
    doc59 = binary_doc("C059: a curve over GF(4) taken over GF(2^22), kic paired with rho", CB, rB, hB, GB,
                       {"known_log": str(kB)}, paired("kic", recipe_b))
    files["C059-subfield-curve-paired.json"] = doc59
    case("C059-subfield-curve-paired",
         "B2b: kic on a curve over GF(4) (k = 2, e = 11), paired with rho-negation on the same point; the "
         "report says rho leaves the q-power Frobenius unused and that the pair is not speedup-eligible",
         {"copy": "{here}/C059-subfield-curve-paired.json"},
         {"exit": 0,
          "json_paths": {"status": "complete", "result.verified": True, "result.known_answer": True,
                         "result.scalar": str(kB), "all_verified": True, "route.ic.pipeline": "kic",
                         "route.rho.pipeline": "rho-negation", "speedup_eligible": False},
          "json_contains": {"disclosures": [{"code": "curve-subfield", "value": {"degree": 2}},
                                            {"code": "rho-frobenius-unused"}],
                            "route.considered": [{"pipeline": "kic", "admitted": True},
                                                 {"pipeline": "rho-koblitz", "admitted": False,
                                                  "gate": "subfield-curve-unsupported"}]}},
         argv_tail=one_rep)

    case("C060-subfield-curve-no-recipe",
         "design §5.3: recipe auto has no rule for kic on a curve over a subfield with k > 1",
         {"copy": "{here}/C059-subfield-curve-paired.json", "set": {"method.index_calculus.recipe": "auto"}},
         refused("no-recipe", "unsupported", 3), argv_tail=one_rep)

    case("C061-subfield-curve-kic-alone",
         "B2b: solve index_calculus runs kic alone; its counts and certificate are the paired run's",
         {"copy": "{here}/C059-subfield-curve-paired.json",
          "set": {"method": {"solve": "index_calculus",
                             "index_calculus": {"pipeline": "kic", "recipe": recipe_b}}}},
         {"exit": 0,
          "json_paths": {"status": "complete", "operation": "ic_single_target", "ic.pipeline": "kic",
                         "result.scalar": str(kB), "result.verified": True, "result.known_answer": True,
                         "speedup_eligible": False}},
         argv_tail=one_rep, extra_files={"paired.json": {"copy": "{here}/C059-subfield-curve-paired.json"}},
         same_outputs_as=ic_alone_matches("paired.json"))

    # C062: kic alone on curve A (k = 1), against B1's paired C010.
    c10 = json.loads((HERE.parent / "v2" / "params" / "C010-curve-a-other-generator.json").read_text())
    case("C062-curve-a-kic-alone",
         "B2b: kic alone on a Koblitz curve gives the paired price's index calculus, count for count",
         {"copy": "{cases}/C010-curve-a-other-generator.json",
          "set": {"method": {"solve": "index_calculus", "index_calculus": c10["method"]["index_calculus"]}}},
         {"exit": 0,
          "json_paths": {"status": "complete", "operation": "ic_single_target", "ic.pipeline": "kic",
                         "result.scalar": c10["target"]["known_log"], "result.verified": True,
                         "speedup_eligible": False}},
         argv_tail=one_rep, extra_files={"paired.json": {"copy": "{cases}/C010-curve-a-other-generator.json"}},
         same_outputs_as=ic_alone_matches("paired.json"))

    # Curve C: over GF(8), taken over GF(2^33) (e = 11), r above 2^29.
    aC, bC = first_subfield_curve(33, 3, lambda r, h: r > h and r >= 1 << 29)
    CC, NC, rC, hC = subfield_curve(33, 3, aC, bC, "ic-conformance-b2b/curve-c")
    assert (rC, hC) == (715829951, 12)
    GC = CC.subgroup_point("ic-conformance-b2b/curve-c/generator", hC, rC)
    kC = b2.known_log("ic-conformance-b2b/curve-c/known_log", rC)
    files["C063-subfield-curve-routed.json"] = binary_doc(
        "C063: a curve over GF(8) taken over GF(2^33), both pipelines routed", CC, rC, hC, GC,
        {"known_log": str(kC)}, paired("auto", recipe(11, rC, 48)))
    case("C063-subfield-curve-routed",
         "design §5.2: with a recipe, auto routes the index calculus to kic on a curve over GF(8) (k = 3, "
         "n = 33 odd), and rho to rho-negation",
         {"copy": "{here}/C063-subfield-curve-routed.json"},
         {"exit": 0,
          "json_paths": {"status": "complete", "result.verified": True, "result.scalar": str(kC),
                         "all_verified": True, "route.ic.pipeline": "kic", "route.rho.pipeline": "rho-negation",
                         "speedup_eligible": False},
          "json_contains": {"disclosures": [{"code": "curve-subfield", "value": {"degree": 3}},
                                            {"code": "rho-frobenius-unused"}],
                            "route.considered": [{"pipeline": "ic-binary-s4", "admitted": False,
                                                  "gate": "recipe-not-taken"},
                                                 {"pipeline": "rho-koblitz", "admitted": False,
                                                  "gate": "subfield-curve-unsupported"}]}},
         argv_tail=one_rep)

    # Curve D: over GF(2^9), taken over GF(2^27) (e = 3): k = 9 > 8.
    aD, bD = first_subfield_curve(27, 9, lambda r, h: r > h, a_values=[0])
    CD, ND, rD, hD = subfield_curve(27, 9, aD, bD, "ic-conformance-b2b/curve-d")
    assert (aD, rD, hD) == (0, 81049, 1656)
    GD = CD.subgroup_point("ic-conformance-b2b/curve-d/generator", hD, rD)
    files["C064-subfield-too-large.json"] = binary_doc(
        "C064: a curve over GF(2^9) taken over GF(2^27), kic named", CD, rD, hD, GD,
        {"known_log": str(b2.known_log("ic-conformance-b2b/curve-d/known_log", rD))},
        paired("kic", recipe(3, rD, 8)))
    case("C064-subfield-too-large", "design §5.2: kic takes a curve over GF(2^k) for k <= 8 only",
         {"copy": "{here}/C064-subfield-too-large.json"}, refused("subfield-too-large", "unsupported", 3),
         argv_tail=one_rep)

    # Curve E: over GF(4), taken over GF(2^20) (e = 10, even).  For even e, #E(GF(q^e)) is
    # #E(GF(q^(e/2))) times its twist's count, so r < h, which validation allows.
    F20 = v2.Field(20, v2.curve_id.find_irreducible_sparse(20))
    w20 = F20.solve_quadratic(1)
    CE, NE, rE, hE = subfield_curve(20, 2, w20, 1, "ic-conformance-b2b/curve-e")
    assert (rE, hE) == (541, 1936)
    GE = CE.subgroup_point("ic-conformance-b2b/curve-e/generator", hE, rE)
    files["C065-subfield-even-extension.json"] = binary_doc(
        "C065: a curve over GF(4) taken over GF(2^20), kic named", CE, rE, hE, GE,
        {"known_log": str(b2.known_log("ic-conformance-b2b/curve-e/known_log", rE))},
        paired("kic", recipe(10, rE, 8)))
    case("C065-subfield-even-extension", "design §5.2: kic needs e = n/k odd",
         {"copy": "{here}/C065-subfield-even-extension.json"},
         refused("even-extension-degree", "unsupported", 3), argv_tail=one_rep)

    # C066-C067: ECC2K-130's r, certified.
    challenge = json.loads((HERE.parent / "v2" / "params" / "ecc2k130-challenge.json").read_text())
    r130 = int(challenge["subgroup"]["order"])
    assert r130 == v2.koblitz_order(0, 131) // 4 and r130 > PSI13
    factors130 = factorise(r130 - 1, "ic-conformance-b2b/ecc2k130")
    assert all(q < PSI13 for q, _ in factors130)
    cert130 = certificate(r130, factors130)
    case("C066-challenge-order-certified",
         "design §3.3: a Pocklington certificate makes r's primality exact at 2^129",
         {"copy": "{cases}/ecc2k130-challenge.json", "set": {"subgroup.order_certificate": cert130}},
         {"exit": 0, "json_paths": {"status": "checks_passed"},
          "json_contains": {"checks": [{"code": "certificate-invalid", "status": "pass", "exact": True},
                                       {"code": "order-composite", "status": "pass", "exact": True}]}},
         command="check", timeout=120)
    short = copy.deepcopy(cert130)
    short["factors"] = [f for f in short["factors"] if int(f["prime"]) != max(q for q, _ in factors130)]
    assert check_certificate(r130, short) == "(F + 1)^2 <= N"
    case("C067-challenge-certificate-invalid",
         "design §4.4: a certificate whose factored part is below sqrt(r) does not verify",
         {"copy": "{cases}/ecc2k130-challenge.json", "set": {"subgroup.order_certificate": short}},
         refused("certificate-invalid", "invalid", 2), command="check", timeout=120)

    # C068-C069: a prime field whose p is certified two levels deep.  q1 < psi_13; q2 = 2 k q1 + 1
    # near 2^150; p = 2 k' q2 + 1 near 2^200 with k' odd, so p = 3 mod 4, and r = (p + 1)/4 prime:
    # y^2 = x^3 + x is supersingular there, #E = p + 1 = 4r.
    label = "ic-conformance-b2b/C068"
    q1 = b2.next_prime((1 << 79) + b2.shake(f"{label}/q1", 16) % (1 << 79))
    assert q1 < PSI13
    k2 = (1 << 69) + b2.shake(f"{label}/q2", 16) % (1 << 68)
    while not probable_prime(2 * k2 * q1 + 1, f"{label}/q2"):
        k2 += 1
    q2 = 2 * k2 * q1 + 1
    cert_q2 = certificate(q2, [(2, 1), (q1, 1)])
    kp = (1 << 49) + b2.shake(f"{label}/p", 16) % (1 << 48) | 1
    while not (probable_prime(2 * kp * q2 + 1, f"{label}/p") and probable_prime((kp * q2 + 1) // 2, f"{label}/r")):
        kp += 2
    p68 = 2 * kp * q2 + 1
    r68 = (p68 + 1) // 4
    assert p68 % 4 == 3 and p68 > PSI13 and q2 > PSI13 and r68 > 4 * math.isqrt(p68)
    cert_p = certificate(p68, [(2, 1), (q2, 1)], {q2: cert_q2})
    E68 = b2.PrimeCurve(p68, 1, 0)
    for i in range(4):
        assert E68.mul(p68 + 1, E68.point(f"{label}/check/{i}")) is None
    G68 = E68.subgroup_point(f"{label}/generator", 4, r68)
    files["C068-prime-field-certified.json"] = b2.prime_doc(
        "C068: y^2 = x^3 + x over a 200-bit p certified two levels deep", p68,
        {"form": "short_weierstrass", "a": "1", "b": "0"}, r68, 4, b2.dec_point(G68),
        {"known_log": str(b2.known_log(f"{label}/known_log", r68))}, b2.rho_alone())
    files["C068-prime-field-certified.json"]["field"]["p_certificate"] = cert_p
    case("C068-prime-field-certified",
         "design §4.2: a nested Pocklington certificate makes p's primality exact at 2^200",
         {"copy": "{here}/C068-prime-field-certified.json"},
         {"exit": 0, "json_paths": {"status": "checks_passed"},
          "json_contains": {"checks": [{"code": "certificate-invalid", "status": "pass", "exact": True},
                                       {"code": "p-composite", "status": "pass", "exact": True}]}},
         command="check", timeout=120)
    broken = copy.deepcopy(cert_p)
    inner = next(f for f in broken["factors"] if "certificate" in f)["certificate"]
    two = next(f for f in inner["factors"] if f["prime"] == "2")
    two["witness"] = "4"  # a square: 4^((q2-1)/2) = 1 (mod q2)
    assert check_certificate(p68, broken).startswith(f"factor {q2}: witness 4 for 2: gcd")
    case("C069-prime-field-certificate-invalid",
         "design §4.4: a nested certificate that fails at its second level does not verify",
         {"copy": "{here}/C068-prime-field-certified.json", "set": {"field.p_certificate": broken}},
         refused("certificate-invalid", "invalid", 2), command="check", timeout=120)

    case("C070-subfield-curve-estimates",
         "design §5.2 step 5: F2 estimates kic on a curve over a subfield and its rho-negation pair",
         {"copy": "{here}/C059-subfield-curve-paired.json", "set": {"method.fidelity": "F2"}},
         {"exit": 0, "json_paths": {"status": "estimated", "estimate.ic.pipeline": "kic",
                                    "estimate.rho.pipeline": "rho-negation"}},
         argv_tail=one_rep)
    return files, cases


def texts() -> dict[str, str]:
    selftest()
    files, cases = build()
    out = {f"params/{name}": json.dumps(doc, indent=1) + "\n" for name, doc in files.items()}
    out["cases.json"] = json.dumps({
        "suite": "ic tool programme conformance suite v2, B2b's cases",
        "design": "research/ic_tool_program/design/schema-v2.md §3.3, §4.2, §4.4, §5.1-§5.3; "
                  "research/ic_tool_program/rounds/B2b-subfield-kic-certificates/PROTOCOL.md",
        "includes": "../v1 (B0), ../v2 (B1), ../v2-b2 (B2), ../v2-b3 (B3); ../run.py runs every step's cases",
        "rules": [
            "The rules of ../v2/cases.json hold. {here} is this directory's params/; {cases} is "
            "../v2/params/.",
            "B2b's expectations are its protocol's, made exact here: the codes, the exit statuses, the "
            "pipelines, and the report keys B2b adds (operation ic_single_target, ic.pipeline, "
            "speedup_eligible, the disclosure rho-frobenius-unused, the check certificate-invalid when it "
            "passes).",
        ],
        "generator": {"path": "research/ic_tool_program/conformance/v2-b2b/make_cases.py",
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
        raise SystemExit("cases.json exists; B2b's cases are frozen (use --check)")
    for rel, t in out.items():
        path = HERE / rel
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, "x") as fh:
            fh.write(t)
    print(f"{len(out)} files written")


if __name__ == "__main__":
    main()
