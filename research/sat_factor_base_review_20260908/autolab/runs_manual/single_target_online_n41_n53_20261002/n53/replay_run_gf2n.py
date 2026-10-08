#!/usr/bin/env python3
"""Independent GF(2^n) elliptic-curve replay for one autolab run (n from logs)."""
import json
import platform
from pathlib import Path

HERE = Path(__file__).resolve().parent
RUN = HERE.parent


def gf_mul(a, b):
    z = 0
    while b:
        if b & 1:
            z ^= a
        b >>= 1
        a <<= 1
        if a & (1 << N):
            a ^= MOD
    return z


def gf_pow(a, e):
    z = 1
    while e:
        if e & 1:
            z = gf_mul(z, a)
        a = gf_mul(a, a)
        e >>= 1
    return z


def gf_inv(a):
    if not a:
        raise ZeroDivisionError
    return gf_pow(a, (1 << N) - 2)


def on_curve(P):
    if P is None:
        return True
    x, y = P
    return gf_mul(y, y) ^ gf_mul(x, y) == gf_mul(gf_mul(x, x), x) ^ 1


def add(P, Q):
    if P is None:
        return Q
    if Q is None:
        return P
    x1, y1 = P
    x2, y2 = Q
    if x1 == x2:
        if y1 ^ y2 == x1:
            return None
        if y1 != y2 or x1 == 0:
            raise ValueError("invalid equal-x group-addition case")
        lam = x1 ^ gf_mul(y1, gf_inv(x1))
        x3 = gf_mul(lam, lam) ^ lam
        y3 = gf_mul(x1, x1) ^ gf_mul(lam ^ 1, x3)
        return (x3, y3)
    lam = gf_mul(y1 ^ y2, gf_inv(x1 ^ x2))
    x3 = gf_mul(lam, lam) ^ lam ^ x1 ^ x2
    y3 = gf_mul(lam, x1 ^ x3) ^ x3 ^ y1
    return (x3, y3)


def mul(k, P):
    R = None
    while k:
        if k & 1:
            R = add(R, P)
        P = add(P, P)
        k >>= 1
    return R


def poly_mod(a, f):
    while a.bit_length() >= f.bit_length():
        a ^= f << (a.bit_length() - f.bit_length())
    return a


def poly_gcd(a, b):
    while b:
        a, b = b, poly_mod(a, b)
    return a


def poly_squaremod(a, f):
    z = 0
    i = 0
    while a:
        if a & 1:
            z |= 1 << (2 * i)
        a >>= 1
        i += 1
    return poly_mod(z, f)


rows = {}
for name in ("direct", "rho"):
    lines = (RUN / f"logs/{name}.stdout.jsonl").read_text().splitlines()
    rows[name] = [json.loads(s) for s in lines if s.strip()]

rho = next(d for d in rows["rho"] if d.get("kind") == "rho_public_fixture")
ic = next(d for d in rows["direct"] if d.get("target_online_wall_ms") is not None)

N = int(rho["n"])
terms = sorted(int(t) for t in rho["field_modulus_low_terms"])
assert terms[-1] < N and 0 in terms, terms
MOD = (1 << N) | sum(1 << t for t in terms)
X = 2

# Rabin irreducibility test for prime degree N.
assert all(N % d for d in range(2, int(N**0.5) + 1)), "degree must be prime for this test"
x2 = poly_squaremod(X, MOD)
field_irreducible = poly_gcd(x2 ^ X, MOD) == 1
xpow = X
for _ in range(N):
    xpow = poly_squaremod(xpow, MOD)
field_irreducible = field_irreducible and xpow == X

G = tuple(map(int, rho["generator"]))
Q = tuple(map(int, ic["published_q"]))
r = int(rho["subgroup_order"])
s_ic = int(ic["recovered_fixture_scalar"])
s_rho = int(rho["recovered_fixture_scalar"])

checks = {
    "field_irreducible": field_irreducible,
    "generator_on_curve": on_curve(G),
    "target_on_curve": on_curve(Q),
    "subgroup_order_annihilates_generator": mul(r, G) is None,
    "ic_scalar_replays_target": mul(s_ic, G) == Q,
    "rho_scalar_replays_target": mul(s_rho, G) == Q,
    "ic_and_rho_target_identical": ic["published_q"] == rho["published_q"],
    "ic_and_rho_scalar_identical": s_ic == s_rho,
}
receipt = {
    "kind": "independent_binary_curve_scalar_replay",
    "run_id": RUN.name,
    "method": "standalone Python polynomial-basis GF(2^%d) arithmetic and binary Weierstrass group law; separate from Rust producers" % N,
    "curve": {
        "field_modulus": " + ".join(["x^%d" % N] + ["x^%d" % t for t in reversed(terms)]),
        "weierstrass_model": "y^2 + x*y = x^3 + 1",
        "a": 0,
        "b": 1,
    },
    "generator": list(G),
    "target": list(Q),
    "subgroup_order": r,
    "ic_recovered_scalar": s_ic,
    "rho_recovered_scalar": s_rho,
    "checks": checks,
    "platform": platform.platform(),
    "python": platform.python_version(),
    "status": "PASS" if all(checks.values()) else "FAIL",
}
out = HERE / "validation.json"
out.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
print(json.dumps(receipt, indent=2, sort_keys=True))
