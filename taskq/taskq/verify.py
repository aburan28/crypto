"""Independent certificate verification: `python -m taskq.verify certificate.json`.

A solver that claims a solve writes `certificate.json` into $TASKQ_OUTPUT_DIR.
The format is crypto-autoresearcher's (docs/claims-and-verification.md), so
one certificate serves both repos:

    {"kind": "discrete_log",
     "curve": {"field": "prime",  "p": P, "a": A, "b": B}                  y^2 = x^3 + ax + b
           or {"field": "binary", "m": M, "modulus": [131, 13, 2, 1, 0],
               "a": A, "b": B}                                           y^2 + xy = x^3 + ax^2 + b
     "statement": {"P": [x, y], "Q": [x, y], "k": K, "n": N (optional)}}

or `{"kind": "none"}` when the run claims nothing. Integers may be JSON ints,
decimal strings or 0x-hex strings; binary-field elements are polynomial-basis
bit vectors and `modulus` lists the exponents of the reduction polynomial.
If `curve` is absent, `statement.curve` is accepted too (the autoresearcher
nests it there in some records).

This file shares no code with any solver, and the worker runs it in a separate
process. It checks that P and Q are on the curve, that P is not the point at
infinity, that k*P == Q, and that n*P == O when an order n is given. Prints
one JSON line and exits 0 verified, 1 refuted, 2 no claim, 3 malformed.
"""
from __future__ import annotations

import json
import sys
from typing import Any

EXIT = {"verified": 0, "refuted": 1, "no_claim": 2, "error": 3}


def _int(v: Any) -> int:
    if isinstance(v, bool):
        raise ValueError("boolean is not an integer")
    if isinstance(v, int):
        return v
    s = str(v).strip().lower()
    return int(s, 16) if s.startswith("0x") else int(s, 10)


# -- prime field -------------------------------------------------------------

class PrimeCurve:
    def __init__(self, p: int, a: int, b: int):
        if p < 3:
            raise ValueError("p must be an odd prime")
        self.p, self.a, self.b = p, a % p, b % p

    def on_curve(self, P) -> bool:
        if P is None:
            return True
        x, y = P
        p = self.p
        return 0 <= x < p and 0 <= y < p and (y * y - (x * x * x + self.a * x + self.b)) % p == 0

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


# -- binary field ------------------------------------------------------------

def _gf2_mul(a: int, b: int, mod: int, m: int) -> int:
    r = 0
    while b:
        if b & 1:
            r ^= a
        b >>= 1
        a <<= 1
        if (a >> m) & 1:
            a ^= mod
    return r


def _gf2_inv(a: int, mod: int, m: int) -> int:
    """a^(2^m - 2) by square-and-multiply: no extended-gcd shortcut to share."""
    if a == 0:
        raise ZeroDivisionError("inverse of 0 in GF(2^m)")
    result, base, e = 1, a, (1 << m) - 2
    while e:
        if e & 1:
            result = _gf2_mul(result, base, mod, m)
        base = _gf2_mul(base, base, mod, m)
        e >>= 1
    return result


class BinaryCurve:
    def __init__(self, m: int, modulus: list[int], a: int, b: int):
        exps = sorted({int(e) for e in modulus}, reverse=True)
        if not exps or exps[0] != m or 0 not in exps:
            raise ValueError("modulus must have degree m and a constant term")
        self.m, self.mod = m, sum(1 << e for e in exps)
        if a >> m or b >> m or b == 0:
            raise ValueError("a, b must be field elements and b != 0")
        self.a, self.b = a, b

    def _mul(self, x: int, y: int) -> int:
        return _gf2_mul(x, y, self.mod, self.m)

    def on_curve(self, P) -> bool:
        if P is None:
            return True
        x, y = P
        if x >> self.m or y >> self.m:
            return False
        mul = self._mul
        x2 = mul(x, x)
        lhs = mul(y, y) ^ mul(x, y)
        rhs = mul(x2, x) ^ mul(self.a, x2) ^ self.b
        return lhs == rhs

    def add(self, P, Q):
        if P is None:
            return Q
        if Q is None:
            return P
        mul = self._mul
        (x1, y1), (x2, y2) = P, Q
        if x1 == x2:
            # -P = (x, x + y): Q == -P gives O; this includes doubling when x == 0.
            if y2 == x1 ^ y1:
                return None
            # otherwise Q == P (a given x has only y and x + y above it)
            lam = x1 ^ mul(y1, _gf2_inv(x1, self.mod, self.m))
            x3 = mul(lam, lam) ^ lam ^ self.a
            y3 = mul(x1, x1) ^ mul(lam ^ 1, x3)
            return x3, y3
        lam = mul(y1 ^ y2, _gf2_inv(x1 ^ x2, self.mod, self.m))
        x3 = mul(lam, lam) ^ lam ^ x1 ^ x2 ^ self.a
        y3 = mul(lam, x1 ^ x3) ^ x3 ^ y1
        return x3, y3


def scalar_mul(curve, k: int, P):
    if k < 0:
        raise ValueError("k must be non-negative")
    R, A = None, P
    while k:
        if k & 1:
            R = curve.add(R, A)
        A = curve.add(A, A)
        k >>= 1
    return R


# -- entry points ------------------------------------------------------------

def _curve(spec: dict[str, Any]):
    field = spec.get("field", "prime")
    if field == "prime":
        return PrimeCurve(_int(spec["p"]), _int(spec["a"]), _int(spec["b"]))
    if field == "binary":
        return BinaryCurve(int(spec["m"]), spec["modulus"], _int(spec["a"]), _int(spec["b"]))
    raise ValueError(f"unknown field {field!r}")


def verify_certificate(cert: dict[str, Any]) -> dict[str, Any]:
    kind = cert.get("kind", "none")
    if kind == "none":
        return {"status": "no_claim", "kind": "none"}
    if kind != "discrete_log":
        return {"status": "error", "kind": kind,
                "detail": f"no built-in verifier for kind {kind!r}"}
    try:
        st = cert["statement"]
        E = _curve(cert.get("curve") or st["curve"])
        P = tuple(_int(c) for c in st["P"])
        Q = tuple(_int(c) for c in st["Q"])
        k = _int(st["k"])
    except (KeyError, TypeError, ValueError) as err:
        return {"status": "error", "kind": kind, "detail": f"malformed certificate: {err}"}
    base = {"kind": kind, "verifier": "taskq.verify/independent-recompute"}
    for name, pt in (("P", P), ("Q", Q)):
        if len(pt) != 2 or not E.on_curve(pt):
            return {**base, "status": "refuted", "detail": f"{name} is not on the curve"}
    if "n" in st:
        n = _int(st["n"])
        if scalar_mul(E, n, P) is not None:
            return {**base, "status": "refuted", "detail": "n*P != O"}
        k %= n
    ok = scalar_mul(E, k, P) == Q
    return {**base, "status": "verified" if ok else "refuted",
            "detail": None if ok else "k*P != Q"}


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if len(argv) != 1:
        print("usage: python -m taskq.verify certificate.json", file=sys.stderr)
        return EXIT["error"]
    try:
        with open(argv[0]) as fh:
            cert = json.load(fh)
        res = verify_certificate(cert)
    except (OSError, ValueError) as err:
        res = {"status": "error", "detail": f"unreadable certificate: {err}"}
    print(json.dumps(res, sort_keys=True))
    return EXIT[res["status"]]


if __name__ == "__main__":
    sys.exit(main())
