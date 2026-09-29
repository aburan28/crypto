#!/usr/bin/env python3
"""ICV1 curve identities: the reference implementation.

Every elliptic curve this repository names is named by its ICV1 identity,
specified in docs/curves/ICV1.md.  This module is the reference: the Rust
module `crypto::cryptanalysis::curve_id` is tested against vectors this
file produces, and `docs/curves/registry.json` is built with it.

    python3 scripts/curve_id.py koblitz 0 41
    python3 scripts/curve_id.py binary 27 --a 1 --b 0x845462 --order 134215648
    python3 scripts/curve_id.py prime 10935329 --a 5320418 --b 8535318 --order 10935433
    python3 scripts/curve_id.py resolve 'K_0 / GF(2^41)'
    python3 scripts/curve_id.py selftest

Only the standard library is used, so the file runs on any CI image.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
from pathlib import Path

VERSION = "1"
REPO = Path(__file__).resolve().parent.parent
REGISTRY = REPO / "docs" / "curves" / "registry.json"

# ---------------------------------------------------------------------------
# GF(2)[x] arithmetic on Python ints: bit i is the coefficient of x^i.
# ---------------------------------------------------------------------------


def _deg(f: int) -> int:
    return f.bit_length() - 1


def _mod(a: int, f: int) -> int:
    df = _deg(f)
    while a and _deg(a) >= df:
        a ^= f << (_deg(a) - df)
    return a


def _mulmod(a: int, b: int, f: int) -> int:
    r = 0
    a = _mod(a, f)
    while b:
        if b & 1:
            r ^= a
        b >>= 1
        a <<= 1
        if a >> _deg(f) & 1:
            a ^= f
    return r


def _gcd2(a: int, b: int) -> int:
    while b:
        a, b = b, _mod(a, b)
    return a


def _inv2(a: int, f: int) -> int:
    """Inverse of a non-zero `a` modulo irreducible `f` (extended Euclid)."""
    r0, r1 = f, _mod(a, f)
    s0, s1 = 0, 1
    if r1 == 0:
        raise ZeroDivisionError("inverse of zero in GF(2^m)")
    while r1 != 1:
        q = 0
        r = r0
        while r and _deg(r) >= _deg(r1):
            shift = _deg(r) - _deg(r1)
            q ^= 1 << shift
            r ^= r1 << shift
        r0, r1 = r1, r
        s0, s1 = s1, s0 ^ _clmul(q, s1)
        if r1 == 0:
            raise ValueError("modulus is not irreducible")
    return _mod(s1, f)


def _clmul(a: int, b: int) -> int:
    r = 0
    while b:
        if b & 1:
            r ^= a
        b >>= 1
        a <<= 1
    return r


def _prime_factors(n: int) -> list[int]:
    out, d = [], 2
    while d * d <= n:
        if n % d == 0:
            out.append(d)
            while n % d == 0:
                n //= d
        d += 1
    if n > 1:
        out.append(n)
    return out


def is_irreducible_f2(f: int) -> bool:
    """Rabin's test over GF(2)."""
    n = _deg(f)
    if n < 1 or not f & 1:
        return n == 1
    x = 2

    def x_pow_2k(k: int) -> int:
        y = x
        for _ in range(k):
            y = _mulmod(y, y, f)
        return y

    if x_pow_2k(n) != _mod(x, f):
        return False
    for p in _prime_factors(n):
        if _gcd2(f, x_pow_2k(n // p) ^ x) != 1:
            return False
    return True


def _sparse_lows(n: int):
    """Odd `low` of weight at most four below x^n, in increasing order."""
    yield 1
    for top in range(1, n):
        below = [0]
        for i in range(1, top):
            below.append(1 << i)
            for j in range(i + 1, top):
                below.append(1 << i | 1 << j)
        for rest in sorted(below):
            yield 1 | rest | 1 << top


def find_irreducible_sparse(n: int) -> int:
    """The modulus `KoblitzCurve::new` and the boundary ladder use.

    Mirrors `koblitz_index_calculus::find_irreducible_sparse`: the
    numerically least `x^n + low` with `low` odd and of weight at most
    four that is irreducible.  The Rust search stops at n = 63 (one
    machine word); this one does not, which is the registry's rule for a
    Koblitz degree with no published polynomial.  Returns the full
    polynomial as an int.
    """
    if n < 1:
        raise ValueError("degree must be positive")
    if n == 1:
        return 0b11
    for low in _sparse_lows(n):
        f = 1 << n | low
        if is_irreducible_f2(f):
            return f
    raise ValueError(f"no sparse irreducible of degree {n}")


# ---------------------------------------------------------------------------
# Group orders that follow from the curve's definition.
# ---------------------------------------------------------------------------


def koblitz_trace(a: int, n: int) -> int:
    """Trace of Frobenius of K_a : y^2 + xy = x^3 + a x^2 + 1 over GF(2^n)."""
    t1 = 1 if a == 1 else -1  # #K_1(F_2) = 2, #K_0(F_2) = 4
    t_prev, t = 2, t1
    for _ in range(n - 1):
        t_prev, t = t, t1 * t - 2 * t_prev
    return t


def koblitz_order(a: int, n: int) -> int:
    return 2**n + 1 - koblitz_trace(a, n)


# ---------------------------------------------------------------------------
# The identity.
# ---------------------------------------------------------------------------


def _hex(v: int) -> str:
    return "0x%x" % v


def _sha256(s: str) -> str:
    return hashlib.sha256(s.encode("ascii")).hexdigest()


def _canonical(obj: dict) -> str:
    return json.dumps(obj, sort_keys=True, separators=(",", ":"), ensure_ascii=True)


def binary_field(m: int, modulus: int) -> str:
    if _deg(modulus) != m:
        raise ValueError("modulus degree must equal m")
    return "f2m-%d-%s" % (m, _sha256("f2m-modulus:" + _hex(modulus))[:8])


def prime_field(p: int) -> str:
    return "fp-%d" % p


def _slug(field_tag: str, trace: int, model: str) -> str:
    t = ("tm%d" % -trace) if trace < 0 else ("t%d" % trace)
    return "icv1-%s-%s-%s" % (field_tag, t, model[:8])


def _check_hasse(q: int, order: int) -> None:
    """|q + 1 - #E| <= 2 sqrt q, or the order is not the whole group's."""
    t = q + 1 - order
    if t * t > 4 * q:
        raise ValueError("order %d violates the Hasse bound over a field of %d elements"
                         % (order, q))


def _identity(field: str, trace: int, order: int, j: str, end: str, model_json: str,
              field_tag: str, level: str = "unk", path: str = "r") -> dict:
    model = _sha256(model_json)
    canonical = "ICV1:%s:%d:%d:%s:%s:%s:%s:%s" % (
        field, trace, order, j, end, level, path, model[:12])
    return {
        "icv1": canonical,
        "slug": _slug(field_tag, trace, model),
        "model_json": model_json,
        "model_sha256": model,
        "field": field,
        "trace": trace,
        "order": order,
        "j": j,
        "end": end,
        "level": level,
        "path": path,
    }


def binary_id(m: int, modulus: int, a: int, b: int, order: int, end: str = "unk") -> dict:
    """y^2 + xy = x^3 + a x^2 + b over GF(2^m) = GF(2)[x]/(modulus).

    `a` and `b` are polynomial-basis field elements (bit i = x^i).  The
    order is #E(GF(2^m)), the whole group.
    """
    if b == 0:
        raise ValueError("b = 0 is singular")
    a, b = _mod(a, modulus), _mod(b, modulus)
    field = binary_field(m, modulus)
    model_json = _canonical({
        "a": _hex(a), "b": _hex(b), "field": field,
        "form": "y^2+xy=x^3+a*x^2+b", "modulus": _hex(modulus), "v": VERSION,
    })
    _check_hasse(2**m, order)
    trace = 2**m + 1 - order
    j = _hex(_inv2(b, modulus))
    return _identity(field, trace, order, j, end, model_json, "f2m%d" % m)


def koblitz_id(a: int, n: int, modulus: int | None = None) -> dict:
    """K_a over GF(2^n).  End(E) = Z[tau] = O_{Q(sqrt -7)}, so end = -7."""
    if a not in (0, 1):
        raise ValueError("Koblitz a is 0 or 1")
    if modulus is None:
        modulus = find_irreducible_sparse(n)
    return binary_id(n, modulus, a, 1, koblitz_order(a, n), end="-7")


def prime_id(p: int, a: int, b: int, order: int, end: str = "unk") -> dict:
    """y^2 = x^3 + a x + b over GF(p), p an odd prime > 3."""
    a, b = a % p, b % p
    disc = (4 * a**3 + 27 * b**2) % p
    if disc == 0:
        raise ValueError("singular curve")
    field = prime_field(p)
    model_json = _canonical({
        "a": str(a), "b": str(b), "field": field,
        "form": "y^2=x^3+a*x+b", "p": str(p), "v": VERSION,
    })
    _check_hasse(p, order)
    trace = p + 1 - order
    j = 1728 * 4 * pow(a, 3, p) * pow(disc, p - 2, p) % p
    return _identity(field, trace, order, str(j), end, model_json, "fp%d" % p.bit_length())


# ---------------------------------------------------------------------------
# Resolving names through the registry.
# ---------------------------------------------------------------------------


def normalise_alias(name: str) -> str:
    """The comparison key for a spelling: case, whitespace, braces,
    underscores, surrounding backticks and subscript digits folded away.
    `crypto::cryptanalysis::curve_id::normalise` is the same function."""
    s = name.strip().strip("`")
    s = s.translate(str.maketrans("₀₁₂₃₄₅₆₇₈₉", "0123456789"))
    s = re.sub(r"[\s{}_]", "", s).lower()
    return s


def load_registry(path: Path = REGISTRY) -> dict:
    return json.loads(path.read_text())


def resolve(name: str, registry: dict | None = None) -> dict | None:
    reg = registry or load_registry()
    key = normalise_alias(name)
    for curve in reg["curves"]:
        if key in (normalise_alias(curve["slug"]), normalise_alias(curve["icv1"])):
            return curve
        if any(key == normalise_alias(x) for x in curve.get("aliases", [])):
            return curve
    return None


# ---------------------------------------------------------------------------
# Self-test: vectors whose facts are independent of this file.
# ---------------------------------------------------------------------------


def _count_prime(p: int, a: int, b: int) -> int:
    sq = [0] * p
    for y in range(p):
        sq[y * y % p] += 1
    return 1 + sum(sq[(x * x * x + a * x + b) % p] for x in range(p))


def _trace2(c: int, m: int, f: int) -> int:
    t, s = 0, c
    for _ in range(m):
        t ^= s
        s = _mulmod(s, s, f)
    return t


def _count_binary(m: int, f: int, a: int, b: int) -> int:
    """#E for y^2 + xy = x^3 + a x^2 + b over GF(2^m), one x at a time.

    x = 0 gives one point (y = sqrt b); for x != 0, y = xz turns the
    equation into z^2 + z = rhs/x^2, with two roots iff that trace is 0.
    """
    total = 2  # infinity and (0, sqrt b)
    for x in range(1, 1 << m):
        rhs = _mulmod(_mulmod(x, x, f), x ^ a, f) ^ b
        c = _mulmod(rhs, _inv2(_mulmod(x, x, f), f), f)
        if _trace2(c, m, f) == 0:
            total += 2
    return total


def selftest() -> None:
    # K_0/GF(2^9): the frozen oracle-pricing cells record #E = 508 and the
    # modulus x^9 + x + 1 (polynomial_low_terms [0, 1]).
    assert koblitz_order(0, 9) == 508
    assert find_irreducible_sparse(9) == (1 << 9) | 0b11
    # sect163k1 (SEC 2): #E = 2 * n.
    n163 = 0x4000000000000000000020108A2E0CC0D99F8A5EF
    assert koblitz_order(1, 163) == 2 * n163
    # ECC2K-130: K_0 over GF(2^131), cofactor 4 times a 129-bit prime.
    assert koblitz_order(0, 131) % 4 == 0
    # Brute-force point counts agree with the orders the frozen Round-5
    # ladder records (docs/ic/runs/ic-boundary-ledger-round5-2026-09-22.json).
    assert _count_prime(827, 1, 15) == 823            # bench-10bit
    assert _count_binary(9, (1 << 9) | 0b11, 0, 1) == 508
    assert _count_binary(15, (1 << 15) | 0b11, 1, 0x524B) == 32638
    # The recorded moduli are the sparse search's (low terms per degree).
    for n, low in ((11, [0, 2]), (13, [0, 1, 3, 4]), (18, [0, 3]), (19, [0, 1, 2, 5]),
                   (24, [0, 1, 3, 4]), (27, [0, 1, 2, 5]), (37, [0, 1, 4, 6]),
                   (39, [0, 4]), (41, [0, 3])):
        assert find_irreducible_sparse(n) == (1 << n) | sum(1 << i for i in low), n
    # Koblitz orders from the recurrence equal the recorded ones.
    for a, n, order in ((1, 11, 1982), (0, 13, 8012), (1, 29, 536911222),
                        (0, 39, 549754332404), (0, 41, 2199025563772)):
        assert koblitz_order(a, n) == order, (a, n)
    # j of y^2 = x^3 + b is 0; of y^2 = x^3 + a x is 1728.
    assert prime_id(10007, 0, 7, 10007 + 1)["j"] == "0"
    assert prime_id(10007, 3, 0, 10007 + 1)["j"] == "1728"
    # Binary j = 1/b: with b = 1 it is 1.
    assert koblitz_id(0, 41)["j"] == "0x1"
    # Inversion round-trips.
    f = find_irreducible_sparse(41)
    for v in (2, 3, 0x845462, (1 << 40) | 5):
        assert _mulmod(v, _inv2(v, f), f) == 1
    # Serialisation is stable.
    k = koblitz_id(0, 41)
    assert k["icv1"].startswith("ICV1:f2m-41-") and k["end"] == "-7"
    assert k["slug"].startswith("icv1-f2m41-t")
    print("curve_id selftest: ok")


def main(argv: list[str]) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    k = sub.add_parser("koblitz")
    k.add_argument("a", type=int)
    k.add_argument("n", type=int)
    b = sub.add_parser("binary")
    b.add_argument("m", type=int)
    b.add_argument("--a", type=lambda s: int(s, 0), required=True)
    b.add_argument("--b", type=lambda s: int(s, 0), required=True)
    b.add_argument("--order", type=int, required=True)
    b.add_argument("--modulus", type=lambda s: int(s, 0))
    p = sub.add_parser("prime")
    p.add_argument("p", type=int)
    p.add_argument("--a", type=int, required=True)
    p.add_argument("--b", type=int, required=True)
    p.add_argument("--order", type=int, required=True)
    r = sub.add_parser("resolve")
    r.add_argument("name")
    sub.add_parser("selftest")
    args = ap.parse_args(argv)
    if args.cmd == "selftest":
        selftest()
        return 0
    if args.cmd == "koblitz":
        out = koblitz_id(args.a, args.n)
    elif args.cmd == "binary":
        mod = args.modulus if args.modulus is not None else find_irreducible_sparse(args.m)
        out = binary_id(args.m, mod, args.a, args.b, args.order)
    elif args.cmd == "prime":
        out = prime_id(args.p, args.a, args.b, args.order)
    else:
        out = resolve(args.name)
        if out is None:
            print("unregistered: %s" % args.name, file=sys.stderr)
            return 1
    print(json.dumps(out, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
