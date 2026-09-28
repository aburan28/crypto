#!/usr/bin/env python3
"""Independent polynomial-product/Euclid replay of the n13 K0-mu4 map."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from audit import HERE, canonical, preflight

O = None


def mul(a: int, b: int, n: int, modulus: int) -> int:
    product = 0
    for i in range(n):
        if (a >> i) & 1:
            product ^= b << i
    for i in range(2 * n - 2, n - 1, -1):
        if (product >> i) & 1:
            product ^= modulus << (i - n)
    assert product < 1 << n
    return product


def square(a: int, n: int, modulus: int) -> int:
    return mul(a, a, n, modulus)


def inverse(a: int, n: int, modulus: int) -> int:
    assert a
    u, v, g, h = a, modulus, 1, 0
    while u != 1:
        assert u != 0
        shift = u.bit_length() - v.bit_length()
        if shift < 0:
            u, v = v, u
            g, h = h, g
            shift = -shift
        u ^= v << shift
        g ^= h << shift
    for i in range(g.bit_length() - 1, n - 1, -1):
        if (g >> i) & 1:
            g ^= modulus << (i - n)
    assert mul(a, g, n, modulus) == 1
    return g


def on_k0(p, n: int, modulus: int) -> bool:
    if p is O:
        return True
    x, y = p
    return square(y, n, modulus) ^ mul(x, y, n, modulus) == mul(square(x, n, modulus), x, n, modulus) ^ 1


def on_mu4(v, n: int, modulus: int) -> bool:
    a, b, c, d = v
    return any(v) and square(a ^ c, n, modulus) == mul(b, d, n, modulus) and square(b ^ d, n, modulus) == mul(a, c, n, modulus)


def image(p, n: int, modulus: int):
    if p is O:
        return (1, 1, 0, 1)
    x, y = p
    x2 = square(x, n, modulus)
    return (x2, x2 ^ y, 1, x2 ^ x ^ y)


def inverse_image(v, n: int, modulus: int):
    assert on_mu4(v, n, modulus)
    a, b, c, d = v
    if c == 0:
        assert a == b == d != 0
        return O
    k = inverse(c, n, modulus)
    p = (mul(b ^ d, k, n, modulus), mul(a ^ b, k, n, modulus))
    assert on_k0(p, n, modulus)
    return p


def add(p, q, n: int, modulus: int):
    if p is O:
        return q
    if q is O:
        return p
    x1, y1 = p
    x2, y2 = q
    if x1 == x2:
        if y1 ^ y2 == x1:
            return O
        assert y1 == y2 and x1
        slope = x1 ^ mul(y1, inverse(x1, n, modulus), n, modulus)
        xr = square(slope, n, modulus) ^ slope
        yr = square(x1, n, modulus) ^ mul(slope ^ 1, xr, n, modulus)
    else:
        slope = mul(y1 ^ y2, inverse(x1 ^ x2, n, modulus), n, modulus)
        xr = square(slope, n, modulus) ^ slope ^ x1 ^ x2
        yr = mul(slope, x1 ^ xr, n, modulus) ^ xr ^ y1
    result = (xr, yr)
    assert on_k0(result, n, modulus)
    return result


def independently_enumerate(n: int, modulus: int):
    # This verifier enumerates the complete Artin-Schreier root table, not
    # producer's trace/half-trace algorithm; square roots are likewise tabulated.
    roots = {}
    square_roots = {}
    for z in range(1 << n):
        roots.setdefault(square(z, n, modulus) ^ z, []).append(z)
        key = square(z, n, modulus)
        assert key not in square_roots
        square_roots[key] = z
    assert len(square_roots) == 1 << n
    src = {(0, 1)}
    for x in range(1, 1 << n):
        xinv = inverse(x, n, modulus)
        rhs = x ^ square(xinv, n, modulus)
        for z in roots.get(rhs, []):
            p = (x, mul(x, z, n, modulus))
            assert on_k0(p, n, modulus)
            src.add(p)
    chart = set()
    for a in range(1 << n):
        x = square_roots[a]
        if x == 0:
            chart.add((0, 1, 1, 1))
            continue
        rhs = mul(square(a ^ 1, n, modulus), inverse(a, n, modulus), n, modulus)
        for z in roots.get(rhs, []):
            b = mul(x, z, n, modulus)
            v = (a, b, 1, b ^ x)
            assert on_mu4(v, n, modulus)
            chart.add(v)
    return src, chart


def verify(producer: dict, inp: dict) -> dict:
    n, modulus = inp["toy"]["n"], inp["toy"]["modulus"]
    src, chart = independently_enumerate(n, modulus)
    images = {image(p, n, modulus) for p in src}
    assert images == chart
    assert all(inverse_image(image(p, n, modulus), n, modulus) == p for p in src)
    assert inverse_image(image(O, n, modulus), n, modulus) is O
    assert not on_mu4((1, 0, 0, 1), n, modulus)
    damaged = list(image((0, 1), n, modulus))
    damaged[1] ^= 1
    assert not on_mu4(damaged, n, modulus)
    try:
        assert image(O, n, modulus)[2] != 0, "O has no affine inverse chart"
    except AssertionError as e:
        assert str(e) == "O has no affine inverse chart"
    else:
        raise AssertionError("O inverse guard failed")
    t = (1, 0)
    chain = []
    p = t
    for _ in range(4):
        p = add(p, t, n, modulus)
        chain.append(p)
    assert chain == [(0, 1), (1, 1), O, t]
    result = {"status": "PASS", "source_affine_points": len(src),
              "mu4_affine_chart_points": len(chart), "projective_points": len(src) + 1,
              "source_points_sha256": hashlib.sha256(canonical(sorted(src))).hexdigest(),
              "mu4_chart_sha256": hashlib.sha256(canonical(sorted(chart))).hexdigest(),
              "controls": {name: "PASS" for name in inp["toy_controls"]},
              "order4_chain": [[0, 1], [1, 1], "O", [1, 0]]}
    assert result == producer["toy_mu4_map"]
    return {"schema": "symbolic-oaware-width-gate-independent-v1",
            "producer_toy_sha256": hashlib.sha256(canonical(producer["toy_mu4_map"])).hexdigest(),
            "toy_mu4_map": result}


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--producer", required=True, type=Path)
    p.add_argument("--output", required=True, type=Path)
    args = p.parse_args()
    _, inp = preflight()
    assert not args.output.exists()
    producer = json.loads(args.producer.read_text())
    assert producer["schema"] == "symbolic-oaware-width-gate-receipt-v1"
    out = verify(producer, inp)
    args.output.write_bytes(canonical(out))
    print("INDEPENDENT_MAP_PASS", hashlib.sha256(args.output.read_bytes()).hexdigest())


if __name__ == "__main__":
    main()
