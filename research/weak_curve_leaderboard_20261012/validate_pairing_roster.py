#!/usr/bin/env python3
"""Replay exact arithmetic checks for the source-pinned pairing roster."""

import json
from pathlib import Path

HERE = Path(__file__).resolve().parent


def read_int(value):
    return int(value, 0) if isinstance(value, str) and value.startswith("0x") else int(value)


def field(p):
    def add(x, y):
        return ((x[0] + y[0]) % p, (x[1] + y[1]) % p)

    def sub(x, y):
        return ((x[0] - y[0]) % p, (x[1] - y[1]) % p)

    def mul(x, y):
        return ((x[0] * y[0] - x[1] * y[1]) % p,
                (x[0] * y[1] + x[1] * y[0]) % p)

    def inv(x):
        den = pow((x[0] * x[0] + x[1] * x[1]) % p, -1, p)
        return (x[0] * den % p, -x[1] * den % p)

    return add, sub, mul, inv


def point_add(u, v, a, p):
    if u is None:
        return v
    if v is None:
        return u
    add, sub, mul, inv = field(p)
    if u[0] == v[0] and add(u[1], v[1]) == (0, 0):
        return None
    if u == v:
        slope = mul(add(mul((3, 0), mul(u[0], u[0])), a), inv(mul((2, 0), u[1])))
    else:
        slope = mul(sub(v[1], u[1]), inv(sub(v[0], u[0])))
    x = sub(sub(mul(slope, slope), u[0]), v[0])
    y = sub(mul(slope, sub(u[0], x)), u[1])
    return x, y


def point_mul(n, point, a, p):
    result = None
    while n:
        if n & 1:
            result = point_add(result, point, a, p)
        point = point_add(point, point, a, p)
        n >>= 1
    return result


def on_curve(point, a, b, p):
    add, _, mul, _ = field(p)
    x, y = point
    return mul(y, y) == add(add(mul(mul(x, x), x), mul(a, x)), b)


def validate(curve):
    p, r = read_int(curve["p"]), read_int(curve["r"])
    k = curve["embedding_degree"]
    assert p > 3 and r > 3 and pow(p, k, r) == 1
    assert all(pow(p, i, r) != 1 for i in range(1, k))
    h1 = read_int(curve["g1_cofactor"])
    expected_trace = (147946756881789318990833708069417712967
                      if curve["name"].startswith("Ethereum") else -15132376222941642751)
    assert h1 * r == p + 1 - expected_trace
    a = (read_int(curve["a"]), 0)
    b = (read_int(curve["b"]), 0)
    g1 = tuple((read_int(x), 0) for x in curve["g1_generator"])
    assert on_curve(g1, a, b, p) and point_mul(r, g1, a, p) is None
    g2 = tuple(tuple(read_int(x) for x in xy) for xy in curve["g2_generator"])
    add, _, mul, inv = field(p)
    if curve["name"].startswith("Ethereum"):
        b2 = mul((3, 0), inv((9, 1)))
    else:
        b2 = (4, 4)
    assert on_curve(g2, a, b2, p) and point_mul(r, g2, a, p) is None
    print(f'{curve["name"]}: p={p.bit_length()} bits, r={r.bit_length()} bits, '
          f'ord_r(p)={k}, G1 and G2 generators satisfy [r]P=O')


def main():
    curves = json.loads((HERE / "pairing_curves.json").read_text())
    assert len(curves) == 2
    for curve in curves:
        validate(curve)


if __name__ == "__main__":
    main()
