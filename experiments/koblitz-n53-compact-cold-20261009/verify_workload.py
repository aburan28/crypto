#!/usr/bin/env python3
"""Independently regenerate the frozen n53 point before any timed run."""

import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
W = json.loads((HERE / "workload.json").read_text())
N = W["n"]
MODULUS = (1 << N) | sum(1 << i for i in W["field_modulus_low_terms"])
MASK = (1 << N) - 1


def mul(a, b):
    result = 0
    while b:
        if b & 1:
            result ^= a
        b >>= 1
        a <<= 1
        if a & (1 << N):
            a ^= MODULUS
    return result & MASK


def inv(a):
    assert a
    u, v, left, right = a, MODULUS, 1, 0
    while u != 1:
        shift = u.bit_length() - v.bit_length()
        if shift < 0:
            u, v, left, right = v, u, right, left
            shift = -shift
        u ^= v << shift
        left ^= right << shift
    while left.bit_length() > N:
        left ^= MODULUS << (left.bit_length() - N - 1)
    return left


def add(p, q):
    if p is None:
        return q
    if q is None:
        return p
    x1, y1 = p
    x2, y2 = q
    if x1 == x2:
        if y1 ^ y2 == x1:
            return None
        assert y1 == y2 and x1
        slope = x1 ^ mul(y1, inv(x1))
        x3 = mul(slope, slope) ^ slope
        y3 = mul(x1, x1) ^ mul(slope ^ 1, x3)
    else:
        slope = mul(y1 ^ y2, inv(x1 ^ x2))
        x3 = mul(slope, slope) ^ slope ^ x1 ^ x2
        y3 = mul(slope, x1 ^ x3) ^ x3 ^ y1
    return x3, y3


def scale(k, point):
    result = None
    while k:
        if k & 1:
            result = add(result, point)
        point = add(point, point)
        k >>= 1
    return result


def on_curve(point):
    x, y = point
    return mul(y, y) ^ mul(x, y) == mul(mul(x, x), x) ^ 1


def main():
    order = W["subgroup_order"]
    seed = W["target_seed_ascii"].encode("ascii")
    scalar = int.from_bytes(hashlib.sha256(seed).digest(), "big") % (order - 1) + 1
    generator = tuple(W["generator"])
    target = scale(scalar, generator)
    assert scalar == W["verification_scalar"]
    assert target == tuple(W["primary_target"])
    assert scale(order, generator) is None
    assert on_curve(generator) and on_curve(target)
    assert (HERE / "target_points.jsonl").read_text() == json.dumps(list(target), separators=(",", ":")) + "\n"
    print(json.dumps({"status": "PASS", "target": target, "scalar": scalar}, sort_keys=True))


if __name__ == "__main__":
    main()
