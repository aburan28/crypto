#!/usr/bin/env python3
"""Freeze the n=83 a=1 public target with standalone GF(2^83) arithmetic.

The scalar is a validation-only sidecar: it is generated here, never sent
to either solver, and only used by the independent replay.  The curve
parameters come from the retained base header + the Koblitz trace
recurrence (pinned by the landed n=61/71/73 fixtures).
"""
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
N = 85
MODULUS_LOW_TERMS = [0, 1, 2, 8]
A = 0
B = 1
R = 14351744671810121
COFACTOR = 2695534732
GENERATOR = (20635810963392012626713602, 3844135993074684155362043)
FIXTURE_SCALAR = 9876543210987654  # fresh; distinct from n41/53/61/71/73/83 sidecars

MODULUS = 1 << N
for term in MODULUS_LOW_TERMS:
    MODULUS |= 1 << term


def gf_mul(a, b):
    out = 0
    while b:
        if b & 1:
            out ^= a
        b >>= 1
        a <<= 1
        if (a >> N) & 1:
            a ^= MODULUS
    return out


def gf_square(a):
    return gf_mul(a, a)


def gf_inverse(value):
    exponent = (1 << N) - 2
    result = 1
    base = value
    while exponent:
        if exponent & 1:
            result = gf_mul(result, base)
        base = gf_square(base)
        exponent >>= 1
    return result


def on_curve(point):
    if point is None:
        return True
    x, y = point
    left = gf_square(y) ^ gf_mul(x, y)
    right = gf_mul(gf_square(x), x) ^ gf_mul(A, gf_square(x)) ^ B
    return left == right


def add(left, right):
    if left is None:
        return right
    if right is None:
        return left
    x, y = left
    u, v = right
    if x == u:
        if y != v or x == 0:
            return None
        slope = x ^ gf_mul(y, gf_inverse(x))
        xx = gf_square(slope) ^ slope ^ A
        yy = gf_square(x) ^ gf_mul(slope ^ 1, xx)
        return xx, yy
    slope = gf_mul(y ^ v, gf_inverse(x ^ u))
    xx = gf_square(slope) ^ slope ^ x ^ u ^ A
    yy = gf_mul(slope, x ^ xx) ^ xx ^ y
    return xx, yy


def scalar_mul(k, point):
    result = None
    addend = point
    while k:
        if k & 1:
            result = add(result, addend)
        addend = add(addend, addend)
        k >>= 1
    return result


def main():
    assert on_curve(GENERATOR), "generator must be on the curve"
    identity = scalar_mul(R, GENERATOR)
    assert identity is None, "[r]G must be infinity"
    q = scalar_mul(FIXTURE_SCALAR, GENERATOR)
    assert q is not None and on_curve(q), "target must be an affine curve point"
    assert scalar_mul(R, q) is None, "[r]Q must be infinity"

    # Trace / order cross-check from the recurrence (t_1 = -1 for a = 0).
    t_prev, t_cur = 2, -1
    for _ in range(2, N + 1):
        t_prev, t_cur = t_cur, t_cur * -1 - 2 * t_prev
    order = (1 << N) + 1 - t_cur
    assert order == R * COFACTOR, f"order mismatch: {order} != {R * COFACTOR}"

    frozen = HERE / "frozen"
    frozen.mkdir(exist_ok=True)
    fixture = {
        "kind": "koblitz_public_fixture_point",
        "n": N,
        "a": A,
        "field_modulus_low_terms": MODULUS_LOW_TERMS,
        "field_value_encoding": "base-10 strings for u128 polynomial-basis coordinates",
        "subgroup_order": str(R),
        "cofactor": str(COFACTOR),
        "curve_order": str(order),
        "curve_trace_over_f2n": str(t_cur),
        "generator": [str(GENERATOR[0]), str(GENERATOR[1])],
        "public_target": [str(q[0]), str(q[1])],
        "fixture_scalar_validation_only": str(FIXTURE_SCALAR),
        "target_count": 1,
        "target_generation_is_outside_both_online_intervals": True,
        "frobenius_order_on_subgroup": N,
    }
    (frozen / "fixture.json").write_text(json.dumps(fixture, sort_keys=True) + "\n")
    (frozen / "target_points.jsonl").write_text(
        json.dumps([str(q[0]), str(q[1])], separators=(",", ":")) + "\n"
    )
    print(json.dumps({
        "public_target_q": [str(q[0]), str(q[1])],
        "fixture_scalar_validation_only": str(FIXTURE_SCALAR),
        "curve_order": str(order),
        "checks": {
            "generator_on_curve": True,
            "subgroup_annihilates_generator": True,
            "target_on_curve": True,
            "subgroup_annihilates_target": True,
        },
    }, indent=1))


if __name__ == "__main__":
    main()
