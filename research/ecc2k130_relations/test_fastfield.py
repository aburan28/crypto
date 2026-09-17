#!/usr/bin/env python3
"""`fastfield` must agree with the reference field, element for element."""

import random
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from ecc2k130_point_decomposition import GF2m, IRR131, N131   # noqa: E402
from fastfield import FastGF2m                                 # noqa: E402

REF = GF2m(N131, IRR131)
FAST = FastGF2m(N131, IRR131)
TRIALS = 400


def _rand(rng):
    return rng.getrandbits(N131)


def test_same_reduction_polynomial():
    assert FAST.irr == REF.irr
    assert FAST.mask == REF.mask
    assert FAST.trace_bits == REF.trace_bits


def test_mul_matches_reference():
    rng = random.Random(20260917)
    for _ in range(TRIALS):
        a, b = _rand(rng), _rand(rng)
        assert FAST.mul(a, b) == REF.mul(a, b)


def test_mul_edge_cases():
    for a in (0, 1, 2, FAST.mask):
        for b in (0, 1, 2, FAST.mask):
            assert FAST.mul(a, b) == REF.mul(a, b)


def test_sqr_matches_reference():
    rng = random.Random(1)
    for _ in range(TRIALS):
        a = _rand(rng)
        assert FAST.sqr(a) == REF.sqr(a)


def test_half_trace_matches_reference():
    rng = random.Random(2)
    for _ in range(TRIALS):
        a = _rand(rng)
        assert FAST.half_trace(a) == REF.half_trace(a)


def test_trace_matches_reference():
    rng = random.Random(3)
    for _ in range(TRIALS):
        a = _rand(rng)
        assert FAST.trace(a) == REF.trace(a)


def test_artin_schreier_solution_is_a_root():
    rng = random.Random(4)
    solved = 0
    for _ in range(TRIALS):
        c = _rand(rng)
        t = FAST.solve_artin_schreier(c)
        if t is None:
            assert FAST.trace(c) == 1
            continue
        solved += 1
        assert FAST.mul(t, t) ^ t == c
        # the other root is t + 1, and there are exactly two
        assert FAST.mul(t ^ 1, t ^ 1) ^ (t ^ 1) == c
    assert solved > 0, "no trace-zero c sampled"


def test_inv_matches_reference():
    rng = random.Random(5)
    for _ in range(80):
        a = _rand(rng) or 1
        assert FAST.inv(a) == REF.inv(a)
        assert FAST.mul(a, FAST.inv(a)) == 1


def test_batch_inv_matches_scalar_inv():
    rng = random.Random(6)
    xs = [_rand(rng) or 1 for _ in range(50)]
    got = FAST.batch_inv(xs)
    for x, g in zip(xs, got):
        assert FAST.mul(x, g) == 1


def test_batch_inv_passes_zero_through():
    rng = random.Random(7)
    xs = [_rand(rng) or 1 for _ in range(20)]
    xs[5] = 0
    xs[19] = 0
    got = FAST.batch_inv(xs)
    assert got[5] == 0 and got[19] == 0
    for i, (x, g) in enumerate(zip(xs, got)):
        if x:
            assert FAST.mul(x, g) == 1


def test_frobenius_is_iterated_squaring():
    rng = random.Random(8)
    for _ in range(40):
        a = _rand(rng)
        acc = a
        for k in range(1, 12):
            acc = FAST.sqr(acc)
            assert FAST.frobenius(a, k) == acc


def test_frobenius_has_order_deg():
    rng = random.Random(9)
    for _ in range(20):
        a = _rand(rng)
        assert FAST.frobenius(a, N131) == a
        assert FAST.frobenius(a, 0) == a


def test_sqrt_via_frobenius_inverts_sqr():
    rng = random.Random(10)
    for _ in range(40):
        a = _rand(rng)
        assert FAST.frobenius(FAST.sqr(a), N131 - 1) == a


if __name__ == "__main__":
    import traceback
    fails = 0
    for name, fn in sorted(globals().items()):
        if not name.startswith("test_") or not callable(fn):
            continue
        try:
            fn()
            print(f"ok   {name}")
        except Exception:
            fails += 1
            print(f"FAIL {name}")
            traceback.print_exc()
    print("failures:", fails)
    sys.exit(1 if fails else 0)
