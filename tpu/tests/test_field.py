"""JAX field arithmetic agrees with the scalar oracle, element by element."""

import random

import numpy as np
import pytest

from ic.field import FieldJax, bits_to_int, int_to_bits
from ic.instances import find_field


def _batch(vals, n):
    return np.stack([int_to_bits(v, n) for v in vals]).astype(np.int32)


def _ints(bits):
    bits = np.asarray(bits)
    return [bits_to_int(bits[i]) for i in range(bits.shape[0])]


@pytest.mark.parametrize("n", [7, 11, 13, 17])
def test_mul_sqr_inv_match_oracle(n):
    f = find_field(n)
    fj = FieldJax(f)
    rng = random.Random(0xC0FFEE ^ n)
    xs = [rng.randrange(1 << n) for _ in range(64)]
    ys = [rng.randrange(1 << n) for _ in range(64)]

    A = _batch(xs, n)
    B = _batch(ys, n)

    mul = _ints(fj.mul(A, B))
    sqr = _ints(fj.sqr(A))
    inv = _ints(fj.inv(A))

    for i, (x, y) in enumerate(zip(xs, ys)):
        assert mul[i] == f.mul(x, y), (n, x, y)
        assert sqr[i] == f.sqr(x), (n, x)
        assert inv[i] == f.inv(x), (n, x)


@pytest.mark.parametrize("n", [7, 11, 13])
def test_inverse_is_inverse(n):
    f = find_field(n)
    fj = FieldJax(f)
    xs = [v for v in range(1, min(1 << n, 200))]
    A = _batch(xs, n)
    prod = _ints(fj.mul(A, fj.inv(A)))
    assert all(p == 1 for p in prod)


def test_add_is_xor():
    f = find_field(11)
    fj = FieldJax(f)
    A = _batch([1, 2, 3, 255], 11)
    B = _batch([1, 3, 3, 1], 11)
    out = _ints(FieldJax.add(A, B))
    assert out == [0, 1, 0, 254]
