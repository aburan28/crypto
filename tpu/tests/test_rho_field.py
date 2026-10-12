"""Limb GF(p) arithmetic agrees with Python integers, including batch inversion."""

import random

import jax
import jax.numpy as jnp
import pytest

from rho import field

SECP_P = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F


@pytest.mark.parametrize("p,rounds", [
    ((1 << 40) - 87, 256),
    (0xFFFFFFFFFFFFFFC5, 128),
    (SECP_P, 64),
])
def test_mul_add_sub_match_python(p, rounds):
    rng = random.Random(1)
    ctx = field.FieldContext(p)
    L = ctx.numLimbs
    aVals = [0, 1, p - 1, p - 2] + [rng.randrange(p) for _ in range(rounds)]
    bVals = [0, p - 1, p - 1, 2] + [rng.randrange(p) for _ in range(rounds)]
    a = jnp.asarray(field.toLimbsBatch(aVals, L))
    b = jnp.asarray(field.toLimbsBatch(bVals, L))

    gotMul = field.fromLimbsBatch(jax.jit(lambda x, y: field.mulMod(ctx, x, y))(a, b))
    gotAdd = field.fromLimbsBatch(jax.jit(lambda x, y: field.addMod(ctx, x, y))(a, b))
    gotSub = field.fromLimbsBatch(jax.jit(lambda x, y: field.subMod(ctx, x, y))(a, b))
    for i in range(len(aVals)):
        assert gotMul[i] == aVals[i] * bVals[i] % p
        assert gotAdd[i] == (aVals[i] + bVals[i]) % p
        assert gotSub[i] == (aVals[i] - bVals[i]) % p


@pytest.mark.parametrize("p", [(1 << 40) - 87, SECP_P])
def test_batch_inverse_matches_python(p):
    rng = random.Random(2)
    ctx = field.FieldContext(p)
    L = ctx.numLimbs
    vals = [rng.randrange(1, p) for _ in range(64)]
    inv = jax.jit(lambda x: field.batchInverse(ctx, x))(jnp.asarray(field.toLimbsBatch(vals, L)))
    got = field.fromLimbsBatch(inv)
    for i in range(64):
        assert got[i] == pow(vals[i], -1, p)


def test_normalize_handles_negative_borrow():
    # (a - b) with a < b must wrap mod B^L and surface a -1 top carry.
    digits, carry = field.normalize(jnp.asarray([[0, 0, 0]]) - jnp.asarray([[1, 0, 0]]))
    assert int(carry[0]) == -1
    assert field.fromLimbs(digits[0]) == field.BASE ** 3 - 1
