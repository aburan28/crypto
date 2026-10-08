"""The device walk step on secp256k1 is the host reference walk, step for step,
and the solver invariant W = c*P + d*Q holds after every step."""

import random

import jax
import jax.numpy as jnp
import numpy as np

from rho import curve as curveLib
from rho import field
from rho import solve
from rho import walk as walkLib


def test_secp256k1_walk_matches_host():
    curve = curveLib.SECP256K1
    rng = random.Random(3)
    k = rng.randrange(1, curve.order)
    Q = curve.mul(k, curve.generator)
    table = solve.buildTable(curve, Q, 32, rng)
    ctx = walkLib.WalkContext(curve, *table, tableSize=32, dpBits=21)
    L = ctx.numLimbs
    Ln = ctx.orderLimbs
    numWalks = 16
    steps = 3

    cs, ds, xs, ys = solve.seedWalks(curve, Q, numWalks, rng)
    x = jnp.asarray(field.toLimbsBatch(xs, L))
    y = jnp.asarray(field.toLimbsBatch(ys, L))
    c = jnp.asarray(field.toLimbsBatch(cs, Ln))
    d = jnp.asarray(field.toLimbsBatch(ds, Ln))
    step = jax.jit(lambda x, y, c, d: walkLib.walkStep(ctx, x, y, c, d))

    hostPoints = list(zip(xs, ys))
    hostC = list(cs)
    hostD = list(ds)
    for _ in range(steps):
        x, y, c, d, bad = step(x, y, c, d)
        assert not bool(np.asarray(bad).any())
        gotX = field.fromLimbsBatch(x)
        gotY = field.fromLimbsBatch(y)
        gotC = field.fromLimbsBatch(c)
        gotD = field.fromLimbsBatch(d)
        for i in range(numWalks):
            hostPoints[i], hostC[i], hostD[i] = solve.hostWalkStep(
                curve, ctx, hostPoints[i], hostC[i], hostD[i], *table)
            assert hostPoints[i] == (gotX[i], gotY[i])
            assert (hostC[i], hostD[i]) == (gotC[i], gotD[i])
            assert curve.add(curve.mul(gotC[i], curve.generator), curve.mul(gotD[i], Q)) == hostPoints[i]


def test_distinguished_point_test_matches_host():
    curve = curveLib.SECP256K1
    rng = random.Random(4)
    table = solve.buildTable(curve, curve.generator, 32, rng)
    ctx = walkLib.WalkContext(curve, *table, tableSize=32, dpBits=5)
    xs = [rng.randrange(curve.p) for _ in range(512)]
    device = np.asarray(walkLib.isDistinguished(ctx, jnp.asarray(field.toLimbsBatch(xs, ctx.numLimbs))))
    host = [solve.hostIsDistinguished(ctx, xv) for xv in xs]
    assert list(map(bool, device)) == host
    assert any(host), "expected some DPs at 5 bits over 512 samples"
