"""End-to-end: generate a small prime-order curve, solve a random ECDLP with the
batched rho harness, and check the recovered scalar. Also exercises the
prime-order curve generator and the collision solver on the +/- branches."""

import random

from rho import curve as curveLib
from rho import solve


def test_generate_prime_order_curve():
    rng = random.Random(7)
    curve = curveLib.generatePrimeOrderCurve(36, rng)
    assert curveLib.isProbablePrime(curve.order)
    assert curve.isOnCurve(curve.generator)
    assert curve.mul(curve.order, curve.generator) is None
    # Hasse bound
    import math
    assert abs(curve.order - (curve.p + 1)) <= 2 * math.isqrt(curve.p) + 1


def test_solve_from_collision_both_signs():
    rng = random.Random(8)
    curve = curveLib.generatePrimeOrderCurve(28, rng)
    n = curve.order
    k = rng.randrange(1, n)
    Q = curve.mul(k, curve.generator)
    c1, d1 = rng.randrange(1, n), rng.randrange(1, n)
    # same point: c2 P + d2 Q == c1 P + d1 Q with d2 != d1
    d2 = (d1 + 12345) % n
    c2 = (c1 - 12345 * k) % n
    assert solve.solveFromCollision(curve, Q, c1, d1, c2, d2) == k
    # negated point: c3 P + d3 Q == -(c1 P + d1 Q)
    d3 = (-(d1) + 777) % n
    c3 = (-(c1) - 777 * k) % n
    assert solve.solveFromCollision(curve, Q, c1, d1, c3, d3) == k


def test_end_to_end_solve_small_curve():
    rng = random.Random(2026)
    curve = curveLib.generatePrimeOrderCurve(32, rng)
    secret = rng.randrange(1, curve.order)
    Q = curve.mul(secret, curve.generator)
    solver = solve.RhoSolver(curve, Q, numWalks=512, stepsPerChunk=16, tableSize=32, dpBits=6, seed=1, verbose=False)
    k = solver.solve(maxChunks=400)
    assert k == secret
    assert curve.mul(k, curve.generator) == Q
