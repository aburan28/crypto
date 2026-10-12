"""
Host driver for parallel collision search (van Oorschot-Wiener) on ECDLP.

  seed N walks at random c*P + d*Q
  loop:
    advance every walk K steps on device
    harvest frozen walks:  record (x -> c, d), reseed
    on a repeated x with different coefficients -> solve for k with Q = kP
"""

import random
import time

import numpy as np
import jax.numpy as jnp

from rho import curve as curveLib
from rho import field
from rho import walk as walkLib


def buildTable(curve, Q, tableSize, rng):
    n = curve.order
    P = curve.generator
    tableAx = []
    tableBx = []
    tableX = []
    tableY = []
    for _ in range(tableSize):
        while True:
            aj = rng.randrange(1, n)
            bj = rng.randrange(1, n)
            point = curve.add(curve.mul(aj, P), curve.mul(bj, Q))
            if point is not None:
                break
        tableAx.append(aj)
        tableBx.append(bj)
        tableX.append(point[0])
        tableY.append(point[1])
    return tableAx, tableBx, tableX, tableY


def hostWalkStep(curve, ctx, point, c, d, tableAx, tableBx, tableX, tableY):
    """Reference implementation of exactly the device walk function."""
    x, y = point
    j = x & (ctx.tableSize - 1)
    nextPoint = curve.add(point, (tableX[j], tableY[j]))
    n = curve.order
    return nextPoint, (c + tableAx[j]) % n, (d + tableBx[j]) % n


def hostIsDistinguished(ctx, x):
    value = (x >> field.LIMB_BITS) & ((1 << (field.LIMB_BITS * walkLib.DP_LIMB_COUNT)) - 1)
    return (value & ctx.dpMask) == 0


def seedWalks(curve, Q, count, rng):
    n = curve.order
    P = curve.generator
    cs = []
    ds = []
    xs = []
    ys = []
    for _ in range(count):
        while True:
            c = rng.randrange(1, n)
            d = rng.randrange(1, n)
            point = curve.add(curve.mul(c, P), curve.mul(d, Q))
            if point is not None:
                break
        cs.append(c)
        ds.append(d)
        xs.append(point[0])
        ys.append(point[1])
    return cs, ds, xs, ys


def solveFromCollision(curve, Q, c1, d1, c2, d2):
    """Two walks landed on the same x: c1 P + d1 Q = +/- (c2 P + d2 Q)."""
    n = curve.order
    P = curve.generator
    candidates = []
    if (d2 - d1) % n != 0:
        candidates.append((c1 - c2) * pow(d2 - d1, -1, n) % n)
    if (d1 + d2) % n != 0:
        candidates.append((-(c1 + c2)) * pow(d1 + d2, -1, n) % n)
    for k in candidates:
        if curve.mul(k, P) == Q:
            return k
    return None


class RhoSolver:
    def __init__(self, curve, Q, numWalks=2048, stepsPerChunk=32, tableSize=32, dpBits=8, seed=0, verbose=True):
        assert numWalks & (numWalks - 1) == 0, "numWalks must be a power of two"
        self.curve = curve
        self.Q = Q
        self.numWalks = numWalks
        self.stepsPerChunk = stepsPerChunk
        self.rng = random.Random(seed)
        self.verbose = verbose
        self.table = buildTable(curve, Q, tableSize, self.rng)
        self.ctx = walkLib.WalkContext(curve, *self.table, tableSize=tableSize, dpBits=dpBits)
        self.runChunk = walkLib.makeChunkRunner(self.ctx, stepsPerChunk)
        self.dpStore = {}
        self.stepsIssued = 0
        self.dpCount = 0
        self.badCount = 0

    def log(self, message):
        if self.verbose:
            print(message, flush=True)

    def solve(self, maxChunks=None):
        curve = self.curve
        L = self.ctx.numLimbs
        Ln = self.ctx.orderLimbs
        cs, ds, xs, ys = seedWalks(curve, self.Q, self.numWalks, self.rng)
        x = jnp.asarray(field.toLimbsBatch(xs, L))
        y = jnp.asarray(field.toLimbsBatch(ys, L))
        c = jnp.asarray(field.toLimbsBatch(cs, Ln))
        d = jnp.asarray(field.toLimbsBatch(ds, Ln))
        frozen = jnp.zeros(self.numWalks, bool)
        bad = jnp.zeros(self.numWalks, bool)

        expected = (3.14159 * curve.order / 2) ** 0.5
        self.log(f"group order {curve.order.bit_length()} bits, expected ~{expected:.3e} steps, "
                 f"{self.numWalks} walks x {self.stepsPerChunk} steps/chunk, dpBits={self.ctx.dpBits}")
        startTime = time.time()
        chunk = 0
        while maxChunks is None or chunk < maxChunks:
            chunk += 1
            x, y, c, d, frozen, bad = self.runChunk(x, y, c, d, frozen, bad)
            self.stepsIssued += self.numWalks * self.stepsPerChunk
            frozenHost = np.asarray(frozen)
            if not frozenHost.any():
                continue
            badHost = np.asarray(bad)
            idx = np.nonzero(frozenHost)[0]
            xHost = np.asarray(x)[idx]
            cHost = np.asarray(c)[idx]
            dHost = np.asarray(d)[idx]

            for row, walkId in enumerate(idx):
                if badHost[walkId]:
                    self.badCount += 1
                    continue
                xVal = field.fromLimbs(xHost[row])
                cVal = field.fromLimbs(cHost[row])
                dVal = field.fromLimbs(dHost[row])
                self.dpCount += 1
                if xVal in self.dpStore:
                    c2, d2 = self.dpStore[xVal]
                    if (cVal, dVal) != (c2, d2):
                        k = solveFromCollision(curve, self.Q, cVal, dVal, c2, d2)
                        if k is not None:
                            elapsed = time.time() - startTime
                            self.log(f"collision after {self.stepsIssued:.3e} issued steps, "
                                     f"{self.dpCount} DPs, {elapsed:.1f}s ({self.stepsIssued/elapsed:.3e} steps/s)")
                            return k
                else:
                    self.dpStore[xVal] = (cVal, dVal)

            # reseed everything that was frozen
            newC, newD, newX, newY = seedWalks(curve, self.Q, len(idx), self.rng)
            x = x.at[idx].set(jnp.asarray(field.toLimbsBatch(newX, L)))
            y = y.at[idx].set(jnp.asarray(field.toLimbsBatch(newY, L)))
            c = c.at[idx].set(jnp.asarray(field.toLimbsBatch(newC, Ln)))
            d = d.at[idx].set(jnp.asarray(field.toLimbsBatch(newD, Ln)))
            frozen = frozen.at[idx].set(False)
            bad = bad.at[idx].set(False)

            if chunk % 10 == 0:
                elapsed = time.time() - startTime
                self.log(f"chunk {chunk}: {self.stepsIssued:.3e} steps, {self.dpCount} DPs, "
                         f"{self.stepsIssued/elapsed:.3e} steps/s")
        return None
