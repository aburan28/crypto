"""
Batched r-adding Pollard rho walk with distinguished points, on device.

State per walk (all little-endian limbs, int8):
  x, y    affine coordinates of the current point W = c*P + d*Q      [N, L]
  c, d    coefficients mod n                                           [N, Ln]
  frozen  walk has reached a distinguished point (or a bad state) and is
          parked until the host harvests and reseeds it                [N] bool
  bad     the parked state is garbage (x == table x: doubling/inverse case),
          so the host must reseed without recording a DP               [N] bool

Walk function: j = x[0] & (r-1) selects a table point R_j = a_j P + b_j Q, and
W <- W + R_j, (c, d) <- (c + a_j, d + b_j) mod n. Distinguished point test
reads bits above the index bits so the two decisions are independent.

One jitted call advances every walk by `stepsPerChunk` steps. The walk is an
affine add with a tree batch inversion across the whole batch, so each step
costs ~6 batched modmuls plus one serial Fermat inversion at the tree root.
"""

import numpy as np
import jax
import jax.numpy as jnp
from jax import lax

from rho import field

DP_LIMB_START = 1          # distinguished-point bits begin at limb 1 (above index bits)
DP_LIMB_COUNT = 3          # up to 21 DP bits read from limbs 1..3


class WalkContext:
    def __init__(self, curve, tableAx, tableBx, tableX, tableY, tableSize, dpBits):
        assert tableSize & (tableSize - 1) == 0 and tableSize <= field.BASE
        assert dpBits <= field.LIMB_BITS * DP_LIMB_COUNT
        self.fieldCtx = field.FieldContext(curve.p)
        self.orderCtx = field.FieldContext(curve.order)
        L = self.fieldCtx.numLimbs
        Ln = self.orderCtx.numLimbs
        assert L >= DP_LIMB_START + DP_LIMB_COUNT, "field too small for the DP window"
        self.numLimbs = L
        self.orderLimbs = Ln
        self.tableSize = tableSize
        self.dpBits = dpBits
        self.dpMask = (1 << dpBits) - 1
        self.tableX = jnp.asarray(field.toLimbsBatch(tableX, L))
        self.tableY = jnp.asarray(field.toLimbsBatch(tableY, L))
        self.tableA = jnp.asarray(field.toLimbsBatch(tableAx, Ln))
        self.tableB = jnp.asarray(field.toLimbsBatch(tableBx, Ln))


def tableIndex(ctx, x):
    return (x[:, 0].astype(jnp.int32)) & (ctx.tableSize - 1)


def dpValue(ctx, x):
    value = jnp.zeros(x.shape[0], jnp.int32)
    for i in range(DP_LIMB_COUNT):
        value = value | (x[:, DP_LIMB_START + i].astype(jnp.int32) << (field.LIMB_BITS * i))
    return value


def isDistinguished(ctx, x):
    return (dpValue(ctx, x) & ctx.dpMask) == 0


def walkStep(ctx, x, y, c, d):
    """One affine r-adding step for every walk. Returns (x, y, c, d, bad)."""
    fc = ctx.fieldCtx
    oc = ctx.orderCtx
    j = tableIndex(ctx, x)
    rx = ctx.tableX[j]
    ry = ctx.tableY[j]

    den = field.subMod(fc, rx, x)
    bad = field.isZero(den)
    den = jnp.where(bad[:, None], jnp.broadcast_to(fc.one, den.shape), den)
    invDen = field.batchInverse(fc, den)

    num = field.subMod(fc, ry, y)
    lam = field.mulMod(fc, num, invDen)
    lamSq = field.squareMod(fc, lam)
    x3 = field.subMod(fc, field.subMod(fc, lamSq, x), rx)
    y3 = field.subMod(fc, field.mulMod(fc, lam, field.subMod(fc, x, x3)), y)

    c3 = field.addMod(oc, c, ctx.tableA[j])
    d3 = field.addMod(oc, d, ctx.tableB[j])
    return x3, y3, c3, d3, bad


def makeChunkRunner(ctx, stepsPerChunk):
    """Build a jitted function advancing the whole batch by stepsPerChunk steps.

    Walks that hit a distinguished point (or a bad denominator) freeze in place
    for the remainder of the chunk; the host harvests them afterwards.
    """

    def body(i, state):
        x, y, c, d, frozen, bad = state
        nx, ny, nc, nd, stepBad = walkStep(ctx, x, y, c, d)
        active = ~frozen
        keep = active[:, None]
        x = jnp.where(keep, nx, x)
        y = jnp.where(keep, ny, y)
        c = jnp.where(keep, nc, c)
        d = jnp.where(keep, nd, d)
        newBad = active & stepBad
        hit = active & (isDistinguished(ctx, nx) | stepBad)
        return x, y, c, d, frozen | hit, bad | newBad

    def run(x, y, c, d, frozen, bad):
        return lax.fori_loop(0, stepsPerChunk, body, (x, y, c, d, frozen, bad))

    return jax.jit(run)
