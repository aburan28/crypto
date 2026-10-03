"""
Multi-limb GF(p) arithmetic for TPU, built around int8 matmuls.

Representation
--------------
A field element is a little-endian vector of LIMB_BITS-bit digits stored in int8.
LIMB_BITS = 7 so that every digit (0..127) is representable as a *signed* int8,
which is what the TPU MXU consumes natively (int8 x int8 -> int32 accumulate).
Intermediate (unnormalised) results live in int32.

Where the MXU earns its keep
----------------------------
A modular multiply a*b mod p is three big-integer products:

  1. a*b           walk-dependent x walk-dependent   -> batched Toeplitz matvec
  2. q1*mu         walk-dependent x constant         -> dense matmul, shared weight
  3. q3*p          walk-dependent x constant         -> dense matmul, shared weight

(2) and (3) are the Barrett reduction. Their right-hand operand is a fixed
Toeplitz matrix built once from mu = floor(B^(2L)/p) and from p, so across a
batch of N walks they are ordinary [N, L+1] @ [L+1, M] int8 matmuls: exactly the
shape the systolic array wants. (1) has no shared operand (every walk multiplies
two different numbers) so it is expressed as a batched Toeplitz product; XLA
lowers it to a batched matmul with N=1 on the inner dimension, which is far less
MXU-efficient. Two of the three big products per modmul are therefore dense
weight-shared MXU work; the third is the known weak spot and the obvious place
for a Pallas kernel later.

Carry resolution is a lax.scan over the limb axis (L sequential vector ops over
the batch). It is exact for arbitrary inputs, including negatives (two's
complement arithmetic shift gives floor division, so a-b normalises to
(a-b) mod B^L with the borrow surfacing as a -1 top carry).
"""

import numpy as np
import jax
import jax.numpy as jnp
from jax import lax

LIMB_BITS = 7
BASE = 1 << LIMB_BITS
MASK = BASE - 1


# ----------------------------------------------------------------------------
# host <-> limb conversion
# ----------------------------------------------------------------------------

def limbCount(modulus):
    return (modulus.bit_length() + LIMB_BITS - 1) // LIMB_BITS


def toLimbs(value, numLimbs):
    out = np.zeros(numLimbs, np.int8)
    for i in range(numLimbs):
        out[i] = value & MASK
        value >>= LIMB_BITS
    assert value == 0, "value does not fit in the given number of limbs"
    return out


def toLimbsBatch(values, numLimbs):
    out = np.zeros((len(values), numLimbs), np.int8)
    for row, value in enumerate(values):
        out[row] = toLimbs(value, numLimbs)
    return out


def fromLimbs(limbs):
    limbs = np.asarray(limbs)
    value = 0
    for i in range(limbs.shape[-1] - 1, -1, -1):
        value = (value << LIMB_BITS) | int(limbs[i])
    return value


def fromLimbsBatch(limbs):
    return [fromLimbs(row) for row in np.asarray(limbs)]


# ----------------------------------------------------------------------------
# field context: all the constant matrices for one modulus
# ----------------------------------------------------------------------------

class FieldContext:
    """Precomputed constants for Barrett reduction modulo p in base 2^LIMB_BITS."""

    def __init__(self, modulus):
        p = int(modulus)
        L = limbCount(p)
        self.modulus = p
        self.numLimbs = L
        self.bitLength = p.bit_length()

        # Barrett constant mu = floor(B^(2L) / p); fits in L+1 limbs.
        mu = (BASE ** (2 * L)) // p
        muLimbs = toLimbs(mu, L + 1)
        pLimbs = toLimbs(p, L)

        # Toeplitz matrix for q1*mu : [L+1, 2L+2], entry (i,k) = mu[k-i]
        tMu = np.zeros((L + 1, 2 * L + 2), np.int8)
        for i in range(L + 1):
            for j in range(L + 1):
                tMu[i, i + j] = muLimbs[j]
        # Truncated Toeplitz for (q3*p) mod B^(L+1): [L+1, L+1], entry (i,k) = p[k-i]
        tP = np.zeros((L + 1, L + 1), np.int8)
        for i in range(L + 1):
            for j in range(L):
                if i + j < L + 1:
                    tP[i, i + j] = pLimbs[j]

        self.tMu = jnp.asarray(tMu)
        self.tP = jnp.asarray(tP)
        self.pLimbs = jnp.asarray(pLimbs)
        self.pLimbsWide = jnp.asarray(toLimbs(p, L + 1))
        self.one = jnp.asarray(toLimbs(1, L))
        self.zero = jnp.zeros(L, jnp.int8)

        # Toeplitz gather indices for the walk-dependent product a*b.
        kIdx = np.arange(2 * L - 1)[:, None]
        jIdx = np.arange(L)[None, :]
        diff = kIdx - jIdx
        self.toeplitzMask = jnp.asarray(((diff >= 0) & (diff < L)).astype(np.int8))
        self.toeplitzIdx = jnp.asarray(np.clip(diff, 0, L - 1))

        # Exponent bits of p-2 (MSB first) for Fermat inversion.
        exponent = p - 2
        bits = [(exponent >> i) & 1 for i in range(exponent.bit_length() - 1, -1, -1)]
        self.inverseExponentBits = jnp.asarray(np.array(bits, np.int32))


# ----------------------------------------------------------------------------
# core limb ops
# ----------------------------------------------------------------------------

def normalize(values):
    """Propagate carries along the last axis. Returns (digits int32, finalCarry int32).

    Exact for any int32 input, positive or negative. The final carry is the
    quotient on division by B^n (so -1 signals a borrow out the top).
    """
    values = values.astype(jnp.int32)
    swapped = jnp.moveaxis(values, -1, 0)

    def body(carry, limb):
        total = limb + carry
        return total >> LIMB_BITS, total & MASK

    finalCarry, digits = lax.scan(body, jnp.zeros(values.shape[:-1], jnp.int32), swapped)
    return jnp.moveaxis(digits, 0, -1), finalCarry


def padLimbs(values, width):
    extra = width - values.shape[-1]
    if extra <= 0:
        return values[..., :width]
    padSpec = [(0, 0)] * (values.ndim - 1) + [(0, extra)]
    return jnp.pad(values, padSpec)


def toeplitzProduct(ctx, a, b):
    """Full product of two L-limb values, unnormalised, shape [..., 2L-1] int32.

    T[b, k, j] = a[b, k-j]; c = T @ b. Both operands vary per batch element so
    this is a batched matvec (the MXU-unfriendly one, see module docstring).
    """
    a8 = a.astype(jnp.int8)
    b8 = b.astype(jnp.int8)
    toeplitz = jnp.take(a8, ctx.toeplitzIdx, axis=-1) * ctx.toeplitzMask  # [..., 2L-1, L]
    batchDims = tuple(range(a.ndim - 1))
    return lax.dot_general(
        toeplitz, b8,
        ((( toeplitz.ndim - 1,), (b8.ndim - 1,)), (batchDims, batchDims)),
        preferred_element_type=jnp.int32,
    )


def constantProduct(values, matrix):
    """[..., K] int8 @ [K, M] int8 -> [..., M] int32. Shared-weight dense matmul."""
    return lax.dot_general(
        values.astype(jnp.int8), matrix,
        (((values.ndim - 1,), (0,)), ((), ())),
        preferred_element_type=jnp.int32,
    )


def conditionalSubtract(values, modulusLimbs):
    """values - modulus if values >= modulus else values. Same width as values."""
    diff, borrow = normalize(values - modulusLimbs.astype(jnp.int32))
    keep = (borrow < 0)[..., None]
    return jnp.where(keep, values, diff)


def barrettReduce(ctx, productDigits):
    """Reduce a normalised 2L-digit product (< p^2) to L digits in [0, p)."""
    L = ctx.numLimbs
    q1 = productDigits[..., L - 1:]                       # L+1 digits
    q2, _ = normalize(constantProduct(q1, ctx.tMu))       # 2L+2 digits, MXU
    q3 = q2[..., L + 1:]                                  # L+1 digits
    r1 = productDigits[..., :L + 1]
    r2, _ = normalize(constantProduct(q3, ctx.tP))        # L+1 digits mod B^(L+1), MXU
    r, _ = normalize(r1 - r2)                             # mod B^(L+1); true value in [0, 3p)
    r = conditionalSubtract(r, ctx.pLimbsWide)
    r = conditionalSubtract(r, ctx.pLimbsWide)
    return r[..., :L].astype(jnp.int8)


def mulMod(ctx, a, b):
    L = ctx.numLimbs
    product = toeplitzProduct(ctx, a, b)
    digits, _ = normalize(padLimbs(product, 2 * L))
    return barrettReduce(ctx, digits)


def squareMod(ctx, a):
    return mulMod(ctx, a, a)


def addMod(ctx, a, b):
    # a+b < 2p < 2*B^L: a carry out of the top limb is possible, keep an extra digit
    wide, carry = normalize(a.astype(jnp.int32) + b.astype(jnp.int32))
    wide = jnp.concatenate([wide, carry[..., None]], axis=-1)
    wide = conditionalSubtract(wide, ctx.pLimbsWide)
    return wide[..., :ctx.numLimbs].astype(jnp.int8)


def subMod(ctx, a, b):
    diff, borrow = normalize(a.astype(jnp.int32) - b.astype(jnp.int32))
    fixed, _ = normalize(diff + ctx.pLimbs.astype(jnp.int32))
    return jnp.where((borrow < 0)[..., None], fixed, diff).astype(jnp.int8)


def negMod(ctx, a):
    return subMod(ctx, jnp.broadcast_to(ctx.zero, a.shape), a)


def isZero(a):
    return jnp.all(a == 0, axis=-1)


def inverseMod(ctx, a):
    """a^(p-2) by left-to-right binary exponentiation. Serial in the exponent bits."""
    bits = ctx.inverseExponentBits
    result = jnp.broadcast_to(ctx.one, a.shape)

    def body(i, result):
        result = squareMod(ctx, result)
        multiplied = mulMod(ctx, result, a)
        return jnp.where(bits[i] == 1, multiplied, result)

    return lax.fori_loop(0, bits.shape[0], body, result)


def batchInverse(ctx, a):
    """Invert a batch [N, L] (N a power of two) with a product tree.

    3N modmuls arranged in 2*log2(N) batched levels, plus one serial Fermat
    inversion at the root. Any zero in the batch poisons the whole tree, so
    callers must substitute a non-zero value for zero inputs first.
    """
    n = a.shape[0]
    assert n & (n - 1) == 0, "batch size must be a power of two"
    levels = [a]
    current = a
    while current.shape[0] > 1:
        current = mulMod(ctx, current[0::2], current[1::2])
        levels.append(current)
    inv = inverseMod(ctx, levels[-1])
    for level in reversed(levels[:-1]):
        left = level[0::2]
        right = level[1::2]
        invLeft = mulMod(ctx, inv, right)
        invRight = mulMod(ctx, inv, left)
        inv = jnp.stack([invLeft, invRight], axis=1).reshape(level.shape)
    return inv
