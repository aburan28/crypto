"""Batched GF(2^m) arithmetic for the TPU, in the bit-matmul formulation.

A field element is a length-``n`` vector of 0/1 ``int32`` coefficients
(coefficient of ``t^i`` at index ``i``), and a batch is ``[B, n]``.  This
is the representation a TPU can actually work in, and it is chosen on
purpose:

* **Addition** is ``XOR`` -- elementwise.
* **Squaring** and **reduction** are ``F_2``-linear, so each is a fixed
  0/1 matrix and the operation is an integer matmul taken ``mod 2`` -- the
  systolic array's native shape.
* **Multiplication** is ``F_2``-bilinear: ``mul(a,b)_k = sum_{i,j} a_i b_j
  T[i,j,k] (mod 2)`` for the field's fixed multiplication tensor ``T``.
  As ``(a outer b) . T`` that is one contraction, i.e. a matmul of the
  ``n^2``-wide outer product against a fixed ``[n^2, n]`` matrix.
* **Inversion** is Fermat ``a^(2^n-2)``: a static square-and-multiply
  chain, hence a fixed sequence of the matmuls above -- which is exactly
  why the batch inversion that dominates the pair-table build becomes a
  batched bit-matmul on a TPU.

Every constant is built from :mod:`tpu.ic.reference`, the scalar oracle,
so this file has no independent notion of the field.  Correctness is "it
agrees with the oracle", checked in ``tpu/tests/test_field.py``.

NOTE ON FIT.  This is faithful, not fast: a single GF(2^m) multiply here
is an ``O(n^2)`` contraction where a CPU with a carry-less-multiply
instruction does it in a few cycles.  The win a TPU could offer is only
in *batch* -- thousands of independent field ops issued as a few large
matmuls -- and whether that batch win ever repays the ``n^2`` blow-up is a
device measurement this repository does not yet have.  See
``tpu/protocol/RESEARCH_TPU_IC.md``.
"""

from __future__ import annotations

from functools import partial

import jax
import jax.numpy as jnp
import numpy as np

from .reference import Field


def int_to_bits(v: int, n: int) -> np.ndarray:
    return np.array([(v >> i) & 1 for i in range(n)], dtype=np.int32)


def bits_to_int(bits) -> int:
    v = 0
    for i, b in enumerate(np.asarray(bits).tolist()):
        if int(b) & 1:
            v |= 1 << i
    return v


class FieldJax:
    """JAX-side constants for a :class:`tpu.ic.reference.Field`."""

    def __init__(self, field: Field):
        self.field = field
        n = field.n
        self.n = n
        # squaring matrix: sqr_mat[k, j] = bit k of sqr(t^j)
        sqr_mat = np.zeros((n, n), dtype=np.int32)
        for j in range(n):
            col = int_to_bits(field.sqr(1 << j), n)
            sqr_mat[:, j] = col
        self.sqr_mat = jnp.asarray(sqr_mat)
        # multiplication tensor: T[i, j, k] = bit k of mul(t^i, t^j)
        mul_t = np.zeros((n, n, n), dtype=np.int32)
        for i in range(n):
            for j in range(n):
                mul_t[i, j, :] = int_to_bits(field.mul(1 << i, 1 << j), n)
        self.mul_tensor = jnp.asarray(mul_t)
        # flattened [n*n, n] form, for the matmul statement of a multiply
        self.mul_mat = jnp.asarray(mul_t.reshape(n * n, n))
        self.one = jnp.asarray(int_to_bits(1, n))
        # Fermat exponent bits, low to high (static; drives the inverse loop)
        self._inv_bits = [((1 << n) - 2) >> i & 1 for i in range(n)]

    # -- elementwise ------------------------------------------------------
    @staticmethod
    def add(a, b):
        return jnp.bitwise_xor(a, b)

    # -- linear maps (matmul mod 2) --------------------------------------
    def sqr(self, a):
        # a: [..., n]  ->  (a @ sqr_mat^T) mod 2
        return jnp.mod(a @ self.sqr_mat.T, 2).astype(jnp.int32)

    def sqr_k(self, a, k: int):
        for _ in range(int(k)):
            a = self.sqr(a)
        return a

    # -- bilinear map (one contraction mod 2) ----------------------------
    def mul(self, a, b):
        # outer product a_i b_j, flattened, matmul against mul_mat.
        outer = (a[..., :, None] * b[..., None, :]).reshape(a.shape[:-1] + (self.n * self.n,))
        return jnp.mod(outer @ self.mul_mat, 2).astype(jnp.int32)

    def inv(self, a):
        # Fermat a^(2^n - 2) by static square-and-multiply; inv(0)=0 falls
        # out because the chain multiplies by a square of 0.
        one = jnp.broadcast_to(self.one, a.shape)
        r = one
        base = a
        for bit in self._inv_bits:
            if bit:
                r = self.mul(r, base)
            base = self.sqr(base)
        return r
