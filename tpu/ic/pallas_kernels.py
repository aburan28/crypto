"""Pallas kernels for the index-calculus TPU backend.

The one operation under everything in this backend is **bit-matrix
multiply mod 2**: GF(2^m) squaring and reduction are such a product
(:mod:`tpu.ic.field`), the GF(2^m) multiply is a pair of them, and the
GF(2) relation algebra is one directly (:mod:`tpu.ic.linalg`).  So there
is one kernel worth writing by hand, and it is here.

:func:`bitmatmul_mod2` tiles the output over a grid and, per tile,
accumulates an integer matmul and takes it ``mod 2``.  On a TPU
(``interpret=False`` on a TPU host) the inner ``jnp.dot`` lowers to the
MXU; here, with ``interpret=True`` (the default), Pallas runs the same
kernel on the CPU so it can be checked without a device -- the same
host-verification discipline ``gpu/ecc2k`` uses for its CUDA.

Integer accumulation stays exact: entries are 0/1 and a tile contracts at
most ``block_k`` of them, so the partial sums are well within ``int32``
before the final ``mod 2``.
"""

from __future__ import annotations

import jax
import jax.numpy as jnp

try:  # Pallas moved modules across JAX versions; tolerate both.
    from jax.experimental import pallas as pl
    _HAVE_PALLAS = True
except Exception:  # pragma: no cover - environment dependent
    pl = None
    _HAVE_PALLAS = False


def _bitmatmul_kernel(x_ref, y_ref, o_ref):
    # Each program instance owns one (row-tile, col-tile) of the output and
    # the full K dimension, accumulating across K blocks.
    acc = jnp.zeros_like(o_ref)
    x = x_ref[...]
    y = y_ref[...]
    acc = acc + jnp.dot(x, y, preferred_element_type=jnp.int32)
    o_ref[...] = jnp.mod(acc, 2).astype(jnp.int32)


def bitmatmul_mod2(X, Y, block_m: int = 64, block_n: int = 64, interpret: bool = True):
    """``(X @ Y) mod 2`` over 0/1 ``int32`` matrices, tiled with Pallas.

    Falls back to a plain ``jnp`` computation when Pallas is unavailable or
    the shapes do not tile evenly, so callers always get the right answer;
    the kernel path is what a TPU run would exercise.
    """
    X = jnp.asarray(X, dtype=jnp.int32)
    Y = jnp.asarray(Y, dtype=jnp.int32)
    M, K = X.shape
    K2, N = Y.shape
    assert K == K2, f"inner dims disagree: {K} vs {K2}"

    bm = min(block_m, M)
    bn = min(block_n, N)
    tiles_ok = _HAVE_PALLAS and M % bm == 0 and N % bn == 0
    if not tiles_ok:
        return jnp.mod(X @ Y, 2).astype(jnp.int32)

    grid = (M // bm, N // bn)
    return pl.pallas_call(
        _bitmatmul_kernel,
        grid=grid,
        in_specs=[
            pl.BlockSpec((bm, K), lambda i, j: (i, 0)),
            pl.BlockSpec((K, bn), lambda i, j: (0, j)),
        ],
        out_specs=pl.BlockSpec((bm, bn), lambda i, j: (i, j)),
        out_shape=jax.ShapeDtypeStruct((M, N), jnp.int32),
        interpret=interpret,
    )(X, Y)
