"""Stage 2 -- GF(2) linear algebra on the array.

The relation matrix of an index calculus is solved by linear algebra, and
over GF(2) that algebra is the one stage of the whole pipeline that is
*natively* the TPU's shape: a GF(2) matrix product is an integer matmul
taken ``mod 2``, which is precisely what the systolic array does.

This module gives:

* :func:`gf2_matmul` -- the MXU op, ``(A @ B) mod 2``.
* :func:`gf2_reduce_matrix` -- Gauss-Jordan elimination over GF(2),
  written swap-free (each column picks a pivot among the not-yet-pivot
  rows and XORs it into every other row carrying that column) so it is a
  static, branch-free XLA graph rather than host control flow.
* :func:`gf2_rank`, :func:`gf2_solve`, :func:`gf2_nullspace_vector`.

Honesty markers (``AGENTS.md`` 5, 8):

* The real ECDLP log solve is over ``F_r`` for the large prime subgroup
  order ``r``, **not** GF(2); modular arithmetic mod a 95--131-bit prime
  is a poor MXU fit and is out of scope here.  This module is the GF(2)
  flavour (the kind the Semaev / WDSat relation systems and parity
  arguments use), offered as the stage whose *shape* a TPU fits.
* The elimination below is the textbook ``O(n^3)`` form, not the blocked
  "four Russians" / recursive-matmul form a TPU would want for size.  It
  is a correctness reference, not a performance claim.  No timing is made.

Checked against a plain NumPy oracle in ``tpu/tests/test_linalg.py``.
"""

from __future__ import annotations

from typing import Optional, Tuple

import jax
import jax.numpy as jnp
import numpy as np


def gf2_matmul(A, B):
    """``(A @ B) mod 2`` -- the MXU-native GF(2) product."""
    return jnp.mod(A.astype(jnp.int32) @ B.astype(jnp.int32), 2).astype(jnp.int32)


def gf2_reduce_matrix(A):
    """Gauss-Jordan over GF(2).

    Returns ``(R, pivot_row_is_pivot, pivot_col_of_row, rank)`` where ``R``
    is in reduced row-echelon form (up to row order): every pivot column
    holds a single 1, in its pivot row.
    """
    A = jnp.asarray(A, dtype=jnp.int32)
    R, C = A.shape
    is_pivot = jnp.zeros((R,), dtype=jnp.int32)
    pivot_col_of_row = jnp.full((R,), -1, dtype=jnp.int32)
    rows = jnp.arange(R)

    for col in range(C):
        colv = A[:, col]
        avail = (colv == 1) & (is_pivot == 0)
        has = jnp.any(avail)
        prow = jnp.argmax(avail.astype(jnp.int32))  # first available row
        prow_onehot = (rows == prow) & has
        pivot_vec = A[prow]
        # rows carrying this column, excluding the pivot row itself
        clear = (colv == 1) & (~prow_onehot)
        add = jnp.where((clear & has)[:, None], pivot_vec[None, :], 0)
        A = jnp.bitwise_xor(A, add.astype(jnp.int32))
        is_pivot = jnp.where(prow_onehot, 1, is_pivot)
        pivot_col_of_row = jnp.where(prow_onehot, col, pivot_col_of_row)

    rank = jnp.sum(is_pivot)
    return A, is_pivot, pivot_col_of_row, rank


def gf2_rank(A) -> int:
    _, _, _, rank = gf2_reduce_matrix(A)
    return int(rank)


def gf2_solve(A, b) -> Optional[np.ndarray]:
    """One solution ``x`` of ``A x = b`` over GF(2), or ``None`` if the
    system is inconsistent.  Free variables are set to 0."""
    A = np.asarray(A, dtype=np.int32)
    b = np.asarray(b, dtype=np.int32).reshape(-1)
    R, C = A.shape
    aug = np.concatenate([A, b[:, None]], axis=1)
    red, is_pivot, pcol, rank = gf2_reduce_matrix(aug)
    red = np.asarray(red)
    is_pivot = np.asarray(is_pivot)
    pcol = np.asarray(pcol)
    x = np.zeros((C,), dtype=np.int32)
    for r in range(R):
        if is_pivot[r]:
            c = int(pcol[r])
            if c == C:  # pivot in the augmented column -> 0 = 1
                return None
            x[c] = red[r, C]
    # consistency: A x must equal b
    if np.any((A @ x) % 2 != b):
        return None
    return x


def gf2_nullspace_vector(A) -> Optional[np.ndarray]:
    """A nonzero ``x`` with ``A x = 0`` over GF(2), or ``None`` if the
    only solution is trivial (columns independent).  This is the shape of
    "a dependency among relations" that an index calculus looks for."""
    A = np.asarray(A, dtype=np.int32)
    R, C = A.shape
    red, is_pivot, pcol, rank = gf2_reduce_matrix(A)
    red = np.asarray(red)
    is_pivot = np.asarray(is_pivot)
    pcol = np.asarray(pcol)
    pivot_cols = set(int(pcol[r]) for r in range(R) if is_pivot[r])
    free = [c for c in range(C) if c not in pivot_cols]
    if not free:
        return None
    f = free[0]
    x = np.zeros((C,), dtype=np.int32)
    x[f] = 1
    for r in range(R):
        if is_pivot[r]:
            c = int(pcol[r])
            x[c] = red[r, f]  # pivot var = coefficient on the chosen free var
    assert np.all((A @ x) % 2 == 0)
    return x
