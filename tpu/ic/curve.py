"""Batched binary-Koblitz point arithmetic for the TPU.

A batch of points is three arrays: ``x`` and ``y`` of shape ``[..., n]``
(bit-vectors, :mod:`tpu.ic.field`) and ``inf`` of shape ``[...]`` (``int32``
0/1).  Every operation computes all branches of the affine law and selects
with :func:`jax.numpy.where`, so there is no data-dependent control flow --
the uniform, lane-parallel shape XLA wants.

The addition formulas are transcribed from ``FastCurve`` in
``src/cryptanalysis/koblitz_fast.rs`` and checked against
:mod:`tpu.ic.reference`.

Why no Montgomery batch-inversion.  The CPU/CUDA pair-table amortizes one
field inversion over a whole row with Montgomery's trick, because an
inversion there is ~``n`` multiplications and the trick trades ``k``
inversions for one inversion and ``3(k-1)`` multiplications.  In this
backend an inversion is a fixed square-and-multiply chain of bit-matmuls
(:meth:`FieldJax.inv`) applied to the *whole batch at once*; that batched
chain already is the vectorized form, and threading Montgomery's sequential
prefix product through it would serialize what the array parallelizes.  So
the pair sum below inverts the batch directly.
"""

from __future__ import annotations

import jax.numpy as jnp
import numpy as np

from .field import FieldJax, bits_to_int, int_to_bits
from .reference import Curve


def _all_eq(u, v):
    return jnp.all(u == v, axis=-1)


def _all_zero(u):
    return jnp.all(u == 0, axis=-1)


class CurveJax:
    def __init__(self, curve: Curve):
        self.curve = curve
        self.f = FieldJax(curve.field)
        self.n = curve.field.n
        self.a = jnp.asarray(int_to_bits(curve.a, self.n))
        self.b = jnp.asarray(int_to_bits(curve.b, self.n))
        self.one = self.f.one

    # -- unary ------------------------------------------------------------
    def neg(self, x, y, inf):
        return x, jnp.bitwise_xor(x, y), inf

    def frobenius(self, x, y, inf):
        return self.f.sqr(x), self.f.sqr(y), inf

    # -- addition (full branch select) -----------------------------------
    def add(self, Px, Py, Pinf, Qx, Qy, Qinf):
        f = self.f
        a = jnp.broadcast_to(self.a, Px.shape)

        # generic branch: distinct abscissae
        dx = f.add(Px, Qx)
        dy = f.add(Py, Qy)
        lam = f.mul(dy, f.inv(dx))
        gx = f.add(f.add(f.add(f.sqr(lam), lam), f.add(Px, Qx)), a)
        gy = f.add(f.add(f.mul(lam, f.add(Px, gx)), gx), Py)

        # doubling branch: equal abscissae, P != -Q, x != 0
        lam_d = f.add(Px, f.mul(Py, f.inv(Px)))
        dx3 = f.add(f.add(f.sqr(lam_d), lam_d), a)
        dy3 = f.add(f.sqr(Px), f.mul(f.add(lam_d, jnp.broadcast_to(self.one, Px.shape)), dx3))

        same_x = _all_eq(Px, Qx)
        is_neg = same_x & _all_eq(dy, Px)        # p.y ^ q.y == p.x
        px_zero = _all_zero(Px)
        do_double = same_x & (~is_neg)
        # 2-torsion (same x, x==0) doubles to O
        double_to_inf = do_double & px_zero

        ddbl = do_double[..., None]

        # choose coordinates: doubling where same_x, else generic
        rx = jnp.where(ddbl, dx3, gx)
        ry = jnp.where(ddbl, dy3, gy)

        # infinity flag of the raw sum (before P/Q infinity handling)
        raw_inf = (is_neg | double_to_inf).astype(jnp.int32)

        # P infinite -> Q ; Q infinite -> P  (checked after the arithmetic)
        pinf = (Pinf != 0)
        qinf = (Qinf != 0)
        bx = jnp.where(qinf[..., None], Px, rx)
        by = jnp.where(qinf[..., None], Py, ry)
        binf = jnp.where(qinf, Pinf, raw_inf)
        ox = jnp.where(pinf[..., None], Qx, bx)
        oy = jnp.where(pinf[..., None], Qy, by)
        oinf = jnp.where(pinf, Qinf, binf)
        return ox, oy, oinf.astype(jnp.int32)


# --------------------------------------------------------------------------
# host-side packing (the key a sum is stored under)
# --------------------------------------------------------------------------


def pack_batch(x, y, inf) -> np.ndarray:
    """``FastPoint::pack`` over a batch, on the host.

    Packing compares ``y`` and ``x ^ y`` as integers and shifts ``x + 1``
    left, which needs full ``n``-bit (<=62) integers; it is a terminal,
    non-arithmetic step, so it is done here in Python rather than on the
    device, where 64-bit integers are second-class.  This is the honest
    seam between what the TPU computes (the sums, as bit-vectors) and what
    the host does with the keys.
    """
    x = np.asarray(x)
    y = np.asarray(y)
    inf = np.asarray(inf).reshape(-1)
    flat_x = x.reshape(-1, x.shape[-1])
    flat_y = y.reshape(-1, y.shape[-1])
    out = np.zeros(flat_x.shape[0], dtype=np.uint64)
    for i in range(flat_x.shape[0]):
        if inf[i]:
            out[i] = 0
            continue
        xi = bits_to_int(flat_x[i])
        yi = bits_to_int(flat_y[i])
        sign = 1 if yi > (xi ^ yi) else 0
        out[i] = np.uint64(((xi + 1) << 1) | sign)
    return out.reshape(x.shape[:-1])
