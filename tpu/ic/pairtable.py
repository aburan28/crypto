"""Stage 1 -- pair-table relation collection.

The meet-in-the-middle half of a decomposition search, ported from
``gpu/ecc2k/pairtable.cuh``: build every sum ``P_i + P_j`` of the factor
base once, then "is ``R`` a sum of two base points" is a lookup and "of
three" is ``|F|`` lookups of ``R - P_k``.  A decomposition
``R = P_i + P_j + P_k`` with ``R = [a]G`` is exactly one relation of the
index calculus.

The division of labour is the honest part:

* **On the array (device, JAX):** all the elliptic-curve arithmetic --
  the ``|F|^2/2`` pair sums and the ``|targets| x |F|`` subtractions
  ``R - P_k`` -- as batched bit-matmuls (:mod:`tpu.ic.curve`).  This is
  the stage's dominant cost and the only part shaped like what a TPU is
  built for.
* **On the host:** packing the sums to 64-bit keys and the sort / binary
  search that answers membership.  A sort is not matmul; on a real TPU it
  would run on the vector unit via XLA ``sort`` over int32-pair keys, or
  come back to the host.  Pretending the MXU does it would be the kind of
  phase mislabelling ``AGENTS.md`` 5 warns against, so it is kept separate
  and named.

Keys are ``FastPoint::pack`` exactly, so a table built here and one built
by the Rust/CUDA backends sort and look up identically -- checked in
``tpu/tests/test_contract.py``.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional, Tuple

import jax.numpy as jnp
import numpy as np

from .curve import CurveJax, pack_batch
from .field import int_to_bits
from .reference import Curve


def _base_to_arrays(base, n: int):
    x = np.stack([int_to_bits(p[0], n) for p in base]).astype(np.int32)
    y = np.stack([int_to_bits(p[1], n) for p in base]).astype(np.int32)
    inf = np.array([1 if p[2] else 0 for p in base], dtype=np.int32)
    return jnp.asarray(x), jnp.asarray(y), jnp.asarray(inf)


@dataclass
class PairTable:
    """Every ``P_i + P_j`` (``i <= j``) keyed by packed identity, sorted."""

    keys: np.ndarray            # uint64 [P], sorted ascending
    pairs: np.ndarray           # int32 [P, 2] (i, j), aligned to keys
    n_base: int

    @staticmethod
    def build(cj: CurveJax, base) -> "PairTable":
        n = cj.n
        bx, by, binf = _base_to_arrays(base, n)
        F = len(base)
        # upper-triangle index lists
        ii, jj = np.triu_indices(F)  # i <= j
        ii = ii.astype(np.int64)
        jj = jj.astype(np.int64)
        # device: P_i + P_j for all upper-triangle pairs, in one batch
        Px, Py, Pinf = bx[ii], by[ii], binf[ii]
        Qx, Qy, Qinf = bx[jj], by[jj], binf[jj]
        sx, sy, sinf = cj.add(Px, Py, Pinf, Qx, Qy, Qinf)
        keys = pack_batch(sx, sy, sinf).astype(np.uint64)
        order = np.argsort(keys, kind="stable")
        keys = keys[order]
        pairs = np.stack([ii[order], jj[order]], axis=1).astype(np.int32)
        return PairTable(keys=keys, pairs=pairs, n_base=F)

    # -- membership (host sort/search) -----------------------------------
    def first_pair(self, key: int) -> Optional[Tuple[int, int]]:
        pos = np.searchsorted(self.keys, np.uint64(key), side="left")
        if pos < len(self.keys) and self.keys[pos] == np.uint64(key):
            i, j = self.pairs[pos]
            return int(i), int(j)
        return None

    def contains(self, key: int) -> bool:
        return self.first_pair(key) is not None


def _neg_base_arrays(cj: CurveJax, base):
    n = cj.n
    bx, by, binf = _base_to_arrays(base, n)
    nx, ny, ninf = cj.neg(bx, by, binf)
    return nx, ny, ninf


def collect_m3(
    cj: CurveJax,
    table: PairTable,
    base,
    targets: List[Tuple[int, int, bool]],
    scalars: Optional[List[int]] = None,
) -> List[dict]:
    """For each target ``R`` (optionally tagged with its scalar ``a`` so
    ``R = [a]G``), find one ``R = P_i + P_j + P_k`` via the table.

    Device work: ``R - P_k`` for every ``(target, k)`` as one batched add
    of ``R`` with the negated base.  Host work: pack and look the
    differences up.  Returns relation dicts ``{"a", "points": [i, j, k]}``
    with ``points`` sorted; the sum is **not** trusted here -- the caller
    re-checks every relation in the group, exactly as
    ``verify_collected_relation`` does on the Rust side.
    """
    n = cj.n
    F = len(base)
    T = len(targets)
    tx = jnp.asarray(np.stack([int_to_bits(t[0], n) for t in targets]).astype(np.int32))
    ty = jnp.asarray(np.stack([int_to_bits(t[1], n) for t in targets]).astype(np.int32))
    tinf = jnp.asarray(np.array([1 if t[2] else 0 for t in targets], dtype=np.int32))

    nx, ny, ninf = _neg_base_arrays(cj, base)  # [F, n]

    # broadcast to [T, F, n]: R_t - P_k
    Rx = tx[:, None, :]
    Ry = ty[:, None, :]
    Rinf = tinf[:, None]
    Nx = nx[None, :, :]
    Ny = ny[None, :, :]
    Ninf = ninf[None, :]
    dx, dy, dinf = cj.add(
        jnp.broadcast_to(Rx, (T, F, n)),
        jnp.broadcast_to(Ry, (T, F, n)),
        jnp.broadcast_to(Rinf, (T, F)),
        jnp.broadcast_to(Nx, (T, F, n)),
        jnp.broadcast_to(Ny, (T, F, n)),
        jnp.broadcast_to(Ninf, (T, F)),
    )
    diff_keys = pack_batch(dx, dy, dinf).astype(np.uint64)  # [T, F]

    out = []
    for t in range(T):
        found = None
        for k in range(F):
            pr = table.first_pair(int(diff_keys[t, k]))
            if pr is not None:
                i, j = pr
                found = sorted([i, j, k])
                break
        if found is not None:
            rec = {"points": found}
            if scalars is not None:
                rec["a"] = int(scalars[t])
            out.append(rec)
    return out


def verify_relation_m3(curve: Curve, base, points: List[int], target) -> bool:
    """Re-check ``target == P_i + P_j + P_k`` in the group (the oracle)."""
    if len(points) != 3 or any(p < 0 or p >= len(base) for p in points):
        return False
    acc = (0, 0, True)
    for p in points:
        acc = curve.add(acc, base[p])
    return acc == target
