"""Batched JAX curve ops agree with the scalar oracle."""

import numpy as np
import pytest

from ic.curve import CurveJax, pack_batch
from ic.field import bits_to_int, int_to_bits
from ic.instances import enumerate_group, make_curve
from ic.reference import INF, pack


def _pt_arrays(pts, n):
    x = np.stack([int_to_bits(p[0], n) for p in pts]).astype(np.int32)
    y = np.stack([int_to_bits(p[1], n) for p in pts]).astype(np.int32)
    inf = np.array([1 if p[2] else 0 for p in pts], dtype=np.int32)
    return x, y, inf


def _to_pts(x, y, inf):
    x = np.asarray(x)
    y = np.asarray(y)
    inf = np.asarray(inf)
    out = []
    for i in range(x.shape[0]):
        if inf[i]:
            out.append(INF)
        else:
            out.append((bits_to_int(x[i]), bits_to_int(y[i]), False))
    return out


@pytest.mark.parametrize("n", [7, 9, 11])
def test_add_matches_oracle_all_pairs(n):
    curve = make_curve(n, a=0)
    cj = CurveJax(curve)
    pts = enumerate_group(curve)[: 60]  # cap the quadratic sweep
    P, Q = zip(*[(p, q) for p in pts for q in pts])
    Px, Py, Pinf = _pt_arrays(P, n)
    Qx, Qy, Qinf = _pt_arrays(Q, n)
    rx, ry, rinf = cj.add(Px, Py, Pinf, Qx, Qy, Qinf)
    got = _to_pts(rx, ry, rinf)
    for idx, (p, q) in enumerate(zip(P, Q)):
        assert got[idx] == curve.add(p, q), (n, p, q, got[idx], curve.add(p, q))


@pytest.mark.parametrize("n", [7, 9, 11])
def test_neg_and_frobenius(n):
    curve = make_curve(n, a=0)
    cj = CurveJax(curve)
    pts = enumerate_group(curve)
    x, y, inf = _pt_arrays(pts, n)
    nx, ny, ninf = cj.neg(x, y, inf)
    assert _to_pts(nx, ny, ninf) == [curve.neg(p) for p in pts]
    fx, fy, finf = cj.frobenius(x, y, inf)
    assert _to_pts(fx, fy, finf) == [curve.frobenius(p) for p in pts]


def test_pack_matches_oracle():
    curve = make_curve(9, a=0)
    pts = enumerate_group(curve)
    x, y, inf = _pt_arrays(pts, 9)
    packed = pack_batch(x, y, inf)
    packed = np.asarray(packed).reshape(-1)
    for i, p in enumerate(pts):
        assert int(packed[i]) == pack(p), (p, int(packed[i]), pack(p))
