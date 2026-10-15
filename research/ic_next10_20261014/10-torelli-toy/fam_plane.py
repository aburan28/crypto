"""Plane curves over F_2: all smooth plane curves of degree d, or (trigonal
genus 5) all plane quintics with a delta = 1 double point at P0 = [0:0:1] and
no other singular point.

Singular points are searched over F_{2^k}, k = 1..6 (an irreducible plane
curve of degree <= 5 has at most 6 singular points, so every Galois orbit of
them is seen).  Point counts N_k, k = 1..g, on the smooth model.

Method: forms are F_2-linear in their coefficient bits, so the value and the
three partials at every point are packed in one 32-bit word per (point,
monomial) and combined by XOR; the coefficient vector is split low/high and a
model's word at P is VL[low] ^ VH[high].
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from gf2k import Field, monomials


def build(d, ks):
    exps = monomials(3, d)
    segs, cols = [], []
    off = 0
    p0_idx = []
    for k in ks:
        F = Field(k)
        pts = F.proj_points(2)
        V = F.monomial_values(pts, exps)
        W = V.astype(np.uint32)
        for var in range(3):
            # partial derivative of monomial wrt var: exponent odd -> x^{e-1}...
            dexps, keep = [], []
            for e in exps:
                if e[var] % 2 == 1:
                    e2 = list(e); e2[var] -= 1
                    dexps.append(tuple(e2)); keep.append(True)
                else:
                    dexps.append(tuple([0, 0, 0])); keep.append(False)
            D = F.monomial_values(pts, dexps)
            D[:, ~np.array(keep)] = 0
            W |= D.astype(np.uint32) << np.uint32(8 * (var + 1))
        cols.append(W)
        segs.append((off, off + len(pts)))
        p0 = np.nonzero((pts[:, 0] == 0) & (pts[:, 1] == 0))[0]
        p0_idx.append(off + int(p0[0]))
        off += len(pts)
    return exps, np.concatenate(cols, axis=0), segs, p0_idx


def span_table(Wsub):
    """T[i] = XOR of columns j with bit j of i set; Wsub is npts x m."""
    m = Wsub.shape[1]
    T = np.zeros((1 << m, Wsub.shape[0]), dtype=np.uint32)
    for j in range(m):
        T[1 << j:1 << (j + 1)] = T[:1 << j] ^ Wsub[:, j][None, :]
    return T


def enumerate_plane(d, g, trigonal=False, ks_sing=(1, 2, 3, 4, 5, 6)):
    ks = sorted(set(range(1, g + 1)) | set(ks_sing))
    exps, W, segs, p0_idx = build(d, ks)
    kpos = {k: i for i, k in enumerate(ks)}
    free = list(range(len(exps)))
    if trigonal:
        forced = [exps.index((0, 0, 5)), exps.index((1, 0, 4)), exps.index((0, 1, 4))]
        free = [i for i in free if i not in forced]
    M = len(free)
    L = M // 2
    lowc, highc = free[:L], free[L:]
    VL = span_table(W[:, lowc])
    VH = span_table(W[:, highc])
    singmask = np.ones(W.shape[0], dtype=bool)
    if trigonal:
        singmask[p0_idx] = False
        iq = [exps.index(e) for e in ((2, 0, 3), (1, 1, 3), (0, 2, 3))]
        i3 = {e: exps.index(e) for e in exps if e[2] == 2}
    hb = np.arange(1 << (M - L))
    rows, nsmooth, ntot = [], 0, 0
    for a in range(1 << L):
        X = VL[a][None, :] ^ VH
        on = (X & np.uint32(0xFF)) == 0
        sing = (X == 0) & singmask[None, :]
        ok = ~sing.any(1)
        bits = np.zeros((len(hb), len(exps)), dtype=np.int64)
        for j, c in enumerate(lowc):
            bits[:, c] = (a >> j) & 1
        for j, c in enumerate(highc):
            bits[:, c] = (hb >> j) & 1
        corr = np.zeros((len(hb), g), dtype=np.int64)
        if trigonal:
            qa, qb, qc = bits[:, iq[0]], bits[:, iq[1]], bits[:, iq[2]]
            ok &= (qa | qb | qc) == 1
            # cusp: delta = 1 iff f3(c, a) != 0
            f3 = np.zeros(len(hb), dtype=np.int64)
            for e, idx in i3.items():
                term = bits[:, idx].copy()
                if e[0]:
                    term &= qc
                if e[1]:
                    term &= qa
                f3 ^= term
            ok &= (qb == 1) | (f3 == 1)
            for k in range(1, g + 1):
                split = (qa == 0) | (qc == 0)
                r_node = np.where(split, 2, 2 if k % 2 == 0 else 0)
                corr[:, k - 1] = np.where(qb == 1, r_node, 1) - 1
        ntot += len(hb)
        if not ok.any():
            continue
        N = np.stack([on[:, segs[kpos[k]][0]:segs[kpos[k]][1]].sum(1) for k in range(1, g + 1)],
                     axis=1) + corr
        rows.append(N[ok].astype(np.int16))
        nsmooth += int(ok.sum())
    allrows = np.concatenate(rows, axis=0)
    u, c = np.unique(allrows, axis=0, return_counts=True)
    out = {tuple(int(x) for x in r): int(n) for r, n in zip(u, c)}
    return out, nsmooth, ntot


if __name__ == "__main__":
    d, g = int(sys.argv[1]), int(sys.argv[2])
    tri = len(sys.argv) > 3
    out, ns, nt = enumerate_plane(d, g, trigonal=tri)
    print(d, g, tri, nt, ns, len(out))
