"""Canonical genus-5 curves: smooth complete intersections of three quadrics
in P^4 over F_2, sampled as uniform random nets (3-dim subspaces of the
15-dim space of quadratic forms).

Stage 1 (all sampled nets): N_1, N_2, N_3 of the scheme V(Q1,Q2,Q3).
Stage 2 (nets passing a target filter on (N_1,N_2,N_3), plus a uniform
subsample): N_4, N_5 and a Jacobian-criterion singularity search over
F_{2^k}, k = 1..5 (k = 6 on request).
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from gf2k import Field, monomials

EXPS = monomials(5, 2)  # 15 quadratic monomials
FIELDS = {}
PTS = {}
VALS = {}


def field(k):
    if k not in FIELDS:
        F = Field(k)
        FIELDS[k] = F
        pts = F.proj_points(4)
        PTS[k] = pts
        VALS[k] = F.monomial_values(pts, EXPS).astype(np.uint8)
    return FIELDS[k]


def rank_lt3(F, J):
    """J: (npts, 3, 5) field elements. True where rank < 3 (all 3x3 minors 0)."""
    mul = F.mul
    from itertools import combinations
    out = np.ones(J.shape[0], dtype=bool)
    for (a, b, c) in combinations(range(5), 3):
        m = (mul[J[:, 0, a], mul[J[:, 1, b], J[:, 2, c]]] ^ mul[J[:, 0, a], mul[J[:, 1, c], J[:, 2, b]]]
             ^ mul[J[:, 0, b], mul[J[:, 1, a], J[:, 2, c]]] ^ mul[J[:, 0, b], mul[J[:, 1, c], J[:, 2, a]]]
             ^ mul[J[:, 0, c], mul[J[:, 1, a], J[:, 2, b]]] ^ mul[J[:, 0, c], mul[J[:, 1, b], J[:, 2, a]]])
        out &= m == 0
    return out


def count_and_sing(net, k, chunk=None):
    """net: 3 ints (15-bit coefficient masks). Returns (#points, has_singular_point)."""
    F = field(k) if k <= 5 else None
    if k <= 5:
        V, pts = VALS[k], PTS[k]
        zero = np.ones(len(pts), dtype=bool)
        for qm in net:
            supp = [m for m in range(15) if (qm >> m) & 1]
            zero &= np.bitwise_xor.reduce(V[:, supp], axis=1) == 0
        cp = pts[zero]
    else:
        F = Field(k)
        cp = curve_points_large(net, F)
    n = len(cp)
    if n == 0:
        return 0, False
    # Jacobian: dQ/dx_l = sum_{j != l} c_{lj} x_j
    J = np.zeros((n, 3, 5), dtype=np.int64)
    for qi, qm in enumerate(net):
        for m, e in enumerate(EXPS):
            if not (qm >> m) & 1:
                continue
            idx = [i for i in range(5) for _ in range(e[i])]
            if idx[0] != idx[1]:
                l, j = idx
                J[:, qi, l] ^= cp[:, j]
                J[:, qi, j] ^= cp[:, l]
    return n, bool(rank_lt3(F, J).any())


def curve_points_large(net, F):
    """Points of V(net) over F_{2^k} for k = 6, enumerated in slices of P^4."""
    q = F.q
    out = []
    import itertools
    for lead in range(5):
        rest = 4 - lead
        if rest == 0:
            arr = np.array([[0] * lead + [1]], dtype=np.int64)
            blocks = [arr]
        else:
            # iterate over first free coordinate, vectorise the others
            tails = np.array(list(itertools.product(range(q), repeat=rest - 1)), dtype=np.int64) \
                if rest > 1 else np.zeros((1, 0), dtype=np.int64)
            blocks = []
            for x in range(q):
                arr = np.concatenate([np.zeros((len(tails), lead), dtype=np.int64),
                                      np.ones((len(tails), 1), dtype=np.int64),
                                      np.full((len(tails), 1), x, dtype=np.int64), tails], axis=1)
                blocks.append(arr)
        for arr in blocks:
            V = F.monomial_values(arr, EXPS)
            zero = np.ones(len(arr), dtype=bool)
            for qm in net:
                supp = [m for m in range(15) if (qm >> m) & 1]
                zero &= np.bitwise_xor.reduce(V[:, supp], axis=1) == 0
            if zero.any():
                out.append(arr[zero])
    return np.concatenate(out, axis=0) if out else np.zeros((0, 5), dtype=np.int64)


def full(net, kmax=5):
    N, sing = [], False
    for k in range(1, kmax + 1):
        n, s = count_and_sing(net, k)
        N.append(n)
        sing |= s
    return N, sing


def stage1_tables(ks=(1, 2, 3)):
    """Span tables for fast cheap counts."""
    T = {}
    for k in ks:
        field(k)
        V = VALS[k]
        lo = np.zeros((128, V.shape[0]), dtype=np.uint8)
        for j in range(7):
            lo[1 << j:1 << (j + 1)] = lo[:1 << j] ^ V[:, j][None, :]
        hi = np.zeros((256, V.shape[0]), dtype=np.uint8)
        for j in range(8):
            hi[1 << j:1 << (j + 1)] = hi[:1 << j] ^ V[:, 7 + j][None, :]
        T[k] = (lo, hi)
    return T


def stage1(nets, T):
    """nets: (B,3) int array. Returns (B, len(T)) counts."""
    out = []
    for k, (lo, hi) in sorted(T.items()):
        z = None
        for i in range(3):
            c = nets[:, i]
            v = lo[c & 127] ^ hi[c >> 7]
            z = v if z is None else (z | v)
        out.append((z == 0).sum(1))
    return np.stack(out, axis=1)


def random_nets(rng, B):
    nets = rng.integers(1, 1 << 15, size=(B, 3))
    # reject dependent triples: a, b, c independent iff a!=b, c not in {a, b, a^b}
    a, b, c = nets[:, 0], nets[:, 1], nets[:, 2]
    ok = (a != b) & (c != a) & (c != b) & (c != (a ^ b))
    return nets[ok]


def gauss_binom(n, k, q=2):
    num = den = 1
    for i in range(k):
        num *= q ** (n - i) - 1
        den *= q ** (i + 1) - 1
    return num // den


def gl_order(n, q=2):
    r = 1
    for i in range(n):
        r *= q ** n - q ** i
    return r
