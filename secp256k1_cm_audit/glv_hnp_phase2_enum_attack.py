"""
GLV-HNP Phase 2, Thread 23c: past the K1 wall by enumeration + box test.

Thread 23b (glv_hnp_phase2_bdd_optimality.py) established, on the two
historical 12-bit curves at m=12:

  E5  d is UNIQUE for every (K1, seed) tested, up to K1=12 (eff=0.236).  The
      wall is therefore not a shortage of data.
  E7  the LLL wall coincides exactly with the rank-1 boundary of the true d in
      the coset-BDD ordering:
        12-bit/2557  K1=4,6,8 -> rank 1 (LLL 5/5);  K1=12 -> ranks 1,1,2,1,3
                                                    (LLL 4/5)
        12-bit/2677  K1=4     -> rank 1 (5/5);  K1=6 -> 1,1,1,1,3 (2/5);
                     K1=8     -> 2,21,19,2,9 (0/5);  K1=12 -> 88..247 (0/5)

So the obstruction is an OBJECTIVE MISMATCH, not hardness: the lattice minimises
an l2 norm, but the true solution is defined by a BOX constraint
(0 <= k1_i < K1, 0 <= k2_i < K2).  Past the wall a wrong d' has a cheaper l2
representative while violating the box.  The l2 optimum is the wrong answer, so
no strengthening of the reduction (BKZ-40 was already tried, 2026-07-26) can
help -- but the true answer is still only a handful of ranks down the list.

This script tests the consequence: enumerate the R shortest vectors of the
Kannan coset instead of taking only the shortest, and filter with the box test
(the `representable` predicate, O(m*K2) per candidate).  The radius is derived
from PUBLIC parameters only (n, K1, K2, m) -- no knowledge of d.

  E8  LLL (1 candidate) vs enumeration with R in {1, 8, 64, 256}, m=12,
      K1 in {4, 6, 8, 12}, 5 seeds, on both historical curves.

Run: python3 glv_hnp_phase2_enum_attack.py
"""

import importlib.util
import math
import random
import time

from fpylll import GSO, IntegerMatrix, LLL, Enumeration, EnumerationError

_spec = importlib.util.spec_from_file_location(
    "_t23", __file__.rsplit("/", 1)[0] + "/glv_hnp_phase2_projected.py")
_t23 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_t23)

find_generator = _t23.find_generator
gen_signatures = _t23.gen_signatures
scales = _t23.scales
build_glv_lattice = _t23.build_glv_lattice
planted_vector = _t23.planted_vector
norm = _t23.norm
modinv = _t23.modinv


def expected_radius_sq(n, m, k1_bound, k2_bound, slack=1.35):
    """Radius from public parameters only.  E[x^2] = X^2/3 for x ~ U[0,X)."""
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    exp_sq = (m * (k1_bound * S_K1) ** 2 / 3.0
              + (n * S_D) ** 2 / 3.0
              + m * (k2_bound * S_K2) ** 2 / 3.0
              + S_KAN ** 2)
    return slack * exp_sq


def box_test(d, sigs, n, lam, k1_bound, k2_bound):
    """Does candidate d admit 0<=k1_i<K1, 0<=k2_i<K2 for EVERY signature?"""
    for s in sigs:
        k = (s['A'] + s['B'] * d) % n
        for k2 in range(k2_bound):
            if (k - lam * k2) % n < k1_bound:
                break
        else:
            return False
    return True


def enum_attack(sigs, n, lam, k1_bound, k2_bound, nr_solutions):
    """Enumerate the nr_solutions shortest coset vectors, box-test each.

    Returns (d or None, n_candidates_examined).  Uses only public data.
    """
    m = len(sigs)
    dim = 2 * m + 2
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    M = build_glv_lattice(sigs, n, lam, k1_bound, k2_bound)
    if nr_solutions == 1:
        # plain LLL, exactly the historical attack: read d off any row whose
        # last coordinate is +-S_KANNAN
        A = IntegerMatrix.from_matrix(M)
        LLL.reduction(A)
        cands = []
        for i in range(A.nrows):
            row = [A[i][j] for j in range(dim)]
            if abs(row[dim - 1]) == S_KAN:
                sign = 1 if row[dim - 1] > 0 else -1
                cands.append(sign * row[m] % n)
    else:
        # CVP-mode enumeration INSIDE the coset {last coordinate = S_KANNAN}.
        # (Plain SVP enumeration on the full lattice is useless here: every
        # short vector of L lies in the last-coordinate-zero sublattice, headed
        # by the trivial n*S_D*e_m of norm n -- Thread 23b, result T5.)
        body = [row[:dim - 1] for row in M[:-1]]
        t = M[-1][:dim - 1]
        B = IntegerMatrix.from_matrix(body)
        LLL.reduction(B)
        g = GSO.Mat(B, update=True)
        radius = expected_radius_sq(n, m, k1_bound, k2_bound) - S_KAN ** 2
        try:
            sols = Enumeration(g, nr_solutions=nr_solutions).enumerate(
                0, B.nrows, radius, 0,
                target=g.from_canonical(tuple(-x for x in t)))
        except EnumerationError:
            sols = []
        cands = []
        for _dist, coeffs in sols:
            v = list(t)
            for i, c in enumerate(coeffs):
                ci = int(c)
                if ci:
                    for j in range(dim - 1):
                        v[j] += ci * B[i][j]
            cands.append(v[m] % n)
    seen = set()
    for d in cands:
        if d in seen or d == 0:
            continue
        seen.add(d)
        if box_test(d, sigs, n, lam, k1_bound, k2_bound):
            return d, len(seen)
    return None, len(seen)


def main():
    print("=" * 78)
    print("GLV-HNP Phase 2, Thread 23c — enumeration + box test past the K1 wall")
    print("=" * 78)

    SEEDS = [42, 1234, 9999, 555, 31337]
    HIST = [("12-bit/2557", 2557, 2, 2659, 1755),
            ("12-bit/2677", 2677, 2, 2647, 185)]
    curves = []
    for label, p, b, n, lam in HIST:
        G = find_generator(p, b, n)
        assert G is not None and (lam * lam + lam + 1) % n == 0
        curves.append((label, (p, b, n, lam, G)))

    M_SIGS = 12
    K1_GRID = [4, 6, 8, 12]
    R_GRID = [1, 8, 64, 256]

    print("\nm=12, 5 seeds.  R = number of enumerated coset vectors")
    print("(R=1 is plain LLL = the historical attack).  Radius from public")
    print("parameters only: 1.35 * E[||v_planted||^2].\n")
    print(f"{'curve':<13} {'K1':>3} {'eff':>6} " +
          " ".join(f"{'R='+str(r):>7}" for r in R_GRID) + f" {'secs':>7}")

    for label, curve in curves:
        p, b, n, lam, G = curve
        k2b = math.isqrt(n) + 1
        for k1 in K1_GRID:
            cells, t0, distinct = [], time.time(), []
            for R in R_GRID:
                wins = 0
                for seed in SEEDS:
                    d_true = random.Random(seed + 7777).randint(1, n - 1)
                    sigs = gen_signatures(G, d_true, M_SIGS, n, lam, p, k1, k2b, seed)
                    d_found, nseen = enum_attack(sigs, n, lam, k1, k2b, R)
                    if R == R_GRID[-1]:
                        distinct.append(nseen)
                    if d_found is not None:
                        # the box test is the only acceptance criterion; scoring
                        # against d_true just records whether it was right
                        wins += (d_found == d_true)
                cells.append(f"{wins}/{len(SEEDS)}".rjust(7))
            print(f"{label:<13} {k1:>3} {k1*k2b/n:>6.3f} " + " ".join(cells) +
                  f" {time.time()-t0:>7.1f}  distinct d at R={R_GRID[-1]}: "
                  f"{distinct}")

    print("\n" + "=" * 78)
    print("done")
    print("=" * 78)


if __name__ == "__main__":
    main()
