"""
GLV-HNP Phase 2, Thread 23b: is the K1 wall information-theoretic or algorithmic?

Thread 23 (glv_hnp_phase2_projected.py) falsified the "make the planted vector
lambda_1" reformulation: deleting the d column raises sv/pv (0.42 -> 0.81 on the
2026-07-26 failure curve) but never above 1, and it makes recovery strictly
WORSE (12-bit/2557 at K1=4: 5/5 -> 2/5; 17-bit eff=0.15 sweep: 2/12 -> 0/12).

That leaves the second branch of the 2026-07-29 falsifier: "the wall is
information-theoretic".  This script settles it with GROUND TRUTH rather than
with lattice heuristics, using the fact that the toy moduli are small enough to
enumerate.

  E5  exhaustive d-uniqueness.  For every candidate d' in [0,n), every
      signature must satisfy  A_i + B_i*d' = k1_i + lam*k2_i (mod n)  with
      0 <= k1_i < K1, 0 <= k2_i < K2.  Representability is decided exactly by
      looping k2 over [0,K2) and testing (k - lam*k2) mod n < K1.  Cost
      O(n*m*K2) ~ 10^6 per instance.
        - count == 1  =>  the instance determines d uniquely; any failure of
                          LLL/BKZ at that point is ALGORITHMIC.
        - count > 1   =>  the instance itself is ambiguous; no algorithm on
                          this lattice can win => INFORMATION-THEORETIC.

  E6  exact BDD optimality.  Even when d is unique, LLL only wins if the
      planted vector is the shortest vector of the coset {last coord = S_KANNAN}.
      Solve that exactly by CVP (fpylll) in the rank-(2m+1) sublattice L_0 of
      vectors with last coordinate 0, target = the Kannan row.
        - CVP optimum == planted  =>  the lattice formulation is sound and the
                                      failure is a reduction-strength failure.
        - CVP optimum != planted  =>  the EMBEDDING is lossy: a shorter wrong
                                      answer exists, so no amount of BKZ helps.

  E7  BDD rank of the true d.  Every coset vector is parameterised by
      (d', k2_1..k2_m): the k1 columns are then forced, and the cheapest
      representative of each is the centred residue.  So
        cost(d') = min over k2 of  sum_i [ (cent(c_i - lam*k2_i)*S_K1)^2
                                           + (k2_i*S_K2)^2 ]
                   + (d'*S_D)^2 + S_KANNAN^2,     c_i = A_i + B_i*d' mod n,
      and the per-i minimisation decouples.  Sorting all n values of cost(d')
      gives the exact rank of the true d in the coset-BDD ordering.  Rank 1 =
      LLL/BKZ can in principle win; rank r = an attack must enumerate r coset
      vectors, which is the concrete price of going past the K1 wall.

Run: python3 glv_hnp_phase2_bdd_optimality.py
"""

import importlib.util
import math
import random
import time

import numpy as np
from fpylll import IntegerMatrix, LLL, CVP

_spec = importlib.util.spec_from_file_location(
    "_t23", __file__.rsplit("/", 1)[0] + "/glv_hnp_phase2_projected.py")
_t23 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_t23)

find_generator = _t23.find_generator
gen_signatures = _t23.gen_signatures
scales = _t23.scales
build_glv_lattice = _t23.build_glv_lattice
build_glv_lattice_projected = _t23.build_glv_lattice_projected
planted_vector = _t23.planted_vector
planted_vector_projected = _t23.planted_vector_projected
recover_d = _t23.recover_d
recover_d_projected = _t23.recover_d_projected
lll_rows = _t23.lll_rows
norm = _t23.norm
lam_star = _t23.lam_star
modinv = _t23.modinv

# ---------------------------------------------------------------------------
# E5 helper: exact representability and exhaustive d enumeration
# ---------------------------------------------------------------------------

def representable(k, n, lam, k1_bound, k2_bound):
    """Is k = k1 + lam*k2 (mod n) with 0<=k1<K1, 0<=k2<K2?  Returns (k1,k2) or None."""
    for k2 in range(k2_bound):
        k1 = (k - lam * k2) % n
        if k1 < k1_bound:
            return k1, k2
    return None

def enumerate_d(sigs, n, lam, k1_bound, k2_bound):
    """All d' in [0,n) for which every signature has a small (k1,k2)."""
    good = []
    for d in range(n):
        ok = True
        for s in sigs:
            k = (s['A'] + s['B'] * d) % n
            if representable(k, n, lam, k1_bound, k2_bound) is None:
                ok = False
                break
        if ok:
            good.append(d)
    return good

def witness_vector_norm(sigs, d, n, lam, k1_bound, k2_bound):
    """||v|| of the lattice vector realising candidate d (original lattice)."""
    m = len(sigs)
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    tot = (d * S_D) ** 2 + S_KAN ** 2
    for s in sigs:
        k = (s['A'] + s['B'] * d) % n
        r = representable(k, n, lam, k1_bound, k2_bound)
        if r is None:
            return None
        k1, k2 = r
        tot += (k1 * S_K1) ** 2 + (k2 * S_K2) ** 2
    return math.sqrt(tot)

# ---------------------------------------------------------------------------
# E6 helper: exact CVP in the last-coordinate-zero sublattice
# ---------------------------------------------------------------------------

def bdd_optimum(sigs, n, lam, k1_bound, k2_bound, projected):
    """Exact shortest vector of the coset {last coordinate = S_KANNAN}.

    Returns (norm, recovered_d or None).  L_0 = rows of the lattice whose last
    coordinate is 0 (all rows but the Kannan row); the Kannan row supplies the
    target t.  We minimise ||t + u|| over u in L_0, i.e. CVP(L_0, -t).
    """
    m = len(sigs)
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    if projected:
        M = build_glv_lattice_projected(sigs, n, lam, k1_bound, k2_bound)
        ncols = 2 * m          # drop the Kannan column
    else:
        M = build_glv_lattice(sigs, n, lam, k1_bound, k2_bound)
        ncols = 2 * m + 1
    body = [row[:ncols] for row in M[:-1]]        # every row but the Kannan row
    t = M[-1][:ncols]                             # the Kannan row's body
    A = IntegerMatrix.from_matrix(body)
    LLL.reduction(A)
    rows = [[A[i][j] for j in range(ncols)] for i in range(A.nrows)]
    rows = [r for r in rows if any(r)]
    B = IntegerMatrix.from_matrix(rows)
    u = CVP.closest_vector(B, tuple(-x for x in t))
    v = [t[j] + u[j] for j in range(ncols)]
    full = math.sqrt(sum(x * x for x in v) + S_KAN ** 2)
    if projected:
        k1s = [v[i] // S_K1 for i in range(m)]
        k2s = [v[m + i] // S_K2 for i in range(m)]
        if any(v[i] % S_K1 for i in range(m)) or any(v[m + i] % S_K2 for i in range(m)):
            return full, None
        if not all(0 <= a < k1_bound for a in k1s) or not all(0 <= b < k2_bound for b in k2s):
            return full, None
        B0 = sigs[0]['B'] % n
        d_c = (k1s[0] + lam * k2s[0] - sigs[0]['A']) * modinv(B0, n) % n
        return full, d_c
    return full, v[m] % n if S_D == 1 else None


# ---------------------------------------------------------------------------
# E7 helper: exact cost of the cheapest coset vector for every candidate d'
# ---------------------------------------------------------------------------

def coset_costs(sigs, n, lam, k1_bound, k2_bound, k2_range=None):
    """cost(d')^2 for every d' in [0,n), exactly (see module docstring, E7)."""
    m = len(sigs)
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    if k2_range is None:
        # |k2_i * S_K2| can never exceed the planted norm ~ n*sqrt(2m/3+1)
        k2_range = int(math.ceil(n * math.sqrt(2 * m / 3.0 + 1) / S_K2)) + 2
    # int64 is exact here: the largest term is (n/2*S_K1)^2 ~ 2e11 for the
    # 12-bit toys, and a sum of 2m+2 such terms stays far below 2^63.
    k2 = np.arange(-k2_range, k2_range + 1, dtype=np.int64)
    lam_k2 = (lam % n) * k2 % n                      # lam*k2 mod n
    cost_k2 = (k2 * S_K2) ** 2
    A = np.array([s['A'] % n for s in sigs], dtype=np.int64)
    B = np.array([s['B'] % n for s in sigs], dtype=np.int64)
    base = S_KAN ** 2
    out = []
    for d in range(n):
        c = (A + B * d) % n                          # per-signature target
        tot = (d * S_D) ** 2 + base
        for ci in c:
            r = (ci - lam_k2) % n
            cent = np.minimum(r, n - r)
            tot += int(((cent * S_K1) ** 2 + cost_k2).min())
        out.append(tot)
    return out


def main():
    print("=" * 78)
    print("GLV-HNP Phase 2, Thread 23b — is the K1 wall information-theoretic?")
    print("=" * 78)

    SEEDS = [42, 1234, 9999, 555, 31337]
    HIST = [
        ("12-bit/2557", 2557, 2, 2659, 1755),
        ("12-bit/2677", 2677, 2, 2647, 185),
    ]
    curves = []
    for label, p, b, n, lam in HIST:
        G = find_generator(p, b, n)
        assert G is not None and (lam * lam + lam + 1) % n == 0
        curves.append((label, (p, b, n, lam, G)))

    M_SIGS = 12
    K1_GRID = [4, 6, 8, 12]

    # -----------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("E5: exhaustive d-uniqueness (ground truth), m=12, 5 seeds")
    print("-" * 78)
    print("#d = number of d' in [0,n) admitting small (k1,k2) for ALL m signatures.")
    print("If #d == 1 on every seed, the instance is unambiguous and any LLL/BKZ")
    print("failure at that (K1,m) is an ALGORITHMIC failure, not a lack of data.\n")
    print(f"{'curve':<13} {'K1':>3} {'eff':>6} {'#d per seed':>20} {'LLL':>6} "
          f"{'pv rank by norm':>16}")

    for label, curve in curves:
        p, b, n, lam, G = curve
        k2b = math.isqrt(n) + 1
        for k1 in K1_GRID:
            counts, wins, ranks = [], 0, []
            for seed in SEEDS:
                d_true = random.Random(seed + 7777).randint(1, n - 1)
                sigs = gen_signatures(G, d_true, M_SIGS, n, lam, p, k1, k2b, seed)
                good = enumerate_d(sigs, n, lam, k1, k2b)
                assert d_true in good, "ground-truth enumeration missed the true d"
                counts.append(len(good))
                # rank of the true d by lattice-vector norm among all admissible d'
                norms = sorted((witness_vector_norm(sigs, d, n, lam, k1, k2b), d)
                               for d in good)
                ranks.append(1 + [d for _, d in norms].index(d_true))
                r = _t23.run_one(curve, M_SIGS, d_true, k1, seed, projected=False)
                wins += bool(r['ok'])
            eff = k1 * k2b / n
            print(f"{label:<13} {k1:>3} {eff:>6.3f} {str(counts):>20} "
                  f"{str(wins)+'/'+str(len(SEEDS)):>6} {str(ranks):>16}")

    # -----------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("E6: exact BDD optimum of the coset {last coord = S_KANNAN}")
    print("-" * 78)
    print("cvp/pv = ||exact BDD optimum|| / ||v_planted||  (<=1 always;")
    print("  == 1 means the planted vector IS the exact optimum).")
    print("d ok  = does the exact optimum decode to the true d?\n")
    print(f"{'curve':<13} {'lat':<5} {'K1':>3} {'seed':>6} {'cvp/pv':>7} "
          f"{'d ok':>5} {'LLL':>5} {'secs':>6}")

    for label, curve in curves:
        p, b, n, lam, G = curve
        k2b = math.isqrt(n) + 1
        for k1 in K1_GRID:
            for seed in SEEDS[:3]:
                d_true = random.Random(seed + 7777).randint(1, n - 1)
                sigs = gen_signatures(G, d_true, M_SIGS, n, lam, p, k1, k2b, seed)
                for projected in (False, True):
                    pv = (planted_vector_projected(sigs, n, k1, k2b) if projected
                          else planted_vector(sigs, d_true, n, k1, k2b))
                    t0 = time.time()
                    cn, d_c = bdd_optimum(sigs, n, lam, k1, k2b, projected)
                    dt = time.time() - t0
                    r = _t23.run_one(curve, M_SIGS, d_true, k1, seed, projected)
                    print(f"{label:<13} {'proj' if projected else 'orig':<5} "
                          f"{k1:>3} {seed:>6} {cn/norm(pv):>7.4f} "
                          f"{str(d_c == d_true):>5} {str(r['ok']):>5} {dt:>6.1f}")

    # -----------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("E7: exact BDD rank of the true d among all n cosets, m=12")
    print("-" * 78)
    print("rank 1 => the planted vector IS the coset optimum (LLL/BKZ can win);")
    print("rank r => an attack must enumerate r coset vectors past the optimum.")
    print("gap = ||coset optimum|| / ||planted||.\n")
    print(f"{'curve':<13} {'K1':>3} {'eff':>6} {'ranks (5 seeds)':>26} "
          f"{'gap':>7} {'LLL':>6} {'secs':>6}")

    for label, curve in curves:
        p, b, n, lam, G = curve
        k2b = math.isqrt(n) + 1
        for k1 in K1_GRID:
            ranks, gaps, wins = [], [], 0
            t0 = time.time()
            for seed in SEEDS:
                d_true = random.Random(seed + 7777).randint(1, n - 1)
                sigs = gen_signatures(G, d_true, M_SIGS, n, lam, p, k1, k2b, seed)
                costs = coset_costs(sigs, n, lam, k1, k2b)
                order = sorted(range(n), key=lambda d: costs[d])
                ranks.append(1 + order.index(d_true))
                gaps.append(math.sqrt(costs[order[0]] / costs[d_true]))
                r = _t23.run_one(curve, M_SIGS, d_true, k1, seed, projected=False)
                wins += bool(r['ok'])
            dt = time.time() - t0
            print(f"{label:<13} {k1:>3} {k1*k2b/n:>6.3f} {str(ranks):>26} "
                  f"{sum(gaps)/len(gaps):>7.4f} "
                  f"{str(wins)+'/'+str(len(SEEDS)):>6} {dt:>6.1f}")

    print("\n" + "=" * 78)
    print("done")
    print("=" * 78)


if __name__ == "__main__":
    main()
