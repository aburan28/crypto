"""
Thread 23 addendum — why do 2 of 10 17-bit curves recover d far past the
gamma_GH threshold, and only in the UNprojected lattice?

E3 of glv_hnp_thread23_projected.py: p=65539 (lam*=0.335) and p=66853
(lam*=0.482) recover 5/5 at eff up to 0.20-0.25 in L while every other curve
collapses by eff=0.15 -- and they recover 0/5 in pi(L) even at eff=0.05, where
all eight others get 5/5.  Excluding these two, L and pi(L) are identical at
every eff, so the entire "projection hurts" effect is these two curves.

Hypothesis H: the recovering row is NOT the planted vector but
    w = v_planted + z,   z in the sublattice with zero d- and zero Kannan
                         coordinate,
so its d-column still reads d_secret (d is only defined mod n) while its
k-coordinates are wrong.  Deleting the d-column destroys that route, which is
why pi(L) fails on exactly these curves.  If true, the full lattice is strictly
stronger than pi(L) and the 2026-07-29 Thread-23 premise is inverted.

Run: python3 glv_hnp_thread23_outliers.py
"""

import math
import os
import random

from fpylll import IntegerMatrix, LLL

_HERE = os.path.dirname(os.path.abspath(__file__))
_SRC = os.path.join(_HERE, "glv_hnp_phase2_lambda_threshold.py")
with open(_SRC) as _f:
    _text = _f.read()
exec(compile(_text.split("# EXPERIMENT T1")[0], _SRC, "exec"), globals())

SEEDS = [42, 1234, 9999, 555, 31337]


def lll_rows(M, ncols):
    A = IntegerMatrix.from_matrix(M)
    LLL.reduction(A)
    return [[A[i][j] for j in range(ncols)] for i in range(A.nrows)]


def glv_small(k, n, lam, K1, K2):
    for k2 in range(K2):
        if (k - lam * k2) % n < K1:
            return True
    return False


print("=" * 78)
print("Thread 23 addendum — anatomy of the two outlier curves")
print("=" * 78)

curves = search_curves(2 ** 16, 2 ** 17, per_bin=1, nbins=10, max_primes=100000)[:10]
M_SIGS = 12

print(f"\n{'p':>8} {'n':>8} {'lam':>8} {'lam*':>6} {'nu_hat':>7} {'n%3':>4} "
      f"{'K2':>5} {'ord(G)=n':>9}")
for (p, b, n, lam, G) in curves:
    k2b = math.isqrt(n) + 1
    S_K1, S_D, S_K2, S_KAN = scales(n, 63, k2b)
    mu, _ = lambda_block_mu(n, lam, S_K1, S_K2)
    nuhat = mu / math.sqrt(n * S_K1 * S_K2)
    ordok = ec_mul(G, n, p) is None and ec_mul(G, 1, p) is not None
    print(f"{p:>8} {n:>8} {lam:>8} {lam_star(lam, n):>6.3f} {nuhat:>7.3f} "
          f"{n % 3:>4} {k2b:>5} {str(ordok):>9}")

print("\n" + "=" * 78)
print("Anatomy of the recovering row at eff=0.25, m=12 (full lattice)")
print("=" * 78)
print(f"{'p':>8} {'seed':>6} {'route':>9} {'||w||/||pv||':>13} {'||z||/||pv||':>13} "
      f"{'z=0':>5} {'z d-col':>9} {'z kan':>7} {'d ok':>5}")

for (p, b, n, lam, G) in curves:
    k2b = math.isqrt(n) + 1
    K1 = max(2, int(0.25 * n / k2b))
    S_K1, S_D, S_K2, S_KAN = scales(n, K1, k2b)
    for seed in SEEDS[:2]:
        rng = random.Random(seed ^ 0xC0FFEE)
        d_secret = rng.randrange(1, n)
        sigs = gen_signatures(G, d_secret, M_SIGS, n, lam, p, K1, k2b, seed)
        if len(sigs) < M_SIGS:
            continue
        pv = planted_vector(sigs, d_secret, n, K1, k2b)
        red = lll_rows(build_glv_lattice(sigs, n, lam, K1, k2b), 2 * M_SIGS + 2)
        hit = None
        for row in red:
            last = row[2 * M_SIGS + 1]
            if abs(last) != S_KAN:
                continue
            w = row if last > 0 else [-x for x in row]
            if w[M_SIGS] % n == d_secret and glv_small(
                    (sigs[0]['A'] + sigs[0]['B'] * (w[M_SIGS] % n)) % n,
                    n, lam, K1, k2b):
                # is it the planted vector itself, or a shifted cousin?
                exact = all(w[i] == pv[i] for i in range(M_SIGS)) and \
                        all(w[M_SIGS + 1 + i] == pv[M_SIGS + 1 + i]
                            for i in range(M_SIGS))
                hit = (w, exact)
                break
        if hit is None:
            print(f"{p:>8} {seed:>6} {'none':>9} {'-':>13} {'-':>13} "
                  f"{'-':>5} {'-':>9} {'-':>7} {str(False):>5}")
            continue
        w, exact = hit
        z = [w[i] - pv[i] for i in range(len(pv))]
        route = "planted" if exact else "d-col"
        print(f"{p:>8} {seed:>6} {route:>9} {norm(w)/norm(pv):>13.3f} "
              f"{norm(z)/norm(pv):>13.3f} {str(not any(z)):>5} "
              f"{z[M_SIGS] % n:>9} {z[2*M_SIGS+1]:>7} {str(True):>5}")

print("\nroute='planted': LLL found v_planted itself (k-coords exact) -- pi(L)")
print("would recover too.  route='d-col': w = v_planted + z with z != 0 and")
print("z's d-column a multiple of n, so only the d-column reading works and")
print("pi(L) necessarily fails.  If the two outliers are all 'd-col' and the")
print("other eight are all 'planted'/'none', hypothesis H is confirmed.")
print("\nDone.")
