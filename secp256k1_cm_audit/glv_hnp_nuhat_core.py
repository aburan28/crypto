"""
RECOVERY NOTE (2026-08-05, autolab)
-----------------------------------
This file is the Thread 20a/20b/20c/20d module that was committed on 2026-07-29
as `secp256k1_cm_audit/glv_hnp_phase2_lambda_threshold.py` by commit e845207.

Two autolab sessions ran independently on 2026-07-29 (d525931 and e845207) and
each wrote a DIFFERENT 500-line file under that same name.  The merge that
brought both onto main (7b5b702 / e31dd19, PR #33) kept d525931's version and
silently discarded this one, which left `glv_hnp_phase2_mu_response.py`,
`glv_hnp_phase2_nuhat_control.py` and `glv_hnp_nuhat_vs_c1c2.py` importing
names (`glv_eigenvalues`, `mu_of`, `identify_twist`, `rival_sublattice_nu`)
that no longer existed -- all three died at import.

Restored here under a non-colliding name; the three dependents now import from
this file.  d525931's file keeps the original name and is unchanged.
Content below is byte-identical to `e845207:.../glv_hnp_phase2_lambda_threshold.py`.
"""

"""
GLV-HNP Phase 2 — Thread 20: lambda/n threshold study for the GLV-aware lattice.

Context
-------
2026-07-26 (run #1) validated the Phase 2 (2m+2)-dimensional column-scaled GLV
lattice and observed:

    lam/n = 0.53, 0.66, 0.34  ->  LLL reaches 3/3
    lam/n = 0.07             ->  LLL and BKZ(40) both fail at every m

and proposed bisecting the threshold in (0.07, 0.34).  This script does that.

Two corrections to the 2026-07-26 framing are applied here:

(1) lam/n is NOT the right coordinate.  For j=0 curves the two GLV eigenvalues
    are lam and lam' = n-1-lam, so lam/n and (n-lam)/n describe the *same*
    curve.  The curve invariant is the symmetric

        mu = min(lam, n - lam) / n   in (0, 1/2].

    Under mu, the 2026-07-26 data reads: mu = 0.467, 0.339, 0.34 succeed,
    mu = 0.070 fails.  The bisection interval is mu in (0.07, 0.34).

(2) mu (equivalently lam/n) was ALREADY falsified as a predictor in the
    *Phase 1* lattice: Exp M/N/O (2026-06-27) found overlapping ranges, and
    Exp S (2026-06-29) recorded n=236209 with lam/n=0.018 succeeding 6/6.
    So the pre-registered expectation here is that the mu-threshold does NOT
    exist in Phase 2 either.

Primary hypothesis (H-nu)
-------------------------
The 2026-06-29 log entry closed with an unexecuted proposal:

    "the separator likely requires a lattice-geometric computation:
     specifically, whether the BV lattice for (n, lam) admits a short
     non-planted vector whose norm is less than sqrt(m) * K1 * S_K1"

In the Phase 2 lattice this is directly computable.  Rows i (0<=i<m) and
row m+1+i are supported on the coordinate pair (i, m+1+i):

    u = (n * S_K1, 0)          [row i]
    v = (-lam * S_K1, S_K2)    [row m+1+i]

They generate a 2-dimensional sublattice L2 whose shortest vector

    nu = lambda_1(L2)

is a *non-planted* short vector (last coordinate 0, so recover_d can never
read d off it).  det(L2) = n * S_K1 * S_K2 does not depend on lam, so

    nu_hat = nu / sqrt(det L2)

isolates exactly the lam-dependence of the lattice geometry.  There are m
independent copies of L2 (one per signature index), so the planted vector
competes against >= 2m vectors of norm ~nu.

H-nu:  LLL recovers d iff nu is not too small relative to the planted norm,
       i.e. the separator is R = nu / ||planted||, not mu.

Falsifier: if success/failure ranges of R overlap as badly as those of mu,
H-nu joins the seven already-falsified invariants (delta/n, kappa(M), q_cf,
max_q_cf, max_a, a_corn/n, mu) and the Phase 2 separator is likewise
non-algebraic.

Run: python3 glv_hnp_phase2_lambda_threshold.py
"""

import math
import random
import sys
import time

import sympy
from fpylll import IntegerMatrix, LLL

# ---------------------------------------------------------------------------
# Experiment configuration
# ---------------------------------------------------------------------------

BITS_LO, BITS_HI = 2 ** 19, 2 ** 20      # 20-bit primes (matches 2026-07-26)
EFF_TARGET = 0.05                        # K1 * K2 / n, gives m_thresh ~ 5
SEEDS = [42, 1234, 9999]                 # 3 seeds, "3/3" convention
M_RANGE = range(4, 15)                   # m = 4..14
MU_TARGETS = [0.03, 0.05, 0.07, 0.10, 0.14, 0.18, 0.22,
              0.26, 0.30, 0.34, 0.38, 0.42, 0.46, 0.49]
MU_TOL = 0.010
CURVES_PER_BUCKET = 2
MAX_PRIMES_SCANNED = 12000

# ---------------------------------------------------------------------------
# Minimal EC arithmetic (short Weierstrass, a=0) -- same core as
# glv_hnp_phase2_20bit.py so results are directly comparable.
# ---------------------------------------------------------------------------

def modinv(a, m):
    return pow(a, -1, m)

def ec_add(P, Q, p):
    if P is None: return Q
    if Q is None: return P
    x1, y1 = P; x2, y2 = Q
    if x1 == x2:
        if (y1 + y2) % p == 0: return None
        s = 3 * x1 * x1 * modinv(2 * y1, p) % p
    else:
        s = (y2 - y1) * modinv(x2 - x1, p) % p
    x3 = (s * s - x1 - x2) % p
    y3 = (s * (x1 - x3) - y1) % p
    return (x3, y3)

def ec_mul(P, k, p):
    if k == 0: return None
    R, Q = None, P
    while k > 0:
        if k & 1: R = ec_add(R, Q, p)
        Q = ec_add(Q, Q, p)
        k >>= 1
    return R

def tonelli_shanks(n, p):
    n %= p
    if n == 0: return 0
    if pow(n, (p - 1) // 2, p) != 1: return None
    if p % 4 == 3: return pow(n, (p + 1) // 4, p)
    q, s = p - 1, 0
    while q % 2 == 0: q //= 2; s += 1
    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1: z += 1
    m, c, t, r = s, pow(z, q, p), pow(n, q, p), pow(n, (q + 1) // 2, p)
    while True:
        if t == 0: return 0
        if t == 1: return r
        i, tmp = 0, t
        while tmp != 1: tmp = tmp * tmp % p; i += 1
        b = pow(c, 1 << (m - i - 1), p)
        m, c, t, r = i, b * b % p, t * b * b % p, r * b % p

def find_point(p, b, seed=42):
    rng = random.Random(seed)
    for _ in range(2000):
        x = rng.randint(0, p - 1)
        rhs = (pow(x, 3, p) + b) % p
        y = tonelli_shanks(rhs, p)
        if y is not None and y != 0:
            return (x, y)
    return None

def find_generator(p, b, n, seed=12345):
    """A point of exact order n (n prime) is a generator of E(F_p) with #E=n."""
    rng = random.Random(seed)
    for _ in range(2000):
        x = rng.randint(0, p - 1)
        rhs = (pow(x, 3, p) + b) % p
        y = tonelli_shanks(rhs, p)
        if y is not None and y != 0:
            P = (x, y)
            if ec_mul(P, n, p) is None:
                return P
    return None

# ---------------------------------------------------------------------------
# CM theory: Eisenstein decomposition for j=0 curves over F_p
# ---------------------------------------------------------------------------

def eisenstein_decompose(p):
    """(a,b) >= 0 with a^2 - a*b + b^2 = p, in O(sqrt(p))."""
    for a in range(1, 2 * math.isqrt(p // 3) + 3):
        disc = 4 * p - 3 * a * a
        if disc < 0:
            break
        s = math.isqrt(disc)
        if s * s != disc:
            continue
        for num in (a + s, a - s):
            if num % 2 == 0:
                b = num // 2
                if b >= 0 and a * a - a * b + b * b == p:
                    return (a, b)
    return None

def j0_traces(a, b):
    """The 6 Frobenius traces of the 6 sextic twists of j=0 over F_p."""
    return [2 * a - b, -2 * a + b, -(a + b), a + b, 2 * b - a, a - 2 * b]

def glv_eigenvalues(n):
    """Both roots of x^2 + x + 1 = 0 mod n (requires n = 1 mod 3)."""
    sq = tonelli_shanks((n - 3) % n, n)
    if sq is None:
        return None
    inv2 = modinv(2, n)
    r1 = (n - 1 + sq) * inv2 % n
    r2 = (n - 1 - sq) * inv2 % n
    if (r1 * r1 + r1 + 1) % n != 0 or (r2 * r2 + r2 + 1) % n != 0:
        return None
    return (min(r1, r2), max(r1, r2))

def mu_of(lam, n):
    """Symmetric curve invariant mu = min(lam, n-lam)/n in (0, 1/2]."""
    return min(lam, n - lam) / n

def identify_twist(p, n):
    """Smallest b with #(y^2 = x^3 + b over F_p) = n."""
    for b in range(1, 400):
        P = find_point(p, b)
        if P is None:
            continue
        if ec_mul(P, n, p) is None:
            return b
    return None

# ---------------------------------------------------------------------------
# Lattice geometry
# ---------------------------------------------------------------------------

def gauss_reduce_2d(u, v):
    """Lagrange-Gauss reduction; returns the shortest nonzero vector norm."""
    def nrm2(w):
        return w[0] * w[0] + w[1] * w[1]
    if nrm2(u) < nrm2(v):
        u, v = v, u
    while True:
        q = round((u[0] * v[0] + u[1] * v[1]) / nrm2(v))
        r = (u[0] - q * v[0], u[1] - q * v[1])
        u, v = v, r
        if nrm2(u) <= nrm2(v):
            return math.sqrt(nrm2(u))

def scales(n, k1_bound, k2_bound):
    """Column-diagonal scaling validated 2026-07-26."""
    return (n // k1_bound, 1, max(1, n // k2_bound), n)

def rival_sublattice_nu(n, lam, k1_bound, k2_bound):
    """
    nu = lambda_1 of the non-planted 2D sublattice L2 spanned by
      u = (n*S_K1, 0)  and  v = (-lam*S_K1, S_K2),
    together with det(L2) and the scale-free nu_hat = nu / sqrt(det).
    """
    s_k1, _, s_k2, _ = scales(n, k1_bound, k2_bound)
    nu = gauss_reduce_2d((n * s_k1, 0), (-lam * s_k1, s_k2))
    det = n * s_k1 * s_k2
    return nu, det, nu / math.sqrt(det)

def planted_norm_expected(n, m, k1_bound, k2_bound):
    """
    E[||planted||] with k1 ~ U[0,K1), k2 ~ U[0,K2), d ~ U[1,n).
    Components: k1_i*S_K1, d*S_D, k2_i*S_K2, S_KANNAN.
    """
    s_k1, s_d, s_k2, s_kan = scales(n, k1_bound, k2_bound)
    e_k1sq = (k1_bound ** 2) / 3.0
    e_k2sq = (k2_bound ** 2) / 3.0
    e_dsq = (n ** 2) / 3.0
    return math.sqrt(m * e_k1sq * s_k1 ** 2 + e_dsq * s_d ** 2
                     + m * e_k2sq * s_k2 ** 2 + s_kan ** 2)

def m_threshold(n, k1_bound, k2_bound):
    """Information-theoretic m at which the planted vector becomes unique."""
    eff = k1_bound * k2_bound / n
    if eff >= 1.0:
        return float('inf')
    return math.ceil(math.log(n) / math.log(1.0 / eff))

# ---------------------------------------------------------------------------
# Phase 2 attack
# ---------------------------------------------------------------------------

def gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed):
    rng = random.Random(seed)
    sigs = []
    attempts = 0
    while len(sigs) < m and attempts < 100000:
        attempts += 1
        k1 = rng.randint(0, k1_bound - 1)
        k2 = rng.randint(0, k2_bound - 1)
        k_full = (k1 + lam * k2) % n
        if k_full == 0: continue
        R = ec_mul(G, k_full, p)
        if R is None: continue
        r = R[0] % n
        if r == 0: continue
        h = rng.randint(0, n - 1)
        s = modinv(k_full, n) * (h + d_secret * r) % n
        if s == 0: continue
        s_inv = modinv(s, n)
        A = h * s_inv % n
        B = r * s_inv % n
        assert (A + B * d_secret) % n == k_full
        sigs.append({'A': A, 'B': B, 'k1': k1, 'k2': k2})
    return sigs

def build_glv_lattice(sigs, n, lam, k1_bound, k2_bound):
    m = len(sigs)
    dim = 2 * m + 2
    s_k1, s_d, s_k2, s_kan = scales(n, k1_bound, k2_bound)
    M = [[0] * dim for _ in range(dim)]
    for i in range(m):
        M[i][i] = n * s_k1
    for i in range(m):
        M[m][i] = sigs[i]['B'] * s_k1
    M[m][m] = s_d
    for i in range(m):
        M[m + 1 + i][i] = -lam * s_k1
        M[m + 1 + i][m + 1 + i] = s_k2
    for i in range(m):
        M[2 * m + 1][i] = sigs[i]['A'] * s_k1
    M[2 * m + 1][dim - 1] = s_kan
    return M, s_kan

def run_experiment(p, n, lam, G, m, d_secret, k1_bound, k2_bound, seed):
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return False
    M, s_kan = build_glv_lattice(sigs, n, lam, k1_bound, k2_bound)
    dim = 2 * m + 2
    A = IntegerMatrix.from_matrix(M)
    LLL.reduction(A)
    for i in range(dim):
        last = A[i][dim - 1]
        if abs(last) != s_kan:
            continue
        sign = 1 if last > 0 else -1
        if (sign * A[i][m]) % n == d_secret:
            return True
    return False

def sweep_curve(p, n, lam, G, k1_bound, k2_bound, m_range=M_RANGE, seeds=SEEDS):
    """Returns (wins_by_m, first_3of3, best_wins)."""
    wins_by_m = {}
    first = None
    best = 0
    for m in m_range:
        w = 0
        for seed in seeds:
            d_trial = random.Random(seed + 7777).randint(1, n - 1)
            w += run_experiment(p, n, lam, G, m, d_trial, k1_bound, k2_bound, seed)
        wins_by_m[m] = w
        best = max(best, w)
        if w == len(seeds) and first is None:
            first = m
    return wins_by_m, first, best

# ---------------------------------------------------------------------------
# Curve harvesting, bucketed by mu
# ---------------------------------------------------------------------------

def harvest_curves():
    buckets = {t: [] for t in MU_TARGETS}
    seen_n = set()
    scanned = 0
    p = int(sympy.nextprime(BITS_LO))
    t0 = time.time()
    while p < BITS_HI and scanned < MAX_PRIMES_SCANNED:
        scanned += 1
        if p % 3 == 1:
            eis = eisenstein_decompose(p)
            if eis is not None:
                for t in j0_traces(*eis):
                    n = p + 1 - t
                    if n < 2 or n in seen_n or not sympy.isprime(n) or n % 3 != 1:
                        continue
                    roots = glv_eigenvalues(n)
                    if roots is None:
                        continue
                    lam = roots[0]
                    mu = mu_of(lam, n)
                    for target in MU_TARGETS:
                        if abs(mu - target) <= MU_TOL and len(buckets[target]) < CURVES_PER_BUCKET:
                            b = identify_twist(p, n)
                            if b is None:
                                break
                            G = find_generator(p, b, n)
                            if G is None:
                                break
                            seen_n.add(n)
                            buckets[target].append(
                                {'p': p, 'b': b, 'n': n, 'lam': lam,
                                 'lam2': roots[1], 'mu': mu, 'target': target, 'G': G})
                            print(f"  [bucket {target:.2f}] p={p} n={n} lam={lam} "
                                  f"mu={mu:.4f}", flush=True)
                            break
        if all(len(v) >= CURVES_PER_BUCKET for v in buckets.values()):
            break
        p = int(sympy.nextprime(p))
    print(f"  scanned {scanned} primes in {time.time()-t0:.1f}s")
    return [c for t in MU_TARGETS for c in buckets[t]]

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("=" * 78)
    print("GLV-HNP Phase 2 — Thread 20: mu = min(lam, n-lam)/n threshold study")
    print("=" * 78)
    print(f"20-bit j=0 CM curves, eff target={EFF_TARGET}, seeds={SEEDS}, "
          f"m in {M_RANGE.start}..{M_RANGE.stop-1}")

    print("\n--- Part A: harvest curves bucketed by mu ---")
    curves = harvest_curves()
    print(f"  harvested {len(curves)} curves")

    print("\n--- Part B: Phase 2 LLL sweep per curve ---")
    rows = []
    t0 = time.time()
    for c in curves:
        n, lam = c['n'], c['lam']
        k2 = math.isqrt(n) + 1
        k1 = max(2, int(EFF_TARGET * math.sqrt(n)))
        nu, det2, nu_hat = rival_sublattice_nu(n, lam, k1, k2)
        wins, first, best = sweep_curve(c['p'], n, lam, c['G'], k1, k2)
        m_ref = first if first is not None else M_RANGE.stop - 1
        planted = planted_norm_expected(n, m_ref, k1, k2)
        row = {**c, 'k1': k1, 'k2': k2, 'nu': nu, 'nu_hat': nu_hat,
               'planted': planted, 'ratio': nu / planted,
               'wins': wins, 'first': first, 'best': best,
               'm_thresh': m_threshold(n, k1, k2)}
        rows.append(row)
        status = f"3/3 at m={first}" if first else f"never 3/3 (best {best}/3)"
        print(f"  n={n:<8} mu={c['mu']:.4f} K1={k1} nu_hat={nu_hat:.4f} "
              f"nu/planted={row['ratio']:.3f}  {status}", flush=True)
    print(f"  sweep wall time {time.time()-t0:.1f}s")

    # ---- Part C: does a mu threshold exist? -------------------------------
    print("\n" + "=" * 78)
    print("Part C: mu threshold")
    print("=" * 78)
    succ = [r for r in rows if r['first'] is not None]
    fail = [r for r in rows if r['first'] is None]
    print(f"  success: {len(succ)}/{len(rows)}   failure: {len(fail)}/{len(rows)}")
    if succ:
        print(f"  mu range (success): [{min(r['mu'] for r in succ):.4f}, "
              f"{max(r['mu'] for r in succ):.4f}]")
    if fail:
        print(f"  mu range (failure): [{min(r['mu'] for r in fail):.4f}, "
              f"{max(r['mu'] for r in fail):.4f}]")
    if succ and fail:
        lo = min(r['mu'] for r in succ)
        hi = max(r['mu'] for r in fail)
        if lo > hi:
            print(f"  CLEAN THRESHOLD: mu* in ({hi:.4f}, {lo:.4f})")
        else:
            n_over = sum(1 for r in succ if r['mu'] <= hi) + \
                     sum(1 for r in fail if r['mu'] >= lo)
            print(f"  NO CLEAN THRESHOLD: ranges overlap on "
                  f"[{lo:.4f}, {hi:.4f}] ({n_over} curves in the overlap)")

    # ---- Part C2: H-nu, the lattice-geometric predictor -------------------
    print("\n" + "=" * 78)
    print("Part C2: H-nu  (nu = lambda_1 of non-planted 2D sublattice)")
    print("=" * 78)
    for key in ('nu_hat', 'ratio'):
        if succ:
            print(f"  {key:7s} success: [{min(r[key] for r in succ):.4f}, "
                  f"{max(r[key] for r in succ):.4f}]")
        if fail:
            print(f"  {key:7s} failure: [{min(r[key] for r in fail):.4f}, "
                  f"{max(r[key] for r in fail):.4f}]")
        if succ and fail:
            lo = min(r[key] for r in succ)
            hi = max(r[key] for r in fail)
            print("    -> CLEAN SEPARATION" if lo > hi else
                  f"    -> OVERLAP on [{lo:.4f}, {hi:.4f}]")

    # best single-cut accuracy for each candidate predictor
    print("\n  best single-cut accuracy (predict success iff x > cut):")
    for key in ('mu', 'nu_hat', 'ratio'):
        cuts = sorted(r[key] for r in rows)
        best_acc, best_cut = 0.0, None
        for cut in cuts:
            acc = sum(1 for r in rows
                      if (r[key] > cut) == (r['first'] is not None)) / len(rows)
            if acc > best_acc:
                best_acc, best_cut = acc, cut
        base = max(len(succ), len(fail)) / len(rows)
        print(f"    {key:7s}: {best_acc*100:5.1f}%  (cut={best_cut:.4f}) "
              f"vs trivial baseline {base*100:.1f}%")

    # ---- Part D: full table ----------------------------------------------
    print("\n" + "=" * 78)
    print("Part D: full table")
    print("=" * 78)
    print(f"{'n':>8} {'mu':>7} {'K1':>4} {'nu_hat':>7} {'nu/pl':>7} "
          f"{'m_th':>5} {'first3/3':>9} {'best':>5}  wins by m")
    for r in sorted(rows, key=lambda x: x['mu']):
        wv = ''.join(str(r['wins'][m]) for m in M_RANGE)
        print(f"{r['n']:>8} {r['mu']:>7.4f} {r['k1']:>4} {r['nu_hat']:>7.4f} "
              f"{r['ratio']:>7.3f} {r['m_thresh']:>5} "
              f"{(str(r['first']) if r['first'] else '-'):>9} {r['best']:>5}  {wv}")

    # ---- Part E: root-choice control -------------------------------------
    print("\n" + "=" * 78)
    print("Part E: root-choice control (lam vs lam' = n-1-lam; mu identical)")
    print("=" * 78)
    for r in rows[:6]:
        n, k1, k2 = r['n'], r['k1'], r['k2']
        _, _, nh2 = rival_sublattice_nu(n, r['lam2'], k1, k2)
        _, f2, b2 = sweep_curve(r['p'], n, r['lam2'], r['G'], k1, k2)
        print(f"  n={n:<8} mu={r['mu']:.4f} | lam : nu_hat={r['nu_hat']:.4f} "
              f"first={r['first']} best={r['best']}/3"
              f"  || lam': nu_hat={nh2:.4f} first={f2} best={b2}/3", flush=True)

    print("\nDone.")

if __name__ == '__main__':
    sys.exit(main())
