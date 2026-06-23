"""
GLV-HNP Phase 2: Per-curve lattice conditioning analysis.

From 2026-06-22 boundary sweep: pairs at the same λ/n show opposite LLL outcomes.
Hypothesis: condition number κ(M) (or GS orthogonality defect) of raw basis predicts
LLL success better than δ/n alone.

Three focal discordant pairs (same λ/n ≈ target, opposite LLL outcomes):
  PAIR 1 (λ/n≈0.20): lam/n=0.2114,δ/n=0.366 (ALL FAIL) vs lam/n=0.2122,δ/n=0.364 (TRIVIALLY SUCCEED)
  PAIR 2 (λ/n≈0.40): lam/n=0.3867,δ/n=0.160 (ALL FAIL) vs lam/n=0.4110,δ/n=0.233 (SUCCEED m=13)
  PAIR 3 (λ/n≈0.43): lam/n=0.4289,δ/n=0.287 (SUCCEED m=13) vs lam/n=0.4257,δ/n=0.277 (ALL FAIL)

For each curve build basis at m=12 (representative), compute:
  (a) κ₂(M) = σ_max/σ_min (numpy.linalg.cond)
  (b) GS orthogonality defect: log2(∏||b_i||) - log2(vol(L))
  (c) log2 of smallest/largest GS norm ratio
  (d) LLL success rate at m=10..20, 3 seeds

Run: python3 glv_hnp_conditioning.py
"""

import math
import random
import sympy
import numpy as np
from fpylll import IntegerMatrix, LLL

# ---------------------------------------------------------------------------
# EC arithmetic (reused from lamn_boundary.py)
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
    m2, c, t, r = s, pow(z, q, p), pow(n, q, p), pow(n, (q + 1) // 2, p)
    while True:
        if t == 0: return 0
        if t == 1: return r
        i, tmp = 0, t
        while tmp != 1: tmp = tmp * tmp % p; i += 1
        b = pow(c, 1 << (m2 - i - 1), p)
        m2, c, t, r = i, b * b % p, t * b * b % p, r * b % p

def find_generator(p, b_param, n):
    rng = random.Random(12345)
    for _ in range(50000):
        x = rng.randint(0, p - 1)
        rhs = (pow(x, 3, p) + b_param) % p
        y = tonelli_shanks(rhs, p)
        if y is not None and y != 0:
            P = (x, y)
            if ec_mul(P, n, p) is None:
                return P
    return None

def eisenstein_decompose(p):
    for a in range(1, 2 * math.isqrt(p // 3) + 3):
        disc = 4 * p - 3 * a * a
        if disc < 0:
            break
        s = math.isqrt(disc)
        if s * s != disc:
            continue
        for num in [a + s, a - s]:
            if num % 2 == 0:
                b2 = num // 2
                if b2 >= 0 and a * a - a * b2 + b2 * b2 == p:
                    return (a, b2)
    return None

def j0_traces(a, b2):
    return [2*a - b2, -2*a + b2, -(a + b2), a + b2, 2*b2 - a, a - 2*b2]

def glv_eigenvalue(n):
    neg3 = (n - 3) % n
    sq = tonelli_shanks(neg3, n)
    if sq is None:
        return None
    inv2 = modinv(2, n)
    r1 = (n - 1 + sq) * inv2 % n
    r2 = (n - 1 + (n - sq)) * inv2 % n
    if (r1 * r1 + r1 + 1) % n != 0:
        r1, r2 = r2, r1
    if (r1 * r1 + r1 + 1) % n != 0:
        return None
    return min(r1, r2)

def delta_ratio(lam, n):
    x = (3 * lam) % n
    return min(x, n - x) / n

def find_b_for_n(p, n):
    for b_try in range(1, min(p, 500)):
        rng2 = random.Random(99)
        for _ in range(20):
            x = rng2.randint(0, p - 1)
            rhs = (pow(x, 3, p) + b_try) % p
            y = tonelli_shanks(rhs, p)
            if y is not None and y != 0:
                P = (x, y)
                if ec_mul(P, n, p) is None:
                    return b_try
                break
    return None

# ---------------------------------------------------------------------------
# Curve search — replicate lamn_boundary.py find_curves_by_lamn
# ---------------------------------------------------------------------------

def find_curves_by_lamn(p_lo, p_hi, targets, tol=0.015, max_per=2):
    found = {t: [] for t in targets}
    p = sympy.nextprime(p_lo - 1)
    while p < p_hi:
        if p % 3 == 1:
            eis = eisenstein_decompose(p)
            if eis is not None:
                a_e, b_e = eis
                for tr in j0_traces(a_e, b_e):
                    n_cand = p + 1 - tr
                    if n_cand < 2: continue
                    if not sympy.isprime(n_cand): continue
                    if n_cand % 3 != 1: continue
                    lam = glv_eigenvalue(n_cand)
                    if lam is None: continue
                    lam_n = lam / n_cand
                    dr = delta_ratio(lam, n_cand)
                    for t in targets:
                        if abs(lam_n - t) <= tol and len(found[t]) < max_per:
                            b_param = find_b_for_n(p, n_cand)
                            if b_param is None: continue
                            found[t].append((p, b_param, n_cand, lam, lam_n, dr))
                            break
        p = sympy.nextprime(p)
        if all(len(found[t]) >= max_per for t in targets):
            break
    return found

# ---------------------------------------------------------------------------
# Lattice construction + conditioning metrics
# ---------------------------------------------------------------------------

def build_lattice(sigs, n, lam, K1):
    m = len(sigs)
    K2 = math.isqrt(n) + 1
    dim = 2 * m + 2
    S_K1 = n // K1
    S_K2 = max(1, n // K2)
    S_KANNAN = n
    M = [[0] * dim for _ in range(dim)]
    for i in range(m):
        M[i][i] = n * S_K1
    for i in range(m):
        M[m][i] = sigs[i]['B'] * S_K1
    M[m][m] = 1
    for i in range(m):
        M[m + 1 + i][i] = -lam * S_K1
        M[m + 1 + i][m + 1 + i] = S_K2
    for i in range(m):
        M[2*m+1][i] = sigs[i]['A'] * S_K1
    M[2*m+1][dim-1] = S_KANNAN
    return M, S_KANNAN

def gram_schmidt_log_norms(M):
    """Compute log2 of GS norms for basis rows of M."""
    basis = [np.array(row, dtype=np.float64) for row in M]
    gs = []
    log_norms = []
    for b in basis:
        v = b.copy()
        for u, u_norm_sq in gs:
            proj = np.dot(b, u) / u_norm_sq
            v = v - proj * u
        n2 = np.dot(v, v)
        if n2 < 1e-30:
            log_norms.append(-1e10)
            gs.append((v, max(n2, 1e-300)))
        else:
            log_norms.append(0.5 * math.log2(max(n2, 1e-300)))
            gs.append((v, n2))
    return log_norms

def lattice_metrics(M_int, n, m):
    """
    Returns dict with:
      kappa2      : 2-norm condition number of M (float, may be large)
      log2_kappa  : log2 of condition number
      gs_defect   : log2 orthogonality defect = sum(log2||b_i||) - log2(vol(L))
      gs_ratio    : log2(gs_norm_max/gs_norm_min) — measures GS telescoping
    All computed in float64; sufficient for 20-bit problem.
    """
    dim = len(M_int)
    A = np.array(M_int, dtype=np.float64)

    # Condition number via SVD (clipped to avoid inf)
    sv = np.linalg.svd(A, compute_uv=False)
    sv_max = sv[0]
    sv_min = sv[-1]
    if sv_min < 1e-10:
        kappa2 = 1e18
        log2_k = 60.0
    else:
        kappa2 = sv_max / sv_min
        log2_k = math.log2(kappa2)

    # log2 of each row norm (||b_i||)
    row_log_norms = [0.5 * math.log2(max(sum(x*x for x in row), 1e-300))
                     for row in M_int]
    sum_row_log = sum(row_log_norms)

    # log2 vol(L) = log2|det(M)| = sum of log2|sv_i|
    sv_nz = sv[sv > 1e-15]
    log2_vol = sum(math.log2(s) for s in sv_nz)

    gs_defect = sum_row_log - log2_vol

    # GS norms
    gs_log = gram_schmidt_log_norms(M_int)
    gs_finite = [x for x in gs_log if x > -1e9]
    gs_ratio = (max(gs_finite) - min(gs_finite)) if gs_finite else 0.0

    return {
        'kappa2': kappa2,
        'log2_kappa': log2_k,
        'gs_defect': gs_defect,
        'gs_ratio': gs_ratio,
    }

# ---------------------------------------------------------------------------
# Signature generation + LLL
# ---------------------------------------------------------------------------

def gen_signatures(p, b_param, n, lam, G, m_sigs, K1, seed):
    K2 = math.isqrt(n) + 1
    rng = random.Random(seed)
    d_secret = rng.randint(1, n - 1)
    sigs = []
    attempts = 0
    while len(sigs) < m_sigs and attempts < 200000:
        attempts += 1
        k1 = rng.randint(0, K1 - 1)
        k2 = rng.randint(0, K2 - 1)
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
        sigs.append({'A': A, 'B': B})
    return d_secret, sigs

def try_recover(M_red, m, n, S_KANNAN, d_secret):
    dim = 2 * m + 2
    for row in M_red:
        last = row[dim - 1]
        if abs(last) != S_KANNAN: continue
        sign = 1 if last > 0 else -1
        d_cand = (sign * row[m]) % n
        if d_cand == d_secret:
            return True
    return False

def run_lll_trial(p, b_param, n, lam, G, m, K1, seed):
    d_secret, sigs = gen_signatures(p, b_param, n, lam, G, m, K1, seed)
    if len(sigs) < m:
        return False
    M_int, S_KANNAN = build_lattice(sigs, n, lam, K1)
    dim = 2 * m + 2
    A = IntegerMatrix.from_matrix(M_int)
    LLL.reduction(A)
    reduced = [[A[i][j] for j in range(dim)] for i in range(dim)]
    return try_recover(reduced, m, n, S_KANNAN, d_secret)

def analyze_curve(label, p, b_param, n, lam, lam_n, dr, K1=72, m_probe=12,
                  m_range=range(10, 23), seeds=(0xDEAD, 0xBEEF, 0xCAFE)):
    """Full analysis: conditioning metrics + LLL sweep."""
    G = find_generator(p, b_param, n)
    if G is None:
        print(f"  [{label}] FAIL: no generator found")
        return None

    # Conditioning at m_probe
    d_secret, sigs = gen_signatures(p, b_param, n, lam, G, m_probe, K1, seed=0xDEAD)
    M_int, _ = build_lattice(sigs, n, lam, K1)
    metrics = lattice_metrics(M_int, n, m_probe)

    # LLL sweep
    first_3of3 = None
    sweep = {}
    for m in m_range:
        wins = 0
        for seed in seeds:
            if run_lll_trial(p, b_param, n, lam, G, m, K1, seed):
                wins += 1
        sweep[m] = wins
        if wins == 3 and first_3of3 is None:
            first_3of3 = m

    print(f"\n  [{label}]  p={p}, n={n}, lam/n={lam_n:.4f}, δ/n={dr:.4f}")
    print(f"    log2_κ(M) = {metrics['log2_kappa']:.2f}  "
          f"GS_defect = {metrics['gs_defect']:.2f}  "
          f"GS_ratio = {metrics['gs_ratio']:.2f}")
    print(f"    LLL sweep (3 seeds): " +
          " ".join(f"m{m}:{sweep[m]}" for m in m_range))
    print(f"    First 3/3: {'m='+str(first_3of3) if first_3of3 else 'NEVER'}")

    return {
        'label': label, 'p': p, 'n': n,
        'lam_n': lam_n, 'dr': dr,
        'log2_kappa': metrics['log2_kappa'],
        'gs_defect': metrics['gs_defect'],
        'gs_ratio': metrics['gs_ratio'],
        'first_3of3': first_3of3,
        'sweep': sweep,
    }

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == '__main__':
    K1 = 72

    # From June 22 log: discordant pairs at λ/n≈0.20, ≈0.40, ≈0.43.
    # We use same targets + tol as lamn_boundary.py to regenerate the same curves.
    # PAIR 1: lam/n≈0.20 → one at 0.2114 (fail) and one at 0.2122 (succeed)
    # PAIR 2: lam/n≈0.40 → one at 0.387 (fail) and one at 0.411 (succeed)
    # PAIR 3: lam/n≈0.43 → one at 0.429 (succeed) and one at 0.426 (fail)
    # We widen search to find at least 2 curves near each discordant λ/n.

    print("=== GLV-HNP Conditioning Analysis 2026-06-23 ===")
    print("Searching for discordant-pair curves (may take ~2 min)...")

    TARGETS = [0.20, 0.40, 0.43]
    curves_by_target = find_curves_by_lamn(2**19, 2**20, TARGETS, tol=0.015, max_per=2)

    results = []
    for t in TARGETS:
        group = curves_by_target[t]
        print(f"\n--- λ/n ≈ {t:.2f}: {len(group)} curve(s) found ---")
        for idx, (p, b_param, n, lam, lam_n, dr) in enumerate(group):
            label = f"λ/n≈{t:.2f} curve_{idx+1}"
            r = analyze_curve(label, p, b_param, n, lam, lam_n, dr,
                              K1=K1, m_probe=12, m_range=range(10, 23))
            if r is not None:
                results.append(r)

    # Summary table
    print("\n\n=== SUMMARY TABLE ===")
    print(f"{'Label':<30} {'lam/n':>7} {'δ/n':>7} {'log2_κ':>9} {'GS_def':>8} {'GS_rat':>8} {'1st_3/3':>9}")
    print("-" * 90)
    for r in results:
        f33 = str(r['first_3of3']) if r['first_3of3'] else 'NEVER'
        print(f"{r['label']:<30} {r['lam_n']:7.4f} {r['dr']:7.4f} "
              f"{r['log2_kappa']:9.2f} {r['gs_defect']:8.2f} {r['gs_ratio']:8.2f} {f33:>9}")

    # Correlation check: within each pair, does log2_κ flip between success/failure?
    print("\n=== PAIR-WISE CONDITIONING COMPARISON ===")
    for t in TARGETS:
        group_r = [r for r in results if abs(r['lam_n'] - t) < 0.025]
        if len(group_r) < 2:
            print(f"λ/n≈{t:.2f}: only {len(group_r)} curve(s), skip.")
            continue
        g0, g1 = group_r[0], group_r[1]
        s0 = "SUCCESS" if g0['first_3of3'] else "FAIL"
        s1 = "SUCCESS" if g1['first_3of3'] else "FAIL"
        k_ratio = g0['log2_kappa'] - g1['log2_kappa']
        print(f"λ/n≈{t:.2f}:")
        print(f"  Curve 1: {s0}, log2_κ={g0['log2_kappa']:.2f}, GS_def={g0['gs_defect']:.2f}")
        print(f"  Curve 2: {s1}, log2_κ={g1['log2_kappa']:.2f}, GS_def={g1['gs_defect']:.2f}")
        print(f"  Δlog2_κ = {k_ratio:.2f} (positive → curve 1 more ill-conditioned)")
        if s0 != s1:
            winner = "curve_1_worse" if k_ratio > 0 else "curve_2_worse"
            consistent = (k_ratio > 0 and s0 == "FAIL") or (k_ratio < 0 and s1 == "FAIL")
            print(f"  DISCORDANT PAIR: κ predicts failure = {winner}. Consistent with LLL? {consistent}")
        else:
            print(f"  Both {s0} — no discordance at this pair (need more seeds or wider m)")

    print("\nDone.")
