"""
Thread 23 — reformulate the GLV-HNP Phase-2 lattice so the planted vector can be
the shortest vector.

Background (RESEARCH_AUTOLAB_LOG.md, 2026-07-29, finding T5):
the Phase-2 lattice always contains the trivial vector  n*S_D*e_m  (norm n*S_D),
while ||v_planted|| ~ n*sqrt(2m/3 + 4/3) > n*S_D.  So the planted vector is NEVER
lambda_1, on any curve, for any m >= 1, and no choice of S_D changes this (both
vectors scale linearly in S_D).  Recovery is therefore a BDD/coset condition, not
an SVP condition.

This script implements the fix proposed as Thread 23: quotient out the trivial
direction.  Since ker(phi) = <n*S_D*e_m> exactly (n prime, B_i != 0), deleting
column m is a lattice homomorphism

    phi : L -> Z^{2m+1},   ker(phi) = <n*S_D*e_m>

whose image is a rank-(2m+1) lattice with det(phi(L)) = det(L)/(n*S_D).
The map is LOSSLESS for the secret: two preimages of phi(v_planted) differ by a
multiple of n*S_D*e_m, i.e. by d -> d + n, so d mod n is still well defined.
It is recovered algebraically from the surviving (k1_i, k2_i) coordinates:

    k_full_0 = k1_0 + lam*k2_0 (mod n),   d = (k_full_0 - A_0) * B_0^{-1} (mod n).

Experiments
  E1  sanity: does the projected attack recover d at all?
  E2  is the planted vector lambda_1 in the projected lattice? (sv/pv)
  E3  the 2026-07-29 T4 K1-grid, original vs projected, 5 seeds, both 12-bit curves
  E4  Gaussian-heuristic prediction of the K1 wall vs observed

Run: python3 glv_hnp_phase2_projected.py
"""

import math
import random
from fpylll import IntegerMatrix, LLL, BKZ

# ---------------------------------------------------------------------------
# EC + number theory helpers (verbatim from glv_hnp_phase2_lambda_threshold.py)
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
        i, t2 = 0, t
        while t2 != 1:
            t2 = t2 * t2 % p; i += 1
        b = pow(c, 1 << (m - i - 1), p)
        m, c = i, b * b % p
        t, r = t * c % p, r * b % p

def find_generator(p, b, n, seed=12345):
    rng = random.Random(seed)
    for _ in range(2000):
        x = rng.randint(0, p - 1)
        rhs = (pow(x, 3, p) + b) % p
        y = tonelli_shanks(rhs, p)
        if y is None or y == 0: continue
        P = (x, y)
        if ec_mul(P, n, p) is None and P is not None:
            return P
    return None

def glv_roots(n):
    """Both roots of x^2+x+1 = 0 mod n, or None."""
    sq = tonelli_shanks((n - 3) % n, n)
    if sq is None: return None
    inv2 = modinv(2, n)
    r1 = (n - 1 + sq) * inv2 % n
    r2 = (n - 1 + (n - sq)) * inv2 % n
    if (r1 * r1 + r1 + 1) % n != 0: r1, r2 = r2, r1
    if (r1 * r1 + r1 + 1) % n != 0: return None
    return (min(r1, r2), max(r1, r2))

def lam_star(lam, n):
    return min(lam, n - lam) / n

# ---------------------------------------------------------------------------
# Phase-2 lattice (ORIGINAL construction, verbatim from glv_hnp_phase2_20bit.py)
# ---------------------------------------------------------------------------

def scales(n, k1_bound, k2_bound):
    return (n // k1_bound, 1, max(1, n // k2_bound), n)

def gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed=42):
    rng = random.Random(seed)
    sigs = []
    attempts = 0
    while len(sigs) < m and attempts < 200000:
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
    """ORIGINAL: dim 2m+2, columns [k1_res | d | k2 | kannan]."""
    m = len(sigs)
    dim = 2 * m + 2
    S_K1, S_D, S_K2, S_KANNAN = scales(n, k1_bound, k2_bound)
    M = [[0] * dim for _ in range(dim)]
    for i in range(m):
        M[i][i] = n * S_K1
    for i in range(m):
        M[m][i] = sigs[i]['B'] * S_K1
    M[m][m] = S_D
    for i in range(m):
        M[m + 1 + i][i] = -lam * S_K1
        M[m + 1 + i][m + 1 + i] = S_K2
    for i in range(m):
        M[2 * m + 1][i] = sigs[i]['A'] * S_K1
    M[2 * m + 1][dim - 1] = S_KANNAN
    return M

def planted_vector(sigs, d_secret, n, k1_bound, k2_bound):
    m = len(sigs)
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    v = [0] * (2 * m + 2)
    for i in range(m):
        v[i] = sigs[i]['k1'] * S_K1
        v[m + 1 + i] = sigs[i]['k2'] * S_K2
    v[m] = d_secret * S_D
    v[2 * m + 1] = S_KAN
    return v

def recover_d(M_reduced, m, n, S_KANNAN, d_secret):
    dim = 2 * m + 2
    for row in M_reduced:
        last = row[dim - 1]
        if abs(last) != S_KANNAN: continue
        sign = 1 if last > 0 else -1
        d_cand = (sign * row[m]) % n
        if d_cand == 0: continue
        if d_cand == d_secret:
            return d_cand
    return None

# ---------------------------------------------------------------------------
# NEW (Thread 23): projected lattice = original with column m deleted
# ---------------------------------------------------------------------------

def build_projected_lattice(sigs, n, lam, k1_bound, k2_bound):
    """
    phi(L) subset Z^{2m+1}, columns [k1_res (m) | k2 (m) | kannan (1)].
    Rows: the same 2m+2 generators of L with coordinate m dropped.  The result
    has rank 2m+1 (one generator becomes redundant); fplll's LLL handles the
    rank deficiency by emitting a leading zero row.
    """
    m = len(sigs)
    S_K1, S_D, S_K2, S_KANNAN = scales(n, k1_bound, k2_bound)
    cols = 2 * m + 1
    rows = []
    for i in range(m):                                   # n*S_K1*e_i
        r = [0] * cols; r[i] = n * S_K1; rows.append(r)
    r = [0] * cols                                       # d-direction generator
    for i in range(m):
        r[i] = sigs[i]['B'] * S_K1
    rows.append(r)
    for i in range(m):                                   # k2 generators
        r = [0] * cols
        r[i] = -lam * S_K1
        r[m + i] = S_K2
        rows.append(r)
    r = [0] * cols                                       # Kannan / target row
    for i in range(m):
        r[i] = sigs[i]['A'] * S_K1
    r[2 * m] = S_KANNAN
    rows.append(r)
    return rows

def planted_vector_projected(sigs, n, k1_bound, k2_bound):
    m = len(sigs)
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    v = [0] * (2 * m + 1)
    for i in range(m):
        v[i] = sigs[i]['k1'] * S_K1
        v[m + i] = sigs[i]['k2'] * S_K2
    v[2 * m] = S_KAN
    return v

def recover_d_projected(reduced, sigs, n, lam, k1_bound, k2_bound):
    """
    Scan reduced rows for the coset marker |last| == S_KANNAN, read off
    (k1_i, k2_i), then solve for d algebraically and verify on ALL signatures.
    """
    m = len(sigs)
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    for row in reduced:
        last = row[2 * m]
        if abs(last) != S_KAN: continue
        sign = 1 if last > 0 else -1
        k1s, k2s, ok = [], [], True
        for i in range(m):
            a, b = sign * row[i], sign * row[m + i]
            if a % S_K1 or b % S_K2:      # not an integral (k1,k2) lift
                ok = False; break
            k1s.append(a // S_K1); k2s.append(b // S_K2)
        if not ok: continue
        # d from signature 0, then verified against every signature
        k_full0 = (k1s[0] + lam * k2s[0]) % n
        try:
            d_cand = (k_full0 - sigs[0]['A']) * modinv(sigs[0]['B'], n) % n
        except ValueError:
            continue
        if d_cand == 0: continue
        if all((sigs[i]['A'] + sigs[i]['B'] * d_cand) % n
               == (k1s[i] + lam * k2s[i]) % n for i in range(m)):
            return d_cand
    return None

# ---------------------------------------------------------------------------
# Runners
# ---------------------------------------------------------------------------

def norm(v):
    return math.sqrt(sum(float(x) * float(x) for x in v))

def run_original(curve, m, d_secret, k1_bound, seed=42, use_bkz=False, beta=20):
    p, b, n, lam, G = curve
    k2 = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2, seed)
    if len(sigs) < m: return None
    M = build_glv_lattice(sigs, n, lam, k1_bound, k2)
    dim = 2 * m + 2
    A = IntegerMatrix.from_matrix(M)
    BKZ.reduction(A, BKZ.Param(beta)) if use_bkz else LLL.reduction(A)
    red = [[A[i][j] for j in range(dim)] for i in range(dim)]
    _, _, _, S_KAN = scales(n, k1_bound, k2)
    ok = recover_d(red, m, n, S_KAN, d_secret) is not None
    pn = norm(planted_vector(sigs, d_secret, n, k1_bound, k2))
    sn = min((norm(r) for r in red if any(r)), default=float('nan'))
    return ok, pn, sn

def run_projected(curve, m, d_secret, k1_bound, seed=42, use_bkz=False, beta=20):
    p, b, n, lam, G = curve
    k2 = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2, seed)
    if len(sigs) < m: return None
    rows = build_projected_lattice(sigs, n, lam, k1_bound, k2)
    cols = 2 * m + 1
    A = IntegerMatrix.from_matrix(rows)
    BKZ.reduction(A, BKZ.Param(beta)) if use_bkz else LLL.reduction(A)
    red = [[A[i][j] for j in range(cols)] for i in range(A.nrows)]
    d_cand = recover_d_projected(red, sigs, n, lam, k1_bound, k2)
    ok = (d_cand == d_secret)
    pn = norm(planted_vector_projected(sigs, n, k1_bound, k2))
    sn = min((norm(r) for r in red if any(r)), default=float('nan'))
    return ok, pn, sn

def rate(runner, curve, m, k1_bound, seeds, **kw):
    wins, ratios = 0, []
    for s in seeds:
        d = random.Random(s + 7777).randint(1, curve[2] - 1)
        res = runner(curve, m, d, k1_bound, s, **kw)
        if res is None: continue
        ok, pn, sn = res
        wins += bool(ok); ratios.append(sn / pn if pn else float('nan'))
    return wins, len(seeds), (sum(ratios) / len(ratios) if ratios else float('nan'))

# ---------------------------------------------------------------------------
# Gaussian-heuristic predictor
# ---------------------------------------------------------------------------

def gh_tau(m, n, k1_bound, k2_bound, projected):
    """
    tau = E||v_planted|| / lambda_1^GH.   Predict recovery when tau < 1.

    det(L)      = (n*S_K1)^m * S_D * S_K2^m * S_KANNAN         [block triangular]
    det(phi(L)) = det(L) / (n*S_D)
    """
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    log_det = (m * math.log(n * S_K1) + math.log(S_D)
               + m * math.log(S_K2) + math.log(S_KAN))
    if projected:
        log_det -= math.log(n * S_D)
        dim = 2 * m + 1
        pn2 = (m * (k1_bound * S_K1) ** 2 / 3.0
               + m * (k2_bound * S_K2) ** 2 / 3.0 + S_KAN ** 2)
    else:
        dim = 2 * m + 2
        pn2 = (m * (k1_bound * S_K1) ** 2 / 3.0 + (n * S_D) ** 2 / 3.0
               + m * (k2_bound * S_K2) ** 2 / 3.0 + S_KAN ** 2)
    lam1 = math.sqrt(dim / (2 * math.pi * math.e)) * math.exp(log_det / dim)
    return math.sqrt(pn2) / lam1


# ===========================================================================

def _main():
    print("=" * 78)
    print("Thread 23 — projected (trivial-vector-free) GLV-HNP Phase-2 lattice")
    print("=" * 78)

    SEEDS = [42, 1234, 9999, 555, 31337]

    # The two 12-bit curves from the 2026-07-29 T4 grid.
    HIST = [
        # label,          p,    b, n,    lam,  m
        ("12-bit/2557",   2557, 2, 2659, 1755, 8),
        ("12-bit/2677",   2677, 2, 2647, 185,  10),
    ]
    curves = []
    for label, p, b, n, lam, m in HIST:
        G = find_generator(p, b, n)
        assert G is not None, f"no generator for {label}"
        assert (lam * lam + lam + 1) % n == 0, f"bad lam for {label}"
        curves.append((label, (p, b, n, lam, G), m))
        print(f"  {label}: n={n} ({n.bit_length()}b) lam={lam} "
              f"lam*={lam_star(lam, n):.4f} m={m}")

    # -----------------------------------------------------------------------
    print("\n" + "=" * 78)
    print("E1/E2 — sanity + is the planted vector lambda_1 now?  (K1=4, LLL)")
    print("=" * 78)
    print(f"{'curve':<14} {'lattice':<10} {'wins':<7} {'sv/pv':<9} {'tau_GH':<8}")
    for label, cur, m in curves:
        n = cur[2]; k2b = math.isqrt(n) + 1
        for name, runner, proj in (("original", run_original, False),
                                   ("projected", run_projected, True)):
            w, t, r = rate(runner, cur, m, 4, SEEDS)
            tau = gh_tau(m, n, 4, k2b, proj)
            print(f"{label:<14} {name:<10} {str(w)+'/'+str(t):<7} {r:<9.3f} {tau:<8.3f}")

    # -----------------------------------------------------------------------
    print("\n" + "=" * 78)
    print("E3 — the 2026-07-29 T4 K1 grid: original vs projected (5 seeds, LLL)")
    print("=" * 78)
    K1_GRID = [2, 3, 4, 6, 8, 12, 16, 24]
    e3 = {}
    for label, cur, m in curves:
        for name, runner in (("original", run_original), ("projected", run_projected)):
            row = []
            for k1 in K1_GRID:
                w, t, _ = rate(runner, cur, m, k1, SEEDS)
                row.append(w)
            e3[(label, name)] = row
            print(f"{label:<14} {name:<10} " +
                  " ".join(f"K1={k}:{w}/5" for k, w in zip(K1_GRID, row)))

    # -----------------------------------------------------------------------
    print("\n" + "=" * 78)
    print("E4 — Gaussian-heuristic tau vs observed success  (tau<1 predicts success)")
    print("=" * 78)
    for label, cur, m in curves:
        n = cur[2]; k2b = math.isqrt(n) + 1
        for name, proj in (("original", False), ("projected", True)):
            taus = [gh_tau(m, n, k1, k2b, proj) for k1 in K1_GRID]
            obs = e3[(label, name)]
            print(f"\n{label} [{name}]  m={m}")
            print("  " + " ".join(f"{'K1='+str(k):>9}" for k in K1_GRID))
            print("  " + " ".join(f"{tau:>9.2f}" for tau in taus))
            print("  " + " ".join(f"{str(w)+'/5':>9}" for w in obs))
            # crossing point
            cross = next((k for k, tau in zip(K1_GRID, taus) if tau >= 1.0), None)
            wall = next((k for k, w in zip(K1_GRID, obs) if w < len(SEEDS)), None)
            print(f"  predicted wall (first tau>=1): K1={cross}   "
                  f"observed wall (first <5/5): K1={wall}")

    print("\n" + "=" * 78)
    print("done")
    print("=" * 78)


if __name__ == "__main__":
    _main()
