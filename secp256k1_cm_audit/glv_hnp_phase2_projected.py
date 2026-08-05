"""
Thread 23 — reformulate the GLV-HNP Phase-2 lattice so the planted vector is
the shortest vector.

Background (RESEARCH_AUTOLAB_LOG.md, 2026-07-29, exp T5):
  The Phase-2 Kannan lattice always contains the trivial vector n*S_D*e_m
  (norm n), while the planted vector has norm ~ n*sqrt(2m/3 + 4/3).  So the
  planted vector is NEVER lambda_1 for m >= 1, and recovery is a BDD/coset
  problem, not an SVP problem.  No choice of S_D fixes this: both vectors
  scale linearly in S_D.

This script tests the proposed fix.  The trivial vector lives entirely in the
d-column (100% of its energy, |sv[m]|/n = 1.0000 exactly), and d is only
defined mod n, so the d-column carries no metric information.  Deleting it is
the orthogonal projection pi along e_m; since ker(pi) ∩ L = n*S_D*e_m*Z has
rank 1, pi(L) is a lattice of rank 2m+1 and the trivial vector maps to 0.

Three methods are compared on identical signature sets:

  BASE  the 2026-06-15 Kannan lattice, dim 2m+2, d read off column m.
  PROJ  the same lattice with column m deleted, dim 2m+1.  d is no longer
        read off directly; it is reconstructed from the recovered (k1_i,k2_i)
        via  d = (k1_i + lam*k2_i - A_i) * B_i^{-1}  mod n.
  CVP   no Kannan embedding at all: Babai nearest-plane (exact rational GS)
        on the LLL-reduced 2m-dimensional lattice, target (-A_i*S_K1, 0).

All three use a SELF-VERIFYING success test (no d_secret oracle): a candidate
d is accepted only if every signature satisfies A_i + B_i*d = k1_i + lam*k2_i
(mod n) with the recovered k1_i in [0,K1) and k2_i in [0,K2).  The test is
cross-checked against d_secret so it cannot be vacuous.

Falsifier (stated 2026-07-29):
  "if sv/pv rises above 1 after the reformulation AND the K1 wall in T4 moves
   outward on the lam*=0.07 curve (currently K1 ~ 4-6), the reformulation is a
   real improvement; if the wall stays at K1 ~ 4-6, then the wall is
   information-theoretic and Phase 2 is at its ceiling."

Run: python3 glv_hnp_phase2_projected.py
Deps: fpylll, sympy   (pip install fpylll cysignals sympy)

!! SUPERSEDED CRITERION -- read before reusing these numbers.
   The success test below (STRICT) additionally demands that the recovered
   k1/k2 blocks be the PLANTED nonce split.  That understates the attack: the
   real attacker reads d off column m and checks d*G == Q against the public
   key, which also succeeds on alternative GLV splits of the same nonce.  Under
   the real criterion BASE reproduces the 2026-07-29 T4 row exactly and PROJ
   ties it -- see glv_hnp_phase2_cmin.py (Thread 23c).  The E2/E3/E4 tables in
   glv_hnp_phase2_projected_output.txt are STRICT numbers; do not compare them
   against log entries from other runs without this adjustment.
"""

import math
import random
from fractions import Fraction

import sympy
from fpylll import IntegerMatrix, LLL, BKZ

# ---------------------------------------------------------------------------
# Minimal EC arithmetic (verbatim from glv_hnp_phase2_lambda_threshold.py)
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

def find_generator(p, b, n, seed=12345):
    rng = random.Random(seed)
    for _ in range(10000):
        x = rng.randint(0, p - 1)
        rhs = (pow(x, 3, p) + b) % p
        y = tonelli_shanks(rhs, p)
        if y is not None and y != 0:
            P = (x, y)
            if ec_mul(P, n, p) is None:
                return P
    return None

# ---------------------------------------------------------------------------
# j=0 CM helpers
# ---------------------------------------------------------------------------

def eisenstein_decompose(p):
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
    return [2*a - b, -2*a + b, -(a + b), a + b, 2*b - a, a - 2*b]

def glv_roots(n):
    sq = tonelli_shanks((n - 3) % n, n)
    if sq is None:
        return None
    inv2 = modinv(2, n)
    r1 = (n - 1 + sq) * inv2 % n
    r2 = (n - 1 + (n - sq)) * inv2 % n
    if (r1 * r1 + r1 + 1) % n != 0:
        r1, r2 = r2, r1
    if (r1 * r1 + r1 + 1) % n != 0:
        return None
    return (min(r1, r2), max(r1, r2))

def lam_star(lam, n):
    return min(lam, n - lam) / n

def build_curve(p, n, seed=12345):
    rng = random.Random(seed)
    for _ in range(400):
        b = rng.randint(1, p - 1)
        x = rng.randint(0, p - 1)
        rhs = (pow(x, 3, p) + b) % p
        y = tonelli_shanks(rhs, p)
        if y is None or y == 0:
            continue
        if ec_mul((x, y), n, p) is None:
            G = find_generator(p, b, n, seed=seed)
            if G is not None:
                return (p, b, n, None, G)
    return None

def search_curves(lo, hi, per_bin=2, nbins=10):
    """j=0 GLV curves with p in [lo,hi), n prime, n = 1 mod 3, bucketed by lam*."""
    bins = {i: [] for i in range(nbins)}
    p = int(sympy.nextprime(lo))
    while p < hi:
        if p % 3 == 1:
            eis = eisenstein_decompose(p)
            if eis is not None:
                a_e, b_e = eis
                for t in j0_traces(a_e, b_e):
                    n_cand = p + 1 - t
                    if n_cand < 2 or n_cand % 3 != 1 or not sympy.isprime(n_cand):
                        continue
                    roots = glv_roots(n_cand)
                    if roots is None:
                        continue
                    lam = roots[0]
                    idx = min(nbins - 1, int(lam_star(lam, n_cand) / (0.5 / nbins)))
                    if len(bins[idx]) >= per_bin:
                        continue
                    cur = build_curve(p, n_cand)
                    if cur is None:
                        continue
                    bins[idx].append((p, cur[1], n_cand, lam, cur[4]))
        if all(len(v) >= per_bin for v in bins.values()):
            break
        p = int(sympy.nextprime(p))
    out = []
    for i in range(nbins):
        out.extend(bins[i])
    return out

# ---------------------------------------------------------------------------
# Signatures + scales (verbatim construction)
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

def norm(v):
    return math.sqrt(sum(float(x) * float(x) for x in v))

# ---------------------------------------------------------------------------
# BASE — the 2026-06-15 Kannan lattice (dim 2m+2)
# ---------------------------------------------------------------------------

def build_base_lattice(sigs, n, lam, k1_bound, k2_bound):
    m = len(sigs)
    dim = 2 * m + 2
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
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
    M[2 * m + 1][dim - 1] = S_KAN
    return M

def base_planted(sigs, d_secret, n, k1_bound, k2_bound):
    m = len(sigs)
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    v = [0] * (2 * m + 2)
    for i in range(m):
        v[i] = sigs[i]['k1'] * S_K1
        v[m + 1 + i] = sigs[i]['k2'] * S_K2
    v[m] = d_secret * S_D
    v[2 * m + 1] = S_KAN
    return v

# ---------------------------------------------------------------------------
# PROJ — the same lattice with the d-column (column m) deleted (dim 2m+1)
#        columns: 0..m-1 = k1 block, m..2m-1 = k2 block, 2m = Kannan
# ---------------------------------------------------------------------------

def build_proj_lattice(sigs, n, lam, k1_bound, k2_bound):
    m = len(sigs)
    ncols = 2 * m + 1
    S_K1, _S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    rows = []
    for i in range(m):                       # modular rows
        r = [0] * ncols
        r[i] = n * S_K1
        rows.append(r)
    r = [0] * ncols                          # d-row: S_D entry is gone
    for i in range(m):
        r[i] = sigs[i]['B'] * S_K1
    rows.append(r)
    for i in range(m):                       # k2 rows
        r = [0] * ncols
        r[i] = -lam * S_K1
        r[m + i] = S_K2
        rows.append(r)
    r = [0] * ncols                          # Kannan row
    for i in range(m):
        r[i] = sigs[i]['A'] * S_K1
    r[2 * m] = S_KAN
    rows.append(r)
    return rows

def proj_planted(sigs, n, k1_bound, k2_bound):
    m = len(sigs)
    S_K1, _S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    v = [0] * (2 * m + 1)
    for i in range(m):
        v[i] = sigs[i]['k1'] * S_K1
        v[m + i] = sigs[i]['k2'] * S_K2
    v[2 * m] = S_KAN
    return v

# ---------------------------------------------------------------------------
# Self-verifying recovery: given candidate (k1_i, k2_i) for all i, rebuild d
# and check it against EVERY signature.  No d_secret is consulted.
# ---------------------------------------------------------------------------

def d_from_k(sigs, k1s, k2s, n, lam, k1_bound, k2_bound):
    m = len(sigs)
    for i in range(m):
        if not (0 <= k1s[i] < k1_bound and 0 <= k2s[i] < k2_bound):
            return None
    for i in range(m):
        B = sigs[i]['B'] % n
        if B == 0 or math.gcd(B, n) != 1:
            continue
        k_full = (k1s[i] + lam * k2s[i]) % n
        d_cand = (k_full - sigs[i]['A']) * modinv(B, n) % n
        if d_cand == 0:
            continue
        ok = True
        for j in range(m):
            if (sigs[j]['A'] + sigs[j]['B'] * d_cand - k1s[j] - lam * k2s[j]) % n != 0:
                ok = False
                break
        if ok:
            return d_cand
    return None

def recover_base(reduced, sigs, n, lam, k1_bound, k2_bound):
    m = len(sigs)
    dim = 2 * m + 2
    S_K1, _S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    for row in reduced:
        if abs(row[dim - 1]) != S_KAN:
            continue
        r = row if row[dim - 1] > 0 else [-x for x in row]
        if any(r[i] % S_K1 for i in range(m)):
            continue
        if any(r[m + 1 + i] % S_K2 for i in range(m)):
            continue
        k1s = [r[i] // S_K1 for i in range(m)]
        k2s = [r[m + 1 + i] // S_K2 for i in range(m)]
        d = d_from_k(sigs, k1s, k2s, n, lam, k1_bound, k2_bound)
        if d is not None:
            return d
    return None

def recover_proj(reduced, sigs, n, lam, k1_bound, k2_bound):
    m = len(sigs)
    S_K1, _S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    for row in reduced:
        if abs(row[2 * m]) != S_KAN:
            continue
        r = row if row[2 * m] > 0 else [-x for x in row]
        if any(r[i] % S_K1 for i in range(m)):
            continue
        if any(r[m + i] % S_K2 for i in range(m)):
            continue
        k1s = [r[i] // S_K1 for i in range(m)]
        k2s = [r[m + i] // S_K2 for i in range(m)]
        d = d_from_k(sigs, k1s, k2s, n, lam, k1_bound, k2_bound)
        if d is not None:
            return d
    return None

# ---------------------------------------------------------------------------
# CVP — Babai nearest plane with exact rational Gram-Schmidt
#       lattice (dim 2m): rows {n*S_K1*e_i}, {B_i*S_K1}, {-lam*S_K1*e_i + S_K2*e_{m+i}}
#       target t = (-A_i*S_K1, 0, ..., 0);  u - t = (k1_i*S_K1, k2_i*S_K2)
# ---------------------------------------------------------------------------

def build_cvp_lattice(sigs, n, lam, k1_bound, k2_bound):
    m = len(sigs)
    ncols = 2 * m
    S_K1, _S_D, S_K2, _S_KAN = scales(n, k1_bound, k2_bound)
    rows = []
    for i in range(m):
        r = [0] * ncols
        r[i] = n * S_K1
        rows.append(r)
    r = [0] * ncols
    for i in range(m):
        r[i] = sigs[i]['B'] * S_K1
    rows.append(r)
    for i in range(m):
        r = [0] * ncols
        r[i] = -lam * S_K1
        r[m + i] = S_K2
        rows.append(r)
    return rows

def gram_schmidt_exact(basis):
    """Exact rational GS.  Returns (bstar, mu) with bstar[i] orthogonal."""
    bstar, mu = [], []
    for v in basis:
        w = [Fraction(x) for x in v]
        coeffs = []
        for j, bs in enumerate(bstar):
            nrm = sum(c * c for c in bs)
            if nrm == 0:
                coeffs.append(Fraction(0))
                continue
            c = sum(Fraction(v[k]) * bs[k] for k in range(len(v))) / nrm
            coeffs.append(c)
            for k in range(len(w)):
                w[k] -= c * bs[k]
        bstar.append(w)
        mu.append(coeffs)
    return bstar, mu

def babai_nearest_plane(basis, target):
    """Exact Babai nearest plane on `basis` (assumed LLL-reduced, no zero rows)."""
    bstar, _ = gram_schmidt_exact(basis)
    b = [Fraction(x) for x in target]
    coeffs = [0] * len(basis)
    for i in range(len(basis) - 1, -1, -1):
        nrm = sum(c * c for c in bstar[i])
        if nrm == 0:
            continue
        c = sum(b[k] * bstar[i][k] for k in range(len(b))) / nrm
        # round-half-to-even on an exact Fraction
        ci = int(c + Fraction(1, 2)) if c >= 0 else -int(-c + Fraction(1, 2))
        coeffs[i] = ci
        for k in range(len(b)):
            b[k] -= ci * Fraction(basis[i][k])
    out = [0] * len(target)
    for i, ci in enumerate(coeffs):
        if ci:
            for k in range(len(out)):
                out[k] += ci * basis[i][k]
    return out

def strip_zero_rows(rows):
    return [r for r in rows if any(r)]

# ---------------------------------------------------------------------------
# One experiment on one (curve, m, K1, seed)
# ---------------------------------------------------------------------------

def run_all_methods(curve, m, d_secret, k1_bound, seed=42, methods=('base', 'proj', 'cvp'),
                    use_bkz=False, bkz_beta=20):
    """Returns dict: method -> (ok, sv_over_pv or None)."""
    p, b, n, lam, G = curve
    d_secret %= n
    k2_bound = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return None
    out = {}

    if 'base' in methods:
        M = build_base_lattice(sigs, n, lam, k1_bound, k2_bound)
        A = IntegerMatrix.from_matrix(M)
        if use_bkz:
            BKZ.reduction(A, BKZ.Param(bkz_beta))
        else:
            LLL.reduction(A)
        red = [[A[i][j] for j in range(2 * m + 2)] for i in range(A.nrows)]
        d = recover_base(red, sigs, n, lam, k1_bound, k2_bound)
        pv = norm(base_planted(sigs, d_secret, n, k1_bound, k2_bound))
        sv = min((norm(r) for r in red if any(r)), default=float('nan'))
        out['base'] = (d == d_secret, sv / pv)
        assert d is None or d == d_secret, "self-verifying test accepted a wrong d (BASE)"

    if 'proj' in methods:
        rows = build_proj_lattice(sigs, n, lam, k1_bound, k2_bound)
        A = IntegerMatrix.from_matrix(rows)
        if use_bkz:
            BKZ.reduction(A, BKZ.Param(bkz_beta))
        else:
            LLL.reduction(A)
        red = [[A[i][j] for j in range(2 * m + 1)] for i in range(A.nrows)]
        d = recover_proj(red, sigs, n, lam, k1_bound, k2_bound)
        pv = norm(proj_planted(sigs, n, k1_bound, k2_bound))
        sv = min((norm(r) for r in red if any(r)), default=float('nan'))
        out['proj'] = (d == d_secret, sv / pv)
        assert d is None or d == d_secret, "self-verifying test accepted a wrong d (PROJ)"

    if 'cvp' in methods:
        S_K1, _S_D, S_K2, _K = scales(n, k1_bound, k2_bound)
        rows = build_cvp_lattice(sigs, n, lam, k1_bound, k2_bound)
        A = IntegerMatrix.from_matrix(rows)
        LLL.reduction(A)
        red = strip_zero_rows([[A[i][j] for j in range(2 * m)] for i in range(A.nrows)])
        t = [0] * (2 * m)
        for i in range(m):
            t[i] = -sigs[i]['A'] * S_K1
        u = babai_nearest_plane(red, t)
        diff = [u[k] - t[k] for k in range(2 * m)]
        d = None
        if not any(diff[i] % S_K1 for i in range(m)) and \
           not any(diff[m + i] % S_K2 for i in range(m)):
            k1s = [diff[i] // S_K1 for i in range(m)]
            k2s = [diff[m + i] // S_K2 for i in range(m)]
            d = d_from_k(sigs, k1s, k2s, n, lam, k1_bound, k2_bound)
        out['cvp'] = (d == d_secret, None)
        assert d is None or d == d_secret, "self-verifying test accepted a wrong d (CVP)"

    return out

def success_rates(curve, m, k1_bound, seeds, methods=('base', 'proj', 'cvp'),
                  use_bkz=False, bkz_beta=20):
    wins = {k: 0 for k in methods}
    ratios = {k: [] for k in methods}
    tot = 0
    for seed in seeds:
        n = curve[2]
        d_trial = random.Random(seed + 7777).randint(1, n - 1)
        res = run_all_methods(curve, m, d_trial, k1_bound, seed, methods,
                              use_bkz=use_bkz, bkz_beta=bkz_beta)
        if res is None:
            continue
        tot += 1
        for k in methods:
            ok, rr = res[k]
            wins[k] += bool(ok)
            if rr is not None:
                ratios[k].append(rr)
    mean = {k: (sum(v) / len(v) if v else float('nan')) for k, v in ratios.items()}
    return wins, tot, mean


# ===========================================================================
print("=" * 78)
print("Thread 23 — projecting out the trivial n*e_m direction (GLV-HNP Phase 2)")
print("=" * 78)

SEEDS = [42, 1234, 9999, 555, 31337]

HIST = [
    # label,              p,    b, n,    lam,  K1, m
    ("8-bit/199",         211,  2, 199,  106,  2,  6),
    ("12-bit/2557",       2557, 2, 2659, 1755, 8,  8),
    ("12-bit/2677 FAIL",  2677, 2, 2647, 185,  8,  10),
]
hist = []
for label, p, b, n, lam, k1, m in HIST:
    G = find_generator(p, b, n)
    assert G is not None, f"no generator for {label}"
    assert (lam * lam + lam + 1) % n == 0, f"bad lam for {label}"
    hist.append((label, (p, b, n, lam, G), k1, m))

# ===========================================================================
# E1 — does the projection make the planted vector the shortest vector?
# ===========================================================================
print("\n" + "-" * 78)
print("E1: shortest/planted norm ratio, BASE vs PROJ  (m=12, K1 as historical)")
print("-" * 78)
print("T5 (2026-07-29) measured sv/pv in [0.337,0.368] for BASE — the trivial")
print("vector n*e_m.  If the projection works, PROJ's ratio must rise to 1.0")
print("whenever the planted vector is found, and >1 is impossible.\n")
print(f"{'curve':<18} {'lam*':>7} {'K1':>3}  {'BASE sv/pv':>11} {'PROJ sv/pv':>11}")
for label, curve, k1, _m in hist:
    p, b, n, lam, G = curve
    w, t, mr = success_rates(curve, 12, k1, SEEDS, methods=('base', 'proj'))
    print(f"{label:<18} {lam_star(lam,n):>7.4f} {k1:>3}  "
          f"{mr['base']:>11.4f} {mr['proj']:>11.4f}")

# ===========================================================================
# E2 — the T4 grid, all three methods.  Does the K1 wall move outward?
# ===========================================================================
print("\n" + "-" * 78)
print("E2: T4 K1-sweep replication (m=12, 5 seeds) — BASE / PROJ / CVP")
print("-" * 78)
print("2026-07-29 T4 BASE numbers to reproduce:")
print("  2557 (lam*=0.340): 5,5,5,5,5,4,1,0  over K1=2,3,4,6,8,12,16,24")
print("  2677 (lam*=0.070): 5,5,5,2,0,0,0,0")
print("Falsifier: the wall on the lam*=0.07 curve must move outward.\n")

K1_GRID = [2, 3, 4, 6, 8, 12, 16, 24]
for label, curve, _k1, _m in hist:
    p, b, n, lam, G = curve
    if n < 2000:
        continue
    print(f"{label}  (lam*={lam_star(lam,n):.4f}, n={n}, K2={math.isqrt(n)+1})")
    lines = {k: [] for k in ('base', 'proj', 'cvp')}
    for k1 in K1_GRID:
        w, t, _ = success_rates(curve, 12, k1, SEEDS)
        for k in lines:
            lines[k].append(f"{w[k]}/{t}")
    hdr = "  " + " ".join(f"K1={k:<4}" for k in K1_GRID)
    print(hdr)
    for k in ('base', 'proj', 'cvp'):
        print(f"  {k.upper():<5}" + " ".join(f"{c:<6}" for c in lines[k]))
    print()

# ===========================================================================
# E3 — at the K1 where BASE dies, does more data rescue PROJ/CVP?
# ===========================================================================
print("-" * 78)
print("E3: more data at the BASE wall (2677 curve, K1=8; BASE was 0/5 for all m)")
print("-" * 78)
fail_curve = [c for _l, c, _k, _m in hist if c[2] == 2647][0]
print(f"{'m':>4}  " + "  ".join(f"{k.upper():>6}" for k in ('base', 'proj', 'cvp')))
for m_try in (8, 12, 16, 24, 32):
    w, t, _ = success_rates(fail_curve, m_try, 8, SEEDS)
    print(f"{m_try:>4}  " + "  ".join(f"{w[k]}/{t:<4}" for k in ('base', 'proj', 'cvp')))

# ===========================================================================
# E4 — 17-bit cross-curve check at the eff levels from T3
# ===========================================================================
print("\n" + "-" * 78)
print("E4: 17-bit curve sweep at eff = K1*K2/n in {0.05, 0.15, 0.25}")
print("-" * 78)
print("T3 (2026-07-29) BASE result: eff=0.05 -> 19/20 curves, eff=0.15 -> 3/20,")
print("eff=0.25 -> 0/20.  Does the reformulation lift the 0.15 / 0.25 rows?\n")
curves17 = search_curves(2**16, 2**17, per_bin=2, nbins=10)
print(f"collected {len(curves17)} 17-bit curves "
      f"(lam* range {min(lam_star(c[3],c[2]) for c in curves17):.4f}"
      f"..{max(lam_star(c[3],c[2]) for c in curves17):.4f})\n")

for eff in (0.05, 0.15, 0.25):
    tally = {k: 0 for k in ('base', 'proj', 'cvp')}
    full = {k: 0 for k in ('base', 'proj', 'cvp')}
    for c in curves17:
        p, b, n, lam, G = c
        k2 = math.isqrt(n) + 1
        k1 = max(2, int(round(eff * n / k2)))
        w, t, _ = success_rates(c, 12, k1, SEEDS)
        for k in tally:
            tally[k] += w[k]
            full[k] += (w[k] == t)
    ncur = len(curves17)
    print(f"eff={eff:.2f}  " + "   ".join(
        f"{k.upper()}: {full[k]}/{ncur} curves 5/5, {tally[k]}/{ncur*len(SEEDS)} seeds"
        for k in ('base', 'proj', 'cvp')))

print("\n" + "=" * 78)
print("done")
print("=" * 78)
