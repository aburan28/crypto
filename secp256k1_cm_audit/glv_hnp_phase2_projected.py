"""
GLV-HNP Phase 2 — Thread 23: quotient out the trivial n*e_m direction.

Background (2026-07-29 log entry, T5)
-------------------------------------
In the Phase-2 Kannan lattice built by `build_glv_lattice` (verbatim from
`glv_hnp_phase2_20bit.py:262`) the shortest vector after LLL is, on every
instance tested, the *trivial* vector

    t = n * S_D * e_m          (||t|| = n,  |t[m]|/n = 1.0000 exactly)

which lies in L because n*row_m - sum_i B_i*row_i = n*S_D*e_m.  It carries no
information (d is only defined mod n) and no choice of S_D removes it: both t
and the planted vector scale linearly in S_D.  Measured sv/pv sat in
[0.337, 0.368] for successes and failures alike.

Thread 23 proposal (verbatim from the 2026-07-29 "Next step proposal"):
  "project the lattice along e_m (quotient out the trivial n*e_m direction)
   and solve BDD in the projection ... Falsifier: if sv/pv rises above 1 after
   the reformulation and the K1 wall in T4 moves outward on the lam*=0.07
   curve (currently K1~4-6), the reformulation is a real improvement; if the
   wall stays at K1~4-6, then the wall is information-theoretic and Phase 2 is
   at its ceiling."

The projection is exact: L cap R*e_m = Z * (n*S_D*e_m), so pi = "delete
coordinate m" is a lattice homomorphism onto a rank-(dim-1) lattice with
covol(pi L) = det(L) / (n*S_D).  The planted vector's k1- and k2-blocks are
untouched by pi, so d is reconstructed afterwards from any single signature:

    d = (k1_i + lam*k2_i - A_i) * B_i^{-1}  (mod n)

Experiments
-----------
E1  Analytic: rho_GH = ||pv|| / GH(L) for the original and projected lattices,
    and its m -> infinity limit.  Prediction to test against T4b.
E2  Correctness: projected lattice + oracle-free reconstruction recovers d
    exactly when the original does.
E3  The falsifier proper: sv/pv in the projected lattice.
E4  T4 K1-wall grid re-run, original vs projected, both anchor curves.
E5  m-sweep at K1=8 on the lam*=0.07 curve (T4b said more data does not help).

Run: python3 glv_hnp_phase2_projected.py
"""

import math
import random
from fpylll import IntegerMatrix, LLL, BKZ

# ---------------------------------------------------------------------------
# Minimal EC arithmetic + signature generation
# (verbatim from glv_hnp_phase2_20bit.py so the comparison to 2026-07-26 /
#  2026-07-29 is exact)
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

def find_generator(p, b, n):
    rng = random.Random(12345)
    for _ in range(10000):
        x = rng.randint(0, p - 1)
        rhs = (pow(x, 3, p) + b) % p
        y = tonelli_shanks(rhs, p)
        if y is not None and y != 0:
            P = (x, y)
            if ec_mul(P, n, p) is None:
                return P
    return None

def gen_signatures(G, d_secret, m, n, lam, p, b, k1_bound, k2_bound, seed=42):
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
        sigs.append({'A': A, 'B': B, 'k1': k1, 'k2': k2, 'k_full': k_full})
    return sigs

# ---------------------------------------------------------------------------
# Lattices
# ---------------------------------------------------------------------------

def build_glv_lattice(sigs, n, lam, k1_bound, k2_bound):
    """VERBATIM from glv_hnp_phase2_20bit.py:262. dim = 2m+2, lower triangular."""
    m = len(sigs)
    dim = 2 * m + 2
    S_K1 = n // k1_bound
    S_D = 1
    S_K2 = max(1, n // k2_bound)
    S_KANNAN = n

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

    return M, S_K1, S_D, S_K2, S_KANNAN


def project_out_d(M, m):
    """
    pi = delete coordinate m.  Exact quotient by L cap R*e_m = Z*(n*S_D*e_m).
    Returns a (2m+2) x (2m+1) generating set of rank 2m+1.

    Column layout of the image:
        0     .. m-1   k1-block   (was 0..m-1)
        m     .. 2m-1  k2-block   (was m+1..2m)
        2m             Kannan     (was 2m+1)
    """
    return [row[:m] + row[m + 1:] for row in M]


def planted_vector(sigs, n, lam, S_K1, S_D, S_K2, S_KANNAN, d_secret):
    """The exact planted lattice point, in the ORIGINAL coordinates."""
    m = len(sigs)
    dim = 2 * m + 2
    v = [0] * dim
    for i, sg in enumerate(sigs):
        v[i] = sg['k1'] * S_K1
        v[m + 1 + i] = sg['k2'] * S_K2
    v[m] = d_secret * S_D
    v[dim - 1] = S_KANNAN
    return v

def norm(v):
    return math.sqrt(sum(float(x) * float(x) for x in v))

# ---------------------------------------------------------------------------
# Recovery
# ---------------------------------------------------------------------------

def recover_d_original_oracle(M_reduced, m, n, S_KANNAN, d_secret):
    """VERBATIM from glv_hnp_phase2_20bit.py — uses d_secret as an oracle."""
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


def _consistent(d_cand, sigs, n, lam, k1_bound, k2_bound):
    """Oracle-free check: does d_cand explain every signature with in-range k1,k2?"""
    if d_cand is None or d_cand == 0:
        return False
    for sg in sigs:
        u = (sg['A'] + d_cand * sg['B']) % n
        # need u = k1 + lam*k2 mod n with 0<=k1<k1_bound, 0<=k2<k2_bound
        ok = False
        for k2 in range(k2_bound):
            if (u - lam * k2) % n < k1_bound:
                ok = True
                break
        if not ok:
            return False
    return True


def recover_d_original_free(M_reduced, sigs, n, lam, S_KANNAN, k1_bound, k2_bound):
    """Oracle-free variant of the original recovery (read d off coordinate m)."""
    m = len(sigs)
    dim = 2 * m + 2
    for row in M_reduced:
        last = row[dim - 1]
        if abs(last) != S_KANNAN: continue
        sign = 1 if last > 0 else -1
        d_cand = (sign * row[m]) % n
        if _consistent(d_cand, sigs, n, lam, k1_bound, k2_bound):
            return d_cand
    return None


def recover_d_projected(M_reduced, sigs, n, lam, S_K1, S_K2, S_KANNAN,
                        k1_bound, k2_bound):
    """
    Oracle-free reconstruction in the projected lattice.  d is gone from the
    coordinates, so it is rebuilt from the k1/k2 blocks of the candidate row:

        d = (k1_i + lam*k2_i - A_i) * B_i^{-1}  (mod n)
    """
    m = len(sigs)
    for row in M_reduced:
        last = row[2 * m]
        if abs(last) != S_KANNAN: continue
        sign = 1 if last > 0 else -1
        # k1 / k2 blocks must be exactly divisible by their scalings
        ok = True
        k1s, k2s = [], []
        for i in range(m):
            a, b = sign * row[i], sign * row[m + i]
            if a % S_K1 != 0 or b % S_K2 != 0:
                ok = False
                break
            k1s.append(a // S_K1)
            k2s.append(b // S_K2)
        if not ok:
            continue
        for i in range(m):
            if sigs[i]['B'] % n == 0:
                continue
            d_cand = (k1s[i] + lam * k2s[i] - sigs[i]['A']) * modinv(sigs[i]['B'], n) % n
            if _consistent(d_cand, sigs, n, lam, k1_bound, k2_bound):
                return d_cand
    return None

# ---------------------------------------------------------------------------
# Experiment drivers
# ---------------------------------------------------------------------------

def reduce_rows(M, nrows, ncols, use_bkz=False, beta=20):
    A = IntegerMatrix.from_matrix(M)
    if use_bkz:
        BKZ.reduction(A, BKZ.Param(beta))
    else:
        LLL.reduction(A)
    return [[A[i][j] for j in range(ncols)] for i in range(A.nrows)]


def run_pair(curve_params, m, d_secret, k1_bound, seed, use_bkz=False, beta=20):
    """
    One instance -> dict with original and projected outcomes + geometry.
    """
    p, b, n, lam, G = curve_params
    k2_bound = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, b, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return None
    M, S_K1, S_D, S_K2, S_KANNAN = build_glv_lattice(sigs, n, lam, k1_bound, k2_bound)
    dim = 2 * m + 2

    pv = planted_vector(sigs, n, lam, S_K1, S_D, S_K2, S_KANNAN, d_secret)
    pv_proj = pv[:m] + pv[m + 1:]

    red_o = reduce_rows(M, dim, dim, use_bkz, beta)
    red_p = reduce_rows(project_out_d(M, m), dim, dim - 1, use_bkz, beta)

    nz_o = [r for r in red_o if any(r)]
    nz_p = [r for r in red_p if any(r)]
    sv_o = min(norm(r) for r in nz_o)
    sv_p = min(norm(r) for r in nz_p)

    # is the projected shortest vector the planted one (up to sign)?
    shortest_p = min(nz_p, key=norm)
    is_planted = (shortest_p == pv_proj or [-x for x in shortest_p] == pv_proj)

    return {
        'ok_o_oracle': recover_d_original_oracle(red_o, m, n, S_KANNAN, d_secret) is not None,
        'ok_o_free': recover_d_original_free(red_o, sigs, n, lam, S_KANNAN,
                                             k1_bound, k2_bound) == d_secret,
        'ok_p': recover_d_projected(red_p, sigs, n, lam, S_K1, S_K2, S_KANNAN,
                                    k1_bound, k2_bound) == d_secret,
        'sv_o': sv_o, 'sv_p': sv_p,
        'pv_o': norm(pv), 'pv_p': norm(pv_proj),
        'shortest_is_planted': is_planted,
        'S': (S_K1, S_D, S_K2, S_KANNAN),
        'k2_bound': k2_bound,
    }


def gh(det_log, dim):
    """Gaussian heuristic given log(det)."""
    return math.sqrt(dim / (2 * math.pi * math.e)) * math.exp(det_log / dim)


def geometry(n, k1_bound, k2_bound, m):
    """Analytic rho_GH for the original and projected lattices."""
    S_K1 = n // k1_bound
    S_K2 = max(1, n // k2_bound)
    S_D, S_KANNAN = 1, n
    dim = 2 * m + 2
    ld_o = (m * math.log(n * S_K1) + math.log(S_D)
            + m * math.log(S_K2) + math.log(S_KANNAN))
    ld_p = ld_o - math.log(n * S_D)
    # E||pv||^2 with k1 ~ U[0,K1), k2 ~ U[0,K2), d ~ U[0,n)
    e_o = (m * (k1_bound ** 2 / 3.0) * S_K1 ** 2
           + (n ** 2 / 3.0) * S_D ** 2
           + m * (k2_bound ** 2 / 3.0) * S_K2 ** 2
           + S_KANNAN ** 2)
    e_p = e_o - (n ** 2 / 3.0) * S_D ** 2
    return {
        'eff': k1_bound * k2_bound / n,
        'rho_o': math.sqrt(e_o) / gh(ld_o, dim),
        'rho_p': math.sqrt(e_p) / gh(ld_p, dim - 1),
    }


# ---------------------------------------------------------------------------
# Anchor curves (same as 2026-06-15 / 2026-07-29)
# ---------------------------------------------------------------------------

print("=" * 74)
print("GLV-HNP Phase 2 — Thread 23: quotient out the trivial n*e_m direction")
print("=" * 74)

CURVES = {}
G1 = find_generator(2557, 2, 2659)
CURVES['12b/2557'] = (2557, 2, 2659, 1755, G1)      # lam* = 904/2659 = 0.340
G2 = find_generator(2677, 2, 2647)
CURVES['12b/2677'] = (2677, 2, 2647, 185, G2)       # lam* = 185/2647 = 0.070
G0 = find_generator(211, 2, 199)
CURVES['8b/199'] = (211, 2, 199, 106, G0)           # lam* = 93/199  = 0.467

for lbl, (p, b, n, lam, G) in CURVES.items():
    assert G is not None, f"no generator for {lbl}"
    lam_star = min(lam, n - lam) / n
    print(f"  {lbl:10s} p={p:5d} n={n:5d} lam={lam:5d} lam*={lam_star:.4f}")

SEEDS = [42, 1234, 9999, 31337, 271828]


# ===========================================================================
# E1 — analytic rho_GH, and the m -> infinity limit
# ===========================================================================
print("\n" + "=" * 74)
print("E1  Gaussian-heuristic ratio rho_GH = ||pv|| / GH(L)")
print("=" * 74)
print("Asymptote as m -> inf:  rho_inf = sqrt(2*pi*e/3) * sqrt(eff)"
      f" = {math.sqrt(2*math.pi*math.e/3):.4f} * sqrt(eff)")
print("  (eff = K1*K2/n).  rho_inf = 1 at eff = "
      f"{3/(2*math.pi*math.e):.4f}, i.e. K1 = {3/(2*math.pi*math.e):.4f}*n/K2.")

print("\n  12b/2677 (n=2647, K2=52):  rho_GH vs m and K1")
print("  " + "-" * 68)
n_ = 2647; k2_ = math.isqrt(n_) + 1
hdr = "  K1   eff    " + "".join(f" m={m:<2d}  " for m in (6, 8, 12, 16, 24, 32)) + "  rho_inf"
print(hdr)
for k1 in (2, 3, 4, 6, 8, 12, 16, 24):
    g = geometry(n_, k1, k2_, 12)
    row = f"  {k1:<4d} {g['eff']:.3f} "
    for m in (6, 8, 12, 16, 24, 32):
        row += f" {geometry(n_, k1, k2_, m)['rho_o']:.3f} "
    row += f"   {math.sqrt(2*math.pi*math.e/3)*math.sqrt(g['eff']):.3f}"
    print(row)

print("\n  projected lattice, same curve (rho_p):")
print("  K1   eff    " + "".join(f" m={m:<2d}  " for m in (6, 8, 12, 16, 24, 32)))
for k1 in (2, 3, 4, 6, 8, 12, 16, 24):
    g = geometry(n_, k1, k2_, 12)
    row = f"  {k1:<4d} {g['eff']:.3f} "
    for m in (6, 8, 12, 16, 24, 32):
        row += f" {geometry(n_, k1, k2_, m)['rho_p']:.3f} "
    print(row)

print("\n  Trivial vector n*S_D relative to GH(L), 12b/2677, m=12:")
for k1 in (2, 4, 8, 16):
    S_K1 = n_ // k1; S_K2 = max(1, n_ // k2_)
    dim = 26
    ld = 12 * math.log(n_ * S_K1) + math.log(1) + 12 * math.log(S_K2) + math.log(n_)
    print(f"    K1={k1:<3d}  ||t||/GH = {n_ / gh(ld, dim):.4f}"
          f"   ||pv||/GH = {geometry(n_, k1, k2_, 12)['rho_o']:.4f}")


# ===========================================================================
# E2/E3 — correctness of the projected reconstruction + the sv/pv falsifier
# ===========================================================================
print("\n" + "=" * 74)
print("E2/E3  projected reconstruction + sv/pv (the Thread-23 falsifier)")
print("=" * 74)
print("  curve      K1  m  seed | orig(orcl) orig(free) proj | sv/pv orig  sv/pv proj  proj-sv-is-pv")
print("  " + "-" * 100)

E23 = []
for lbl, k1 in (('8b/199', 2), ('12b/2557', 8), ('12b/2677', 8), ('12b/2677', 4)):
    cp = CURVES[lbl]
    n = cp[2]
    for seed in SEEDS[:3]:
        d_trial = random.Random(seed + 7777).randint(1, n - 1)
        r = run_pair(cp, 12, d_trial, k1, seed)
        if r is None:
            continue
        E23.append((lbl, k1, r))
        print(f"  {lbl:10s} {k1:<3d} 12 {seed:<5d}| "
              f"{str(r['ok_o_oracle']):<10s} {str(r['ok_o_free']):<10s} {str(r['ok_p']):<5s}| "
              f"{r['sv_o']/r['pv_o']:.4f}      {r['sv_p']/r['pv_p']:.4f}      "
              f"{r['shortest_is_planted']}")

agree = sum(1 for _, _, r in E23 if r['ok_o_free'] == r['ok_p'])
print(f"\n  original(oracle-free) and projected agree on {agree}/{len(E23)} instances")
oracle_gap = sum(1 for _, _, r in E23 if r['ok_o_oracle'] != r['ok_o_free'])
print(f"  oracle vs oracle-free disagreement on the original: {oracle_gap}/{len(E23)}")
if E23:
    print(f"  sv/pv projected: min={min(r['sv_p']/r['pv_p'] for _,_,r in E23):.4f} "
          f"max={max(r['sv_p']/r['pv_p'] for _,_,r in E23):.4f}")


# ===========================================================================
# E4 — the T4 K1-wall grid, original vs projected
# ===========================================================================
print("\n" + "=" * 74)
print("E4  K1 wall, m=12, 5 seeds — original vs projected")
print("=" * 74)

K1_GRID = (2, 3, 4, 6, 8, 12, 16, 24)
walls = {}
for lbl in ('12b/2557', '12b/2677'):
    cp = CURVES[lbl]
    n = cp[2]
    lam_star = min(cp[3], n - cp[3]) / n
    k2b = math.isqrt(n) + 1
    print(f"\n  {lbl}  (lam*={lam_star:.3f}, K2={k2b})")
    print("    K1  |  eff   rho_o rho_p | orig(free)  proj  | T4 (2026-07-29, oracle)")
    print("    " + "-" * 74)
    T4 = {'12b/2557': {2:5,3:5,4:5,6:5,8:5,12:4,16:1,24:0},
          '12b/2677': {2:5,3:5,4:5,6:2,8:0,12:0,16:0,24:0}}[lbl]
    row_o, row_p = {}, {}
    for k1 in K1_GRID:
        wo = wp = 0
        for seed in SEEDS:
            d_trial = random.Random(seed + 7777).randint(1, n - 1)
            r = run_pair(cp, 12, d_trial, k1, seed)
            if r is None: continue
            wo += r['ok_o_free']
            wp += r['ok_p']
        row_o[k1], row_p[k1] = wo, wp
        g = geometry(n, k1, k2b, 12)
        print(f"    {k1:<4d}|  {g['eff']:.3f} {g['rho_o']:.3f} {g['rho_p']:.3f} |"
              f"   {wo}/5       {wp}/5  |   {T4[k1]}/5")
    walls[lbl] = (row_o, row_p)

def wall_of(row):
    """largest K1 with 5/5."""
    good = [k for k, v in row.items() if v == 5]
    return max(good) if good else None

print("\n  wall (largest K1 achieving 5/5):")
for lbl, (ro, rp) in walls.items():
    print(f"    {lbl}: original K1<={wall_of(ro)}   projected K1<={wall_of(rp)}")


# ===========================================================================
# E5 — m-sweep at K1=8 on the lam*=0.07 curve (T4b: 0,0,1,0,1 of 5)
# ===========================================================================
print("\n" + "=" * 74)
print("E5  m-sweep at K1=8 on 12b/2677 (lam*=0.070) — does more data rescue it?")
print("=" * 74)
cp = CURVES['12b/2677']; n = cp[2]; k2b = math.isqrt(n) + 1
print("    m  | rho_o rho_p | orig(free)  proj | T4b (orig, 2026-07-29)")
print("    " + "-" * 60)
T4B = {8: 0, 12: 0, 16: 1, 24: 0, 32: 1}
for m in (8, 12, 16, 24, 32):
    wo = wp = 0
    for seed in SEEDS:
        d_trial = random.Random(seed + 7777).randint(1, n - 1)
        r = run_pair(cp, m, d_trial, 8, seed)
        if r is None: continue
        wo += r['ok_o_free']
        wp += r['ok_p']
    g = geometry(n, 8, k2b, m)
    print(f"    {m:<3d}| {g['rho_o']:.3f} {g['rho_p']:.3f} |   {wo}/5       {wp}/5 |   {T4B[m]}/5")


# ===========================================================================
# E6 — reconcile the residual with nu_hat (Thread 20c, 2026-07-29)
# ===========================================================================

def gauss_reduce(u, v):
    """Exact Lagrange-Gauss reduction of a rank-2 integer lattice; returns lambda_1."""
    def dot(a, b): return a[0] * b[0] + a[1] * b[1]
    def n2(a): return dot(a, a)
    if n2(u) > n2(v):
        u, v = v, u
    while True:
        q = (2 * dot(u, v) + n2(u)) // (2 * n2(u))     # round(<u,v>/<u,u>)
        v = (v[0] - q * u[0], v[1] - q * u[1])
        if n2(v) >= n2(u):
            return math.sqrt(n2(u))
        u, v = v, u


def nu_hat(n, lam, k1_bound, k2_bound):
    """nu_hat = lambda_1(L2)/sqrt(det L2) for the non-planted 2-D block L2."""
    S_K1 = n // k1_bound
    S_K2 = max(1, n // k2_bound)
    l1 = gauss_reduce((n * S_K1, 0), (-lam * S_K1, S_K2))
    return l1 / math.sqrt(float(n) * S_K1 * S_K2)


print("\n" + "=" * 74)
print("E6  rho_GH is lam-blind; nu_hat carries the residual")
print("=" * 74)
print("  rho_GH depends only on (n, K1, K2, m) — it cannot see lam. The two")
print("  anchor curves have near-identical rho at every K1 yet walls 8 vs 3.")
print("  nu_hat = lambda_1(L2)/sqrt(det L2) is the lam-dependent part.\n")
print("    curve      K1  |  rho_o(m=12)  nu_hat  | obs 5-seed")
print("    " + "-" * 56)
obs = {'12b/2557': {2:5,3:5,4:5,6:5,8:5,12:3,16:1,24:0},
       '12b/2677': {2:5,3:5,4:4,6:2,8:0,12:0,16:0,24:0}}
for lbl in ('12b/2557', '12b/2677'):
    _, _, n, lam, _ = CURVES[lbl]
    k2b = math.isqrt(n) + 1
    for k1 in K1_GRID:
        print(f"    {lbl:10s} {k1:<3d} |  {geometry(n, k1, k2b, 12)['rho_o']:.3f}"
              f"        {nu_hat(n, lam, k1, k2b):.4f}  |  {obs[lbl][k1]}/5")
    print()

print("  lam* vs nu_hat (K1=8, m=12) — note the ordering INVERTS:")
for lbl in ('12b/2557', '12b/2677', '8b/199'):
    _, _, n, lam, _ = CURVES[lbl]
    k2b = math.isqrt(n) + 1
    k1 = 2 if lbl == '8b/199' else 8
    print(f"    {lbl:10s} lam*={min(lam, n-lam)/n:.4f}  nu_hat(K1={k1})={nu_hat(n, lam, k1, k2b):.4f}")


# ===========================================================================
# E7 — combined predictor  chi = rho_GH * nu_hat
# ===========================================================================
def spearman(xs, ys):
    def rank(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
                j += 1
            avg = (i + j) / 2.0 + 1
            for k in range(i, j + 1):
                r[order[k]] = avg
            i = j + 1
        return r
    rx, ry = rank(xs), rank(ys)
    mx, my = sum(rx) / len(rx), sum(ry) / len(ry)
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    den = math.sqrt(sum((a - mx) ** 2 for a in rx) * sum((b - my) ** 2 for b in ry))
    return num / den if den else float('nan')


print("\n" + "=" * 74)
print("E7  combined predictor  chi = rho_GH * nu_hat   (16 cells, 2 curves)")
print("=" * 74)
rows = []
for lbl in ('12b/2557', '12b/2677'):
    _, _, n, lam, _ = CURVES[lbl]
    k2b = math.isqrt(n) + 1
    for k1 in K1_GRID:
        rho = geometry(n, k1, k2b, 12)['rho_o']
        nh = nu_hat(n, lam, k1, k2b)
        rows.append((rho * nh, rho, nh, lbl, k1, obs[lbl][k1]))
rows.sort()
print("    chi     rho    nu_hat | curve      K1  | wins/5")
print("    " + "-" * 58)
for chi, rho, nh, lbl, k1, w in rows:
    print(f"    {chi:.3f}   {rho:.3f}  {nh:.4f} | {lbl:10s} {k1:<3d} | {w}/5")
w = [r[5] for r in rows]
print(f"\n    Spearman(chi,    wins) = {spearman([r[0] for r in rows], w):+.3f}")
print(f"    Spearman(rho,    wins) = {spearman([r[1] for r in rows], w):+.3f}")
print(f"    Spearman(nu_hat, wins) = {spearman([r[2] for r in rows], w):+.3f}")
print("\n    NOTE: 16 cells over 2 curves only. This is a HYPOTHESIS, not a result.")
print("    Falsifier for the next run: recompute chi on the 20 fresh 17-bit")
print("    curves of T3 (2026-07-29) at eff in {0.05, 0.15, 0.25}. If chi does")
print("    not beat rho and nu_hat individually there, drop the product form.")

print("\nDone.")
