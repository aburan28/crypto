"""
GLV-HNP Thread 23 — reformulating the Phase-2 lattice so the planted vector wins.

Background (RESEARCH_AUTOLAB_LOG.md 2026-07-29, Thread 20 / experiments T1-T5)
------------------------------------------------------------------------------
T5 established that the shortest vector after LLL is, on every instance tested,
the trivial lattice vector  n*S_D*e_m  (norm n), never the planted vector
(norm ~ n*sqrt(2m/3 + 4/3)).  The 2026-07-29 entry proposed Thread 23:
"reformulate the Phase-2 lattice so the target is lambda_1", by projecting along
e_m or by an explicit CVP call.

This script tests a different — and much cheaper — reformulation first, plus the
quantitative model that explains every Phase-2 success/failure observed so far.

The model (C1)
--------------
Recovery is a unique-SVP condition on the *coset* of vectors with last
coordinate +-S_KANNAN.  The natural scale-free statistic is

    r = ||v_planted|| / GH(L),      GH(L) = sqrt(dim/(2*pi*e)) * det(L)^(1/dim)

With the validated column scaling S_K1 = n/K1, S_D = 1, S_K2 = n/K2,
S_KANNAN = n and dim = 2m+2:

    det(L)   = (n*S_K1)^m * S_D * S_K2^m * S_KANNAN = n^(3m+1) / (K1*K2)^m
    ||v||^2  ~ n^2 * (2m/3 + 4/3)                      [uncentered]

so for large m,  det^(1/dim) -> n / sqrt(eff)  with eff = K1*K2/n, and

    r  ->  sqrt(2*pi*e/3) * sqrt(eff)  =  2.386 * sqrt(eff).

Prediction: the wall sits at r = 1, i.e. at eff = 1/2.386^2 = 0.176 — a
*lattice* quantity, not a curve invariant.  That is exactly why the six
curve-level invariants tried 2026-06-21..06-29 (delta/n, kappa(M), q_cf,
max_q_cf, max_a, a_corn/n) all failed to separate the C1/C2 classes.

The reformulation (C2)
----------------------
The Phase-2 lattice never centred its unknowns.  Signatures are generated with
k1 ~ U[0,K1), k2 ~ U[0,K2), and the lattice targets those raw values, so each
scaled coordinate has E[x^2] = n^2/3.  Substituting

    k1 = k1' + K1//2,   k2 = k2' + K2//2,   k1',k2' centred on 0

is a pure change of representation of the SAME instance: it only shifts the
A_i column by  c = K1//2 + lam*(K2//2)  (mod n).  Now E[x^2] = n^2/12, so

    ||v||^2  ~ n^2 * (m/6 + 4/3)      [centred]

which is a factor ~2 shorter for large m, hence  r -> 1.193 * sqrt(eff)  and the
predicted wall moves from eff ~ 0.176 to eff ~ 0.70 — a 4x widening of the
attackable bias regime, at zero cost.

Falsifier: if the centred lattice does NOT move the K1 wall outward on the
lam*=0.070 curve (currently K1 ~ 4-6), the GH model is wrong and the wall is
information-theoretic.

Experiments
-----------
  C1  closed-form vs measured r; reproduction of the T4 grid (uncentred)
  C2  T4 grid, centred vs uncentred, on both historical 12-bit curves
  C3  fresh 17-bit sweep across eff, centred vs uncentred; is r a separator?
  C4  is the planted vector lambda_1 after centring?  (+ S_KANNAN = n/2 variant)

Run: python3 glv_hnp_phase2_centered.py
"""

import math
import random

import sympy
from fpylll import IntegerMatrix, LLL, BKZ

TWO_PI_E = 2.0 * math.pi * math.e

# ---------------------------------------------------------------------------
# Minimal EC arithmetic (short Weierstrass, a=0) — same as Thread 20 scripts
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

def tonelli_shanks(a, p):
    a %= p
    if a == 0: return 0
    if pow(a, (p - 1) // 2, p) != 1: return None
    if p % 4 == 3: return pow(a, (p + 1) // 4, p)
    q, s = p - 1, 0
    while q % 2 == 0: q //= 2; s += 1
    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1: z += 1
    mm, c, t, r = s, pow(z, q, p), pow(a, q, p), pow(a, (q + 1) // 2, p)
    while True:
        if t == 0: return 0
        if t == 1: return r
        i, tmp = 0, t
        while tmp != 1: tmp = tmp * tmp % p; i += 1
        b = pow(c, 1 << (mm - i - 1), p)
        mm, c, t, r = i, b * b % p, t * b * b % p, r * b % p

def find_generator(p, b, n, seed=12345):
    rng = random.Random(seed)
    for _ in range(10000):
        x = rng.randint(0, p - 1)
        y = tonelli_shanks((pow(x, 3, p) + b) % p, p)
        if y is not None and y != 0 and ec_mul((x, y), n, p) is None:
            return (x, y)
    return None

# ---------------------------------------------------------------------------
# j=0 CM curve construction
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
    return min(lam % n, n - (lam % n)) / n

def build_curve(p, n, seed=12345):
    rng = random.Random(seed)
    for _ in range(400):
        b = rng.randint(1, p - 1)
        x = rng.randint(0, p - 1)
        y = tonelli_shanks((pow(x, 3, p) + b) % p, p)
        if y is None or y == 0:
            continue
        if ec_mul((x, y), n, p) is None:
            G = find_generator(p, b, n, seed=seed)
            if G is not None:
                return (p, b, n, None, G)
    return None

def search_curves(lo, hi, per_bin=2, nbins=10, max_primes=100000):
    """j=0 GLV curves with p in [lo,hi), n prime, n = 1 mod 3, bucketed by lam*."""
    bins = {i: [] for i in range(nbins)}
    p = int(sympy.nextprime(lo))
    seen = 0
    while p < hi and seen < max_primes:
        seen += 1
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
# Phase-2 lattice, with the centring reformulation as a flag
# ---------------------------------------------------------------------------

def scales(n, k1_bound, k2_bound, kan_frac=1.0):
    """Column-diagonal scaling validated 2026-07-26.

    kan_frac generalises S_KANNAN = kan_frac * n (the 2026-07-26 builder is
    kan_frac = 1.0).  C5 shows 1.0 is not optimal once the unknowns are centred.
    """
    return (n // k1_bound, 1, max(1, n // k2_bound), max(1, int(kan_frac * n)))

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

def centre_shift(lam, k1_bound, k2_bound, n, centred):
    """c such that A_i - c is the constant term for the centred unknowns."""
    if not centred:
        return 0
    return (k1_bound // 2 + lam * (k2_bound // 2)) % n

def build_glv_lattice(sigs, n, lam, k1_bound, k2_bound, centred=False,
                      kan_frac=1.0):
    """Identical to the 2026-07-26 builder except for the A-column shift."""
    m = len(sigs)
    dim = 2 * m + 2
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound, kan_frac)
    c = centre_shift(lam, k1_bound, k2_bound, n, centred)
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
        M[2 * m + 1][i] = ((sigs[i]['A'] - c) % n) * S_K1
    M[2 * m + 1][dim - 1] = S_KAN
    return M

def planted_vector(sigs, d_secret, n, lam, k1_bound, k2_bound, centred=False,
                   kan_frac=1.0):
    m = len(sigs)
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound, kan_frac)
    o1 = k1_bound // 2 if centred else 0
    o2 = k2_bound // 2 if centred else 0
    c = centre_shift(lam, k1_bound, k2_bound, n, centred)
    v = [0] * (2 * m + 2)
    for i in range(m):
        k1c = sigs[i]['k1'] - o1
        k2c = sigs[i]['k2'] - o2
        # membership identity: A_i - c + B_i*d - lam*k2c = k1c  (mod n)
        assert ((sigs[i]['A'] - c) + sigs[i]['B'] * d_secret
                - lam * k2c - k1c) % n == 0
        v[i] = k1c * S_K1
        v[m + 1 + i] = k2c * S_K2
    v[m] = d_secret * S_D
    v[2 * m + 1] = S_KAN
    return v

def norm(v):
    return math.sqrt(sum(x * x for x in v))

def gaussian_heuristic(n, m, k1_bound, k2_bound, kan_frac=1.0):
    """GH(L) = sqrt(dim/(2*pi*e)) * det(L)^(1/dim), computed in logs."""
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound, kan_frac)
    dim = 2 * m + 2
    log_det = (m * math.log(n * S_K1) + math.log(S_D)
               + m * math.log(S_K2) + math.log(S_KAN))
    return math.sqrt(dim / TWO_PI_E) * math.exp(log_det / dim)

def planted_norm_expected(n, m, k1_bound, k2_bound, centred=False,
                          kan_frac=1.0):
    """E[||v||] over k1 ~ U[0,K1), k2 ~ U[0,K2), d ~ U[1,n)."""
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound, kan_frac)
    div = 12.0 if centred else 3.0
    e_k1sq = k1_bound ** 2 / div
    e_k2sq = k2_bound ** 2 / div
    e_dsq = n ** 2 / 3.0            # d is never centred (it is the unknown)
    return math.sqrt(m * e_k1sq * S_K1 ** 2 + e_dsq * S_D ** 2
                     + m * e_k2sq * S_K2 ** 2 + S_KAN ** 2)

def r_predicted(n, m, k1_bound, k2_bound, centred=False, kan_frac=1.0):
    return (planted_norm_expected(n, m, k1_bound, k2_bound, centred, kan_frac)
            / gaussian_heuristic(n, m, k1_bound, k2_bound, kan_frac))

def recover_d(reduced, m, n, S_KAN, d_secret):
    dim = 2 * m + 2
    for row in reduced:
        last = row[dim - 1]
        if abs(last) != S_KAN: continue
        sign = 1 if last > 0 else -1
        d_cand = (sign * row[m]) % n
        if d_cand != 0 and d_cand == d_secret:
            return d_cand
    return None

def run_experiment(curve, m, d_secret, k1_bound, seed=42, centred=False,
                   kan_frac=1.0, use_bkz=False, bkz_beta=20):
    """Returns (recovered, ||v_planted||, ||shortest after reduction||, GH)."""
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return None
    M = build_glv_lattice(sigs, n, lam, k1_bound, k2_bound, centred, kan_frac)
    dim = 2 * m + 2
    A = IntegerMatrix.from_matrix(M)
    if use_bkz:
        BKZ.reduction(A, BKZ.Param(bkz_beta))
    else:
        LLL.reduction(A)
    reduced = [[A[i][j] for j in range(dim)] for i in range(dim)]
    _, _, _, S_KAN = scales(n, k1_bound, k2_bound, kan_frac)
    ok = recover_d(reduced, m, n, S_KAN, d_secret) is not None
    pn = norm(planted_vector(sigs, d_secret, n, lam, k1_bound, k2_bound,
                             centred, kan_frac))
    sn = min(norm(r) for r in reduced)
    return (ok, pn, sn, gaussian_heuristic(n, m, k1_bound, k2_bound, kan_frac))

def success_rate(curve, m, k1_bound, seeds, centred=False, kan_frac=1.0,
                 use_bkz=False, bkz_beta=20):
    """Returns (wins, total, mean r_measured, mean sv/pv)."""
    wins, rs, svpv = 0, [], []
    n = curve[2]
    for seed in seeds:
        d_trial = random.Random(seed + 7777).randint(1, n - 1)
        res = run_experiment(curve, m, d_trial, k1_bound, seed, centred,
                             kan_frac, use_bkz, bkz_beta)
        if res is None:
            continue
        ok, pn, sn, gh = res
        wins += bool(ok)
        rs.append(pn / gh)
        svpv.append(sn / pn if pn else float('nan'))
    mr = sum(rs) / len(rs) if rs else float('nan')
    ms = sum(svpv) / len(svpv) if svpv else float('nan')
    return wins, len(seeds), mr, ms

def auc(rows):
    """rows = [(score, label)]; AUC for 'lower score => positive label'."""
    pos = [s for s, y in rows if y]
    neg = [s for s, y in rows if not y]
    if not pos or not neg:
        return float('nan')
    wins = sum((1.0 if a < b else 0.5 if a == b else 0.0)
               for a in pos for b in neg)
    return wins / (len(pos) * len(neg))


# ===========================================================================
print("=" * 78)
print("Thread 23 — centred Phase-2 lattice and the Gaussian-heuristic wall")
print("=" * 78)

SEEDS = [42, 1234, 9999, 555, 31337]

HIST = [
    # label,             p,    b, n,    lam,  m
    ("8-bit/199",        211,  2, 199,  106,  6),
    ("12-bit/2557",      2557, 2, 2659, 1755, 12),
    ("12-bit/2677 FAIL", 2677, 2, 2647, 185,  12),
]

hist = []
for label, p, b, n, lam, m in HIST:
    G = find_generator(p, b, n)
    assert G is not None, f"no generator for {label}"
    assert (lam * lam + lam + 1) % n == 0, f"bad lam for {label}"
    hist.append((label, (p, b, n, lam, G), m))


# ===========================================================================
# C1 — the closed form
# ===========================================================================
print("\n" + "-" * 78)
print("C1: r = ||v_planted|| / GH(L) — closed form vs exact per-instance value")
print("-" * 78)
print("Asymptotic prediction (m -> inf):  r = sqrt(2*pi*e/3)*sqrt(eff) = "
      f"{math.sqrt(TWO_PI_E / 3):.3f}*sqrt(eff)   [uncentred]")
print("                                   r = sqrt(2*pi*e/12)*sqrt(eff) = "
      f"{math.sqrt(TWO_PI_E / 12):.3f}*sqrt(eff)  [centred]\n")

print(f"{'curve':<18} {'m':>3} {'K1':>4} {'K2':>5} {'eff':>6} "
      f"{'r_un':>7} {'2.386v':>7} {'r_ctr':>7} {'1.193v':>7}")
for label, curve, m in hist:
    n = curve[2]
    k2 = math.isqrt(n) + 1
    for k1 in (2, 8):
        eff = k1 * k2 / n
        ru = r_predicted(n, m, k1, k2, centred=False)
        rc = r_predicted(n, m, k1, k2, centred=True)
        print(f"{label:<18} {m:>3} {k1:>4} {k2:>5} {eff:>6.3f} "
              f"{ru:>7.3f} {math.sqrt(TWO_PI_E/3*eff):>7.3f} "
              f"{rc:>7.3f} {math.sqrt(TWO_PI_E/12*eff):>7.3f}")
print("\n(the closed form is the m -> inf limit; the small-m values carry the")
print(" +4/3*n^2 from the d and Kannan coordinates, so r_exact > r_asymptotic)")


# ===========================================================================
# C2 — the T4 grid, centred vs uncentred (the falsifier)
# ===========================================================================
print("\n" + "-" * 78)
print("C2: K1 grid — does centring move the wall?   (m=12, 5 seeds)")
print("-" * 78)
print("Falsifier: if the centred row does not extend past the uncentred row on")
print("the lam*=0.070 curve, the GH model is wrong.\n")

K1_GRID = [2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64]
c2_rows = []
for label, curve, m in hist[1:]:            # the two 12-bit curves from T4
    p, b, n, lam, G = curve
    k2 = math.isqrt(n) + 1
    print(f"\n{label}   lam*={lam_star(lam,n):.3f}  n={n}  K2={k2}  m={m}")
    print(f"{'K1':>5} {'eff':>6} | {'un: LLL':>8} {'r':>6} {'sv/pv':>6} "
          f"| {'ctr: LLL':>9} {'r':>6} {'sv/pv':>6}")
    for k1 in K1_GRID:
        eff = k1 * k2 / n
        if eff >= 1.0:
            continue
        wu, tu, ru, su = success_rate(curve, m, k1, SEEDS, centred=False)
        wc, tc, rc, sc = success_rate(curve, m, k1, SEEDS, centred=True)
        c2_rows.append((label, k1, eff, wu, tu, ru, wc, tc, rc))
        print(f"{k1:>5} {eff:>6.3f} | {str(wu)+'/'+str(tu):>8} {ru:>6.3f} "
              f"{su:>6.3f} | {str(wc)+'/'+str(tc):>9} {rc:>6.3f} {sc:>6.3f}")

print("\nLargest eff with a full 5/5 sweep:")
for label, _, _ in hist[1:]:
    best_u = max([e for (l, k, e, wu, tu, ru, wc, tc, rc) in c2_rows
                  if l == label and wu == tu] or [0.0])
    best_c = max([e for (l, k, e, wu, tu, ru, wc, tc, rc) in c2_rows
                  if l == label and wc == tc] or [0.0])
    gain = best_c / best_u if best_u else float('inf')
    print(f"  {label:<18} uncentred eff={best_u:.3f}   centred eff={best_c:.3f}"
          f"   gain x{gain:.2f}")


# ===========================================================================
# C3 — fresh 17-bit sweep across eff, both variants
# ===========================================================================
print("\n" + "-" * 78)
print("C3: fresh 17-bit curves, eff sweep, centred vs uncentred (m=12)")
print("-" * 78)

LO, HI = 2 ** 16, 2 ** 17
M_SIGS = 12
print(f"Searching j=0 GLV curves with p in [{LO}, {HI}) ...")
curves = search_curves(LO, HI, per_bin=2, nbins=10)
print(f"Found {len(curves)} curves (lam* spread over (0, 0.5)).\n")

EFF_TARGETS = [0.05, 0.15, 0.25, 0.40, 0.60, 0.80]
c3_rows = []
print(f"{'eff*':>6} | {'uncentred wins':>15} {'mean r':>7} "
      f"| {'centred wins':>13} {'mean r':>7}")
for eff_t in EFF_TARGETS:
    wu_tot = wc_tot = tot = 0
    ru_all, rc_all = [], []
    for (p, b, n, lam, G) in curves:
        k2 = math.isqrt(n) + 1
        k1 = max(2, int(eff_t * n / k2))
        if k1 * k2 / n >= 1.0:
            continue
        curve = (p, b, n, lam, G)
        wu, tu, ru, _ = success_rate(curve, M_SIGS, k1, SEEDS, centred=False)
        wc, tc, rc, _ = success_rate(curve, M_SIGS, k1, SEEDS, centred=True)
        wu_tot += wu; wc_tot += wc; tot += tu
        ru_all.append(ru); rc_all.append(rc)
        c3_rows.append((lam_star(lam, n), k1 * k2 / n, ru, wu == tu,
                        rc, wc == tc))
    mru = sum(ru_all) / len(ru_all) if ru_all else float('nan')
    mrc = sum(rc_all) / len(rc_all) if rc_all else float('nan')
    print(f"{eff_t:>6.2f} | {str(wu_tot)+'/'+str(tot):>15} {mru:>7.3f} "
          f"| {str(wc_tot)+'/'+str(tot):>13} {mrc:>7.3f}")

# Is r a single separator across BOTH variants and all eff?
pooled = ([(r, ok) for (_, _, r, ok, _, _) in c3_rows]
          + [(r, ok) for (_, _, _, _, r, ok) in c3_rows])
pooled_eff = ([(e, ok) for (_, e, _, ok, _, _) in c3_rows]
              + [(e, ok) for (_, e, _, _, _, ok) in c3_rows])
pooled_lam = ([(l, ok) for (l, _, _, ok, _, _) in c3_rows]
              + [(l, ok) for (l, _, _, _, _, ok) in c3_rows])
npos = sum(1 for _, y in pooled if y)
print(f"\nPooled instances (both variants): {len(pooled)}, "
      f"full-sweep successes: {npos}")
print(f"  AUC(r)     = {auc(pooled):.3f}      <- lattice statistic")
print(f"  AUC(eff)   = {auc(pooled_eff):.3f}      <- ignores centring")
print(f"  AUC(lam*)  = {auc(pooled_lam):.3f}      <- curve invariant (T3: none work)")

succ_r = [r for r, y in pooled if y]
fail_r = [r for r, y in pooled if not y]
if succ_r and fail_r:
    print(f"  r | success: max = {max(succ_r):.3f}, mean = "
          f"{sum(succ_r)/len(succ_r):.3f}")
    print(f"  r | failure: min = {min(fail_r):.3f}, mean = "
          f"{sum(fail_r)/len(fail_r):.3f}")


# ===========================================================================
# C4 — is the planted vector lambda_1 now?  And does S_KANNAN = n/2 help?
# ===========================================================================
print("\n" + "-" * 78)
print("C4: sv/pv after centring, and the S_KANNAN = n/2 variant")
print("-" * 78)
print("T5 (2026-07-29): uncentred, the shortest vector is always the trivial")
print("n*S_D*e_m, so sv/pv < 1 always.  Centring shrinks ||v|| while leaving the")
print("trivial vector untouched, so sv/pv rises toward 1 but cannot reach it:")
print("the planted vector is still NOT lambda_1.  Recovery does not need it to")
print("be — it needs the planted vector to be the shortest in its Kannan coset.\n")

print(f"{'curve':<18} {'K1':>4} {'variant':<22} {'LLL':>6} {'r':>6} {'sv/pv':>6}")
for label, curve, m in hist:
    n = curve[2]
    k2 = math.isqrt(n) + 1
    k1 = 8 if n > 1000 else 2
    for name, kw in (("uncentred, S_KAN=n", dict(centred=False)),
                     ("centred,   S_KAN=n", dict(centred=True)),
                     ("centred,   S_KAN=n/2", dict(centred=True, kan_frac=0.5))):
        w, t, r, s = success_rate(curve, m, k1, SEEDS, **kw)
        print(f"{label:<18} {k1:>4} {name:<22} {str(w)+'/'+str(t):>6} "
              f"{r:>6.3f} {s:>6.3f}")


# ===========================================================================
# C5 — the optimal Kannan scale
# ===========================================================================
print("\n" + "-" * 78)
print("C5: optimal S_KANNAN = t*n, and the wall with the fully tuned lattice")
print("-" * 78)
print("Minimising r(t) = sqrt(C + (t*n)^2) / (const * (t*n)^(1/dim)) over t gives")
print("t*^2*n^2 = C/(dim-1) with C = m*(E[k1^2]*S_K1^2 + E[k2^2]*S_K2^2) + n^2/3.")
print("Centred, large m:  C -> n^2*(m/6 + 1/3), dim-1 = 2m+1, so t* -> 1/sqrt(12)")
print(f"= {1/math.sqrt(12):.3f}.  The 2026-07-26 builder uses t = 1.0.\n")


def t_star(n, m, k1_bound, k2_bound, centred):
    S_K1, S_D, S_K2, _ = scales(n, k1_bound, k2_bound, 1.0)
    div = 12.0 if centred else 3.0
    C = (m * (k1_bound ** 2 / div) * S_K1 ** 2
         + m * (k2_bound ** 2 / div) * S_K2 ** 2 + (n ** 2 / 3.0) * S_D ** 2)
    return math.sqrt(C / (2 * m + 1)) / n


T_GRID = [0.15, 0.20, 0.25, 0.29, 0.35, 0.50, 0.75, 1.00, 1.50]
for label, curve, m in hist[1:]:
    p, b, n, lam, G = curve
    k2 = math.isqrt(n) + 1
    k1 = 16
    print(f"\n{label}  K1={k1}  eff={k1*k2/n:.3f}  m={m}  "
          f"t* (predicted) = {t_star(n, m, k1, k2, True):.3f}")
    print(f"{'t':>6} {'r':>7} {'LLL':>6} {'sv/pv':>7}")
    for t in T_GRID:
        w, tt, r, s = success_rate(curve, m, k1, SEEDS, centred=True, kan_frac=t)
        print(f"{t:>6.2f} {r:>7.3f} {str(w)+'/'+str(tt):>6} {s:>7.3f}")

print("\nK1 grid with the tuned lattice (centred, t = 0.29):")
print(f"{'curve':<18} {'K1':>4} {'eff':>6} {'r':>6} {'LLL':>6}")
for label, curve, m in hist[1:]:
    p, b, n, lam, G = curve
    k2 = math.isqrt(n) + 1
    for k1 in (16, 24, 32, 48, 64):
        eff = k1 * k2 / n
        if eff >= 1.0:
            continue
        w, t, r, s = success_rate(curve, m, k1, SEEDS, centred=True,
                                  kan_frac=0.29)
        print(f"{label:<18} {k1:>4} {eff:>6.3f} {r:>6.3f} "
              f"{str(w)+'/'+str(t):>6}")


# ===========================================================================
# C6 — the June C1/C2 split, re-run on the centred lattice
# ===========================================================================
print("\n" + "-" * 78)
print("C6: does the C1/C2 split survive centring?   (Exp S protocol regime)")
print("-" * 78)
print("Exp S (2026-06-29) fixed K1=72, m=12, 20-bit curves, and called a curve")
print("C2 if LLL recovers d on every seed, C1 otherwise.  74% were C1, and six")
print("curve invariants failed to separate them; nu_hat reached AUC 0.935.")
print("The GH model says Exp S sits *exactly at the wall* (r ~ 1), which is why a")
print("second-order statistic could separate at all.  Centring drops r well below")
print("1, so the prediction is that the C1 class disappears entirely.\n")

EXPS_K1, EXPS_M = 72, 12
EXPS_SEEDS = [42, 1234, 9999, 555, 31337, 271828]
print(f"Harvesting 20-bit j=0 GLV curves ...")
big = search_curves(2 ** 19, 2 ** 20, per_bin=10, nbins=10)
print(f"Found {len(big)} curves.\n")

c2_un = c2_ct = 0
wins_un = wins_ct = trials = 0
r_un_all, r_ct_all = [], []
for (p, b, n, lam, G) in big:
    k2 = math.isqrt(n) + 1
    curve = (p, b, n, lam, G)
    wu, tu, ru, _ = success_rate(curve, EXPS_M, EXPS_K1, EXPS_SEEDS, centred=False)
    wc, tc, rc, _ = success_rate(curve, EXPS_M, EXPS_K1, EXPS_SEEDS, centred=True,
                                 kan_frac=0.29)
    c2_un += (wu == tu); c2_ct += (wc == tc)
    wins_un += wu; wins_ct += wc; trials += tu
    r_un_all.append(ru); r_ct_all.append(rc)

nc = len(big)
print(f"{'lattice':<28} {'C2 curves':>10} {'seed wins':>12} {'mean r':>8}")
print(f"{'uncentred, t=1 (Exp S)':<28} {str(c2_un)+'/'+str(nc):>10} "
      f"{str(wins_un)+'/'+str(trials):>12} "
      f"{sum(r_un_all)/len(r_un_all):>8.3f}")
print(f"{'centred,   t=0.29':<28} {str(c2_ct)+'/'+str(nc):>10} "
      f"{str(wins_ct)+'/'+str(trials):>12} "
      f"{sum(r_ct_all)/len(r_ct_all):>8.3f}")
print(f"\nK1={EXPS_K1}, m={EXPS_M}, K2 ~ sqrt(n), eff ~ {EXPS_K1/2**10:.3f}")

print("\nDone.")
