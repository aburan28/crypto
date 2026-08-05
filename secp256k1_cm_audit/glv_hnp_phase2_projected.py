"""
GLV-HNP Phase 2 — Thread 23: reformulate the lattice so the planted vector
is lambda_1 (or a genuine CVP target), and test whether that moves the K1 wall.

Background (RESEARCH_AUTOLAB_LOG.md, 2026-07-29, EXP T5)
--------------------------------------------------------
The Phase-2 lattice L subset Z^{2m+2} (glv_hnp_phase2_20bit.py:262) always
contains the *trivial* vector

    w0 = n * S_D * e_m           (n * row_m  -  sum_i B_i * row_i)

of norm n*S_D = n, while the planted vector has expected norm
n * sqrt(2m/3 + 4/3) > n.  T5 measured ||sv||/||v_planted|| in [0.34, 0.61] on
every curve, success and failure alike, with 100% of the shortest vector's
energy in the d-column.  So the planted vector is NEVER lambda_1 of L, and
recovery is a BDD/coset condition, not an SVP condition.

w0 carries no information (d is only defined mod n) and no choice of S_D
removes it: both ||w0|| and ||v_planted|| scale linearly in S_D.

Thread 23 (proposed 2026-07-29)
-------------------------------
Quotient the trivial direction out and re-pose the problem.  Two variants,
both named in the proposal:

  P1  "projected SVP":  pi = delete coordinate m.  pi(L) subset Z^{2m+1} has
      rank 2m+1, det pi(L) = det(L) / (n*S_D) = n^m * S_K1^m * S_K2^m, and
      pi(w0) = 0.  The image of the 2m+2 basis rows generates pi(L) (the
      d-row does NOT become redundant: it contributes B_i*S_K1 with B_i < n).
      pi(v_planted) = (k1_i*S_K1 | k2_i*S_K2 | S_KANNAN), norm ~ n*sqrt(2m/3+1).
      d is still recoverable: d = (k1_0 + lam*k2_0 - A_0) * B_0^{-1} mod n.

  P2  "explicit CVP":  drop the Kannan column as well.  L'' subset Z^{2m}
      spanned by {n*S_K1*e_i}, {(B_i*S_K1 | 0)}, {(-lam*S_K1*e_i | S_K2*e_i)}.
      Target t = (-A_i*S_K1 | 0).  The closest vector w satisfies
      w - t = (k1_i*S_K1 | k2_i*S_K2).  Solved by Babai nearest-plane on the
      LLL-reduced basis (GSO.Mat.babai), no embedding.

Falsifier (stated verbatim in the 2026-07-29 log)
-------------------------------------------------
  "if sv/pv rises above 1 after the reformulation AND the K1 wall in T4 moves
   outward on the lam*=0.07 curve (currently K1 ~ 4-6), the reformulation is a
   real improvement; if the wall stays at K1 ~ 4-6, then the wall is
   information-theoretic and Phase 2 is at its ceiling."

Experiments
-----------
  E1  Structure: is pi(v_planted) = lambda_1 of pi(L)?  Report ||sv||/||pi(v)||
      and the energy split, mirroring T5's table so the two are comparable.
  E2  THE FALSIFIER: the T4 K1 grid (K1 in {2,3,4,6,8,12,16,24}, m=12, 5 seeds)
      re-run for all three methods (ORIG / P1 / P2) on the historical pair.
  E2b More data at the wall (T4b replication) for P1 and P2.
  E3  Gaussian-heuristic prediction: does ||pi(v)|| / GH(pi(L)) predict the eff
      wall (T3 measured eff=0.05 works, 0.15 marginal, 0.25 fails)?
  E4  Fresh 17-bit curves at eff in {0.05, 0.15, 0.25}: ORIG vs P1 vs P2.

RESULTS (2026-08-05 run; full output in glv_hnp_phase2_projected_output.txt)
---------------------------------------------------------------------------
  E2  THREAD 23 ANSWERED NEGATIVELY.  The K1 wall does not move.

        curve              ORIG  P1   P2
        12-bit/2557 (0.34)   12   12    8
        12-bit/2677 (0.07)    4    4    4     <- wall = largest K1 with >=3/5

      P1 reproduces ORIG cell for cell across the whole 8x2 grid; P2 (Babai,
      no embedding) is strictly worse on the easy curve and equal on the hard
      one.  E2b: more data does not rescue either (m=8..32 at K1=8 -> 0-1/5).

  E1  The planted vector is STILL not lambda_1 after projection: sv/pv =
      0.843 / 0.532 / 0.813 (vs 0.603 / 0.517 / 0.422 unprojected).  The new
      shortest vector again has Kannan-energy 0.000 -- it is a different
      parasite, split across the k1/k2 blocks instead of sitting in the
      d-column.  Removing w0 promotes the next non-Kannan vector; the
      obstruction is the whole 2m-dimensional kan=0 sublattice, not one vector.

  E5  ORIG and P1 agree on 119/120 individual (curve, K1, seed) trials
      (99.2%); the single disagreement is a P1 *win* at K1=16 on the 8-bit
      curve, i.e. past the wall.  So the projection is a no-op for LLL.
      This does not contradict T5 -- T5's observation that the planted vector
      is never lambda_1 is correct -- but it shows that fact is NOT CAUSAL
      for failure.

  E3  ||pi(v)||/GH(pi(L)) -> sqrt(2*pi*e*eff/3) asymptotically, giving a
      closed-form wall at eff* = 3/(2*pi*e) = 0.1947.  That brackets the
      measured band (T3: 0.05 works, 0.15 marginal, 0.25 dead; E4 here:
      6/6, 2/6, 0/6 curves at 3/3).  But it is NOT a separator: 12-bit/2557
      and 12-bit/2677 at K1=8 have the SAME eff (0.156 / 0.157) and the SAME
      exact ratio (1.251) yet score 5/5 and 0/5.  eff fixes where the wall
      is; something lam-dependent decides the marginal band.

  E6  That something is nu_hat (2026-07-29, commit e845207), replicated
      OUT OF SAMPLE on 14 fresh 17-bit curves held at eff = 0.149:
        AUC(nu_hat, low = success) = 0.958  (23/24 concordant pairs)
        AUC(lam*,   low = success) = 0.333  (worse than chance)
        successes mean nu_hat 0.423 vs failures 0.800.
      Compare the 2026-07-29 in-sample AUC of 0.935.

  Consequence: Phase 2 is at its reduction ceiling.  The remaining lever is
  eff = K1*K2/n (an information-theoretic quantity fixed by the signature
  bias), not the lattice's presentation.

Everything below the helper block is new; the helper block is copied verbatim
from glv_hnp_phase2_lambda_threshold.py:87-350 (which itself copied
glv_hnp_phase2_20bit.py:262) so the ORIG column reproduces T4 exactly.
That file executes its whole T1-T5 suite at import time, so it is copied
rather than imported.

Run: python3 glv_hnp_phase2_projected.py
"""

import math
import random

import sympy
from fpylll import IntegerMatrix, LLL, GSO

# ---------------------------------------------------------------------------
# Helper block — verbatim from glv_hnp_phase2_lambda_threshold.py
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
    for _ in range(20000):
        x = rng.randint(0, p - 1)
        rhs = (pow(x, 3, p) + b) % p
        y = tonelli_shanks(rhs, p)
        if y is not None and y != 0:
            P = (x, y)
            if ec_mul(P, n, p) is None:
                return P
    return None

def eisenstein_decompose(p):
    """(a,b) with a^2 - a*b + b^2 = p, a,b >= 0.  O(sqrt p)."""
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

def norm(v):
    return math.sqrt(sum(x * x for x in v))

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
# NEW: the two reformulations
# ---------------------------------------------------------------------------

def lll_rows(M, ncols):
    """LLL-reduce a (possibly rank-deficient) generating set; drop zero rows."""
    A = IntegerMatrix.from_matrix(M)
    LLL.reduction(A)
    rows = [[A[i][j] for j in range(ncols)] for i in range(A.nrows)]
    return [r for r in rows if any(r)]


def project_out_d(M, m):
    """pi: delete coordinate m.  Returns the (2m+2) x (2m+1) generating set of
    pi(L).  pi(w0) = 0, so the image has rank 2m+1."""
    return [row[:m] + row[m + 1:] for row in M]


def planted_projected(sigs, n, k1_bound, k2_bound):
    """pi(v_planted) = (k1_i*S_K1 | k2_i*S_K2 | S_KANNAN), dim 2m+1."""
    m = len(sigs)
    S_K1, _S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    v = [0] * (2 * m + 1)
    for i in range(m):
        v[i] = sigs[i]['k1'] * S_K1
        v[m + i] = sigs[i]['k2'] * S_K2
    v[2 * m] = S_KAN
    return v


def d_from_k(k1, k2, sig, n, lam):
    """A + B*d = k1 + lam*k2 (mod n)  =>  d = (k1 + lam*k2 - A) * B^{-1}."""
    try:
        return (k1 + lam * k2 - sig['A']) * modinv(sig['B'], n) % n
    except ValueError:
        return None


def recover_d_projected(rows, sigs, n, k1_bound, k2_bound, lam, d_secret):
    """P1 recovery: any reduced row whose last coordinate is +-S_KANNAN is a
    candidate for +-pi(v_planted); read (k1_i, k2_i) off it and solve for d."""
    m = len(sigs)
    S_K1, _S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    for row in rows:
        last = row[2 * m]
        if abs(last) != S_KAN:
            continue
        sgn = 1 if last > 0 else -1
        for i in range(m):
            k1 = sgn * row[i]
            k2 = sgn * row[m + i]
            if k1 % S_K1 or k2 % S_K2:
                continue
            d_cand = d_from_k(k1 // S_K1, k2 // S_K2, sigs[i], n, lam)
            if d_cand is not None and d_cand == d_secret:
                return d_cand
    return None


def build_cvp_lattice(sigs, n, lam, k1_bound, k2_bound):
    """P2: L'' subset Z^{2m}, columns = (k1-block | k2-block).  No d column,
    no Kannan column.  Generators: m modular rows, the d-row, m lam-rows."""
    m = len(sigs)
    S_K1, _S_D, S_K2, _S_KAN = scales(n, k1_bound, k2_bound)
    gens = []
    for i in range(m):
        r = [0] * (2 * m)
        r[i] = n * S_K1
        gens.append(r)
    r = [0] * (2 * m)
    for i in range(m):
        r[i] = sigs[i]['B'] * S_K1
    gens.append(r)
    for i in range(m):
        r = [0] * (2 * m)
        r[i] = -lam * S_K1
        r[m + i] = S_K2
        gens.append(r)
    return gens


def cvp_target(sigs, n, k1_bound, k2_bound):
    m = len(sigs)
    S_K1, _S_D, _S_K2, _S_KAN = scales(n, k1_bound, k2_bound)
    return [-sigs[i]['A'] * S_K1 for i in range(m)] + [0] * m


def babai_nearest_plane(basis_rows, target):
    """Babai nearest-plane on an already LLL-reduced basis.  Returns the
    lattice vector.  Uses fpylll's GSO.Mat.babai for the coefficients."""
    dim = len(target)
    A = IntegerMatrix.from_matrix(basis_rows)
    Mg = GSO.Mat(A, float_type="mpfr")
    Mg.update_gso()
    coeffs = Mg.babai(list(target))
    out = [0] * dim
    for i, c in enumerate(coeffs):
        if c:
            for j in range(dim):
                out[j] += c * basis_rows[i][j]
    return out


def recover_d_cvp(sigs, n, lam, k1_bound, k2_bound, d_secret):
    """P2 recovery: Babai CVP, then read (k1_i,k2_i) off (w - t)."""
    m = len(sigs)
    S_K1, _S_D, S_K2, _S_KAN = scales(n, k1_bound, k2_bound)
    gens = build_cvp_lattice(sigs, n, lam, k1_bound, k2_bound)
    basis = lll_rows(gens, 2 * m)
    if len(basis) != 2 * m:
        return None, None
    t = cvp_target(sigs, n, k1_bound, k2_bound)
    w = babai_nearest_plane(basis, t)
    diff = [w[j] - t[j] for j in range(2 * m)]
    for i in range(m):
        if diff[i] % S_K1 or diff[m + i] % S_K2:
            continue
        d_cand = d_from_k(diff[i] // S_K1, diff[m + i] // S_K2, sigs[i], n, lam)
        if d_cand is not None and d_cand == d_secret:
            return d_cand, diff
    return None, diff


# ---------------------------------------------------------------------------
# Unified trial driver: ORIG / P1 / P2 on the same signature set
# ---------------------------------------------------------------------------

def trial(curve, m, d_secret, k1_bound, seed, method):
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return None
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)

    if method == "ORIG":
        M = build_glv_lattice(sigs, n, lam, k1_bound, k2_bound)
        rows = lll_rows(M, 2 * m + 2)
        ok = recover_d(rows, m, n, S_KAN, d_secret) is not None
        pv = norm(planted_vector(sigs, d_secret, n, k1_bound, k2_bound))
        sv = min(norm(r) for r in rows)
        return ok, sv, pv, rows, sigs

    if method == "P1":
        M = build_glv_lattice(sigs, n, lam, k1_bound, k2_bound)
        rows = lll_rows(project_out_d(M, m), 2 * m + 1)
        ok = recover_d_projected(rows, sigs, n, k1_bound, k2_bound,
                                 lam, d_secret) is not None
        pv = norm(planted_projected(sigs, n, k1_bound, k2_bound))
        sv = min(norm(r) for r in rows)
        return ok, sv, pv, rows, sigs

    if method == "P2":
        d_cand, diff = recover_d_cvp(sigs, n, lam, k1_bound, k2_bound, d_secret)
        pvv = planted_projected(sigs, n, k1_bound, k2_bound)
        pv = norm(pvv[:2 * m])          # same vector minus the Kannan slot
        sv = norm(diff) if diff else float('nan')
        return d_cand is not None, sv, pv, None, sigs

    raise ValueError(method)


SEEDS = [42, 1234, 9999, 555, 31337]


def rate(curve, m, k1_bound, method, seeds=SEEDS):
    n = curve[2]
    wins, ratios = 0, []
    for seed in seeds:
        d_trial = random.Random(seed + 7777).randint(1, n - 1)
        res = trial(curve, m, d_trial, k1_bound, seed, method)
        if res is None:
            continue
        ok, sv, pv, _rows, _sigs = res
        wins += bool(ok)
        if pv:
            ratios.append(sv / pv)
    mr = sum(ratios) / len(ratios) if ratios else float('nan')
    return wins, len(seeds), mr


def gh_log(det_log, dim):
    """log(Gaussian heuristic) = log sqrt(dim/(2 pi e)) + det_log/dim."""
    return 0.5 * math.log(dim / (2 * math.pi * math.e)) + det_log / dim


# ===========================================================================
print("=" * 78)
print("Thread 23 — project out the trivial direction; does the K1 wall move?")
print("=" * 78)

HIST = [
    # label,             p,    b, n,    lam,  K1, m      (2026-07-29 EXP T4/T5)
    ("8-bit/199",        211,  2, 199,  106,  2,  6),
    ("12-bit/2557",      2557, 2, 2659, 1755, 8,  8),
    ("12-bit/2677 FAIL", 2677, 2, 2647, 185,  8,  10),
]
hist = []
for label, p, b, n, lam, k1, m in HIST:
    G = find_generator(p, b, n)
    assert G is not None, f"no generator for {label}"
    assert (lam * lam + lam + 1) % n == 0, f"bad lam for {label}"
    hist.append((label, (p, b, n, lam, G), k1, m))


# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("EXP E1: is pi(v_planted) the shortest vector of pi(L)?")
print("-" * 78)
print("T5 (2026-07-29) measured sv/pv in [0.34, 0.61] in the UNPROJECTED")
print("lattice, with the shortest vector 100% in the d-column.  Same curves,")
print("same seeds, same K1 -- now in pi(L).  sv/pv >= 1 means the planted")
print("vector is (tied for) lambda_1.\n")
print(f"{'curve':<18} {'K1':>3} {'m':>3} "
      f"{'ORIG sv/pv':>11} {'P1 sv/pv':>10} {'P1 planted=lam1?':>18} "
      f"{'||pi(v)||/GH':>13}")
for label, curve, k1, m in hist:
    p, b, n, lam, G = curve
    k2b = math.isqrt(n) + 1
    d0 = random.Random(42 + 7777).randint(1, n - 1)
    sigs = gen_signatures(G, d0, m, n, lam, p, k1, k2b, 42)
    S_K1, S_D, S_K2, S_KAN = scales(n, k1, k2b)

    M = build_glv_lattice(sigs, n, lam, k1, k2b)
    rows_o = lll_rows(M, 2 * m + 2)
    sv_o = min(norm(r) for r in rows_o)
    pv_o = norm(planted_vector(sigs, d0, n, k1, k2b))

    rows_p = lll_rows(project_out_d(M, m), 2 * m + 1)
    pvv = planted_projected(sigs, n, k1, k2b)
    pv_p = norm(pvv)
    sv_p = min(norm(r) for r in rows_p)

    det_log = m * math.log(n) + m * math.log(S_K1) + m * math.log(S_K2)
    ratio_gh = math.exp(math.log(pv_p) - gh_log(det_log, 2 * m + 1))

    print(f"{label:<18} {k1:>3} {m:>3} {sv_o/pv_o:>11.3f} {sv_p/pv_p:>10.3f} "
          f"{str(sv_p/pv_p >= 0.999):>18} {ratio_gh:>13.3f}")

print("\nEnergy split of the shortest vector of pi(L) "
      "(k1-blk = cols 0..m-1, k2-blk = m..2m-1, kan = col 2m):")
print(f"{'curve':<18} {'k1-blk':>8} {'k2-blk':>8} {'kan':>8} {'kan=S?':>7}")
for label, curve, k1, m in hist:
    p, b, n, lam, G = curve
    k2b = math.isqrt(n) + 1
    d0 = random.Random(42 + 7777).randint(1, n - 1)
    sigs = gen_signatures(G, d0, m, n, lam, p, k1, k2b, 42)
    _S_K1, _S_D, _S_K2, S_KAN = scales(n, k1, k2b)
    M = build_glv_lattice(sigs, n, lam, k1, k2b)
    rows_p = lll_rows(project_out_d(M, m), 2 * m + 1)
    sv = min(rows_p, key=norm)
    tot = sum(x * x for x in sv) or 1
    e_k1 = sum(sv[i] ** 2 for i in range(m)) / tot
    e_k2 = sum(sv[m + i] ** 2 for i in range(m)) / tot
    e_kn = sv[2 * m] ** 2 / tot
    print(f"{label:<18} {e_k1:>8.3f} {e_k2:>8.3f} {e_kn:>8.3f} "
          f"{str(abs(sv[2*m]) == S_KAN):>7}")


# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("EXP E2 (THE FALSIFIER): T4's K1 grid, all three methods, m=12, 5 seeds")
print("-" * 78)
print("T4 (2026-07-29) found the wall at K1 ~ 12-16 for lam*=0.34 and")
print("K1 ~ 4-6 for lam*=0.07.  If P1/P2 push those outward, the")
print("reformulation is a real improvement; if not, the wall is")
print("information-theoretic.\n")

K1_GRID = [2, 3, 4, 6, 8, 12, 16, 24]
M_SIGS = 12
print(f"{'curve':<18} {'lam*':>7} {'method':<7} "
      + " ".join(f"K1={k:<4}" for k in K1_GRID))
e2 = {}
for label, curve, _k1, _m in hist:
    p, b, n, lam, G = curve
    if n < 2000:
        continue                      # 8-bit curve: K2=15, grid not comparable
    for method in ("ORIG", "P1", "P2"):
        cells = []
        for k1 in K1_GRID:
            w, t, _ = rate(curve, M_SIGS, k1, method)
            e2[(label, method, k1)] = w
            cells.append(f"{w}/{t} ")
        print(f"{label:<18} {lam_star(lam,n):>7.4f} {method:<7} "
              + " ".join(cells))
    print()

print("wall(K1) = largest K1 with >=3/5 recovery:")
for label, curve, _k1, _m in hist:
    if curve[2] < 2000:
        continue
    for method in ("ORIG", "P1", "P2"):
        ks = [k for k in K1_GRID if e2.get((label, method, k), 0) >= 3]
        print(f"  {label:<18} {method:<5} wall = {max(ks) if ks else '<2'}")


# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("EXP E2b: does more data rescue P1/P2 at the wall?  (T4b replication)")
print("-" * 78)
fail_curve = [c for lbl, c, _, _ in hist if c[2] == 2647][0]
for method in ("ORIG", "P1", "P2"):
    cells = []
    for m_try in (8, 12, 16, 24, 32):
        w, t, _ = rate(fail_curve, m_try, 8, method)
        cells.append(f"m={m_try}: {w}/{t}")
    print(f"  K1=8  {method:<5} " + "   ".join(cells))


# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("EXP E3: does ||pi(v)|| / GH(pi(L)) predict the eff wall?")
print("-" * 78)
print("Asymptotically ||pi(v)||/GH -> sqrt(2 pi e * eff / 3) = 2.266*sqrt(eff),")
print("so the heuristic predicts success iff eff < 3/(2 pi e) = 0.1947.")
print("T3 (2026-07-29) measured: eff=0.05 works, 0.15 marginal, 0.25 fails.\n")
print(f"{'eff':>6} {'pred ||pi(v)||/GH':>18} {'heuristic':>12}")
for eff in (0.02, 0.05, 0.10, 0.15, 0.1947, 0.25, 0.40):
    r = math.sqrt(2 * math.pi * math.e * eff / 3.0)
    print(f"{eff:>6.4f} {r:>18.3f} {'findable' if r < 1 else 'buried':>12}")

print("\nMeasured on the historical curves (exact, not asymptotic):")
print(f"{'curve':<18} {'K1':>3} {'m':>3} {'eff':>7} {'||pi(v)||/GH':>13} "
      f"{'P1 5-seed':>10}")
for label, curve, k1, m in hist:
    p, b, n, lam, G = curve
    k2b = math.isqrt(n) + 1
    for k1v in (2, 8):
        eff = k1v * k2b / n
        S_K1, _S_D, S_K2, _S_KAN = scales(n, k1v, k2b)
        d0 = random.Random(42 + 7777).randint(1, n - 1)
        sigs = gen_signatures(G, d0, M_SIGS, n, lam, p, k1v, k2b, 42)
        if len(sigs) < M_SIGS:
            continue
        pv = norm(planted_projected(sigs, n, k1v, k2b))
        det_log = (M_SIGS * math.log(n) + M_SIGS * math.log(S_K1)
                   + M_SIGS * math.log(S_K2))
        r = math.exp(math.log(pv) - gh_log(det_log, 2 * M_SIGS + 1))
        w, t, _ = rate(curve, M_SIGS, k1v, "P1")
        print(f"{label:<18} {k1v:>3} {M_SIGS:>3} {eff:>7.4f} {r:>13.3f} "
              f"{w}/{t:>9}")


# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("EXP E4: fresh 17-bit curves, eff in {0.05, 0.15, 0.25}, ORIG vs P1 vs P2")
print("-" * 78)

def search_curves_17bit(lo, hi, want=8):
    """j=0 GLV curves with n prime, n = 1 mod 3, spread over lam*."""
    out = []
    pr = sympy.nextprime(lo)
    while pr < hi and len(out) < want:
        if pr % 3 == 1:
            eis = eisenstein_decompose(pr)
            if eis:
                a_e, b_e = eis
                for t in j0_traces(a_e, b_e):
                    nc = pr + 1 - t
                    if nc < 2 or not sympy.isprime(nc) or nc % 3 != 1:
                        continue
                    roots = glv_roots(nc)
                    if roots is None:
                        continue
                    lam = roots[0]
                    b = None
                    rng = random.Random(7 * pr)
                    for _ in range(300):
                        bt = rng.randint(1, pr - 1)
                        x = rng.randint(0, pr - 1)
                        rhs = (pow(x, 3, pr) + bt) % pr
                        y = tonelli_shanks(rhs, pr)
                        if y is None or y == 0:
                            continue
                        if ec_mul((x, y), nc, pr) is None:
                            b = bt
                            break
                    if b is None:
                        continue
                    G = find_generator(pr, b, nc, seed=pr)
                    if G is None:
                        continue
                    out.append((pr, b, nc, lam, G))
                    break
        pr = sympy.nextprime(pr)
    return out

curves17 = search_curves_17bit(2 ** 16, 2 ** 17, want=6)
print(f"harvested {len(curves17)} 17-bit curves\n")
E4_SEEDS = [42, 1234, 9999]
print(f"{'n':>8} {'lam*':>7} {'eff':>6} {'ORIG':>6} {'P1':>6} {'P2':>6}")
for eff_t in (0.05, 0.15, 0.25):
    for cv in curves17:
        p, b, n, lam, G = cv
        k2b = math.isqrt(n) + 1
        k1v = max(2, int(round(eff_t * n / k2b)))
        row = []
        for method in ("ORIG", "P1", "P2"):
            w, t, _ = rate(cv, 10, k1v, method, seeds=E4_SEEDS)
            row.append(f"{w}/{t}")
        print(f"{n:>8} {lam_star(lam,n):>7.4f} {k1v*k2b/n:>6.3f} "
              f"{row[0]:>6} {row[1]:>6} {row[2]:>6}")
    print()


# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("EXP E5: per-seed agreement -- is P1 the SAME attack as ORIG?")
print("-" * 78)
print("E2 shows identical 5-seed counts.  That could still be two different")
print("attacks with equal power.  Here we compare outcomes trial by trial:")
print("if ORIG and P1 succeed on exactly the same (curve, K1, seed) triples,")
print("the projection is a lattice reformulation with no effect on LLL at all.\n")

agree = {"ORIG-P1": [0, 0], "ORIG-P2": [0, 0]}
disagree_examples = []
for label, curve, _k1, _m in hist:
    p, b, n, lam, G = curve
    for k1 in K1_GRID:
        for seed in SEEDS:
            d_trial = random.Random(seed + 7777).randint(1, n - 1)
            outs = {}
            for method in ("ORIG", "P1", "P2"):
                res = trial(curve, M_SIGS, d_trial, k1, seed, method)
                outs[method] = bool(res[0]) if res else None
            for pair in ("ORIG-P1", "ORIG-P2"):
                a, c = pair.split("-")
                if outs[a] is None or outs[c] is None:
                    continue
                agree[pair][1] += 1
                if outs[a] == outs[c]:
                    agree[pair][0] += 1
                elif pair == "ORIG-P1" and len(disagree_examples) < 8:
                    disagree_examples.append(
                        (label, k1, seed, outs["ORIG"], outs["P1"]))

for pair, (ok, tot) in agree.items():
    print(f"  {pair}: {ok}/{tot} trials agree ({100.0*ok/tot:.1f}%)")
if disagree_examples:
    print("\n  ORIG/P1 disagreements (label, K1, seed, ORIG, P1):")
    for row in disagree_examples:
        print(f"    {row}")
else:
    print("\n  ORIG and P1 agree on EVERY trial.")


# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("EXP E6: out-of-sample check of nu_hat (2026-07-29, commit e845207)")
print("-" * 78)
print("E4 at eff~0.149 splits the 17-bit curves into successes and failures")
print("that lam* does NOT explain (lam*=0.335 and 0.318 succeed; 0.388,")
print("0.358, 0.410, 0.027 fail).  nu_hat = lambda_1(L2)/sqrt(det L2) with")
print("L2 = <(n*S_K1,0), (-lam*S_K1,S_K2)> claims low nu_hat => easier.")
print("Prediction: the successes sit at LOWER nu_hat than the failures.\n")

def gauss_reduce_2d(u, v):
    """Lagrange/Gauss reduction; returns the exact shortest nonzero vector.
    Verbatim from glv_hnp_phase2_lambda_threshold.py:188."""
    def nrm2(w):
        return w[0] * w[0] + w[1] * w[1]
    def dot(w, z):
        return w[0] * z[0] + w[1] * z[1]
    if nrm2(u) > nrm2(v):
        u, v = v, u
    while True:
        num, den = dot(v, u), nrm2(u)
        if den == 0:
            break
        q = (2*num + den) // (2*den) if num >= 0 else -((-2*num + den) // (2*den))
        v = (v[0] - q * u[0], v[1] - q * u[1])
        if nrm2(v) >= nrm2(u):
            break
        u, v = v, u
    return u

def nu_hat(n, lam, S_K1, S_K2):
    w = gauss_reduce_2d((n * S_K1, 0), (-(lam % n) * S_K1, S_K2))
    return math.sqrt(w[0]**2 + w[1]**2) / math.sqrt(n * S_K1 * S_K2)

curves17b = search_curves_17bit(2 ** 16, 2 ** 17, want=14)
EFF_BAND = 0.149
print(f"{'n':>8} {'lam*':>7} {'nu_hat':>8} {'eff':>6} {'ORIG 5-seed':>12}")
rows6 = []
for cv in curves17b:
    p, b, n, lam, G = cv
    k2b = math.isqrt(n) + 1
    k1v = max(2, int(round(EFF_BAND * n / k2b)))
    S_K1, _S_D, S_K2, _S_KAN = scales(n, k1v, k2b)
    nh = nu_hat(n, lam, S_K1, S_K2)
    w, t, _ = rate(cv, 10, k1v, "ORIG")
    rows6.append((n, lam_star(lam, n), nh, k1v * k2b / n, w, t))
    print(f"{n:>8} {lam_star(lam,n):>7.4f} {nh:>8.4f} {k1v*k2b/n:>6.3f} "
          f"{w}/{t:>10}")

succ = [r for r in rows6 if r[4] >= 3]
fail = [r for r in rows6 if r[4] < 3]
def mean(xs): return sum(xs) / len(xs) if xs else float('nan')
print(f"\n  successes (>=3/5): {len(succ)}  mean nu_hat = {mean([r[2] for r in succ]):.4f}"
      f"  mean lam* = {mean([r[1] for r in succ]):.4f}")
print(f"  failures   (<3/5): {len(fail)}  mean nu_hat = {mean([r[2] for r in fail]):.4f}"
      f"  mean lam* = {mean([r[1] for r in fail]):.4f}")
if succ and fail:
    for key, idx in (("nu_hat", 2), ("lam*", 1)):
        pairs = sum(1 for a in succ for c in fail if a[idx] < c[idx])
        tot = len(succ) * len(fail)
        print(f"  AUC({key}, low=success) = {pairs/tot:.3f}   "
              f"({pairs}/{tot} concordant pairs)")

print("\nDone.")
