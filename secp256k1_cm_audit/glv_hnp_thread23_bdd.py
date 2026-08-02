"""
GLV-HNP Phase 2, Thread 23: is the K1 wall a reduction failure or a BDD wall?

Background (RESEARCH_AUTOLAB_LOG.md, 2026-07-29 entry, EXP T5):
  The Phase-2 Kannan lattice (dim 2m+2) always has lambda_1 = the trivial
  vector n*S_D*e_m, never the planted vector (|sv[m]|/n = 1.0000 on every
  curve tested, sv/pv in [0.34, 0.61]).  That entry concluded recovery is a
  BDD/coset condition rather than an SVP condition, and proposed Thread 23:

    "project the lattice along e_m ... or replace the Kannan embedding with
     an explicit CVP call targeting (A_i*S_K1, 0, ...).  Falsifier: if sv/pv
     rises above 1 after the reformulation AND the K1 wall in T4 moves
     outward on the lam*=0.07 curve (currently K1 ~ 4-6), the reformulation
     is a real improvement; if the wall stays at K1 ~ 4-6, the wall is
     information-theoretic and Phase 2 is at its ceiling."

  That entry also asserted, parenthetically, that "no choice of S_D removes
  [the trivial vector] -- both vectors scale linearly in S_D".  The REASON is
  wrong (only the d-coordinate of the planted vector scales with S_D, so
  S_D > sqrt(m + 3/2) does beat the trivial vector) but the CONCLUSION holds
  for a different reason, measured in EXP B: raising S_D just hands lambda_1
  to a different rival, the 2-D block L2 = <(n*S_K1,0), (-lam*S_K1,S_K2)>,
  whose length is independent of S_D.  sv/pv peaks at 0.855 (S_D = 3) and
  falls again; the planted vector is never lambda_1 at any S_D.

RESULTS (2026-08-02 run, output in glv_hnp_thread23_bdd_output.txt)
  A   tau < 1 is sufficient but not necessary: 6/6 configurations with tau < 1
      recover, and successes reach tau = 1.362.  Failures start at tau = 1.161.
  B   sv/pv never crosses 1 (max 0.855); success degrades monotonically once
      S_D > 4 as tau grows.  S_D = 1 is optimal.
  C   the K1 wall MOVES on the hard curve under exact CVP (4 -> 6) and not on
      the easy one (8 -> 8).  BKZ(20) on the Kannan lattice does not move it.
  D   on 6 fresh 17-bit curves the aggregate wall sits at eff ~ 0.10 (tau ~ 1.0
      to 1.3), consistent with tau = 2.386*sqrt(eff).

This script runs the falsifier.

  L_kan  (dim 2m+2)  the historical Kannan lattice, S_D parameterised
  L_cvp  (dim 2m+1)  the same lattice with the Kannan row/column deleted;
                     recovery becomes CVP with target t_i = -A_i*S_K1.

  tau  = ||e|| / GH(L_cvp),  GH(L) = sqrt(dim/(2*pi*e)) * det(L)^(1/dim),
         e = (k1_i*S_K1, d*S_D, k2_i*S_K2)  the true CVP error vector.
  tau < 1 is the unique-decoding (BDD) regime; tau > 1 means the target is
  not the closest lattice point and NO reduction algorithm can recover d.

Experiments
  A   tau vs the measured T4 K1-wall on the two historical 12-bit curves.
      Closed form: det(L_cvp) = (n*S_K1)^m * S_D * S_K2^m, so with
      S_K1 = n/K1, S_K2 = n/K2, S_D = 1 and eff = K1*K2/n,
          tau  ->  sqrt(2*pi*e/3) * sqrt(eff)  =  2.386 * sqrt(eff)
      as m -> infinity.  tau = 1 at eff = 0.1756.
  B   S_D sweep {1,2,3,4,6,8,16,32} on L_kan: does sv/pv cross 1, and does
      success follow?
  C   L_cvp with Babai nearest-plane and with exact CVP enumeration, vs
      L_kan, on the full T4 K1 grid.  Does the wall move?
  D   eff sweep on fresh curves: locate the empirical wall and compare to
      the tau = 1 prediction.

Run: python3 glv_hnp_thread23_bdd.py
"""

import math
import random
import time

import sympy
from fpylll import CVP, GSO, LLL, BKZ, IntegerMatrix

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
        i, t2 = 0, t
        while t2 != 1:
            t2 = t2 * t2 % p; i += 1
        b = pow(c, 1 << (mm - i - 1), p)
        mm, c, t, r = i, b * b % p, t * b * b % p, r * b % p

def find_generator(p, b, n, seed=12345):
    rng = random.Random(seed)
    for _ in range(4000):
        x = rng.randint(0, p - 1)
        y = tonelli_shanks((pow(x, 3, p) + b) % p, p)
        if y is None or y == 0: continue
        P = (x, y)
        if ec_mul(P, n, p) is None and P is not None:
            return P
    return None

def eisenstein_decompose(p):
    """p = a^2 - a*b + b^2 (p = 1 mod 3).  Returns (a, b) or None."""
    if p % 3 != 1: return None
    for a in range(1, int(math.isqrt(4 * p // 3)) + 2):
        disc = 4 * p - 3 * a * a
        if disc < 0: break
        r = math.isqrt(disc)
        if r * r == disc and (a + r) % 2 == 0:
            return (a, (a + r) // 2)
    return None

def j0_traces(a, b):
    """The six traces of the sextic-twist family of a j=0 curve."""
    t = 2 * a - b
    return [t, -t, a + b, -(a + b), a - 2 * b, -(a - 2 * b)]

def glv_roots(n):
    """Both roots of x^2 + x + 1 = 0 mod n, or None."""
    if n % 3 != 1: return None
    for _ in range(200):
        g = random.Random(n).randint(2, n - 2) if False else None
        break
    # deterministic: use a small search for a cube root of unity
    for g in range(2, 200):
        r = pow(g, (n - 1) // 3, n)
        if r != 1 and (r * r + r + 1) % n == 0:
            return (min(r, n - 1 - r), max(r, n - 1 - r))
    return None

def lam_star(lam, n):
    return min(lam % n, n - (lam % n)) / n

# ---------------------------------------------------------------------------
# Scales, signatures (verbatim, S_D now a parameter)
# ---------------------------------------------------------------------------

def scales(n, k1_bound, k2_bound, S_D=1):
    return (n // k1_bound, S_D, max(1, n // k2_bound), n)

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

# ---------------------------------------------------------------------------
# The two lattices
# ---------------------------------------------------------------------------

def build_kannan(sigs, n, lam, k1_bound, k2_bound, S_D=1):
    """Historical Phase-2 Kannan lattice, dim 2m+2."""
    m = len(sigs)
    dim = 2 * m + 2
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound, S_D)
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

def build_cvp(sigs, n, lam, k1_bound, k2_bound, S_D=1):
    """Same lattice with the Kannan row and column deleted, dim 2m+1,
    plus the CVP target t_i = -A_i*S_K1 (0 elsewhere)."""
    m = len(sigs)
    dim = 2 * m + 1
    S_K1, S_D, S_K2, _ = scales(n, k1_bound, k2_bound, S_D)
    M = [[0] * dim for _ in range(dim)]
    for i in range(m):
        M[i][i] = n * S_K1
    for i in range(m):
        M[m][i] = sigs[i]['B'] * S_K1
    M[m][m] = S_D
    for i in range(m):
        M[m + 1 + i][i] = -lam * S_K1
        M[m + 1 + i][m + 1 + i] = S_K2
    t = [0] * dim
    for i in range(m):
        t[i] = -sigs[i]['A'] * S_K1
    return M, t

def planted_kannan(sigs, d, n, k1_bound, k2_bound, S_D=1):
    m = len(sigs)
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound, S_D)
    v = [0] * (2 * m + 2)
    for i in range(m):
        v[i] = sigs[i]['k1'] * S_K1
        v[m + 1 + i] = sigs[i]['k2'] * S_K2
    v[m] = d * S_D
    v[2 * m + 1] = S_KAN
    return v

def error_vector(sigs, d, n, k1_bound, k2_bound, S_D=1):
    """The true CVP error e = v_closest - t, dim 2m+1."""
    m = len(sigs)
    S_K1, S_D, S_K2, _ = scales(n, k1_bound, k2_bound, S_D)
    e = [0] * (2 * m + 1)
    for i in range(m):
        e[i] = sigs[i]['k1'] * S_K1
        e[m + 1 + i] = sigs[i]['k2'] * S_K2
    e[m] = d * S_D
    return e

def norm(v):
    return math.sqrt(sum(float(x) * float(x) for x in v))

# ---------------------------------------------------------------------------
# tau: BDD radius of the CVP formulation
# ---------------------------------------------------------------------------

def gh_cvp(m, n, k1_bound, k2_bound, S_D=1):
    """Gaussian-heuristic lambda_1 of L_cvp.  det = (n*S_K1)^m * S_D * S_K2^m."""
    S_K1, S_D, S_K2, _ = scales(n, k1_bound, k2_bound, S_D)
    dim = 2 * m + 1
    log_det = m * math.log(n * S_K1) + math.log(S_D) + m * math.log(S_K2)
    return math.sqrt(dim / (2 * math.pi * math.e)) * math.exp(log_det / dim)

def tau_expected(m, n, k1_bound, k2_bound, S_D=1):
    """E[||e||] / GH, with k1~U[0,K1), k2~U[0,K2), d~U[0,n)."""
    S_K1, S_D, S_K2, _ = scales(n, k1_bound, k2_bound, S_D)
    e_norm = math.sqrt(m * (k1_bound * S_K1) ** 2 / 3.0
                       + (n * S_D) ** 2 / 3.0
                       + m * (k2_bound * S_K2) ** 2 / 3.0)
    return e_norm / gh_cvp(m, n, k1_bound, k2_bound, S_D)

TAU_ASYMPTOTIC = math.sqrt(2 * math.pi * math.e / 3.0)   # 2.3861...

# ---------------------------------------------------------------------------
# Recovery routines
# ---------------------------------------------------------------------------

def recover_kannan(reduced, m, n, S_KAN, S_D, d_secret):
    dim = 2 * m + 2
    for row in reduced:
        last = row[dim - 1]
        if abs(last) != S_KAN: continue
        sign = 1 if last > 0 else -1
        val = sign * row[m]
        if val % S_D != 0: continue
        d_cand = (val // S_D) % n
        if d_cand and d_cand == d_secret:
            return True
    return False

def run_kannan(sigs, n, lam, k1_bound, k2_bound, d, S_D=1, beta=0):
    m = len(sigs)
    dim = 2 * m + 2
    M = build_kannan(sigs, n, lam, k1_bound, k2_bound, S_D)
    A = IntegerMatrix.from_matrix(M)
    if beta:
        BKZ.reduction(A, BKZ.Param(beta))
    else:
        LLL.reduction(A)
    reduced = [[A[i][j] for j in range(dim)] for i in range(dim)]
    _, S_Dv, _, S_KAN = scales(n, k1_bound, k2_bound, S_D)
    ok = recover_kannan(reduced, m, n, S_KAN, S_Dv, d)
    pv = norm(planted_kannan(sigs, d, n, k1_bound, k2_bound, S_D))
    sv = min(norm(r) for r in reduced)
    return ok, pv, sv

def run_cvp(sigs, n, lam, k1_bound, k2_bound, d, S_D=1, exact=False, beta=0):
    """Babai nearest-plane (exact=False) or CVP enumeration (exact=True)."""
    m = len(sigs)
    dim = 2 * m + 1
    M, t = build_cvp(sigs, n, lam, k1_bound, k2_bound, S_D)
    A = IntegerMatrix.from_matrix(M)
    LLL.reduction(A)
    if beta:
        BKZ.reduction(A, BKZ.Param(beta))
    _, S_Dv, _, _ = scales(n, k1_bound, k2_bound, S_D)
    if exact:
        v = CVP.closest_vector(A, tuple(t))
    else:
        G = GSO.Mat(A, float_type="d")
        G.update_gso()
        coeffs = G.babai(list(t))
        v = [sum(coeffs[i] * A[i][j] for i in range(dim)) for j in range(dim)]
    e = [v[j] - t[j] for j in range(dim)]
    if e[m] % S_Dv != 0:
        return False, e
    return (e[m] // S_Dv) % n == d % n, e

# ---------------------------------------------------------------------------
# Curves
# ---------------------------------------------------------------------------

HIST = [
    # label,             p,    b, n,    lam,  K1_hist, m
    ("8-bit/199",        211,  2, 199,  106,  2,  6),
    ("12-bit/2557",      2557, 2, 2659, 1755, 8,  8),
    ("12-bit/2677 FAIL", 2677, 2, 2647, 185,  8,  10),
]

hist = []
for label, p, b, n, lam, k1h, m in HIST:
    G = find_generator(p, b, n)
    assert G is not None, label
    assert (lam * lam + lam + 1) % n == 0, label
    hist.append((label, (p, b, n, lam, G), k1h, m))

SEEDS = [42, 1234, 9999, 555, 31337]

def trial_d(seed, n):
    return random.Random(seed + 7777).randint(1, n - 1)

def rate(curve, m, k1_bound, seeds, mode="kannan", S_D=1, exact=False, beta=0):
    """Returns (wins, total, mean sv/pv or mean ||e||/GH)."""
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    wins, ratios, taus = 0, [], []
    for seed in seeds:
        d = trial_d(seed, n)
        sigs = gen_signatures(G, d, m, n, lam, p, k1_bound, k2_bound, seed)
        if len(sigs) < m:
            continue
        taus.append(norm(error_vector(sigs, d, n, k1_bound, k2_bound, S_D))
                    / gh_cvp(m, n, k1_bound, k2_bound, S_D))
        if mode == "kannan":
            ok, pv, sv = run_kannan(sigs, n, lam, k1_bound, k2_bound, d, S_D, beta)
            ratios.append(sv / pv)
        else:
            ok, _ = run_cvp(sigs, n, lam, k1_bound, k2_bound, d, S_D, exact, beta)
        wins += bool(ok)
    mr = sum(ratios) / len(ratios) if ratios else float('nan')
    mt = sum(taus) / len(taus) if taus else float('nan')
    return wins, len(seeds), mr, mt


def find_curves(lo, hi, want=6):
    out = []
    p = int(sympy.nextprime(lo))
    while p < hi and len(out) < want:
        if p % 3 == 1:
            eis = eisenstein_decompose(p)
            if eis is not None:
                for t in j0_traces(*eis):
                    nc = p + 1 - t
                    if nc < 2 or nc % 3 != 1 or not sympy.isprime(nc):
                        continue
                    roots = glv_roots(nc)
                    if roots is None:
                        continue
                    rng = random.Random(nc)
                    for _ in range(200):
                        b = rng.randint(1, p - 1)
                        x = rng.randint(0, p - 1)
                        y = tonelli_shanks((pow(x, 3, p) + b) % p, p)
                        if y is None or y == 0:
                            continue
                        if ec_mul((x, y), nc, p) is None:
                            G = find_generator(p, b, nc)
                            if G is not None:
                                out.append((p, b, nc, roots[0], G))
                            break
                    break
        p = int(sympy.nextprime(p))
    return out


def main():
    print("=" * 78)
    print("Thread 23 — BDD reformulation of the GLV-HNP Phase-2 lattice")
    print("=" * 78)

    # ===========================================================================
    # EXP A — tau vs the measured T4 K1 wall
    # ===========================================================================
    print("\n" + "-" * 78)
    print("EXP A: tau = ||e|| / GH(L_cvp) against the measured T4 K1 wall")
    print("-" * 78)
    print(f"asymptotic tau -> {TAU_ASYMPTOTIC:.4f} * sqrt(eff);  tau = 1 at "
          f"eff = {1.0 / TAU_ASYMPTOTIC**2:.4f}")

    K1_GRID = [2, 3, 4, 6, 8, 12, 16, 24]
    T4_CURVES = [c for c in hist if c[0].startswith("12-bit")]

    print(f"\n{'curve':<20}{'K1':>4}{'eff':>8}{'tau_e':>8}{'tau_pred':>9}"
          f"{'LLL':>7}{'sv/pv':>8}")
    expA = []
    for label, curve, _k1h, m in T4_CURVES:
        p, b, n, lam, G = curve
        k2_bound = math.isqrt(n) + 1
        for k1 in K1_GRID:
            eff = k1 * k2_bound / n
            w, tot, mr, mt = rate(curve, m, k1, SEEDS, mode="kannan")
            expA.append((label, m, k1, eff, mt, w, tot, mr))
            print(f"{label:<20}{k1:>4}{eff:>8.3f}{mt:>8.3f}"
                  f"{TAU_ASYMPTOTIC*math.sqrt(eff):>9.3f}{f'{w}/{tot}':>7}{mr:>8.3f}")

    print("\nConfusion of the rule 'predict success iff tau < 1':")
    tp = sum(1 for r in expA if r[4] < 1 and r[5] >= 3)
    fp = sum(1 for r in expA if r[4] < 1 and r[5] < 3)
    fn = sum(1 for r in expA if r[4] >= 1 and r[5] >= 3)
    tn = sum(1 for r in expA if r[4] >= 1 and r[5] < 3)
    print(f"  tau<1 & success {tp}   tau<1 & fail {fp}   "
          f"tau>=1 & success {fn}   tau>=1 & fail {tn}")
    succ = [r[4] for r in expA if r[5] >= 3]
    fail = [r[4] for r in expA if r[5] < 3]
    if succ and fail:
        print(f"  max tau among successes = {max(succ):.3f}   "
              f"min tau among failures = {min(fail):.3f}")

    # ===========================================================================
    # EXP B — S_D sweep: does making the planted vector lambda_1 help?
    # ===========================================================================
    print("\n" + "-" * 78)
    print("EXP B: S_D sweep on L_kan.  Trivial rival = n*S_D*e_m (norm n*S_D);")
    print("       ||v_planted||^2 = n^2(2m/3 + 1) + (n*S_D)^2/3, so the planted")
    print("       vector becomes shorter once S_D > sqrt(m + 3/2).")
    print("-" * 78)

    SD_GRID = [1, 2, 3, 4, 6, 8, 16, 32]
    for label, curve, _k1h, m in T4_CURVES:
        p, b, n, lam, G = curve
        k1 = 4                      # inside the wall for both curves (T4)
        print(f"\n{label}  m={m}  K1={k1}  (crossover predicted at "
              f"S_D > {math.sqrt(m + 1.5):.2f})")
        print(f"{'S_D':>5}{'sv/pv':>9}{'tau_e':>9}{'LLL':>8}")
        for S_D in SD_GRID:
            w, tot, mr, mt = rate(curve, m, k1, SEEDS, mode="kannan", S_D=S_D)
            print(f"{S_D:>5}{mr:>9.3f}{mt:>9.3f}{f'{w}/{tot}':>8}")

    # ===========================================================================
    # EXP C — L_cvp (Babai / exact CVP) vs L_kan on the full K1 grid
    # ===========================================================================
    print("\n" + "-" * 78)
    print("EXP C: does the CVP reformulation move the K1 wall?")
    print("-" * 78)
    print(f"{'curve':<20}{'K1':>4}{'tau_e':>8}{'kannan':>9}{'babai':>8}"
          f"{'cvp_ex':>8}{'kan+BKZ20':>11}")
    expC = []
    for label, curve, _k1h, m in T4_CURVES:
        for k1 in K1_GRID:
            wk, _, _, mt = rate(curve, m, k1, SEEDS, mode="kannan")
            wb, _, _, _ = rate(curve, m, k1, SEEDS, mode="cvp", exact=False)
            t0 = time.time()
            try:
                wx, _, _, _ = rate(curve, m, k1, SEEDS, mode="cvp", exact=True)
                wx_s = f"{wx}/5"
            except Exception as exc:                       # enumeration blow-up
                wx_s = f"ERR({type(exc).__name__})"
            el = time.time() - t0
            wz, _, _, _ = rate(curve, m, k1, SEEDS, mode="kannan", beta=20)
            expC.append((label, k1, mt, wk, wb, wx_s, wz))
            print(f"{label:<20}{k1:>4}{mt:>8.3f}{f'{wk}/5':>9}{f'{wb}/5':>8}"
                  f"{wx_s:>8}{f'{wz}/5':>11}   [cvp_ex {el:.1f}s]")

    def wall(rows, idx):
        """Largest K1 with >=3/5 recoveries, per curve."""
        out = {}
        for r in rows:
            lab, k1 = r[0], r[1]
            v = r[idx]
            v = int(v.split('/')[0]) if isinstance(v, str) and '/' in v else v
            if isinstance(v, int) and v >= 3:
                out[lab] = max(out.get(lab, 0), k1)
        return out

    print("\nK1 wall (largest K1 with >= 3/5 recoveries):")
    for name, idx in (("kannan", 3), ("babai", 4), ("cvp_exact", 5), ("kannan+BKZ20", 6)):
        print(f"  {name:<14}{wall(expC, idx)}")

    # ===========================================================================
    # EXP D — eff sweep on fresh curves: where is the empirical wall?
    # ===========================================================================
    print("\n" + "-" * 78)
    print("EXP D: eff sweep on fresh 17-bit curves; empirical wall vs tau = 1")
    print("-" * 78)

    curves_d = find_curves(1 << 16, 1 << 17, want=6)
    print(f"found {len(curves_d)} fresh 17-bit j=0 GLV curves")
    m_d = 12
    print(f"\n{'eff_target':>11}{'K1':>5}{'tau_e':>8}{'kannan':>9}{'babai':>8}")
    for eff_t in [0.03, 0.05, 0.08, 0.12, 0.15, 0.18, 0.22, 0.28]:
        tot_k = tot_b = tot_n = 0
        taus = []
        for cur in curves_d:
            p, b, n, lam, G = cur
            k2 = math.isqrt(n) + 1
            k1 = max(2, round(eff_t * n / k2))
            wk, t_, _, mt = rate(cur, m_d, k1, SEEDS[:3], mode="kannan")
            wb, _, _, _ = rate(cur, m_d, k1, SEEDS[:3], mode="cvp")
            tot_k += wk; tot_b += wb; tot_n += 3
            taus.append(mt)
        mt = sum(taus) / len(taus)
        print(f"{eff_t:>11.2f}{'~':>5}{mt:>8.3f}{f'{tot_k}/{tot_n}':>9}"
              f"{f'{tot_b}/{tot_n}':>8}")

    print("\n" + "=" * 78)
    print("done")
    print("=" * 78)


if __name__ == "__main__":
    main()
