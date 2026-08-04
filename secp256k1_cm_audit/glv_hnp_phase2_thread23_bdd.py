"""
GLV-HNP Phase 2, Thread 23: is the K1 wall algorithmic or information-theoretic?

Background (RESEARCH_AUTOLAB_LOG.md, 2026-07-29 run, exp T5):
  In the Phase-2 Kannan lattice the shortest vector after LLL is ALWAYS the
  trivial vector n*S_D*e_m (100% of its energy in the d-column, Kannan
  coordinate 0).  It carries no information because d is only defined mod n.
  Consequently the planted vector is never lambda_1, and recovery is not an
  SVP condition but a BDD/coset condition: v_planted must be shortest among
  lattice vectors whose Kannan coordinate is +-S_KANNAN.

  The 2026-07-29 entry proposed Thread 23: reformulate so that the target is
  actually the object the algorithm optimises, then re-measure the "K1 wall"
  (T4: both 12-bit historical curves recover for K1 <= 4; the lam*=0.07 curve
  dies at K1 ~ 6-8, the lam*=0.34 curve at K1 ~ 12-16).

Three lattice/algorithm variants are compared on the identical signature sets:

  A  BASELINE   dim 2m+2 Kannan lattice + LLL.       (verbatim 2026-07-26/29)
  B  PROJECTED  d-column dropped + LLL.
       The d-column is the only place the trivial vector lives: the lattice
       satisfies  L cap R*e_m = n*S_D*Z*e_m, so dropping coordinate m is
       exactly the orthogonal projection that quotients the trivial direction
       out.  covol drops by the factor n*S_D = n, dim by 1.  d is recovered
       afterwards from the k1_i,k2_i via  d = (k1_i + lam*k2_i - A_i)/B_i.
  C  EXACT CVP  no Kannan row/column; exact closest-vector (method="proved")
       to the target  -a,  a = (A_i*S_K1, 0, ..., 0),  in the dim-(2m+1)
       lattice L0 = < n*S_K1*e_i ; (B_i*S_K1 | S_D) ; (-lam*S_K1 | S_K2) >.
       By construction  v_planted = a + w  with w in L0, so v_planted is a
       *candidate* solution of this CVP instance and nothing else.

Why C is decisive.  Variants A and B can fail for two different reasons:
the planted vector is not the optimum (information-theoretic), or it is the
optimum but LLL is too weak to find it (algorithmic).  Exact CVP separates
them.  Define the decoding margin

    tau = dist(-a, L0) / ||v_planted||   in (0, 1].

  tau == 1  <=>  v_planted IS the closest vector.  Then any failure of A/B is
                 purely a reduction-strength problem and better lattice
                 reduction is the right next move.
  tau <  1  <=>  some other lattice point is strictly closer.  The instance is
                 then unsolvable by ANY algorithm in this lattice formulation,
                 no matter how strong -- the wall is information-theoretic and
                 Phase 2 is at its ceiling.

Falsifier stated in advance (2026-07-29 log):
  "if sv/pv rises above 1 after the reformulation and the K1 wall moves
   outward on the lam*=0.07 curve (currently K1~4-6), the reformulation is a
   real improvement; if the wall stays at K1~4-6, then the wall is
   information-theoretic and Phase 2 is at its ceiling."

Experiments:
  E1  sanity: all three variants on the historical 8-bit curve.
  E2  K1-wall grid (both historical 12-bit curves x K1 in {2,3,4,6,8,12,16,24}
      x 5 seeds) for variants A and B.  Does the wall move?
  E3  exact-CVP decoding margin tau on the same grid.  Algorithmic or
      information-theoretic?
  E4  sv/pv and energy split for variant B: is the planted vector lambda_1
      once the trivial direction is gone?

Run:  python3 secp256k1_cm_audit/glv_hnp_phase2_thread23_bdd.py
Deps: fpylll, sympy   (pip install fpylll cysignals sympy)
"""

import math
import os
import random
import time

from fpylll import IntegerMatrix, LLL, BKZ, CVP

# ---------------------------------------------------------------------------
# Reuse the Thread-20 helpers VERBATIM so the comparison with the 2026-07-29
# numbers is exact.  That file runs its experiments at import time, so only the
# helper prefix (everything before the first experiment block) is exec'd.
# ---------------------------------------------------------------------------

_HERE = os.path.dirname(os.path.abspath(__file__))
_SRC = os.path.join(_HERE, "glv_hnp_phase2_lambda_threshold.py")
_MARKER = "# EXPERIMENT T1"

with open(_SRC) as fh:
    _text = fh.read()
_prefix = _text[: _text.index(_MARKER)]
exec(compile(_prefix, _SRC, "exec"), globals())

# Names now in scope from the Thread-20 helper prefix:
#   modinv ec_add ec_mul tonelli_shanks find_generator eisenstein_decompose
#   j0_traces glv_roots lam_star gauss_reduce_2d lambda_block_mu
#   planted_norm_expected scales gen_signatures build_glv_lattice
#   planted_vector norm recover_d build_curve search_curves

SEEDS = [42, 1234, 9999, 555, 31337]

# Historical Phase-2 curves (RESEARCH_AUTOLAB_LOG.md 2026-06-15 / 07-26 / 07-29)
HIST = [
    # label,          p,    b, n,    lam,  K1_hist, m
    ("8-bit/199",     211,  2, 199,  106,  2,  6),
    ("12-bit/2557",   2557, 2, 2659, 1755, 8,  12),
    ("12-bit/2677",   2677, 2, 2647, 185,  8,  12),
]

CURVES = {}
for _label, _p, _b, _n, _lam, _k1, _m in HIST:
    _G = find_generator(_p, _b, _n)
    assert _G is not None, f"no generator for {_label}"
    assert (_lam * _lam + _lam + 1) % _n == 0, f"bad lam for {_label}"
    CURVES[_label] = (_p, _b, _n, _lam, _G)


# ===========================================================================
# Variant B — d-column dropped (orthogonal projection along e_m)
# ===========================================================================

def build_lattice_projected(sigs, n, lam, k1_bound, k2_bound):
    """Columns: [k1_0..k1_{m-1} | k2_0..k2_{m-1} | kannan].  dim = 2m+1.

    Same generators as build_glv_lattice with column m deleted.  There are
    2m+2 rows spanning a rank-(2m+1) lattice: n*rowB is an integer combination
    of the n*S_K1*e_i rows once the d-coordinate is gone, which is precisely
    the trivial vector being quotiented out.  LLL handles the dependency and
    emits a leading zero row.
    """
    m = len(sigs)
    ncol = 2 * m + 1
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    M = []
    for i in range(m):                       # n*S_K1*e_i
        row = [0] * ncol
        row[i] = n * S_K1
        M.append(row)
    row = [0] * ncol                         # B-row (d multiplier)
    for i in range(m):
        row[i] = sigs[i]['B'] * S_K1
    M.append(row)
    for j in range(m):                       # lambda rows
        row = [0] * ncol
        for i in range(m):
            if i == j:
                row[i] = -lam * S_K1
        row[m + j] = S_K2
        M.append(row)
    row = [0] * ncol                         # Kannan row
    for i in range(m):
        row[i] = sigs[i]['A'] * S_K1
    row[2 * m] = S_KAN
    M.append(row)
    return M


def planted_vector_projected(sigs, n, k1_bound, k2_bound):
    m = len(sigs)
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    v = [0] * (2 * m + 1)
    for i in range(m):
        v[i] = sigs[i]['k1'] * S_K1
        v[m + i] = sigs[i]['k2'] * S_K2
    v[2 * m] = S_KAN
    return v


def recover_d_projected(rows, sigs, n, lam, k1_bound, k2_bound, d_secret):
    """Read k1_i,k2_i off a row with |kannan| == S_KANNAN, then solve for d.

    d = (k1_i + lam*k2_i - A_i) * B_i^{-1}  (mod n),  cross-checked on a
    second signature so a spurious row cannot produce a false positive.
    """
    m = len(sigs)
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    for row in rows:
        last = row[2 * m]
        if abs(last) != S_KAN:
            continue
        sgn = 1 if last > 0 else -1
        k1s, k2s = [], []
        ok = True
        for i in range(m):
            a, b = sgn * row[i], sgn * row[m + i]
            if a % S_K1 or b % S_K2:
                ok = False
                break
            k1s.append(a // S_K1)
            k2s.append(b // S_K2)
        if not ok:
            continue
        for i in range(m):
            B = sigs[i]['B'] % n
            if B == 0:
                continue
            d_cand = (k1s[i] + lam * k2s[i] - sigs[i]['A']) * modinv(B, n) % n
            if d_cand == 0:
                continue
            j = (i + 1) % m
            lhs = (sigs[j]['A'] + sigs[j]['B'] * d_cand - lam * k2s[j]) % n
            if lhs != k1s[j] % n:
                continue
            if d_cand == d_secret:
                return d_cand
            break
    return None


# ===========================================================================
# Variant C — exact CVP (no Kannan embedding)
# ===========================================================================

def build_lattice_cvp(sigs, n, lam, k1_bound, k2_bound):
    """L0, columns [k1_0..k1_{m-1} | d | k2_0..k2_{m-1}], dim 2m+1, full rank.

    This is the Kannan lattice with the Kannan row AND column removed, i.e.
    the sublattice of vectors with Kannan coordinate 0.
    det(L0) = (n*S_K1)^m * S_D * S_K2^m = n^(2m) / eff^m.
    """
    m = len(sigs)
    ncol = 2 * m + 1
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    M = []
    for i in range(m):
        row = [0] * ncol
        row[i] = n * S_K1
        M.append(row)
    row = [0] * ncol
    for i in range(m):
        row[i] = sigs[i]['B'] * S_K1
    row[m] = S_D
    M.append(row)
    for j in range(m):
        row = [0] * ncol
        row[j] = -lam * S_K1
        row[m + 1 + j] = S_K2
        M.append(row)
    return M


def cvp_target(sigs, n, k1_bound, k2_bound):
    """a = (A_i*S_K1, 0, ..., 0).  v_planted = a + w with w in L0, so the CVP
    target is -a and the recovered vector is a + w."""
    m = len(sigs)
    S_K1, _, _, _ = scales(n, k1_bound, k2_bound)
    a = [0] * (2 * m + 1)
    for i in range(m):
        a[i] = sigs[i]['A'] * S_K1
    return a


def planted_vector_cvp(sigs, d_secret, n, k1_bound, k2_bound, centered=True):
    """The planted vector in L0-coset form.

    IMPORTANT (bug found in the first run of this script, 2026-08-04):
    the whole coset  v_planted + c*n*S_D*e_m,  c in Z,  decodes to the SAME d
    (d is only defined mod n), so ||v_planted|| with the raw representative
    d in [0,n) is NOT the quantity a CVP solver competes against.  The minimal
    representative uses the centered residue d_c = ((d + n/2) mod n) - n/2.
    With the raw representative one measures tau < 1 on instances that are in
    fact decoded perfectly, which is exactly the `rec > opt` inconsistency the
    first run showed.  Default is the centered (correct) representative.
    """
    m = len(sigs)
    S_K1, S_D, S_K2, _ = scales(n, k1_bound, k2_bound)
    v = [0] * (2 * m + 1)
    for i in range(m):
        v[i] = sigs[i]['k1'] * S_K1
        v[m + 1 + i] = sigs[i]['k2'] * S_K2
    d_rep = ((d_secret + n // 2) % n) - n // 2 if centered else d_secret
    v[m] = d_rep * S_D
    return v


def run_cvp(sigs, d_secret, n, lam, k1_bound, k2_bound, method="proved"):
    """Exact CVP.  Returns (d_recovered_or_None, tau, tau_raw, dist, pn)."""
    m = len(sigs)
    M = build_lattice_cvp(sigs, n, lam, k1_bound, k2_bound)
    A = IntegerMatrix.from_matrix(M)
    LLL.reduction(A)
    a = cvp_target(sigs, n, k1_bound, k2_bound)
    t = tuple(-x for x in a)
    w = CVP.closest_vector(A, t, method=method)
    v = [a[i] + w[i] for i in range(2 * m + 1)]
    dist = norm(v)
    pn = norm(planted_vector_cvp(sigs, d_secret, n, k1_bound, k2_bound,
                                 centered=True))
    pn_raw = norm(planted_vector_cvp(sigs, d_secret, n, k1_bound, k2_bound,
                                     centered=False))
    d_cand = v[m] % n                        # S_D == 1
    got = d_cand if d_cand == d_secret else None
    tau = dist / pn if pn else float('nan')
    tau_raw = dist / pn_raw if pn_raw else float('nan')
    return got, tau, tau_raw, dist, pn


# ===========================================================================
# Shared drivers
# ===========================================================================

def make_sigs(curve, m, d_secret, k1_bound, seed):
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed)
    return sigs, k2_bound


def run_variant_A(curve, m, d_secret, k1_bound, seed, sigs, k2_bound,
                  use_bkz=False, beta=20):
    p, b, n, lam, G = curve
    dim = 2 * m + 2
    M = build_glv_lattice(sigs, n, lam, k1_bound, k2_bound)
    A = IntegerMatrix.from_matrix(M)
    if use_bkz:
        BKZ.reduction(A, BKZ.Param(beta))
    else:
        LLL.reduction(A)
    rows = [[A[i][j] for j in range(dim)] for i in range(dim)]
    _, _, _, S_KAN = scales(n, k1_bound, k2_bound)
    ok = recover_d(rows, m, n, S_KAN, d_secret) is not None
    pv = planted_vector(sigs, d_secret, n, k1_bound, k2_bound)
    nz = [r for r in rows if any(r)]
    sv = min(nz, key=norm)
    return ok, norm(sv) / norm(pv), sv, pv


def run_variant_B(curve, m, d_secret, k1_bound, seed, sigs, k2_bound,
                  use_bkz=False, beta=20):
    p, b, n, lam, G = curve
    ncol = 2 * m + 1
    M = build_lattice_projected(sigs, n, lam, k1_bound, k2_bound)
    A = IntegerMatrix.from_matrix(M)
    if use_bkz:
        BKZ.reduction(A, BKZ.Param(beta))
    else:
        LLL.reduction(A)
    rows = [[A[i][j] for j in range(ncol)] for i in range(A.nrows)]
    ok = recover_d_projected(rows, sigs, n, lam, k1_bound, k2_bound,
                             d_secret) is not None
    pv = planted_vector_projected(sigs, n, k1_bound, k2_bound)
    nz = [r for r in rows if any(r)]
    sv = min(nz, key=norm)
    return ok, norm(sv) / norm(pv), sv, pv


def gaussian_heuristic(log_det, dim):
    """GH(L) = sqrt(dim/(2*pi*e)) * det^(1/dim), computed in logs."""
    return math.sqrt(dim / (2 * math.pi * math.e)) * math.exp(log_det / dim)


def banner(title):
    print()
    print("=" * 78)
    print(title)
    print("=" * 78)


# ===========================================================================
# E1 — sanity check on the historical 8-bit curve
# ===========================================================================

banner("Thread 23 — BDD/coset reformulation of the Phase-2 GLV lattice")

print("""
Variants
  A  baseline  dim 2m+2 Kannan lattice, LLL           (2026-07-26/29 verbatim)
  B  projected d-column dropped (trivial vector quotiented out), LLL
  C  exact CVP dim 2m+1, no Kannan embedding, method="proved"

Decoding margin  tau = dist(-a, L0) / ||v_planted||  in (0,1].
  tau == 1  planted vector IS the closest vector -> failures are algorithmic.
  tau <  1  something else is closer -> the instance is unsolvable in this
            lattice formulation by any algorithm (information-theoretic).
""")

print("-" * 78)
print("E1: sanity — all three variants, 8-bit/199, K1=2, m=6")
print("-" * 78)
print(f"{'seed':>7} {'A':>4} {'B':>4} {'C':>4} {'tau':>8} "
      f"{'svA/pv':>8} {'svB/pv':>8}")

curve = CURVES["8-bit/199"]
p, b, n, lam, G = curve
for seed in SEEDS:
    d0 = random.Random(seed + 7777).randint(1, n - 1)
    sigs, k2b = make_sigs(curve, 6, d0, 2, seed)
    if len(sigs) < 6:
        continue
    okA, rA, _, _ = run_variant_A(curve, 6, d0, 2, seed, sigs, k2b)
    okB, rB, _, _ = run_variant_B(curve, 6, d0, 2, seed, sigs, k2b)
    gotC, tau, tau_raw, _, _ = run_cvp(sigs, d0, n, lam, 2, k2b)
    print(f"{seed:>7} {str(bool(okA)):>4} {str(bool(okB)):>4} "
          f"{str(gotC is not None):>4} {tau:>8.4f} {rA:>8.4f} {rB:>8.4f}")


# ===========================================================================
# E2 — the K1 wall, variants A and B
# ===========================================================================

print()
print("-" * 78)
print("E2: K1 wall — does dropping the d-column move it?  m=12, 5 seeds")
print("-" * 78)
print("(2026-07-29 T4 baseline for A: 2557 -> 5,5,5,5,5,4,1,0 ; "
      "2677 -> 5,5,5,2,0,0,0,0)")

K1_GRID = [2, 3, 4, 6, 8, 12, 16, 24]
M_SIGS = 12

wall_rows = {}
print(f"\n{'curve':<14} {'var':>4} " + " ".join(f"K1={k:<3}" for k in K1_GRID))
for label in ("12-bit/2557", "12-bit/2677"):
    curve = CURVES[label]
    p, b, n, lam, G = curve
    for variant, fn in (("A", run_variant_A), ("B", run_variant_B)):
        cells, ratios = [], []
        for k1 in K1_GRID:
            wins = 0
            rr = []
            for seed in SEEDS:
                d0 = random.Random(seed + 7777).randint(1, n - 1)
                sigs, k2b = make_sigs(curve, M_SIGS, d0, k1, seed)
                if len(sigs) < M_SIGS:
                    continue
                ok, r, _, _ = fn(curve, M_SIGS, d0, k1, seed, sigs, k2b)
                wins += bool(ok)
                rr.append(r)
            cells.append(wins)
            ratios.append(sum(rr) / len(rr) if rr else float('nan'))
        wall_rows[(label, variant)] = (cells, ratios)
        print(f"{label:<14} {variant:>4} " +
              " ".join(f"{c}/5{'':<2}" for c in cells))

print(f"\n{'curve':<14} {'var':>4} mean sv/pv (shortest nonzero / planted)")
for label in ("12-bit/2557", "12-bit/2677"):
    for variant in ("A", "B"):
        _, ratios = wall_rows[(label, variant)]
        print(f"{label:<14} {variant:>4} " +
              " ".join(f"{r:>6.3f}" for r in ratios))


# ===========================================================================
# E3 — exact CVP: algorithmic wall or information-theoretic wall?
# ===========================================================================

print()
print("-" * 78)
print("E3: exact CVP on the same grid — is the planted vector the optimum?")
print("-" * 78)
print("opt = #seeds with tau == 1 (planted vector IS the closest vector)")
print("rec = #seeds exact CVP returns the true d;  lll = variant-A LLL wins")
print("nu  = ||v_planted|| / GH(L0), the Gaussian-heuristic decodability index")

print(f"\n{'curve':<14} {'K1':>4} {'opt':>5} {'rec':>5} {'lll':>5} "
      f"{'tau_min':>9} {'tau_mean':>9} {'tauraw_m':>9} {'nu=pv/GH':>9} "
      f"{'sec':>7}")
cvp_summary = {}
for label in ("12-bit/2557", "12-bit/2677"):
    curve = CURVES[label]
    p, b, n, lam, G = curve
    for ki, k1 in enumerate(K1_GRID):
        t0 = time.time()
        taus, taus_raw, opt, rec = [], [], 0, 0
        nus = []
        for seed in SEEDS:
            d0 = random.Random(seed + 7777).randint(1, n - 1)
            sigs, k2b = make_sigs(curve, M_SIGS, d0, k1, seed)
            if len(sigs) < M_SIGS:
                continue
            got, tau, tau_raw, dist, pn = run_cvp(sigs, d0, n, lam, k1, k2b)
            taus.append(tau)
            taus_raw.append(tau_raw)
            opt += (tau > 1.0 - 1e-9)
            rec += (got is not None)
            S_K1, S_D, S_K2, _ = scales(n, k1, k2b)
            dim0 = 2 * M_SIGS + 1
            log_det = (M_SIGS * math.log(n * S_K1) + math.log(S_D)
                       + M_SIGS * math.log(S_K2))
            nus.append(pn / gaussian_heuristic(log_det, dim0))
        el = time.time() - t0
        lll_wins = wall_rows[(label, "A")][0][ki]
        cvp_summary[(label, k1)] = (opt, rec, lll_wins, taus, nus)
        print(f"{label:<14} {k1:>4} {opt:>5} {rec:>5} {lll_wins:>5} "
              f"{min(taus):>9.4f} {sum(taus)/len(taus):>9.4f} "
              f"{sum(taus_raw)/len(taus_raw):>9.4f} "
              f"{sum(nus)/len(nus):>9.4f} {el:>7.1f}")

print("\nTotals over the 80 (curve,K1,seed) cells of the grid:")
_o = sum(v[0] for v in cvp_summary.values())
_r = sum(v[1] for v in cvp_summary.values())
_l = sum(v[2] for v in cvp_summary.values())
print(f"  planted vector is the CVP optimum : {_o}/80")
print(f"  exact CVP recovers d              : {_r}/80")
print(f"  variant-A LLL recovers d          : {_l}/80")


# ===========================================================================
# E4 — energy split of the shortest vector, variant B
# ===========================================================================

print()
print("-" * 78)
print("E4: variant B — where does the shortest vector's energy sit now?")
print("-" * 78)
print("(2026-07-29 T5 for variant A: e_d = 1.000 exactly, |sv[m]|/n = 1.0000)")
print(f"\n{'curve':<14} {'K1':>4} {'sv/pv':>7} {'e_k1':>7} {'e_k2':>7} "
      f"{'e_kan':>7} {'kan=S?':>7}")
for label in ("12-bit/2557", "12-bit/2677"):
    curve = CURVES[label]
    p, b, n, lam, G = curve
    for k1 in (2, 8, 24):
        d0 = random.Random(42 + 7777).randint(1, n - 1)
        sigs, k2b = make_sigs(curve, M_SIGS, d0, k1, 42)
        if len(sigs) < M_SIGS:
            continue
        ok, r, sv, pv = run_variant_B(curve, M_SIGS, d0, k1, 42, sigs, k2b)
        m = M_SIGS
        tot = sum(x * x for x in sv) or 1
        e_k1 = sum(sv[i] ** 2 for i in range(m)) / tot
        e_k2 = sum(sv[m + i] ** 2 for i in range(m)) / tot
        e_kan = sv[2 * m] ** 2 / tot
        _, _, _, S_KAN = scales(n, k1, k2b)
        print(f"{label:<14} {k1:>4} {r:>7.3f} {e_k1:>7.3f} {e_k2:>7.3f} "
              f"{e_kan:>7.3f} {str(abs(sv[2*m]) == S_KAN):>7}")


# ===========================================================================
# E5 — the nu law, and the m-regime it predicts
# ===========================================================================

print()
print("-" * 78)
print("E5: nu = ||v_planted|| / GH(L0) as a closed form, and its prediction")
print("-" * 78)
print("""
E3 shows nu is a near-perfect gate: every one of the 30 grid cells with
nu < 1 recovered d, and every one of the 10 cells with nu > 2.1 failed.
nu depends only on (n, K1, K2, m) -- NOT on the curve -- which is exactly the
kind of invariant 2026-07-29 T5 argued had to exist (no curve-level invariant
can separate C1 from C2, because the condition is a property of the projected
GS profile and the signature set).

Closed form.  With k1 ~ U{0..K1-1}, k2 ~ U{0..K2-1}, d centered in
[-n/2, n/2), and S_K1 = floor(n/K1), S_K2 = max(1, floor(n/K2)), S_D = 1:

  ||v_p||^2 = m*S_K1^2*(K1-1)(2K1-1)/6 + n^2/12 + m*S_K2^2*(K2-1)(2K2-1)/6
  det(L0)   = (n*S_K1)^m * S_D * S_K2^m = n^(2m) / eff^m,   eff = K1*K2/n
  GH(L0)    = sqrt((2m+1)/(2*pi*e)) * det(L0)^(1/(2m+1))

Taking S_K1 -> n/K1, S_K2 -> n/K2 and K1,K2 -> infinity gives the asymptotic

  nu ~ (n * eff^m)^(1/(2m+1)) * sqrt((2m/3 + 1/12) * 2*pi*e / (2m+1))
     -> (n * eff^m)^(1/(2m+1)) * sqrt(2*pi*e/3)      as m -> infinity.

sqrt(2*pi*e/3) = 2.3864, so nu < 1 needs m*ln(1/eff) > ln n + (2m+1)*ln(2.3864),
i.e. asymptotically  eff < 1/2.3864^2 = 0.1756  REGARDLESS of m.  Above that
bias strength no number of signatures ever makes the planted vector decodable
in this lattice -- a hard ceiling that the information-theoretic bound
m >= ln n / ln(1/eff) (m_thresh in glv_hnp_phase2_20bit.py) completely misses.
""")


def nu_exact(n, m, k1_bound, k2_bound):
    S_K1, S_D, S_K2, _ = scales(n, k1_bound, k2_bound)
    e1 = (k1_bound - 1) * (2 * k1_bound - 1) / 6.0
    e2 = (k2_bound - 1) * (2 * k2_bound - 1) / 6.0
    pn2 = m * S_K1 ** 2 * e1 + (n * S_D) ** 2 / 12.0 + m * S_K2 ** 2 * e2
    dim0 = 2 * m + 1
    log_det = (m * math.log(n * S_K1) + math.log(S_D) + m * math.log(S_K2))
    return math.sqrt(pn2) / gaussian_heuristic(log_det, dim0)


print("Predicted-vs-measured nu at m=12 (measured column from E3):")
print(f"{'curve':<14} {'K1':>4} {'eff':>7} {'nu_pred':>9} {'nu_meas':>9} "
      f"{'ratio':>7}")
for label in ("12-bit/2557", "12-bit/2677"):
    n = CURVES[label][2]
    k2b = math.isqrt(n) + 1
    for k1 in K1_GRID:
        nus = cvp_summary[(label, k1)][4]
        meas = sum(nus) / len(nus)
        pred = nu_exact(n, M_SIGS, k1, k2b)
        print(f"{label:<14} {k1:>4} {k1*k2b/n:>7.4f} {pred:>9.4f} "
              f"{meas:>9.4f} {pred/meas:>7.4f}")

print("\nm* = smallest m with nu_exact < 1  (m_it = information-theoretic "
      "bound ln n / ln(1/eff))")
print(f"{'curve':<14} {'K1':>4} {'eff':>7} {'m_it':>6} {'m*':>6} "
      f"{'nu(12)':>8} {'nu(m*)':>8}")
MSTAR = {}
for label in ("12-bit/2557", "12-bit/2677"):
    n = CURVES[label][2]
    k2b = math.isqrt(n) + 1
    for k1 in K1_GRID:
        eff = k1 * k2b / n
        m_it = math.ceil(math.log(n) / math.log(1.0 / eff)) if eff < 1 else -1
        mstar = None
        for mm in range(4, 4001):
            if nu_exact(n, mm, k1, k2b) < 1.0:
                mstar = mm
                break
        MSTAR[(label, k1)] = mstar
        ms = str(mstar) if mstar else ">4000"
        nv = f"{nu_exact(n, mstar, k1, k2b):8.4f}" if mstar else "     n/a"
        print(f"{label:<14} {k1:>4} {eff:>7.4f} {m_it:>6} {ms:>6} "
              f"{nu_exact(n, 12, k1, k2b):>8.4f} {nv}")

print("""
The 2026-07-29 T4b experiment ("at K1=8 more data does not rescue the lam*=0.07
curve": m = 8,12,16,24,32 -> 0,0,1,0,1 of 5) stopped at m=32.  The nu law says
m* for that cell is far beyond 32.  So T4b did not test its own hypothesis --
it sampled entirely inside the nu > 1 region.  Direct test below.
""")

print("-" * 78)
print("E5b: FALSIFIER — sweep m at the cells T4b declared hopeless")
print("-" * 78)
print("Prediction: recovery switches ON at m ~ m*, not at the m_it bound.")

E5_CELLS = [("12-bit/2677", 8), ("12-bit/2557", 16), ("12-bit/2677", 6)]
E5_SEEDS = [42, 1234, 9999]
E5_M = [12, 24, 40, 60, 80, 100, 120, 160]

print(f"\n{'curve':<14} {'K1':>4} {'m*':>5} " +
      " ".join(f"m={mm:<4}" for mm in E5_M))
for label, k1 in E5_CELLS:
    curve = CURVES[label]
    p, b, n, lam, G = curve
    k2b = math.isqrt(n) + 1
    cells = []
    for mm in E5_M:
        wins = 0
        for seed in E5_SEEDS:
            d0 = random.Random(seed + 7777).randint(1, n - 1)
            sigs, _ = make_sigs(curve, mm, d0, k1, seed)
            if len(sigs) < mm:
                cells.append("skip")
                break
            ok, _, _, _ = run_variant_A(curve, mm, d0, k1, seed, sigs, k2b)
            wins += bool(ok)
        cells.append(f"{wins}/3")
    ms = MSTAR[(label, k1)]
    print(f"{label:<14} {k1:>4} {str(ms):>5} " +
          " ".join(f"{c:<6}" for c in cells))

print(f"\n{'curve':<14} {'K1':>4} " + " ".join(f"m={mm:<4}" for mm in E5_M) +
      "   (nu_exact)")
for label, k1 in E5_CELLS:
    n = CURVES[label][2]
    k2b = math.isqrt(n) + 1
    print(f"{label:<14} {k1:>4} " +
          " ".join(f"{nu_exact(n, mm, k1, k2b):<6.3f}" for mm in E5_M))


# ===========================================================================
# E6 — the coset factorises into m independent 2-D blocks
# ===========================================================================

print()
print("-" * 78)
print("E6: does the recoverable coset factorise?  (explains E5b's flat curves)")
print("-" * 78)
print("""
E5b kills nu as a sufficient condition: for (2677, K1=8) nu drops below 1 at
m=40 and reaches 0.919 at m=160, yet LLL recovery stays flat at ~1/8 for every
m in 12..160.  For (2557, K1=16) nu is ABOVE 1 at every m yet recovery is a
flat ~46%.  Success does not move with m at all -- in either direction.  So
neither nu < 1 nor any dimension-dependent LLL gap factor delta_0^dim explains
the data: both vary by more than an order of magnitude over that sweep.

The structural explanation.  Look at which lattice vectors recover_d() accepts:
a reduced row is accepted iff its Kannan coordinate is +-S_KANNAN and its
d-coordinate is ~ d (mod n).  The d-coordinate of a coset element is d + c_B
where c_B is the coefficient of the B-row, so acceptance forces c_B = 0 for any
row shorter than n.  The set of accepted vectors is therefore the coset

    v_planted + L_lam,   L_lam = < n*S_K1*e_i ; (-lam*S_K1 | S_K2) >_{i<m}

and L_lam is the ORTHOGONAL DIRECT SUM of m copies of one 2-D lattice

    B = < (n*S_K1, 0), (-lam*S_K1, S_K2) >      (same B for every i).

So the minimum over the accepted coset separates:

    min ||.||^2 = S_KANNAN^2 + d_c^2 + sum_i dist( t_i(d), B )^2,
    t_i(d) = ( ((A_i + d*B_i) mod n) * S_K1,  0 ).

Every block is 2-dimensional, and LLL is EXACT in dimension 2 (Gauss/Lagrange).
The lattice problem therefore has NO dimension-2m+2 content whatsoever: its
difficulty is independent of m -- exactly the flatness E5b measured -- and it
depends on lam only through the single 2-D lattice B, which is why lam* shifts
the K1 wall (2026-07-29 T4) without ever creating a hard obstruction.

Consequence, if true: the only genuinely unknown quantity is d, and f(d) below
is a separable score computable in O(m) per candidate.  Phase 2's (2m+2)-dim
lattice would then buy no dimension reduction over exhaustive search on d.

TEST: compute f(d) = sum_i dist(t_i(d), B)^2 for EVERY d in [1,n) and report
the rank of the true d.  rank 1 = the score function identifies d outright.
""")


def gauss_basis_2d(u, v):
    """Lagrange-reduced basis of the 2-D lattice <u, v>."""
    def n2(w):
        return w[0] * w[0] + w[1] * w[1]
    def dt(w, z):
        return w[0] * z[0] + w[1] * z[1]
    if n2(u) > n2(v):
        u, v = v, u
    while True:
        den = n2(u)
        if den == 0:
            break
        num = dt(v, u)
        q = (2 * num + den) // (2 * den) if num >= 0 else -((-2 * num + den) // (2 * den))
        v = (v[0] - q * u[0], v[1] - q * u[1])
        if n2(v) >= n2(u):
            break
        u, v = v, u
    return u, v


def dist2_2d(t, b1, b2):
    """Exact squared distance from t to the 2-D lattice <b1,b2> (reduced).

    Babai rounding on a Lagrange-reduced 2-D basis plus a +-1 neighbourhood is
    exact: the rounding error per coordinate is <= 1/2 and the basis is within
    a factor sqrt(4/3) of orthogonal, so the optimum lies in the 3x3 window.
    """
    d11 = b1[0] * b1[0] + b1[1] * b1[1]
    d22 = b2[0] * b2[0] + b2[1] * b2[1]
    d12 = b1[0] * b2[0] + b1[1] * b2[1]
    det = d11 * d22 - d12 * d12
    if det == 0:
        return float('inf')
    t1 = t[0] * b1[0] + t[1] * b1[1]
    t2 = t[0] * b2[0] + t[1] * b2[1]
    # solve [d11 d12; d12 d22] [x1;x2] = [t1;t2]
    x1 = (t1 * d22 - t2 * d12) / det
    x2 = (t2 * d11 - t1 * d12) / det
    best = None
    r1, r2 = round(x1), round(x2)
    for a in (r1 - 1, r1, r1 + 1):
        for b in (r2 - 1, r2, r2 + 1):
            wx = t[0] - a * b1[0] - b * b2[0]
            wy = t[1] - a * b1[1] - b * b2[1]
            c = wx * wx + wy * wy
            if best is None or c < best:
                best = c
    return best


def score_f(d, sigs, n, lam, S_K1, b1, b2):
    tot = 0
    for s in sigs:
        r = (s['A'] + d * s['B']) % n
        tot += dist2_2d((r * S_K1, 0), b1, b2)
    return tot


print(f"{'curve':<14} {'K1':>4} {'seed':>7} {'LLL':>5} {'rank(d)':>9} "
      f"{'f(d)':>12} {'f_min':>12} {'gap':>7}")
E6_CELLS = [("12-bit/2677", 8), ("12-bit/2677", 6), ("12-bit/2557", 16),
            ("12-bit/2557", 8), ("12-bit/2557", 24), ("12-bit/2677", 24)]
e6_ranks = {}
for label, k1 in E6_CELLS:
    curve = CURVES[label]
    p, b, n, lam, G = curve
    k2b = math.isqrt(n) + 1
    S_K1, S_D, S_K2, S_KAN = scales(n, k1, k2b)
    b1, b2 = gauss_basis_2d((n * S_K1, 0), (-(lam % n) * S_K1, S_K2))
    ranks = []
    for seed in SEEDS:
        d0 = random.Random(seed + 7777).randint(1, n - 1)
        sigs, _ = make_sigs(curve, M_SIGS, d0, k1, seed)
        if len(sigs) < M_SIGS:
            continue
        okA, _, _, _ = run_variant_A(curve, M_SIGS, d0, k1, seed, sigs, k2b)
        f_true = score_f(d0, sigs, n, lam, S_K1, b1, b2)
        better = 0
        f_min = f_true
        for dd in range(1, n):
            if dd == d0:
                continue
            fv = score_f(dd, sigs, n, lam, S_K1, b1, b2)
            if fv < f_min:
                f_min = fv
            if fv < f_true:
                better += 1
        rank = better + 1
        ranks.append(rank)
        gap = (f_true / f_min) if f_min else float('nan')
        print(f"{label:<14} {k1:>4} {seed:>7} {str(bool(okA)):>5} "
              f"{rank:>9} {f_true:>12.4g} {f_min:>12.4g} {gap:>7.3f}")
    e6_ranks[(label, k1)] = ranks

print(f"\n{'curve':<14} {'K1':>4} {'rank(d) over 5 seeds':>34} "
      f"{'#rank==1':>9}")
for label, k1 in E6_CELLS:
    rk = e6_ranks[(label, k1)]
    print(f"{label:<14} {k1:>4} {str(rk):>34} "
          f"{sum(1 for r in rk if r == 1):>9}")

print("""
Read-off: wherever rank(d) == 1 the separable score f already pins down d
outright, so all the information the (2m+2)-dimensional lattice can hold is
present in m independent 2-D problems -- LLL's failure there is purely a
failure to reach the coset minimum, not a missing-information failure.
Wherever rank(d) >> 1 the information genuinely is not there and no reduction
strength can help.
""")

print()
print("Done.")
