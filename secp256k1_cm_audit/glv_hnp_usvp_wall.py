"""
GLV-HNP — Thread 23: is the K1/eff wall a uSVP (Gaussian-heuristic) wall?

Context
-------
2026-07-29 (Thread 20, `d525931`) established three things about the Phase-2
(2m+2)-dimensional lattice built by `glv_hnp_phase2_lambda_threshold.py:254`:

  T3  success is governed by eff = K1*K2/n, not by lam:  eff=0.05 -> 19/20
      curves recover, eff=0.15 -> 3/20, eff=0.25 -> 0/20.
  T4  on a fixed curve the K1 wall sits at K1 ~ 4..16.
  T4b at K1=8 more signatures do NOT rescue it (m = 8,12,16,24,32 -> 0,0,1,0,1).

None of those entries says *why* the wall is where it is, and T4b (data does not
help) was recorded as an unexplained brute fact.  This script tests a single
closed-form explanation.

Hypothesis (T23)
----------------
Recovery is a unique-SVP event in the Kannan-embedded lattice L, so it is
governed by the ratio of the planted vector to the Gaussian heuristic:

    nu(m, K1, K2, n) = E||v_planted|| / GH(L),
    GH(L) = sqrt(dim / (2*pi*e)) * det(L)^(1/dim),      dim = 2m+2.

L is triangular, so det(L) = (n*S_K1)^m * S_D * S_K2^m * S_KANNAN exactly.
With S_K1 = n/K1, S_D = 1, S_K2 = n/K2, S_KANNAN = n and r = eff = K1*K2/n:

    E||v_planted||^2 = n^2 * (2m + 4) / 3
    det(L)^(1/dim)   = n^((2m+1)/(2m+2)) * r^(-m/(2m+2))

    ==>  nu(m, r, n) = sqrt( 2*pi*e * (2m+4) / (3*(2m+2)) )
                       * n^(1/(2m+2)) * r^(m/(2m+2))                      (*)

Two consequences, both falsifiable:

  (H1) COLLAPSE.  Success probability is a function of nu alone.  Curves with
       different (n, m, K1, K2) but equal nu succeed equally often, and nu
       out-predicts eff on a pooled sample that mixes bit sizes and m.

  (H2) FLOOR.  nu(m) is decreasing in m but bounded below by

           nu_inf(r) = sqrt(2*pi*e/3) * sqrt(r) = 2.3866 * sqrt(r),

       which is INDEPENDENT of m and of n.  So if 2.3866*sqrt(r) > nu_crit, no
       number of signatures can ever work; and if it is below nu_crit, the
       attack must switch on at the specific m solving nu(m) = nu_crit.
       (H2) is what would explain T4b.

Design: nu_crit is FIT on the 14- and 17-bit strata (T23-B) and then used to
make a blind, pre-registered prediction of the switch-on m at 20 bits (T23-C).

Run: python3 glv_hnp_usvp_wall.py
"""

import math
import os
import random
import sys
import time

# ---------------------------------------------------------------------------
# Reuse the Thread-20 harness verbatim (functions only -- the module runs its
# own experiments at import time, so cut the file at the first experiment).
# ---------------------------------------------------------------------------
_HERE = os.path.dirname(os.path.abspath(__file__))
_SRC = open(os.path.join(_HERE, "glv_hnp_phase2_lambda_threshold.py")).read()
_CUT = _SRC.index("# EXPERIMENT T1")
_lib = {"__name__": "_t20lib", "__file__": os.path.join(_HERE, "_t20lib.py")}
exec(compile(_SRC[:_CUT], "glv_hnp_phase2_lambda_threshold.py(prefix)", "exec"), _lib)

scales = _lib["scales"]
lam_star = _lib["lam_star"]
lambda_block_mu = _lib["lambda_block_mu"]
planted_norm_expected = _lib["planted_norm_expected"]
planted_vector = _lib["planted_vector"]
gen_signatures = _lib["gen_signatures"]
build_glv_lattice = _lib["build_glv_lattice"]
recover_d = _lib["recover_d"]
norm = _lib["norm"]
search_curves = _lib["search_curves"]
glv_roots = _lib["glv_roots"]

from fpylll import IntegerMatrix, LLL, BKZ

TWO_PI_E = 2 * math.pi * math.e
NU_INF_CONST = math.sqrt(TWO_PI_E / 3.0)          # 2.38658...

SEEDS = [42, 1234, 9999, 555, 31337]
T0 = time.time()


def stamp():
    return f"[{time.time() - T0:6.1f}s]"


# ---------------------------------------------------------------------------
# The predictor
# ---------------------------------------------------------------------------

def log_det(m, n, k1, k2):
    """log det(L) for the triangular Phase-2 lattice, exactly as built."""
    S_K1, S_D, S_K2, S_KAN = scales(n, k1, k2)
    return (m * math.log(n * S_K1) + math.log(S_D)
            + m * math.log(S_K2) + math.log(S_KAN))


def gh(m, n, k1, k2):
    dim = 2 * m + 2
    return math.exp(0.5 * math.log(dim / TWO_PI_E) + log_det(m, n, k1, k2) / dim)


def nu_of(m, n, k1, k2):
    """nu = E||v_planted|| / GH(L), computed from the exact scales (not (*))."""
    S_K1, S_D, S_K2, S_KAN = scales(n, k1, k2)
    pv = planted_norm_expected(m, n, k1, k2, S_K1, S_D, S_K2, S_KAN)
    return pv / gh(m, n, k1, k2)


def nu_closed(m, n, k1, k2):
    """The asymptotic closed form (*), for checking the algebra."""
    r = k1 * k2 / n
    return (math.sqrt(TWO_PI_E * (2 * m + 4) / (3.0 * (2 * m + 2)))
            * n ** (1.0 / (2 * m + 2)) * r ** (m / (2.0 * m + 2)))


def nu_inf(k1, k2, n):
    return NU_INF_CONST * math.sqrt(k1 * k2 / n)


def m_for_nu(n, k1, k2, nu_target, m_hi=400):
    """Smallest m with nu(m) <= nu_target, or None if the floor is above it."""
    for m in range(2, m_hi + 1):
        if nu_of(m, n, k1, k2) <= nu_target:
            return m
    return None


# ---------------------------------------------------------------------------
# Trials
# ---------------------------------------------------------------------------

def trial(curve, m, k1, seed, use_bkz=False, beta=20):
    """One LLL run.  Returns (recovered, nu_actual) where nu_actual uses the
    true planted vector rather than its expectation."""
    p, b, n, lam, G = curve
    k2 = math.isqrt(n) + 1
    rng = random.Random(seed ^ 0x5EED)
    d_secret = rng.randrange(1, n)
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1, k2, seed)
    if len(sigs) < m:
        return None, None
    M = build_glv_lattice(sigs, n, lam, k1, k2)
    dim = 2 * m + 2
    A = IntegerMatrix.from_matrix(M)
    if use_bkz:
        BKZ.reduction(A, BKZ.Param(beta))
    else:
        LLL.reduction(A)
    red = [[A[i][j] for j in range(dim)] for i in range(dim)]
    S_KAN = scales(n, k1, k2)[3]
    ok = recover_d(red, m, n, S_KAN, d_secret) is not None
    nu_act = norm(planted_vector(sigs, d_secret, n, k1, k2)) / gh(m, n, k1, k2)
    return ok, nu_act


def rate(curve, m, k1, seeds, use_bkz=False, beta=20):
    w = t = 0
    nus = []
    for s in seeds:
        ok, nua = trial(curve, m, k1, s, use_bkz, beta)
        if ok is None:
            continue
        t += 1
        w += 1 if ok else 0
        nus.append(nua)
    return w, t, (sum(nus) / len(nus) if nus else float("nan"))


# ---------------------------------------------------------------------------
# Metrics
# ---------------------------------------------------------------------------

def auc(rows, key, positive="ok", lower_is_better=True):
    """AUC of `key` for predicting rows[positive] (bool).  Ties count 0.5."""
    pos = [r[key] for r in rows if r[positive]]
    neg = [r[key] for r in rows if not r[positive]]
    if not pos or not neg:
        return float("nan"), len(pos), len(neg)
    good = 0.0
    for a in pos:
        for b in neg:
            if a == b:
                good += 0.5
            elif (a < b) if lower_is_better else (a > b):
                good += 1.0
    return good / (len(pos) * len(neg)), len(pos), len(neg)


def best_threshold(rows, key, lower_is_better=True):
    """Best single cut; returns (accuracy, threshold, baseline)."""
    vals = sorted({r[key] for r in rows})
    n_tot = len(rows)
    base = max(sum(1 for r in rows if r["ok"]),
               sum(1 for r in rows if not r["ok"])) / n_tot
    best = (0.0, None)
    for thr in vals:
        acc = sum(1 for r in rows
                  if r["ok"] == ((r[key] <= thr) if lower_is_better
                                 else (r[key] >= thr))) / n_tot
        if acc > best[0]:
            best = (acc, thr)
    return best[0], best[1], base


print("=" * 78)
print("Thread 23 — is the GLV-HNP Phase-2 wall a uSVP/Gaussian-heuristic wall?")
print("=" * 78)
print(f"nu_inf constant sqrt(2*pi*e/3) = {NU_INF_CONST:.5f}")


# ===========================================================================
# T23-A — algebra check + the 2026-07-29 T4 anchor grid
# ===========================================================================
print("\n" + "=" * 78)
print("T23-A  algebra check and replication of the 2026-07-29 T4 K1 grid")
print("=" * 78)

print("\nExact nu vs closed form (*)  [n=2557-class, K2=isqrt(n)+1]")
print(f"{'n':>7} {'K1':>4} {'m':>4} {'r':>7} {'nu_exact':>9} {'nu_(*)':>9} {'rel':>8}")
for n_ in (2647, 65269, 524269):
    k2_ = math.isqrt(n_) + 1
    for k1_ in (4, 8, 24):
        for m_ in (6, 12, 32):
            e_, c_ = nu_of(m_, n_, k1_, k2_), nu_closed(m_, n_, k1_, k2_)
            print(f"{n_:>7} {k1_:>4} {m_:>4} {k1_*k2_/n_:>7.4f} "
                  f"{e_:>9.4f} {c_:>9.4f} {abs(e_-c_)/c_:>8.1e}")

# The two historical 12-bit curves used by T4 (2026-07-29 / 2026-07-26).
HIST = [
    ("12-bit/2557", 2557, 2, 2659, 1755),
    ("12-bit/2677", 2677, 2, 2647, 185),
]
K1_GRID = [2, 3, 4, 6, 8, 12, 16, 24]
M_ANCHOR = 12

anchor_rows = []
for label, p, b, n, lam in HIST:
    curve = _lib["build_curve"](p, n)
    if curve is None:
        print(f"  !! could not rebuild {label}")
        continue
    curve = (p, curve[1], n, lam, curve[4])
    k2 = math.isqrt(n) + 1
    print(f"\n{label}: p={p} n={n} lam={lam} lam*={lam_star(lam,n):.3f} "
          f"K2={k2}  (m={M_ANCHOR})")
    print(f"{'K1':>4} {'r=eff':>7} {'nu':>7} {'nu_act':>7} {'nu_inf':>7} {'LLL':>6}")
    for k1 in K1_GRID:
        w, t, nua = rate(curve, M_ANCHOR, k1, SEEDS)
        r_ = k1 * k2 / n
        row = {"label": label, "n": n, "m": M_ANCHOR, "k1": k1, "k2": k2,
               "eff": r_, "nu": nu_of(M_ANCHOR, n, k1, k2), "nu_act": nua,
               "wins": w, "total": t, "ok": w == t, "any": w > 0,
               "lam_star": lam_star(lam, n)}
        anchor_rows.append(row)
        print(f"{k1:>4} {r_:>7.4f} {row['nu']:>7.4f} {nua:>7.4f} "
              f"{nu_inf(k1,k2,n):>7.4f} {str(w)+'/'+str(t):>6}")
    print(f"  {stamp()}")

if anchor_rows:
    a_full = [r for r in anchor_rows if r["ok"]]
    a_none = [r for r in anchor_rows if not r["any"]]
    if a_full and a_none:
        print(f"\n  max nu over full-success cells : {max(r['nu'] for r in a_full):.4f}")
        print(f"  min nu over zero-success cells : {min(r['nu'] for r in a_none):.4f}")


# ===========================================================================
# T23-B — collapse test: pool bit sizes, m and eff; does nu alone predict?
# ===========================================================================
print("\n" + "=" * 78)
print("T23-B  collapse test (fit stratum: 14-bit and 17-bit)")
print("=" * 78)

FIT_BITS = [14, 17]
EFF_GRID = [0.02, 0.05, 0.09, 0.15, 0.22, 0.32]
M_GRID = [6, 10, 16]
PER_BIN = 1
NBINS = 6

fit_rows = []
for nbits in FIT_BITS:
    print(f"\n{stamp()} searching j=0 GLV curves, p in [2^{nbits-1}, 2^{nbits}) ...")
    curves = search_curves(2 ** (nbits - 1), 2 ** nbits,
                           per_bin=PER_BIN, nbins=NBINS)
    print(f"  {len(curves)} curves")
    for (p, b, n, lam, G) in curves:
        k2 = math.isqrt(n) + 1
        for eff_t in EFF_GRID:
            k1 = max(2, int(round(eff_t * n / k2)))
            r_ = k1 * k2 / n
            for m in M_GRID:
                w, t, nua = rate((p, b, n, lam, G), m, k1, SEEDS)
                if t == 0:
                    continue
                fit_rows.append({"bits": nbits, "p": p, "n": n, "m": m,
                                 "k1": k1, "k2": k2, "eff": r_,
                                 "nu": nu_of(m, n, k1, k2), "nu_act": nua,
                                 "nu_inf": nu_inf(k1, k2, n),
                                 "wins": w, "total": t, "ok": w == t,
                                 "rate": w / t,
                                 "lam_star": lam_star(lam, n)})
        print(f"  {stamp()} p={p} n={n} done ({len(fit_rows)} cells)")

print(f"\nPooled cells: {len(fit_rows)}  "
      f"(success {sum(1 for r in fit_rows if r['ok'])})")

print("\nPredictor comparison (cell = full success on all 5 seeds):")
print(f"{'predictor':>10} {'AUC':>7} {'best acc':>9} {'thr':>8} {'baseline':>9}")
for key, lo in (("nu", True), ("eff", True), ("m", False), ("nu_act", True)):
    a, npos, nneg = auc(fit_rows, key, lower_is_better=lo)
    acc, thr, base = best_threshold(fit_rows, key, lower_is_better=lo)
    thr_s = f"{thr:.4f}" if isinstance(thr, float) else str(thr)
    print(f"{key:>10} {a:>7.3f} {acc:>9.3f} {thr_s:>8} {base:>9.3f}")

# nu-binned success rate, split by stratum, to see whether the strata collapse.
print("\nSuccess rate vs nu, by stratum (H1: identical columns):")
BINS = [0.0, 0.4, 0.55, 0.7, 0.85, 1.0, 1.2, 1.5, 2.0, 99]
hdr = "  ".join(f"{BINS[i]:.2f}-{BINS[i+1]:.2f}" for i in range(len(BINS) - 1))
print(f"{'stratum':>14}  {hdr}")
for strat in FIT_BITS + ["ALL"]:
    sub = fit_rows if strat == "ALL" else [r for r in fit_rows if r["bits"] == strat]
    cells = []
    for i in range(len(BINS) - 1):
        b = [r for r in sub if BINS[i] <= r["nu"] < BINS[i + 1]]
        cells.append(f"{sum(r['wins'] for r in b)}/{sum(r['total'] for r in b)}"
                     if b else "-")
    print(f"{str(strat):>14}  " + "  ".join(f"{c:>9}" for c in cells))

# also stratify by m, to check nu absorbs the m-dependence
print("\nSuccess rate vs nu, by m (H1: identical columns):")
print(f"{'m':>14}  {hdr}")
for mm in M_GRID:
    sub = [r for r in fit_rows if r["m"] == mm]
    cells = []
    for i in range(len(BINS) - 1):
        b = [r for r in sub if BINS[i] <= r["nu"] < BINS[i + 1]]
        cells.append(f"{sum(r['wins'] for r in b)}/{sum(r['total'] for r in b)}"
                     if b else "-")
    print(f"{mm:>14}  " + "  ".join(f"{c:>9}" for c in cells))

# nu_crit = midpoint between the largest all-success nu and smallest all-fail nu
succ_nu = [r["nu"] for r in fit_rows if r["ok"]]
fail_nu = [r["nu"] for r in fit_rows if r["wins"] == 0]
NU_CRIT = best_threshold(fit_rows, "nu", lower_is_better=True)[1]
print(f"\n  max nu with 5/5      : {max(succ_nu):.4f}" if succ_nu else "")
print(f"  min nu with 0/5      : {min(fail_nu):.4f}" if fail_nu else "")
print(f"  NU_CRIT (best cut)   : {NU_CRIT:.4f}")


# ===========================================================================
# T23-C — blind prediction at 20 bits: where does the attack switch on in m?
# ===========================================================================
print("\n" + "=" * 78)
print("T23-C  held-out prediction at 20 bits using NU_CRIT from T23-B")
print("=" * 78)

M_SWEEP = [4, 6, 8, 10, 12, 16, 20, 26, 32, 40]
pred_rows = []

print(f"\n{stamp()} searching 20-bit curves ...")
c20 = search_curves(2 ** 19, 2 ** 20, per_bin=1, nbins=4)
print(f"  {len(c20)} curves")

# Pick eff values that straddle the floor: nu_inf(r) vs NU_CRIT.
#   r_floor = (NU_CRIT / 2.3866)^2  -- above it, H2 says NO m ever works.
R_FLOOR = (NU_CRIT / NU_INF_CONST) ** 2
print(f"  predicted floor  r_floor = (NU_CRIT/{NU_INF_CONST:.4f})^2 = {R_FLOOR:.4f}")
EFF_C = sorted({round(R_FLOOR * f, 4) for f in (0.45, 0.8, 1.25, 2.0)})
print(f"  eff grid (as multiples 0.45/0.8/1.25/2.0 of r_floor): {EFF_C}")

for (p, b, n, lam, G) in c20[:3]:
    k2 = math.isqrt(n) + 1
    for eff_t in EFF_C:
        k1 = max(2, int(round(eff_t * n / k2)))
        r_ = k1 * k2 / n
        ninf = nu_inf(k1, k2, n)
        m_pred = m_for_nu(n, k1, k2, NU_CRIT)
        print(f"\n  p={p} n={n} K1={k1} eff={r_:.4f} nu_inf={ninf:.4f}  "
              f"PREDICTION: " +
              (f"switch on at m={m_pred}" if m_pred
               else "never (floor above NU_CRIT)"))
        print(f"    {'m':>4} {'nu':>7} {'pred':>5} {'LLL':>6}")
        first_ok = None
        for m in M_SWEEP:
            w, t, nua = rate((p, b, n, lam, G), m, k1, SEEDS)
            nuv = nu_of(m, n, k1, k2)
            pr = "OK" if nuv <= NU_CRIT else "-"
            if w == t and first_ok is None:
                first_ok = m
            pred_rows.append({"p": p, "n": n, "m": m, "k1": k1, "eff": r_,
                              "nu": nuv, "nu_inf": ninf, "wins": w, "total": t,
                              "ok": w == t, "pred": nuv <= NU_CRIT})
            print(f"    {m:>4} {nuv:>7.4f} {pr:>5} {str(w)+'/'+str(t):>6}")
        print(f"    observed switch-on m = {first_ok}   {stamp()}")

if pred_rows:
    tp = sum(1 for r in pred_rows if r["pred"] and r["ok"])
    fp = sum(1 for r in pred_rows if r["pred"] and not r["ok"])
    fn = sum(1 for r in pred_rows if not r["pred"] and r["ok"])
    tn = sum(1 for r in pred_rows if not r["pred"] and not r["ok"])
    print(f"\n  held-out confusion (predict 5/5 iff nu <= {NU_CRIT:.4f}):")
    print(f"    TP={tp}  FP={fp}  FN={fn}  TN={tn}   "
          f"acc={(tp+tn)/len(pred_rows):.3f}  n={len(pred_rows)}")
    a, _, _ = auc(pred_rows, "nu")
    print(f"    AUC(nu) held out = {a:.3f}")


# ===========================================================================
# T23-D — what the floor means at 256 bits
# ===========================================================================
print("\n" + "=" * 78)
print("T23-D  extrapolation of the floor to cryptographic size")
print("=" * 78)
print(f"""
H2 floor:   nu_inf(r) = {NU_INF_CONST:.4f} * sqrt(r),   r = K1*K2/n.
No number of signatures can push nu below nu_inf, so the attack is capped at

    r  <  r_floor = (NU_CRIT / {NU_INF_CONST:.4f})^2 = {R_FLOOR:.4f}
    log2(K1) + log2(K2)  <  log2(n) - {(-math.log2(R_FLOOR)):.2f}

i.e. the two GLV nonce halves together must lose ~{(-math.log2(R_FLOOR)):.1f} bits
relative to the modulus, independent of how many signatures are collected.
""")
for nb in (20, 64, 128, 256):
    n_ = 2 ** nb
    k2_ = math.isqrt(n_)
    k1_max = R_FLOOR * n_ / k2_
    print(f"  n ~ 2^{nb:<4}  K2 = sqrt(n) = 2^{nb/2:<5.1f}  "
          f"=> K1 < 2^{math.log2(k1_max):.2f}  "
          f"(k1 must be {nb/2 - math.log2(k1_max):.2f} bits short of sqrt(n))")

print(f"\n{stamp()} done.")
