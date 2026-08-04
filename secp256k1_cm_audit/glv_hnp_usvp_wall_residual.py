"""
GLV-HNP — Thread 23, part 2: what is left over after nu?

`glv_hnp_usvp_wall.py` (T23-A..D) showed that

    nu = E||v_planted|| / GH(L)          [GH = Gaussian heuristic of the
                                          (2m+2)-dim Kannan lattice]

orders success far better than eff (AUC 0.955 vs 0.894 on a pooled 14/17-bit
sample; 0.938 held out at 20 bits), and that the m-independent floor

    nu_inf(r) = sqrt(2*pi*e/3) * sqrt(r) = 2.386 * sqrt(r),   r = K1*K2/n

correctly predicts which (curve, eff) cells can never be rescued by more data.

But two residuals survived:

  R1  the strata do NOT collapse.  At the same nu, 17-bit cells succeed less
      often than 14-bit cells (bin nu in [1.2,1.5): 86/90 at 14 bits vs 30/60
      at 17 bits).  So p_success is NOT a function of nu alone.
  R2  one of the three held-out 20-bit curves (p=524347, n=523969) recovered
      nothing at ANY m up to 40, at eff=0.141 where nu falls to 1.09 -- well
      inside the regime where the other two curves recover 5/5.

This script tests one explanation for each.

  T23-E (for R1): hold nu and m FIXED by construction and vary n across four
        bit strata.  If p_success drops with n at fixed (nu, m), nu is an
        incomplete statistic and the residual is an n-effect, not an artefact
        of the eff/m grid used in T23-B.

  T23-F (for R2): test whether the residual is the 2026-07-29 rival-sublattice
        invariant nu_hat = lambda_1(L2)/sqrt(det L2), L2 = <(n*S_K1,0),
        (-lam*S_K1, S_K2)>, which that run found separates the June C1/C2
        classes at AUC 0.935.  Prediction from that run's sign convention:
        LOW nu_hat (skew L2) => easier.  Here we ask whether nu_hat explains
        the residual of nu, i.e. whether it separates cells that nu ranks
        equally.

Run: python3 glv_hnp_usvp_wall_residual.py
"""

import math
import os
import time

_HERE = os.path.dirname(os.path.abspath(__file__))
_SRC = open(os.path.join(_HERE, "glv_hnp_usvp_wall.py")).read()
_CUT = _SRC.index('print("=" * 78)')
_lib = {"__name__": "_t23lib", "__file__": os.path.join(_HERE, "glv_hnp_usvp_wall.py")}
exec(compile(_SRC[:_CUT], "glv_hnp_usvp_wall.py(prefix)", "exec"), _lib)

nu_of = _lib["nu_of"]
nu_inf = _lib["nu_inf"]
gh = _lib["gh"]
rate = _lib["rate"]
auc = _lib["auc"]
best_threshold = _lib["best_threshold"]
scales = _lib["scales"]
lam_star = _lib["lam_star"]
lambda_block_mu = _lib["lambda_block_mu"]
search_curves = _lib["search_curves"]
build_curve = _lib["_lib"]["build_curve"] if "_lib" in _lib else None
_t20 = _lib["_lib"]
build_curve = _t20["build_curve"]

T0 = time.time()
def stamp():
    return f"[{time.time() - T0:6.1f}s]"

SEEDS8 = [42, 1234, 9999, 555, 31337, 271828, 141421, 161803]


def nu_hat(n, lam, k1, k2):
    """lambda_1(L2)/sqrt(det L2) for the rival 2D sublattice (2026-07-29)."""
    S_K1, _S_D, S_K2, _S_KAN = scales(n, k1, k2)
    l1, _w = lambda_block_mu(n, lam, S_K1, S_K2)
    return l1 / math.sqrt(n * S_K1 * S_K2)


def k1_for_nu(n, m, nu_target):
    """Largest K1 with nu(m,n,K1,K2) <= nu_target (K2 = isqrt(n)+1)."""
    k2 = math.isqrt(n) + 1
    lo, hi = 2, max(3, n // k2)
    if nu_of(m, n, lo, k2) > nu_target:
        return None, k2
    while lo < hi:
        mid = (lo + hi + 1) // 2
        if nu_of(m, n, mid, k2) <= nu_target:
            lo = mid
        else:
            hi = mid - 1
    return lo, k2


print("=" * 78)
print("Thread 23 part 2 — the residual after nu")
print("=" * 78)


# ===========================================================================
# T23-E — is p_success a function of nu alone?  (fixed nu and m, varying n)
# ===========================================================================
print("\n" + "=" * 78)
print("T23-E  fixed nu, fixed m, varying n  (H1 collapse test, by construction)")
print("=" * 78)

BIT_STRATA = [14, 17, 20, 22]
NU_TARGETS = [0.90, 1.20, 1.50]
M_FIXED = 12

strata_curves = {}
for nb in BIT_STRATA:
    print(f"{stamp()} searching {nb}-bit curves ...", flush=True)
    cs = search_curves(2 ** (nb - 1), 2 ** nb, per_bin=1, nbins=5)
    strata_curves[nb] = cs
    print(f"   {len(cs)} curves")

e_rows = []
for nu_t in NU_TARGETS:
    print(f"\n--- target nu = {nu_t:.2f}, m = {M_FIXED} ---")
    print(f"{'bits':>5} {'n':>9} {'K1':>7} {'eff':>7} {'nu':>7} "
          f"{'nu_inf':>7} {'nu_hat':>7} {'lam*':>6} {'LLL':>7}")
    for nb in BIT_STRATA:
        agg_w = agg_t = 0
        for (p, b, n, lam, G) in strata_curves[nb]:
            k1, k2 = k1_for_nu(n, M_FIXED, nu_t)
            if k1 is None:
                continue
            nuv = nu_of(M_FIXED, n, k1, k2)
            w, t, _ = rate((p, b, n, lam, G), M_FIXED, k1, SEEDS8)
            if t == 0:
                continue
            agg_w += w
            agg_t += t
            nh = nu_hat(n, lam, k1, k2)
            e_rows.append({"bits": nb, "n": n, "lam": lam, "m": M_FIXED,
                           "k1": k1, "k2": k2, "eff": k1 * k2 / n, "nu": nuv,
                           "nu_target": nu_t, "nu_hat": nh,
                           "lam_star": lam_star(lam, n),
                           "wins": w, "total": t, "rate": w / t, "ok": w == t})
            print(f"{nb:>5} {n:>9} {k1:>7} {k1*k2/n:>7.4f} {nuv:>7.4f} "
                  f"{nu_inf(k1,k2,n):>7.4f} {nh:>7.4f} "
                  f"{lam_star(lam,n):>6.3f} {str(w)+'/'+str(t):>7}")
        if agg_t:
            print(f"{nb:>5} {'POOLED':>9} {'':>7} {'':>7} {'':>7} {'':>7} "
                  f"{'':>7} {'':>6} {str(agg_w)+'/'+str(agg_t):>7}  "
                  f"= {agg_w/agg_t:.3f}   {stamp()}")

print("\nT23-E summary: pooled success rate at fixed (nu, m), by stratum")
print(f"{'nu_target':>10} " + " ".join(f"{nb:>10}b" for nb in BIT_STRATA))
for nu_t in NU_TARGETS:
    cells = []
    for nb in BIT_STRATA:
        sub = [r for r in e_rows if r["nu_target"] == nu_t and r["bits"] == nb]
        if sub:
            w = sum(r["wins"] for r in sub); t = sum(r["total"] for r in sub)
            cells.append(f"{w}/{t}={w/t:.2f}")
        else:
            cells.append("-")
    print(f"{nu_t:>10.2f} " + " ".join(f"{c:>11}" for c in cells))


# ===========================================================================
# T23-F — does nu_hat explain the residual?
# ===========================================================================
print("\n" + "=" * 78)
print("T23-F  nu_hat as the residual predictor (cells matched on nu and m)")
print("=" * 78)

# Within each (nu_target, stratum) cell all rows share nu and m by
# construction, so any within-cell variation in outcome is residual.
print("\nWithin-cell spread (rows sharing nu_target and stratum):")
print(f"{'nu_t':>6} {'bits':>5} {'cells':>6} {'rates':>28} "
      f"{'nu_hat(win)':>12} {'nu_hat(lose)':>13}")
resid = []
for nu_t in NU_TARGETS:
    for nb in BIT_STRATA:
        sub = [r for r in e_rows if r["nu_target"] == nu_t and r["bits"] == nb]
        if len(sub) < 2:
            continue
        rates = sorted(r["rate"] for r in sub)
        if rates[0] == rates[-1]:
            spread = "(no spread)"
            nhw = nhl = float("nan")
        else:
            spread = ""
            win = [r["nu_hat"] for r in sub if r["rate"] >= 0.99]
            lose = [r["nu_hat"] for r in sub if r["rate"] <= 0.01]
            nhw = sum(win) / len(win) if win else float("nan")
            nhl = sum(lose) / len(lose) if lose else float("nan")
            resid.extend(sub)
        print(f"{nu_t:>6.2f} {nb:>5} {len(sub):>6} "
              f"{str([f'{x:.2f}' for x in rates]):>28} "
              f"{nhw:>12.4f} {nhl:>13.4f}  {spread}")

if resid:
    print(f"\nResidual sample (only cells with within-cell spread): {len(resid)} rows")
    print(f"{'predictor':>10} {'AUC(low=win)':>13} {'AUC(high=win)':>14} "
          f"{'best acc':>9} {'baseline':>9}")
    for key in ("nu_hat", "lam_star", "eff", "nu"):
        a_lo, np_, nn_ = auc(resid, key, lower_is_better=True)
        a_hi, _, _ = auc(resid, key, lower_is_better=False)
        acc, thr, base = best_threshold(resid, key, lower_is_better=True)
        acc2, _, _ = best_threshold(resid, key, lower_is_better=False)
        print(f"{key:>10} {a_lo:>13.3f} {a_hi:>14.3f} "
              f"{max(acc,acc2):>9.3f} {base:>9.3f}")
    print(f"  (positives = {np_} full-success rows, negatives = {nn_})")

# Direct check on the three held-out 20-bit T23-C curves plus the two 12-bit
# anchors, at the exact parameters where they disagreed.
print("\nDiagnostic: the cells where nu failed to predict (from T23-A/T23-C)")
CASES = [
    # label,               p,      n,      lam,   K1,  m,  observed
    ("12b/2557 K1=6",      2557,   2659,   1755,  6,   12, "5/5"),
    ("12b/2677 K1=6",      2677,   2647,   185,   6,   12, "1/5"),
    ("12b/2557 K1=8",      2557,   2659,   1755,  8,   12, "5/5"),
    ("12b/2677 K1=8",      2677,   2647,   185,   8,   12, "0/5"),
    ("20b/523717 K1=102",  524599, 523717, None,  102, 26, "5/5"),
    ("20b/525583 K1=102",  524341, 525583, None,  102, 26, "5/5"),
    ("20b/523969 K1=102",  524347, 523969, None,  102, 26, "0/5  <-- R2"),
]
glv_roots = _t20["glv_roots"]
print(f"{'case':>20} {'lam*':>6} {'nu':>7} {'nu_hat':>7} {'observed':>12}")
for label, p, n, lam, k1, m, obs in CASES:
    if lam is None:
        rts = glv_roots(n)
        lam = rts[0] if rts else None
        if lam is None:
            print(f"{label:>20}  (no GLV eigenvalue)")
            continue
    k2 = math.isqrt(n) + 1
    print(f"{label:>20} {lam_star(lam,n):>6.3f} {nu_of(m,n,k1,k2):>7.4f} "
          f"{nu_hat(n,lam,k1,k2):>7.4f} {obs:>12}")

print(f"\n{stamp()} done.")
