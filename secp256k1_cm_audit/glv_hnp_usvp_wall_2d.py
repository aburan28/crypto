"""
GLV-HNP — Thread 23, part 3: the (nu, nu_hat) plane.

Parts 1 and 2 established:
  * nu = E||v_planted||/GH(L) orders success (AUC 0.955 pooled, 0.938 held out)
    and its m-independent floor 2.386*sqrt(eff) says which cells no amount of
    data can rescue;
  * at nu and m held fixed by construction, the residual is explained by the
    2026-07-29 rival-sublattice invariant nu_hat = lambda_1(L2)/sqrt(det L2)
    at AUC 0.970, while lam*, eff and nu itself sit near the baseline;
  * BUT at fixed (nu, m) the success rate still fell with n
    (nu=1.20: 1.00 / 0.75 / 0.62 / 0.62 at 14 / 17 / 20 / 22 bits).

The last point is the open question this script settles: is that n-dependence a
real third effect, or is it a nu_hat confound?  Part 2 used only 5 curves per
stratum, and the 17-bit sample happened to carry the highest nu_hat values
(mean nu_hat of the total-failure rows was 0.813).

Design: 4 strata x 10 curves x 6 nu targets x 8 seeds at m=12, all with nu
fixed by construction.  Then
  (G1) the pooled 2-D success table over (nu, nu_hat);
  (G2) success rate by stratum WITHIN nu_hat terciles -- if the n-effect is a
       confound it should vanish here;
  (G3) nu_hat distribution by stratum, to see whether the confound is real.

Run: python3 glv_hnp_usvp_wall_2d.py
"""

import math
import os
import time

_HERE = os.path.dirname(os.path.abspath(__file__))
_SRC = open(os.path.join(_HERE, "glv_hnp_usvp_wall_residual.py")).read()
_CUT = _SRC.index('\nprint("=" * 78)\nprint("Thread 23 part 2')
_lib = {"__name__": "_t23p2", "__file__": os.path.join(_HERE, "glv_hnp_usvp_wall_residual.py")}
exec(compile(_SRC[:_CUT], "glv_hnp_usvp_wall_residual.py(prefix)", "exec"), _lib)

nu_of = _lib["nu_of"]
nu_inf = _lib["nu_inf"]
nu_hat = _lib["nu_hat"]
k1_for_nu = _lib["k1_for_nu"]
rate = _lib["rate"]
auc = _lib["auc"]
best_threshold = _lib["best_threshold"]
lam_star = _lib["lam_star"]
search_curves = _lib["search_curves"]

T0 = time.time()
def stamp():
    return f"[{time.time() - T0:6.1f}s]"

SEEDS8 = [42, 1234, 9999, 555, 31337, 271828, 141421, 161803]
BIT_STRATA = [14, 17, 20, 22]
NU_TARGETS = [0.90, 1.05, 1.20, 1.35, 1.50, 1.70]
M_FIXED = 12
PER_BIN, NBINS = 2, 5           # 10 curves per stratum

print("=" * 78)
print("Thread 23 part 3 — the (nu, nu_hat) plane")
print("=" * 78)

rows = []
for nb in BIT_STRATA:
    print(f"{stamp()} {nb}-bit: searching curves ...", flush=True)
    cs = search_curves(2 ** (nb - 1), 2 ** nb, per_bin=PER_BIN, nbins=NBINS)
    print(f"   {len(cs)} curves", flush=True)
    for nu_t in NU_TARGETS:
        for (p, b, n, lam, G) in cs:
            k1, k2 = k1_for_nu(n, M_FIXED, nu_t)
            if k1 is None:
                continue
            w, t, _ = rate((p, b, n, lam, G), M_FIXED, k1, SEEDS8)
            if t == 0:
                continue
            rows.append({"bits": nb, "n": n, "lam": lam, "k1": k1, "k2": k2,
                         "eff": k1 * k2 / n, "nu": nu_of(M_FIXED, n, k1, k2),
                         "nu_target": nu_t, "nu_hat": nu_hat(n, lam, k1, k2),
                         "lam_star": lam_star(lam, n), "wins": w, "total": t,
                         "rate": w / t, "ok": w == t})
    print(f"   {stamp()} {nb}-bit done ({len(rows)} rows)", flush=True)

print(f"\nTotal rows: {len(rows)}  trials: {sum(r['total'] for r in rows)}")


def tab(rs):
    w = sum(r["wins"] for r in rs); t = sum(r["total"] for r in rs)
    return f"{w/t:.2f}({t})" if t else "-"


# ---------------------------------------------------------------------------
# G1 — the 2-D table
# ---------------------------------------------------------------------------
print("\n" + "=" * 78)
print("G1  pooled success rate over (nu, nu_hat)   [all strata, m=12]")
print("=" * 78)

NH_BINS = [0.0, 0.35, 0.50, 0.65, 0.80, 1.01]
print(f"\n{'nu':>7} | " + " ".join(
    f"{'nh<'+format(NH_BINS[i+1],'.2f'):>11}" for i in range(len(NH_BINS) - 1))
    + f" | {'row':>11}")
print("-" * 78)
for nu_t in NU_TARGETS:
    sub = [r for r in rows if r["nu_target"] == nu_t]
    cells = []
    for i in range(len(NH_BINS) - 1):
        cells.append(tab([r for r in sub
                          if NH_BINS[i] <= r["nu_hat"] < NH_BINS[i + 1]]))
    print(f"{nu_t:>7.2f} | " + " ".join(f"{c:>11}" for c in cells)
          + f" | {tab(sub):>11}")
print("-" * 78)
cells = []
for i in range(len(NH_BINS) - 1):
    cells.append(tab([r for r in rows
                      if NH_BINS[i] <= r["nu_hat"] < NH_BINS[i + 1]]))
print(f"{'col':>7} | " + " ".join(f"{c:>11}" for c in cells)
      + f" | {tab(rows):>11}")

print("\nPredictor AUC on all rows (positive = 8/8 recovery):")
print(f"{'predictor':>10} {'AUC':>7} {'best acc':>9} {'thr':>8} {'baseline':>9}")
for key, lo in (("nu", True), ("nu_hat", True), ("eff", True),
                ("bits", True), ("lam_star", True)):
    a, npos, nneg = auc(rows, key, lower_is_better=lo)
    acc, thr, base = best_threshold(rows, key, lower_is_better=lo)
    acc2, thr2, _ = best_threshold(rows, key, lower_is_better=False)
    if acc2 > acc:
        acc, thr = acc2, thr2
    print(f"{key:>10} {a:>7.3f} {acc:>9.3f} "
          f"{(format(thr,'.4f') if isinstance(thr,float) else str(thr)):>8} "
          f"{base:>9.3f}")
print(f"  (positives {npos}, negatives {nneg})")

# joint cut
print("\nJoint rule (nu <= a AND nu_hat <= b): best over a grid")
best = (0.0, None, None)
for a_ in sorted({round(r["nu"], 3) for r in rows}):
    for b_ in sorted({round(r["nu_hat"], 3) for r in rows}):
        acc = sum(1 for r in rows
                  if r["ok"] == (r["nu"] <= a_ and r["nu_hat"] <= b_)) / len(rows)
        if acc > best[0]:
            best = (acc, a_, b_)
print(f"  best accuracy {best[0]:.3f} at nu <= {best[1]:.3f} and "
      f"nu_hat <= {best[2]:.3f}   (single-predictor best was "
      f"{max(best_threshold(rows,'nu_hat')[0], best_threshold(rows,'nu')[0]):.3f})")


# ---------------------------------------------------------------------------
# G2 — does the n-effect survive conditioning on nu_hat?
# ---------------------------------------------------------------------------
print("\n" + "=" * 78)
print("G2  success by stratum WITHIN nu_hat terciles (nu fixed by construction)")
print("=" * 78)

nhs = sorted(r["nu_hat"] for r in rows)
q1, q2 = nhs[len(nhs) // 3], nhs[2 * len(nhs) // 3]
print(f"nu_hat terciles: low < {q1:.3f} <= mid < {q2:.3f} <= high")

for name, lo_, hi_ in (("low", 0.0, q1), ("mid", q1, q2), ("high", q2, 1.01)):
    print(f"\n  -- nu_hat {name} [{lo_:.3f}, {hi_:.3f}) --")
    print(f"{'nu':>7} | " + " ".join(f"{nb:>11}b" for nb in BIT_STRATA))
    for nu_t in NU_TARGETS:
        cells = []
        for nb in BIT_STRATA:
            cells.append(tab([r for r in rows
                              if r["nu_target"] == nu_t and r["bits"] == nb
                              and lo_ <= r["nu_hat"] < hi_]))
        print(f"{nu_t:>7.2f} | " + " ".join(f"{c:>12}" for c in cells))
    print(f"{'ALL':>7} | " + " ".join(
        f"{tab([r for r in rows if r['bits']==nb and lo_<=r['nu_hat']<hi_]):>12}"
        for nb in BIT_STRATA))

# AUC of `bits` restricted to each tercile
print("\n  AUC(bits, lower=win) within each nu_hat tercile "
      "(0.5 = no residual n-effect):")
for name, lo_, hi_ in (("low", 0.0, q1), ("mid", q1, q2), ("high", q2, 1.01)):
    sub = [r for r in rows if lo_ <= r["nu_hat"] < hi_]
    a, np_, nn_ = auc(sub, "bits", lower_is_better=True)
    print(f"    {name:>5}: AUC={a:.3f}  (pos {np_}, neg {nn_})")
a, np_, nn_ = auc(rows, "bits", lower_is_better=True)
print(f"    {'all':>5}: AUC={a:.3f}  (pos {np_}, neg {nn_})")


# ---------------------------------------------------------------------------
# G3 — is nu_hat distributed differently across strata?
# ---------------------------------------------------------------------------
print("\n" + "=" * 78)
print("G3  nu_hat distribution by stratum (the suspected confound)")
print("=" * 78)
print(f"{'bits':>6} {'rows':>6} {'min':>7} {'median':>7} {'mean':>7} {'max':>7}")
for nb in BIT_STRATA:
    v = sorted(r["nu_hat"] for r in rows if r["bits"] == nb)
    if not v:
        continue
    print(f"{nb:>6} {len(v):>6} {v[0]:>7.3f} {v[len(v)//2]:>7.3f} "
          f"{sum(v)/len(v):>7.3f} {v[-1]:>7.3f}")

print("\nnu_hat vs lam*, by stratum (is nu_hat itself curve-arithmetic?):")
print(f"{'bits':>6} {'spearman(nu_hat, lam*)':>24}")
def spearman(xs, ys):
    def rk(v):
        s = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        i = 0
        while i < len(s):
            j = i
            while j + 1 < len(s) and v[s[j + 1]] == v[s[i]]:
                j += 1
            avg = (i + j) / 2.0 + 1
            for k in range(i, j + 1):
                r[s[k]] = avg
            i = j + 1
        return r
    a, b = rk(xs), rk(ys)
    ma, mb = sum(a) / len(a), sum(b) / len(b)
    num = sum((x - ma) * (y - mb) for x, y in zip(a, b))
    den = math.sqrt(sum((x - ma) ** 2 for x in a) * sum((y - mb) ** 2 for y in b))
    return num / den if den else float("nan")

for nb in BIT_STRATA:
    sub = [r for r in rows if r["bits"] == nb]
    if len(sub) < 3:
        continue
    print(f"{nb:>6} {spearman([r['nu_hat'] for r in sub], [r['lam_star'] for r in sub]):>24.3f}")
print(f"{'all':>6} {spearman([r['nu_hat'] for r in rows], [r['lam_star'] for r in rows]):>24.3f}")


# ---------------------------------------------------------------------------
# G4 — is there a single combined statistic?
# ---------------------------------------------------------------------------
print("\n" + "=" * 78)
print("G4  combined statistic  Psi = nu * nu_hat^gamma")
print("=" * 78)

best_g = (0.0, None)
print(f"{'gamma':>7} {'AUC':>7} {'best acc':>9}")
g = 0.0
while g <= 3.001:
    for r in rows:
        r["psi"] = r["nu"] * (r["nu_hat"] ** g)
    a, _, _ = auc(rows, "psi", lower_is_better=True)
    acc, thr, base = best_threshold(rows, "psi", lower_is_better=True)
    if round(g, 2) in (0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0):
        print(f"{g:>7.2f} {a:>7.3f} {acc:>9.3f}")
    if a > best_g[0]:
        best_g = (a, g, acc, thr)
    g += 0.05

a_best, g_best, acc_best, thr_best = best_g
print(f"\n  best gamma = {g_best:.2f}:  AUC = {a_best:.3f}, "
      f"accuracy = {acc_best:.3f} at Psi <= {thr_best:.4f}  "
      f"(baseline {base:.3f})")
print(f"  compare: nu alone AUC {auc(rows,'nu')[0]:.3f}, "
      f"nu_hat alone AUC {auc(rows,'nu_hat')[0]:.3f}, "
      f"best rectangular (nu,nu_hat) cut accuracy {best[0]:.3f}")

for r in rows:
    r["psi"] = r["nu"] * (r["nu_hat"] ** g_best)
ps = sorted(r["psi"] for r in rows)
edges = [ps[0]] + [ps[int(len(ps) * f)] for f in (0.15, 0.3, 0.45, 0.6, 0.75, 0.9)] + [ps[-1] * 1.01]
print(f"\n  success rate vs Psi (deciles), split by stratum "
      f"(collapse check for Psi):")
hdr = " ".join(f"{edges[i]:.2f}-{edges[i+1]:.2f}".rjust(12)
               for i in range(len(edges) - 1))
print(f"{'stratum':>8} {hdr}")
for nb in BIT_STRATA + ["ALL"]:
    sub = rows if nb == "ALL" else [r for r in rows if r["bits"] == nb]
    cells = [tab([r for r in sub if edges[i] <= r["psi"] < edges[i + 1]])
             for i in range(len(edges) - 1)]
    print(f"{str(nb):>8} " + " ".join(f"{c:>12}" for c in cells))

print(f"\n{stamp()} done.")
