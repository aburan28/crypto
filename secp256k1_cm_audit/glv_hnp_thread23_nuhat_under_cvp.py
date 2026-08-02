"""
GLV-HNP Phase 2, Thread 23b: is nu_hat an information-theoretic invariant or a
proxy for LLL's reduction quality?

Background
----------
2026-07-29 (commit e845207) found the first working Phase-2 separator:

    nu_hat = lambda_1(L2) / sqrt(det L2),
    L2 = < (n*S_K1, 0), (-lam*S_K1, S_K2) >   (det = n*S_K1*S_K2, lam-free)

with Spearman(nu_hat, p_hat) = -0.34 .. -0.57 for m = 8..11 and AUC 0.935
against the June C1/C2 classes.  Low nu_hat => higher recovery probability.

Thread 23a (glv_hnp_thread23_bdd.py, EXP C) then showed that on the historical
"hard" curve 12-bit/2677 (lam* = 0.070) the K1 wall MOVES when the Kannan
embedding is replaced by exact CVP enumeration on the (2m+1)-dim lattice:

    K1 wall (>= 3/5 recoveries)     2557 (easy)   2677 (hard)
      Kannan + LLL                       8             4
      Kannan + BKZ(20)                   8             4
      L_cvp + Babai                      8             3
      L_cvp + exact CVP                  8             6      <-- wall moved

Exact CVP returns THE closest lattice point, so a CVP failure means the planted
vector genuinely is not closest.  If nu_hat were an information-theoretic
invariant it would correlate with p_hat_cvp too.  If it is a reduction-quality
proxy, the correlation survives for LLL and collapses for exact CVP.

  H23b : Spearman(nu_hat, p_hat_cvp) is null, while Spearman(nu_hat, p_hat_lll)
         reproduces the 2026-07-29 negative correlation on the same sample.
  Falsifier: if nu_hat correlates equally under exact CVP, nu_hat is
         information-theoretic and H23b is dead.

RESULT (2026-08-02): H23b FALSIFIED.  nu_hat correlates MORE strongly under
exact CVP than under LLL, and not at all with the LLL-vs-CVP gap:

    Spearman(nu_hat, p_hat_LLL)      = -0.567  perm p = 0.0014
    Spearman(nu_hat, p_hat_BKZ30)    = -0.409  perm p = 0.0265
    Spearman(nu_hat, p_hat_exactCVP) = -0.744  perm p < 1e-4   AUC 1.000
    Spearman(nu_hat, p_cvp - p_lll)  = +0.142  perm p = 0.45   (null)
    controls on the same sample: lam* and tau are both null (p >= 0.44).

nu_hat is therefore a property of the lattice geometry, not of the reduction
algorithm.  See glv_hnp_thread23_nuhat_causal.py for the controlled version.

CAVEAT recorded so no later run repeats it: the premise "p_hat_cvp is the
ceiling" is WRONG.  The Kannan arm scans every reduced-basis row with
|last coord| = S_KANNAN, so it is a LIST decoder and can succeed when the
planted point is not the closest one; exact CVP is a UNIQUE decoder.  Neither
dominates -- here p_lll ~ p_cvp (0.267 vs 0.269), in the causal arm
p_lll > p_cvp (0.267 vs 0.198), and on 12-bit/2677 at K1=6 p_cvp >> p_lll
(8/12 vs 1/12).

Design
------
NC fresh 17-bit j=0 GLV curves, m = 12, SEEDS seeds, eff pinned per curve to
EFF_TARGET (the discriminating band tau ~ 1.3 found in Thread 23a).  Three arms
per (curve, seed): Kannan+LLL, Kannan+BKZ(30), L_cvp + exact CVP.
Significance by permutation test on the Spearman statistic (no scipy).

Run: python3 glv_hnp_thread23_nuhat_under_cvp.py
"""

import math
import random
import time

from glv_hnp_thread23_bdd import (
    ec_mul, find_curves, find_generator, gen_signatures, gh_cvp, error_vector,
    hist, norm, run_cvp, run_kannan, scales, lam_star, trial_d,
)

NC = 30
M_SIGS = 12
SEEDS = [42, 1234, 9999, 555, 31337, 7, 20260802, 4242, 31, 99991, 271828, 161803]
EFF_TARGET = 0.13
BKZ_BETA = 30

# ---------------------------------------------------------------------------
# nu_hat
# ---------------------------------------------------------------------------

def gauss_reduce_2d(u, v):
    """Exact shortest nonzero vector of the 2D lattice <u, v> (Lagrange)."""
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
        q = (2 * num + den) // (2 * den) if num >= 0 else -((-2 * num + den) // (2 * den))
        v = (v[0] - q * u[0], v[1] - q * u[1])
        if nrm2(v) >= nrm2(u):
            break
        u, v = v, u
    return u

def nu_hat(n, lam, k1_bound, k2_bound):
    """lambda_1(L2)/sqrt(det L2), det L2 = n*S_K1*S_K2 (independent of lam)."""
    S_K1, _, S_K2, _ = scales(n, k1_bound, k2_bound)
    w = gauss_reduce_2d((n * S_K1, 0), (-(lam % n) * S_K1, S_K2))
    return math.sqrt(w[0] ** 2 + w[1] ** 2) / math.sqrt(n * S_K1 * S_K2)

# ---------------------------------------------------------------------------
# Statistics (no scipy)
# ---------------------------------------------------------------------------

def ranks(xs):
    order = sorted(range(len(xs)), key=lambda i: xs[i])
    r = [0.0] * len(xs)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and xs[order[j + 1]] == xs[order[i]]:
            j += 1
        avg = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            r[order[k]] = avg
        i = j + 1
    return r

def pearson(xs, ys):
    nn = len(xs)
    mx, my = sum(xs) / nn, sum(ys) / nn
    num = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    dx = math.sqrt(sum((x - mx) ** 2 for x in xs))
    dy = math.sqrt(sum((y - my) ** 2 for y in ys))
    return num / (dx * dy) if dx and dy else float('nan')

def spearman(xs, ys):
    return pearson(ranks(xs), ranks(ys))

def perm_pvalue(xs, ys, trials=20000, seed=2026):
    obs = spearman(xs, ys)
    if math.isnan(obs):
        return obs, float('nan')
    rng = random.Random(seed)
    ys2 = list(ys)
    hits = 0
    for _ in range(trials):
        rng.shuffle(ys2)
        if abs(spearman(xs, ys2)) >= abs(obs) - 1e-15:
            hits += 1
    return obs, (hits + 1) / (trials + 1)

def auc(labels, scores):
    """P(score of a positive > score of a negative), ties counted 0.5."""
    pos = [s for l, s in zip(labels, scores) if l]
    neg = [s for l, s in zip(labels, scores) if not l]
    if not pos or not neg:
        return float('nan')
    tot = 0.0
    for a in pos:
        for b in neg:
            tot += 1.0 if a > b else (0.5 if a == b else 0.0)
    return tot / (len(pos) * len(neg))


def main():
    print("=" * 78)
    print("Thread 23b — nu_hat under exact CVP vs under LLL")
    print("=" * 78)
    print(f"NC={NC} fresh 17-bit curves, m={M_SIGS}, {len(SEEDS)} seeds, "
          f"eff_target={EFF_TARGET}, BKZ beta={BKZ_BETA}")

    t0 = time.time()
    curves = find_curves(1 << 16, 1 << 17, want=NC)
    print(f"harvested {len(curves)} curves in {time.time()-t0:.1f}s")

    rows = []
    t0 = time.time()
    for ci, cur in enumerate(curves):
        p, b, n, lam, G = cur
        k2 = math.isqrt(n) + 1
        k1 = max(2, round(EFF_TARGET * n / k2))
        eff = k1 * k2 / n
        nh = nu_hat(n, lam, k1, k2)
        ls = lam_star(lam, n)
        w_lll = w_bkz = w_cvp = 0
        taus = []
        for seed in SEEDS:
            d = trial_d(seed, n)
            sigs = gen_signatures(G, d, M_SIGS, n, lam, p, k1, k2, seed)
            if len(sigs) < M_SIGS:
                continue
            taus.append(norm(error_vector(sigs, d, n, k1, k2)) / gh_cvp(M_SIGS, n, k1, k2))
            w_lll += bool(run_kannan(sigs, n, lam, k1, k2, d)[0])
            w_bkz += bool(run_kannan(sigs, n, lam, k1, k2, d, beta=BKZ_BETA)[0])
            w_cvp += bool(run_cvp(sigs, n, lam, k1, k2, d, exact=True)[0])
        ns = len(taus)
        rows.append(dict(n=n, lam=lam, k1=k1, eff=eff, nu=nh, ls=ls,
                         tau=sum(taus) / ns, ns=ns,
                         p_lll=w_lll / ns, p_bkz=w_bkz / ns, p_cvp=w_cvp / ns))
        print(f"  [{ci+1:>2}/{len(curves)}] n={n:<7} lam*={ls:.3f} nu={nh:.3f} "
              f"K1={k1:<4} eff={eff:.3f} tau={rows[-1]['tau']:.3f}  "
              f"LLL {w_lll:>2}/{ns}  BKZ{BKZ_BETA} {w_bkz:>2}/{ns}  CVP {w_cvp:>2}/{ns}")
    print(f"trials done in {time.time()-t0:.1f}s")

    # ---------------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("Aggregate success rates")
    print("-" * 78)
    for key, name in (("p_lll", "Kannan + LLL"), ("p_bkz", f"Kannan + BKZ({BKZ_BETA})"),
                      ("p_cvp", "L_cvp + exact CVP")):
        mean = sum(r[key] for r in rows) / len(rows)
        print(f"  {name:<22} mean p_hat = {mean:.3f}")
    gap = sum(r["p_cvp"] - r["p_lll"] for r in rows) / len(rows)
    print(f"  reduction loss (CVP - LLL)  = {gap:+.3f}")
    print(f"  curves where CVP > LLL: "
          f"{sum(1 for r in rows if r['p_cvp'] > r['p_lll'])}/{len(rows)};  "
          f"CVP < LLL: {sum(1 for r in rows if r['p_cvp'] < r['p_lll'])}")

    print("\n" + "-" * 78)
    print("H23b: does nu_hat still predict once the reduction is exact?")
    print("-" * 78)
    nus = [r["nu"] for r in rows]
    print(f"{'target':<24}{'spearman':>10}{'perm p':>10}")
    for key, name in (("p_lll", "p_hat (Kannan+LLL)"), ("p_bkz", f"p_hat (Kannan+BKZ{BKZ_BETA})"),
                      ("p_cvp", "p_hat (exact CVP)")):
        s, pv = perm_pvalue(nus, [r[key] for r in rows])
        print(f"{name:<24}{s:>10.3f}{pv:>10.4f}")
    s, pv = perm_pvalue(nus, [r["p_cvp"] - r["p_lll"] for r in rows])
    print(f"{'reduction loss CVP-LLL':<24}{s:>10.3f}{pv:>10.4f}")

    print("\nControls (same sample):")
    for var, name in (("ls", "lam*"), ("tau", "tau")):
        for key in ("p_lll", "p_cvp"):
            s, pv = perm_pvalue([r[var] for r in rows], [r[key] for r in rows])
            print(f"  spearman({name:<5}, {key}) = {s:>6.3f}  perm p = {pv:.4f}")

    print("\nAUC of nu_hat for 'curve recovers on a majority of seeds':")
    for key in ("p_lll", "p_bkz", "p_cvp"):
        lab = [r[key] >= 0.5 for r in rows]
        # nu_hat predicts success with a NEGATIVE sign, so score = -nu_hat
        a = auc(lab, [-r["nu"] for r in rows])
        print(f"  {key}: {sum(lab)}/{len(rows)} majority-success, AUC(-nu_hat) = {a:.3f}")

    print("\n" + "-" * 78)
    print("Historical anchors under exact CVP (12 seeds, K1 grid)")
    print("-" * 78)
    print(f"{'curve':<20}{'K1':>4}{'nu_hat':>9}{'tau':>8}{'LLL':>8}{'CVP':>8}")
    for label, curve, _k1h, m in hist:
        p, b, n, lam, G = curve
        k2 = math.isqrt(n) + 1
        for k1 in (2, 3, 4, 6, 8, 12):
            wl = wc = 0
            taus = []
            for seed in SEEDS:
                d = trial_d(seed, n)
                sigs = gen_signatures(G, d, m, n, lam, p, k1, k2, seed)
                if len(sigs) < m:
                    continue
                taus.append(norm(error_vector(sigs, d, n, k1, k2)) / gh_cvp(m, n, k1, k2))
                wl += bool(run_kannan(sigs, n, lam, k1, k2, d)[0])
                wc += bool(run_cvp(sigs, n, lam, k1, k2, d, exact=True)[0])
            ns = len(taus)
            print(f"{label:<20}{k1:>4}{nu_hat(n, lam, k1, k2):>9.3f}"
                  f"{sum(taus)/ns:>8.3f}{f'{wl}/{ns}':>8}{f'{wc}/{ns}':>8}")

    print("\n" + "=" * 78)
    print("done")
    print("=" * 78)


if __name__ == "__main__":
    main()
