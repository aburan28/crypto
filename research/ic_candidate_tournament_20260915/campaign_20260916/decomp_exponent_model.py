#!/usr/bin/env python3
"""Exponent audit for decomposition strategies on binary Koblitz curves (model only).

Nothing here is a measurement.  It evaluates closed-form cost laws, each stated below with the
assumption it rests on, so that DECOMPOSITION-SURVEY.md can say which strategy could move the
exponent against a matched (signed-Frobenius) rho and under what condition.

Part A  generic j-table decomposition (the tournament's current architecture)
Part B  algebraic subspace decomposition: the exact threshold on the per-trial solve cost
Part C  algebraic subspace decomposition priced under four solve-cost laws, against matched rho

Units: log2 of "operations".  Group operations, table probes and Macaulay-matrix F_2 operations
are NOT converted into one another (the repo's ledgers do that at measured ratios); a conversion
moves every algebraic row by O(log n) bits and moves no exponent.  Constants below 2x are
dropped.  Matched rho = sqrt(pi*r/(4n)) with r = 2^(n-1) (cofactor 2), the signed-Frobenius walk.
"""
import math

LOG2 = math.log2


def log2_binom_upto(N, D):
    """log2 of sum_{i<=D} C(N, i): the number of multilinear monomials of degree <= D."""
    D = max(0, min(int(D), N))
    s, term = 1.0, 1.0
    for i in range(1, D + 1):
        term = term * (N - i + 1) / i
        s += term
    return LOG2(s)


def log2_add(a, b):
    hi, lo = max(a, b), min(a, b)
    return hi + LOG2(1 + 2 ** (lo - hi))


def rho_log2(n):
    r_log2 = n - 1
    return 0.5 * (LOG2(math.pi) + r_log2 - LOG2(4 * n))


# --------------------------------------------------------------------------------------------
print("=" * 100)
print("Part A. Generic j-table decomposition: table of all j-sums of an F-point base, (m-j) summands")
print("looped per trial, one relation per column.  Cost F^j + r/F^(j-1)  =>  r^(j/(2j-1)).")
print("Holds for any m >= j; the summand count only moves constants.  Rho is r^(1/2).")
print("=" * 100)
print(f"{'j':>3} {'table':>14} {'exponent':>9} {'ratio IC/rho grows as':>24}  memory at optimum")
for j, name in ((1, "base only"), (2, "pair table"), (3, "triple table"), (4, "4-sum table"),
                (6, "6-sum table")):
    e = j / (2 * j - 1)
    print(f"{j:>3} {name:>14} {e:>9.4f} {'r^' + format(e - 0.5, '.4f'):>24}  r^{j / (2 * j - 1):.3f} entries")
print("  j -> inf: exponent -> 1/2 from above, never below.  Every such method is generic (it uses the")
print("  group law, equality of encodings and the Frobenius, which is [lambda] and costs O(log r) group")
print("  operations to simulate), so Shoup's Omega(sqrt r) generic lower bound (Eurocrypt 1997) forbids")
print("  an exponent below 1/2 for the whole family, whatever the table, fold, walk or summand count.")
print("  Measured two-point rate of the pair arm vs the unmatched rho (RESULTS.md round 0021): r^0.153;")
print("  model r^(1/6) = r^0.1667.")

# --------------------------------------------------------------------------------------------
print()
print("=" * 100)
print("Part B. Algebraic subspace decomposition, free-form solve cost.")
print("Base: points with x in an l-dim F_2-subspace V (|F| ~ 2^l).  H1: a target decomposes with")
print("probability min(1, 2^(m l)/(m! 2^n)) (the repo's yield law, confirmed to 1.09x on 12 toy cells).")
print("Relations needed ~ |F|; linear algebra ~ |F|^2.  Per-trial solve cost T = 2^(c n).")
print("With x = l/n the log-cost/n is max(1 - (m-1)x + c, 2x) for x <= 1/m, so:")
print("   E(c, m) = 2(1+c)/(m+1)      if c <= 1/m      (balanced optimum x = (1+c)/(m+1))")
print("   E(c, m) = 1/m + c           if c >  1/m      (saturated optimum x = 1/m)")
print("E < 1/2 (an exponent below rho) iff c < c*(m):")
print("=" * 100)


def E(c, m):
    return 2 * (1 + c) / (m + 1) if c <= 1 / m else 1 / m + c


def c_star(m):
    # largest c with E(c, m) < 1/2 (supremum)
    if m <= 3:
        return 0.0
    bal = (m - 3) / 4
    return bal if bal <= 1 / m else 0.5 - 1 / m


print(f"{'m':>3} {'c*(m)':>8} {'E(0,m)':>8}   gamma* chained S3 (N=(m-1)n at x=1/m)   gamma* direct S_(m+1) (N=n)")
for m in range(3, 11):
    cs = c_star(m)
    # at the saturated optimum x = 1/m: chained unknowns N = m*l + (m-2)*n = (m-1) n; direct N = m*l = n
    g_ch = cs / (m - 1)
    g_di = cs / 1.0
    print(f"{m:>3} {cs:>8.4f} {E(0, m):>8.4f}   {g_ch:>36.4f}   {g_di:>26.4f}")
print("  m = 3 cannot go below 1/2 even with a free oracle (E(0,3) = 1/2): it ties rho at best.")
print("  gamma* is the largest per-unknown exponent (T = 2^(gamma N)) a solver may have on the")
print("  descended system and still beat rho asymptotically.  Exhaustive search is gamma = 1; the")
print("  generic j-table of Part A is what 'solving' by the group law costs.  So the algebraic route")
print("  needs the descended systems to be far easier than generic Boolean systems of their size.")

# --------------------------------------------------------------------------------------------
print()
print("=" * 100)
print("Part C. Algebraic subspace decomposition priced under four solve laws, fold = n (stable V),")
print("chained S3 unknowns N = m l + (m-2) n, Macaulay cost T = (monomials of degree <= D)^omega.")
print("  L_free : T = 1                         (the repo's free-oracle floor)")
print("  L_const: D = m^2 - m + 1, omega = 2    (Kousidis-Wiemers first-fall bound used AS a solving")
print("                                          degree: the first-fall-degree assumption, disputed)")
print("  L_sem  : D = 4 for every m, omega = 2 (Semaev ePrint 2015/310 Assumption 1 on chained S3;")
print("                                          the assumption Huang-Kosters-Yeo turn into a reductio)")
print("  L_meas : D = 3.5 + l/2,  omega = 2     (fit to this repo's m = 3 refutation degrees 5,6,6,7,7")
print("                                          over l = 3..7; carried to m >= 4 unmeasured, which is")
print("                                          optimistic if m = 4 degrees are higher)")
print("Direct-S_(m+1) variant (N = m l, same D laws) shown as 'dir' for L_const: the optimistic reading.")
print("=" * 100)


def best_cost(n, m, law, direct=False, fold=True):
    best = None
    for l in range(3, n):
        if m * l > n + 40:
            break
        N = m * l if direct else m * l + (m - 2) * n
        if law == "free":
            logT = 0.0
        else:
            D = {"const": m * m - m + 1, "sem": 4}.get(law, 3.5 + l / 2)
            logT = 2 * log2_binom_upto(N, D)
        logF = max(0.0, l - (LOG2(2 * n) if fold else 1))  # columns up to sign (and Frobenius)
        log_trials_per_rel = max(0.0, n + LOG2(math.factorial(m)) - m * l)
        log_rel = logF + log_trials_per_rel + max(logT, 0.0)
        log_la = LOG2(m) + 2 * logF
        tot = log2_add(log_rel, log_la)
        if best is None or tot < best[0]:
            best = (tot, l, logT)
    return best


ns = [31, 61, 131, 163, 283, 571, 1000, 2000, 4000]
MEAS_NOTE = ("  L_meas optimum sits at a tiny base (l = 4-6) because a Macaulay degree growing with l\n"
             "  makes each extra dimension of V cost more than the 2^(m-1) trials it saves: the model then\n"
             "  degenerates to ~2^n (slope 1.0), i.e. worse than the product law's enumeration oracle.")
for law, direct in (("free", False), ("const", False), ("const", True), ("sem", False), ("meas", False)):
    tag = {"free": "L_free", "const": "L_const", "meas": "L_meas", "sem": "L_sem"}[law] + (" dir" if direct else "")
    print(f"\n{tag}: log2(IC cost) - log2(matched rho) at the best l  [best l in brackets]")
    print(f"{'m':>3} " + " ".join(f"{'n=' + str(n):>11}" for n in ns) + "   slope(2000->4000)/n")
    for m in ((3, 4, 5, 6, 8, 12) if law == "sem" else (3, 4, 5, 6)):
        cells, pts = [], {}
        for n in ns:
            tot, l, _ = best_cost(n, m, law, direct)
            pts[n] = tot
            cells.append(f"{tot - rho_log2(n):>+7.1f}[{l:>2}]")
        slope = (pts[4000] - pts[2000]) / 2000
        print(f"{m:>3} " + " ".join(cells) + f"   {slope:.3f}")

print(MEAS_NOTE)
print()
print("Crossover with matched rho, first n on a step-5 grid from 60 (model; assumptions as labelled):")
for law, cases in (("const", [(m, d) for m in (4, 5, 6) for d in (False, True)]),
                   ("sem", [(m, False) for m in (4, 5, 6, 8, 12)])):
    for m, direct in cases:
        hit = None
        for n in range(60, 4001, 5):
            if best_cost(n, m, law, direct)[0] < rho_log2(n):
                hit = n
                break
        kind = "direct S_(m+1)" if direct else "chained S3    "
        print(f"  {'L_const' if law == 'const' else 'L_sem  '} m = {m:>2} {kind}: n ~ {hit}")
print()
print("Reading Part C: negative = cheaper than matched rho.  'slope' is the empirical exponent of the")
print("model's IC cost per unit n at large n (rho's is 0.5).  L_free is a floor nothing reaches;")
print("L_meas is the only law with measured support, and it is m = 3 support only.  Finite-n slopes")
print("need not equal the asymptote E(0, m) = 2/(m+1) of the constant-degree laws.")

# --------------------------------------------------------------------------------------------
print()
print("=" * 100)
print("Part D. The oracle budget at tournament sizes: the most a per-target decomposition solve may")
print("cost (group-operation units; every phase free except the linear algebra m*K^2) for an m-summand")
print("subspace arm to stay under matched rho, maximised over l, fold = n (stable or orbit-union base,")
print("the optimistic case).  Compare the repo's measured per-target oracle prices")
print("(RESEARCH_IC_BOUNDARY_LEDGER.md sec. 5 / BOUNDARY_TARGETS.md Regime B): F4 at n = 23, m = 2")
print("(22 unknowns): 8.8e4 GAE/target; F4 at n = 15, m = 3 (30 unknowns): 2.4e6 GAE/target;")
print("pair-table probe 0.5-25 GAE.")
print("=" * 100)
print(f"{'n':>4} {'m':>2} {'l':>4} {'unknowns':>9} {'log2 trials':>12} {'log2 rho':>9} {'max budget/trial':>17}")
for n in (23, 31, 37, 43, 53, 61):
    for m in (3, 4, 5):
        rho = rho_log2(n)
        best = None
        for l in range(3, n):
            logF = max(0.0, l - LOG2(2 * n))
            log_trials = logF + max(0.0, n + LOG2(math.factorial(m)) - m * l)
            la = LOG2(m) + 2 * logF
            room = 2 ** rho - 2 ** la
            if room <= 0:
                break
            budget = room / 2 ** log_trials
            if best is None or budget > best[0]:
                best = (budget, l, log_trials)
        budget, l, log_trials = best
        unknowns = m * l + (m - 2) * n
        print(f"{n:>4} {m:>2} {l:>4} {unknowns:>9} {log_trials:>12.1f} {rho:>9.1f} {budget:>17.3g}")
print("  'unknowns' is the chained-S3 system each trial must refute at that l.  The measured F4 prices")
print("  above are for systems of 22-30 unknowns; every row here needs a larger system, and the m = 3")
print("  price at 30 unknowns (2.4e6) already exceeds every budget in the table.")
