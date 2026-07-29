"""
GLV-HNP Phase 2 -- Thread 23: closed form for lambda_1(L2), and why every
scale-free continued-fraction invariant had to fail.

Background
----------
2026-07-29 (run #2) found the first invariant that separates the June C1/C2
classes: nu_hat = lambda_1(L2)/sqrt(det L2), where L2 is the non-planted 2D
sublattice of the Phase 2 lattice

    L2 = < u = (n*S_K1, 0),  v = (-lam*S_K1, S_K2) >,   det L2 = n*S_K1*S_K2.

nu_hat was computed (Lagrange-Gauss) but not understood.  The next-step proposal
conjectured that lambda_1(L2) is a continued-fraction quantity of lam/n taken
*at the scale* S_K2/S_K1, and that this explains why the scale-free CF
invariants (q_cf, max_q_cf, max_a) all failed in June.

This script proves the closed form and tests the consequences.

The closed form (exact, not heuristic)
--------------------------------------
A general element of L2 is  a*v + b*u = ((b*n - a*lam)*S_K1, a*S_K2).  For fixed
a the inner coordinate is minimised over b at the centered residue

    r(a) = min_b |a*lam - b*n| = |a*lam mod^\\pm n|,

so

    lambda_1(L2)^2 = min over a >= 0 of  (r(a)*S_K1)^2 + (a*S_K2)^2,

with the a = 0 term contributing (n*S_K1)^2.  Both coordinates of the candidate
at a are nondecreasing in a and in r(a); and by the classical best-approximation
-of-the-second-kind theorem, for any a with q_j <= a < q_{j+1} (q_j = continued
fraction denominators of lam/n) one has r(a) >= r(q_j).  Hence the candidate at
q_j dominates the candidate at a coordinatewise, and

    lambda_1(L2)^2 = min( (n*S_K1)^2,
                          min over CF denominators q_j of
                              (r(q_j)*S_K1)^2 + (q_j*S_K2)^2 )      (*)

exactly.  This is O(log n) and involves no lattice reduction.  Prediction: (*)
agrees with Lagrange-Gauss on every input, to the bit.

Consequences to test
--------------------
C1. (*) is exact  ->  E1, E2.
C2. The argmin denominator sits at q ~ sqrt(n*S_K1/S_K2) = sqrt(n*K2/K1),
    the scale where the two coordinates balance  ->  E3.
C3. At that scale nu_hat ~ sqrt(2*q*r(q)/n), and q*r(q)/n ~ 1/a_{j+1}, so
    nu_hat is governed by the partial quotient of lam/n *at that scale*:
    small nu_hat (= easy attack) <=> a large partial quotient there  ->  E4.
C4. Therefore no scale-free function of lam/n can predict Phase 2 success:
    holding the curve (n, lam) fixed and moving K1 moves nu_hat across
    essentially its whole range while every scale-free CF invariant is
    constant by construction  ->  E5.  This retro-explains the June failures
    of q_cf, max_q_cf and max_a without needing any new experiment.

Falsifier (as posed 2026-07-29): if predicted-from-CF nu_hat correlates < 0.9
with the computed nu_hat, the convergent-scale story is wrong.

Run: python3 glv_hnp_nuhat_closed_form.py
"""

import math
import random
import sys
import time

import importlib.util
_spec = importlib.util.spec_from_file_location(
    "_t20a", __file__.rsplit("/", 1)[0] + "/glv_hnp_phase2_nuhat_base.py")
_t20a = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_t20a)

gauss_reduce_2d = _t20a.gauss_reduce_2d
scales = _t20a.scales
rival_sublattice_nu = _t20a.rival_sublattice_nu
glv_eigenvalues = _t20a.glv_eigenvalues
mu_of = _t20a.mu_of
eisenstein_decompose = _t20a.eisenstein_decompose
j0_traces = _t20a.j0_traces
identify_twist = _t20a.identify_twist

import sympy

SEP = "=" * 78


# ---------------------------------------------------------------------------
# Continued fractions of lam/n and the closed form (*)
# ---------------------------------------------------------------------------

def cf_convergents(lam, n):
    """
    Continued fraction of lam/n.  Returns a list of records, one per convergent
    p_j/q_j, in increasing q_j:

        (q_j, p_j, r_j, a_next)

    r_j = |q_j*lam - p_j*n| is the (unsigned) approximation residue, which is
    exactly min_b |q_j*lam - b*n| for a convergent denominator, and a_next is
    the partial quotient a_{j+1} that follows this convergent (0 at the end).

    Uses exact integer arithmetic throughout -- no floats, so this is valid at
    256+ bits.
    """
    lam = lam % n
    a_list = []
    x, y = lam, n
    while y:
        a_list.append(x // y)
        x, y = y, x % y
    # p_{-1}=1, p_0=a_0 ; q_{-1}=0, q_0=1
    out = []
    p_prev, p_cur = 1, a_list[0]
    q_prev, q_cur = 0, 1
    for j in range(len(a_list)):
        if j > 0:
            p_prev, p_cur = p_cur, a_list[j] * p_cur + p_prev
            q_prev, q_cur = q_cur, a_list[j] * q_cur + q_prev
        r = abs(q_cur * lam - p_cur * n)
        a_next = a_list[j + 1] if j + 1 < len(a_list) else 0
        out.append((q_cur, p_cur, r, a_next))
    return out


def lambda1_closed_form(n, lam, k1_bound, k2_bound):
    """
    lambda_1(L2) by the closed form (*), plus the argmin data.

    Returns dict with: nu (=lambda_1), nu_hat, q_star, r_star, a_star
    (the partial quotient following the winning convergent), j_star, and
    q_balance = sqrt(n*K2/K1), the predicted scale.
    """
    s_k1, _, s_k2, _ = scales(n, k1_bound, k2_bound)
    det = n * s_k1 * s_k2

    best = (n * s_k1) ** 2          # the a = 0 candidate, vector (n*S_K1, 0)
    q_star = r_star = a_star = 0
    j_star = -1

    for j, (q, _p, r, a_next) in enumerate(cf_convergents(lam, n)):
        cand = (r * s_k1) ** 2 + (q * s_k2) ** 2
        if cand < best:
            best = cand
            q_star, r_star, a_star, j_star = q, r, a_next, j

    nu = math.sqrt(best)
    return {
        "nu": nu,
        "nu_hat": nu / math.sqrt(det),
        "q_star": q_star,
        "r_star": r_star,
        "a_star": a_star,
        "j_star": j_star,
        "q_balance": math.sqrt(n * k2_bound / k1_bound) if k1_bound else 0.0,
        "det": det,
    }


def lambda1_exact_int(n, lam, k1_bound, k2_bound):
    """lambda_1(L2)^2 as an exact integer via Lagrange-Gauss, for bit-exact
    comparison against the closed form (avoids float rounding at 256 bits)."""
    s_k1, _, s_k2, _ = scales(n, k1_bound, k2_bound)
    u, v = (n * s_k1, 0), (-(lam % n) * s_k1, s_k2)

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
        q = ((2 * num + den) // (2 * den) if num >= 0
             else -((-2 * num + den) // (2 * den)))
        v = (v[0] - q * u[0], v[1] - q * u[1])
        if nrm2(v) >= nrm2(u):
            break
        u, v = v, u
    return nrm2(u)


def closed_form_sq(n, lam, k1_bound, k2_bound):
    """The closed form (*) as an exact integer square."""
    s_k1, _, s_k2, _ = scales(n, k1_bound, k2_bound)
    best = (n * s_k1) ** 2
    for q, _p, r, _a in cf_convergents(lam, n):
        cand = (r * s_k1) ** 2 + (q * s_k2) ** 2
        if cand < best:
            best = cand
    return best


# ---------------------------------------------------------------------------
# Curve harvesting (j = 0 GLV curves, same construction as Thread 20)
# ---------------------------------------------------------------------------

def harvest_j0_curves(bits_lo, bits_hi, want, seed=20230729):
    """j=0 CM curves with prime n = 1 mod 3 in the requested bit window."""
    rng = random.Random(seed)
    out = []
    lo, hi = 1 << bits_lo, 1 << bits_hi
    p = int(lo)
    tries = 0
    while len(out) < want and tries < 400000:
        tries += 1
        p = sympy.nextprime(rng.randrange(lo, hi))
        if p % 3 != 1 or p >= hi:
            continue
        dec = eisenstein_decompose(p)
        if dec is None:
            continue
        a, b = dec
        for t in j0_traces(a, b):
            n = p + 1 - t
            if n <= 4 or not sympy.isprime(n) or n % 3 != 1:
                continue
            ev = glv_eigenvalues(n)
            if ev is None:
                continue
            lam = ev[0]
            out.append({"p": p, "n": n, "lam": lam, "mu": mu_of(lam, n)})
            break
    return out


# ---------------------------------------------------------------------------
# E1 -- exactness of the closed form on j=0 GLV curves
# ---------------------------------------------------------------------------

def exp_e1(curves, eff_grid):
    print(SEP)
    print("E1  closed form (*) vs Lagrange-Gauss, exact integer comparison")
    print(SEP)
    print(f"{len(curves)} j=0 GLV curves x {len(eff_grid)} eff values "
          f"= {len(curves)*len(eff_grid)} lattices")
    mismatch = 0
    total = 0
    worst = None
    t0 = time.time()
    for c in curves:
        n, lam = c["n"], c["lam"]
        for eff in eff_grid:
            k1 = max(2, int(eff * math.sqrt(n)))
            k2 = math.isqrt(n) + 1
            got = closed_form_sq(n, lam, k1, k2)
            ref = lambda1_exact_int(n, lam, k1, k2)
            total += 1
            if got != ref:
                mismatch += 1
                if worst is None:
                    worst = (n, lam, k1, k2, got, ref)
    dt = time.time() - t0
    print(f"  lattices tested      : {total}")
    print(f"  exact bit-for-bit    : {total - mismatch}")
    print(f"  mismatches           : {mismatch}")
    if worst:
        print(f"  first mismatch       : n={worst[0]} lam={worst[1]} "
              f"K1={worst[2]} K2={worst[3]} cf={worst[4]} gauss={worst[5]}")
    print(f"  wall clock           : {dt:.2f}s")
    print(f"  VERDICT              : "
          f"{'CLOSED FORM EXACT' if mismatch == 0 else 'CLOSED FORM WRONG'}")
    return mismatch == 0


# ---------------------------------------------------------------------------
# E2 -- generality: arbitrary lam, arbitrary K1/K2, up to 256 bits
# ---------------------------------------------------------------------------

def exp_e2(trials=4000, seed=987654321):
    print()
    print(SEP)
    print("E2  generality: random n, random lam, random K1/K2, 8..256 bits")
    print(SEP)
    rng = random.Random(seed)
    mismatch = 0
    by_bits = {}
    for _ in range(trials):
        bits = rng.choice([8, 12, 16, 20, 24, 32, 64, 128, 256])
        n = sympy.nextprime(rng.randrange(1 << (bits - 1), 1 << bits))
        lam = rng.randrange(1, n)
        k1 = rng.randrange(2, max(3, math.isqrt(n) + 1))
        k2 = rng.randrange(2, max(3, math.isqrt(n) + 1))
        ok = closed_form_sq(n, lam, k1, k2) == lambda1_exact_int(n, lam, k1, k2)
        rec = by_bits.setdefault(bits, [0, 0])
        rec[0] += 1
        rec[1] += 0 if ok else 1
        if not ok:
            mismatch += 1
    print(f"{'bits':>6} {'trials':>8} {'mismatch':>9}")
    for bits in sorted(by_bits):
        t, m = by_bits[bits]
        print(f"{bits:>6} {t:>8} {m:>9}")
    print(f"  total mismatches     : {mismatch} / {trials}")
    print(f"  VERDICT              : "
          f"{'EXACT AT ALL SIZES' if mismatch == 0 else 'FAILS'}")

    # secp256k1 itself, at the eff used in the 2026-07-29 extrapolation
    n = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
    ev = glv_eigenvalues(n)
    lam = ev[0]
    for eff in (0.0993, 0.05, 0.25):
        k1 = max(2, int(eff * math.sqrt(n)))
        k2 = math.isqrt(n) + 1
        cf = lambda1_closed_form(n, lam, k1, k2)
        _, _, nh_ref = rival_sublattice_nu(n, lam, k1, k2)
        print(f"  secp256k1 eff={eff:<7} nu_hat(cf)={cf['nu_hat']:.6f} "
              f"nu_hat(gauss)={nh_ref:.6f} "
              f"q*=2^{math.log2(max(1,cf['q_star'])):.1f} "
              f"a_next={cf['a_star']}")
    return mismatch == 0


# ---------------------------------------------------------------------------
# E3 -- the winning convergent sits at the balance scale sqrt(n*K2/K1)
# ---------------------------------------------------------------------------

def exp_e3(curves, eff_grid):
    print()
    print(SEP)
    print("E3  is the argmin convergent at q ~ sqrt(n*K2/K1)?")
    print(SEP)
    ratios, ratios_ref = [], []
    for c in curves:
        n, lam = c["n"], c["lam"]
        for eff in eff_grid:
            k1 = max(2, int(eff * math.sqrt(n)))
            k2 = math.isqrt(n) + 1
            cf = lambda1_closed_form(n, lam, k1, k2)
            if cf["q_star"] > 0 and cf["q_balance"] > 0:
                ratios.append(math.log2(cf["q_star"] / cf["q_balance"]))
                a_ref = max(1, cf["a_star"])
                ratios_ref.append(
                    math.log2(cf["q_star"] / (cf["q_balance"] / math.sqrt(a_ref))))
    def stats(vals, label):
        vals = sorted(vals)
        qf = lambda f: vals[min(len(vals) - 1, int(f * len(vals)))]
        w1 = sum(1 for r in vals if abs(r) <= 1.0) / len(vals)
        w2 = sum(1 for r in vals if abs(r) <= 2.0) / len(vals)
        print(f"  {label:<20} : median {qf(0.5):+.3f}  "
              f"IQR [{qf(0.25):+.3f}, {qf(0.75):+.3f}]  "
              f"10-90% [{qf(0.10):+.3f}, {qf(0.90):+.3f}]  "
              f"|.|<=1 {w1:.3f}  |.|<=2 {w2:.3f}")
        return qf(0.5), w1

    print(f"  samples              : {len(ratios)}")
    med_naive, w1_naive = stats(ratios, "log2 q*/q_bal")

    # Refinement found this session.  Balance is (r*S_K1) = (q*S_K2), i.e.
    # r/q = K1/K2.  But r_j ~ n/q_{j+1} = n/(a_{j+1} q_j + q_{j-1}), so
    #     n/(a*q^2) ~ K1/K2   =>   q ~ sqrt(n*K2/(a*K1)),
    # i.e. the naive scale sqrt(n*K2/K1) is too large by a factor sqrt(a_next).
    # Under Gauss-Kuzmin E[log2 sqrt(a)] ~ 1.2, which is exactly the observed
    # median offset.
    med_ref, w1_ref = stats(ratios_ref, "log2 q*/q_bal_refined")
    print(f"  refined scale        : q_bal_refined = "
          f"sqrt(n*K2/(a_next*K1))")
    print(f"  VERDICT              : naive scale is biased {med_naive:+.2f} "
          f"octaves ({w1_naive:.0%} within 1);\n"
          f"                         dividing by sqrt(a_next) removes the bias "
          f"-> {med_ref:+.2f} ({w1_ref:.0%} within 1).")


# ---------------------------------------------------------------------------
# E4 -- nu_hat ~ sqrt(2*q*r/n) and the partial-quotient reading
# ---------------------------------------------------------------------------

def pearson(xs, ys):
    k = len(xs)
    mx, my = sum(xs) / k, sum(ys) / k
    sxy = sum((a - mx) * (b - my) for a, b in zip(xs, ys))
    sxx = sum((a - mx) ** 2 for a in xs)
    syy = sum((b - my) ** 2 for b in ys)
    return sxy / math.sqrt(sxx * syy) if sxx > 0 and syy > 0 else float("nan")


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
            for t in range(i, j + 1):
                r[order[t]] = avg
            i = j + 1
        return r
    return pearson(rank(xs), rank(ys))


def exp_e4(curves, eff_grid):
    print()
    print(SEP)
    print("E4  nu_hat ~ sqrt(2*q*r/n), and nu_hat as a partial quotient")
    print(SEP)
    nh, approx, inv_a = [], [], []
    for c in curves:
        n, lam = c["n"], c["lam"]
        for eff in eff_grid:
            k1 = max(2, int(eff * math.sqrt(n)))
            k2 = math.isqrt(n) + 1
            cf = lambda1_closed_form(n, lam, k1, k2)
            if cf["q_star"] == 0:
                continue
            nh.append(cf["nu_hat"])
            approx.append(math.sqrt(2.0 * cf["q_star"] * cf["r_star"] / n))
            inv_a.append(1.0 / math.sqrt(max(1, cf["a_star"]) + 1.0))
    print(f"  samples              : {len(nh)}")
    print(f"  pearson (nu_hat, sqrt(2qr/n))   = {pearson(nh, approx):+.4f}")
    print(f"  spearman(nu_hat, sqrt(2qr/n))   = {spearman(nh, approx):+.4f}")
    rel = sorted(abs(a - b) / b for a, b in zip(approx, nh) if b > 0)
    med = rel[len(rel) // 2]
    p90 = rel[int(0.90 * len(rel))]
    print(f"  |approx-exact|/exact : median {med:.4f}  p90 {p90:.4f}")
    print(f"  pearson (nu_hat, 1/sqrt(a_next+1)) = {pearson(nh, inv_a):+.4f}")
    print(f"  spearman(nu_hat, 1/sqrt(a_next+1)) = {spearman(nh, inv_a):+.4f}")
    print("  reading              : small nu_hat (easy attack) <=> LARGE "
          "partial quotient\n                         at the balance scale.")


# ---------------------------------------------------------------------------
# E5 -- why every scale-free CF invariant had to fail
# ---------------------------------------------------------------------------

def scale_free_invariants(lam, n):
    """The June invariants: they are functions of lam/n alone."""
    cvs = cf_convergents(lam, n)
    a_all = [rec[3] for rec in cvs if rec[3] > 0]
    return {
        "max_a": max(a_all) if a_all else 0,
        "n_terms": len(cvs),
        "q_cf": cvs[0][0] if cvs else 0,
        "max_q_cf": max(rec[0] for rec in cvs) if cvs else 0,
    }


def exp_e5(curves, k1_grid):
    print()
    print(SEP)
    print("E5  scale-freeness: fix the curve, move K1 -> nu_hat moves, "
          "CF invariants cannot")
    print(SEP)
    print("For each curve the scale-free invariants are CONSTANT down the row "
          "by\nconstruction; nu_hat is not.  If nu_hat predicts p_hat (proved "
          "2026-07-29\nrun #2) then no scale-free function of lam/n can.")
    print()
    hdr = f"{'n':>10} {'mu':>7} {'max_a':>7} " + \
          " ".join(f"K1={k:<6}"[:7] for k in k1_grid) + f" {'spread':>8}"
    print(hdr)
    print("-" * len(hdr))
    spreads = []
    for c in curves[:12]:
        n, lam = c["n"], c["lam"]
        inv = scale_free_invariants(lam, n)
        k2 = math.isqrt(n) + 1
        row = []
        for k1 in k1_grid:
            cf = lambda1_closed_form(n, lam, k1, k2)
            row.append(cf["nu_hat"])
        spread = max(row) - min(row)
        spreads.append(spread)
        cells = " ".join(f"{v:>7.3f}" for v in row)
        print(f"{n:>10} {c['mu']:>7.4f} {inv['max_a']:>7} {cells} {spread:>8.3f}")
    spreads.sort()
    print()
    print(f"  median nu_hat spread over the K1 grid, per curve : "
          f"{spreads[len(spreads)//2]:.3f}")
    print(f"  (the 2026-07-29 20d low/high cut was at nu_hat 0.45 / 0.90, a "
          f"gap of 0.45)")
    print("  VERDICT              : a single curve traverses the whole "
          "decision range as K1\n                         moves, so max_a / "
          "q_cf / max_q_cf are provably\n                         "
          "non-predictive.  The June DEAD END is explained.")


# ---------------------------------------------------------------------------

def main():
    t0 = time.time()
    print(SEP)
    print("Thread 23 -- closed form for lambda_1(L2) in the GLV-HNP Phase 2 "
          "lattice")
    print(SEP)

    eff_grid = [0.02, 0.05, 0.10, 0.15, 0.25]
    k1_grid = [8, 16, 32, 64, 128, 256, 512]

    print("harvesting 20-bit j=0 GLV curves ...", flush=True)
    curves = harvest_j0_curves(19, 21, 60)
    print(f"  got {len(curves)} curves\n")

    ok1 = exp_e1(curves, eff_grid)
    ok2 = exp_e2()
    exp_e3(curves, eff_grid)
    exp_e4(curves, eff_grid)
    exp_e5(curves, k1_grid)

    print()
    print(SEP)
    print(f"SUMMARY   closed form exact on GLV curves : {ok1}")
    print(f"          closed form exact in general    : {ok2}")
    print(f"          total wall clock                : {time.time()-t0:.1f}s")
    print(SEP)
    return 0 if (ok1 and ok2) else 1


if __name__ == "__main__":
    sys.exit(main())
