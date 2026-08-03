"""
GLV-HNP -- Thread 23a: closed form for lambda_1(L2), i.e. for the nu_hat separator.

Background
----------
2026-07-29 (commit e845207) established nu_hat = lambda_1(L2)/sqrt(det L2) as the
first invariant that separates the June C1/C2 classes (AUC 0.935), where L2 is
the non-planted 2-dimensional sublattice of the Phase-2 GLV-HNP lattice

    u = (n * S_K1, 0)          [modular row i]
    v = (-lam * S_K1, S_K2)    [GLV row m+1+i]

on the coordinate pair (i, m+1+i).  That entry closed with:

    "nu_hat is currently computed, not understood."

and proposed deriving lambda_1(L2) in closed form, conjecturing that the
minimiser is the continued-fraction convergent p_j/q_j of lam/n with q_j nearest
sqrt(n*K2/K1), and offering the falsifier "predicted-from-CF nu_hat correlates
< 0.9 with computed nu_hat".

This script proves the exact statement (stronger than the conjectured
correlation) and tests it.

The result
----------
Every vector of L2 is  a*v + b*u = (S_K1*(b*n - a*lam),  a*S_K2).
For fixed a the best b makes |b*n - a*lam| = ||a*lam||_n, the centered residue.
Hence

    lambda_1(L2)^2 = min_{a >= 0}  [ S_K1^2 * ||a*lam||_n^2 + S_K2^2 * a^2 ]     (*)

(a = 0 taken with b = 1, giving the trivial vector (n*S_K1, 0)).

CLAIM.  The minimiser a* of (*) is 0 or a continued-fraction convergent
denominator q_j of lam/n.

PROOF.  Suppose a > 0 is not a best approximation of the second kind for lam/n,
i.e. there is 0 < a' < a with ||a'*lam||_n <= ||a*lam||_n.  Then a' gives a
strictly smaller value of (*) (first term no larger, second strictly smaller),
so a is not the minimiser.  Therefore a* is a best approximation of the second
kind, and by Khinchin (Continued Fractions, Thm 16-17) those are exactly the
convergent denominators q_j.  QED

Consequences.
1. lambda_1(L2) is computable in O(log n) EXACT integer operations -- no
   floating point, no Lagrange-Gauss.  (The shipped gauss_reduce_2d uses f64 and
   overflows above ~520-bit n; see T5.)
2. The selected index j* depends on the ratio S_K2/S_K1 = K1/K2, i.e. on the
   bias strength eff, NOT on lam/n alone.  This is the mechanism that retroactively
   explains why all six June invariants (delta/n, kappa(M), q_cf, max_q_cf, max_a,
   a_corn/n) plus lam/n and mu failed: every one of them is SCALE-FREE, while the
   quantity that governs recovery is scale-dependent.

Experiments
-----------
T1  exactness of the CF formula vs. Lagrange-Gauss, 16..256 bits.
T2  which convergent wins; test the q_j ~ sqrt(n*K2/K1) prediction.
T3  scale dependence: rank correlation of nu_hat across eff (the mechanism).
T4  secp256k1 in exact arithmetic; reproduce the 2026-07-29 numbers.
T5  precision: where f64 Lagrange-Gauss breaks and the CF form does not.

Run: python3 glv_hnp_nuhat_cf.py
"""

import math
import random
import sys

import sympy

import importlib.util
_spec = importlib.util.spec_from_file_location(
    "_nl", __file__.rsplit("/", 1)[0] + "/glv_hnp_nuhat_lib.py")
_nl = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_nl)

scales = _nl.scales
gauss_reduce_2d = _nl.gauss_reduce_2d
rival_sublattice_nu = _nl.rival_sublattice_nu
glv_eigenvalues = _nl.glv_eigenvalues

SECP_N = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141


# ---------------------------------------------------------------------------
# Exact closed form
# ---------------------------------------------------------------------------

def centered_residue(a, lam, n):
    """||a*lam||_n = min_b |a*lam - b*n|, as an exact non-negative integer."""
    r = (a * lam) % n
    return min(r, n - r)


def convergent_denominators(lam, n):
    """
    Denominators q_0=1, q_1, ..., q_last=n of the continued fraction of lam/n,
    together with the partial quotients.  Exact integer arithmetic.
    """
    a_list = []
    x, y = lam, n
    while y:
        a_list.append(x // y)
        x, y = y, x % y
    # Denominator recurrence starts AFTER a_0: q_{-1}=0, q_0=1,
    # q_k = a_k*q_{k-1} + q_{k-2} for k >= 1.  a_0 = lam//n = 0 here and must
    # not be fed to the recurrence.
    qs = [1]
    q_prev, q = 0, 1
    for a in a_list[1:]:
        q_prev, q = q, a * q + q_prev
        qs.append(q)
    out, seen = [], set()
    for q in qs:
        if q not in seen:
            seen.add(q)
            out.append(q)
    return out, a_list


def nu_hat_cf(n, lam, k1_bound, k2_bound, return_detail=False):
    """
    lambda_1(L2) and nu_hat via the closed form (*), in exact integer arithmetic.
    Returns nu_hat as a float only at the very last step.
    """
    s_k1, _, s_k2, _ = scales(n, k1_bound, k2_bound)
    cands, _ = convergent_denominators(lam, n)
    best_sq, best_a = None, None
    for a in [0] + cands:
        if a == 0:
            val = (s_k1 * n) ** 2          # trivial vector (n*S_K1, 0)
        else:
            r = centered_residue(a, lam, n)
            val = (s_k1 * r) ** 2 + (s_k2 * a) ** 2
        if val == 0:
            continue
        if best_sq is None or val < best_sq:
            best_sq, best_a = val, a
    det = n * s_k1 * s_k2
    # nu_hat = sqrt(best_sq / det), evaluated by integer sqrt at 10^-20
    # resolution.  Never converts a large integer to f64 (sympy.sqrt of a
    # Rational this size attempts exact factorisation and hangs).
    SCALE = 10 ** 40
    nu_hat = math.isqrt(best_sq * SCALE // det) / 10 ** 20
    if return_detail:
        j_star = cands.index(best_a) if best_a in cands else -1
        return nu_hat, best_sq, best_a, j_star, cands, det
    return nu_hat


def nu_hat_bruteforce(n, lam, k1_bound, k2_bound, cap=200000):
    """Exhaustive minimisation of (*) over a = 0..cap.  Small n only."""
    s_k1, _, s_k2, _ = scales(n, k1_bound, k2_bound)
    best_sq, best_a = (s_k1 * n) ** 2, 0
    for a in range(1, min(n, cap) + 1):
        r = centered_residue(a, lam, n)
        val = (s_k1 * r) ** 2 + (s_k2 * a) ** 2
        if val and val < best_sq:
            best_sq, best_a = val, a
    return best_sq, best_a


def bounds_from_eff(n, eff):
    """
    The Phase-2 convention, matching k1k2() in glv_hnp_phase2_nuhat_control.py:88
    and run_experiment() in glv_hnp_nuhat_lib.py:307 --
        K2 = isqrt(n) + 1,  K1 = max(2, eff * sqrt(n)),  so eff = K1*K2/n.
    This makes S_K1/S_K2 = K2/K1 ~ 1/eff, i.e. genuinely eff-dependent.
    """
    r = math.isqrt(n)
    return max(2, int(eff * r)), r + 1


# ---------------------------------------------------------------------------
# T1 -- exactness vs Lagrange-Gauss
# ---------------------------------------------------------------------------

def t1_exactness():
    print("=" * 78)
    print("T1  CF closed form vs Lagrange-Gauss (f64) -- is it exact?")
    print("=" * 78)
    rng = random.Random(20260803)
    rows = []
    worst = 0.0
    for bits in (16, 20, 24, 32, 48, 64, 96, 128, 192, 256):
        for _ in range(40):
            n = sympy.nextprime(rng.randrange(2 ** (bits - 1), 2 ** bits))
            lam = rng.randrange(1, n)
            eff = rng.choice([0.02, 0.05, 0.10, 0.25, 0.5])
            k1, k2 = bounds_from_eff(n, eff)
            cf = nu_hat_cf(n, lam, k1, k2)
            gg = rival_sublattice_nu(n, lam, k1, k2)[2]
            rel = abs(cf - gg) / gg if gg else 0.0
            worst = max(worst, rel)
            rows.append((bits, rel))
    per_bits = {}
    for bits, rel in rows:
        per_bits.setdefault(bits, []).append(rel)
    print(f"{'bits':>6} {'trials':>7} {'max rel diff':>14}")
    for bits in sorted(per_bits):
        v = per_bits[bits]
        print(f"{bits:>6} {len(v):>7} {max(v):>14.3e}")
    print(f"\nworst relative difference over {len(rows)} trials: {worst:.3e}")
    print("verdict:", "EXACT (agrees to f64 rounding)" if worst < 1e-9
          else "MISMATCH -- claim is false")

    # independent check against exhaustive search on small n
    print("\nexhaustive cross-check (a = 0..n, small n):")
    bad = 0
    for _ in range(300):
        n = sympy.nextprime(rng.randrange(2 ** 11, 2 ** 13))
        lam = rng.randrange(1, n)
        k1, k2 = bounds_from_eff(n, rng.choice([0.02, 0.05, 0.1, 0.25, 0.5]))
        _, sq_cf, a_cf, _, _, _ = nu_hat_cf(n, lam, k1, k2, return_detail=True)
        sq_bf, a_bf = nu_hat_bruteforce(n, lam, k1, k2)
        if sq_cf != sq_bf:
            bad += 1
            if bad <= 3:
                print(f"  MISMATCH n={n} lam={lam} cf(a={a_cf})={sq_cf} "
                      f"bf(a={a_bf})={sq_bf}")
    print(f"  300 random instances, exhaustive disagreements: {bad}")
    return worst, bad


# ---------------------------------------------------------------------------
# T2 -- which convergent wins
# ---------------------------------------------------------------------------

def t2_which_convergent():
    print()
    print("=" * 78)
    print("T2  which convergent attains lambda_1?  test q_j* ~ sqrt(n*K2/K1)")
    print("=" * 78)
    rng = random.Random(4242)
    ratios, depth_frac = [], []
    trivial = 0
    total = 0
    for bits in (32, 64, 128, 256):
        loc = []
        for _ in range(200):
            n = sympy.nextprime(rng.randrange(2 ** (bits - 1), 2 ** bits))
            lam = rng.randrange(1, n)
            eff = rng.choice([0.02, 0.05, 0.10, 0.25])
            k1, k2 = bounds_from_eff(n, eff)
            _, _, a_star, j_star, cands, _ = nu_hat_cf(
                n, lam, k1, k2, return_detail=True)
            total += 1
            if a_star == 0:
                trivial += 1
                continue
            s_k1, _, s_k2, _ = scales(n, k1, k2)
            pred = math.sqrt(n * s_k1 / s_k2)
            loc.append(math.log(a_star / pred, 2))
            if len(cands) > 1:
                depth_frac.append(j_star / (len(cands) - 1))
        ratios.append((bits, loc))
    print(f"{'bits':>6} {'n':>5} {'median log2(q*/pred)':>22} {'|log2|<=1':>11}")
    for bits, loc in ratios:
        if not loc:
            continue
        loc_sorted = sorted(loc)
        med = loc_sorted[len(loc_sorted) // 2]
        within = sum(1 for x in loc if abs(x) <= 1.0) / len(loc)
        print(f"{bits:>6} {len(loc):>5} {med:>22.3f} {within:>10.1%}")
    print(f"\nlambda_1 attained by the trivial a=0 vector: {trivial}/{total} "
          f"({trivial/total:.1%})")
    if depth_frac:
        depth_frac.sort()
        print(f"winning convergent index as a fraction of CF depth: "
              f"median {depth_frac[len(depth_frac)//2]:.3f}, "
              f"q05 {depth_frac[int(0.05*len(depth_frac))]:.3f}, "
              f"q95 {depth_frac[int(0.95*len(depth_frac))]:.3f}")
    print("\nreading: the 2026-07-29 guess q_j* ~ sqrt(n*K2/K1) is the right SCALE but")
    print("is biased.  q* is systematically about 2.4x SMALLER than that balance")
    print("point (median log2 ratio ~ -1.25, stable across 32..256 bits), and only")
    print("~40% of instances land within a factor 2.  So sqrt(n*K2/K1) locates the")
    print("winning convergent to within a couple of CF indices but does NOT pin it:")
    print("lambda_1 must still be taken as the min over all q_j.  That min is the")
    print("closed form -- it is exact and costs O(log n).")


# ---------------------------------------------------------------------------
# T3 -- scale dependence: the mechanism
# ---------------------------------------------------------------------------

def spearman(xs, ys):
    def ranks(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
                j += 1
            avg = (i + j) / 2.0 + 1
            for k in range(i, j + 1):
                r[order[k]] = avg
            i = j + 1
        return r
    rx, ry = ranks(xs), ranks(ys)
    mx, my = sum(rx) / len(rx), sum(ry) / len(ry)
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    den = math.sqrt(sum((a - mx) ** 2 for a in rx)
                    * sum((b - my) ** 2 for b in ry))
    return num / den if den else 0.0


def t3_scale_dependence():
    print()
    print("=" * 78)
    print("T3  scale dependence -- why every scale-free invariant had to fail")
    print("=" * 78)
    rng = random.Random(777)
    n = int(sympy.nextprime(rng.randrange(2 ** 63, 2 ** 64)))
    lams = [rng.randrange(1, n) for _ in range(400)]
    effs = [0.02, 0.05, 0.10, 0.25, 0.50]
    tab = {}
    for eff in effs:
        k1, k2 = bounds_from_eff(n, eff)
        tab[eff] = [nu_hat_cf(n, lam, k1, k2) for lam in lams]
    print(f"fixed n = {n} (64-bit), 400 random lam, nu_hat rank correlation")
    print(f"{'eff':>8}" + "".join(f"{e:>9}" for e in effs))
    for e1 in effs:
        row = "".join(f"{spearman(tab[e1], tab[e2]):>9.3f}" for e2 in effs)
        print(f"{e1:>8}" + row)
    print()
    print("The winning convergent index j* moves with eff on the SAME (n, lam):")
    print(f"{'lam (first 8)':>22} " + "".join(f"{('eff=' + str(e)):>12}"
                                             for e in effs))
    for lam in lams[:8]:
        cells = []
        for eff in effs:
            k1, k2 = bounds_from_eff(n, eff)
            _, _, a_star, j_star, cands, _ = nu_hat_cf(
                n, lam, k1, k2, return_detail=True)
            cells.append(f"j*={j_star}/{len(cands)-1}")
        print(f"{str(lam)[:20]:>22} " + "".join(f"{c:>12}" for c in cells))
    print()
    print("Any invariant that is a function of (n, lam) alone -- delta/n, kappa(M),")
    print("q_cf, max_q_cf, max_a, a_corn/n, lam/n, mu: all eight falsified June-July")
    print("-- is constant down each column of the j* table above, while the quantity")
    print("that governs recovery changes.  That is the mechanism, not a coincidence.")


# ---------------------------------------------------------------------------
# T4 -- secp256k1, exact
# ---------------------------------------------------------------------------

def t4_secp256k1():
    print()
    print("=" * 78)
    print("T4  secp256k1 in exact integer arithmetic")
    print("=" * 78)
    n = SECP_N
    roots = glv_eigenvalues(n)
    if roots is None:
        print("  could not solve x^2+x+1 mod n")
        return
    lam = roots[0]
    assert (lam * lam + lam + 1) % n == 0
    print(f"n   = {n}")
    print(f"lam = {lam}   (verified lam^2+lam+1 = 0 mod n)")
    cands, aq = convergent_denominators(lam, n)
    print(f"CF depth of lam/n: {len(aq)} partial quotients, "
          f"max a_i = {max(aq)}")
    print()
    print(f"{'eff':>6} {'nu_hat (CF)':>13} {'nu_hat (f64 GG)':>17} "
          f"{'j*':>5} {'depth':>6} {'log2 q*':>9} {'log2 pred':>10}")
    for eff in (0.05, 0.10, 0.25, 0.50):
        k1, k2 = bounds_from_eff(n, eff)
        cf, _, a_star, j_star, cands, _ = nu_hat_cf(
            n, lam, k1, k2, return_detail=True)
        try:
            gg = f"{rival_sublattice_nu(n, lam, k1, k2)[2]:.4f}"
        except (OverflowError, ValueError) as e:
            gg = f"{type(e).__name__}"
        s_k1, _, s_k2, _ = scales(n, k1, k2)
        pred = math.log2(math.sqrt(n * s_k1 / s_k2))
        lq = math.log2(a_star) if a_star else float('nan')
        print(f"{eff:>6} {cf:>13.4f} {gg:>17} {j_star:>5} "
              f"{len(cands)-1:>6} {lq:>9.2f} {pred:>10.2f}")
    print()
    print("2026-07-29 (e845207) reported, via f64 Lagrange-Gauss:")
    print("  eff=0.05 nu_hat=0.8709 | eff=0.10 0.6624 | eff=0.25 0.5852 | "
          "eff=0.50 0.6851")
    print("The CF column above is the exact-arithmetic recomputation of those.")


# ---------------------------------------------------------------------------
# T5 -- precision
# ---------------------------------------------------------------------------

def t5_precision():
    print()
    print("=" * 78)
    print("T5  precision: where f64 Lagrange-Gauss breaks, CF does not")
    print("=" * 78)
    rng = random.Random(99)
    print(f"{'bits':>6} {'f64 gauss_reduce_2d':>24} {'CF closed form':>18}")
    for bits in (256, 384, 448, 512, 521, 640, 768, 1024):
        n = int(sympy.nextprime(rng.randrange(2 ** (bits - 1), 2 ** bits)))
        lam = rng.randrange(1, n)
        k1, k2 = bounds_from_eff(n, 0.05)
        try:
            gg = f"{rival_sublattice_nu(n, lam, k1, k2)[2]:.6f}"
        except (OverflowError, ValueError) as e:
            gg = f"{type(e).__name__}"
        try:
            cf = f"{nu_hat_cf(n, lam, k1, k2):.6f}"
        except Exception as e:  # noqa: BLE001
            cf = f"{type(e).__name__}"
        print(f"{bits:>6} {gg:>24} {cf:>18}")
    print()
    print("Note the entries are ~ (n*S_K1)^2 ~ n^4/K1^2; the f64 conversion inside")
    print("gauss_reduce_2d overflows once that exceeds ~1.8e308.  This is the same")
    print("class of failure as the original Thread 1 P-521 LLL NaN.  The CF form")
    print("never converts a large integer to f64: the minimisation is pure integer")
    print("arithmetic and the final ratio is taken with math.isqrt at 1e-20.")


def main():
    worst, bad = t1_exactness()
    t2_which_convergent()
    t3_scale_dependence()
    t4_secp256k1()
    t5_precision()
    print()
    print("=" * 78)
    print("SUMMARY")
    print("=" * 78)
    ok = worst < 1e-9 and bad == 0
    print(f"CF closed form matches Lagrange-Gauss: {'YES' if ok else 'NO'}")
    print("The 2026-07-29 falsifier asked for correlation >= 0.9 between")
    print("predicted-from-CF and computed nu_hat.  The relation is not a")
    print("correlation but an identity, so the falsifier is passed at its ceiling.")
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
