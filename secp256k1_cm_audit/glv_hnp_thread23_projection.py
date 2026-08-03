"""
GLV-HNP Phase 2, Thread 23: does removing the trivial vector help?

Background (RESEARCH_AUTOLAB_LOG.md, 2026-07-29 run, exp T5):
  The Phase-2 Kannan lattice L (dim 2m+2, built by build_glv_lattice) always
  contains the "trivial" vector

        n * S_D * e_m      (norm n, since S_D = 1)

  because n*r_d - sum_i B_i*r_i kills every k1-column.  The planted vector has
  expected norm n*sqrt(2m/3 + 4/3) >= n*sqrt(2), so the planted vector is NEVER
  lambda_1.  T5 measured |sv[m]|/n = 1.0000 on every curve tested.  Recovery is
  therefore a BDD/coset condition, not an SVP condition, which retro-explains
  the six curve-level separators falsified 2026-06-21..06-29.

  The 2026-07-29 entry proposed Thread 23: reformulate so that the target IS
  lambda_1 -- "project the lattice along e_m ... or replace the Kannan
  embedding with an explicit CVP call".  Falsifier as stated there:

      if sv/pv rises above 1 after the reformulation AND the K1 wall on the
      lam*=0.07 curve (currently K1 ~ 4-6) moves outward, the reformulation is
      a real improvement; if the wall stays at K1 ~ 4-6, the wall is
      information-theoretic and Phase 2 is at its ceiling.

This script executes exactly that.

Four variants of the same instance, all sharing gen_signatures / scales /
build_glv_lattice verbatim from glv_hnp_phase2_lambda_threshold.py (which in
turn are verbatim from glv_hnp_phase2_20bit.py:262), so every number here is
directly comparable with the 2026-07-26 and 2026-07-29 tables.

  A  KANNAN-FULL   dim 2m+2.  Original.  LLL, scan rows for |last| = S_KANNAN.
  B  KANNAN-PROJ   dim 2m+1.  Column m (the d-column) deleted, so the trivial
                   vector n*S_D*e_m is quotiented out.  LLL on the 2m+2
                   generators (rank 2m+1), scan for |last| = S_KANNAN, then
                   recover d from the recovered (k1_i, k2_i) by
                   d = (k1_i + lam*k2_i - A_i) * B_i^-1 mod n.
  C  BDD-BABAI     dim 2m.  Kannan row dropped AND column m deleted.  Babai
                   nearest-plane (exact rational GSO) against the target
                   t = (-A_i*S_K1, 0).  Answer offset is (k1_i*S_K1, k2_i*S_K2).
  D  BDD-CVP       same lattice/target as C, but fplll's exact CVP enumeration.
                   D is an ORACLE: it returns the true closest vector, so a D
                   failure means the planted vector is not the closest vector
                   at all -- an instance-level wall no lattice reduction can
                   cross.  D succeeding where A fails means the wall was
                   reduction quality.

Every coordinate of every vector of the projected lattice is divisible by S_K1
in the k1-columns and by S_K2 in the k2-columns (all four generator families
have that property), so the (k1_i, k2_i) read-off in B/C/D is exact, never
rounded.

Experiments:
  E0  Gaussian-heuristic ratio ||target|| / GH(L) for A, B, C, and the
      information-theoretic bound, as a function of eff = K1*K2/n.  Predicts
      where each variant's wall should sit.
  E1  sv/pv in the projected lattice: is the target lambda_1 now?
  E2  the T4 K1-wall grid on the two historical 12-bit curves, all 4 variants.
      This is the falsifier.
  E3  17-bit fresh-curve sweep at eff = 0.05 / 0.15 / 0.25, variants A/B/D.

Run: python3 glv_hnp_thread23_projection.py
Runtime ~5 min (E3's exact CVP at dim 24 dominates).
"""

import math
import os
import random
import sys
from fractions import Fraction

from fpylll import IntegerMatrix, LLL, BKZ, CVP

HERE = os.path.dirname(os.path.abspath(__file__))

# ---------------------------------------------------------------------------
# Reuse the Thread-20 helpers verbatim (everything above its "EXPERIMENT T1"
# banner is pure definitions; the experiments below it run at import time, so
# we exec only the prefix).  This guarantees the lattice, the scaling and the
# signature generator are bit-identical to the 2026-07-29 run.
# ---------------------------------------------------------------------------
_SRC_PATH = os.path.join(HERE, "glv_hnp_phase2_lambda_threshold.py")
with open(_SRC_PATH) as _f:
    _src = _f.read()
_marker = "# EXPERIMENT T1"
assert _marker in _src, "helper source layout changed; update the marker"
_helpers = _src.split(_marker)[0]
exec(compile(_helpers, _SRC_PATH + " (helpers)", "exec"), globals())

# Names pulled in from the exec above (listed for the reader / linters):
#   modinv ec_add ec_mul tonelli_shanks find_generator eisenstein_decompose
#   j0_traces glv_roots lam_star gauss_reduce_2d lambda_block_mu
#   planted_norm_expected scales gen_signatures build_glv_lattice
#   planted_vector norm recover_d build_curve search_curves


# ===========================================================================
# Variant B — projected Kannan lattice (d-column deleted)
# ===========================================================================

def build_projected_kannan(sigs, n, lam, K1, K2):
    """The dim-(2m+2) Kannan lattice with column m removed.

    Returns a (2m+2) x (2m+1) generator matrix of rank 2m+1.  Column layout:
      0 .. m-1     k1 columns   (entries always divisible by S_K1)
      m .. 2m-1    k2 columns   (entries always divisible by S_K2)
      2m           Kannan column
    """
    m = len(sigs)
    M = build_glv_lattice(sigs, n, lam, K1, K2)
    return [row[:m] + row[m + 1:] for row in M]


def projected_planted(sigs, n, K1, K2):
    """pi(v_planted) = (k1_i*S_K1, k2_i*S_K2, S_KANNAN)."""
    m = len(sigs)
    S_K1, _S_D, S_K2, S_KAN = scales(n, K1, K2)
    v = [0] * (2 * m + 1)
    for i in range(m):
        v[i] = sigs[i]['k1'] * S_K1
        v[m + i] = sigs[i]['k2'] * S_K2
    v[2 * m] = S_KAN
    return v


def _d_from_k(k1s, k2s, sigs, n, lam):
    """Recover d from recovered nonce splits.  Requires agreement across all
    signatures (strictly stronger than the baseline recover_d, which only
    compares one candidate against the known d)."""
    d0 = None
    for i in range(len(sigs)):
        kf = (k1s[i] + lam * k2s[i]) % n
        B = sigs[i]['B'] % n
        if B == 0:
            continue
        d = (kf - sigs[i]['A']) * modinv(B, n) % n
        if d0 is None:
            d0 = d
        elif d != d0:
            return None
    return d0


def recover_d_projected(rows, sigs, n, lam, K1, K2, d_secret):
    """Scan a reduced basis of the projected Kannan lattice."""
    m = len(sigs)
    S_K1, _S_D, S_K2, S_KAN = scales(n, K1, K2)
    for row in rows:
        last = row[2 * m]
        if abs(last) != S_KAN:
            continue
        r = row if last > 0 else [-x for x in row]
        if any(r[i] % S_K1 for i in range(m)):
            continue
        if any(r[m + i] % S_K2 for i in range(m)):
            continue
        k1s = [r[i] // S_K1 for i in range(m)]
        k2s = [r[m + i] // S_K2 for i in range(m)]
        d = _d_from_k(k1s, k2s, sigs, n, lam)
        if d is not None and d == d_secret:
            return d
    return None


# ===========================================================================
# Variants C / D — BDD in the projected lattice (no Kannan row)
# ===========================================================================

def build_bdd_lattice(sigs, n, lam, K1, K2):
    """Generators of pi(L_0): the Phase-2 lattice minus its Kannan row, with
    the d-column deleted.  (2m+1) generators spanning a rank-2m lattice."""
    m = len(sigs)
    M = build_glv_lattice(sigs, n, lam, K1, K2)
    gens = [row[:m] + row[m + 1:2 * m + 1] for row in M[:2 * m + 1]]
    return gens


def bdd_target(sigs, n, K1, K2):
    """t = (-A_i*S_K1, 0, ..., 0); the planted offset is v - t."""
    m = len(sigs)
    S_K1, _S_D, _S_K2, _S_KAN = scales(n, K1, K2)
    t = [0] * (2 * m)
    for i in range(m):
        t[i] = -sigs[i]['A'] * S_K1
    return t


def recover_d_bdd(w, t, sigs, n, lam, K1, K2, d_secret):
    """w is the (claimed) closest lattice vector to t."""
    m = len(sigs)
    S_K1, _S_D, S_K2, _S_KAN = scales(n, K1, K2)
    off = [w[j] - t[j] for j in range(2 * m)]
    if any(off[i] % S_K1 for i in range(m)):
        return None
    if any(off[m + i] % S_K2 for i in range(m)):
        return None
    k1s = [off[i] // S_K1 for i in range(m)]
    k2s = [off[m + i] // S_K2 for i in range(m)]
    d = _d_from_k(k1s, k2s, sigs, n, lam)
    if d is not None and d == d_secret:
        return d
    return None


def lll_basis(gens):
    """LLL-reduce a (possibly rank-deficient) generator set and drop the zero
    rows fplll leaves at the top.  Returns a list of int rows."""
    A = IntegerMatrix.from_matrix(gens)
    LLL.reduction(A)
    ncols = A.ncols
    out = []
    for i in range(A.nrows):
        row = [A[i][j] for j in range(ncols)]
        if any(row):
            out.append(row)
    return out


def gso_exact(B):
    """Exact rational Gram-Schmidt.  dim <= ~26 here, so Fractions are cheap
    and immune to the f64 mantissa problems logged for the P-521 work."""
    k = len(B)
    Bs, mu, nrm = [], [[Fraction(0)] * k for _ in range(k)], []
    for i in range(k):
        v = [Fraction(x) for x in B[i]]
        for j in range(i):
            if nrm[j] == 0:
                continue
            num = sum(Fraction(B[i][t]) * Bs[j][t] for t in range(len(v)))
            mu[i][j] = num / nrm[j]
            v = [v[t] - mu[i][j] * Bs[j][t] for t in range(len(v))]
        Bs.append(v)
        nrm.append(sum(x * x for x in v))
    return Bs, nrm


def babai_nearest_plane(B, t):
    """Babai's nearest-plane on an already-reduced basis B, exact arithmetic."""
    Bs, nrm = gso_exact(B)
    b = [Fraction(x) for x in t]
    coeffs = [0] * len(B)
    for i in reversed(range(len(B))):
        if nrm[i] == 0:
            continue
        c = sum(b[j] * Bs[i][j] for j in range(len(b))) / nrm[i]
        ci = int(c + Fraction(1, 2)) if c >= 0 else -int(-c + Fraction(1, 2))
        coeffs[i] = ci
        b = [b[j] - ci * Fraction(B[i][j]) for j in range(len(b))]
    return [sum(coeffs[i] * B[i][j] for i in range(len(B))) for j in range(len(t))]


def exact_cvp(B, t):
    A = IntegerMatrix.from_matrix(B)
    return list(CVP.closest_vector(A, tuple(t)))


# ===========================================================================
# One instance, four variants
# ===========================================================================

def run_variants(curve, m, d_secret, K1, seed, variants=("A", "B", "C", "D"),
                 bkz_beta=None):
    """Returns {variant: bool}.  All four attack the SAME signature set."""
    p, b, n, lam, G = curve
    K2 = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, K1, K2, seed)
    if len(sigs) < m:
        return None
    S_K1, S_D, S_K2, S_KAN = scales(n, K1, K2)
    res = {}

    if "A" in variants:
        M = build_glv_lattice(sigs, n, lam, K1, K2)
        A = IntegerMatrix.from_matrix(M)
        if bkz_beta:
            BKZ.reduction(A, BKZ.Param(bkz_beta))
        else:
            LLL.reduction(A)
        dim = 2 * m + 2
        rows = [[A[i][j] for j in range(dim)] for i in range(dim)]
        res["A"] = recover_d(rows, m, n, S_KAN, d_secret) is not None

    if "B" in variants:
        gens = build_projected_kannan(sigs, n, lam, K1, K2)
        Ab = IntegerMatrix.from_matrix(gens)
        if bkz_beta:
            BKZ.reduction(Ab, BKZ.Param(bkz_beta))
        else:
            LLL.reduction(Ab)
        rows = [[Ab[i][j] for j in range(2 * m + 1)] for i in range(Ab.nrows)]
        res["B"] = recover_d_projected(rows, sigs, n, lam, K1, K2, d_secret) is not None

    if "C" in variants or "D" in variants:
        gens = build_bdd_lattice(sigs, n, lam, K1, K2)
        Bmat = lll_basis(gens)
        t = bdd_target(sigs, n, K1, K2)
        if "C" in variants:
            w = babai_nearest_plane(Bmat, t)
            res["C"] = recover_d_bdd(w, t, sigs, n, lam, K1, K2, d_secret) is not None
        if "D" in variants:
            try:
                w = exact_cvp(Bmat, t)
                res["D"] = recover_d_bdd(w, t, sigs, n, lam, K1, K2, d_secret) is not None
            except Exception as exc:            # enumeration failure / overflow
                res["D"] = None
                res["D_err"] = str(exc)[:60]
    return res


def rate(curve, m, K1, seeds, variants=("A", "B", "C", "D"), bkz_beta=None):
    acc = {v: 0 for v in variants}
    tot = 0
    p, b, n, lam, G = curve
    for s in seeds:
        d_secret = random.Random(s ^ 0xC0FFEE).randint(1, n - 1)
        r = run_variants(curve, m, d_secret, K1, s, variants, bkz_beta)
        if r is None:
            continue
        tot += 1
        for v in variants:
            if r.get(v):
                acc[v] += 1
    return acc, tot


# ===========================================================================
print("=" * 78)
print("Thread 23 — quotient out the trivial vector: does the wall move?")
print("=" * 78)

SEEDS = [42, 1234, 9999, 555, 31337]

HIST = [
    # label,             p,    b, n,    lam,  K1, m
    ("12-bit/2557",      2557, 2, 2659, 1755, 8,  8),
    ("12-bit/2677 FAIL", 2677, 2, 2647, 185,  8,  10),
]
hist = []
for label, p, b, n, lam, k1, m in HIST:
    G = find_generator(p, b, n)
    assert G is not None, f"no generator for {label}"
    assert (lam * lam + lam + 1) % n == 0
    hist.append((label, (p, b, n, lam, G), k1, m))


# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("E0: Gaussian-heuristic ratio ||target|| / GH  for each variant")
print("-" * 78)
print("A  dim 2m+2  det = (n*S_K1)^m * S_D * S_K2^m * S_KAN     |t|^2 = n^2(2m/3+4/3)")
print("B  dim 2m+1  det = det(A)/(n*S_D)                        |t|^2 = n^2(2m/3+1)")
print("C  dim 2m    det = (n*S_K1)^m * S_K2^m / n               |t|^2 = n^2(2m/3)")
print("IT bound: (m-1)*log n >= m*log(K1*K2), i.e. eff <= n^(-1/m).")
print("Ratio < 1 => the target is (heuristically) the shortest / closest vector.\n")


def gh_ratios(n, m, K1):
    K2 = math.isqrt(n) + 1
    S_K1, S_D, S_K2, S_KAN = scales(n, K1, K2)
    eff = K1 * K2 / n
    logs = math.log

    def gh(dim, logdet):
        return math.sqrt(dim / (2 * math.pi * math.e)) * math.exp(logdet / dim)

    logdet_A = m * logs(n * S_K1) + logs(S_D) + m * logs(S_K2) + logs(S_KAN)
    logdet_B = logdet_A - logs(n * S_D)
    logdet_C = m * logs(n * S_K1) + m * logs(S_K2) - logs(n)

    tA = n * math.sqrt(2 * m / 3 + 4 / 3)
    tB = n * math.sqrt(2 * m / 3 + 1)
    tC = n * math.sqrt(2 * m / 3)
    return (eff,
            tA / gh(2 * m + 2, logdet_A),
            tB / gh(2 * m + 1, logdet_B),
            tC / gh(2 * m, logdet_C),
            eff / n ** (-1.0 / m))


print(f"{'n bits':>7} {'m':>3} {'K1':>5} {'eff':>7} {'A':>7} {'B':>7} {'C':>7} {'eff/IT':>7}")
for nb, m in ((12, 10), (17, 10), (17, 16), (20, 12)):
    n = 2 ** nb - 1
    for K1 in (2, 4, 8, 16, 32, 64):
        eff, ra, rb, rc, rit = gh_ratios(n, m, K1)
        if eff > 1.2:
            continue
        print(f"{nb:>7} {m:>3} {K1:>5} {eff:>7.3f} {ra:>7.3f} {rb:>7.3f} "
              f"{rc:>7.3f} {rit:>7.3f}")

print("\nReading: B improves on A by ~4%, C on A by ~11%, in the ratio.  Since the")
print("ratio scales as (K1*K2)^-1/2 near the wall, an 11% ratio gain predicts the")
print("K1 wall moving out by ~1.11^2 = 1.23x -- i.e. K1 ~ 4-6 becomes K1 ~ 5-7.")
print("A large jump (2x+) would falsify the 'wall is information-theoretic' reading;")
print("no movement at all would confirm it.")


# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("E1: is the target lambda_1 after the projection?")
print("-" * 78)
print("sv/pv = ||shortest reduced row|| / ||planted||.  T5 (2026-07-29) measured")
print("0.42-0.60 in variant A, with 100% of the shortest vector's energy in the")
print("d-column.  If the projection works, B's sv/pv should be ~1.\n")

print(f"{'curve':<18} {'K1':>3} {'m':>3} {'A sv/pv':>8} {'A |sv_m|/n':>11} "
      f"{'B sv/pv':>8} {'B is planted':>13}")
for label, curve, K1, m in hist:
    p, b, n, lam, G = curve
    K2 = math.isqrt(n) + 1
    d_secret = random.Random(42 ^ 0xC0FFEE).randint(1, n - 1)
    sigs = gen_signatures(G, d_secret, m, n, lam, p, K1, K2, 42)
    S_K1, S_D, S_K2, S_KAN = scales(n, K1, K2)

    M = build_glv_lattice(sigs, n, lam, K1, K2)
    A = IntegerMatrix.from_matrix(M)
    LLL.reduction(A)
    dim = 2 * m + 2
    rowsA = [[A[i][j] for j in range(dim)] for i in range(dim)]
    pv = norm(planted_vector(sigs, d_secret, n, K1, K2))
    svA = min((norm(r) for r in rowsA if any(r)), default=0.0)
    svA_row = min((r for r in rowsA if any(r)), key=norm)
    frac_m = abs(svA_row[m]) / n

    gens = build_projected_kannan(sigs, n, lam, K1, K2)
    Ab = IntegerMatrix.from_matrix(gens)
    LLL.reduction(Ab)
    rowsB = [[Ab[i][j] for j in range(2 * m + 1)] for i in range(Ab.nrows)]
    pvB = norm(projected_planted(sigs, n, K1, K2))
    svB_row = min((r for r in rowsB if any(r)), key=norm)
    svB = norm(svB_row)
    is_planted = (svB_row == projected_planted(sigs, n, K1, K2) or
                  [-x for x in svB_row] == projected_planted(sigs, n, K1, K2))

    print(f"{label:<18} {K1:>3} {m:>3} {svA/pv:>8.3f} {frac_m:>11.4f} "
          f"{svB/pvB:>8.3f} {str(is_planted):>13}")


# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("E2: the K1 wall, all four variants (THE FALSIFIER)")
print("-" * 78)
print("Replicates the 2026-07-29 T4 grid.  T4 found the wall at K1~12-16 for")
print("lam*=0.34 and K1~4-6 for lam*=0.07, both in variant A.\n")

K1_GRID = [2, 3, 4, 6, 8, 12, 16, 24, 32]
for label, curve, _K1hist, m in hist:
    p, b, n, lam, G = curve
    K2 = math.isqrt(n) + 1
    print(f"\n{label}: p={p} n={n} lam={lam} lam*={lam_star(lam,n):.4f} m={m} "
          f"K2={K2}")
    print(f"{'K1':>4} {'eff':>7} | {'A':>5} {'B':>5} {'C':>5} {'D':>5}   "
          f"(x/{len(SEEDS)} seeds)")
    for K1 in K1_GRID:
        acc, tot = rate(curve, m, K1, SEEDS)
        eff = K1 * K2 / n
        cells = []
        for v in ("A", "B", "C", "D"):
            cells.append(f"{acc[v]}/{tot}")
        print(f"{K1:>4} {eff:>7.3f} | " + " ".join(f"{c:>5}" for c in cells))
    sys.stdout.flush()


# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("E3: 17-bit fresh curves, eff = 0.05 / 0.15 / 0.25")
print("-" * 78)
print("2026-07-29 T3 (variant A): 19/20, 3/20, 0/20 curves recovering.\n")

curves17 = search_curves(2 ** 16, 2 ** 17, per_bin=1, nbins=8)
print(f"found {len(curves17)} 17-bit j=0 GLV curves\n")

M_SIGS = 10
SEEDS3 = [42, 1234, 9999]
for eff_t in (0.05, 0.15, 0.25):
    tallies = {v: 0 for v in ("A", "B", "C", "D")}
    trials = 0
    per_curve = []
    for (p, b, n, lam, G) in curves17:
        K2 = math.isqrt(n) + 1
        K1 = max(2, int(eff_t * n / K2))
        acc, tot = rate((p, b, n, lam, G), M_SIGS, K1, SEEDS3)
        trials += tot
        for v in tallies:
            tallies[v] += acc[v]
        per_curve.append((n, lam_star(lam, n), K1 * K2 / n, acc, tot))
    print(f"eff_target={eff_t}:  trials={trials}")
    print(f"{'n':>8} {'lam*':>7} {'eff':>7} | {'A':>5} {'B':>5} {'C':>5} {'D':>5}")
    for n, ls, eff, acc, tot in per_curve:
        print(f"{n:>8} {ls:>7.4f} {eff:>7.3f} | " +
              " ".join(f"{acc[v]}/{tot:<3}" for v in ("A", "B", "C", "D")))
    print("  TOTAL: " + "  ".join(f"{v}={tallies[v]}/{trials}"
                                  for v in ("A", "B", "C", "D")))
    sys.stdout.flush()

print("\n" + "=" * 78)
print("done")
print("=" * 78)
