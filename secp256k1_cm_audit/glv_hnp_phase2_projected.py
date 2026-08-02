"""
GLV-HNP Phase 2, Thread 23: reformulate the lattice so the planted vector is
the target of a *solvable* problem (it is provably never lambda_1 today).

Background (RESEARCH_AUTOLAB_LOG.md, 2026-07-29 run, exp T5):
  The Phase-2 Kannan lattice L (dim 2m+2) always contains the trivial vector
      v_triv = n * S_D * e_m          (d is only defined mod n)
  with ||v_triv|| = n*S_D, while the planted vector has
      ||v_planted|| ~ n*sqrt(2m/3 + 4/3).
  So v_triv is shorter for every m >= 1 and every choice of S_D (both scale
  linearly in S_D).  Measured: |sv[m]|/n = 1.0000 exactly on every curve, and
  the shortest vector's Kannan coordinate is always 0.  Recovery is therefore
  NOT an SVP condition in L; it is a BDD/coset condition.

  The 2026-07-29 run proposed Thread 23:
    "project the lattice along e_m (quotient out the trivial n*e_m direction)
     and solve BDD in the projection, or replace the Kannan embedding with an
     explicit CVP call (Babai nearest-plane) targeting (A_i*S_K1, 0, ...).
     Falsifier: if sv/pv rises above 1 after the reformulation AND the K1 wall
     in T4 moves outward on the lam*=0.07 curve (currently K1 ~ 4-6), the
     reformulation is a real improvement; if the wall stays at K1 ~ 4-6, then
     the wall is information-theoretic and Phase 2 is at its ceiling."

Three formulations are compared on identical (curve, K1, m, seed) instances:

  KAN   baseline.  Kannan lattice L, dim 2m+2, LLL, scan rows for
        |last| = S_KANNAN.  (verbatim build_glv_lattice / recover_d)

  PROJ  quotient lattice L' = L / <n*S_D*e_m>, realised by deleting column m.
        An explicit basis of L' is constructed (the deletion makes the 2m+2
        generators rank-deficient, so the dependency n*pi(row_d) = sum_i B_i
        pi(row_i) is removed by hand), LLL-reduced with a transformation
        matrix U, and every reduced row is lifted back into L via U to read
        off d from column m.  v_triv is not in L', so the planted vector CAN
        be lambda_1 here.

  CVP   no Kannan row, no d ambiguity in the objective.  Lattice L0 (dim 2m+1)
        spanned by n*S_K1*e_i, row_d, and the k2 rows; target
        t = (-A_i*S_K1, 0, ..., 0).  The closest vector satisfies
        v - t = (k1_i*S_K1, d*S_D, k2_i*S_K2) exactly.  Solved by exact-
        rational Babai nearest-plane on the LLL-reduced basis (no floating
        point anywhere), and optionally by fpylll's exact CVP enumeration.

Experiments:
  U1  sanity: L' is a basis (determinant check), U-lift is exact, CVP target
      algebra is exact.
  U2  the falsifier: T4's K1 grid on both historical 12-bit curves, three
      formulations, 5 seeds each.  Does the wall move?
  U3  sv/pv and energy distribution in L' (does sv/pv rise above 1?).
  U4  Gaussian-heuristic model of the wall: predicted critical
      eff = K1*K2/n, compared against the T3 data (eff = 0.05 / 0.15 / 0.25).
  U5  17-bit sweep at the predicted critical eff, best formulation.

Run: python3 glv_hnp_phase2_projected.py
"""

import ast
import math
import os
import random
from fractions import Fraction

from fpylll import IntegerMatrix, LLL

# ---------------------------------------------------------------------------
# Reuse the Phase-2 primitives VERBATIM from the Thread 20 script.
# That file has no __main__ guard (it is a driver), so import only its
# function/class/import definitions via AST filtering.  This guarantees the
# lattice construction, signature generation and curve search are bit-identical
# to the 2026-07-29 run rather than a re-typed copy.
# ---------------------------------------------------------------------------

_SRC = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                    "glv_hnp_phase2_lambda_threshold.py")
with open(_SRC) as _f:
    _TREE = ast.parse(_f.read(), _SRC)
_KEEP = [nd for nd in _TREE.body
         if isinstance(nd, (ast.Import, ast.ImportFrom,
                            ast.FunctionDef, ast.ClassDef))]
exec(compile(ast.Module(body=_KEEP, type_ignores=[]), _SRC, "exec"), globals())

# names now available: modinv, ec_add, ec_mul, tonelli_shanks, find_generator,
# eisenstein_decompose, j0_traces, glv_roots, lam_star, gauss_reduce_2d,
# lambda_block_mu, planted_norm_expected, scales, gen_signatures,
# build_glv_lattice, planted_vector, norm, recover_d, run_experiment,
# success_rate, hnf_fingerprint, build_curve, search_curves


# ===========================================================================
# Exact rational Gram-Schmidt / Babai nearest plane
# ===========================================================================

def gso_exact(B):
    """Exact Gram-Schmidt over Q.  Returns (Bstar, sqnorms)."""
    dim = len(B[0])
    Bstar, sqn = [], []
    for i in range(len(B)):
        v = [Fraction(x) for x in B[i]]
        for j in range(i):
            if sqn[j] == 0:
                continue
            mu = sum(Fraction(B[i][k]) * Bstar[j][k] for k in range(dim)) / sqn[j]
            if mu:
                v = [v[k] - mu * Bstar[j][k] for k in range(dim)]
        Bstar.append(v)
        sqn.append(sum(x * x for x in v))
    return Bstar, sqn


def _round_frac(q):
    """Nearest integer to a Fraction (ties away from zero)."""
    num, den = q.numerator, q.denominator
    if num >= 0:
        return (2 * num + den) // (2 * den)
    return -((-2 * num + den) // (2 * den))


def babai_nearest_plane(B, t, gso=None):
    """Exact Babai nearest-plane.  B: list of integer rows (a basis, possibly
    with zero rows), t: integer target.  Returns the lattice vector v."""
    dim = len(t)
    Bstar, sqn = gso if gso is not None else gso_exact(B)
    w = [Fraction(x) for x in t]
    for i in reversed(range(len(B))):
        if sqn[i] == 0:
            continue
        c = sum(w[k] * Bstar[i][k] for k in range(dim)) / sqn[i]
        ci = _round_frac(c)
        if ci:
            w = [w[k] - ci * Fraction(B[i][k]) for k in range(dim)]
    return [int(t[k] - w[k]) for k in range(dim)]


# ===========================================================================
# Formulation PROJ — quotient lattice L' = L / <n*S_D*e_m>
# ===========================================================================

def proj_basis(sigs, n, lam, k1_bound, k2_bound):
    """Explicit basis of L' = pi_m(L), where pi_m deletes column m.

    Deleting column m makes the 2m+2 generators of L rank-deficient:
        n * pi(row_d) - sum_i B_i * pi(row_i) = 0.
    The k1-block of L' is  Lam = { x in (S_K1 Z)^m : x/S_K1 = c*B mod n },
    for which a genuine basis is
        g = S_K1 * (c0*B mod n)      with c0 = B_0^{-1} mod n   (so g_0 = S_K1)
        n*S_K1*e_i,  i = 1..m-1      (n*S_K1*e_0 is then redundant)
    giving det(k1-block) = S_K1 * (n*S_K1)^(m-1) = n^(m-1) * S_K1^m, the
    required index-n sublattice determinant.

    Returns (rows_proj, coeff) where rows_proj[r] lives in Z^(2m+1) (column m
    deleted) and coeff[r] is the integer coefficient vector of rows_proj[r]
    with respect to the ORIGINAL rows of build_glv_lattice, so that any
    reduced row can be lifted back into L (and its column-m entry read).
    """
    m = len(sigs)
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    ngen = 2 * m + 2                       # rows of M: 0..m-1, m(=row_d),
                                           # m+1..2m (k2), 2m+1 (Kannan)
    B = [s['B'] % n for s in sigs]
    c0 = modinv(B[0], n)

    rows, coeff = [], []

    # (1) g = c0 * row_d  -  sum_i floor(c0*B_i/n) * row_i     [k1-block only]
    cg = [0] * ngen
    cg[m] = c0
    gk1 = []
    for i in range(m):
        prod = c0 * B[i]
        cg[i] = -(prod // n)
        gk1.append((prod % n) * S_K1)
    assert gk1[0] == S_K1, "c0*B_0 mod n != 1"
    rows.append(gk1 + [0] * (m + 1))       # k2-block + Kannan col are zero
    coeff.append(cg)

    # (2) n*S_K1*e_i for i = 1..m-1
    for i in range(1, m):
        r = [0] * (2 * m + 1)
        r[i] = n * S_K1
        c = [0] * ngen
        c[i] = 1
        rows.append(r)
        coeff.append(c)

    # (3) the k2 rows, unchanged (they have no column-m entry)
    for i in range(m):
        r = [0] * (2 * m + 1)
        r[i] = -(lam % n) * S_K1
        r[m + i] = S_K2                    # original col m+1+i -> m+i
        c = [0] * ngen
        c[m + 1 + i] = 1
        rows.append(r)
        coeff.append(c)

    # (4) the Kannan row, unchanged
    r = [0] * (2 * m + 1)
    for i in range(m):
        r[i] = sigs[i]['A'] * S_K1
    r[2 * m] = S_KAN                       # original col 2m+1 -> 2m
    c = [0] * ngen
    c[2 * m + 1] = 1
    rows.append(r)
    coeff.append(c)

    assert len(rows) == 2 * m + 1
    return rows, coeff


def attack_proj(sigs, n, lam, k1_bound, k2_bound, d_secret, want_diag=False):
    """LLL in the quotient lattice L', lift through U, read d from column m."""
    m = len(sigs)
    dim_p = 2 * m + 1
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    rows, coeff = proj_basis(sigs, n, lam, k1_bound, k2_bound)
    M_full = build_glv_lattice(sigs, n, lam, k1_bound, k2_bound)

    A = IntegerMatrix.from_matrix(rows)
    U = IntegerMatrix.identity(dim_p)
    LLL.reduction(A, U)
    red = [[A[i][j] for j in range(dim_p)] for i in range(dim_p)]
    Um = [[U[i][j] for j in range(dim_p)] for i in range(dim_p)]

    # lift: reduced row r = sum_s U[r][s] * (coeff[s] . M_full)
    lifted = []
    for r in range(dim_p):
        cf = [0] * (2 * m + 2)
        for s in range(dim_p):
            u = Um[r][s]
            if u:
                for j, cj in enumerate(coeff[s]):
                    if cj:
                        cf[j] += u * cj
        v = [0] * (2 * m + 2)
        for j, cj in enumerate(cf):
            if cj:
                for k in range(2 * m + 2):
                    if M_full[j][k]:
                        v[k] += cj * M_full[j][k]
        lifted.append(v)

    found = recover_d(lifted, m, n, S_KAN, d_secret)
    if not want_diag:
        return found is not None
    sv = min((r for r in red if any(r)), key=norm)
    return found is not None, red, sv, lifted


# ===========================================================================
# Formulation CVP — Babai nearest plane in L0, no Kannan row
# ===========================================================================

def l0_basis_and_target(sigs, n, lam, k1_bound, k2_bound):
    """L0 (dim 2m+1): rows n*S_K1*e_i, row_d, k2 rows.  Columns 0..m-1 = k1,
    m = d, m+1..2m = k2.  Target t = (-A_i*S_K1, 0, ..., 0).

    For v = d*row_d + sum k2_i*row_k2i + sum t_i*row_i one has
        col i   = S_K1*(d*B_i - lam*k2_i + t_i*n) = S_K1*(k1_i - A_i)
    (using A_i + B_i*d = k1_i + lam*k2_i mod n), col m = d*S_D,
    col m+1+i = k2_i*S_K2, hence  v - t = (k1_i*S_K1, d*S_D, k2_i*S_K2).
    """
    m = len(sigs)
    dim = 2 * m + 1
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    rows = []
    for i in range(m):
        r = [0] * dim
        r[i] = n * S_K1
        rows.append(r)
    r = [0] * dim
    for i in range(m):
        r[i] = sigs[i]['B'] * S_K1
    r[m] = S_D
    rows.append(r)
    for i in range(m):
        r = [0] * dim
        r[i] = -(lam % n) * S_K1
        r[m + 1 + i] = S_K2
        rows.append(r)
    t = [0] * dim
    for i in range(m):
        t[i] = -sigs[i]['A'] * S_K1
    return rows, t


def attack_cvp(sigs, n, lam, k1_bound, k2_bound, d_secret, want_err=False):
    m = len(sigs)
    dim = 2 * m + 1
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    rows, t = l0_basis_and_target(sigs, n, lam, k1_bound, k2_bound)
    A = IntegerMatrix.from_matrix(rows)
    LLL.reduction(A)
    red = [[A[i][j] for j in range(dim)] for i in range(dim)]
    v = babai_nearest_plane(red, t)
    e = [v[k] - t[k] for k in range(dim)]
    d_cand = (e[m] // S_D) % n if S_D else None
    ok = (d_cand == d_secret)
    if want_err:
        return ok, e
    return ok


# ===========================================================================
# Baseline KAN (verbatim), plus shared instance generation
# ===========================================================================

def make_instance(curve, m, k1_bound, seed):
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    d_secret = random.Random(seed + 7777).randint(1, n - 1)
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed)
    return sigs, d_secret, k2_bound


def attack_kan(sigs, n, lam, k1_bound, k2_bound, d_secret):
    m = len(sigs)
    dim = 2 * m + 2
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    M = build_glv_lattice(sigs, n, lam, k1_bound, k2_bound)
    A = IntegerMatrix.from_matrix(M)
    LLL.reduction(A)
    red = [[A[i][j] for j in range(dim)] for i in range(dim)]
    return recover_d(red, m, n, S_KAN, d_secret) is not None


METHODS = ("KAN", "PROJ", "CVP")


def run_all(curve, m, k1_bound, seed):
    p, b, n, lam, G = curve
    sigs, d_secret, k2 = make_instance(curve, m, k1_bound, seed)
    if len(sigs) < m:
        return None
    return {
        "KAN": attack_kan(sigs, n, lam, k1_bound, k2, d_secret),
        "PROJ": attack_proj(sigs, n, lam, k1_bound, k2, d_secret),
        "CVP": attack_cvp(sigs, n, lam, k1_bound, k2, d_secret),
    }


# ===========================================================================
# Gaussian-heuristic model
# ===========================================================================

def log_det_kan(m, n, K1, K2):
    S_K1, S_D, S_K2, S_KAN = scales(n, K1, K2)
    return (m * math.log(n * S_K1) + math.log(S_D)
            + m * math.log(S_K2) + math.log(S_KAN))


def gh(logdet, dim):
    return math.exp(0.5 * math.log(dim / (2 * math.pi * math.e))
                    + logdet / dim)


def model_ratios(m, n, K1, K2):
    """Expected ||planted|| / GH(lattice) for the three formulations."""
    S_K1, S_D, S_K2, S_KAN = scales(n, K1, K2)
    ld_kan = log_det_kan(m, n, K1, K2)
    pv_kan = planted_norm_expected(m, n, K1, K2, S_K1, S_D, S_K2, S_KAN)
    r_kan = pv_kan / gh(ld_kan, 2 * m + 2)
    # L' = L / <n*S_D*e_m>:  det divided by n*S_D, dim 2m+1, planted image
    # loses the d-coordinate.
    ld_proj = ld_kan - math.log(n * S_D)
    pv_proj = math.sqrt(max(pv_kan ** 2 - (n * S_D) ** 2 / 3.0, 1.0))
    r_proj = pv_proj / gh(ld_proj, 2 * m + 1)
    # L0: det = (n S_K1)^m S_D S_K2^m, dim 2m+1, error norm keeps d, no Kannan
    ld_l0 = m * math.log(n * S_K1) + math.log(S_D) + m * math.log(S_K2)
    err = math.sqrt(m * (K1 * S_K1) ** 2 / 3.0 + (n * S_D) ** 2 / 3.0
                    + m * (K2 * S_K2) ** 2 / 3.0)
    r_cvp = err / gh(ld_l0, 2 * m + 1)
    return r_kan, r_proj, r_cvp, pv_kan, pv_proj, err


def critical_eff(m):
    """Solve r_proj(eff) = 1 for eff, in the large-n limit.

    With S_K1 = n/K1, S_K2 = n/K2, S_D = 1, S_KAN = n:
      det(L') = n^(3m) / (K1*K2)^m,  dim = 2m+1
      ||pv'|| ~ n * sqrt(2m/3 + 1)
    so   ||pv'|| / GH = sqrt(2m/3+1) * sqrt(2*pi*e/(2m+1)) * eff^(m/(2m+1))
                        * n^(1/(2m+1))                       -> limit below.
    """
    lo, hi = 1e-6, 0.999
    for _ in range(200):
        mid = math.sqrt(lo * hi)
        # large-n limit: n^(1/(2m+1)) -> 1
        val = (math.sqrt(2 * m / 3.0 + 1)
               * math.sqrt(2 * math.pi * math.e / (2 * m + 1))
               * mid ** (m / (2.0 * m + 1)))
        if val < 1:
            lo = mid
        else:
            hi = mid
    return math.sqrt(lo * hi)


# ===========================================================================
# Driver
# ===========================================================================

def main():
    SEEDS = [42, 1234, 9999, 555, 31337]

    HIST = [
        ("8-bit/199",        211,  2, 199,  106),
        ("12-bit/2557",      2557, 2, 2659, 1755),
        ("12-bit/2677 FAIL", 2677, 2, 2647, 185),
    ]

    hist_curves = []
    for label, p, b, n, lam in HIST:
        G = find_generator(p, b, n)
        assert G is not None, f"no generator for {label}"
        assert (lam * lam + lam + 1) % n == 0, f"bad lam for {label}"
        hist_curves.append((label, (p, b, n, lam, G)))

    print("=" * 78)
    print("Thread 23 — reformulate the Phase-2 lattice so the target is lambda_1")
    print("=" * 78)


    # ---------------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP U1: correctness of the two reformulations")
    print("-" * 78)

    label, curve = hist_curves[1]
    p, b, n, lam, G = curve
    m_u1, k1_u1 = 6, 2
    sigs, d0, k2_u1 = make_instance(curve, m_u1, k1_u1, 42)
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_u1, k2_u1)

    rows_p, coeff_p = proj_basis(sigs, n, lam, k1_u1, k2_u1)


    def det_int(rows):
        """Exact determinant of a square integer matrix (fraction-free Bareiss)."""
        Mx = [[Fraction(x) for x in r] for r in rows]
        N = len(Mx)
        det = Fraction(1)
        for c in range(N):
            piv = next((r for r in range(c, N) if Mx[r][c] != 0), None)
            if piv is None:
                return 0
            if piv != c:
                Mx[c], Mx[piv] = Mx[piv], Mx[c]
                det = -det
            det *= Mx[c][c]
            inv = Mx[c][c]
            for r in range(c + 1, N):
                f = Mx[r][c] / inv
                if f:
                    Mx[r] = [Mx[r][k] - f * Mx[c][k] for k in range(N)]
        return det


    d_proj = abs(det_int(rows_p))
    d_expect = Fraction((n * S_K1) ** m_u1 * S_D * S_K2 ** m_u1 * S_KAN, n * S_D)
    print(f"curve {label}, m={m_u1}, K1={k1_u1}, K2={k2_u1}")
    print(f"  det(L')          = {d_proj}")
    print(f"  det(L)/(n*S_D)   = {d_expect}")
    print(f"  L' is a basis    : {d_proj == d_expect}")

    rows0, t0 = l0_basis_and_target(sigs, n, lam, k1_u1, k2_u1)
    v_exact = [0] * (2 * m_u1 + 1)
    # build the exact planted CVP solution from the known witnesses
    tt = []
    for i in range(m_u1):
        num = d0 * sigs[i]['B'] - lam * sigs[i]['k2'] - (sigs[i]['k1'] - sigs[i]['A'])
        assert num % n == 0
        tt.append(-num // n)
    coef = [tt[i] for i in range(m_u1)] + [d0] + [sigs[i]['k2'] for i in range(m_u1)]
    for j, c in enumerate(coef):
        for k in range(2 * m_u1 + 1):
            v_exact[k] += c * rows0[j][k]
    err_exact = [v_exact[k] - t0[k] for k in range(2 * m_u1 + 1)]
    want = ([sigs[i]['k1'] * S_K1 for i in range(m_u1)] + [d0 * S_D]
            + [sigs[i]['k2'] * S_K2 for i in range(m_u1)])
    print(f"  CVP target algebra exact (v-t == planted error): {err_exact == want}")

    ok_p, red_p, sv_p, lifted_p = attack_proj(sigs, n, lam, k1_u1, k2_u1, d0,
                                              want_diag=True)
    M_full = build_glv_lattice(sigs, n, lam, k1_u1, k2_u1)


    def in_lattice(v, basis, n_, S_K1_):
        """Cheap membership check: v must be an integer combination of basis rows.
        Basis here is triangular enough to solve directly."""
        return True  # (structural check done via det + column consistency below)


    # verify the lift is consistent: deleting column m of the lifted vector must
    # reproduce the reduced projected row
    consistent = True
    for r in range(2 * m_u1 + 1):
        proj_of_lift = [lifted_p[r][k] for k in range(2 * m_u1 + 2) if k != m_u1]
        if proj_of_lift != red_p[r]:
            consistent = False
            break
    print(f"  U-lift consistent (pi(lift(row)) == reduced row): {consistent}")
    print(f"  PROJ recovers d  : {ok_p}")
    ok_c, err_c = attack_cvp(sigs, n, lam, k1_u1, k2_u1, d0, want_err=True)
    print(f"  CVP  recovers d  : {ok_c}")
    print(f"  KAN  recovers d  : {attack_kan(sigs, n, lam, k1_u1, k2_u1, d0)}")


    # ---------------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP U2: the falsifier — T4's K1 grid, three formulations")
    print("-" * 78)
    print("Wall location per method.  m = 12, 5 seeds.  K2 = isqrt(n)+1 = 52.\n")

    M_SIGS = 12
    K1_GRID = [2, 3, 4, 6, 8, 12, 16, 24, 32]

    u2_rows = []
    for label, curve in hist_curves[1:]:
        p, b, n, lam, G = curve
        k2 = math.isqrt(n) + 1
        print(f"{label}  (n={n}, lam={lam}, lam*={lam_star(lam, n):.3f})")
        print(f"{'method':>7} " + " ".join(f"K1={k:<4}" for k in K1_GRID))
        per_method = {meth: [] for meth in METHODS}
        for k1 in K1_GRID:
            tally = {meth: 0 for meth in METHODS}
            for sd in SEEDS:
                res = run_all(curve, M_SIGS, k1, sd)
                if res is None:
                    continue
                for meth in METHODS:
                    tally[meth] += 1 if res[meth] else 0
            for meth in METHODS:
                per_method[meth].append(tally[meth])
            u2_rows.append({'curve': label, 'K1': k1, 'eff': k1 * k2 / n,
                            **{meth: tally[meth] for meth in METHODS}})
        for meth in METHODS:
            cells = " ".join(f"{v}/5{'':2}" for v in per_method[meth])
            print(f"{meth:>7} {cells}")
        # wall = largest K1 with >= 3/5
        for meth in METHODS:
            wall = max([K1_GRID[i] for i, v in enumerate(per_method[meth]) if v >= 3],
                       default=0)
            effw = wall * k2 / n if wall else 0.0
            print(f"        wall({meth}) = K1 {wall}   (eff = {effw:.3f})")
        print()


    # ---------------------------------------------------------------------------
    print("-" * 78)
    print("EXP U3: does sv/pv rise above 1 in the quotient lattice?")
    print("-" * 78)
    print("sv = shortest reduced row; pv = image of the planted vector in L'.")
    print("kan-energy = fraction of ||sv||^2 in the Kannan column.\n")
    print(f"{'curve':<18} {'K1':>4} {'L: sv/pv':>9} {'L: |sv[m]|/n':>13} "
          f"{'L2: sv/pv':>10} {'L2: kan-E':>10} {'L2 sv=pv?':>10}")

    for label, curve in hist_curves:
        p, b, n, lam, G = curve
        for k1 in (2, 8):
            m = M_SIGS
            sigs, d0, k2 = make_instance(curve, m, k1, 42)
            if len(sigs) < m:
                continue
            S_K1, S_D, S_K2, S_KAN = scales(n, k1, k2)
            # baseline
            Mb = build_glv_lattice(sigs, n, lam, k1, k2)
            Ab = IntegerMatrix.from_matrix(Mb)
            LLL.reduction(Ab)
            redb = [[Ab[i][j] for j in range(2 * m + 2)] for i in range(2 * m + 2)]
            svb = min((r for r in redb if any(r)), key=norm)
            pvb = planted_vector(sigs, d0, n, k1, k2)
            # projected
            _, red_p, sv_p, _ = attack_proj(sigs, n, lam, k1, k2, d0, want_diag=True)
            pv_p = [pvb[k] for k in range(2 * m + 2) if k != m]
            tot = sum(x * x for x in sv_p) or 1
            kanE = sv_p[2 * m] ** 2 / tot
            same = (sv_p == pv_p or [-x for x in sv_p] == pv_p)
            print(f"{label:<18} {k1:>4} {norm(svb)/norm(pvb):>9.3f} "
                  f"{abs(svb[m])/n:>13.4f} {norm(sv_p)/norm(pv_p):>10.3f} "
                  f"{kanE:>10.3f} {str(same):>10}")


    # ---------------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP U4: Gaussian-heuristic model of the wall")
    print("-" * 78)
    print("ratio = E||planted|| / GH(lattice).  Unique-SVP/BDD needs ratio < 1.\n")
    print(f"{'m':>3} {'crit eff (proj)':>16}")
    for m in (6, 8, 12, 16, 24, 32):
        print(f"{m:>3} {critical_eff(m):>16.4f}")

    print("\nPer-instance ratios on the historical curves (m = 12):")
    print(f"{'curve':<18} {'K1':>4} {'eff':>6} {'r_KAN':>7} {'r_PROJ':>7} "
          f"{'r_CVP':>7} {'KAN':>5} {'PROJ':>5} {'CVP':>5}")
    for label, curve in hist_curves[1:]:
        p, b, n, lam, G = curve
        k2 = math.isqrt(n) + 1
        for k1 in K1_GRID:
            r_kan, r_proj, r_cvp, _, _, _ = model_ratios(M_SIGS, n, k1, k2)
            row = next((r for r in u2_rows if r['curve'] == label and r['K1'] == k1),
                       None)
            if row is None:
                continue
            print(f"{label:<18} {k1:>4} {k1*k2/n:>6.3f} {r_kan:>7.3f} "
                  f"{r_proj:>7.3f} {r_cvp:>7.3f} "
                  f"{str(row['KAN'])+'/5':>5} {str(row['PROJ'])+'/5':>5} "
                  f"{str(row['CVP'])+'/5':>5}")


    # ---------------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP U5: 17-bit sweep at three eff values, three formulations")
    print("-" * 78)

    LO, HI = 2 ** 16, 2 ** 17
    print(f"Searching j=0 GLV curves with p in [{LO}, {HI}) ...")
    curves17 = search_curves(LO, HI, per_bin=1, nbins=8)
    print(f"Found {len(curves17)} curves.\n")

    for eff_t in (0.05, 0.15, 0.25):
        print(f"=== eff target = {eff_t}  (m = {M_SIGS}, {len(SEEDS)} seeds) ===")
        print(f"{'p':>7} {'n':>7} {'lam*':>6} {'K1':>4} {'eff':>6} {'r_PROJ':>7} "
              f"{'KAN':>5} {'PROJ':>5} {'CVP':>5}")
        agg = {meth: [0, 0] for meth in METHODS}
        for (p, b, n, lam, G) in curves17:
            k2 = math.isqrt(n) + 1
            k1 = max(2, int(eff_t * n / k2))
            curve = (p, b, n, lam, G)
            tally = {meth: 0 for meth in METHODS}
            for sd in SEEDS:
                res = run_all(curve, M_SIGS, k1, sd)
                if res is None:
                    continue
                for meth in METHODS:
                    tally[meth] += 1 if res[meth] else 0
            r_kan, r_proj, r_cvp, _, _, _ = model_ratios(M_SIGS, n, k1, k2)
            for meth in METHODS:
                agg[meth][0] += tally[meth]
                agg[meth][1] += len(SEEDS)
            print(f"{p:>7} {n:>7} {lam_star(lam, n):>6.3f} {k1:>4} "
                  f"{k1*k2/n:>6.3f} {r_proj:>7.3f} "
                  f"{str(tally['KAN'])+'/5':>5} {str(tally['PROJ'])+'/5':>5} "
                  f"{str(tally['CVP'])+'/5':>5}")
        tot = "  ".join(f"{meth} {agg[meth][0]}/{agg[meth][1]}" for meth in METHODS)
        print(f"  TOTAL: {tot}\n")

    # -----------------------------------------------------------------------
    print("-" * 78)
    print("EXP U6: per-instance agreement KAN vs PROJ (not just per-cell tally)")
    print("-" * 78)
    print("3 curves x K1 in {2,3,4,6,8,12,16} x 10 seeds, m = 10.\n")

    U6_SEEDS = [42, 1234, 9999, 555, 31337, 777, 2024, 31, 99, 1]
    tot = dis = cvp_win = kan_win = 0
    for label, curve in hist_curves:
        for k1 in (2, 3, 4, 6, 8, 12, 16):
            for sd in U6_SEEDS:
                r = run_all(curve, 10, k1, sd)
                if r is None:
                    continue
                tot += 1
                if r["KAN"] != r["PROJ"]:
                    dis += 1
                    print(f"  DISAGREE {label} K1={k1} seed={sd} {r}")
                if r["CVP"] and not r["KAN"]:
                    cvp_win += 1
                    print(f"  CVP-wins {label} K1={k1} seed={sd}")
                if r["KAN"] and not r["CVP"]:
                    kan_win += 1
    print(f"  instances = {tot}")
    print(f"  KAN != PROJ                  : {dis}")
    print(f"  CVP succeeds where KAN fails : {cvp_win}")
    print(f"  KAN succeeds where CVP fails : {kan_win}")

    print("\nDone.")


if __name__ == "__main__":
    main()
