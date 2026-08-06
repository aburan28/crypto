"""
Thread 23 — reformulate the GLV-HNP Phase-2 lattice so the planted vector can be
lambda_1.

Motivation (2026-07-29 log, finding T5): in the Phase-2 lattice L of dimension
2m+2 the shortest vector is ALWAYS the trivial vector n*S_D*e_m, because
    n*row_m - sum_i B_i*row_i = (0,...,0, n*S_D, 0,...,0),
and ||n*S_D*e_m|| = n while ||v_planted|| ~ n*sqrt(2m/3 + 4/3).  So the planted
vector is never lambda_1 and recovery is a BDD/coset condition, not SVP.

This script executes the proposed fix: project L along e_m.  Since
L cap span(e_m) = n*S_D*Z*e_m is exactly the trivial-vector line, the orthogonal
projection pi = pi_{e_m^perp} maps L onto a lattice pi(L) of rank 2m+1, and the
trivial vector maps to 0.  Concretely: DELETE COLUMN m from the basis.  The 2m+2
rows then have rank 2m+1 and LLL emits exactly one zero row -- the trivial vector,
annihilated.

d is no longer read off directly (its column is gone), so recovery instead reads
k1_i and k2_i off the surviving coordinates and solves
    d = (k1_0 + lam*k2_0 - A_0) * B_0^{-1}  (mod n)
verifying against all m congruences.  That is an ORACLE-FREE recovery; for a fair
comparison the baseline lattice is scored with an oracle-free recovery too (the
2026-07-29 harness used d_cand == d_secret, which is strictly easier).

Experiments
  E0  structural sanity: rank, zero row, trivial vector gone
  E1  sv/pv on the three historical anchor curves, L vs pi(L)
  E2  the K1 wall (2026-07-29 T4 grid) re-run on pi(L)
  E3  17-bit sweep across eff = K1*K2/n, L vs pi(L)
  E4  Gaussian-heuristic predictor gamma = ||pi(v)|| / lambda_1^GH(pi L)

Falsifier (stated 2026-07-29): "if sv/pv rises above 1 after the reformulation
AND the K1 wall in T4 moves outward on the lam*=0.07 curve (currently K1~4-6),
the reformulation is a real improvement; if the wall stays at K1~4-6, the wall is
information-theoretic and Phase 2 is at its ceiling."

Run: python3 glv_hnp_thread23_projected.py
"""

import math
import os
import random
import sys

from fpylll import IntegerMatrix, LLL

# ---------------------------------------------------------------------------
# Reuse the 2026-07-29 helper prefix VERBATIM so the comparison is exact.
# glv_hnp_phase2_lambda_threshold.py runs its experiments at module level, so we
# exec only the part above the first experiment banner.
# ---------------------------------------------------------------------------

_HERE = os.path.dirname(os.path.abspath(__file__))
_SRC = os.path.join(_HERE, "glv_hnp_phase2_lambda_threshold.py")
with open(_SRC) as _f:
    _text = _f.read()
_MARKER = "# EXPERIMENT T1"
assert _MARKER in _text, "helper-prefix marker moved; update _MARKER"
_prefix = _text.split(_MARKER)[0]
# drop the trailing '# ====' banner line of the prefix
exec(compile(_prefix, _SRC, "exec"), globals())

# names now in scope: modinv, ec_add, ec_mul, tonelli_shanks, find_generator,
# eisenstein_decompose, j0_traces, glv_roots, lam_star, gauss_reduce_2d,
# lambda_block_mu, planted_norm_expected, scales, gen_signatures,
# build_glv_lattice, planted_vector, norm, recover_d, build_curve, search_curves

SEEDS = [42, 1234, 9999, 555, 31337]


# ---------------------------------------------------------------------------
# Projected lattice: delete column m
# ---------------------------------------------------------------------------

def build_glv_lattice_proj(sigs, n, lam, k1_bound, k2_bound):
    """pi_{e_m^perp}(L): the Phase-2 basis with the d-column removed.
    2m+2 rows spanning a rank-(2m+1) lattice in Z^{2m+1}.
    Column layout: [0..m-1] k1-cols, [m..2m-1] k2-cols, [2m] Kannan col."""
    M = build_glv_lattice(sigs, n, lam, k1_bound, k2_bound)
    m = len(sigs)
    return [row[:m] + row[m + 1:] for row in M]


def planted_vector_proj(sigs, n, k1_bound, k2_bound):
    m = len(sigs)
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    v = [0] * (2 * m + 1)
    for i in range(m):
        v[i] = sigs[i]['k1'] * S_K1
        v[m + i] = sigs[i]['k2'] * S_K2
    v[2 * m] = S_KAN
    return v


def det_proj(m, n, k1_bound, k2_bound):
    """det pi(L) = det(L) / (n*S_D) = (n*S_K1)^m * S_K2^m."""
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    return (n * S_K1) ** m * S_K2 ** m


def gh_lambda1(det, dim):
    """Gaussian-heuristic lambda_1 = sqrt(dim/(2*pi*e)) * det^(1/dim)."""
    return math.sqrt(dim / (2 * math.pi * math.e)) * math.exp(math.log(det) / dim)


# ---------------------------------------------------------------------------
# Oracle-free recovery
# ---------------------------------------------------------------------------

def _d_from_row(w, sigs, n, lam, k1_bound, k2_bound, S_K1, S_K2):
    """Given a candidate (sign-normalised) row w in the projected coords,
    read k1_i,k2_i and solve for d.  Returns d or None.  No oracle."""
    m = len(sigs)
    k1s, k2s = [], []
    for i in range(m):
        a, b = w[i], w[m + i]
        if a % S_K1 or b % S_K2:
            return None
        k1, k2 = a // S_K1, b // S_K2
        if not (0 <= k1 < k1_bound) or not (0 <= k2 < k2_bound):
            return None
        k1s.append(k1)
        k2s.append(k2)
    B0 = sigs[0]['B'] % n
    if B0 == 0:
        return None
    d = (k1s[0] + lam * k2s[0] - sigs[0]['A']) * modinv(B0, n) % n
    for i in range(m):
        if (sigs[i]['A'] + sigs[i]['B'] * d - k1s[i] - lam * k2s[i]) % n:
            return None
    return d


def recover_d_proj(reduced, sigs, n, lam, k1_bound, k2_bound):
    m = len(sigs)
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    for row in reduced:
        last = row[2 * m]
        if abs(last) != S_KAN:
            continue
        w = row if last > 0 else [-x for x in row]
        d = _d_from_row(w, sigs, n, lam, k1_bound, k2_bound, S_K1, S_K2)
        if d is not None:
            return d
    return None


def glv_small(k, n, lam, k1_bound, k2_bound):
    """Oracle-free test: is k GLV-decomposable as k1 + lam*k2 (mod n) with
    0 <= k1 < k1_bound, 0 <= k2 < k2_bound?  O(k2_bound) ~ O(sqrt n)."""
    for k2 in range(k2_bound):
        if (k - lam * k2) % n < k1_bound:
            return True
    return False


def glv_verify(d, sigs, n, lam, k1_bound, k2_bound):
    """Oracle-free verification of a candidate d: every nonce implied by d must
    be GLV-small.  Uses only (A_i, B_i), never sigs[i]['k1'/'k2']."""
    for s in sigs:
        if not glv_small((s['A'] + s['B'] * d) % n, n, lam, k1_bound, k2_bound):
            return False
    return True


def recover_d_full_blind(reduced, sigs, n, lam, k1_bound, k2_bound):
    """Oracle-free recovery on the UNprojected lattice, so both columns of every
    table are scored by the same rule.  The full lattice keeps one extra route
    the projection destroys: reading d straight off the d-column."""
    m = len(sigs)
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    for row in reduced:
        last = row[2 * m + 1]
        if abs(last) != S_KAN:
            continue
        w = row if last > 0 else [-x for x in row]
        wp = w[:m] + w[m + 1:]
        d = _d_from_row(wp, sigs, n, lam, k1_bound, k2_bound, S_K1, S_K2)
        if d is not None:
            return d
        d2 = w[m] % n
        if d2 and glv_verify(d2, sigs, n, lam, k1_bound, k2_bound):
            return d2
    return None


# ---------------------------------------------------------------------------
# Experiment drivers
# ---------------------------------------------------------------------------

def lll_rows(M, ncols):
    A = IntegerMatrix.from_matrix(M)
    LLL.reduction(A)
    return [[A[i][j] for j in range(ncols)] for i in range(A.nrows)]


def run_once(curve, m, d_secret, k1_bound, seed, mode):
    """mode in {'full','proj'}.  Returns (ok, sigs, reduced_rows)."""
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return None, None, None
    if mode == 'proj':
        M = build_glv_lattice_proj(sigs, n, lam, k1_bound, k2_bound)
        red = lll_rows(M, 2 * m + 1)
        d = recover_d_proj(red, sigs, n, lam, k1_bound, k2_bound)
    else:
        M = build_glv_lattice(sigs, n, lam, k1_bound, k2_bound)
        red = lll_rows(M, 2 * m + 2)
        d = recover_d_full_blind(red, sigs, n, lam, k1_bound, k2_bound)
    return (d == d_secret), sigs, red


def success(curve, m, k1_bound, seeds, mode):
    n = curve[2]
    ok = 0
    for s in seeds:
        rng = random.Random(s ^ 0xC0FFEE)
        d_secret = rng.randrange(1, n)
        r, _, _ = run_once(curve, m, d_secret, k1_bound, s, mode)
        if r:
            ok += 1
    return ok


def shortest_nonzero(rows):
    best, bv = None, None
    for r in rows:
        nn = sum(x * x for x in r)
        if nn == 0:
            continue
        if best is None or nn < best:
            best, bv = nn, r
    return math.sqrt(best), bv


def hline(t=""):
    print("=" * 78)
    if t:
        print(t)
        print("=" * 78)


# ===========================================================================
hline("Thread 23 — project the Phase-2 lattice along e_m")
print(f"seeds = {SEEDS}")

HIST = [
    # label,               p,    b, n,    lam,  K1, m
    ("8-bit/199",          211,  2, 199,  106,  2,  6),
    ("12-bit/2557",        2557, 2, 2659, 1755, 8,  8),
    ("12-bit/2677 FAIL",   2677, 2, 2647, 185,  8,  10),
]
hist = []
for label, p, b, n, lam, k1, m in HIST:
    G = find_generator(p, b, n)
    assert G is not None and (lam * lam + lam + 1) % n == 0, label
    hist.append((label, (p, b, n, lam, G), k1, m))

# ---------------------------------------------------------------------------
hline("E0 — structural sanity: is the trivial vector actually annihilated?")
print(f"{'curve':<20} {'dim L':>6} {'dim piL':>8} {'zero rows':>10} "
      f"{'|triv|/n':>9} {'triv in piL':>12}")
for label, curve, k1, m in hist:
    p, b, n, lam, G = curve
    k2b = math.isqrt(n) + 1
    rng = random.Random(1)
    d_secret = rng.randrange(1, n)
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1, k2b, 42)
    Mf = build_glv_lattice(sigs, n, lam, k1, k2b)
    Mp = build_glv_lattice_proj(sigs, n, lam, k1, k2b)
    redf = lll_rows(Mf, 2 * m + 2)
    redp = lll_rows(Mp, 2 * m + 1)
    zf = sum(1 for r in redf if not any(r))
    zp = sum(1 for r in redp if not any(r))
    svf, _ = shortest_nonzero(redf)
    # is any nonzero projected row equal to the image of the trivial vector? (=0)
    print(f"{label:<20} {2*m+2:>6} {2*m+1:>8} {str(zf)+'/'+str(zp):>10} "
          f"{svf/n:>9.4f} {str(zp == 1):>12}")
print("\n'zero rows' is full/projected.  |triv|/n == 1.0000 confirms the")
print("2026-07-29 T5 result on L; zp == 1 confirms exactly one dependency in")
print("pi(L), i.e. the trivial vector and nothing else was removed.")

# ---------------------------------------------------------------------------
hline("E1 — sv/pv (shortest reduced vector / planted vector), L vs pi(L)")
print(f"{'curve':<20} {'K1':>4} {'m':>3} {'sv/pv L':>9} {'sv/pv piL':>10} "
      f"{'gamma_GH':>9} {'kan?':>6}")
for label, curve, k1, m in hist:
    p, b, n, lam, G = curve
    k2b = math.isqrt(n) + 1
    rng = random.Random(1)
    d_secret = rng.randrange(1, n)
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1, k2b, 42)
    S_K1, S_D, S_K2, S_KAN = scales(n, k1, k2b)
    pvf = planted_vector(sigs, d_secret, n, k1, k2b)
    pvp = planted_vector_proj(sigs, n, k1, k2b)
    redf = lll_rows(build_glv_lattice(sigs, n, lam, k1, k2b), 2 * m + 2)
    redp = lll_rows(build_glv_lattice_proj(sigs, n, lam, k1, k2b), 2 * m + 1)
    svf, _ = shortest_nonzero(redf)
    svp, bvp = shortest_nonzero(redp)
    gh = gh_lambda1(det_proj(m, n, k1, k2b), 2 * m + 1)
    print(f"{label:<20} {k1:>4} {m:>3} {svf/norm(pvf):>9.3f} "
          f"{svp/norm(pvp):>10.3f} {norm(pvp)/gh:>9.3f} "
          f"{str(abs(bvp[2*m]) == S_KAN):>6}")
print("\nsv/pv >= 1 means the planted vector IS the shortest vector of the")
print("reduced basis.  gamma_GH = ||pi(v)||/lambda_1^GH(pi L): < 1 predicts the")
print("planted vector is shorter than a random lattice vector of that covolume.")

# ---------------------------------------------------------------------------
hline("E2 — the K1 wall (2026-07-29 T4 grid) on L vs pi(L)")
K1_GRID = [2, 3, 4, 6, 8, 12, 16, 24]
print(f"{'curve':<20} {'lat':<5} " + " ".join(f"{k:>4}" for k in K1_GRID))
walls = {}
for label, curve, _k1, _m in hist[1:]:
    p, b, n, lam, G = curve
    ls = lam_star(lam, n)
    for mode in ('full', 'proj'):
        row = []
        for K1 in K1_GRID:
            row.append(success(curve, 8, K1, SEEDS, mode))
        walls[(label, mode)] = row
        print(f"{label:<20} {mode:<5} " + " ".join(f"{c:>4}" for c in row)
              + f"   lam*={ls:.3f}")
print("\ncell = #seeds (of 5) recovering d at m=8.  The 2026-07-29 T4 wall was")
print("K1<=12-16 (2557) and K1<=4-6 (2677).  If the proj row extends further")
print("right, the reformulation is a real improvement.")


# ---------------------------------------------------------------------------
hline("E3 — 17-bit sweep across eff = K1*K2/n, L vs pi(L)")
print("harvesting 17-bit j=0 GLV curves ...")
sys.stdout.flush()
curves = search_curves(2 ** 16, 2 ** 17, per_bin=1, nbins=10, max_primes=100000)
curves = curves[:10]
print(f"  {len(curves)} curves: " +
      ", ".join(f"p={c[0]}/lam*={lam_star(c[3], c[2]):.3f}" for c in curves))

M_SIGS = 12
EFFS = [0.05, 0.10, 0.15, 0.20, 0.25]
tot = {e: [0, 0] for e in EFFS}
gamma_rows = []
print(f"\n{'eff':>6} {'p':>8} {'lam*':>6} {'K1':>5} {'L':>4} {'piL':>5} "
      f"{'gamma_GH':>9}")
for eff in EFFS:
    for (p, b, n, lam, G) in curves:
        k2b = math.isqrt(n) + 1
        K1 = max(2, int(eff * n / k2b))
        curve = (p, b, n, lam, G)
        sf = success(curve, M_SIGS, K1, SEEDS, 'full')
        sp = success(curve, M_SIGS, K1, SEEDS, 'proj')
        tot[eff][0] += sf
        tot[eff][1] += sp
        S_K1, S_D, S_K2, S_KAN = scales(n, K1, k2b)
        pv_exp = math.sqrt(M_SIGS * (K1 * S_K1) ** 2 / 3.0
                           + M_SIGS * (k2b * S_K2) ** 2 / 3.0 + S_KAN ** 2)
        g = pv_exp / gh_lambda1(det_proj(M_SIGS, n, K1, k2b), 2 * M_SIGS + 1)
        gamma_rows.append((eff, p, sp, sf, g))
        print(f"{eff:>6.2f} {p:>8} {lam_star(lam, n):>6.3f} {K1:>5} "
              f"{sf:>4} {sp:>5} {g:>9.3f}")
den = len(curves) * len(SEEDS)
print(f"\n{'eff':>6} {'L total':>10} {'piL total':>11} {'gamma_GH (mean)':>16}")
for e in EFFS:
    gs = [g for (ee, _p, _sp, _sf, g) in gamma_rows if ee == e]
    print(f"{e:>6.2f} {str(tot[e][0])+'/'+str(den):>10} "
          f"{str(tot[e][1])+'/'+str(den):>11} {sum(gs)/len(gs):>16.3f}")

# ---------------------------------------------------------------------------
hline("E4 — is gamma_GH a predictor?  (threshold scan on the E3 data)")
pts = [(g, sp) for (_e, _p, sp, _sf, g) in gamma_rows]
ns = len(SEEDS)
best = None
for thr in [x / 100.0 for x in range(50, 260, 5)]:
    lo = [s for g, s in pts if g < thr]
    hi = [s for g, s in pts if g >= thr]
    if not lo or not hi:
        continue
    acc = (sum(lo) + sum(ns - s for s in hi)) / (ns * len(pts))
    if best is None or acc > best[1]:
        best = (thr, acc, sum(lo) / (ns * len(lo)), sum(hi) / (ns * len(hi)))
maj = max(sum(s for _g, s in pts), sum(ns - s for _g, s in pts)) / (ns * len(pts))
if best:
    print(f"best split gamma_GH < {best[0]:.2f}: accuracy {best[1]:.3f} "
          f"(majority baseline {maj:.3f})")
    print(f"  p_hat below threshold = {best[2]:.3f}, above = {best[3]:.3f}")
print("\nasymptotic prediction: with S_K1*K1 ~ S_K2*K2 ~ S_KAN ~ n and")
print("det pi(L) = n^{3m}/(K1 K2)^m, gamma_GH < 1 reduces (m -> inf) to")
print(f"  eff < 3/(2*pi*e) = {3/(2*math.pi*math.e):.4f}")
print("which brackets the 2026-07-29 T3 observation (19/20 at eff=0.05,")
print("3/20 at 0.15, 0/20 at 0.25).")

# ---------------------------------------------------------------------------
hline("E5 — what is the SECOND parasite?  (sv/pv stayed < 1 in E1)")
print("Hypothesis: the lambda-block sublattice L2_i = <(n*S_K1,0),(-lam*S_K1,S_K2)>")
print("supported on the coordinate pair (i, m+i).  It survives the projection")
print("(it never touched the d-column).  mu = lambda_1(L2) is the same quantity")
print("as the nu_hat separator of commit e845207 (nu_hat = mu/sqrt(det L2)).\n")
print(f"{'curve':<20} {'K1':>4} {'m':>3} {'||sv piL||':>11} {'mu':>10} "
      f"{'sv/mu':>7} {'pv/mu':>7} {'pairs':>6} {'lam-blk?':>9} {'nu_hat':>7}")


def block_support(v, m):
    """number of coordinate pairs (i, m+i) carrying nonzero mass, and whether
    the Kannan coordinate is zero."""
    pairs = sum(1 for i in range(m) if v[i] or v[m + i])
    return pairs, v[2 * m] == 0


e5_curves = [(lbl, c, k1, mm) for lbl, c, k1, mm in hist]
for (p, b, n, lam, G) in curves[:4]:
    e5_curves.append((f"17-bit/{p}", (p, b, n, lam, G),
                      max(2, int(0.05 * n / (math.isqrt(n) + 1))), 12))
for label, curve, k1, m in e5_curves:
    p, b, n, lam, G = curve
    k2b = math.isqrt(n) + 1
    S_K1, S_D, S_K2, S_KAN = scales(n, k1, k2b)
    rng = random.Random(1)
    d_secret = rng.randrange(1, n)
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1, k2b, 42)
    if len(sigs) < m:
        continue
    redp = lll_rows(build_glv_lattice_proj(sigs, n, lam, k1, k2b), 2 * m + 1)
    svn, sv = shortest_nonzero(redp)
    pvn = norm(planted_vector_proj(sigs, n, k1, k2b))
    mu, _w = lambda_block_mu(n, lam, S_K1, S_K2)
    nuhat = mu / math.sqrt(n * S_K1 * S_K2)
    pairs, kan0 = block_support(sv, m)
    print(f"{label:<20} {k1:>4} {m:>3} {svn:>11.1f} {mu:>10.1f} "
          f"{svn/mu:>7.3f} {pvn/mu:>7.3f} {pairs:>6} "
          f"{str(pairs == 1 and kan0):>9} {nuhat:>7.3f}")
print("\nsv/mu ~ 1.0 with pairs == 1 confirms the shortest vector of pi(L) is")
print("exactly a lambda-block vector.  pv/mu < 1 is then the per-instance")
print("condition for the planted vector to become lambda_1.")

# ---------------------------------------------------------------------------
hline("E6 — can rescaling S_K2 push ||pv|| below mu?")
print("mu ~ nu_hat*sqrt(n*S_K1*S_K2) scales as sqrt(c), ||pv||'s k2-part scales")
print("as c.  So shrinking S_K2 by c should shrink ||pv|| faster than mu.")


def scales_c(n, k1_bound, k2_bound, c2):
    return (n // k1_bound, 1, max(1, int(c2 * (n // k2_bound))), n)


def build_proj_c(sigs, n, lam, k1_bound, k2_bound, c2):
    m = len(sigs)
    S_K1, S_D, S_K2, S_KAN = scales_c(n, k1_bound, k2_bound, c2)
    dim = 2 * m + 1
    M = [[0] * dim for _ in range(2 * m + 2)]
    for i in range(m):
        M[i][i] = n * S_K1
    for i in range(m):
        M[m][i] = sigs[i]['B'] * S_K1
    for i in range(m):
        M[m + 1 + i][i] = -lam * S_K1
        M[m + 1 + i][m + i] = S_K2
    for i in range(m):
        M[2 * m + 1][i] = sigs[i]['A'] * S_K1
    M[2 * m + 1][dim - 1] = S_KAN
    return M


def run_once_c(curve, m, d_secret, k1_bound, seed, c2):
    p, b, n, lam, G = curve
    k2b = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2b, seed)
    if len(sigs) < m:
        return None
    S_K1, S_D, S_K2, S_KAN = scales_c(n, k1_bound, k2b, c2)
    red = lll_rows(build_proj_c(sigs, n, lam, k1_bound, k2b, c2), 2 * m + 1)
    for row in red:
        if abs(row[2 * m]) != S_KAN:
            continue
        w = row if row[2 * m] > 0 else [-x for x in row]
        d = _d_from_row(w, sigs, n, lam, k1_bound, k2b, S_K1, S_K2)
        if d is not None:
            return d == d_secret
    return False


C_GRID = [0.03125, 0.0625, 0.125, 0.25, 0.5, 1.0, 2.0]
print(f"\n{'curve':<20} {'K1':>4} " + " ".join(f"{c:>7}" for c in C_GRID)
      + f" {'| pv/mu at c':>14}")
for label, curve, _k1, _m in hist[1:]:
    p, b, n, lam, G = curve
    k2b = math.isqrt(n) + 1
    for K1 in (6, 8, 12):
        row, ratios = [], []
        for c2 in C_GRID:
            ok = 0
            for s in SEEDS:
                rng = random.Random(s ^ 0xC0FFEE)
                ds = rng.randrange(1, n)
                if run_once_c(curve, 8, ds, K1, s, c2):
                    ok += 1
            row.append(ok)
            S_K1, S_D, S_K2, S_KAN = scales_c(n, K1, k2b, c2)
            mu, _w = lambda_block_mu(n, lam, S_K1, S_K2)
            pv = math.sqrt(8 * (K1 * S_K1) ** 2 / 3.0
                           + 8 * (k2b * S_K2) ** 2 / 3.0 + S_KAN ** 2)
            ratios.append(pv / mu)
        print(f"{label:<20} {K1:>4} " + " ".join(f"{v:>7}" for v in row)
              + "  | " + " ".join(f"{r:.2f}" for r in ratios))
print("\ncell = #seeds (of 5) at m=8.  If success survives (or improves) while")
print("pv/mu drops below 1, the lambda-block parasite was the binding constraint.")

print("\nDone.")
