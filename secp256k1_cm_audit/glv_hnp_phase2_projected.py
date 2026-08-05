"""
GLV-HNP Phase 2 -- Thread 23: reformulate the lattice so the planted vector is
lambda_1.

Background (2026-07-29 run, RESEARCH_AUTOLAB_LOG.md "T5")
--------------------------------------------------------
The Phase-2 lattice L (dim 2m+2, `glv_hnp_phase2_20bit.py:263`) always contains
the vector

    t = n * S_D * e_m        (d-column axis, norm n*S_D)

because n*row_m reduces to it modulo the modular rows n*S_K1*e_i.  Meanwhile

    ||v_planted||^2 ~ n^2 * (2m/3 + 4/3)      (with S_D = 1)

so t is shorter than v_planted for every m >= 1, on every curve.  Measured
sv/pv sat in [0.337, 0.368] for successes and failures alike, and 100% of the
shortest vector's energy sat in the d-column.  t is information-free: d is only
defined mod n, and no choice of S_D removes t because both vectors scale
linearly in S_D.

So Phase-2 recovery was never an SVP condition.  It is a coset/BDD condition.

This script tests the fix proposed on 2026-07-29:

  (P) PROJECTED  -- quotient out Z*t, i.e. delete the d-column.  d is recovered
      afterwards from any single signature equation.  dim 2m+1.
  (C) CVP        -- drop the Kannan embedding too and call Babai nearest-plane
      against an explicitly centred target.  dim 2m.

Falsifier (stated 2026-07-29, verbatim)
---------------------------------------
  "if sv/pv rises above 1 after the reformulation AND the K1 wall in T4 moves
   outward on the lam*=0.07 curve (currently K1~4-6), the reformulation is a
   real improvement; if the wall stays at K1~4-6, then the wall is
   information-theoretic and Phase 2 is at its ceiling."

Arms
----
A0  closed-form norm / determinant / Gaussian-heuristic table, ORIG vs PROJ.
A1  correctness: PROJ and CVP recover d on the known-good anchors.
A2  sv/pv and the rank of v_planted in the LLL-reduced basis, ORIG vs PROJ.
A3  the T4 K1 wall, ORIG vs PROJ vs CVP, on both anchor curves, 5 seeds.
A4  does the Gaussian-heuristic ratio ||v||/GH = 1 predict where the wall is?

Run: python3 glv_hnp_phase2_projected.py
"""

import math
import random
import sys  # noqa: F401

from fpylll import IntegerMatrix, LLL, GSO

import importlib.util

_spec = importlib.util.spec_from_file_location(
    "_t20a", __file__.rsplit("/", 1)[0] + "/glv_hnp_phase2_lambda_threshold.py")
_t20a = importlib.util.module_from_spec(_spec)
try:
    _spec.loader.exec_module(_t20a)   # helpers only; the T1-T5 suite self-skips
except SystemExit:
    pass

modinv = _t20a.modinv
ec_mul = _t20a.ec_mul
find_generator = _t20a.find_generator
gen_signatures = _t20a.gen_signatures
scales = _t20a.scales
build_glv_lattice = _t20a.build_glv_lattice
planted_vector = _t20a.planted_vector
recover_d = _t20a.recover_d
norm = _t20a.norm
lam_star = _t20a.lam_star

SEEDS = [42, 1234, 9999, 555, 31337]
A5_M = 12   # signatures per curve in the A5 fresh-curve sweep

# Historical Phase-2 anchors (RESEARCH_AUTOLAB_LOG.md 2026-06-15 / 07-26 / 07-29)
HIST = [
    # label,              p,    b, n,    lam,  K1, m
    ("8-bit/199",         211,  2, 199,  106,  2,  6),
    ("12-bit/2557",       2557, 2, 2659, 1755, 8,  8),
    ("12-bit/2677 FAIL",  2677, 2, 2647, 185,  8,  10),
]


# ---------------------------------------------------------------------------
# centring: an orthogonal, one-line change to ANY of the formulations
# ---------------------------------------------------------------------------

def centre_sigs(sigs, n, lam, k1_bound, k2_bound):
    """Shift the inhomogeneous term so the unknowns are centred on 0.

    A_i + B_i*d = k1_i + lam*k2_i  (mod n)  with k1_i in [0,K1), k2_i in [0,K2).
    Put c1 = K1//2, c2 = K2//2 and A_i' = A_i - c1 - lam*c2 (mod n); then

        A_i' + B_i*d = (k1_i - c1) + lam*(k2_i - c2)   (mod n)

    with the unknowns now in [-K1/2, K1/2) and [-K2/2, K2/2).  E[x^2] falls
    from X^2/3 to X^2/12, so ||v_planted||^2 drops from ~n^2(2m/3 + 1) to
    ~n^2(m/6 + 1) -- a factor ~2 in norm, hence ~4 in the eff threshold.
    Nothing else about the lattice changes.
    """
    c1 = k1_bound // 2
    c2 = k2_bound // 2
    out = []
    for s in sigs:
        t = dict(s)
        t['A'] = (s['A'] - c1 - lam * c2) % n
        t['k1'] = s['k1'] - c1
        t['k2'] = s['k2'] - c2
        out.append(t)
    return out


# ---------------------------------------------------------------------------
# (P) projected lattice: same lattice, d-column quotiented out
# ---------------------------------------------------------------------------

def build_projected_lattice(sigs, n, lam, k1_bound, k2_bound):
    """dim 2m+1: [k1_0..k1_{m-1} | k2_0..k2_{m-1} | kannan].

    The d-column of the Phase-2 lattice is deleted.  The rank-m sublattice
    Lambda = {x in Z^m : x = d*B mod n for some d} that the d-row + modular
    rows generate is re-expressed by the triangular basis

        g_0 = (1, B_1/B_0, ..., B_{m-1}/B_0) mod n      [the d direction]
        g_i = n*e_i                        i = 1..m-1

    which has det n^{m-1} -- exactly det(orig)/(n*S_D), i.e. the projection
    along the trivial vector t = n*S_D*e_m.  S_D drops out entirely.
    """
    m = len(sigs)
    dim = 2 * m + 1
    S_K1, _S_D, S_K2, S_KANNAN = scales(n, k1_bound, k2_bound)
    B0inv = modinv(sigs[0]['B'], n)

    M = [[0] * dim for _ in range(dim)]
    # row 0: the d direction, normalised so column 0 has coefficient 1
    for i in range(m):
        M[0][i] = (sigs[i]['B'] * B0inv % n) * S_K1
    # rows 1..m-1: modular
    for i in range(1, m):
        M[i][i] = n * S_K1
    # rows m..2m-1: GLV rows
    for i in range(m):
        M[m + i][i] = -lam * S_K1
        M[m + i][m + i] = S_K2
    # row 2m: Kannan row
    for i in range(m):
        M[2 * m][i] = sigs[i]['A'] * S_K1
    M[2 * m][dim - 1] = S_KANNAN
    return M


def projected_planted_vector(sigs, n, k1_bound, k2_bound):
    """(k1_i*S_K1, k2_i*S_K2, S_KANNAN) -- no d entry."""
    m = len(sigs)
    S_K1, _S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    v = [0] * (2 * m + 1)
    for i in range(m):
        v[i] = sigs[i]['k1'] * S_K1
        v[m + i] = sigs[i]['k2'] * S_K2
    v[2 * m] = S_KAN
    return v


def recover_d_projected(reduced, sigs, n, lam, k1_bound, k2_bound, d_secret):
    """Read (k1, k2) off a reduced row, then solve one congruence for d and
    verify it against all m signatures."""
    m = len(sigs)
    dim = 2 * m + 1
    S_K1, _S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    for row in reduced:
        if abs(row[dim - 1]) != S_KAN:
            continue
        sgn = 1 if row[dim - 1] > 0 else -1
        k1 = [sgn * row[i] for i in range(m)]
        k2 = [sgn * row[m + i] for i in range(m)]
        if any(x % S_K1 for x in k1) or any(x % S_K2 for x in k2):
            continue
        k1 = [x // S_K1 for x in k1]
        k2 = [x // S_K2 for x in k2]
        d = (k1[0] + lam * k2[0] - sigs[0]['A']) * modinv(sigs[0]['B'], n) % n
        if d == 0:
            continue
        if all((sigs[i]['A'] + sigs[i]['B'] * d - lam * k2[i] - k1[i]) % n == 0
               for i in range(m)):
            if d == d_secret:
                return d
    return None


# ---------------------------------------------------------------------------
# (C) CVP formulation: no Kannan column, explicit centred Babai target
# ---------------------------------------------------------------------------

def build_cvp_lattice(sigs, n, lam, k1_bound, k2_bound):
    """dim 2m: [k1 slots | k2 slots].  Same as (P) with the Kannan row/col
    removed; the inhomogeneous part A_i moves into the CVP target."""
    m = len(sigs)
    dim = 2 * m
    S_K1, _S_D, S_K2, _S_KAN = scales(n, k1_bound, k2_bound)
    B0inv = modinv(sigs[0]['B'], n)
    M = [[0] * dim for _ in range(dim)]
    for i in range(m):
        M[0][i] = (sigs[i]['B'] * B0inv % n) * S_K1
    for i in range(1, m):
        M[i][i] = n * S_K1
    for i in range(m):
        M[m + i][i] = -lam * S_K1
        M[m + i][m + i] = S_K2
    return M


def run_cvp(sigs, n, lam, k1_bound, k2_bound, d_secret, centred=True):
    """Babai nearest-plane on the LLL-reduced (C) basis.

    w - T = (k1_i*S_K1, k2_i*S_K2) for the true solution, where
    T = (-A_i*S_K1, 0).  Centring shifts the target by half a box so the
    residual is centred on 0 instead of living in the positive orthant.
    """
    m = len(sigs)
    dim = 2 * m
    S_K1, _S_D, S_K2, _S_KAN = scales(n, k1_bound, k2_bound)
    M = build_cvp_lattice(sigs, n, lam, k1_bound, k2_bound)
    A = IntegerMatrix.from_matrix(M)
    LLL.reduction(A)
    G = GSO.Mat(A)
    G.update_gso()

    T = [-sigs[i]['A'] * S_K1 for i in range(m)] + [0] * m
    if centred:
        c1 = (k1_bound // 2) * S_K1
        c2 = (k2_bound // 2) * S_K2
        tgt = [T[i] + c1 for i in range(m)] + [T[m + i] + c2 for i in range(m)]
    else:
        tgt = list(T)

    coeffs = G.babai([float(x) for x in tgt])
    w = [sum(coeffs[r] * A[r][c] for r in range(dim)) for c in range(dim)]

    k1 = [w[i] - T[i] for i in range(m)]
    k2 = [w[m + i] - T[m + i] for i in range(m)]
    if any(x % S_K1 for x in k1) or any(x % S_K2 for x in k2):
        return None
    k1 = [x // S_K1 for x in k1]
    k2 = [x // S_K2 for x in k2]
    d = (k1[0] + lam * k2[0] - sigs[0]['A']) * modinv(sigs[0]['B'], n) % n
    if d == 0:
        return None
    if all((sigs[i]['A'] + sigs[i]['B'] * d - lam * k2[i] - k1[i]) % n == 0
           for i in range(m)):
        return d if d == d_secret else None
    return None


# ---------------------------------------------------------------------------
# geometry: determinants and the Gaussian heuristic
# ---------------------------------------------------------------------------

def log_det_orig(m, n, k1_bound, k2_bound):
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    return (m * math.log(n * S_K1) + math.log(S_D)
            + m * math.log(S_K2) + math.log(S_KAN))


def log_det_proj(m, n, k1_bound, k2_bound):
    S_K1, _S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    return (m * math.log(S_K1) + (m - 1) * math.log(n)
            + m * math.log(S_K2) + math.log(S_KAN))


def gh(log_det, dim):
    """Gaussian heuristic for lambda_1."""
    return math.exp(0.5 * math.log(dim / (2 * math.pi * math.e))
                    + log_det / dim)


def expected_norm_orig(m, n, k1_bound, k2_bound, centred=False):
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    q = 12.0 if centred else 3.0
    return math.sqrt(m * (k1_bound * S_K1) ** 2 / q
                     + (n * S_D) ** 2 / 3.0
                     + m * (k2_bound * S_K2) ** 2 / q
                     + S_KAN ** 2)


def expected_norm_proj(m, n, k1_bound, k2_bound, centred=False):
    S_K1, _S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    q = 12.0 if centred else 3.0
    return math.sqrt(m * (k1_bound * S_K1) ** 2 / q
                     + m * (k2_bound * S_K2) ** 2 / q
                     + S_KAN ** 2)


# ---------------------------------------------------------------------------
# experiment drivers
# ---------------------------------------------------------------------------

def run_orig(curve, m, d_secret, k1_bound, seed, centred=False):
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return None
    if centred:
        sigs = centre_sigs(sigs, n, lam, k1_bound, k2_bound)
    M = build_glv_lattice(sigs, n, lam, k1_bound, k2_bound)
    dim = 2 * m + 2
    A = IntegerMatrix.from_matrix(M)
    LLL.reduction(A)
    reduced = [[A[i][j] for j in range(dim)] for i in range(dim)]
    _, _, _, S_KAN = scales(n, k1_bound, k2_bound)
    ok = recover_d(reduced, m, n, S_KAN, d_secret) is not None
    pv = norm(planted_vector(sigs, d_secret, n, k1_bound, k2_bound))
    norms = sorted(norm(r) for r in reduced if any(r))
    rank = sum(1 for x in norms if x < pv)
    return ok, pv, norms[0], rank


def run_proj(curve, m, d_secret, k1_bound, seed, centred=False):
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return None
    if centred:
        sigs = centre_sigs(sigs, n, lam, k1_bound, k2_bound)
    M = build_projected_lattice(sigs, n, lam, k1_bound, k2_bound)
    dim = 2 * m + 1
    A = IntegerMatrix.from_matrix(M)
    LLL.reduction(A)
    reduced = [[A[i][j] for j in range(dim)] for i in range(dim)]
    ok = recover_d_projected(reduced, sigs, n, lam, k1_bound, k2_bound,
                             d_secret) is not None
    pv = norm(projected_planted_vector(sigs, n, k1_bound, k2_bound))
    norms = sorted(norm(r) for r in reduced if any(r))
    rank = sum(1 for x in norms if x < pv)
    return ok, pv, norms[0], rank


def run_cvp_arm(curve, m, d_secret, k1_bound, seed, centred=True):
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return None
    return run_cvp(sigs, n, lam, k1_bound, k2_bound, d_secret,
                   centred=centred) is not None


def rate(fn, curve, m, k1_bound, seeds, **kw):
    wins, ratios, ranks = 0, [], []
    for seed in seeds:
        d_trial = random.Random(seed + 7777).randint(1, curve[2] - 1)
        res = fn(curve, m, d_trial, k1_bound, seed, **kw)
        if res is None:
            continue
        if isinstance(res, bool):
            wins += res
            continue
        ok, pv, sv, rk = res
        wins += bool(ok)
        ratios.append(sv / pv)
        ranks.append(rk)
    mr = sum(ratios) / len(ratios) if ratios else float('nan')
    mk = sum(ranks) / len(ranks) if ranks else float('nan')
    return wins, len(seeds), mr, mk


# ---------------------------------------------------------------------------

def main():
    print("=" * 78)
    print("Thread 23 -- make the planted vector lambda_1 (GLV-HNP Phase 2)")
    print("=" * 78)

    curves = []
    for label, p, b, n, lam, k1, m in HIST:
        G = find_generator(p, b, n)
        assert G is not None, f"no generator for {label}"
        assert (lam * lam + lam + 1) % n == 0, f"bad lam for {label}"
        curves.append((label, (p, b, n, lam, G), k1, m))

    # ---------------- A0: closed-form geometry ----------------
    print("\n" + "-" * 78)
    print("A0: geometry of ORIG (dim 2m+2) vs PROJ (dim 2m+1)")
    print("-" * 78)
    print("ORIG always contains the trivial vector t = n*S_D*e_m (norm n*S_D),")
    print("so lambda_1(ORIG) <= n and ||v||/lambda_1 > 1 by construction.")
    print("PROJ quotients Z*t out: det drops by exactly n*S_D, dim by 1,")
    print("and ||v||^2 drops by (n*S_D)^2/3.  S_D cancels from det(PROJ).\n")
    hdr = (f"{'curve':<18} {'K1':>3} {'m':>3} {'||v||o/GHo':>10} "
           f"{'||v||p/GHp':>10} {'l1(orig)=n':>11} {'||v||o/n':>9}")
    print(hdr)
    for label, cv, k1, m in curves:
        p, b, n, lam, G = cv
        k2 = math.isqrt(n) + 1
        vo = expected_norm_orig(m, n, k1, k2)
        vp = expected_norm_proj(m, n, k1, k2)
        gho = gh(log_det_orig(m, n, k1, k2), 2 * m + 2)
        ghp = gh(log_det_proj(m, n, k1, k2), 2 * m + 1)
        print(f"{label:<18} {k1:>3} {m:>3} {vo/gho:>10.3f} {vp/ghp:>10.3f} "
              f"{n:>11} {vo/n:>9.3f}")

    # ---------------- A1: correctness ----------------
    print("\n" + "-" * 78)
    print("A1: correctness of the two reformulations on the anchors")
    print("-" * 78)
    print("ORIG/PROJ = Kannan embedding, d-column kept / quotiented out.")
    print("suffix c = centred inhomogeneous term.  CVP = Babai nearest-plane.\n")
    print(f"{'curve':<18} {'K1':>3} {'m':>3} {'ORIG':>6} {'ORIGc':>6} "
          f"{'PROJ':>6} {'PROJc':>6} {'CVPc':>6} {'CVPu':>6}")
    for label, cv, k1, m in curves:
        ns = len(SEEDS)
        wo = rate(run_orig, cv, m, k1, SEEDS)[0]
        woc = rate(run_orig, cv, m, k1, SEEDS, centred=True)[0]
        wp = rate(run_proj, cv, m, k1, SEEDS)[0]
        wpc = rate(run_proj, cv, m, k1, SEEDS, centred=True)[0]
        wc = rate(run_cvp_arm, cv, m, k1, SEEDS, centred=True)[0]
        wu = rate(run_cvp_arm, cv, m, k1, SEEDS, centred=False)[0]
        print(f"{label:<18} {k1:>3} {m:>3} {wo:>4}/{ns} {woc:>4}/{ns} "
              f"{wp:>4}/{ns} {wpc:>4}/{ns} {wc:>4}/{ns} {wu:>4}/{ns}")

    # ---------------- A2: is the planted vector lambda_1 now? ----------------
    print("\n" + "-" * 78)
    print("A2: sv/pv and rank of v_planted among reduced basis vectors")
    print("-" * 78)
    print("sv/pv = ||shortest reduced row|| / ||v_planted||.  rank = how many")
    print("reduced rows are strictly shorter than v_planted (0 => v is lambda_1).\n")
    print(f"{'curve':<18} {'K1':>3} {'sv/pv ORIG':>11} {'rk ORIG':>8} "
          f"{'sv/pv PROJ':>11} {'rk PROJ':>8} {'sv/pv PROJc':>12} "
          f"{'rk PROJc':>9}")
    for label, cv, k1, m in curves:
        _, _, ro, ko = rate(run_orig, cv, m, k1, SEEDS)
        _, _, rp, kp = rate(run_proj, cv, m, k1, SEEDS)
        _, _, rpc, kpc = rate(run_proj, cv, m, k1, SEEDS, centred=True)
        print(f"{label:<18} {k1:>3} {ro:>11.3f} {ko:>8.1f} "
              f"{rp:>11.3f} {kp:>8.1f} {rpc:>12.3f} {kpc:>9.1f}")

    # ---------------- A3: the K1 wall (the falsifier) ----------------
    print("\n" + "-" * 78)
    print("A3: the T4 K1 wall -- ORIG vs PROJ vs CVP (5 seeds each)")
    print("-" * 78)
    print("2026-07-29 T4 measured, for ORIG: 12-bit/2557 walls at K1~12-16,")
    print("12-bit/2677 walls at K1~4-6.  Falsifier: does PROJ/CVP move it out?\n")
    K1S = [2, 3, 4, 6, 8, 12, 16, 24, 32, 48]
    for label, cv, _k1, m in curves[1:]:
        p, b, n, lam, G = cv
        print(f"\n{label}   (lam* = {lam_star(lam, n):.3f}, m = {m}, "
              f"n = {n})")
        head = "  " + "".join(f"{k:>6}" for k in K1S)
        print(f"{'K1':<8}{head}")
        for name, fn, kw in (("ORIG", run_orig, {}),
                             ("ORIGc", run_orig, {'centred': True}),
                             ("PROJ", run_proj, {}),
                             ("PROJc", run_proj, {'centred': True}),
                             ("CVPc", run_cvp_arm, {'centred': True})):
            cells = []
            for k1 in K1S:
                w, t, _, _ = rate(fn, cv, m, k1, SEEDS, **kw)
                cells.append(f"{w}/{t}")
            print(f"{name:<8}  " + "".join(f"{c:>6}" for c in cells))

    # ---------------- A4: does GH ratio = 1 predict the wall? ----------------
    print("\n" + "-" * 78)
    print("A4: Gaussian-heuristic prediction of the wall location")
    print("-" * 78)
    print("PROJ:  ||v||^2 ~ n^2(2m/3 + 1),  det = S_K1^m n^(m-1) S_K2^m S_KAN.")
    print("Asymptotically ||v||/GH -> sqrt(2*pi*e/3) * sqrt(eff) = 2.38*sqrt(eff)")
    print("with eff = K1*K2/n, so the predicted uncentred wall is")
    print(f"  eff = 3/(2*pi*e)  = {3/(2*math.pi*math.e):.4f}")
    print("and centring (E[x^2]: X^2/3 -> X^2/12) moves it out by 4x:")
    print(f"  eff = 12/(2*pi*e) = {12/(2*math.pi*math.e):.4f}\n")
    print(f"{'curve':<18} {'K1':>4} {'eff':>7} {'||v||/GH':>9} "
          f"{'PROJ obs':>9} {'||v||c/GH':>10} {'PROJc obs':>10}")
    for label, cv, _k1, m in curves[1:]:
        p, b, n, lam, G = cv
        k2 = math.isqrt(n) + 1
        for k1 in K1S:
            eff = k1 * k2 / n
            ghp = gh(log_det_proj(m, n, k1, k2), 2 * m + 1)
            vp = expected_norm_proj(m, n, k1, k2)
            vc = expected_norm_proj(m, n, k1, k2, centred=True)
            w, t, _, _ = rate(run_proj, cv, m, k1, SEEDS)
            wc, tc, _, _ = rate(run_proj, cv, m, k1, SEEDS, centred=True)
            print(f"{label:<18} {k1:>4} {eff:>7.3f} {vp/ghp:>9.3f} "
                  f"{str(w)+'/'+str(t):>9} {vc/ghp:>10.3f} "
                  f"{str(wc)+'/'+str(tc):>10}")

    # ---------------- A5: does centring move the eff ceiling generally? ------
    print("\n" + "-" * 78)
    print("A5: fresh 17-bit curves -- the eff ceiling, uncentred vs centred")
    print("-" * 78)
    print("2026-07-29 T3 measured (uncentred, ORIG): eff 0.05 -> 19/20 curves")
    print("recovering, eff 0.15 -> 3/20, eff 0.25 -> 0/20.  If centring is the")
    print("real variable the ceiling should move out by ~4x in eff.\n")
    print(f"Searching j=0 GLV curves with p in [2^16, 2^17) ...")
    fresh = _t20a.search_curves(2 ** 16, 2 ** 17, per_bin=1, nbins=10)
    print(f"Found {len(fresh)} curves.  m = {A5_M}, {len(SEEDS)} seeds each.\n")
    print(f"{'eff':>6} {'PROJ (uncentred)':>18} {'PROJc (centred)':>17} "
          f"{'lam* of PROJc wins':>20}")
    for eff_t in (0.05, 0.15, 0.25, 0.40, 0.60, 0.80):
        nu = nc = 0
        ls_win = []
        for (p, b, n, lam, G) in fresh:
            k2 = math.isqrt(n) + 1
            k1 = max(2, int(eff_t * n / k2))
            cv = (p, b, n, lam, G)
            wu, t, _, _ = rate(run_proj, cv, A5_M, k1, SEEDS)
            wc, _, _, _ = rate(run_proj, cv, A5_M, k1, SEEDS, centred=True)
            nu += (wu == t)
            nc += (wc == t)
            if wc == t:
                ls_win.append(lam_star(lam, n))
        rng = (f"[{min(ls_win):.3f}, {max(ls_win):.3f}]" if ls_win else "-")
        print(f"{eff_t:>6.2f} {str(nu)+'/'+str(len(fresh)):>18} "
              f"{str(nc)+'/'+str(len(fresh)):>17} {rng:>20}")

    print("\n" + "=" * 78)
    print("done")
    print("=" * 78)


if __name__ == "__main__":
    main()
