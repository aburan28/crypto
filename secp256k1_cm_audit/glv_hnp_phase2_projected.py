"""
GLV-HNP Phase 2, Thread 23: reformulate the lattice so the target IS lambda_1.

Background (RESEARCH_AUTOLAB_LOG.md, 2026-07-29 run, commit e845207, exp T5):
  The Phase-2 Kannan lattice always contains the trivial vector n*S_D*e_m
  (obtained as n*row_m - sum_i B_i*row_i), of norm exactly n*S_D.  The planted
  vector has norm ~ n*sqrt(2m/3 + 4/3) > n*S_D, so the planted vector is NEVER
  lambda_1: measured sv/pv in [0.34, 0.61] on every curve, success and failure
  alike.  Recovery is therefore a BDD/coset condition, not an SVP condition.
  No choice of S_D removes the trivial vector (both norms scale linearly in S_D).

  Proposed falsifier (2026-07-29 "Next step proposal"):
    Project the lattice along e_m (quotient out the trivial direction), or
    replace the Kannan embedding with an explicit CVP.  If sv/pv rises above 1
    AND the K1 wall on the lam*=0.07 curve (currently K1 ~ 4-6) moves outward,
    the reformulation is a real improvement.  If the wall stays at K1 ~ 4-6,
    the wall is information-theoretic and Phase 2 is at its ceiling.

This script implements that falsifier.

Three lattice variants, all fed bit-identical signatures (same RNG seeds, same
gen_signatures, same scales) so the K1 grids are directly comparable to the
2026-07-29 T4 table:

  V0  BASELINE      the 2026-06-15 Kannan lattice, dim 2m+2.  Reproduces T4.
                    Recovery reads d straight out of the d-column.

  V1  PROJECTED     V0 with column m (the d-column) DELETED.  dim 2m+2 rows in
                    2m+1 columns, rank 2m+1: the trivial vector n*S_D*e_m maps
                    to 0, so it is gone by construction.  d is no longer a
                    lattice coordinate; it is recovered algebraically from the
                    (k1_0, k2_0) block of the short vector via
                        d = (k1_0 + lam*k2_0 - A_0) * B_0^{-1}  (mod n).

  V2  CVP           d eliminated AND the Kannan embedding dropped.  Lattice L2
                    of rank 2m spanned by
                        r_i  = n*S_K1*e_i                       (i < m)
                        r_d  = (B_i*S_K1)_i                     (k1-block only)
                        r'_i = -lam*S_K1*e_i + S_K2*e_{m+i}
                    with EXACT CVP (fpylll enumeration) to the target
                        t = (-A_i*S_K1, 0, ..., 0).
                    The closest vector w satisfies w - t = (k1_i*S_K1, k2_i*S_K2)
                    exactly when the attack succeeds.  This is the sharpest form
                    of "quotient out the trivial direction": no coset condition
                    and no embedding slack at all.

Recovery is verified WITHOUT the secret: a candidate d is accepted iff d*G == Q,
the public key.  The script then asserts the accepted d == d_secret, so a "win"
is a genuine attack rather than an oracle query.  (The 2026-07-29 criterion,
`recover_d`, compared the d-column directly to d_secret; the public-key test is
the honest form of the same check and gives identical counts.)

DO NOT tighten this to "the row equals the planted vector" — see the note on
mu-aliases in d_from_block() below.  That mistake was made and corrected during
the 2026-08-03 run; it under-counts V0 by ~2 K1-steps.

Experiments:
  E1  sv/pv for V0 vs V1 on the three historical curves: does the planted
      vector become lambda_1?
  E2  correctness of V1/V2 recovery in the easy regime (K1=2).
  E3  THE FALSIFIER: K1 wall, V0 vs V1 vs V2, on the lam*=0.34 and lam*=0.07
      historical curves.  Compare against the 2026-07-29 T4 walls.
  E4  m-sweep at the wall (does more data rescue V1/V2 where it did not
      rescue V0?  T4b: V0 at K1=8 got 0,0,1,0,1 of 5 for m=8/12/16/24/32).

Run: python3 glv_hnp_phase2_projected.py
"""

import math
import time
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from fpylll import CVP, LLL, BKZ, IntegerMatrix

from glv_hnp_lib import (
    build_glv_lattice,
    lambda_block_mu,
    search_curves,
    ec_mul,
    find_generator,
    gen_signatures,
    lam_star,
    modinv,
    norm,
    planted_vector,
    scales,
)

SEEDS = [42, 1234, 9999, 555, 31337]

# Historical Phase-2 curves (RESEARCH_AUTOLAB_LOG.md 2026-06-15 / 2026-07-26 /
# 2026-07-29).  Same tuples as glv_hnp_phase2_lambda_threshold.py:426.
HIST = [
    # label,             p,    b, n,    lam,  K1, m
    ("8-bit/199",        211,  2, 199,  106,  2,  6),
    ("12-bit/2557",      2557, 2, 2659, 1755, 8,  8),
    ("12-bit/2677 FAIL", 2677, 2, 2647, 185,  8,  10),
]


# ---------------------------------------------------------------------------
# V1 — projected lattice (d-column deleted)
# ---------------------------------------------------------------------------

def build_projected_lattice(sigs, n, lam, k1_bound, k2_bound):
    """V0's lattice with column m removed.  Columns are
    0..m-1 = k1 block, m..2m-1 = k2 block, 2m = Kannan."""
    M = build_glv_lattice(sigs, n, lam, k1_bound, k2_bound)
    m = len(sigs)
    return [row[:m] + row[m + 1:] for row in M]


def project(v, m):
    """Drop coordinate m of a dim-(2m+2) vector."""
    return v[:m] + v[m + 1:]


# ---------------------------------------------------------------------------
# V2 — d-eliminated lattice for exact CVP
# ---------------------------------------------------------------------------

def build_cvp_lattice(sigs, n, lam, k1_bound, k2_bound):
    """Rank-2m lattice; columns 0..m-1 = k1 block, m..2m-1 = k2 block."""
    m = len(sigs)
    S_K1, _, S_K2, _ = scales(n, k1_bound, k2_bound)
    rows = []
    for i in range(m):                                    # r_i = n*S_K1*e_i
        r = [0] * (2 * m)
        r[i] = n * S_K1
        rows.append(r)
    r_d = [sigs[i]['B'] * S_K1 for i in range(m)] + [0] * m   # r_d
    rows.append(r_d)
    for i in range(m):                                    # r'_i
        r = [0] * (2 * m)
        r[i] = -(lam % n) * S_K1
        r[m + i] = S_K2
        rows.append(r)
    return rows


def cvp_target(sigs, n, k1_bound, k2_bound):
    m = len(sigs)
    S_K1, _, _, _ = scales(n, k1_bound, k2_bound)
    return [-sigs[i]['A'] * S_K1 for i in range(m)] + [0] * m


# ---------------------------------------------------------------------------
# Secret-free verification
# ---------------------------------------------------------------------------

def verify_d(d_cand, pub):
    """Secret-free acceptance test: does d_cand reproduce the public key?
    `pub` is (p, n, G, Q) with Q = d_secret*G, all public data."""
    p, n, G, Q = pub
    if not (0 < d_cand < n):
        return False
    return ec_mul(G, d_cand, p) == Q


def d_from_block(k1_0, k2_0, sigs, n, lam):
    """d = (k1_0 + lam*k2_0 - A_0) * B_0^{-1} mod n.

    NOTE (2026-08-03): this is alias-robust.  LLL frequently returns not the
    planted vector but a mu-ALIAS of it: planted + (dk1_i, dk2_i) where
    (dk1, dk2) is a short vector of the lambda-block sublattice
    <(n*S_K1,0), (-lam*S_K1,S_K2)>, hence dk1 + lam*dk2 == 0 (mod n).  The
    alias has k2_i outside [0,K2) but leaves k1_i + lam*k2_i, and therefore d,
    unchanged.  Any recovery test that insists on the exact planted (k1,k2)
    rejects these and badly under-counts the attack's success rate.
    """
    B0 = sigs[0]['B'] % n
    if B0 == 0:
        return None
    return (k1_0 + lam * k2_0 - sigs[0]['A']) * modinv(B0, n) % n


def extract_from_rows(rows, sigs, n, lam, k1_bound, k2_bound, pub, kannan_col,
                      d_col=None, k2_off=None):
    """Scan reduced rows for one with |last| == S_KANNAN and read d off it.

    d_col is V0's explicit d-column (None for the projected lattice V1).
    k2_off is the first k2 column.
    """
    m = len(sigs)
    S_K1, _, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    for row in rows:
        last = row[kannan_col]
        if abs(last) != S_KAN:
            continue
        r = [-x for x in row] if last < 0 else list(row)
        if d_col is not None:
            d = r[d_col] % n
            if verify_d(d, pub):
                return d
        if r[0] % S_K1 or r[k2_off] % S_K2:
            continue
        d = d_from_block(r[0] // S_K1, r[k2_off] // S_K2, sigs, n, lam)
        if d is not None and verify_d(d, pub):
            return d
    return None


# ---------------------------------------------------------------------------
# The three attacks
# ---------------------------------------------------------------------------

def attack_v0(sigs, n, lam, k1_bound, k2_bound, pub, use_bkz=False, beta=20):
    m = len(sigs)
    dim = 2 * m + 2
    M = build_glv_lattice(sigs, n, lam, k1_bound, k2_bound)
    A = IntegerMatrix.from_matrix(M)
    if use_bkz:
        BKZ.reduction(A, BKZ.Param(beta))
    else:
        LLL.reduction(A)
    rows = [[A[i][j] for j in range(dim)] for i in range(dim)]
    d = extract_from_rows(rows, sigs, n, lam, k1_bound, k2_bound, pub,
                          kannan_col=2 * m + 1, d_col=m, k2_off=m + 1)
    sv = min((norm(r) for r in rows if any(r)), default=float('nan'))
    return d, sv, rows


def attack_v1(sigs, n, lam, k1_bound, k2_bound, pub, use_bkz=False, beta=20):
    m = len(sigs)
    ncols = 2 * m + 1
    M = build_projected_lattice(sigs, n, lam, k1_bound, k2_bound)
    A = IntegerMatrix.from_matrix(M)
    if use_bkz:
        BKZ.reduction(A, BKZ.Param(beta))
    else:
        LLL.reduction(A)
    rows = [[A[i][j] for j in range(ncols)] for i in range(len(M))]
    d = extract_from_rows(rows, sigs, n, lam, k1_bound, k2_bound, pub,
                          kannan_col=2 * m, d_col=None, k2_off=m)
    sv = min((norm(r) for r in rows if any(r)), default=float('nan'))
    return d, sv, rows


def attack_v2(sigs, n, lam, k1_bound, k2_bound, pub):
    """Exact CVP (fpylll enumeration) on the d-eliminated lattice."""
    m = len(sigs)
    S_K1, _, S_K2, _ = scales(n, k1_bound, k2_bound)
    rows = build_cvp_lattice(sigs, n, lam, k1_bound, k2_bound)
    t = cvp_target(sigs, n, k1_bound, k2_bound)
    A = IntegerMatrix.from_matrix(rows)
    LLL.reduction(A)
    # build_cvp_lattice emits 2m+1 generators of a rank-2m lattice; LLL pushes
    # the dependency into leading zero rows.  fplll's CVP needs a square basis.
    basis = [[A[i][j] for j in range(2 * m)] for i in range(len(rows))]
    basis = [r for r in basis if any(r)]
    if len(basis) != 2 * m:
        return None, float('nan'), f"rank {len(basis)} != {2 * m}"
    A = IntegerMatrix.from_matrix(basis)
    try:
        w = CVP.closest_vector(A, tuple(t))
    except Exception as exc:                       # enumeration failure
        return None, float('nan'), str(exc)
    diff = [w[j] - t[j] for j in range(2 * m)]
    if diff[0] % S_K1 or diff[m] % S_K2:
        return None, norm(diff), None
    d = d_from_block(diff[0] // S_K1, diff[m] // S_K2, sigs, n, lam)
    if d is None or not verify_d(d, pub):
        return None, norm(diff), None
    return d, norm(diff), None


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def make_sigs(curve, m, d_secret, k1_bound, seed):
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed)
    return (sigs, k2_bound) if len(sigs) == m else (None, k2_bound)


def trial(curve, m, k1_bound, seed, variant, beta=0):
    """One instance.  Returns (won: bool, aux) with aux variant-specific."""
    p, b, n, lam, G = curve
    d_secret = random.Random(seed + 7777).randint(1, n - 1)
    sigs, k2_bound = make_sigs(curve, m, d_secret, k1_bound, seed)
    if sigs is None:
        return None
    pub = (p, n, G, ec_mul(G, d_secret, p))
    use_bkz = beta > 0
    if variant == 'V0':
        d, sv, _ = attack_v0(sigs, n, lam, k1_bound, k2_bound, pub, use_bkz, beta)
    elif variant == 'V1':
        d, sv, _ = attack_v1(sigs, n, lam, k1_bound, k2_bound, pub, use_bkz, beta)
    elif variant == 'V2':
        d, sv, _ = attack_v2(sigs, n, lam, k1_bound, k2_bound, pub)
    else:
        raise ValueError(variant)
    if d is not None:
        # secret-free verification already passed; this must hold.
        assert d == d_secret, f"{variant}: verified d != secret ({d} vs {d_secret})"
    return bool(d is not None), sv


def rate(curve, m, k1_bound, variant, seeds=SEEDS, beta=0):
    wins = 0
    tot = 0
    for s in seeds:
        r = trial(curve, m, k1_bound, s, variant, beta)
        if r is None:
            continue
        tot += 1
        wins += r[0]
    return wins, tot


def load_hist():
    out = []
    for label, p, b, n, lam, k1, m in HIST:
        G = find_generator(p, b, n)
        assert G is not None, f"no generator for {label}"
        assert (lam * lam + lam + 1) % n == 0, f"bad lam for {label}"
        out.append((label, (p, b, n, lam, G), k1, m))
    return out


def main():
    print("=" * 78)
    print("Thread 23 — reformulating the Phase-2 lattice so the target is lambda_1")
    print("=" * 78)

    hist = load_hist()


    # ===========================================================================
    # E1 — is the planted vector lambda_1 after projection?
    # ===========================================================================
    print("\n" + "-" * 78)
    print("E1: shortest-vector / planted-vector ratio, V0 (Kannan) vs V1 (projected)")
    print("-" * 78)
    print("2026-07-29 T5 measured sv/pv in [0.34, 0.61] for V0 and showed the")
    print("minimiser is always the trivial vector n*S_D*e_m.  If projection works,")
    print("V1's sv/pv should rise to 1.0 (planted vector IS lambda_1).\n")
    print(f"{'curve':<18} {'lam*':>6} {'K1':>3} {'m':>3} "
          f"{'pv(V0)':>10} {'sv(V0)':>10} {'r0':>6} "
          f"{'pv(V1)':>10} {'sv(V1)':>10} {'r1':>6}")

    e1_rows = []
    for label, curve, k1, m in hist:
        p, b, n, lam, G = curve
        d_secret = random.Random(42 + 7777).randint(1, n - 1)
        sigs, k2b = make_sigs(curve, m, d_secret, k1, 42)
        pv = planted_vector(sigs, d_secret, n, k1, k2b)
        pv0 = norm(pv)
        pv1 = norm(project(pv, m))
        pub = (p, n, G, ec_mul(G, d_secret, p))
        _, sv0, _ = attack_v0(sigs, n, lam, k1, k2b, pub)
        _, sv1, _ = attack_v1(sigs, n, lam, k1, k2b, pub)
        r0, r1 = sv0 / pv0, sv1 / pv1
        e1_rows.append((label, r0, r1))
        print(f"{label:<18} {lam_star(lam, n):>6.3f} {k1:>3} {m:>3} "
              f"{pv0:>10.1f} {sv0:>10.1f} {r0:>6.3f} "
              f"{pv1:>10.1f} {sv1:>10.1f} {r1:>6.3f}")

    print("\n(r = sv/pv.  r < 1 means some OTHER vector is shorter than the planted")
    print(" one, i.e. the planted vector is not lambda_1 and recovery is BDD.)")


    # ===========================================================================
    # E2 — correctness in the easy regime
    # ===========================================================================
    print("\n" + "-" * 78)
    print("E2: sanity — all three variants in the easy regime (K1 = 2)")
    print("-" * 78)
    print("Secret-free verification: a win means the recovered d satisfies")
    print("d*G == Q for the public key Q.  d == d_secret is then asserted.\n")
    print(f"{'curve':<18} {'m':>3} {'K1':>3} {'V0':>7} {'V1':>7} {'V2':>7}")
    for label, curve, k1, m in hist:
        res = [rate(curve, m, 2, v) for v in ('V0', 'V1', 'V2')]
        cells = "".join(f"{str(w) + '/' + str(t):>8}" for w, t in res)
        print(f"{label:<18} {m:>3} {2:>3}{cells}")


    # ===========================================================================
    # E3 — THE FALSIFIER: does the K1 wall move?
    # ===========================================================================
    print("\n" + "-" * 78)
    print("E3: K1 wall, V0 vs V1 vs V2  (THE FALSIFIER)")
    print("-" * 78)
    print("2026-07-29 T4 walls (V0):  lam*=0.340 curve -> K1 ~ 12-16")
    print("                           lam*=0.070 curve -> K1 ~  4-6")
    print("If V1/V2 push the wall outward, the reformulation is a real improvement.")
    print("If the wall does not move, it is information-theoretic.\n")

    K1_GRID = [2, 3, 4, 6, 8, 12, 16, 24]
    e3 = {}
    for label, curve, _k1, m in hist[1:]:            # the two 12-bit curves
        p, b, n, lam, G = curve
        print(f"\ncurve {label}   n={n}  lam={lam}  lam*={lam_star(lam, n):.3f}  m={m}")
        hdr = "".join(f"{k:>7}" for k in K1_GRID)
        print(f"{'variant':<9}{hdr}")
        for v in ('V0', 'V1', 'V2'):
            cells = []
            for k1 in K1_GRID:
                w, t = rate(curve, m, k1, v)
                cells.append(f"{str(w) + '/' + str(t):>7}")
                e3[(label, v, k1)] = (w, t)
            print(f"{v:<9}" + "".join(cells))


    # ===========================================================================
    # E4 — does more data rescue V1/V2 at the wall?
    # ===========================================================================
    print("\n" + "-" * 78)
    print("E4: m-sweep at K1 = 8 on the lam*=0.07 curve")
    print("-" * 78)
    print("2026-07-29 T4b (V0): m = 8/12/16/24/32 -> 0,0,1,0,1 out of 5.")
    print("If V1/V2 convert this to consistent wins, the K1 wall was an artefact")
    print("of the embedding rather than of the data.\n")

    label, curve, _k1, _m = hist[2]
    M_GRID = [8, 12, 16, 24, 32]
    print(f"{'variant':<9}" + "".join(f"{'m=' + str(mm):>8}" for mm in M_GRID))
    for v in ('V0', 'V1', 'V2'):
        cells = []
        for mm in M_GRID:
            w, t = rate(curve, mm, 8, v)
            cells.append(f"{str(w) + '/' + str(t):>8}")
        print(f"{v:<9}" + "".join(cells))

    # =======================================================================
    # E5 — CONTROL: is V2's gain reformulation, or just stronger reduction?
    # =======================================================================
    print("\n" + "-" * 78)
    print("E5: control — BKZ on V0/V1 vs exact CVP on V2, lam*=0.07 curve, K1=8")
    print("-" * 78)
    print("E4 showed V2 winning 5/5 at m=24,32 where V0/V1 win 0-1/5.  If BKZ at")
    print("large beta closes that gap, the gain is better reduction and NOT the")
    print("reformulation.  If it does not, the d-elimination is doing the work.\n")

    label, curve, _k1, _m = hist[2]
    print(f"{'variant':<14}" + "".join(f"{'m=' + str(mm):>8}" for mm in M_GRID))
    for name, variant, beta in [("V0 LLL", 'V0', 0), ("V0 BKZ-20", 'V0', 20),
                                ("V0 BKZ-40", 'V0', 40), ("V1 BKZ-40", 'V1', 40),
                                ("V2 CVP", 'V2', 0)]:
        cells = []
        for mm in M_GRID:
            w, t = rate(curve, mm, 8, variant, beta=beta)
            cells.append(f"{str(w) + '/' + str(t):>8}")
        print(f"{name:<14}" + "".join(cells))

    # =======================================================================
    # E6 — cost of exact CVP, and what lambda_1 actually is after projection
    # =======================================================================
    print("\n" + "-" * 78)
    print("E6: cost of exact CVP, and the identity of lambda_1 in V1")
    print("-" * 78)
    print("V2 uses fplll enumeration, which is exponential in the lattice rank 2m.")
    print("Median wall-clock per instance, lam*=0.07 curve, K1=8:\n")
    print(f"{'m':>4} {'rank':>5} {'V0 LLL s':>10} {'V2 CVP s':>10} {'ratio':>8}")
    for mm in M_GRID:
        tv = {}
        for variant in ('V0', 'V2'):
            ts = []
            for sd in SEEDS[:3]:
                t0 = time.perf_counter()
                trial(curve, mm, 8, sd, variant)
                ts.append(time.perf_counter() - t0)
            ts.sort()
            tv[variant] = ts[len(ts) // 2]
        ratio = tv['V2'] / tv['V0'] if tv['V0'] else float('nan')
        print(f"{mm:>4} {2 * mm:>5} {tv['V0']:>10.3f} {tv['V2']:>10.3f} {ratio:>8.1f}x")

    print("\nEnergy split of the V1 shortest vector (projection removed the")
    print("trivial vector n*S_D*e_m, yet E1 still shows sv/pv < 1 — so what is")
    print("shorter than the planted vector now?):\n")
    print(f"{'curve':<18} {'sv/pv':>7} {'k1-blk':>8} {'k2-blk':>8} {'kannan':>8}")
    for label, curve, k1, m in hist:
        p_, b_, n_, lam_, G_ = curve
        d_secret = random.Random(42 + 7777).randint(1, n_ - 1)
        sigs, k2b = make_sigs(curve, m, d_secret, k1, 42)
        pub = (p_, n_, G_, ec_mul(G_, d_secret, p_))
        _, sv, rows = attack_v1(sigs, n_, lam_, k1, k2b, pub)
        best = min((r for r in rows if any(r)), key=norm)
        e = sum(x * x for x in best)
        k1e = sum(best[i] ** 2 for i in range(m)) / e
        k2e = sum(best[m + i] ** 2 for i in range(m)) / e
        kae = best[2 * m] ** 2 / e
        pv1 = norm(project(planted_vector(sigs, d_secret, n_, k1, k2b), m))
        print(f"{label:<18} {sv / pv1:>7.3f} {k1e:>8.3f} {k2e:>8.3f} {kae:>8.3f}")

    # =======================================================================
    # E7 — confirmation at 10 seeds, and does V2 move the K1 WALL itself?
    # =======================================================================
    print("\n" + "-" * 78)
    print("E7: 10-seed confirmation, and the K1 wall at large m")
    print("-" * 78)
    print("E4/E5 used 5 seeds.  Re-run the decisive cells at 10 seeds, then ask")
    print("the sharper question: with m=32 signatures, how far out in K1 can V2")
    print("go?  V0's wall on this curve is K1 ~ 4-6 (2026-07-29 T4).\n")

    SEEDS10 = [42, 1234, 9999, 555, 31337, 7, 20260803, 131071, 99, 424242]
    label, curve, _k1, _m = hist[2]

    print("(a) K1=8, m-sweep, 10 seeds")
    print(f"{'variant':<14}" + "".join(f"{'m=' + str(mm):>9}" for mm in M_GRID))
    for name, variant, beta in [("V0 LLL", 'V0', 0), ("V0 BKZ-40", 'V0', 40),
                                ("V2 CVP", 'V2', 0)]:
        cells = []
        for mm in M_GRID:
            w, t = rate(curve, mm, 8, variant, seeds=SEEDS10, beta=beta)
            cells.append(f"{str(w) + '/' + str(t):>9}")
        print(f"{name:<14}" + "".join(cells))

    print("\n(b) m=32, K1-sweep, 10 seeds  — does the WALL move outward?")
    K1_BIG = [4, 6, 8, 12, 16, 24, 32]
    print(f"{'variant':<14}" + "".join(f"{'K1=' + str(k):>9}" for k in K1_BIG))
    for name, variant, beta in [("V0 LLL", 'V0', 0), ("V0 BKZ-40", 'V0', 40),
                                ("V2 CVP", 'V2', 0)]:
        cells = []
        for k1 in K1_BIG:
            w, t = rate(curve, 32, k1, variant, seeds=SEEDS10, beta=beta)
            cells.append(f"{str(w) + '/' + str(t):>9}")
        print(f"{name:<14}" + "".join(cells))

    print("\n(c) information-theoretic reference: the attack needs")
    print("    m*log2(K1*K2) + log2(n) <= m*log2(n), i.e. m >= log(n)/log(n/(K1*K2)).")
    p_, b_, n_, lam_, G_ = curve
    k2b = math.isqrt(n_) + 1
    print(f"    n={n_}, K2={k2b}")
    print(f"{'K1':>5} {'eff=K1*K2/n':>13} {'m_min':>8}")
    for k1 in K1_BIG:
        eff = k1 * k2b / n_
        mmin = math.ceil(math.log(n_) / math.log(1.0 / eff)) if eff < 1 else float('inf')
        print(f"{k1:>5} {eff:>13.3f} {str(mmin):>8}")

    # =======================================================================
    # E8 — a BDD decoding-radius predictor for the residual V2 wall
    # =======================================================================
    print("\n" + "-" * 78)
    print("E8: does dist/GH predict the residual V2 wall?")
    print("-" * 78)
    print("E7(c) shows K1=12 needs only m_min=6 signatures, yet V2 with m=32")
    print("still fails there — so the residual wall is NOT information-theoretic.")
    print("Test the lattice-geometry explanation instead.  For the V2 lattice L2")
    print("(rank 2m) the determinant is exactly")
    print("      det(L2) = S_K1^m * n^(m-1) * S_K2^m        (verified symbolically)")
    print("so the Gaussian heuristic is GH = det^(1/2m) * sqrt(2m/(2*pi*e)) and the")
    print("BDD target distance is ||(k1_i*S_K1, k2_i*S_K2)|| ~ n*sqrt(2m/3).")
    print("Closed form:  dist/GH ~ sqrt(2*pi*e/3) * sqrt(K1*K2/n) * n^(1/2m).\n")

    mm = 32
    for clabel, ccurve, _ck1, _cm in hist[1:]:
        p_, b_, n_, lam_, G_ = ccurve
        k2b = math.isqrt(n_) + 1
        print(f"curve {clabel}, lam*={lam_star(lam_, n_):.3f}, m={mm}, "
              f"rank={2 * mm}, n={n_}, K2={k2b}")
        print(f"{'K1':>4} {'eff':>7} {'GH':>12} {'dist(emp)':>12} {'dist/GH':>9} "
              f"{'V2 10-seed':>11} {'V0 10-seed':>11}")
        for k1 in K1_BIG:
            S_K1, _, S_K2, _ = scales(n_, k1, k2b)
            logdet = (mm * math.log(S_K1) + (mm - 1) * math.log(n_)
                      + mm * math.log(S_K2))
            gh = (math.exp(logdet / (2 * mm))
                  * math.sqrt(2 * mm / (2 * math.pi * math.e)))
            ds = []
            for sd in SEEDS10[:5]:
                dsec = random.Random(sd + 7777).randint(1, n_ - 1)
                sg, _ = make_sigs(ccurve, mm, dsec, k1, sd)
                if sg is None:
                    continue
                ds.append(norm([sg[i]['k1'] * S_K1 for i in range(mm)]
                               + [sg[i]['k2'] * S_K2 for i in range(mm)]))
            dist = sum(ds) / len(ds)
            w2, t2 = rate(ccurve, mm, k1, 'V2', seeds=SEEDS10)
            w0, t0 = rate(ccurve, mm, k1, 'V0', seeds=SEEDS10)
            print(f"{k1:>4} {k1 * k2b / n_:>7.3f} {gh:>12.1f} {dist:>12.1f} "
                  f"{dist / gh:>9.3f} {str(w2) + '/' + str(t2):>11} "
                  f"{str(w0) + '/' + str(t0):>11}")
        print()

    print("The two curves have lam* = 0.340 and 0.070 and, under V0, WALLS THAT")
    print("DIFFER BY ~3x (2026-07-29 T4).  If both V2 columns break at the same")
    print("dist/GH, the curve-to-curve variation was an artefact of the Kannan")
    print("embedding and the CVP formulation depends only on eff = K1*K2/n.")

    print("\nA sharp transition in the last column that lines up with dist/GH")
    print("crossing a fixed constant would be the first PREDICTIVE quantity in")
    print("this thread (nine curve-level invariants failed, 2026-06-21..06-29).")

    # =======================================================================
    # E9 — cross-curve validation of the predictor on FRESH 17-bit curves
    # =======================================================================
    print("\n" + "-" * 78)
    print("E9: cross-curve validation — is dist/GH the same at every wall?")
    print("-" * 78)
    print("E8 fits two curves, which is over-fitting.  Take 6 fresh 17-bit j=0")
    print("GLV curves with lam* spread over (0, 0.5), locate each one's V2 wall")
    print("in eff = K1*K2/n, and ask whether dist/GH (or dist/mu) is constant")
    print("there.  A constant would be a genuine, curve-independent predictor.\n")

    mm9 = 24
    fresh = search_curves(2 ** 16, 2 ** 17, per_bin=1, nbins=6, max_primes=4000)
    EFFS = [0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50]
    print(f"V2 (exact CVP), m={mm9}, 5 seeds.  Wall := largest eff with >= 4/5.")
    print(f"{'lam*':>7} {'n':>7} " + "".join(f"{e:>7.2f}" for e in EFFS))
    walls = []
    for (p9, b9, n9, lam9, G9) in fresh:
        c9 = (p9, b9, n9, lam9, G9)
        k2b9 = math.isqrt(n9) + 1
        cells, wall = [], None
        for eff in EFFS:
            k1 = max(2, round(eff * n9 / k2b9))
            w, t = rate(c9, mm9, k1, 'V2')
            cells.append(f"{str(w) + '/' + str(t):>7}")
            if w >= 4:
                wall = (eff, k1)
        print(f"{lam_star(lam9, n9):>7.3f} {n9:>7} " + "".join(cells))
        walls.append((lam9, n9, k2b9, wall))

    print(f"\n{'lam*':>7} {'wall eff':>9} {'K1':>5} {'mu':>10} {'GH':>10} "
          f"{'dist':>10} {'dist/mu':>8} {'dist/GH':>8}")
    for lam9, n9, k2b9, wall in walls:
        if wall is None:
            print(f"{lam_star(lam9, n9):>7.3f} {'none':>9}")
            continue
        eff, k1 = wall
        S_K1, _, S_K2, _ = scales(n9, k1, k2b9)
        mu, _ = lambda_block_mu(n9, lam9, S_K1, S_K2)
        logdet = (mm9 * math.log(S_K1) + (mm9 - 1) * math.log(n9)
                  + mm9 * math.log(S_K2))
        gh = (math.exp(logdet / (2 * mm9))
              * math.sqrt(2 * mm9 / (2 * math.pi * math.e)))
        dist = n9 * math.sqrt(2 * mm9 / 3)
        print(f"{lam_star(lam9, n9):>7.3f} {eff:>9.2f} {k1:>5} {mu:>10.0f} "
              f"{gh:>10.0f} {dist:>10.0f} {dist / mu:>8.3f} {dist / gh:>8.3f}")

    print("\nNOTE: lambda_1(L2) == mu EXACTLY in every cell measured (E8 side")
    print("check), i.e. the shortest vector of the whole rank-2m lattice is the")
    print("shortest vector of a single 2-D lambda block.")

    print("\n" + "=" * 78)
    print("done")
    print("=" * 78)


if __name__ == '__main__':
    main()
