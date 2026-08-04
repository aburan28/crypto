"""
Thread 23, part 3 — replication + per-instance predictor.

Part 2 (glv_hnp_thread23_cvp_wall.py) found, on the three historical curves:

  * an INFORMATION wall, located by exact CVP: beyond it a decoy vector is
    strictly closer to the target than the planted error, so no decoder can
    win.  12-bit/2557: K1 ~ 8-12.  12-bit/2677: K1 ~ 8.
  * an ALGORITHMIC wall, located by Kannan+LLL/BKZ: 12-bit/2557 K1 ~ 12
    (coincides with its information wall), 12-bit/2677 K1 ~ 4 (well inside).

i.e. the lam*=0.070 curve has a REDUCTION GAP (CVP 4/5 vs LLL 0/5 at K1=6)
while the lam*=0.340 curve has none.  That would mean lam* moves only the
algorithmic wall, not the information wall — refining EXP T4 of 2026-07-29.
Two curves is not evidence, so:

  U8  per-instance check of the nearest-plane criterion
        Babai returns the correct lattice point  iff  |<e,b*_i>|/||b*_i||^2 < 1/2
        for every i.  Confusion matrix, per instance (not averaged as in U7).
  U9  replication on fresh 17-bit j=0 GLV curves with lam* spread over (0,0.5).
        For each curve measure the LLL wall and the CVP wall on a common
        eff = K1*K2/n grid, and correlate gap = eff_CVP - eff_LLL with lam*
        and with nu_hat = lambda_1(L2)/sqrt(det L2)  (the 2026-07-29 separator,
        commit e845207).

Run: python3 glv_hnp_thread23_replication.py
"""

import math
import random

import sympy
from fpylll import IntegerMatrix, LLL, BKZ, CVP

from glv_hnp_thread23_bdd import (
    ec_mul, tonelli_shanks, find_generator, lam_star,
    gen_signatures, scales, build_kannan, recover_kannan,
    build_bdd, error_vector, recover_bdd, babai_nearest_plane,
    gram_schmidt, norm,
)

SEEDS = [42, 1234, 9999, 555, 31337]

# ---------------------------------------------------------------------------
# curve search helpers (same maths as glv_hnp_phase2_lambda_threshold.py)
# ---------------------------------------------------------------------------

def eisenstein_decompose(p):
    for a in range(1, 2 * math.isqrt(p // 3) + 3):
        disc = 4 * p - 3 * a * a
        if disc < 0:
            break
        s = math.isqrt(disc)
        if s * s != disc:
            continue
        for num in (a + s, a - s):
            if num % 2 == 0:
                b = num // 2
                if b >= 0 and a * a - a * b + b * b == p:
                    return (a, b)
    return None

def j0_traces(a, b):
    return [2*a - b, -2*a + b, -(a + b), a + b, 2*b - a, a - 2*b]

def glv_roots(n):
    sq = tonelli_shanks((n - 3) % n, n)
    if sq is None:
        return None
    inv2 = pow(2, -1, n)
    r1 = (n - 1 + sq) * inv2 % n
    r2 = (n - 1 + (n - sq)) * inv2 % n
    if (r1 * r1 + r1 + 1) % n != 0:
        r1, r2 = r2, r1
    if (r1 * r1 + r1 + 1) % n != 0:
        return None
    return (min(r1, r2), max(r1, r2))

def build_curve(p, n, seed=12345):
    rng = random.Random(seed)
    for _ in range(400):
        b = rng.randint(1, p - 1)
        x = rng.randint(0, p - 1)
        rhs = (pow(x, 3, p) + b) % p
        y = tonelli_shanks(rhs, p)
        if y is None or y == 0:
            continue
        if ec_mul((x, y), n, p) is None:
            G = find_generator(p, b, n, seed=seed)
            if G is not None:
                return (p, b, n, None, G)
    return None

def gauss_reduce_2d(u, v):
    def nrm2(w): return w[0] * w[0] + w[1] * w[1]
    def dot(w, z): return w[0] * z[0] + w[1] * z[1]
    if nrm2(u) > nrm2(v):
        u, v = v, u
    while True:
        num, den = dot(v, u), nrm2(u)
        if den == 0:
            break
        q = (2*num + den) // (2*den) if num >= 0 else -((-2*num + den) // (2*den))
        v = (v[0] - q*u[0], v[1] - q*u[1])
        if nrm2(v) >= nrm2(u):
            break
        u, v = v, u
    return u

def nu_hat(n, lam, S_K1, S_K2):
    """lambda_1(L2)/sqrt(det L2) for L2 = <(n*S_K1,0), (-lam*S_K1,S_K2)>."""
    w = gauss_reduce_2d((n * S_K1, 0), (-(lam % n) * S_K1, S_K2))
    l1 = math.sqrt(w[0]*w[0] + w[1]*w[1])
    return l1 / math.sqrt(float(n) * S_K1 * S_K2)

def search_curves(lo, hi, nbins=6, per_bin=1):
    bins = {i: [] for i in range(nbins)}
    p = int(sympy.nextprime(lo))
    while p < hi:
        if p % 3 == 1:
            eis = eisenstein_decompose(p)
            if eis is not None:
                for t in j0_traces(*eis):
                    nc = p + 1 - t
                    if nc < 2 or nc % 3 != 1 or not sympy.isprime(nc):
                        continue
                    roots = glv_roots(nc)
                    if roots is None:
                        continue
                    lam = roots[0]
                    ls = lam_star(lam, nc)
                    idx = min(nbins - 1, int(ls / (0.5 / nbins)))
                    if len(bins[idx]) >= per_bin:
                        continue
                    cur = build_curve(p, nc)
                    if cur is None:
                        continue
                    bins[idx].append((p, cur[1], nc, lam, cur[4]))
        if all(len(v) >= per_bin for v in bins.values()):
            break
        p = int(sympy.nextprime(p))
    out = []
    for i in range(nbins):
        out.extend(bins[i])
    return out


def main():
    # ===========================================================================
    print("=" * 78)
    print("Thread 23 part 3 — per-instance predictor + replication")
    print("=" * 78)

    HIST = [
        ("8-bit/199",     211,  2, 199,  106,  6),
        ("12-bit/2557",   2557, 2, 2659, 1755, 8),
        ("12-bit/2677",   2677, 2, 2647, 185,  10),
    ]
    CURVES = []
    for label, p, b, n, lam, m in HIST:
        G = find_generator(p, b, n)
        CURVES.append((label, (p, b, n, lam, G), m))

    # ---------------------------------------------------------------------------
    print("\n" + "=" * 78)
    print("EXP U8 — per-instance nearest-plane criterion")
    print("=" * 78)
    print("Babai returns the exact planted point iff max_i |<e,b*_i>|/||b*_i||^2 < 1/2.")
    print("Recovery is weaker: only the d-coordinate has to come out right, so")
    print("'criterion holds' must imply 'recovered', but not conversely.\n")

    cells = {(True, True): 0, (True, False): 0, (False, True): 0, (False, False): 0}
    rows_dbg = []
    for label, curve, m in CURVES:
        p, b, n, lam, G = curve
        for k1 in (2, 3, 4, 6, 8, 12, 16, 24):
            k2b = math.isqrt(n) + 1
            for s in SEEDS:
                d = random.Random(s + 7777).randint(1, n - 1)
                sigs = gen_signatures(G, d, m, n, lam, p, k1, k2b, s)
                if len(sigs) < m:
                    continue
                M, t = build_bdd(sigs, n, lam, k1, k2b, 1)
                dim = 2 * m + 1
                A = IntegerMatrix.from_matrix(M)
                LLL.reduction(A)
                rows = [[A[i][j] for j in range(dim)] for i in range(dim)]
                Bs, _ = gram_schmidt(rows)
                e = error_vector(sigs, d, n, k1, k2b, 1, reduced=True)
                e2 = error_vector(sigs, d, n, k1, k2b, 1, reduced=False)
                if norm(e2) < norm(e):
                    e = e2
                mx = 0.0
                for i in range(dim):
                    nb = sum(y * y for y in Bs[i])
                    if nb == 0:
                        continue
                    mx = max(mx, abs(sum(float(e[j]) * Bs[i][j] for j in range(dim))) / nb)
                crit = mx < 0.5
                w = babai_nearest_plane(rows, t)
                S_K1, S_D, S_K2, _ = scales(n, k1, k2b, 1)
                got = recover_bdd(w, t, m, n, S_D, d) is not None
                cells[(crit, got)] += 1
                if crit and not got:
                    rows_dbg.append((label, k1, s, mx))

    tot = sum(cells.values())
    print(f"  instances: {tot}")
    print("                        recovered   not recovered")
    print(f"  criterion  < 1/2       {cells[(True,True)]:5d}        {cells[(True,False)]:5d}")
    print(f"  criterion >= 1/2       {cells[(False,True)]:5d}        {cells[(False,False)]:5d}")
    if rows_dbg:
        print(f"  !! {len(rows_dbg)} violations of the implication (crit => recovered):")
        for r in rows_dbg[:10]:
            print(f"     {r}")
    else:
        print("  implication (criterion => recovered) holds on every instance.")
    slack = cells[(False, True)]
    print(f"  slack (recovered despite criterion failing): {slack}"
          f"  ({100.0*slack/tot:.1f}% of all instances)")

    # ---------------------------------------------------------------------------
    print("\n" + "=" * 78)
    print("EXP U9 — replication on fresh 17-bit curves, lam* spread")
    print("=" * 78)
    print("Searching 17-bit j=0 GLV curves binned by lam* ...")
    FRESH = search_curves(2 ** 16, 2 ** 17, nbins=6, per_bin=1)
    print(f"  found {len(FRESH)} curves")

    M_SIGS = 10
    EFF_GRID = [0.01, 0.02, 0.04, 0.08, 0.12, 0.16, 0.24]


    def k1_for_eff(n, eff):
        k2b = math.isqrt(n) + 1
        return max(2, int(round(eff * n / k2b)))


    def lll_rate(curve, m, k1, seeds):
        p, b, n, lam, G = curve
        k2b = math.isqrt(n) + 1
        wins = 0
        for s in seeds:
            d = random.Random(s + 7777).randint(1, n - 1)
            sigs = gen_signatures(G, d, m, n, lam, p, k1, k2b, s)
            if len(sigs) < m:
                continue
            M = build_kannan(sigs, n, lam, k1, k2b, 1)
            dim = 2 * m + 2
            A = IntegerMatrix.from_matrix(M)
            LLL.reduction(A)
            rows = [[A[i][j] for j in range(dim)] for i in range(dim)]
            S_K1, S_D, S_K2, S_KAN = scales(n, k1, k2b, 1)
            wins += recover_kannan(rows, m, n, S_KAN, S_D, d) is not None
        return wins


    def cvp_rate(curve, m, k1, seeds):
        p, b, n, lam, G = curve
        k2b = math.isqrt(n) + 1
        wins = 0
        for s in seeds:
            d = random.Random(s + 7777).randint(1, n - 1)
            sigs = gen_signatures(G, d, m, n, lam, p, k1, k2b, s)
            if len(sigs) < m:
                continue
            M, t = build_bdd(sigs, n, lam, k1, k2b, 1)
            dim = 2 * m + 1
            A = IntegerMatrix.from_matrix(M)
            LLL.reduction(A)
            w = list(CVP.closest_vector(A, tuple(t)))
            S_K1, S_D, S_K2, _ = scales(n, k1, k2b, 1)
            wins += recover_bdd(w, t, m, n, S_D, d) is not None
        return wins


    def wall(rates, grid, thresh=3):
        """largest grid value still scoring >= thresh out of 5; 0 if none."""
        best = 0.0
        for g, r in zip(grid, rates):
            if r >= thresh:
                best = g
        return best


    print(f"\nm = {M_SIGS} signatures, 5 seeds, eff grid = {EFF_GRID}")
    print("\ncurve                lam*    nu_hat |  arm  | " +
          " ".join(f"{g:<5.2f}" for g in EFF_GRID) + " | wall")
    print("-" * 100)
    SUMMARY = []
    for (p, b, n, lam, G) in FRESH:
        curve = (p, b, n, lam, G)
        k2b = math.isqrt(n) + 1
        S_K1, S_D, S_K2, _ = scales(n, 2, k2b, 1)
        nh = nu_hat(n, lam, n // 8, S_K2)          # K1=8 reference scaling
        ls = lam_star(lam, n)
        lll_r, cvp_r = [], []
        for eff in EFF_GRID:
            k1 = k1_for_eff(n, eff)
            lll_r.append(lll_rate(curve, M_SIGS, k1, SEEDS))
            cvp_r.append(cvp_rate(curve, M_SIGS, k1, SEEDS))
        wl, wc = wall(lll_r, EFF_GRID), wall(cvp_r, EFF_GRID)
        lbl = f"p={p}"
        print(f"{lbl:20s} {ls:.4f}  {nh:.4f} |  LLL | " +
              " ".join(f"{r}/5  " for r in lll_r) + f" | {wl:.2f}")
        print(f"{'':20s} {'':6s}  {'':6s} |  CVP | " +
              " ".join(f"{r}/5  " for r in cvp_r) + f" | {wc:.2f}")
        SUMMARY.append((p, ls, nh, wl, wc, wc - wl))

    print("\n" + "-" * 78)
    print("gap = eff_wall(CVP) - eff_wall(LLL)   [>0 means a reduction gap exists]")
    print("  p          lam*     nu_hat    wall_LLL  wall_CVP   gap")
    for p, ls, nh, wl, wc, gap in sorted(SUMMARY, key=lambda r: r[1]):
        print(f"  {p:<10d} {ls:.4f}   {nh:.4f}    {wl:.2f}      {wc:.2f}     {gap:+.2f}")


    def spearman(xs, ys):
        def rank(v):
            order = sorted(range(len(v)), key=lambda i: v[i])
            r = [0.0] * len(v)
            for pos, i in enumerate(order):
                r[i] = pos
            return r
        rx, ry = rank(xs), rank(ys)
        k = len(xs)
        mx, my = sum(rx) / k, sum(ry) / k
        num = sum((rx[i] - mx) * (ry[i] - my) for i in range(k))
        dx = math.sqrt(sum((rx[i] - mx) ** 2 for i in range(k)))
        dy = math.sqrt(sum((ry[i] - my) ** 2 for i in range(k)))
        return num / (dx * dy) if dx and dy else float('nan')

    if len(SUMMARY) >= 4:
        gaps = [r[5] for r in SUMMARY]
        print(f"\n  spearman(lam*,   gap) = {spearman([r[1] for r in SUMMARY], gaps):+.3f}")
        print(f"  spearman(nu_hat, gap) = {spearman([r[2] for r in SUMMARY], gaps):+.3f}")
        print(f"  spearman(lam*,   wall_CVP) = "
              f"{spearman([r[1] for r in SUMMARY], [r[4] for r in SUMMARY]):+.3f}")
        print(f"  spearman(lam*,   wall_LLL) = "
              f"{spearman([r[1] for r in SUMMARY], [r[3] for r in SUMMARY]):+.3f}")

    print("\n" + "=" * 78)
    print("done")
    print("=" * 78)


if __name__ == "__main__":
    main()
