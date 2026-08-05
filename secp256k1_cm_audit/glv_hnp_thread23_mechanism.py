"""
GLV-HNP Phase 2, Thread 23 (part 2): what actually decides recovery?

Companion to glv_hnp_thread23_gh_ceiling.py, which established (E1-E3):
  * the Gaussian-heuristic ratio R = ||v_planted||/GH(L) tracks eff and nothing
    more, and R < 1 is NOT sufficient -- at eff = 0.157 the 12-bit/2677 curve
    stays at 0-1 of 3 all the way to m = 160 (dim 322) even though R < 1 from
    m ~ 104 on;
  * Thread 23's proposed reformulation (delete the d column, killing the trivial
    vector n*S_D*e_m) gives BIT-IDENTICAL success rates to the original lattice
    on all 12 configurations tested.  The trivial vector is not the obstruction.

This script isolates the mechanism.

  E4  Reduction strength.  LLL vs BKZ-beta for beta in {20, 26, 32, 40} at
      m = 12 (dim 26), where BKZ-26 is exact HKZ reduction of the full lattice.
      If the wall were a reduction-strength wall, beta would move it.

  E5  T5's stated criterion.  The 2026-07-29 entry concluded that recovery is
      "the planted vector is shortest among vectors with last coordinate
      +-S_KANNAN".  That coset is RK + L0 (L0 = the last-coord-0 sublattice), so
      the criterion is exactly a CVP.  Measure gamma = coset_min/||v_planted||;
      the criterion predicts gamma = 1 whenever recovery succeeds.

  E6  Mechanism.  recover_d() scans EVERY row of the reduced basis and accepts
      any row whose m-th coordinate is = d (mod n).  Measure how often the
      planted vector is literally a basis row (hit), and how many rows land in
      the +-S_KANNAN coset (nkan).

  E7  nkan as a separator, on the E1 stratum eff ~ 0.12 where outcomes are
      mixed and R / lam* are both at chance.

  E8  Head-to-head against nu_hat = lambda_1(L2)/sqrt(det L2) (Thread 20c,
      2026-07-29, reported AUC 0.935), the existing a-priori separator.
      24 fresh 17-bit curves x 2 eff levels x 6 seeds.

RESULTS (2026-08-05 run; full output in glv_hnp_thread23_mechanism_output.txt)

  E4  The wall is NOT a reduction-strength wall.  LLL, BKZ-20, BKZ-26, BKZ-32
      and BKZ-40 give byte-identical success counts on all 8 configurations,
      both curves.  At m = 12 the lattice is 26-dimensional, so BKZ-26 is exact
      HKZ reduction of the whole lattice: no amount of reduction moves it.
      It also exposes the sharpest form of the C1/C2 split -- at essentially
      the same eff (0.1572 vs 0.1564) and the same R (1.428 vs 1.425),
      12-bit/2677 is 0/5 and 12-bit/2557 is 5/5, at every beta.

  E5  T5's CRITERION IS FALSIFIED.  gamma = coset_min/||v_planted|| < 1 in
      56 of the 60 instances measured, INCLUDING every 5/5 success:
        2557 K1=4  gamma 0.962  5/5      2677 K1=4  gamma 0.978  5/5
        2557 K1=8  gamma 0.893  5/5      2677 K1=8  gamma 0.885  0/5
      A shorter coset vector almost always exists, and recovery happens anyway.
      Note also 0.893 (5/5) vs 0.885 (0/5): gamma does not separate either.

  E6  MECHANISM.  The planted vector is essentially never a row of the reduced
      basis (hit = 1/5 in one config, 0/5 in the other seven) even where
      recovery is 5/5.  recover_d() succeeds because it accepts ANY row of the
      +-S_KANNAN coset whose m-th coordinate is = d (mod n), and the trivial
      vector n*S_D*e_m generates the kernel of that map -- so a whole
      1-parameter family of rows decodes correctly.  Recovery is therefore a
      coupon-collector event over the coset rows, not an SVP event and not a
      coset-minimum event.  The count of such rows, nkan, orders all 8
      configurations monotonically:
        nkan  7.2  7.2  7.0  5.4  5.4  4.4  3.8  2.8
        win   5/5  5/5  5/5  5/5  4/5  2/5  0/5  0/5

  E7  On the mixed stratum eff ~ 0.12, where R (AUC 0.407) and lam* (0.481) are
      both at chance, nkan separates with AUC 0.9630
      (success [5.0, 7.2], failure [2.6, 5.6]).

  E8  48 configs, 24 fresh 17-bit curves x 2 eff levels x 6 seeds:
        AUC(nkan)   = 0.9843     AUC(nu_hat) = 0.9617
        AUC(R)      = 0.5679     AUC(lam*)   = 0.5314
      nu_hat REPLICATES on a fresh sample (Thread 20c reported 0.935).
      Spearman(nu_hat, nkan) = -0.6205, the sign predicted by Thread 20c's
      reading: a skew L2 (low nu_hat) leaves the planted vector uncrowded, so
      more rows survive in the Kannan coset (high nkan), so recovery fires.
      nkan is a-posteriori -- it is free once you have run the reduction, so it
      is a diagnostic, not a screening test.  nu_hat remains the a-priori
      predictor of record.

Run: python3 glv_hnp_thread23_mechanism.py
"""

import math
import random
import sys

import sympy
from fpylll import IntegerMatrix, LLL, CVP

import importlib.util
_here = __file__.rsplit("/", 1)[0]
_spec = importlib.util.spec_from_file_location(
    "_t23a", _here + "/glv_hnp_thread23_gh_ceiling.py")
# the ceiling script runs its experiments at import time, so pull the
# definitions out of its source instead of executing the whole module.
_src = open(_here + "/glv_hnp_thread23_gh_ceiling.py").read()
_g = {"__name__": "_t23a_defs"}
exec(_src.split("def banner(t):")[0], _g)

find_generator = _g["find_generator"]
gen_signatures = _g["gen_signatures"]
scales = _g["scales"]
planted_vector = _g["planted_vector"]
norm = _g["norm"]
build_glv_lattice = _g["build_glv_lattice"]
run_once = _g["run_once"]
sweep = _g["sweep"]
R_ratio = _g["R_ratio"]
auc = _g["auc"]
glv_eigenvalue = _g["glv_eigenvalue"]
eisenstein_decompose = _g["eisenstein_decompose"]
j0_traces = _g["j0_traces"]
identify_b_for_n = _g["identify_b_for_n"]

HIST = [("2677", 2677, 2, 2647, 185), ("2557", 2557, 2, 2659, 1755)]


def banner(t):
    print("\n" + "=" * 78)
    print(t)
    print("=" * 78)
    sys.stdout.flush()


def find_curves(bits, want):
    out = []
    p = sympy.nextprime(2 ** (bits - 1))
    while p < 2 ** bits and len(out) < want:
        if p % 3 == 1:
            dec = eisenstein_decompose(p)
            if dec:
                a, bb = dec
                for t in j0_traces(a, bb):
                    n = p + 1 - t
                    if n > 4 and sympy.isprime(n) and n % 3 == 1:
                        lam = glv_eigenvalue(n)
                        if lam is None:
                            continue
                        bc = identify_b_for_n(p, n)
                        if bc is None:
                            continue
                        G = find_generator(p, bc, n)
                        if G is None:
                            continue
                        out.append((p, bc, n, lam, G))
                        break
        p = sympy.nextprime(p + 1)
    return out


# ---------------------------------------------------------------------------
# E5 helper: exact minimum over the coset {last coord = S_KANNAN}
# ---------------------------------------------------------------------------

def coset_min(sigs, n, lam, k1, k2):
    """min_{u in L0} ||RK + u||, computed as a CVP in the (2m+1)-dim L0'."""
    m = len(sigs)
    S_K1, S_D, S_K2, S_KAN = scales(n, k1, k2)
    d0 = 2 * m + 1
    rows = []
    for i in range(m):
        r = [0] * d0; r[i] = n * S_K1; rows.append(r)
    r = [0] * d0
    for i in range(m):
        r[i] = sigs[i]['B'] * S_K1
    r[m] = S_D
    rows.append(r)
    for i in range(m):
        r = [0] * d0
        r[i] = -(lam % n) * S_K1
        r[m + 1 + i] = S_K2
        rows.append(r)
    B = IntegerMatrix.from_matrix(rows)
    LLL.reduction(B)
    t = [-sigs[i]['A'] * S_K1 for i in range(m)] + [0] * (m + 1)
    c = CVP.closest_vector(B, t)
    diff = [c[i] - t[i] for i in range(d0)]
    return math.sqrt(sum(x * x for x in diff) + S_KAN * S_KAN)


# ---------------------------------------------------------------------------
# E8 helper: nu_hat (Thread 20c, glv_hnp_phase2_nuhat_control.py:11)
# ---------------------------------------------------------------------------

def gauss_reduce_2d(u, v):
    def nrm2(w):
        return w[0] * w[0] + w[1] * w[1]
    if nrm2(u) > nrm2(v):
        u, v = v, u
    while True:
        num = u[0] * v[0] + u[1] * v[1]
        den = nrm2(u)
        if den == 0:
            break
        q = (2 * num + den) // (2 * den) if num >= 0 else -((-2 * num + den) // (2 * den))
        v = (v[0] - q * u[0], v[1] - q * u[1])
        if nrm2(v) >= nrm2(u):
            break
        u, v = v, u
    return u


def nu_hat(n, lam, k1, k2):
    S_K1, S_D, S_K2, S_KAN = scales(n, k1, k2)
    w = gauss_reduce_2d((n * S_K1, 0), (-(lam % n) * S_K1, S_K2))
    return math.sqrt(w[0] ** 2 + w[1] ** 2) / math.sqrt(n * S_K1 * S_K2)


# ---------------------------------------------------------------------------
# shared inner loop: LLL once, report nkan / hit / recovery
# ---------------------------------------------------------------------------

def probe(G, p, n, lam, k1, k2, m, seed):
    d = random.Random(seed + 7777).randint(1, n - 1)
    sigs = gen_signatures(G, d, m, n, lam, p, k1, k2, seed)
    if len(sigs) < m:
        return None
    dim = 2 * m + 2
    S_K1, S_D, S_K2, S_KAN = scales(n, k1, k2)
    v = planted_vector(sigs, d, n, k1, k2)
    A = IntegerMatrix.from_matrix(build_glv_lattice(sigs, n, lam, k1, k2))
    LLL.reduction(A)
    rows = [[A[i][j] for j in range(dim)] for i in range(dim)]
    hit = any(r == v or [-x for x in r] == v for r in rows)
    nkan = sum(1 for r in rows if abs(r[dim - 1]) == S_KAN)
    ok = False
    for r in rows:
        if abs(r[dim - 1]) != S_KAN:
            continue
        sg = 1 if r[dim - 1] > 0 else -1
        if (sg * r[m]) % n == d:
            ok = True
            break
    norms = sorted(norm(r) for r in rows)
    pv = norm(v)
    return {'ok': ok, 'hit': hit, 'nkan': nkan, 'pv': pv, 'sigs': sigs, 'd': d,
            'rank': sum(1 for x in norms if x < pv) + 1, 'pv_over_r1': pv / norms[0]}


print("=" * 78)
print("Thread 23 (part 2) -- mechanism of the GLV-HNP Phase-2 recovery event")
print("=" * 78)

SEEDS5 = [42, 1234, 9999, 555, 31337]
SEEDS6 = [42, 1234, 9999, 555, 31337, 271828]

hist = []
for label, p, b, n, lam in HIST:
    G = find_generator(p, b, n)
    assert G is not None and (lam * lam + lam + 1) % n == 0, label
    hist.append((label, (p, b, n, lam, G)))

# ===========================================================================
banner("E4: is the K1 wall a reduction-STRENGTH wall?  (dim 26 => BKZ-26 = HKZ)")
for label, cv in hist:
    p, b, n, lam, G = cv
    k2 = math.isqrt(n) + 1
    print(f"\n=== 12-bit/{label}  n={n} K2={k2} lam*={min(lam,n-lam)/n:.4f} ===")
    print(f"{'K1':>4}{'m':>4}{'dim':>5}{'eff':>8}{'R':>8}  "
          + "".join(f"{('b=' + str(x) if x else 'LLL'):>8}" for x in [0, 20, 26, 32, 40]))
    for k1, m in [(6, 12), (8, 12), (8, 16), (12, 12)]:
        line = f"{k1:>4}{m:>4}{2*m+2:>5}{k1*k2/n:>8.4f}{R_ratio(n,m,k1,k2):>8.4f}  "
        for beta in [0, 20, 26, 32, 40]:
            w, t, _, _ = sweep(cv, m, k1, SEEDS5, use_bkz=bool(beta),
                               bkz_beta=max(beta, 2))
            line += f"{f'{w}/{t}':>8}"
            sys.stdout.flush()
        print(line)

# ===========================================================================
banner("E5: T5's criterion -- is v_planted the COSET minimum?  (gamma = 1 iff yes)")
print(f"{'curve':>6}{'lam*':>8}{'K1':>4}{'m':>4}{'eff':>8}{'R':>8}"
      f"{'gamma':>8}{'gam<1':>7}{'win':>7}")
for label, cv in hist:
    p, b, n, lam, G = cv
    k2 = math.isqrt(n) + 1
    for k1, m in [(4, 12), (6, 12), (8, 12), (8, 16), (12, 12), (16, 12)]:
        gams, wins, strict = [], 0, 0
        for seed in SEEDS5:
            d = random.Random(seed + 7777).randint(1, n - 1)
            sigs = gen_signatures(G, d, m, n, lam, p, k1, k2, seed)
            if len(sigs) < m:
                continue
            pv = norm(planted_vector(sigs, d, n, k1, k2))
            cm = coset_min(sigs, n, lam, k1, k2)
            gams.append(cm / pv)
            strict += (cm < pv * (1 - 1e-12))
            wins += run_once(cv, m, d, k1, seed)['ok']
        print(f"{label:>6}{min(lam,n-lam)/n:>8.4f}{k1:>4}{m:>4}{k1*k2/n:>8.4f}"
              f"{R_ratio(n,m,k1,k2):>8.4f}{sum(gams)/len(gams):>8.4f}"
              f"{f'{strict}/{len(gams)}':>7}{f'{wins}/{len(SEEDS5)}':>7}")
        sys.stdout.flush()

# ===========================================================================
banner("E6: mechanism -- is v_planted even a basis row?  how many rows in the coset?")
print(f"{'curve':>6}{'K1':>4}{'m':>4}{'eff':>8}{'hit':>6}{'nkan':>6}"
      f"{'rank':>7}{'pv/r1':>8}{'win':>7}")
for label, cv in hist:
    p, b, n, lam, G = cv
    k2 = math.isqrt(n) + 1
    for k1, m in [(4, 12), (6, 12), (8, 12), (12, 12)]:
        rs = [probe(G, p, n, lam, k1, k2, m, s) for s in SEEDS5]
        rs = [r for r in rs if r]
        ns = len(rs)
        hits = sum(r['hit'] for r in rs)
        wins = sum(r['ok'] for r in rs)
        print(f"{label:>6}{k1:>4}{m:>4}{k1*k2/n:>8.4f}"
              f"{f'{hits}/{ns}':>6}"
              f"{sum(r['nkan'] for r in rs)/ns:>6.1f}"
              f"{sum(r['rank'] for r in rs)/ns:>7.1f}"
              f"{sum(r['pv_over_r1'] for r in rs)/ns:>8.3f}"
              f"{f'{wins}/{ns}':>7}")
        sys.stdout.flush()

# ===========================================================================
banner("E7: nkan on the mixed stratum eff ~ 0.12 (12 curves), where R and lam* fail")
C12 = find_curves(17, 12)
recs7 = []
m = 12
print(f"{'p':>7}{'n':>7}{'lam*':>8}{'eff':>8}{'R':>7}{'nkan':>7}{'win':>7}")
for (p, b, n, lam, G) in C12:
    k2 = math.isqrt(n) + 1
    k1 = max(2, int(round(0.12 * n / k2)))
    rs = [probe(G, p, n, lam, k1, k2, m, s) for s in SEEDS5]
    rs = [r for r in rs if r]
    mn = sum(r['nkan'] for r in rs) / len(rs)
    wins = sum(r['ok'] for r in rs)
    recs7.append({'nkan': mn, 'win': wins, 'tot': len(rs),
                  'lam': min(lam, n - lam) / n, 'R': R_ratio(n, m, k1, k2)})
    print(f"{p:>7}{n:>7}{min(lam,n-lam)/n:>8.4f}{k1*k2/n:>8.4f}"
          f"{R_ratio(n,m,k1,k2):>7.3f}{mn:>7.1f}{f'{wins}/{len(rs)}':>7}")
    sys.stdout.flush()
lab7 = [r['win'] * 2 >= r['tot'] for r in recs7]
print(f"\nstratum eff~0.12: {sum(lab7)}/{len(lab7)} majority-successes")
print(f"  AUC(nkan, high=success) = {auc([-r['nkan'] for r in recs7], lab7):.4f}")
print(f"  AUC(R,    low=success)  = {auc([r['R'] for r in recs7], lab7):.4f}")
print(f"  AUC(lam*, high=success) = {auc([-r['lam'] for r in recs7], lab7):.4f}")

# ===========================================================================
banner("E8: nkan (a-posteriori) vs nu_hat (a-priori, Thread 20c) -- 24 curves x 2 eff")
C24 = find_curves(17, 24)
recs8 = []
print(f"{'p':>7}{'n':>7}{'eff':>7}{'lam*':>7}{'nu_hat':>8}{'nkan':>6}{'win':>7}")
for (p, b, n, lam, G) in C24:
    k2 = math.isqrt(n) + 1
    for teff in (0.12, 0.15):
        k1 = max(2, int(round(teff * n / k2)))
        rs = [probe(G, p, n, lam, k1, k2, m, s) for s in SEEDS6]
        rs = [r for r in rs if r]
        if not rs:
            continue
        rec = {'nkan': sum(r['nkan'] for r in rs) / len(rs),
               'win': sum(r['ok'] for r in rs), 'tot': len(rs),
               'nu': nu_hat(n, lam, k1, k2), 'lam': min(lam, n - lam) / n,
               'eff': k1 * k2 / n, 'R': R_ratio(n, m, k1, k2)}
        recs8.append(rec)
        wl = f"{rec['win']}/{len(rs)}"
        print(f"{p:>7}{n:>7}{rec['eff']:>7.3f}{rec['lam']:>7.3f}"
              f"{rec['nu']:>8.4f}{rec['nkan']:>6.1f}{wl:>7}")
        sys.stdout.flush()

lab8 = [r['win'] * 2 >= r['tot'] for r in recs8]
print(f"\nN={len(recs8)} configs, {sum(lab8)} majority-successes")
for name, sc in [('nu_hat (low=success) ', [r['nu'] for r in recs8]),
                 ('nkan   (high=success)', [-r['nkan'] for r in recs8]),
                 ('R      (low=success) ', [r['R'] for r in recs8]),
                 ('lam*   (high=success)', [-r['lam'] for r in recs8])]:
    print(f"  AUC({name}) = {auc(sc, lab8):.4f}")

# rank correlation between the a-priori and a-posteriori observables
def spearman(x, y):
    def rk(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        for pos, i in enumerate(order):
            r[i] = pos + 1.0
        return r
    rx, ry = rk(x), rk(y)
    mx, my = sum(rx) / len(rx), sum(ry) / len(ry)
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    den = math.sqrt(sum((a - mx) ** 2 for a in rx) * sum((b - my) ** 2 for b in ry))
    return num / den if den else float('nan')

print(f"\n  Spearman(nu_hat, nkan) = "
      f"{spearman([r['nu'] for r in recs8], [r['nkan'] for r in recs8]):.4f}")

banner("DONE")
