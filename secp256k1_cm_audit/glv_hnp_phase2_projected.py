"""
GLV-HNP Phase 2 — Thread 23: does projecting out the trivial vector help?

Background (RESEARCH_AUTOLAB_LOG.md, 2026-07-29, T5)
----------------------------------------------------
The Phase-2 Kannan lattice K (dim 2m+2) contains the parasitic vector
n*S_D*e_m, of norm n, which is shorter than the planted vector
(||v_planted|| ~ n*sqrt(2m/3 + 4/3)) for every m >= 1.  So the planted vector
is never lambda_1 and recovery is a BDD/coset condition, not SVP.

The 2026-07-29 next-step proposal was: quotient out that direction and see
whether the attack gets stronger.  Since the parasitic vector is axis-aligned,
the projection is literally "delete column m"; `build_projected_lattice` in
glv_hnp_phase2_core.py gives an explicit basis of pi(L_K) (dim 2m+1,
det = det(L_K)/(n*S_D)).  d is recovered from the k's instead of read off a
coordinate, so no information is lost.

Falsifier as stated on 2026-07-29
---------------------------------
  "if sv/pv rises above 1 after the reformulation AND the K1 wall in T4 moves
   outward on the lam*=0.07 curve (currently K1 ~ 4-6), the reformulation is a
   real improvement; if the wall stays at K1 ~ 4-6, then the wall is
   information-theoretic and Phase 2 is at its ceiling."

Experiments
-----------
E0  correctness of lattice P: determinants, and planted vector membership.
E1  anchors: sv/pv and the energy profile of lambda_1, both lattices.
E2  the T4 K1-wall grid, both lattices (the falsifier's second clause).
E3  m-sweep at the wall.
E4  fresh 17-bit curves at three bias strengths; and a test of which
    lattice-geometric quantity predicts success in P.

Also reports the ORACLE-FREE success criterion alongside Thread 20a's oracle
criterion (which counted a trial as a win iff some row's d-entry equalled the
secret -- not a criterion an attacker can evaluate, and at 8-bit n it fires
by chance).

Run: python3 glv_hnp_phase2_projected.py
"""

import math
import random
import sys
import time

import sympy

import importlib.util
_spec = importlib.util.spec_from_file_location(
    "_core", __file__.rsplit("/", 1)[0] + "/glv_hnp_phase2_core.py")
core = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(core)

from sympy import Matrix

SEEDS = [42, 1234, 9999, 555, 31337]

print("=" * 78)
print("Thread 23 — projecting out the parasitic vector of the Phase-2 lattice")
print("=" * 78)

# Anchor curves as recorded in the 2026-06-15 / 2026-07-29 log entries.
def anchor(p, b, n, lam, k1, label):
    G = core.find_generator(p, b, n)
    assert G is not None, label
    assert core.ec_mul(G, n, p) is None
    assert (lam * lam + lam + 1) % n == 0, f"{label}: lam is not a GLV eigenvalue"
    return {'label': label, 'p': p, 'b': b, 'n': n, 'lam': lam, 'G': G,
            'k1': k1, 'k2': math.isqrt(n) + 1}

ANCHORS = [
    anchor(211,  2,  199, 106, 2, "8-bit/199"),
    anchor(2557, 2, 2659, 1755, 8, "12-bit/2557"),
    anchor(2677, 2, 2647,  185, 8, "12-bit/2677"),
]


# ===========================================================================
# E0 — is lattice P correct?
# ===========================================================================
print("\n" + "=" * 78)
print("E0 — correctness of the projected lattice P")
print("=" * 78)
print("Checks: det(K) and det(P) against the closed forms, det(K)/det(P) =")
print("n*S_D, and integrality of the planted vector's coordinates in each.")
print()

def coords_in_lattice(M, v):
    """Integer coordinates of v in the row lattice of M, or None."""
    A = Matrix(M).T
    try:
        x = A.solve(Matrix(v))
    except Exception:
        return None
    if all(c.q == 1 for c in x):
        return [int(c) for c in x]
    return None

print(f"{'curve':<14} {'m':>2} {'det(K) ok':>10} {'det(P) ok':>10} "
      f"{'det ratio ok':>13} {'plant in K':>11} {'plant in P':>11} {'triv in K':>10}")
e0_all_ok = True
for a in ANCHORS:
    for m in (3, 4):
        n, lam, k1, k2 = a['n'], a['lam'], a['k1'], a['k2']
        d = random.Random(11).randint(1, n - 1)
        sigs = core.gen_signatures(a['G'], d, m, n, lam, a['p'], k1, k2, seed=7)
        if len(sigs) < m:
            continue
        S_K1, S_D, S_K2, S_KAN = core.scales(n, k1, k2)
        K = core.build_glv_lattice(sigs, n, lam, k1, k2)
        P = core.build_projected_lattice(sigs, n, lam, k1, k2)
        detK = abs(Matrix(K).det())
        detP = abs(Matrix(P).det())
        expK = (n * S_K1) ** m * S_D * S_K2 ** m * S_KAN
        expP = n ** (m - 1) * S_K1 ** m * S_K2 ** m * S_KAN
        okK = (detK == expK)
        okP = (detP == expP)
        okR = (detK == detP * n * S_D)
        pK = coords_in_lattice(K, core.planted_vector(sigs, d, n, k1, k2))
        pP = coords_in_lattice(P, core.planted_vector_proj(sigs, n, k1, k2))
        tK = coords_in_lattice(K, core.trivial_vector(n, m, k1, k2))
        row_ok = okK and okP and okR and pK is not None and pP is not None \
            and tK is not None
        e0_all_ok &= row_ok
        print(f"{a['label']:<14} {m:>2} {str(okK):>10} {str(okP):>10} "
              f"{str(okR):>13} {str(pK is not None):>11} "
              f"{str(pP is not None):>11} {str(tK is not None):>10}")
print(f"\nE0 verdict: {'ALL CHECKS PASS' if e0_all_ok else 'FAILURE'}")
print("So P is exactly pi(K), the planted vector survives the projection, and")
print("the parasitic n*S_D*e_m really is in K.")


# ===========================================================================
# E1 — where does lambda_1 sit, in each lattice?
# ===========================================================================
print("\n" + "=" * 78)
print("E1 — energy profile of the shortest reduced vector")
print("=" * 78)
print("blocks: k1 = cols 0..m-1, d = col m (K only), k2 = k2-block,")
print("kan = Kannan column.  Fractions of ||sv||^2.")
print()

def energy_profile(row, m, lattice):
    tot = sum(x * x for x in row) or 1
    if lattice == 'K':
        k1 = sum(row[i] ** 2 for i in range(m))
        dd = row[m] ** 2
        k2 = sum(row[m + 1 + i] ** 2 for i in range(m))
        kan = row[2 * m + 1] ** 2
    else:
        k1 = sum(row[i] ** 2 for i in range(m))
        dd = 0
        k2 = sum(row[m + i] ** 2 for i in range(m))
        kan = row[2 * m] ** 2
    return k1 / tot, dd / tot, k2 / tot, kan / tot

M_E1 = 8
print(f"{'curve':<14} {'lat':>4} {'K1':>3} {'sv/pv':>7} {'k1':>6} {'d':>6} "
      f"{'k2':>6} {'kan':>6} {'kan=S?':>7}")
for a in ANCHORS:
    n, lam, k1, k2 = a['n'], a['lam'], a['k1'], a['k2']
    d = random.Random(42 + 7777).randint(1, n - 1)
    sigs = core.gen_signatures(a['G'], d, M_E1, n, lam, a['p'], k1, k2, seed=42)
    S_K1, S_D, S_K2, S_KAN = core.scales(n, k1, k2)
    for lat in ('K', 'P'):
        if lat == 'K':
            dim = 2 * M_E1 + 2
            M = core.build_glv_lattice(sigs, n, lam, k1, k2)
            pv = core.planted_vector(sigs, d, n, k1, k2)
        else:
            dim = 2 * M_E1 + 1
            M = core.build_projected_lattice(sigs, n, lam, k1, k2)
            pv = core.planted_vector_proj(sigs, n, k1, k2)
        red = core._reduce(M, dim)
        nz = [r for r in red if any(r)]
        sv = min(nz, key=core.norm)
        f1, fd, f2, fk = energy_profile(sv, M_E1, lat)
        print(f"{a['label']:<14} {lat:>4} {k1:>3} "
              f"{core.norm(sv)/core.norm(pv):>7.3f} {f1:>6.3f} {fd:>6.3f} "
              f"{f2:>6.3f} {fk:>6.3f} {str(abs(sv[dim-1])==S_KAN):>7}")


# ===========================================================================
# E2 — the T4 K1-wall grid, both lattices
# ===========================================================================
print("\n" + "=" * 78)
print("E2 — K1 wall (replicates T4 of 2026-07-29), both lattices, 5 seeds")
print("=" * 78)
print("Each cell: verified-wins/exact-wins out of 5, at m = 12.")
print("  verified = d recovered and confirmed by d*G == Q  (attacker-executable)")
print("  exact    = the in-box planted decomposition was the vector found")
print()

K1_GRID = [2, 3, 4, 6, 8, 12, 16, 24]
M_E2 = 12
print(f"{'curve':<14} {'lat':>4} " + " ".join(f"{k:>7}" for k in K1_GRID))
e2 = {}
for a in ANCHORS[1:]:                     # the two 12-bit curves, as in T4
    n, k2, lam = a['n'], a['k2'], a['lam']
    for lat in ('K', 'P'):
        cells, row = [], []
        for k1 in K1_GRID:
            wo, wv, we, tot, ratio, kf = core.success(
                a['p'], n, lam, a['G'], M_E2, k1, k2, SEEDS, lattice=lat)
            cells.append(f"{wv}/{we}")
            row.append((k1, wo, wv, we, tot, ratio))
        e2[(a['label'], lat)] = row
        print(f"{a['label']:<14} {lat:>4} " + " ".join(f"{c:>7}" for c in cells))

def wall(row, idx):
    """Largest K1 whose win count at tuple position idx equals the total."""
    best = None
    for k1, wo, wv, we, tot, ratio in row:
        if row_win(k1, wo, wv, we, idx) == tot:
            best = k1
    return best

def row_win(k1, wo, wv, we, idx):
    return {2: wv, 3: we}[idx]

print("\nWall (largest K1 with 5/5), by criterion:")
for a in ANCHORS[1:]:
    vK = wall(e2[(a['label'], 'K')], 2)
    vP = wall(e2[(a['label'], 'P')], 2)
    eK = wall(e2[(a['label'], 'K')], 3)
    eP = wall(e2[(a['label'], 'P')], 3)
    print(f"  {a['label']:<14} verified  K: {str(vK):>4}  P: {str(vP):>4}   "
          f"exact  K: {str(eK):>4}  P: {str(eP):>4}")


# ===========================================================================
# E3 — does more data move the wall in P?
# ===========================================================================
print("\n" + "=" * 78)
print("E3 — m-sweep at the wall (12-bit/2677, K1 = 8, the T4b failure point)")
print("=" * 78)

a = ANCHORS[2]
print(f"{'m':>4} {'K ver':>7} {'K exact':>8} {'P ver':>7} {'P exact':>8} "
      f"{'K sv/pv':>9} {'P sv/pv':>9}")
for m in (8, 12, 16, 24, 32):
    _, wvK, weK, totK, rK, _ = core.success(a['p'], a['n'], a['lam'], a['G'], m,
                                            8, a['k2'], SEEDS, lattice='K')
    _, wvP, weP, totP, rP, _ = core.success(a['p'], a['n'], a['lam'], a['G'], m,
                                            8, a['k2'], SEEDS, lattice='P')
    print(f"{m:>4} {str(wvK)+'/'+str(totK):>7} {str(weK)+'/'+str(totK):>8} "
          f"{str(wvP)+'/'+str(totP):>7} {str(weP)+'/'+str(totP):>8} "
          f"{rK:>9.3f} {rP:>9.3f}")


# ===========================================================================
# E4 — fresh curves at three bias strengths, and what predicts success in P
# ===========================================================================
print("\n" + "=" * 78)
print("E4 — fresh 17-bit curves at eff = 0.05 / 0.15 / 0.25 (replicates T3)")
print("=" * 78)

def harvest(lo, hi, want):
    out, seen = [], set()
    p = int(sympy.nextprime(lo))
    while p < hi and len(out) < want:
        if p % 3 == 1:
            eis = core.eisenstein_decompose(p)
            if eis is not None:
                for t in core.j0_traces(*eis):
                    n = p + 1 - t
                    if n < 2 or n in seen or n % 3 != 1 or not sympy.isprime(n):
                        continue
                    roots = core.glv_roots(n)
                    if roots is None:
                        continue
                    b = core.identify_twist(p, n)
                    if b is None:
                        continue
                    G = core.find_generator(p, b, n)
                    if G is None:
                        continue
                    seen.add(n)
                    out.append({'p': p, 'b': b, 'n': n, 'lam': roots[0], 'G': G})
                    break
        p = int(sympy.nextprime(p))
    return out

t0 = time.time()
CURVES = harvest(2 ** 16, 2 ** 17, 20)
print(f"harvested {len(CURVES)} 17-bit j=0 GLV curves in {time.time()-t0:.1f}s\n")

M_E4 = 12
e4_rows = []
for eff_t in (0.05, 0.15, 0.25):
    nK = nP = 0
    print(f"--- eff target {eff_t} ---")
    print(f"{'n':>7} {'lam*':>6} {'K1':>4} {'eff':>6} {'nu_hat':>7} "
          f"{'pv/l1L2':>8} {'K ver':>7} {'P ver':>7} {'K exa':>6} {'P exa':>6}")
    for c in CURVES:
        n, lam = c['n'], c['lam']
        k2 = math.isqrt(n) + 1
        k1 = max(2, int(eff_t * n / k2))
        eff = k1 * k2 / n
        S_K1, S_D, S_K2, S_KAN = core.scales(n, k1, k2)
        l1L2, rootdet, nu = core.rival_sublattice_nu(n, lam, k1, k2)
        pv_exp = math.sqrt(M_E4 * (k1 * S_K1) ** 2 / 3.0
                           + M_E4 * (k2 * S_K2) ** 2 / 3.0 + S_KAN ** 2)
        _, wvK, weK, tot, rK, _ = core.success(c['p'], n, lam, c['G'], M_E4, k1,
                                               k2, SEEDS, lattice='K')
        _, wvP, weP, totP, rP, _ = core.success(c['p'], n, lam, c['G'], M_E4, k1,
                                                k2, SEEDS, lattice='P')
        nK += (wvK == tot); nP += (wvP == totP)
        e4_rows.append({'eff': eff_t, 'n': n, 'lam_star': core.lam_star(lam, n),
                        'nu': nu, 'ratio': pv_exp / l1L2,
                        'okK': wvK == tot, 'okP': wvP == totP,
                        'wfK': wvK, 'wfP': wvP, 'weK': weK, 'weP': weP,
                        'tot': tot})
        print(f"{n:>7} {core.lam_star(lam,n):>6.3f} {k1:>4} {eff:>6.3f} "
              f"{nu:>7.3f} {pv_exp/l1L2:>8.3f} "
              f"{str(wvK)+'/'+str(tot):>7} {str(wvP)+'/'+str(totP):>7} "
              f"{str(weK)+'/'+str(tot):>6} {str(weP)+'/'+str(totP):>6}")
    print(f"  => 5/5 verified: K {nK}/{len(CURVES)}   P {nP}/{len(CURVES)}\n")

# --- predictor check on the pooled E4 sample --------------------------------
print("-" * 78)
print("Predictors of oracle-free success in the PROJECTED lattice P")
print("-" * 78)

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
    rx, ry = rank(xs), rank(ys)
    mx, my = sum(rx) / len(rx), sum(ry) / len(ry)
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    dx = math.sqrt(sum((a - mx) ** 2 for a in rx))
    dy = math.sqrt(sum((b - my) ** 2 for b in ry))
    return num / (dx * dy) if dx and dy else float('nan')

pooled = [r for r in e4_rows if 0 < r['wfP'] < r['tot'] or True]
phatP = [r['wfP'] / r['tot'] for r in pooled]
phatK = [r['wfK'] / r['tot'] for r in pooled]
for key, label in (('nu', 'nu_hat (Thread 20c)'),
                   ('ratio', 'pv / lambda_1(L2)'),
                   ('lam_star', 'lam* (falsified 20a)')):
    xs = [r[key] for r in pooled]
    print(f"  spearman({label:<22}, p_hat in P) = {spearman(xs, phatP):+.3f}"
          f"   [in K: {spearman(xs, phatK):+.3f}]")

# within-eff correlations, since eff dominates
for eff_t in (0.05, 0.15, 0.25):
    sub = [r for r in e4_rows if r['eff'] == eff_t]
    ph = [r['wfP'] / r['tot'] for r in sub]
    if len(set(ph)) < 2:
        print(f"  eff={eff_t}: p_hat constant in P ({ph[0]:.2f}) — no correlation")
        continue
    print(f"  eff={eff_t}: spearman(nu_hat, p_hat_P) = "
          f"{spearman([r['nu'] for r in sub], ph):+.3f}, "
          f"spearman(pv/l1L2, p_hat_P) = "
          f"{spearman([r['ratio'] for r in sub], ph):+.3f}")

print("\nDone.")
