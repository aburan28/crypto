"""
GLV-HNP Phase 2 — Thread 24: the recovered vector is planted mod the rival
sublattice, and its reduced norm is the predictor.

Motivation
----------
Thread 23 (glv_hnp_phase2_projected.py) found that projecting out the trivial
vector n*S_D*e_m changes nothing: lattices K and P agree cell-for-cell.  It also
found a large gap between two success criteria:

    verified  d recovered and confirmed by d*G == Q          (12-bit/2557, K1=6: 5/5)
    exact     the in-box planted decomposition was found     (same cell:      0/5)

So d is being recovered from a lattice vector that is NOT the planted vector.
This script identifies that vector.

Hypothesis H-coset
------------------
Let L2_i = < (n*S_K1, 0), (-lam*S_K1, S_K2) > be the rival 2D sublattice on the
coordinate pair (k1-col i, k2-col i) -- the same L2 whose normalised first
minimum is Thread 20c's nu_hat.  Every element of L2_i has zero d-coordinate
and zero Kannan coordinate, so adding one to the planted vector preserves both
d and the Kannan marker.  The trivial vector T = n*S_D*e_m does too: it shifts
the d-coordinate by n, and d is only ever read mod n.  Recovery therefore only
needs SOME vector of the coset

    v_planted + Z*T + (L2_0 + ... + L2_{m-1})

and LLL will find the shortest such, whose length is not ||v_planted|| but

    ||v_planted mod L2^m||  =: pv_red.

Predictions
  (P1) whenever `verified` succeeds, the winning row minus the planted vector
       lies in L2_0 + ... + L2_{m-1}  (exactly, as an integer membership test);
  (P2) pv_red predicts success better than pv;
  (P3) pv_red explains the sign of nu_hat: skew L2 (small nu_hat) => the coset
       representative is shorter => attack easier, which is why 20c found low
       nu_hat easier despite the rival vector itself being short.

Falsifier: if (P1) fails on any winning trial, H-coset is wrong outright.  If
(P1) holds but pv_red is no better than pv as a predictor, (P2)/(P3) fail and
the mechanism is incomplete.

Run: python3 glv_hnp_phase2_coset.py
"""

import math
import random
import time

import sympy

import importlib.util
_spec = importlib.util.spec_from_file_location(
    "_core", __file__.rsplit("/", 1)[0] + "/glv_hnp_phase2_core.py")
core = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(core)

SEEDS = [42, 1234, 9999, 555, 31337]

print("=" * 78)
print("Thread 24 — the coset mechanism behind Phase-2 recovery")
print("=" * 78)


# ---------------------------------------------------------------------------
# The rival sublattice L2: membership and reduction
# ---------------------------------------------------------------------------

def rival_basis(n, lam, S_K1, S_K2):
    return (n * S_K1, 0), (-(lam % n) * S_K1, S_K2)

def in_rival(w, n, lam, S_K1, S_K2):
    """Is w = (w0, w1) in L2 = < (n*S_K1,0), (-lam*S_K1,S_K2) > ?

    w1 = t*S_K2 forces t, then w0 + t*lam*S_K1 must be a multiple of n*S_K1."""
    w0, w1 = w
    if w1 % S_K2:
        return False
    t = w1 // S_K2
    return (w0 + t * (lam % n) * S_K1) % (n * S_K1) == 0

def gauss_reduce_pair(u, v):
    """Full Lagrange-Gauss reduction; returns the reduced BASIS (u, v)."""
    nrm2 = lambda w: w[0] * w[0] + w[1] * w[1]
    dot = lambda w, z: w[0] * z[0] + w[1] * z[1]
    if nrm2(u) > nrm2(v):
        u, v = v, u
    while True:
        num, den = dot(v, u), nrm2(u)
        if den == 0:
            break
        q = (2 * num + den) // (2 * den) if num >= 0 else -((-2 * num + den) // (2 * den))
        v = (v[0] - q * u[0], v[1] - q * u[1])
        if nrm2(v) >= nrm2(u):
            break
        u, v = v, u
    return u, v

def reduce_mod_rival(w, n, lam, S_K1, S_K2):
    """Coset representative of w modulo L2, via Babai rounding on the
    Gauss-reduced basis (exact enough in 2D: repeat until stable)."""
    b1, b2 = gauss_reduce_pair(*rival_basis(n, lam, S_K1, S_K2))
    det = b1[0] * b2[1] - b1[1] * b2[0]
    if det == 0:
        return w
    cur = w
    for _ in range(8):
        x = (cur[0] * b2[1] - cur[1] * b2[0]) / det
        y = (b1[0] * cur[1] - b1[1] * cur[0]) / det
        cx, cy = round(x), round(y)
        if cx == 0 and cy == 0:
            break
        cur = (cur[0] - cx * b1[0] - cy * b2[0],
               cur[1] - cx * b1[1] - cy * b2[1])
    # local search over the 8 neighbours, cheap and removes rounding slop
    best = cur
    for a in (-1, 0, 1):
        for b in (-1, 0, 1):
            cand = (cur[0] - a * b1[0] - b * b2[0], cur[1] - a * b1[1] - b * b2[1])
            if cand[0] ** 2 + cand[1] ** 2 < best[0] ** 2 + best[1] ** 2:
                best = cand
    return best

def planted_reduced_norm(sigs, d_secret, n, lam, k1_bound, k2_bound):
    """||v_planted mod L2^m|| in lattice K (d and Kannan coords untouched)."""
    m = len(sigs)
    S_K1, S_D, S_K2, S_KAN = core.scales(n, k1_bound, k2_bound)
    tot = (d_secret * S_D) ** 2 + S_KAN ** 2
    for i in range(m):
        w = reduce_mod_rival((sigs[i]['k1'] * S_K1, sigs[i]['k2'] * S_K2),
                             n, lam, S_K1, S_K2)
        tot += w[0] ** 2 + w[1] ** 2
    return math.sqrt(tot)


# ===========================================================================
# E5 — (P1): is the winning row  planted + rival ?
# ===========================================================================
print("\n" + "=" * 78)
print("E5 — decomposition of the winning row  (tests P1)")
print("=" * 78)
print("For every trial where `verified` succeeds, take the winning row w,")
print("sign-normalise it, and test w - v_planted:")
print("  d,kan ok : the Kannan coordinate of the difference is 0 and its")
print("             d-coordinate is a multiple of n*S_D (i.e. a multiple of")
print("             the trivial vector T, which does not change d mod n)")
print("  in L2^m  : each (k1-col i, k2-col i) block lies in L2_i (exact test)")
print()

def anchor(p, b, n, lam, k1, label):
    G = core.find_generator(p, b, n)
    assert G is not None and core.ec_mul(G, n, p) is None
    assert (lam * lam + lam + 1) % n == 0
    return {'label': label, 'p': p, 'b': b, 'n': n, 'lam': lam, 'G': G,
            'k1': k1, 'k2': math.isqrt(n) + 1}

ANCHORS = [
    anchor(211,  2,  199, 106, 2, "8-bit/199"),
    anchor(2557, 2, 2659, 1755, 8, "12-bit/2557"),
    anchor(2677, 2, 2647,  185, 8, "12-bit/2647"),
]

M_E5 = 12
print(f"{'curve':<14} {'K1':>3} {'trials':>7} {'verified':>9} {'exact':>7} "
      f"{'d,kan ok':>10} {'in L2^m':>9} {'pv_red/pv':>10}")

e5_fail = 0
for a in ANCHORS[1:]:
    n, lam, k2 = a['n'], a['lam'], a['k2']
    for k1 in (2, 4, 6, 8):
        S_K1, S_D, S_K2, S_KAN = core.scales(n, k1, k2)
        dim = 2 * M_E5 + 2
        nver = nexact = nzero = ncoset = 0
        ratios = []
        for s in SEEDS:
            d = random.Random(s + 7777).randint(1, n - 1)
            sigs = core.gen_signatures(a['G'], d, M_E5, n, lam, a['p'], k1, k2, s)
            if len(sigs) < M_E5:
                continue
            Q = core.ec_mul(a['G'], d, a['p'])
            M = core.build_glv_lattice(sigs, n, lam, k1, k2)
            red = core._reduce(M, dim)
            pv = core.planted_vector(sigs, d, n, k1, k2)
            ratios.append(planted_reduced_norm(sigs, d, n, lam, k1, k2)
                          / core.norm(pv))
            got = core.recover_verified_K(red, sigs, n, lam, k1, k2, a['G'], Q,
                                          a['p'])
            if got is None:
                continue
            nver += 1
            if core.recover_d_free(red, sigs, n, lam, k1, k2) is not None:
                nexact += 1
            # decompose every winning row; the trial counts if any of them
            # is  planted + (multiple of T) + (element of L2^m)
            hit_zero = hit_coset = False
            for row in red:
                if abs(row[dim - 1]) != S_KAN:
                    continue
                sgn = 1 if row[dim - 1] > 0 else -1
                u = [sgn * x for x in row]
                if u[M_E5] % n != d:          # not a winning row
                    continue
                w = [u[j] - pv[j] for j in range(dim)]
                if w[dim - 1] != 0 or w[M_E5] % (n * S_D) != 0:
                    continue
                hit_zero = True
                if all(in_rival((w[i], w[M_E5 + 1 + i]), n, lam, S_K1, S_K2)
                       for i in range(M_E5)):
                    hit_coset = True
                    break
            nzero += hit_zero
            ncoset += hit_coset
        if nver and ncoset != nver:
            e5_fail += 1
        print(f"{a['label']:<14} {k1:>3} {len(SEEDS):>7} {nver:>9} {nexact:>7} "
              f"{nzero:>10} {ncoset:>9} "
              f"{sum(ratios)/len(ratios) if ratios else float('nan'):>10.3f}")

print(f"\nE5 verdict: {'P1 HOLDS on every winning trial' if e5_fail == 0 else 'P1 FAILS'}")


# ===========================================================================
# E6 — (P2)/(P3): does pv_red predict better than pv?
# ===========================================================================
print("\n" + "=" * 78)
print("E6 — pv_red vs pv as a predictor, on fresh 17-bit curves")
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
print(f"harvested {len(CURVES)} curves in {time.time()-t0:.1f}s\n")

M_E6 = 12
rows = []
for eff_t in (0.05, 0.10, 0.15, 0.20, 0.25):
    for c in CURVES:
        n, lam = c['n'], c['lam']
        k2 = math.isqrt(n) + 1
        k1 = max(2, int(eff_t * n / k2))
        S_K1, S_D, S_K2, S_KAN = core.scales(n, k1, k2)
        dim = 2 * M_E6 + 2
        det = (n * S_K1) ** M_E6 * S_D * S_K2 ** M_E6 * S_KAN
        gh = math.sqrt(dim / (2 * math.pi * math.e)) * det ** (1.0 / dim)
        wins, pvs, pvrs = 0, [], []
        for s in SEEDS:
            d = random.Random(s + 7777).randint(1, n - 1)
            sigs = core.gen_signatures(c['G'], d, M_E6, n, lam, c['p'], k1, k2, s)
            if len(sigs) < M_E6:
                continue
            Q = core.ec_mul(c['G'], d, c['p'])
            red = core._reduce(core.build_glv_lattice(sigs, n, lam, k1, k2), dim)
            wins += core.recover_verified_K(red, sigs, n, lam, k1, k2, c['G'],
                                            Q, c['p']) is not None
            pvs.append(core.norm(core.planted_vector(sigs, d, n, k1, k2)))
            pvrs.append(planted_reduced_norm(sigs, d, n, lam, k1, k2))
        avg = lambda xs: sum(xs) / len(xs)
        _, _, nu = core.rival_sublattice_nu(n, lam, k1, k2)
        rows.append({'eff': eff_t, 'n': n, 'nu': nu,
                     'p_hat': wins / len(SEEDS),
                     'pv_gh': avg(pvs) / gh, 'pvred_gh': avg(pvrs) / gh,
                     'shrink': avg(pvrs) / avg(pvs)})

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

def auc(xs, labels):
    """P(x_pos < x_neg) — predictors here are 'smaller is easier'."""
    pos = [x for x, l in zip(xs, labels) if l]
    neg = [x for x, l in zip(xs, labels) if not l]
    if not pos or not neg:
        return float('nan')
    c = sum((a < b) + 0.5 * (a == b) for a in pos for b in neg)
    return c / (len(pos) * len(neg))

print(f"{'eff':>5} {'n_curves':>9} {'mean p_hat':>11} {'mean pv/GH':>11} "
      f"{'mean pvred/GH':>14} {'mean shrink':>12}")
for eff_t in (0.05, 0.10, 0.15, 0.20, 0.25):
    sub = [r for r in rows if r['eff'] == eff_t]
    a = lambda k: sum(r[k] for r in sub) / len(sub)
    print(f"{eff_t:>5} {len(sub):>9} {a('p_hat'):>11.3f} {a('pv_gh'):>11.3f} "
          f"{a('pvred_gh'):>14.3f} {a('shrink'):>12.3f}")

lab = [r['p_hat'] == 1.0 for r in rows]
print(f"\npooled over all {len(rows)} (curve, eff) cells, "
      f"{sum(lab)} of them 5/5:")
for key, name in (('pv_gh', 'pv/GH'), ('pvred_gh', 'pv_red/GH'),
                  ('nu', 'nu_hat')):
    xs = [r[key] for r in rows]
    print(f"  {name:<12} spearman(., p_hat) = {spearman(xs, [r['p_hat'] for r in rows]):+.3f}"
          f"   AUC(5/5 vs not) = {auc(xs, lab):.3f}")

print("\nwithin-eff (removes the dominant eff effect):")
for eff_t in (0.05, 0.10, 0.15, 0.20, 0.25):
    sub = [r for r in rows if r['eff'] == eff_t]
    ph = [r['p_hat'] for r in sub]
    if len(set(ph)) < 2:
        print(f"  eff={eff_t}: p_hat constant ({ph[0]:.2f})")
        continue
    print(f"  eff={eff_t}: pv/GH {spearman([r['pv_gh'] for r in sub], ph):+.3f}   "
          f"pv_red/GH {spearman([r['pvred_gh'] for r in sub], ph):+.3f}   "
          f"nu_hat {spearman([r['nu'] for r in sub], ph):+.3f}")

print("\nDone.")
