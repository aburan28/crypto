"""
Thread 23c — the K1 wall is set by a curve-level invariant: cmin.

Chain of reasoning from Thread 23 / 23b (same day):

 * 23  Projecting the trivial vector n*e_m out of the Phase-2 lattice does not
       move the K1 wall.  The projection works (the trivial vector is gone) but
       buys nothing.
 * 23b Reason: the attack never needed the planted vector.  A reduced row with
       last coordinate S_KANNAN and d-coefficient delta = d (mod n) yields d
       directly, and such rows exist for EVERY gamma in the coset
            w[i]/S_K1 = k1_i + lam*(k2_i - gamma_i)  (mod n),
       i.e. every alternative GLV split of the SAME nonce
            k_full,i = k1_i + lam*k2_i = k1'_i + lam*gamma_i.
       Measured: on 12-bit/2557 at K1=8 the recovered rows have
       gamma - k2 = +-50 with the 2-D lam-lattice minimum at (-3,50):
       shifting k2 by 50 costs only 3 in k1, and 3 < K1 = 8, so the alternative
       split is admissible.  On 12-bit/2677 the minimum is (-14,43): a shift
       costs 14 > 8, no alternative exists, and the attack dies at K1 ~ 4-6.

Hypothesis (H-cmin).  Define the curve-level invariant

       cmin(n, lam, K2) = min over 0 < |e| <= K2 of |lam*e mod n|_centered

which is the cheapest k1-cost of an alternative GLV split.  Then the attack
survives past the point where the planted vector itself is findable exactly
when K1 > cmin, so the K1 wall is an increasing function of 1/cmin.

This CONTRADICTS the 2026-07-29 conclusion "no curve-level invariant can
separate C1 from C2".  That conclusion was drawn after six invariants
(delta/n, kappa(M), q_cf, max_q_cf, max_a, a_corn/n) failed, and after T2
tested the 2-D minimum mu only as the RATIO rho = mu/||v_planted||, which
divides out the very quantity that matters.  cmin is mu's k1-coordinate, not
its length.

Success criterion here is the REAL attack, not an oracle: every candidate d is
checked by testing d*G == Q against the public key (~dim scalar mults).

Run: python3 glv_hnp_phase2_cmin.py
"""

import math
import random

from fpylll import IntegerMatrix, LLL

_src = open(__file__.replace('cmin', 'projected')).read()
_head = _src.split('# ' + '=' * 75)[0]
_G = {}
exec(compile(_head, 'glv_hnp_phase2_projected.py', 'exec'), _G)

find_generator = _G['find_generator']
gen_signatures = _G['gen_signatures']
scales = _G['scales']
build_base_lattice = _G['build_base_lattice']
build_proj_lattice = _G['build_proj_lattice']
lam_star = _G['lam_star']
search_curves = _G['search_curves']
ec_mul = _G['ec_mul']
d_from_k = _G['d_from_k']

# ---------------------------------------------------------------------------
# The curve-level invariant
# ---------------------------------------------------------------------------

def cmin(n, lam, k2_bound):
    """min over 0<|e|<=K2 of the centered residue |lam*e mod n|."""
    best = n
    best_e = 0
    for e in range(1, k2_bound + 1):
        r = lam * e % n
        r = min(r, n - r)
        if r < best:
            best, best_e = r, e
    return best, best_e

# ---------------------------------------------------------------------------
# Real attack criterion: candidate d verified against the public key
# ---------------------------------------------------------------------------

def attack_base(red, sigs, n, lam, k1_bound, k2_bound, G, Q, p):
    """HIST read-out: d off column m of any Kannan row, verified by d*G == Q."""
    m = len(sigs)
    dim = 2 * m + 2
    _S_K1, _S_D, _S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    for row in red:
        if abs(row[dim - 1]) != S_KAN:
            continue
        sign = 1 if row[dim - 1] > 0 else -1
        d_cand = (sign * row[m]) % n
        if d_cand == 0:
            continue
        if ec_mul(G, d_cand, p) == Q:
            return d_cand
    return None

def attack_proj(red, sigs, n, lam, k1_bound, k2_bound, G, Q, p):
    """PROJ has no d-column: d must come from the recovered (k1,gamma) split.
    Accept any split, not just the planted one, then verify d*G == Q."""
    m = len(sigs)
    S_K1, _S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    for row in red:
        if abs(row[2 * m]) != S_KAN:
            continue
        r = row if row[2 * m] > 0 else [-x for x in row]
        if any(r[i] % S_K1 for i in range(m)) or any(r[m + i] % S_K2 for i in range(m)):
            continue
        k1s = [r[i] // S_K1 for i in range(m)]
        k2s = [r[m + i] // S_K2 for i in range(m)]
        for i in range(m):
            B = sigs[i]['B'] % n
            if B == 0 or math.gcd(B, n) != 1:
                continue
            k_full = (k1s[i] + lam * k2s[i]) % n
            d_cand = (k_full - sigs[i]['A']) * pow(B, -1, n) % n
            if d_cand and ec_mul(G, d_cand, p) == Q:
                return d_cand
    return None

def trial(curve, m, k1_bound, seed, methods=('base', 'proj')):
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    d = random.Random(seed + 7777).randint(1, n - 1)
    Q = ec_mul(G, d, p)
    sigs = gen_signatures(G, d, m, n, lam, p, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return None
    out = {}
    if 'base' in methods:
        A = IntegerMatrix.from_matrix(build_base_lattice(sigs, n, lam, k1_bound, k2_bound))
        LLL.reduction(A)
        red = [[A[i][j] for j in range(2 * m + 2)] for i in range(A.nrows)]
        out['base'] = attack_base(red, sigs, n, lam, k1_bound, k2_bound, G, Q, p) == d
    if 'proj' in methods:
        A = IntegerMatrix.from_matrix(build_proj_lattice(sigs, n, lam, k1_bound, k2_bound))
        LLL.reduction(A)
        red = [[A[i][j] for j in range(2 * m + 1)] for i in range(A.nrows)]
        out['proj'] = attack_proj(red, sigs, n, lam, k1_bound, k2_bound, G, Q, p) == d
    return out

def rate(curve, m, k1_bound, seeds, methods=('base', 'proj')):
    w = {k: 0 for k in methods}
    t = 0
    for s in seeds:
        r = trial(curve, m, k1_bound, s, methods)
        if r is None:
            continue
        t += 1
        for k in methods:
            w[k] += bool(r[k])
    return w, t


print("=" * 78)
print("Thread 23c — cmin as the curve-level predictor of the K1 wall")
print("=" * 78)
print("Success = real attack: candidate d verified by d*G == Q (no oracle).\n")

SEEDS = [42, 1234, 9999, 555, 31337]

# ---------------------------------------------------------------------------
# C1 — the historical pair, with cmin annotated
# ---------------------------------------------------------------------------
print("-" * 78)
print("C1: historical pair, K1 sweep under the real criterion (m=12)")
print("-" * 78)
HIST = [("12-bit/2557", 2557, 2, 2659, 1755), ("12-bit/2677", 2677, 2, 2647, 185)]
K1_GRID = [2, 3, 4, 6, 8, 12, 16, 24]
for label, p, b, n, lam in HIST:
    G = find_generator(p, b, n)
    curve = (p, b, n, lam, G)
    k2 = math.isqrt(n) + 1
    c, e = cmin(n, lam, k2)
    print(f"\n{label}: n={n}, lam*={lam_star(lam,n):.4f}, K2={k2}, "
          f"cmin={c} (at e={e})")
    print("  " + " ".join(f"K1={k:<4}" for k in K1_GRID))
    rows = {'base': [], 'proj': []}
    for k1 in K1_GRID:
        w, t = rate(curve, 12, k1, SEEDS)
        for k in rows:
            rows[k].append(f"{w[k]}/{t}")
    for k in ('base', 'proj'):
        print(f"  {k.upper():<5}" + " ".join(f"{c2:<6}" for c2 in rows[k]))
    print(f"  H-cmin predicts the wall at K1 ~ {c}"
          f" (alternative splits admissible once K1 > {c})")

# ---------------------------------------------------------------------------
# C2 — 17-bit sweep: does cmin rank-order the K1 wall?
# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("C2: 17-bit curves — cmin vs measured K1 wall (m=12, 5 seeds, BASE)")
print("-" * 78)
curves = search_curves(2 ** 16, 2 ** 17, per_bin=2, nbins=10)
print(f"collected {len(curves)} curves\n")

K1_SWEEP = [2, 4, 8, 16, 32, 64, 128]
recs = []
for c in curves:
    p, b, n, lam, G = c
    k2 = math.isqrt(n) + 1
    cm, e_star = cmin(n, lam, k2)
    wall = 0
    succ = []
    for k1 in K1_SWEEP:
        w, t = rate(c, 12, k1, SEEDS, methods=('base',))
        succ.append(w['base'])
        if w['base'] >= 4:
            wall = k1
    recs.append({'n': n, 'lam': lam, 'ls': lam_star(lam, n), 'cmin': cm,
                 'e': e_star, 'wall': wall, 'succ': succ, 'k2': k2})

print(f"{'n':>7} {'lam*':>7} {'cmin':>6} {'e*':>5} {'wall':>5}   "
      + " ".join(f"K1={k:<4}" for k in K1_SWEEP))
for r in sorted(recs, key=lambda z: z['cmin']):
    print(f"{r['n']:>7} {r['ls']:>7.4f} {r['cmin']:>6} {r['e']:>5} {r['wall']:>5}   "
          + " ".join(f"{s}/5   " for s in r['succ']))

# rank correlation between cmin and wall
def spearman(xs, ys):
    def ranks(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        rk = [0.0] * len(v)
        i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
                j += 1
            avg = (i + j) / 2.0 + 1
            for k in range(i, j + 1):
                rk[order[k]] = avg
            i = j + 1
        return rk
    rx, ry = ranks(xs), ranks(ys)
    nn = len(xs)
    mx, my = sum(rx) / nn, sum(ry) / nn
    num = sum((rx[i] - mx) * (ry[i] - my) for i in range(nn))
    dx = math.sqrt(sum((rx[i] - mx) ** 2 for i in range(nn)))
    dy = math.sqrt(sum((ry[i] - my) ** 2 for i in range(nn)))
    return num / (dx * dy) if dx and dy else float('nan')

cs = [r['cmin'] for r in recs]
ws = [r['wall'] for r in recs]
ls = [r['ls'] for r in recs]
print(f"\nSpearman(cmin, wall)  = {spearman(cs, ws):+.4f}   (H-cmin predicts negative)")
print(f"Spearman(lam*, wall)  = {spearman(ls, ws):+.4f}   (2026-07-29: lam* is not a predictor)")

print("\n" + "=" * 78)
print("done")
print("=" * 78)
