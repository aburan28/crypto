"""
Thread 23d — scoping H-cmin: the alternative-split mechanism only operates
when cmin < K1, which is a rare fluke of the 2-D lam-lattice.

Thread 23c measured, over 20 17-bit curves,
    Spearman(cmin, K1-wall) = +0.4708
where H-cmin predicted a NEGATIVE correlation.  But on every one of those 20
curves cmin >= 49 while the measured wall was <= 32, i.e. cmin > K1 throughout:
no alternative GLV split is admissible anywhere in the swept range, so the
mechanism is switched OFF on the whole sample and the +0.47 cannot be it.

The mechanism was observed directly (Thread 23b) only on 12-bit/2557, which has
cmin = 3 against K2 = 52.  The 2-D lattice <(n,0),(-lam,1)> has covolume n, so
its minimum is generically ~sqrt(n) ~ 51 there: cmin = 3 is a ~17x downward
fluke.  12-bit/2677 has cmin = 14, no fluke, and dies at K1 ~ 4-6.

This script tests the scoped claim on matched curves:

  H1 (distribution)  r = cmin/sqrt(n) has P(r < c) ~ c^2 for small c
                     (area heuristic for a random 2-D lattice of covolume n).

  H2 (matched test)  Among 12-bit curves of the same size, those with small r
                     have a systematically LARGER K1 wall than those with
                     typical r.  Falsifier: if the low-r and typical-r groups
                     have the same wall distribution, the 23b mechanism does
                     not generalise beyond the single curve it was seen on.

Run: python3 glv_hnp_phase2_cmin_matched.py
"""

import math
import random

import sympy
from fpylll import IntegerMatrix, LLL

_src = open(__file__.replace('cmin_matched', 'projected')).read()
_head = _src.split('# ' + '=' * 75)[0]
_G = {}
exec(compile(_head, 'glv_hnp_phase2_projected.py', 'exec'), _G)

find_generator = _G['find_generator']
gen_signatures = _G['gen_signatures']
scales = _G['scales']
build_base_lattice = _G['build_base_lattice']
lam_star = _G['lam_star']
eisenstein_decompose = _G['eisenstein_decompose']
j0_traces = _G['j0_traces']
glv_roots = _G['glv_roots']
build_curve = _G['build_curve']
ec_mul = _G['ec_mul']

_c = open(__file__.replace('cmin_matched', 'cmin')).read()
exec(compile(_c.split("print(\"=\" * 78)")[0].split('_src = ')[0], 'x', 'exec'), {})


def cmin(n, lam, k2_bound):
    best, best_e = n, 0
    for e in range(1, k2_bound + 1):
        r = lam * e % n
        r = min(r, n - r)
        if r < best:
            best, best_e = r, e
    return best, best_e


def attack_base(red, n, m, k1_bound, k2_bound, G, Q, p):
    dim = 2 * m + 2
    _S1, _SD, _S2, S_KAN = scales(n, k1_bound, k2_bound)
    for row in red:
        if abs(row[dim - 1]) != S_KAN:
            continue
        sign = 1 if row[dim - 1] > 0 else -1
        d_cand = (sign * row[m]) % n
        if d_cand and ec_mul(G, d_cand, p) == Q:
            return d_cand
    return None


def rate(curve, m, k1_bound, seeds):
    p, b, n, lam, G = curve
    k2 = math.isqrt(n) + 1
    w = 0
    for seed in seeds:
        d = random.Random(seed + 7777).randint(1, n - 1)
        Q = ec_mul(G, d, p)
        sigs = gen_signatures(G, d, m, n, lam, p, k1_bound, k2, seed)
        if len(sigs) < m:
            continue
        A = IntegerMatrix.from_matrix(build_base_lattice(sigs, n, lam, k1_bound, k2))
        LLL.reduction(A)
        red = [[A[i][j] for j in range(2 * m + 2)] for i in range(A.nrows)]
        w += (attack_base(red, n, m, k1_bound, k2, G, Q, p) == d)
    return w


def enumerate_12bit(lo, hi):
    """All j=0 GLV curves with p in [lo,hi): n prime, n = 1 mod 3."""
    out = []
    p = int(sympy.nextprime(lo))
    while p < hi:
        if p % 3 == 1:
            eis = eisenstein_decompose(p)
            if eis is not None:
                a_e, b_e = eis
                for t in j0_traces(a_e, b_e):
                    n_c = p + 1 - t
                    if n_c < 2 or n_c % 3 != 1 or not sympy.isprime(n_c):
                        continue
                    roots = glv_roots(n_c)
                    if roots is None:
                        continue
                    lam = roots[0]
                    k2 = math.isqrt(n_c) + 1
                    cm, e = cmin(n_c, lam, k2)
                    out.append({'p': p, 'n': n_c, 'lam': lam, 'k2': k2,
                                'cmin': cm, 'e': e,
                                'r': cm / math.sqrt(n_c),
                                'ls': lam_star(lam, n_c)})
        p = int(sympy.nextprime(p))
    return out


print("=" * 78)
print("Thread 23d — scoping H-cmin (matched low-r vs typical-r 12-bit curves)")
print("=" * 78)

SEEDS = [42, 1234, 9999, 555, 31337]
K1_GRID = [2, 3, 4, 6, 8, 12, 16, 24, 32]

# ---------------------------------------------------------------------------
# H1 — distribution of r = cmin/sqrt(n)
# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("H1: distribution of r = cmin/sqrt(n) over all 12-bit j=0 GLV curves")
print("-" * 78)
pool = enumerate_12bit(2 ** 11, 2 ** 12)
print(f"enumerated {len(pool)} curves with p in [2048, 4096)\n")
print(f"{'c':>6} {'#(r<c)':>8} {'emp P':>8} {'c^2':>8}")
for c in (0.1, 0.15, 0.2, 0.3, 0.5, 0.75, 1.0):
    k = sum(1 for x in pool if x['r'] < c)
    print(f"{c:>6.2f} {k:>8} {k/len(pool):>8.4f} {min(1.0,c*c):>8.4f}")

rs = sorted(x['r'] for x in pool)
print(f"\nr: min={rs[0]:.4f}  median={rs[len(rs)//2]:.4f}  max={rs[-1]:.4f}")

# ---------------------------------------------------------------------------
# H2 — matched groups
# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("H2: K1 wall, low-r group vs typical-r group (m=12, 5 seeds, real criterion)")
print("-" * 78)

low = sorted([x for x in pool if x['r'] < 0.20], key=lambda z: z['r'])[:8]
typ = sorted([x for x in pool if 0.45 <= x['r'] <= 0.85],
             key=lambda z: abs(z['n'] - 2650))[:8]

print(f"low-r group: {len(low)} curves,  typical-r group: {len(typ)} curves")
print(f"\n{'group':>7} {'n':>6} {'lam*':>7} {'cmin':>5} {'e*':>5} {'r':>6} {'wall':>5}   "
      + " ".join(f"K1={k:<3}" for k in K1_GRID))

results = {'low': [], 'typ': []}
for gname, grp in (('low', low), ('typ', typ)):
    for x in grp:
        cur = build_curve(x['p'], x['n'])
        if cur is None:
            continue
        curve = (x['p'], cur[1], x['n'], x['lam'], cur[4])
        succ, wall = [], 0
        for k1 in K1_GRID:
            w = rate(curve, 12, k1, SEEDS)
            succ.append(w)
            if w >= 4:
                wall = k1
        results[gname].append({'x': x, 'wall': wall, 'succ': succ})
        print(f"{gname:>7} {x['n']:>6} {x['ls']:>7.4f} {x['cmin']:>5} {x['e']:>5} "
              f"{x['r']:>6.3f} {wall:>5}   " + " ".join(f"{s}/5 " for s in succ))

for g in ('low', 'typ'):
    walls = [z['wall'] for z in results[g]]
    if walls:
        walls_s = sorted(walls)
        print(f"\n{g}-r walls: {walls}  median={walls_s[len(walls_s)//2]}  "
              f"mean={sum(walls)/len(walls):.1f}")

lw = [z['wall'] for z in results['low']]
tw = [z['wall'] for z in results['typ']]
if lw and tw:
    # Mann-Whitney U (small samples, exact enough for a log entry)
    u = sum(1 for a in lw for b in tw if a > b) + 0.5 * sum(1 for a in lw for b in tw if a == b)
    auc = u / (len(lw) * len(tw))
    print(f"\nAUC(low-r wall > typical-r wall) = {auc:.4f}   "
          f"(0.5 = no effect; H2 predicts >> 0.5)")

print("\n" + "=" * 78)
print("done")
print("=" * 78)
