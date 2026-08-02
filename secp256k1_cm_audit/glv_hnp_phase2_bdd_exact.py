"""
GLV-HNP Phase 2 — Thread 23b: is the K1 wall information-theoretic or reduction-quality?

Why this script exists
----------------------
The 2026-07-29 log entry proposed Thread 23 with a two-branch falsifier:

    "if sv/pv rises above 1 after the reformulation and the K1 wall in T4 moves
     outward on the lam*=0.07 curve (currently K1~4-6), the reformulation is a real
     improvement; if the wall stays at K1~4-6, then the wall is information-theoretic
     and Phase 2 is at its ceiling."

`glv_hnp_phase2_projected.py` (Thread 23a, today) settled the first branch: the wall
does NOT move.  Variant B (d-column quotiented out) has exactly the same wall as the
2026-07-29 baseline -- K1=8 on lam*=0.34, K1=4 on lam*=0.07 -- even though sv/pv rose
from 0.37 to as high as 1.000, i.e. the trivial vector really was removed.

But the second branch of that disjunction does not follow, and this script shows it is
false.  "The wall stays" and "the wall is information-theoretic" are different claims,
and the entropy count separates them:

    m signatures give m*log2(n) bits of constraint.
    The unknowns are m*log2(K1) + m*log2(K2) + log2(n) bits.
    The solution is unique (whp) once  m*log2(K1*K2/n) < -log2(n),  i.e.

        eff := K1*K2/n  <  n^(-1/m).

    At m=12, n~2650 that is eff < 0.518 -- while the OBSERVED wall sits at
    eff = 8*52/2659 = 0.156 (lam*=0.34) and eff = 4*52/2647 = 0.079 (lam*=0.07).

So the observed wall is a factor 3.3x-6.6x in eff BELOW the unicity bound.  There is a
wide band in which d is information-theoretically determined but the lattice attack
does not find it.  That is a reduction-quality / BDD-radius wall, not an information
wall -- unless the planted vector genuinely stops being the closest lattice point well
before unicity fails, which is exactly what the experiments below test directly.

Experiments
-----------
  Y1  unicity bound vs the observed wall, tabulated in eff.
  Y2  EXACT CVP (fpylll enumeration, not Babai) on the rank-2m homogeneous lattice L0.
      If the exact closest vector to -t IS the planted vector at a K1 where LLL/Babai
      fail, the wall is reduction-quality.  If the exact closest vector is some OTHER
      lattice point, the BDD instance itself is ambiguous at that radius.
  Y3  BKZ ladder (beta = 2/10/20/30/40) on variant B at the failure point, to see how
      far stronger reduction pushes the wall.

Run: python3 glv_hnp_phase2_bdd_exact.py
"""

import math
import os
import random

from fpylll import CVP, IntegerMatrix, LLL, BKZ

_HERE = os.path.dirname(os.path.abspath(__file__))
_SRC = os.path.join(_HERE, "glv_hnp_phase2_projected.py")

# Same pattern as Thread 23a: exec only the helper prefix (everything above that
# module's "EXPERIMENT X1" banner) so its experiment suite does not re-run here.
with open(_SRC) as _f:
    _text = _f.read()
_P = {'__file__': _SRC, '__name__': '_t23a'}
exec(compile(_text[:_text.index("# EXPERIMENT X1")].rsplit("# ====", 1)[0],
             _SRC, "exec"), _P)

scales = _P['scales']
gen_signatures = _P['gen_signatures']
find_generator = _P['find_generator']
lam_star = _P['lam_star']
norm = _P['norm']
d_from_k = _P['d_from_k']
build_homogeneous_lattice = _P['build_homogeneous_lattice']
build_projected_lattice = _P['build_projected_lattice']
recover_d_projected = _P['recover_d_projected']
projected_planted = _P['projected_planted']

SEEDS = [42, 1234, 9999, 555, 31337]
HIST = [
    ("12-bit/2557", 2557, 2, 2659, 1755),
    ("12-bit/2677", 2677, 2, 2647, 185),
]
curves = []
for label, p, b, n, lam in HIST:
    G = find_generator(p, b, n)
    assert G is not None and (lam * lam + lam + 1) % n == 0, label
    curves.append((label, (p, b, n, lam, G)))

print("=" * 78)
print("Thread 23b — is the Phase-2 K1 wall information-theoretic, or reduction-quality?")
print("=" * 78)

# --- Y1: unicity bound -------------------------------------------------------
print("\n" + "-" * 78)
print("Y1: unicity bound (eff < n^(-1/m)) vs the observed wall")
print("-" * 78)
print(f"{'curve':<14} {'lam*':>7} {'n':>6} {'K2':>4} {'m':>3} "
      f"{'wall K1':>8} {'eff@wall':>9} {'eff_unicity':>12} {'gap':>6}")
WALL = {"12-bit/2557": 8, "12-bit/2677": 4}     # from Thread 23a X3, variants A and B
for label, curve in curves:
    p, b, n, lam, G = curve
    K2 = math.isqrt(n) + 1
    m = 12
    k1w = WALL[label]
    eff = k1w * K2 / n
    eff_u = n ** (-1.0 / m)
    print(f"{label:<14} {lam_star(lam,n):>7.4f} {n:>6} {K2:>4} {m:>3} "
          f"{k1w:>8} {eff:>9.4f} {eff_u:>12.4f} {eff_u/eff:>5.2f}x")
print("\n  Both walls sit well below the unicity bound: d is information-theoretically")
print("  determined at the wall, so 'information-theoretic' cannot be the explanation")
print("  unless the BDD instance is itself ambiguous.  Y2 tests that directly.")

# --- Y2: exact CVP -----------------------------------------------------------
print("\n" + "-" * 78)
print("Y2: EXACT CVP (enumeration) on L0, rank 2m, m=8 -> dim 16")
print("-" * 78)
print("planted = the exact closest vector IS (k1_i*S_K1, k2_i*S_K2)?")
print("If planted==True where Babai/LLL fail  -> reduction-quality wall.")
print("If planted==False                      -> the BDD instance is ambiguous.\n")

M_EXACT = 8                       # dim 2m = 16; exact enumeration is comfortable here
K1_EXACT = [2, 4, 6, 8, 12, 16]

def exact_cvp_trial(curve, m, k1_bound, seed):
    """Returns (exact_is_planted, babai_ok, dist_exact/dist_planted)."""
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    d_secret = random.Random(seed + 7777).randint(1, n - 1)
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return None
    S_K1, _S_D, S_K2, _S = scales(n, k1_bound, k2_bound)

    M = build_homogeneous_lattice(sigs, n, lam, k1_bound, k2_bound)
    A = IntegerMatrix.from_matrix(M)
    LLL.reduction(A)
    rows = [[A[i][j] for j in range(2 * m)] for i in range(A.nrows)]
    rows = [r for r in rows if any(r)]
    B = IntegerMatrix.from_matrix(rows)

    t_pos = [0] * (2 * m)
    for i in range(m):
        t_pos[i] = (sigs[i]['A'] % n) * S_K1
    tgt = tuple(-x for x in t_pos)                       # u ~ CVP(L0, -t_pos)

    u_exact = list(CVP.closest_vector(B, tgt))
    v_exact = [t_pos[j] + u_exact[j] for j in range(2 * m)]
    v_planted = [0] * (2 * m)
    for i in range(m):
        v_planted[i] = sigs[i]['k1'] * S_K1
        v_planted[m + i] = sigs[i]['k2'] * S_K2

    is_planted = (v_exact == v_planted)
    d_exact = None
    if all(v_exact[i] % S_K1 == 0 and v_exact[m + i] % S_K2 == 0 for i in range(m)):
        cands = {d_from_k(v_exact[i] // S_K1, v_exact[m + i] // S_K2,
                          sigs[i], n, lam) for i in range(m)}
        if len(cands) == 1:
            d_exact = cands.pop()
    ratio = norm(v_exact) / norm(v_planted) if norm(v_planted) else float('nan')

    # variant-B (LLL + Kannan) result on the SAME instance, for the comparison
    Mb = build_projected_lattice(sigs, n, lam, k1_bound, k2_bound)
    Ab = IntegerMatrix.from_matrix(Mb)
    LLL.reduction(Ab)
    rb = [[Ab[i][j] for j in range(2 * m + 1)] for i in range(Ab.nrows)]
    lll_ok = (recover_d_projected(rb, sigs, n, lam, k1_bound, k2_bound) == d_secret)

    return is_planted, (d_exact == d_secret), lll_ok, ratio

print(f"{'curve':<14} {'K1':>3} {'eff':>6} {'exact=planted':>14} "
      f"{'exactCVP d ok':>14} {'LLL(B) d ok':>12} {'|ve|/|vp|':>10}")
for label, curve in curves:
    p, b, n, lam, G = curve
    K2 = math.isqrt(n) + 1
    for k1 in K1_EXACT:
        pl = cv = ll = 0
        rats, tot = [], 0
        for s in SEEDS:
            r = exact_cvp_trial(curve, M_EXACT, k1, s)
            if r is None:
                continue
            tot += 1
            pl += r[0]; cv += r[1]; ll += r[2]; rats.append(r[3])
        mr = sum(rats) / len(rats) if rats else float('nan')
        print(f"{label:<14} {k1:>3} {k1*K2/n:>6.3f} {str(pl)+'/'+str(tot):>14} "
              f"{str(cv)+'/'+str(tot):>14} {str(ll)+'/'+str(tot):>12} {mr:>10.4f}")

print(f"\n  (unicity bound at m={M_EXACT}: eff < n^(-1/{M_EXACT}) = "
      f"{2647 ** (-1.0/M_EXACT):.4f})")

# --- Y3: BKZ ladder ----------------------------------------------------------
print("\n" + "-" * 78)
print("Y3: BKZ ladder on variant B, m=12, the lam*=0.07 curve")
print("-" * 78)

def rate_bkz(curve, m, k1_bound, seeds, beta):
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    wins = 0
    for seed in seeds:
        d_secret = random.Random(seed + 7777).randint(1, n - 1)
        sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed)
        if len(sigs) < m:
            continue
        M = build_projected_lattice(sigs, n, lam, k1_bound, k2_bound)
        A = IntegerMatrix.from_matrix(M)
        LLL.reduction(A)
        if beta > 2:
            BKZ.reduction(A, BKZ.Param(beta))
        rows = [[A[i][j] for j in range(2 * m + 1)] for i in range(A.nrows)]
        wins += (recover_d_projected(rows, sigs, n, lam, k1_bound, k2_bound)
                 == d_secret)
    return wins, len(seeds)

fail = [c for lbl, c in curves if c[2] == 2647][0]
good = [c for lbl, c in curves if c[2] == 2659][0]
BETAS = [2, 10, 20, 30, 40]
for lbl, cur, grid in (("12-bit/2677 (lam*=0.07)", fail, [4, 6, 8, 12]),
                       ("12-bit/2557 (lam*=0.34)", good, [8, 12, 16, 24])):
    print(f"\n  {lbl}   (beta=2 means plain LLL)")
    print(f"  {'beta':>5} " + " ".join(f"K1={k:<5}" for k in grid))
    for beta in BETAS:
        cells = []
        for k1 in grid:
            w, t = rate_bkz(cur, 12, k1, SEEDS, beta)
            cells.append(f"{w}/{t}  ")
        print(f"  {beta:>5} " + " ".join(cells))

print("\n" + "=" * 78)
print("done")
print("=" * 78)
