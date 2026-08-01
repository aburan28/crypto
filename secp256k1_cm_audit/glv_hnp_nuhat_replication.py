"""
GLV-HNP Thread 23A, addendum: blind replication of nu_hat, and the (K1, m)
generalisation grid that the 2026-07-29 run asked for.

Why this file exists
--------------------
Two autolab runs happened on 2026-07-29 (commits d525931 and e845207).  The
merge at 7b5b702 / e31dd19 kept only the d525931 log entry, so the e845207
entry -- which introduced

    nu_hat = lambda_1(L2) / sqrt(det L2),
    L2 = < (n*S_K1, 0), (-lam*S_K1, S_K2) >,   det L2 = n*S_K1*S_K2

as a Phase-2 separator with AUC 0.935 -- is present in git but ABSENT from
RESEARCH_AUTOLAB_LOG.md at HEAD.  Today's Thread 23A session rediscovered the
same quantity blind (as mu/||v_planted|| in glv_hnp_phase2_projected.py, EXP
E5) and, crucially, the same counter-intuitive SIGN: a SHORT rival vector makes
the attack EASIER.

That makes today's E5 an independent replication rather than a new result, and
this script turns it into a clean one, plus runs the follow-up that e845207
explicitly proposed and never executed:

  Part A  nu_hat on the 17-bit E5 grid (10 curves x 4 eff levels).  Does the
          "C2 ceiling" nu_hat <= 0.645 measured at 20 bits hold at 17 bits?
  Part B  the generalisation grid: K1 in {36, 72, 144} x m in {10, 12, 14} at
          20 bits.  e845207's own words: "if the C2 ceiling stays near
          nu_hat ~ 0.65 across all nine cells, the cut is a property of the
          lattice family, not of the sample."
          FALSIFIER: if the ceiling drifts by more than ~0.1 across the nine
          cells, nu_hat is a per-(K1,m) statistic, not a family invariant.

Run: python3 glv_hnp_nuhat_replication.py
"""

import math
import random
import sympy
from fpylll import IntegerMatrix, LLL


# glv_hnp_phase2_projected.py runs its experiments at import time, so pull the
# helpers in by exec'ing only the part above the experiment marker.
_src = open(__file__.rsplit("/", 1)[0] + "/glv_hnp_phase2_projected.py").read()
_head = _src.split("# " + "=" * 75)[0]
exec(compile(_head, "glv_hnp_phase2_projected.py[helpers]", "exec"))


def nu_hat(n, lam, k1_bound, k2_bound):
    """lambda_1(L2)/sqrt(det L2); det L2 = n*S_K1*S_K2 is independent of lam."""
    S_K1, _S_D, S_K2, _ = scales(n, k1_bound, k2_bound)
    mu = lambda_block_mu(n, lam, S_K1, S_K2)
    return mu / math.sqrt(n * S_K1 * S_K2)


def auc(pos, neg):
    """Mann-Whitney AUC for 'small value => positive'; 1.0 = perfect."""
    if not pos or not neg:
        return float('nan')
    wins = 0.0
    for a in pos:
        for b in neg:
            wins += 1.0 if a < b else (0.5 if a == b else 0.0)
    return wins / (len(pos) * len(neg))


def best_cut(rows, key, label):
    """Best 'predict positive iff value <= thr' cut."""
    vals = sorted({r[key] for r in rows})
    best = (0, None)
    for thr in vals:
        acc = sum(1 for r in rows if (r[key] <= thr) == r[label])
        if acc > best[0]:
            best = (acc, thr)
    return best


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


def collect_curves(lo, hi, want, seed=1):
    """Random j=0 GLV curves with p in [lo,hi), n prime, n = 1 mod 3."""
    out = []
    rng = random.Random(seed)
    tries = 0
    while len(out) < want and tries < 20000:
        tries += 1
        p = int(sympy.nextprime(rng.randint(lo, hi)))
        if p >= hi or p % 3 != 1:
            continue
        eis = eisenstein_decompose(p)
        if eis is None:
            continue
        a_e, b_e = eis
        for t in j0_traces(a_e, b_e):
            n_c = p + 1 - t
            if n_c < 2 or n_c % 3 != 1 or not sympy.isprime(n_c):
                continue
            roots = glv_roots(n_c)
            if roots is None:
                continue
            cur = build_curve(p, n_c)
            if cur is None:
                continue
            out.append((p, cur[1], n_c, roots[0], cur[4]))
            break
    return out


def kan_wins(curve, m, k1_bound, seeds):
    """Success count of the baseline (2m+2) Kannan attack."""
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    w = 0
    for s in seeds:
        d = random.Random(s + 7777).randint(1, n - 1)
        sigs = gen_signatures(G, d, m, n, lam, p, k1_bound, k2_bound, s)
        if len(sigs) < m:
            continue
        ok, _, _ = attack_kannan(sigs, n, lam, k1_bound, k2_bound, d)
        w += bool(ok)
    return w


print("=" * 82)
print("Thread 23A addendum — nu_hat replication + (K1, m) generalisation grid")
print("=" * 82)
print("nu_hat = lambda_1(L2)/sqrt(det L2),  L2 = <(n S_K1,0), (-lam S_K1, S_K2)>")
print("Reference (2026-07-29, commit e845207, 20-bit, K1=72, m=12, 100 curves):")
print("  AUC 0.935, best acc 89.0%, C2 ceiling nu_hat <= 0.645 (no curve above")
print("  0.645 ever reached full recovery).\n")

SEEDS = [42, 1234, 9999, 555, 31337]

# ---------------------------------------------------------------------------
print("-" * 82)
print("PART A — 17-bit replication on today's EXP E5 grid")
print("-" * 82)

LO, HI = 2 ** 16, 2 ** 17

def search_curves_deciles(lo, hi, per_bin=1, nbins=10):
    bins = {i: [] for i in range(nbins)}
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
                    idx = min(nbins - 1, int(lam_star(lam, n_c) / (0.5 / nbins)))
                    if len(bins[idx]) >= per_bin:
                        continue
                    cur = build_curve(p, n_c)
                    if cur is None:
                        continue
                    bins[idx].append((p, cur[1], n_c, lam, cur[4]))
        if all(len(v) >= per_bin for v in bins.values()):
            break
        p = int(sympy.nextprime(p))
    return [c for i in range(nbins) for c in bins[i]]

curves17 = search_curves_deciles(LO, HI, per_bin=1, nbins=10)
print(f"{len(curves17)} 17-bit curves (one per lam* decile), m=12, "
      f"{len(SEEDS)} seeds\n")

M_A = 12
rowsA = []
for eff_t in (0.05, 0.10, 0.15, 0.20):
    print(f"  === eff target {eff_t} ===")
    print("    " + f"{'p':>7} {'n':>7} {'lam*':>6} {'K1':>4} {'eff':>6} "
          f"{'nu_hat':>7} {'KAN':>6}")
    cells = []
    for c in curves17:
        p, b, n, lam, G = c
        k2 = math.isqrt(n) + 1
        k1 = max(2, int(eff_t * n / k2))
        nh = nu_hat(n, lam, k1, k2)
        w = kan_wins(c, M_A, k1, SEEDS)
        r = {'p': p, 'n': n, 'lam_star': lam_star(lam, n), 'k1': k1,
             'eff': k1 * k2 / n, 'nu_hat': nh, 'wins': w,
             'C2': w == len(SEEDS), 'eff_t': eff_t}
        rowsA.append(r); cells.append(r)
        print("    " + f"{p:>7} {n:>7} {r['lam_star']:>6.3f} {k1:>4} "
              f"{r['eff']:>6.3f} {nh:>7.3f} {str(w)+'/'+str(len(SEEDS)):>6}")
    pos = [r['nu_hat'] for r in cells if r['C2']]
    neg = [r['nu_hat'] for r in cells if not r['C2']]
    if pos and neg:
        print(f"    -> AUC {auc(pos, neg):.3f}   C2 max nu_hat {max(pos):.3f}   "
              f"non-C2 min nu_hat {min(neg):.3f}   "
              f"separated={max(pos) < min(neg)}")
    else:
        print(f"    -> degenerate: {len(pos)}/{len(cells)} C2")
    print()

posA = [r['nu_hat'] for r in rowsA if r['C2']]
negA = [r['nu_hat'] for r in rowsA if not r['C2']]
accA, thrA = best_cut(rowsA, 'nu_hat', 'C2')
acc_ls, thr_ls = best_cut([dict(r, neg_ls=-r['lam_star']) for r in rowsA],
                          'neg_ls', 'C2')
print(f"POOLED PART A ({len(rowsA)} cells, {sum(r['C2'] for r in rowsA)} C2):")
print(f"  AUC(nu_hat)      = {auc(posA, negA):.3f}")
print(f"  best cut nu_hat <= {thrA:.3f} -> {accA}/{len(rowsA)} = "
      f"{accA/len(rowsA):.1%}   (baseline "
      f"{max(sum(r['C2'] for r in rowsA), len(rowsA)-sum(r['C2'] for r in rowsA))}"
      f"/{len(rowsA)})")
print(f"  C2 nu_hat range     [{min(posA):.3f}, {max(posA):.3f}]")
print(f"  non-C2 nu_hat range [{min(negA):.3f}, {max(negA):.3f}]")
print(f"  20-bit C2 ceiling was 0.645; here it is {max(posA):.3f}")

# ---------------------------------------------------------------------------
print("\n" + "-" * 82)
print("PART B — (K1, m) generalisation grid at 20 bits")
print("-" * 82)
print("e845207's unexecuted sub-task: does the C2 ceiling stay near 0.65")
print("across K1 in {36,72,144} x m in {10,12,14}?\n")

N_CURVES_B = 30
SEEDS_B = [42, 1234, 9999, 555, 31337, 8675309]
curves20 = collect_curves(2 ** 19, 2 ** 20, N_CURVES_B, seed=20260801)
print(f"{len(curves20)} fresh 20-bit j=0 GLV curves, {len(SEEDS_B)} seeds/cell\n")

print(f"{'K1':>5} {'m':>4} {'eff':>6} {'C2':>7} {'AUC':>6} "
      f"{'C2 ceil':>8} {'nonC2 min':>10} {'best cut':>9} {'acc':>7}")
ceilings = []
for k1 in (36, 72, 144):
    for m in (10, 12, 14):
        cells = []
        effs = []
        for c in curves20:
            p, b, n, lam, G = c
            k2 = math.isqrt(n) + 1
            effs.append(k1 * k2 / n)
            w = kan_wins(c, m, k1, SEEDS_B)
            cells.append({'nu_hat': nu_hat(n, lam, k1, k2), 'wins': w,
                          'C2': w == len(SEEDS_B)})
        pos = [r['nu_hat'] for r in cells if r['C2']]
        neg = [r['nu_hat'] for r in cells if not r['C2']]
        if pos and neg:
            a = auc(pos, neg)
            acc, thr = best_cut(cells, 'nu_hat', 'C2')
            ceilings.append((k1, m, max(pos)))
            print(f"{k1:>5} {m:>4} {sum(effs)/len(effs):>6.3f} "
                  f"{str(len(pos))+'/'+str(len(cells)):>7} {a:>6.3f} "
                  f"{max(pos):>8.3f} {min(neg):>10.3f} {thr:>9.3f} "
                  f"{acc/len(cells):>7.1%}")
        else:
            print(f"{k1:>5} {m:>4} {sum(effs)/len(effs):>6.3f} "
                  f"{str(len(pos))+'/'+str(len(cells)):>7} "
                  f"{'degen':>6} {'-':>8} {'-':>10} {'-':>9} {'-':>7}")

if ceilings:
    cv = [c for _, _, c in ceilings]
    print(f"\nC2 ceiling across {len(ceilings)} non-degenerate cells: "
          f"min {min(cv):.3f}, max {max(cv):.3f}, spread {max(cv)-min(cv):.3f}")
    print("FALSIFIER: spread > ~0.10 => nu_hat is a per-(K1,m) statistic, "
          "not a family invariant.")

print("\nDone.")
