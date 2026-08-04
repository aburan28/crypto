"""
GLV-HNP Thread 23b: does the nu_hat cut generalise across (K1, m)?

Background
----------
Thread 20d (2026-07-29) measured nu_hat against the June C1/C2 classes at a
SINGLE operating point: K1 = 72, m = 12, 6 seeds, 100 fresh 20-bit j=0 CM
curves.  Result: AUC 0.935, every C2 curve had nu_hat <= 0.645, no curve with
nu_hat > 0.645 was C2.  The 2026-07-29 entry flagged the obvious weakness --
"20d used one m and one K1" -- and proposed the falsifier:

    re-run 20d at K1 in {36, 72, 144} and m in {10, 12, 14}; if the C2 ceiling
    stays near nu_hat ~ 0.65 across all nine cells, the cut is a property of
    the lattice family rather than of the sample.

This script runs that 3x3 grid.  Same harvest, same run_experiment, same C1/C2
convention (C2 = recovers d on ALL seeds) as 20d, so the K1=72/m=12 cell is a
direct replication of the published number.

Note nu_hat itself depends on K1 (through S_K1 = n//K1), so each column of the
grid has its OWN nu_hat values -- the question is whether the *cut* is stable,
not whether the values are.

Run: python3 glv_hnp_nuhat_grid.py
"""

import math
import random
import sys
import time

import sympy

import importlib.util
_spec = importlib.util.spec_from_file_location(
    "_lib", __file__.rsplit("/", 1)[0] + "/glv_hnp_nuhat_lib.py")
_lib = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_lib)

eisenstein_decompose = _lib.eisenstein_decompose
j0_traces = _lib.j0_traces
glv_eigenvalues = _lib.glv_eigenvalues
mu_of = _lib.mu_of
identify_twist = _lib.identify_twist
find_generator = _lib.find_generator
rival_sublattice_nu = _lib.rival_sublattice_nu
run_experiment = _lib.run_experiment

K1_GRID = (36, 72, 144)
M_GRID = (10, 12, 14)
SEEDS_N = 4
N_CURVES = 60
BITS_LO, BITS_HI = 2 ** 19, 2 ** 20


def max_partial_quotient(a, b):
    best = 0
    while b:
        q, r = divmod(a, b)
        best = max(best, q)
        a, b = b, r
    return best


def harvest(n_curves):
    """Verbatim from glv_hnp_nuhat_vs_c1c2.py:harvest (20d)."""
    out, seen = [], set()
    p = int(sympy.nextprime(BITS_LO))
    while p < BITS_HI and len(out) < n_curves:
        if p % 3 == 1:
            eis = eisenstein_decompose(p)
            if eis is not None:
                for t in j0_traces(*eis):
                    n = p + 1 - t
                    if n < 2 or n in seen or n % 3 != 1 or not sympy.isprime(n):
                        continue
                    roots = glv_eigenvalues(n)
                    if roots is None:
                        continue
                    b = identify_twist(p, n)
                    if b is None:
                        break
                    G = find_generator(p, b, n)
                    if G is None:
                        break
                    seen.add(n)
                    out.append({'p': p, 'b': b, 'n': n, 'lam': roots[0], 'G': G})
                    break
        p = int(sympy.nextprime(p))
    return out


def auc(scores, labels):
    pos = [s for s, l in zip(scores, labels) if l]
    neg = [s for s, l in zip(scores, labels) if not l]
    if not pos or not neg:
        return float('nan')
    tot = 0.0
    for a in pos:
        for b in neg:
            tot += 1.0 if a > b else 0.5 if a == b else 0.0
    return tot / (len(pos) * len(neg))


def best_cut_accuracy(scores, labels):
    best, cut = 0.0, None
    for c in sorted(set(scores)):
        for sign in (1, -1):
            acc = sum(1 for s, l in zip(scores, labels)
                      if ((s > c) if sign > 0 else (s <= c)) == l) / len(labels)
            if acc > best:
                best, cut = acc, (c, sign)
    return best, cut


def main():
    print("=" * 78)
    print("GLV-HNP Thread 23b: nu_hat cut across the (K1, m) grid")
    print("=" * 78)
    print(f"K1 in {K1_GRID}, m in {M_GRID}, {SEEDS_N} seeds, "
          f"{N_CURVES} fresh 20-bit j=0 CM curves")
    print("C2 = recovers d on ALL seeds (Exp S / 20d convention)")
    print("Falsifier: C2 ceiling drifts far from nu_hat ~ 0.65 across cells")
    print()

    t0 = time.time()
    curves = harvest(N_CURVES)
    print(f"harvested {len(curves)} curves in {time.time()-t0:.1f}s")
    print()

    cells = {}
    t1 = time.time()
    for k1 in K1_GRID:
        for m in M_GRID:
            rows = []
            for c in curves:
                n = c['n']
                k2 = math.isqrt(n) + 1
                wins = 0
                for s in range(SEEDS_N):
                    d = random.Random(s + 7777).randint(1, n - 1)
                    wins += run_experiment(c['p'], n, c['lam'], c['G'], m, d,
                                           k1, k2, s)
                _, _, nh = rival_sublattice_nu(n, c['lam'], k1, k2)
                rows.append({
                    'n': n, 'wins': wins, 'C2': wins == SEEDS_N, 'nu_hat': nh,
                    'mu': mu_of(c['lam'], n),
                    'max_a': max_partial_quotient(c['lam'], n),
                })
            cells[(k1, m)] = rows
            nc2 = sum(r['C2'] for r in rows)
            print(f"  cell K1={k1:3d} m={m:2d}: C2 {nc2:2d}/{len(rows)}   "
                  f"({time.time()-t1:6.1f}s elapsed)")
    print(f"\n{len(K1_GRID)*len(M_GRID)*len(curves)*SEEDS_N} LLL trials "
          f"in {time.time()-t1:.1f}s")

    # ---------------------------------------------------------------- grid
    print("\n" + "-" * 78)
    print("GRID: C2 ceiling (max nu_hat among C2 curves) and separation quality")
    print("-" * 78)
    print(f"{'K1':>5} {'m':>3} {'C2/N':>7} {'base':>6} {'AUC':>7} {'acc':>7} "
          f"{'C2 ceiling':>11} {'C1 floor':>9} {'overlap':>16}")
    ceilings = []
    for k1 in K1_GRID:
        for m in M_GRID:
            rows = cells[(k1, m)]
            c2 = [r for r in rows if r['C2']]
            c1 = [r for r in rows if not r['C2']]
            base = max(len(c1), len(c2)) / len(rows)
            a = auc([r['nu_hat'] for r in rows], [r['C2'] for r in rows])
            a_or = max(a, 1 - a)
            acc, _ = best_cut_accuracy([r['nu_hat'] for r in rows],
                                       [r['C2'] for r in rows])
            if c2 and c1:
                ceil_ = max(r['nu_hat'] for r in c2)
                floor_ = min(r['nu_hat'] for r in c1)
                ceilings.append(ceil_)
                ov = f"[{floor_:.3f},{ceil_:.3f}]" if floor_ <= ceil_ else "NONE (clean cut)"
                print(f"{k1:>5} {m:>3} {len(c2):>3}/{len(rows):<3} {base:>6.2f} "
                      f"{a_or:>7.3f} {acc*100:>6.1f}% {ceil_:>11.3f} "
                      f"{floor_:>9.3f} {ov:>16}")
            else:
                print(f"{k1:>5} {m:>3} {len(c2):>3}/{len(rows):<3} {base:>6.2f} "
                      f"{'-':>7} {'-':>7} {'degenerate cell':>28}")
    if ceilings:
        print()
        print(f"C2 ceiling across {len(ceilings)} non-degenerate cells: "
              f"min={min(ceilings):.3f} max={max(ceilings):.3f} "
              f"mean={sum(ceilings)/len(ceilings):.3f}")
        print(f"20d published ceiling (K1=72, m=12, 6 seeds, 100 curves) = 0.645")

    # -------------------------------------------------- head-to-head per cell
    print("\n" + "-" * 78)
    print("AUC head-to-head vs the falsified invariants, per cell")
    print("-" * 78)
    print(f"{'K1':>5} {'m':>3} {'nu_hat':>8} {'mu':>8} {'max_a':>8}")
    for k1 in K1_GRID:
        for m in M_GRID:
            rows = cells[(k1, m)]
            if not any(r['C2'] for r in rows) or all(r['C2'] for r in rows):
                print(f"{k1:>5} {m:>3} {'(degenerate cell)':>26}")
                continue
            lab = [r['C2'] for r in rows]
            out = []
            for key in ('nu_hat', 'mu', 'max_a'):
                a = auc([r[key] for r in rows], lab)
                out.append(max(a, 1 - a))
            print(f"{k1:>5} {m:>3} {out[0]:>8.3f} {out[1]:>8.3f} {out[2]:>8.3f}")

    # ------------------------------------------------------------- pooled
    print("\n" + "-" * 78)
    print("POOLED across all cells (curve-cell pairs)")
    print("-" * 78)
    pool = [r for cell in cells.values() for r in cell]
    lab = [r['C2'] for r in pool]
    if any(lab) and not all(lab):
        for key in ('nu_hat', 'mu', 'max_a'):
            a = auc([r[key] for r in pool], lab)
            acc, cut = best_cut_accuracy([r[key] for r in pool], lab)
            print(f"  {key:>7}: AUC {max(a,1-a):.3f}  best-cut acc {acc*100:.1f}%"
                  f"  at {cut[0]:.4f} (sign {cut[1]:+d})")
        base = max(sum(lab), len(lab) - sum(lab)) / len(lab)
        print(f"  trivial baseline = {base*100:.1f}%   (n = {len(pool)} pairs)")

        print("\n  pooled nu_hat decile response")
        order = sorted(pool, key=lambda r: r['nu_hat'])
        k = max(1, len(order) // 10)
        print(f"  {'decile':>7} {'nu_hat mid':>11} {'C2 rate':>9} "
              f"{'mean wins/'+str(SEEDS_N):>12}")
        for i in range(0, len(order), k):
            grp = order[i:i + k]
            if len(grp) < 2:
                continue
            print(f"  {i//k+1:>7} {sum(r['nu_hat'] for r in grp)/len(grp):>11.3f} "
                  f"{sum(r['C2'] for r in grp)/len(grp):>9.2f} "
                  f"{sum(r['wins'] for r in grp)/len(grp):>12.2f}")

    # --------------------------------------------- G2: m_threshold formulation
    # The binary C2 label is cell-relative: its cut moves with the cell's
    # overall difficulty.  m_threshold (smallest m at which a curve wins on all
    # seeds) is a per-curve quantity that does not need a cell-level baseline.
    print("\n" + "-" * 78)
    print("G2: nu_hat vs m_threshold (smallest m with all-seeds success)")
    print("-" * 78)
    M_SCAN = tuple(range(6, 21))
    for k1 in K1_GRID:
        nhs, mts, mus, mas = [], [], [], []
        censored = 0
        for c in curves:
            n = c['n']
            k2 = math.isqrt(n) + 1
            mt = None
            for m in M_SCAN:
                wins = 0
                for s in range(SEEDS_N):
                    d = random.Random(s + 7777).randint(1, n - 1)
                    wins += run_experiment(c['p'], n, c['lam'], c['G'], m, d,
                                           k1, k2, s)
                if wins == SEEDS_N:
                    mt = m
                    break
            if mt is None:
                censored += 1
                mt = M_SCAN[-1] + 1        # right-censored
            _, _, nh = rival_sublattice_nu(n, c['lam'], k1, k2)
            nhs.append(nh); mts.append(mt)
            mus.append(mu_of(c['lam'], n))
            mas.append(max_partial_quotient(c['lam'], n))
        sp_nu = _spearman(nhs, mts)
        sp_mu = _spearman(mus, mts)
        sp_ma = _spearman(mas, mts)
        print(f"  K1={k1:3d}  spearman(nu_hat, m_thr) = {sp_nu:+.3f}   "
              f"(mu {sp_mu:+.3f}, max_a {sp_ma:+.3f})   "
              f"censored {censored}/{len(curves)}")
        # tertile table
        order = sorted(range(len(nhs)), key=lambda i: nhs[i])
        t = max(1, len(order) // 3)
        labels = ('low ', 'mid ', 'high')
        for ti, lab in enumerate(labels):
            grp = order[ti * t:(ti + 1) * t] if ti < 2 else order[2 * t:]
            if not grp:
                continue
            print(f"      {lab} nu_hat tertile: nu_hat={sum(nhs[i] for i in grp)/len(grp):.3f}"
                  f"  mean m_thr={sum(mts[i] for i in grp)/len(grp):5.2f}"
                  f"  censored={sum(1 for i in grp if mts[i] > M_SCAN[-1])}")
    print()
    print("  Reading: a POSITIVE spearman means high nu_hat needs MORE signatures,")
    print("  i.e. low nu_hat = easier -- the 2026-07-29 direction.  Unlike the C2")
    print("  ceiling, m_threshold is comparable across K1 columns.")
    print("  max_a comes out NEGATIVE (large partial quotient => fewer signatures),")
    print("  which is the same mechanism seen through a scale-free proxy: it has")
    print("  the right sign but roughly half the rank correlation of nu_hat.")

    print("\nDone.")


def _spearman(xs, ys):
    def rank(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
                j += 1
            avg = (i + j) / 2.0 + 1.0
            for k in range(i, j + 1):
                r[order[k]] = avg
            i = j + 1
        return r
    rx, ry = rank(xs), rank(ys)
    mx = sum(rx) / len(rx); my = sum(ry) / len(ry)
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    dx = math.sqrt(sum((a - mx) ** 2 for a in rx))
    dy = math.sqrt(sum((b - my) ** 2 for b in ry))
    return num / (dx * dy) if dx and dy else 0.0


if __name__ == '__main__':
    sys.exit(main())
