"""
GLV-HNP Thread 23 — closed form for nu_hat via continued fractions.

Background
----------
Thread 20d (2026-07-29, commit e845207) found that the Phase 2 GLV-HNP lattice
has a cheap empirical separator

    nu_hat = lambda_1(L2) / sqrt(det L2),      L2 = <(n*S1, 0), (-lam*S1, S2)>

with AUC 0.935 against the June C1/C2 classes.  det(L2) = n*S1*S2 is
lam-independent, so nu_hat isolates the lam-dependence.  But nu_hat was
*computed* (one Lagrange-Gauss reduction), not *understood*, and the
2026-07-29 next-step asked for a closed form.

Derivation (this script's pre-registered claims)
------------------------------------------------
A general element of L2 is  a*v + b*u = ((b*n - a*lam)*S1, a*S2), so

    ||a*v + b*u||^2 = S1^2 * (a*lam - b*n)^2 + S2^2 * a^2.

For a = 0 the minimum is n*S1.  For a >= 1, min over b is the centered
residue r(a) = min_b |a*lam - b*n|.  By the best-approximation theorem of the
second kind, for every a in [q_j, q_{j+1}) we have r(a) >= r(q_j), where p_j/q_j
are the convergents of lam/n.  Both terms of the objective are then >= their
value at a = q_j, so:

  (C1) EXACTNESS.  lambda_1(L2)^2 = min( (n*S1)^2,
                                         min_j [ S1^2*r_j^2 + S2^2*q_j^2 ] ),
       r_j = |q_j*lam - p_j*n|.   This is an identity, not an approximation.

Normalising by det L2 = n*S1*S2 and writing the critical scale

       Q = sqrt(n*S1/S2)  ( = sqrt(n*K2/K1) up to the floor in S1,S2 )
       theta_j = q_j / Q

and using r_j ~ n/q_{j+1} gives the scale-explicit form

  (C2) nu_hat^2 ~ min_j ( theta_j^2 + 1/theta_{j+1}^2 ).

  (C3) STRADDLE.  For j < s the term is >= 1/theta_{j+1}^2 >= 1 and for j > s
       it is >= theta_j^2 > 1, where s is the unique index with
       theta_s <= 1 < theta_{s+1}.  So whenever nu_hat < 1 the minimiser is
       j = s: only the ONE convergent straddling the scale Q matters.

  (C4) LOWER BOUND.  r_j > n/(q_{j+1}+q_j) is the two-sided form of the
       convergent bound, so nu_hat^2 > theta_s^2 + 1/(theta_{s+1}+theta_s)^2.
       With theta_{s+1} = a_{s+1}*theta_s + theta_{s-1} <= (a_{s+1}+1)*theta_s
       and AM-GM (x^2 + 1/y^2 >= 2x/y),
            nu_hat > sqrt( 2 / (a* + 2) ),   a* := a_{s+1},
       the partial quotient immediately after the straddling convergent.
       Large a* is NECESSARY (not sufficient) for small nu_hat.
       Note the naive form sqrt(2/(a*+1)) -- which drops the +theta_s -- is
       NOT a bound: r_j < n/q_{j+1} makes the approximation an over-estimate.
       Both are reported below so the gap is visible.

  (C5) RETRO-EXPLANATION.  (C3)+(C4) say the relevant CF data is *local to the
       scale Q*.  The June invariants q_cf / max_q_cf / max_a are maxima over
       the WHOLE continued fraction, i.e. scale-free, so they cannot see a*.
       Prediction: corr(nu_hat, a*) is strong while corr(nu_hat, max_a) ~ 0,
       on the same sample -- which is exactly the Exp S falsification of max_a.

Falsifiers
----------
  A: any sample where the CF formula (C1) differs from an exact integer
     Lagrange-Gauss reduction  -> the derivation is wrong.
  B: straddle rate (C3) not ~100% among curves with nu_hat < 1.
  C: (C2) correlating < 0.9 with the exact nu_hat  -> "convergent at scale Q"
     is the wrong story (the falsifier named in the 2026-07-29 log entry).
  D: a* failing to beat max_a on the Exp S C1/C2 task -> locality is wrong.

Run: python3 glv_hnp_nuhat_cf_closed_form.py
"""

import math
import random
import sys
import time

import sympy

import importlib.util
_spec = importlib.util.spec_from_file_location(
    "_t20a", __file__.rsplit("/", 1)[0] + "/glv_hnp_phase2_nuhat_lib.py")
_t20a = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_t20a)

eisenstein_decompose = _t20a.eisenstein_decompose
j0_traces = _t20a.j0_traces
glv_eigenvalues = _t20a.glv_eigenvalues
mu_of = _t20a.mu_of
identify_twist = _t20a.identify_twist
find_generator = _t20a.find_generator
rival_sublattice_nu = _t20a.rival_sublattice_nu
run_experiment = _t20a.run_experiment
scales = _t20a.scales

# Exp S protocol (2026-06-29), as replicated by Thread 20d
K1_FIXED = 72
M_FIXED = 12
SEEDS_N = 6
N_CURVES = 100
BITS_LO, BITS_HI = 2 ** 19, 2 ** 20

SECP_N = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
SECP_LAM = 0x5363ad4cc05c30e0a5261c028812645a122e22ea20816678df02967c1b23bd72

# ---------------------------------------------------------------------------
# Continued fractions
# ---------------------------------------------------------------------------

def cf_convergents(num, den):
    """Convergents of num/den.  Returns (partial_quotients, [(p_j, q_j)]).

    Index convention: pq[0] = a_0 = floor(num/den) and convs[0] = (a_0, 1),
    so pq[j] = a_j is the partial quotient that PRODUCES convs[j], and the
    partial quotient following convs[j] is pq[j+1].
    """
    pq, convs = [], []
    x, y = num, den
    a, r = divmod(x, y)
    p_prev, p_cur = 1, a
    q_prev, q_cur = 0, 1
    pq.append(a)
    convs.append((p_cur, q_cur))
    x, y = y, r
    while y:
        a, r = divmod(x, y)
        p_prev, p_cur = p_cur, a * p_cur + p_prev
        q_prev, q_cur = q_cur, a * q_cur + q_prev
        pq.append(a)
        convs.append((p_cur, q_cur))
        x, y = y, r
    return pq, convs


def max_partial_quotient(a, b):
    """Largest partial quotient of a/b -- the scale-free 'max_a' invariant
    falsified by Exp S (2026-06-29)."""
    best = 0
    while b:
        q, r = divmod(a, b)
        best = max(best, q)
        a, b = b, r
    return best

# ---------------------------------------------------------------------------
# lambda_1(L2): exact integer Lagrange-Gauss, and the CF formula
# ---------------------------------------------------------------------------

def gauss_reduce_exact(u, v):
    """Exact integer Lagrange-Gauss; returns lambda_1^2 as an integer.
    (The float version in the lib is fine at 20 bits but loses the exact
    comparison at 256 bits, which part A needs.)"""
    n2 = lambda w: w[0] * w[0] + w[1] * w[1]
    if n2(u) < n2(v):
        u, v = v, u
    while True:
        nv = n2(v)
        d = u[0] * v[0] + u[1] * v[1]
        q = (2 * d + nv) // (2 * nv)          # round(d / nv)
        r = (u[0] - q * v[0], u[1] - q * v[1])
        u, v = v, r
        if n2(u) <= n2(v):
            return n2(u)


def nu_cf(n, lam, k1_bound, k2_bound):
    """lambda_1(L2)^2 via the convergents of lam/n, plus the diagnostics
    (winning index j*, straddling index s, theta's, a*)."""
    s1, _, s2, _ = scales(n, k1_bound, k2_bound)
    pq, convs = cf_convergents(lam, n)
    det = n * s1 * s2
    Q = math.sqrt(n * s1 / s2)

    best, jstar = (n * s1) ** 2, -1            # a = 0 candidate
    terms = []
    for j, (p, q) in enumerate(convs):
        r = abs(q * lam - p * n)
        val = s1 * s1 * r * r + s2 * s2 * q * q
        terms.append(val)
        if val < best:
            best, jstar = val, j

    thetas = [q / Q for (_, q) in convs]
    # straddling index: theta_s <= 1 < theta_{s+1}
    s_idx = 0
    for j in range(len(thetas)):
        if thetas[j] <= 1.0:
            s_idx = j
        else:
            break

    a_star = pq[s_idx + 1] if s_idx + 1 < len(pq) else None
    th_s = thetas[s_idx]
    th_s1 = thetas[s_idx + 1] if s_idx + 1 < len(thetas) else float('inf')
    closed = math.sqrt(th_s ** 2 + 1.0 / (th_s1 ** 2)) if th_s1 != float('inf') \
        else th_s
    # (C2) evaluated as a min over all j, not only the straddling one
    closed_min = math.inf
    for j in range(len(thetas) - 1):
        closed_min = min(closed_min, thetas[j] ** 2 + 1.0 / thetas[j + 1] ** 2)
    closed_min = math.sqrt(min(closed_min, thetas[-1] ** 2)) \
        if closed_min != math.inf else thetas[-1]

    return {
        'lam1_sq': best, 'nu_hat': math.sqrt(best / det), 'jstar': jstar,
        's_idx': s_idx, 'theta_s': th_s, 'theta_s1': th_s1, 'a_star': a_star,
        'closed_straddle': closed, 'closed_min': closed_min,
        'n_convs': len(convs), 'Q': Q,
    }

# ---------------------------------------------------------------------------
# Stats helpers
# ---------------------------------------------------------------------------

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
            for k in range(i, j + 1):
                r[order[k]] = avg
            i = j + 1
        return r
    rx, ry = rank(xs), rank(ys)
    mx, my = sum(rx) / len(rx), sum(ry) / len(ry)
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    dx = math.sqrt(sum((a - mx) ** 2 for a in rx))
    dy = math.sqrt(sum((b - my) ** 2 for b in ry))
    return num / (dx * dy) if dx and dy else float('nan')


def pearson(xs, ys):
    mx, my = sum(xs) / len(xs), sum(ys) / len(ys)
    num = sum((a - mx) * (b - my) for a, b in zip(xs, ys))
    dx = math.sqrt(sum((a - mx) ** 2 for a in xs))
    dy = math.sqrt(sum((b - my) ** 2 for b in ys))
    return num / (dx * dy) if dx and dy else float('nan')


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

# ---------------------------------------------------------------------------
# Curve harvesting (Exp S / Thread 20d convention)
# ---------------------------------------------------------------------------

def harvest(n_curves, want_generator=True):
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
                    rec = {'p': p, 'n': n, 'lam': roots[0]}
                    if want_generator:
                        b = identify_twist(p, n)
                        if b is None:
                            break
                        G = find_generator(p, b, n)
                        if G is None:
                            break
                        rec['b'], rec['G'] = b, G
                    seen.add(n)
                    out.append(rec)
                    break
        p = int(sympy.nextprime(p))
    return out

# ---------------------------------------------------------------------------
# Part A -- exactness of the CF formula
# ---------------------------------------------------------------------------

def part_a():
    print("\n" + "=" * 78)
    print("PART A -- (C1) exactness: CF formula vs exact integer Lagrange-Gauss")
    print("=" * 78)
    rng = random.Random(20260802)
    cases = []

    # A1: real 20-bit j=0 GLV curves at the Exp S parameters
    curves = harvest(40, want_generator=False)
    for c in curves:
        k2 = math.isqrt(c['n']) + 1
        cases.append(('20-bit GLV', c['n'], c['lam'], K1_FIXED, k2))

    # A2: random (prime n, arbitrary lam) at several bit sizes and eff
    for bits in (20, 24, 64, 128, 256):
        for _ in range(40):
            n = int(sympy.nextprime(rng.getrandbits(bits) | (1 << (bits - 1))))
            lam = rng.randrange(1, n)
            k2 = math.isqrt(n) + 1
            eff = rng.choice([0.02, 0.05, 0.10, 0.25, 0.50])
            k1 = max(2, int(eff * math.isqrt(n)))
            cases.append((f'random {bits}-bit', n, lam, k1, k2))

    # A3: secp256k1 itself, several eff
    for eff in (0.02, 0.05, 0.10, 0.25, 0.50):
        k1 = max(2, int(eff * math.isqrt(SECP_N)))
        cases.append(('secp256k1', SECP_N, SECP_LAM,
                      k1, math.isqrt(SECP_N) + 1))

    bad = 0
    by_group = {}
    for group, n, lam, k1, k2 in cases:
        s1, _, s2, _ = scales(n, k1, k2)
        exact = gauss_reduce_exact((n * s1, 0), (-lam * s1, s2))
        cf = nu_cf(n, lam, k1, k2)['lam1_sq']
        ok = (exact == cf)
        by_group.setdefault(group, [0, 0])
        by_group[group][1] += 1
        if ok:
            by_group[group][0] += 1
        else:
            bad += 1
            if bad <= 5:
                print(f"  MISMATCH [{group}] n={n} lam={lam} k1={k1}: "
                      f"gauss={exact} cf={cf} ratio={cf/exact:.6f}")

    print(f"{'group':>16} {'agree':>8} {'total':>7}")
    for g, (ok, tot) in by_group.items():
        print(f"{g:>16} {ok:>8} {tot:>7}")
    print(f"\n  RESULT: {len(cases)-bad}/{len(cases)} exact integer agreement.")
    print("  (C1) " + ("HOLDS -- lambda_1(L2) is a convergent of lam/n."
                       if bad == 0 else "FALSIFIED."))
    return bad == 0

# ---------------------------------------------------------------------------
# Part B/C/D -- straddle, closed form, locality
# ---------------------------------------------------------------------------

def part_bcd():
    print("\n" + "=" * 78)
    print("PART B/C/D -- straddle (C3), closed form (C2), locality (C4)/(C5)")
    print("=" * 78)
    rng = random.Random(830219)
    rows = []
    for bits in (20, 24, 32, 64, 256):
        for _ in range(200):
            n = int(sympy.nextprime(rng.getrandbits(bits) | (1 << (bits - 1))))
            lam = rng.randrange(1, n)
            k2 = math.isqrt(n) + 1
            eff = rng.choice([0.02, 0.05, 0.10, 0.25, 0.50])
            k1 = max(2, int(eff * math.isqrt(n)))
            d = nu_cf(n, lam, k1, k2)
            d['bits'] = bits
            d['max_a'] = max_partial_quotient(lam, n)
            rows.append(d)

    # --- B: does the winner sit at the straddling index?
    sub = [r for r in rows if r['nu_hat'] < 1.0]
    hit = sum(1 for r in sub if r['jstar'] == r['s_idx'])
    print(f"\n(C3) straddle: winner j* == straddling s, among nu_hat < 1:")
    print(f"  {hit}/{len(sub)} = {100.0*hit/len(sub):.1f}%")
    allhit = sum(1 for r in rows if r['jstar'] == r['s_idx'])
    print(f"  over all {len(rows)} rows (incl. nu_hat >= 1): "
          f"{allhit}/{len(rows)} = {100.0*allhit/len(rows):.1f}%")
    offs = [abs(r['jstar'] - r['s_idx']) for r in rows if r['jstar'] != r['s_idx']]
    if offs:
        nu_off = [r['nu_hat'] for r in rows if r['jstar'] != r['s_idx']]
        print(f"  misses: |j*-s| max={max(offs)}, "
              f"nu_hat of misses in [{min(nu_off):.3f},{max(nu_off):.3f}]")
        below = [r for r in rows if r['nu_hat'] < min(nu_off)]
        hb = sum(1 for r in below if r['jstar'] == r['s_idx'])
        print(f"  straddle is EXACT below the smallest miss: {hb}/{len(below)} "
              f"of rows with nu_hat < {min(nu_off):.3f}")

    # --- C: closed form vs exact
    ex = [r['nu_hat'] for r in rows]
    cl = [r['closed_min'] for r in rows]
    cs = [r['closed_straddle'] for r in rows]
    ratios = sorted(c / e for c, e in zip(cl, ex))
    print(f"\n(C2) closed form  nu_hat^2 ~ min_j(theta_j^2 + 1/theta_{{j+1}}^2):")
    print(f"  pearson  = {pearson(ex, cl):.4f}   spearman = {spearman(ex, cl):.4f}")
    print(f"  predicted/exact ratio: min={ratios[0]:.4f} "
          f"q25={ratios[len(ratios)//4]:.4f} median={ratios[len(ratios)//2]:.4f} "
          f"q75={ratios[3*len(ratios)//4]:.4f} max={ratios[-1]:.4f}")
    print(f"  straddle-only variant: pearson = {pearson(ex, cs):.4f}")
    print("  falsifier C (corr < 0.9): " +
          ("NOT triggered" if pearson(ex, cl) >= 0.9 else "TRIGGERED"))

    # --- D: lower bound nu_hat > sqrt(2/(a*+2)), and the naive sqrt(2/(a*+1))
    withstar = [r for r in rows if r['a_star'] is not None]
    print(f"\n(C4) lower bounds on nu_hat in terms of a*:")
    for off, lab in ((2, 'sqrt(2/(a*+2))  [rigorous]'),
                     (1, 'sqrt(2/(a*+1))  [naive]   ')):
        viol = [r for r in withstar
                if r['nu_hat'] < math.sqrt(2.0 / (r['a_star'] + off)) - 1e-12]
        slack = sorted(r['nu_hat'] / math.sqrt(2.0 / (r['a_star'] + off))
                       for r in withstar)
        line = (f"  {lab}: violations {len(viol):>4}/{len(withstar)}  "
                f"slack median={slack[len(slack)//2]:.3f} "
                f"q90={slack[int(0.9*len(slack))]:.3f} max={slack[-1]:.1f}")
        if viol:
            w = min(viol, key=lambda r: r['nu_hat']
                    / math.sqrt(2.0 / (r['a_star'] + off)))
            line += (f"  worst={w['nu_hat']/math.sqrt(2.0/(w['a_star']+off)):.4f}"
                     f" (a*={w['a_star']})")
        print(line)

    # --- E: locality -- a* vs the scale-free max_a
    astar = [float(r['a_star']) for r in withstar]
    maxa = [float(r['max_a']) for r in withstar]
    nus = [r['nu_hat'] for r in withstar]
    print(f"\n(C5) locality -- correlation with nu_hat, same sample:")
    print(f"  spearman(nu_hat, a*)     = {spearman(nus, astar):+.4f}   "
          f"[local: partial quotient at the scale Q]")
    print(f"  spearman(nu_hat, max_a)  = {spearman(nus, maxa):+.4f}   "
          f"[scale-free: the June invariant]")
    print(f"  spearman(a*, max_a)      = {spearman(astar, maxa):+.4f}")

    # a* bins
    print(f"\n  nu_hat by a* bin:")
    print(f"{'a*':>6} {'count':>6} {'mean nu_hat':>12} {'min nu_hat':>11} "
          f"{'bound sqrt(2/(a*+2))':>21}")
    for lo, hi, lab in [(1, 1, '1'), (2, 2, '2'), (3, 4, '3-4'),
                        (5, 8, '5-8'), (9, 16, '9-16'), (17, 10 ** 9, '>=17')]:
        grp = [r for r in withstar if lo <= r['a_star'] <= hi]
        if not grp:
            continue
        nn = [r['nu_hat'] for r in grp]
        bd = min(math.sqrt(2.0 / (r['a_star'] + 2)) for r in grp)
        print(f"{lab:>6} {len(grp):>6} {sum(nn)/len(nn):>12.4f} "
              f"{min(nn):>11.4f} {bd:>21.4f}")
    return rows

# ---------------------------------------------------------------------------
# Part E -- does a* separate the C1/C2 recovery classes?
# ---------------------------------------------------------------------------

def part_e():
    print("\n" + "=" * 78)
    print("PART E -- a* vs nu_hat vs max_a on the Exp S C1/C2 task")
    print("=" * 78)
    print(f"  K1={K1_FIXED}, m={M_FIXED}, {SEEDS_N} seeds, {N_CURVES} "
          f"fresh 20-bit j=0 CM curves (C2 = recovers on all seeds)")
    t0 = time.time()
    curves = harvest(N_CURVES)
    print(f"  harvested {len(curves)} curves in {time.time()-t0:.1f}s")

    t1 = time.time()
    for c in curves:
        n = c['n']
        k2 = math.isqrt(n) + 1
        wins = 0
        for s in range(SEEDS_N):
            d = random.Random(s + 7777).randint(1, n - 1)
            wins += run_experiment(c['p'], n, c['lam'], c['G'], M_FIXED, d,
                                   K1_FIXED, k2, s)
        c['wins'] = wins
        c['C2'] = (wins == SEEDS_N)
        cf = nu_cf(n, c['lam'], K1_FIXED, k2)
        c['nu_hat'] = cf['nu_hat']
        c['a_star'] = float(cf['a_star'] if cf['a_star'] is not None else 0)
        c['closed'] = cf['closed_min']
        c['max_a'] = float(max_partial_quotient(c['lam'], n))
        c['mu'] = mu_of(c['lam'], n)
    print(f"  {len(curves)*SEEDS_N} LLL trials in {time.time()-t1:.1f}s")

    lab = [c['C2'] for c in curves]
    base = max(sum(lab), len(lab) - sum(lab)) / len(lab)
    print(f"\n  C1 (fails): {len(lab)-sum(lab)}/{len(lab)}   "
          f"C2 (6/6): {sum(lab)}/{len(lab)}   baseline = {base*100:.1f}%")

    print(f"\n{'invariant':>10} {'C1 range':>18} {'C2 range':>18} "
          f"{'AUC':>7} {'best acc':>9}")
    for key in ('nu_hat', 'closed', 'a_star', 'max_a', 'mu'):
        sc = [c[key] for c in curves]
        s1 = [c[key] for c in curves if not c['C2']]
        s2 = [c[key] for c in curves if c['C2']]
        a = auc(sc, lab)
        acc, _ = best_cut_accuracy(sc, lab)
        r1 = f"[{min(s1):.3f},{max(s1):.3f}]" if s1 else "-"
        r2 = f"[{min(s2):.3f},{max(s2):.3f}]" if s2 else "-"
        print(f"{key:>10} {r1:>18} {r2:>18} {max(a,1-a):>7.3f} {acc*100:>8.1f}%")

    print(f"\n  C2 rate by a*:")
    print(f"{'a*':>6} {'count':>6} {'C2 rate':>9} {'mean wins/6':>12} "
          f"{'mean nu_hat':>12}")
    for lo, hi, l in [(1, 1, '1'), (2, 2, '2'), (3, 4, '3-4'), (5, 8, '5-8'),
                      (9, 10 ** 9, '>=9')]:
        grp = [c for c in curves if lo <= c['a_star'] <= hi]
        if not grp:
            continue
        print(f"{l:>6} {len(grp):>6} {sum(c['C2'] for c in grp)/len(grp):>9.2f} "
              f"{sum(c['wins'] for c in grp)/len(grp):>12.2f} "
              f"{sum(c['nu_hat'] for c in grp)/len(grp):>12.3f}")

    c2 = [c for c in curves if c['C2']]
    if c2:
        print(f"\n  every C2 curve: nu_hat <= {max(c['nu_hat'] for c in c2):.3f}, "
              f"a* >= {int(min(c['a_star'] for c in c2))}")
    return curves

# ---------------------------------------------------------------------------
# Part F -- secp256k1 in closed form
# ---------------------------------------------------------------------------

def part_f():
    print("\n" + "=" * 78)
    print("PART F -- secp256k1 (n, lam): which convergent, which partial quotient")
    print("=" * 78)
    assert (SECP_LAM * SECP_LAM + SECP_LAM + 1) % SECP_N == 0
    pq, convs = cf_convergents(SECP_LAM, SECP_N)
    print(f"  CF of lam/n has {len(pq)} partial quotients; "
          f"max_a = {max(pq)} (at index {pq.index(max(pq))})")
    print(f"  geometric-mean partial quotient = "
          f"{math.exp(sum(math.log(a) for a in pq if a > 0)/sum(1 for a in pq if a>0)):.3f} "
          f"(Khinchin 2.685)")
    print(f"\n{'eff':>6} {'K1':>22} {'s':>4} {'a*':>6} {'theta_s':>9} "
          f"{'theta_s1':>9} {'nu_hat':>8} {'closed':>8} {'bound':>8}")
    for eff in (0.02, 0.05, 0.0993, 0.10, 0.25, 0.50):
        k1 = max(2, int(eff * math.isqrt(SECP_N)))
        k2 = math.isqrt(SECP_N) + 1
        d = nu_cf(SECP_N, SECP_LAM, k1, k2)
        bd = math.sqrt(2.0 / (d['a_star'] + 1)) if d['a_star'] else float('nan')
        print(f"{eff:>6.4f} {k1:>22} {d['s_idx']:>4} {d['a_star']:>6} "
              f"{d['theta_s']:>9.4f} {d['theta_s1']:>9.4f} {d['nu_hat']:>8.4f} "
              f"{d['closed_min']:>8.4f} {bd:>8.4f}")
    print("\n  Caveat (unchanged from 2026-07-29): this thread is conditional on a")
    print("  NON-STANDARD nonce generator k = k1 + lam*k2 with k1 bounded.  It is")
    print("  not an attack on secp256k1 with correctly generated nonces.")

# ---------------------------------------------------------------------------

def main():
    print("=" * 78)
    print("GLV-HNP Thread 23 -- closed form for nu_hat via continued fractions")
    print("=" * 78)
    t0 = time.time()
    part_a()
    part_bcd()
    part_e()
    part_f()
    print(f"\nTotal {time.time()-t0:.1f}s.  Done.")
    return 0


if __name__ == '__main__':
    sys.exit(main())
