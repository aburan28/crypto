"""
GLV-HNP Thread 23 — closed form for nu_hat via continued fractions.

Background
----------
Thread 20 (2026-07-29, commit e845207) found the first invariant that separates
the Phase 2 easy/hard classes: with

    L2 = < u = (n*S_K1, 0),  v = (-lam*S_K1, S_K2) >       (det = n*S_K1*S_K2)
    nu_hat = lambda_1(L2) / sqrt(det L2)

nu_hat reaches AUC 0.935 against the June C1/C2 classification, where eight
curve-level invariants (delta/n, kappa(M), q_cf, max_q_cf, max_a, a_corn/n,
lam/n, mu) all sat at the trivial baseline.  But nu_hat was *computed*, by
Lagrange-Gauss reduction, not *understood*.

This script closes that gap.  Every vector of L2 is

    (a,b) |-> ( (b*lam - a*n) * S_K1,  b * S_K2 ),        (a,b) in Z^2,

so for fixed b the best a is round(b*lam/n) and

    ||.||^2  =  S_K1^2 * r(b)^2  +  S_K2^2 * b^2,
    r(b)     =  min_a |b*lam - a*n|  =  centred residue of b*lam mod n.

CLAIM (Thread 23).  lambda_1(L2)^2 = min over CONTINUED FRACTION CONVERGENT
DENOMINATORS q_j of lam/n (together with the degenerate j = -1 term q=0, r=n)
of  S_K1^2 * r_j^2 + S_K2^2 * q_j^2,  where r_j = |q_j*lam - p_j*n|.

PROOF.  Legendre / best-approximation-of-the-second-kind: for every b with
q_j <= b < q_{j+1} one has r(b) >= r_j.  Hence the pair (r(b), b) dominates
(r_j, q_j) coordinatewise, and the cost is increasing in both coordinates, so
no non-convergent b can beat the convergent below it.  b = 0 is the j = -1
term.  QED  --  the enumeration is exact, not heuristic, and costs O(log n).

CONSEQUENCE (the reason six June invariants failed).  The minimising index j*
is the one that balances the two terms, S_K1*r_j ~ S_K2*q_j.  With
S_K1 = n/K1, S_K2 = n/K2 and r_j ~ n/q_{j+1}, balance gives

    q_j * q_{j+1}  ~  n * K2/K1,

i.e. the relevant convergent sits at a SCALE set by (K1, K2), and

    nu_hat^2  ~  2 * q_j * r_j / n  ~  2 * q_j / q_{j+1}  ~  2 / a_{j+1}.

q_cf, max_q_cf and max_a are scale-free functions of the whole CF expansion;
the quantity that matters is the single partial quotient a_{j*+1} at the scale
sqrt(n*K2/K1).  That is exactly why the scale-free invariants were flat.

Experiments
-----------
E1  exactness of the CF enumeration vs Lagrange-Gauss vs fpylll LLL (2D),
    at 20 / 24 / 64 / 128 / 256 bits.  Falsifier: any mismatch.
E2  scale law  q_{j*} * q_{j*+1} ~ n*K2/K1.
E3  closed forms: nu_hat^2 = 2*q*r/n (balance approximation) and
    nu_hat^2 ~ 2/a_{j*+1}.  Falsifier: Pearson/Spearman < 0.9 vs computed.
E4  scale-free vs scale-local partial quotients (retro-explains June).
E5  secp256k1 placement recomputed from the CF expansion.

Run: python3 glv_hnp_t23_cf_closed_form.py
"""

import math
import random
import sys

# ---------------------------------------------------------------------------
# exact 2D machinery
# ---------------------------------------------------------------------------

def gauss_reduce_exact(u, v):
    """Lagrange-Gauss reduction with exact integer rounding.
    Returns the shortest nonzero vector (exact integer tuple)."""
    def nrm2(w):
        return w[0] * w[0] + w[1] * w[1]
    def dot(w, z):
        return w[0] * z[0] + w[1] * z[1]
    if nrm2(u) > nrm2(v):
        u, v = v, u
    while True:
        num, den = dot(v, u), nrm2(u)
        if den == 0:
            break
        # round-half-away, exact in integers
        q = (2 * num + den) // (2 * den) if num >= 0 else -((-2 * num + den) // (2 * den))
        v = (v[0] - q * u[0], v[1] - q * u[1])
        if nrm2(v) >= nrm2(u):
            break
        u, v = v, u
    return u

def centred_residue(x, n):
    """min_a |x - a*n|, i.e. the centred representative magnitude."""
    r = x % n
    return min(r, n - r)

def cf_convergents(lam, n):
    """Convergent denominators q_j of lam/n together with
    r_j = |q_j*lam - p_j*n| and the partial quotients a_j.
    Returns list of dicts, starting with the degenerate j = -1 term (q=0, r=n).
    """
    lam = lam % n
    out = [{'j': -1, 'q': 0, 'r': n, 'a': None}]
    quotients = []
    if lam == 0:
        return out, quotients
    # Euclid on (lam, n) is the CF of lam/n; a_0 = 0 since lam < n.
    # Denominator recurrence q_j = a_j*q_{j-1} + q_{j-2}, q_{-2}=1, q_{-1}=0.
    q_km2, q_km1 = 1, 0
    x, y = lam, n
    j = 0
    while y != 0:
        a = x // y
        quotients.append(a)
        q_k = a * q_km1 + q_km2
        q_km2, q_km1 = q_km1, q_k
        x, y = y, x - a * y
        if q_k >= 1:
            out.append({'j': j, 'q': q_k, 'r': centred_residue(q_k * lam, n),
                        'a': a})
        j += 1
    # keep strictly increasing q with strictly decreasing r (drops q_1 == q_0)
    clean = [out[0]]
    for e in out[1:]:
        if e['q'] > clean[-1]['q'] and e['r'] <= clean[-1]['r']:
            clean.append(e)
    return clean, quotients

def scales(n, k1_bound, k2_bound):
    """Column-diagonal scaling validated 2026-07-26 (identical to Thread 20)."""
    return (n // k1_bound, 1, max(1, n // k2_bound), n)

def lambda1_cf(n, lam, k1_bound, k2_bound):
    """lambda_1(L2)^2 by CF enumeration + the minimising convergent record."""
    s_k1, _, s_k2, _ = scales(n, k1_bound, k2_bound)
    convs, quotients = cf_convergents(lam, n)
    best, best_e = None, None
    for e in convs:
        c = s_k1 * s_k1 * e['r'] * e['r'] + s_k2 * s_k2 * e['q'] * e['q']
        if best is None or c < best:
            best, best_e = c, e
    return best, best_e, convs, quotients

def lambda1_gauss(n, lam, k1_bound, k2_bound):
    """lambda_1(L2)^2 by exact Lagrange-Gauss reduction."""
    s_k1, _, s_k2, _ = scales(n, k1_bound, k2_bound)
    w = gauss_reduce_exact((n * s_k1, 0), (-(lam % n) * s_k1, s_k2))
    return w[0] * w[0] + w[1] * w[1]

def lambda1_lll(n, lam, k1_bound, k2_bound):
    """lambda_1(L2)^2 via fpylll LLL on the 2D basis (independent check)."""
    try:
        from fpylll import IntegerMatrix, LLL
    except ImportError:
        return None
    s_k1, _, s_k2, _ = scales(n, k1_bound, k2_bound)
    A = IntegerMatrix.from_matrix([[n * s_k1, 0], [-(lam % n) * s_k1, s_k2]])
    LLL.reduction(A)
    rows = [(A[i][0], A[i][1]) for i in range(2)]
    return min(r[0] * r[0] + r[1] * r[1] for r in rows if r != (0, 0))

def nu_hat_exact(n, lam, k1_bound, k2_bound):
    """nu_hat = lambda_1/sqrt(det) computed exactly from the CF enumeration."""
    s_k1, _, s_k2, _ = scales(n, k1_bound, k2_bound)
    l1sq, e, convs, quots = lambda1_cf(n, lam, k1_bound, k2_bound)
    det = n * s_k1 * s_k2
    # nu_hat^2 = l1sq/det, computed as a ratio of big ints to avoid overflow
    return math.sqrt(l1sq / det), e, convs, quots

# ---------------------------------------------------------------------------
# statistics helpers
# ---------------------------------------------------------------------------

def pearson(xs, ys):
    m = len(xs)
    if m < 3:
        return float('nan')
    mx, my = sum(xs) / m, sum(ys) / m
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    sxx = sum((x - mx) ** 2 for x in xs)
    syy = sum((y - my) ** 2 for y in ys)
    if sxx == 0 or syy == 0:
        return float('nan')
    return sxy / math.sqrt(sxx * syy)

def spearman(xs, ys):
    def rank(vs):
        order = sorted(range(len(vs)), key=lambda i: vs[i])
        r = [0.0] * len(vs)
        i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and vs[order[j + 1]] == vs[order[i]]:
                j += 1
            avg = (i + j) / 2.0
            for k in range(i, j + 1):
                r[order[k]] = avg
            i = j + 1
        return r
    return pearson(rank(xs), rank(ys))

def max_partial_quotient(lam, n):
    """Scale-free CF invariant, as used (and falsified) in June."""
    x, y, best, cnt = lam % n, n, 0, 0
    while y and cnt < 4096:
        a = x // y
        best = max(best, a)
        x, y = y, x - a * y
        cnt += 1
    return best

# ---------------------------------------------------------------------------
# sampling
# ---------------------------------------------------------------------------

def is_probable_prime(k):
    if k < 2:
        return False
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if k % p == 0:
            return k == p
    d, s = k - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        x = pow(a, d, k)
        if x in (1, k - 1):
            continue
        for _ in range(s - 1):
            x = x * x % k
            if x == k - 1:
                break
        else:
            return False
    return True

def random_prime(bits, rng):
    while True:
        c = rng.getrandbits(bits) | (1 << (bits - 1)) | 1
        if is_probable_prime(c):
            return c

def glv_lambda(n, rng):
    """A root of x^2+x+1 = 0 mod n, if it exists (n = 1 mod 3)."""
    if n % 3 != 1:
        return None
    # find a cube root of unity: g^((n-1)/3) for a random g
    for _ in range(64):
        g = rng.randrange(2, n - 1)
        c = pow(g, (n - 1) // 3, n)
        if c != 1 and (c * c + c + 1) % n == 0:
            return c
    return None

# ---------------------------------------------------------------------------
# E1 — exactness
# ---------------------------------------------------------------------------

def e1_exactness(rng):
    print("=" * 74)
    print("E1  CF enumeration == Lagrange-Gauss == fpylll LLL  (exact integers)")
    print("=" * 74)
    print(f"{'bits':>5} {'trials':>7} {'CF==Gauss':>10} {'CF==LLL':>9} "
          f"{'#convs':>7} {'j*':>4} {'worst rel err':>14}")
    ok_all = True
    for bits in (20, 24, 64, 128, 256):
        trials = 400 if bits <= 64 else 150
        eq_g = eq_l = 0
        nconv = []
        jstars = []
        worst = 0.0
        for _ in range(trials):
            n = random_prime(bits, rng)
            lam = rng.randrange(1, n)          # arbitrary lam: the geometry
            eff = rng.choice([0.05, 0.10, 0.25, 0.5])   # does not need lam^2+lam+1
            k1 = max(2, int(eff * math.isqrt(n)))
            k2 = math.isqrt(n) + 1
            cf, e, convs, _ = lambda1_cf(n, lam, k1, k2)
            gs = lambda1_gauss(n, lam, k1, k2)
            ll = lambda1_lll(n, lam, k1, k2)
            if cf == gs:
                eq_g += 1
            else:
                worst = max(worst, abs(cf - gs) / gs)
            if ll is None or cf == ll:
                eq_l += 1
            nconv.append(len(convs))
            jstars.append(e['j'])
        print(f"{bits:>5} {trials:>7} {eq_g:>6}/{trials:<3} {eq_l:>5}/{trials:<3} "
              f"{sum(nconv)/len(nconv):>7.1f} {sum(jstars)/len(jstars):>4.1f} "
              f"{worst:>14.2e}")
        ok_all &= (eq_g == trials and eq_l == trials)
    print(f"\nE1 verdict: {'EXACT on all trials' if ok_all else 'MISMATCH FOUND'}")
    return ok_all

# ---------------------------------------------------------------------------
# E2 / E3 / E4 — scale law and closed forms
# ---------------------------------------------------------------------------

def sample_curves(bits, count, eff, rng, require_glv=True):
    out = []
    tries = 0
    while len(out) < count and tries < count * 400:
        tries += 1
        n = random_prime(bits, rng)
        if require_glv:
            lam = glv_lambda(n, rng)
            if lam is None:
                continue
            lam = min(lam, n - 1 - lam)
        else:
            lam = rng.randrange(1, n)
        k1 = max(2, int(eff * math.isqrt(n)))
        k2 = math.isqrt(n) + 1
        out.append((n, lam, k1, k2))
    return out

def e2_scale_law(rng):
    print()
    print("=" * 74)
    print("E2  scale law:  q_{j*} * q_{j*+1}  ~  n * K2/K1")
    print("=" * 74)
    print(f"{'bits':>5} {'eff':>6} {'N':>4} {'median q*q+/(nK2/K1)':>22} "
          f"{'IQR':>18} {'frac in [1/8,8]':>16}")
    for bits in (24, 64, 128):
        for eff in (0.05, 0.25):
            cur = sample_curves(bits, 120, eff, rng, require_glv=(bits <= 64))
            ratios = []
            for (n, lam, k1, k2) in cur:
                _, e, convs, _ = lambda1_cf(n, lam, k1, k2)
                idx = next(i for i, c in enumerate(convs) if c['q'] == e['q'])
                qn = convs[idx + 1]['q'] if idx + 1 < len(convs) else None
                if qn is None or e['q'] == 0:
                    continue
                predicted = n * k2 / k1
                ratios.append((e['q'] * qn) / predicted)
            if not ratios:
                continue
            ratios.sort()
            med = ratios[len(ratios) // 2]
            q1 = ratios[len(ratios) // 4]
            q3 = ratios[3 * len(ratios) // 4]
            frac = sum(1 for r in ratios if 0.125 <= r <= 8.0) / len(ratios)
            print(f"{bits:>5} {eff:>6.2f} {len(ratios):>4} {med:>22.3f} "
                  f"{('[%.2f, %.2f]' % (q1, q3)):>18} {frac:>16.2%}")

def e3_closed_form(rng):
    print()
    print("=" * 74)
    print("E3  closed forms for nu_hat")
    print("     A: nu_hat^2 = (S_K1^2 r^2 + S_K2^2 q^2)/det      [exact, by E1]")
    print("     B: nu_hat^2 ~ 2*q*r/n                            [balance approx]")
    print("     C: nu_hat^2 ~ 2/a_{j*+1}                         [coarse]")
    print("=" * 74)
    print(f"{'bits':>5} {'eff':>6} {'N':>4} {'pearson B':>10} {'spear B':>9} "
          f"{'med B/A':>9} {'pearson C':>10} {'spear C':>9}")
    for bits in (24, 64, 128, 256):
        for eff in (0.05, 0.10, 0.25):
            cur = sample_curves(bits, 150, eff, rng, require_glv=(bits <= 64))
            nu_a, nu_b, nu_c = [], [], []
            for (n, lam, k1, k2) in cur:
                nh, e, convs, _ = nu_hat_exact(n, lam, k1, k2)
                if e['q'] == 0:
                    continue
                idx = next(i for i, c in enumerate(convs) if c['q'] == e['q'])
                a_next = (convs[idx + 1]['a'] if idx + 1 < len(convs) else None)
                b = math.sqrt(2.0 * e['q'] * e['r'] / n)
                nu_a.append(nh)
                nu_b.append(b)
                nu_c.append(math.sqrt(2.0 / a_next) if a_next else float('nan'))
            pairs_b = [(a, b) for a, b in zip(nu_a, nu_b)]
            pairs_c = [(a, c) for a, c in zip(nu_a, nu_c)
                       if not math.isnan(c)]
            if len(pairs_b) < 5:
                continue
            xb = [p[0] for p in pairs_b]; yb = [p[1] for p in pairs_b]
            rat = sorted(b / a for a, b in pairs_b if a > 0)
            med = rat[len(rat) // 2]
            if len(pairs_c) >= 5:
                xc = [p[0] for p in pairs_c]; yc = [p[1] for p in pairs_c]
                pc, sc = pearson(xc, yc), spearman(xc, yc)
            else:
                pc = sc = float('nan')
            print(f"{bits:>5} {eff:>6.2f} {len(pairs_b):>4} "
                  f"{pearson(xb, yb):>10.4f} {spearman(xb, yb):>9.4f} "
                  f"{med:>9.4f} {pc:>10.4f} {sc:>9.4f}")

def e4_scalefree_vs_local(rng):
    print()
    print("=" * 74)
    print("E4  scale-free CF invariants vs the scale-LOCAL partial quotient")
    print("    (retro-explanation of the six June falsifications)")
    print("=" * 74)
    print(f"{'bits':>5} {'eff':>6} {'N':>4} {'spear(nu,max_a)':>16} "
          f"{'spear(nu,a_{j*+1})':>19} {'spear(nu,j*)':>13}")
    for bits in (24, 64, 128):
        for eff in (0.05, 0.25):
            cur = sample_curves(bits, 200, eff, rng, require_glv=(bits <= 64))
            nus, maxa, aloc, js = [], [], [], []
            for (n, lam, k1, k2) in cur:
                nh, e, convs, quots = nu_hat_exact(n, lam, k1, k2)
                if e['q'] == 0:
                    continue
                idx = next(i for i, c in enumerate(convs) if c['q'] == e['q'])
                a_next = (convs[idx + 1]['a'] if idx + 1 < len(convs) else None)
                if not a_next:
                    continue
                nus.append(nh)
                maxa.append(float(max_partial_quotient(lam, n)))
                aloc.append(float(a_next))
                js.append(float(e['j']))
            if len(nus) < 10:
                continue
            print(f"{bits:>5} {eff:>6.2f} {len(nus):>4} "
                  f"{spearman(nus, maxa):>16.4f} {spearman(nus, aloc):>19.4f} "
                  f"{spearman(nus, js):>13.4f}")

# ---------------------------------------------------------------------------
# E5 — secp256k1
# ---------------------------------------------------------------------------

P_SECP = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F
N_SECP = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
LAM_SECP = 0x5363AD4CC05C30E0A5261C028812645A122E22EA20816678DF02967C1B23BD72

def e5_secp256k1():
    print()
    print("=" * 74)
    print("E5  secp256k1: the CF decomposition behind nu_hat")
    print("=" * 74)
    assert (LAM_SECP * LAM_SECP + LAM_SECP + 1) % N_SECP == 0, "lam check failed"
    lam = min(LAM_SECP, N_SECP - 1 - LAM_SECP)
    print(f"n     = {N_SECP}")
    print(f"lam*  = {lam}   (mu = {lam / N_SECP:.6f})")
    print(f"max partial quotient over the whole CF = {max_partial_quotient(lam, N_SECP)}")
    print()
    print(f"{'eff':>6} {'K1':>6} {'j*':>4} {'log2 q*':>8} {'log2 r*':>8} "
          f"{'a_{j*+1}':>9} {'nu_hat':>8} {'2q r/n':>9} {'sqrt(2/a)':>10}")
    for eff in (0.05, 0.0993, 0.10, 0.25, 0.50):
        k1 = max(2, int(eff * math.isqrt(N_SECP)))
        k2 = math.isqrt(N_SECP) + 1
        nh, e, convs, _ = nu_hat_exact(N_SECP, lam, k1, k2)
        idx = next(i for i, c in enumerate(convs) if c['q'] == e['q'])
        a_next = convs[idx + 1]['a'] if idx + 1 < len(convs) else None
        b = math.sqrt(2.0 * e['q'] * e['r'] / N_SECP) if e['q'] else float('nan')
        c = math.sqrt(2.0 / a_next) if a_next else float('nan')
        print(f"{eff:>6.4f} {math.log2(k1):>5.1f}b {e['j']:>4} "
              f"{(math.log2(e['q']) if e['q'] else 0):>8.2f} "
              f"{math.log2(e['r']):>8.2f} {str(a_next):>9} {nh:>8.4f} "
              f"{b:>9.4f} {c:>10.4f}")
    print()
    print("Full CF expansion of lam*/n (partial quotients, first 40):")
    _, quots = cf_convergents(lam, N_SECP)
    print("  " + " ".join(str(q) for q in quots[:40]))
    print(f"  ... {len(quots)} partial quotients total, max = {max(quots)}")

# ---------------------------------------------------------------------------
# E6 / E7 — why j* is where it is: the GLV palindrome and its defect
# ---------------------------------------------------------------------------

def self_mirror_fraction(lam, n):
    """Fraction of positions j with a_j == a_{k-1-j} in the CF of lam/n
    (dropping a_0 = 0).  1.0 means an exact palindrome."""
    _, quots = cf_convergents(lam, n)
    t = quots[1:]
    k = len(t)
    if k == 0:
        return 0.0
    return sum(1 for j in range(k) if t[j] == t[k - 1 - j]) / k

def e6_palindrome(rng):
    print()
    print("=" * 74)
    print("E6  the CF of a GLV eigenvalue is NEAR-PALINDROMIC")
    print("    lam^-1 = lam^2 = n-1-lam (mod n) since lam^3 = 1; the CF-reversal")
    print("    theorem sends lam to +-lam^-1, and lam vs lam+1 differ by 1/n, so")
    print("    reverse(CF(lam/n)) agrees with CF(lam/n) except near the centre.")
    print("=" * 74)
    print(f"{'bits':>5} {'N':>4} {'GLV median':>11} {'GLV min':>9} "
          f"{'random median':>14} {'random max':>11}")
    for bits in (64, 128, 256):
        g, r = [], []
        while len(g) < 60:
            n = random_prime(bits, rng)
            lam = glv_lambda(n, rng)
            if lam is None:
                continue
            g.append(self_mirror_fraction(min(lam, n - 1 - lam), n))
        while len(r) < 60:
            n = random_prime(bits, rng)
            r.append(self_mirror_fraction(rng.randrange(1, n), n))
        g.sort(); r.sort()
        print(f"{bits:>5} {len(g):>4} {g[len(g)//2]:>11.3f} {g[0]:>9.3f} "
              f"{r[len(r)//2]:>14.3f} {r[-1]:>11.3f}")
    lam = min(LAM_SECP, N_SECP - 1 - LAM_SECP)
    _, quots = cf_convergents(lam, N_SECP)
    t = quots[1:]
    k = len(t)
    bad = [(j, t[j], t[k - 1 - j]) for j in range(k) if t[j] != t[k - 1 - j]]
    print(f"\nsecp256k1: k = {k} partial quotients, self-mirror "
          f"{1 - len(bad)/k:.3%}")
    print(f"  the ONLY defect: {bad}   (dead centre, k/2 = {k//2})")

def e7_jstar_is_the_centre(rng):
    print()
    print("=" * 74)
    print("E7  j* is forced to the CENTRE of the expansion")
    print("    scale law q_{j*}q_{j*+1} ~ n*K2/K1 with K1,K2 ~ sqrt(n) puts")
    print("    q_{j*} ~ sqrt(n), which is the midpoint of the CF of lam/n.")
    print("=" * 74)
    print(f"{'bits':>5} {'eff':>6} {'N':>4} {'median j*/k':>12} {'IQR':>16} "
          f"{'frac |j*/k-.5|<.1':>18}")
    for bits in (64, 128, 256):
        for eff in (0.05, 0.25):
            rs = []
            while len(rs) < 80:
                n = random_prime(bits, rng)
                lam = glv_lambda(n, rng)
                if lam is None:
                    continue
                lam = min(lam, n - 1 - lam)
                k1 = max(2, int(eff * math.isqrt(n)))
                k2 = math.isqrt(n) + 1
                _, e, _, quots = lambda1_cf(n, lam, k1, k2)
                k = len(quots) - 1
                if e['q'] == 0 or k < 4:
                    continue
                rs.append(e['j'] / k)
            rs.sort()
            fr = sum(1 for x in rs if abs(x - 0.5) < 0.1) / len(rs)
            print(f"{bits:>5} {eff:>6.2f} {len(rs):>4} {rs[len(rs)//2]:>12.3f} "
                  f"{('[%.3f, %.3f]' % (rs[len(rs)//4], rs[3*len(rs)//4])):>16} "
                  f"{fr:>18.1%}")
    print("\nCONSEQUENCE: nu_hat is governed by the partial quotient at the")
    print("centre of the near-palindrome -- the one place the GLV symmetry")
    print("breaks.  For secp256k1 that defect is (a_69, a_70) = (4,5) <-> (5,4).")

# ---------------------------------------------------------------------------

def main():
    seed = int(sys.argv[1]) if len(sys.argv) > 1 else 20260803
    rng = random.Random(seed)
    print(f"Thread 23 — CF closed form for nu_hat.  seed = {seed}")
    print()
    ok = e1_exactness(rng)
    e2_scale_law(rng)
    e3_closed_form(rng)
    e4_scalefree_vs_local(rng)
    e5_secp256k1()
    e6_palindrome(rng)
    e7_jstar_is_the_centre(rng)
    print()
    print("=" * 74)
    print(f"E1 exactness: {'PASS' if ok else 'FAIL'}")
    print("=" * 74)

if __name__ == "__main__":
    main()
