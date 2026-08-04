"""
GLV-HNP -- Thread 23: closed form for lambda_1(L2), i.e. for nu_hat.

Background (2026-07-29 log entry, sections 3-5)
----------------------------------------------
The (2m+2)-dimensional column-scaled Phase 2 lattice contains m independent
copies of a *non-planted* 2-dimensional sublattice

    L2 = < (n*S_K1, 0), (-lam*S_K1, S_K2) >,      det L2 = n*S_K1*S_K2,

and the scale-free

    nu_hat = lambda_1(L2) / sqrt(det L2)

was measured to separate the June C1/C2 classes at AUC 0.935 while eight
algebraic curve invariants sat at the trivial baseline.  nu_hat was *computed*
(one Lagrange-Gauss reduction) but not *understood*.

What this script establishes
----------------------------
A. lambda_1(L2) has an exact closed form in the continued fraction of lam/n:

       lambda_1(L2)^2 = min_j [ (r_j*S_K1)^2 + (q_j*S_K2)^2 ]

   where p_j/q_j runs over the convergents of lam/n (including the seed
   convergent 1/0, which contributes the vector (n*S_K1, 0)) and
   r_j = |q_j*lam - p_j*n|.  The claim is that the minimiser b of
   |b*lam mod+- n| against b is always a convergent DENOMINATOR -- a best
   approximation of the second kind -- so the minimisation over the whole
   lattice collapses to a scan over ~log(n) convergents.  Verified exactly
   (bit-for-bit, integer arithmetic) below, not statistically.

B. The relevant convergent index j* is the one at the SCALE
   q_j ~ sqrt(n*S_K1/S_K2) = sqrt(n*K2/K1), confirming the 2026-07-29
   conjecture.

C. Writing W = S_K1/S_K2 and, per convergent,

       theta_j = r_j*q_j/n   (in (0,1], ~ 1/a_{j+1})
       t_j     = r_j*W/q_j   (the balance ratio)

   the closed form becomes

       nu_hat^2 = min_j theta_j * (t_j + 1/t_j)

   whose floor 2*theta_{j*} is attained when the scale balances.  Hence
   nu_hat is small -- attack easy -- exactly when lam/n has a LARGE partial
   quotient a_{j*+1} at the scale sqrt(n*K2/K1), and

       nu_hat >= sqrt(2*min_j theta_j) ~ sqrt(2/max_j a_{j+1}).

D. This retro-explains the June falsifications.  max_a = max_j a_j is
   scale-free: it reports the largest partial quotient anywhere in the CF,
   while nu_hat is governed by a_{j+1} at ONE index j* selected by K1/K2.
   Changing K1 moves j* and therefore changes nu_hat for a FIXED curve --
   so no scale-free invariant of (n, lam) can predict it, which is precisely
   what six June sweeps found the hard way.

Falsifier (from the 2026-07-29 next-step proposal)
--------------------------------------------------
"if predicted-from-CF nu_hat correlates <0.9 with computed nu_hat, the
convergent-scale story is wrong."  Test A is much stronger than a 0.9
correlation: it asserts exact equality of two independently computed integers.

Run: python3 glv_hnp_nuhat_cf.py
No third-party dependencies (sympy optional, only for the GLV curve harvest).
"""

import math
import random

# ---------------------------------------------------------------------------
# Scaling, verbatim from glv_hnp_phase2_lambda_threshold.py (e845207) so the
# nu_hat computed here is the same number as the 2026-07-29 run's.
# ---------------------------------------------------------------------------

def scales(n, k1_bound, k2_bound):
    """Column-diagonal scaling validated 2026-07-26."""
    return (n // k1_bound, 1, max(1, n // k2_bound), n)


def gauss_reduce_2d_exact(u, v):
    """
    Lagrange-Gauss reduction over the integers; returns lambda_1^2 as an int.

    The e845207 version (`glv_hnp_nuhat_lib.rival_sublattice_nu`) rounds via
    float division.  Measured head to head against this one, 200 draws per size
    at eff=0.10: it agrees to f64 epsilon (max rel err 3.1e-16) at 20, 24, 64,
    128, 192, 256 and 384 bits, and raises OverflowError at 521 bits, where
    nrm2 ~ 2^1080 exceeds the f64 range.  So the published 20/24-bit and
    secp256k1 nu_hat numbers stand unchanged; only P-521-sized inputs need this
    exact version.  Use this one anyway -- it is not measurably slower and it
    removes the size caveat entirely.
    """
    def nrm2(w):
        return w[0] * w[0] + w[1] * w[1]

    if nrm2(u) < nrm2(v):
        u, v = v, u
    while True:
        nv = nrm2(v)
        if nv == 0:
            return nrm2(u)
        d = u[0] * v[0] + u[1] * v[1]
        # round-half-to-even-free nearest integer of d/nv, exactly:
        q = (2 * d + nv) // (2 * nv)
        r = (u[0] - q * v[0], u[1] - q * v[1])
        u, v = v, r
        if nrm2(u) <= nrm2(v):
            return nrm2(u)


def nu_hat_gauss(n, lam, k1_bound, k2_bound):
    """nu_hat via Lagrange-Gauss (the 2026-07-29 definition)."""
    s_k1, _, s_k2, _ = scales(n, k1_bound, k2_bound)
    l1sq = gauss_reduce_2d_exact((n * s_k1, 0), (-lam * s_k1, s_k2))
    det = n * s_k1 * s_k2
    return math.sqrt(l1sq / det), l1sq, det


# ---------------------------------------------------------------------------
# A. the continued-fraction closed form
# ---------------------------------------------------------------------------

def convergents(a, b):
    """
    Convergents p_j/q_j of a/b, seeded with 1/0 so that the b=0 lattice vector
    (n*S_K1, 0) appears as the j=-1 term.  Returns (q_j, r_j, a_{j+1}) with
    r_j = |q_j*a - p_j*b|; a_{j+1} is the next partial quotient (0 at the end).
    """
    quotients = []
    x, y = a, b
    while y:
        qt, rem = divmod(x, y)
        quotients.append(qt)
        x, y = y, rem

    # h_{-2}=0, h_{-1}=1 ; k_{-2}=1, k_{-1}=0
    h2, h1 = 0, 1
    k2, k1 = 1, 0
    # seed term k_{-1}=0: the lattice vector (n*S_K1, 0), r = |0*a - 1*b| = b
    out = [(0, b)]
    for qt in quotients:
        h = qt * h1 + h2
        k = qt * k1 + k2
        out.append((k, abs(k * a - h * b)))
        h2, h1 = h1, h
        k2, k1 = k1, k
    # attach the partial quotient a_{j+1} that FOLLOWS each convergent.
    # out[0] is the seed (index j = -1), out[i] is convergent j = i-1,
    # whose following quotient is quotients[i] when it exists.
    res = []
    for idx, (q, r) in enumerate(out):
        nxt = quotients[idx] if idx < len(quotients) else 0
        res.append((q, r, nxt))
    return res


def nu_hat_cf(n, lam, k1_bound, k2_bound):
    """
    nu_hat by scanning the convergents of lam/n.  Returns
    (nu_hat, l1sq, best_index, best_q, best_r, next_partial_quotient, table).
    """
    s_k1, _, s_k2, _ = scales(n, k1_bound, k2_bound)
    det = n * s_k1 * s_k2
    table = convergents(lam, n)
    best = None
    for j, (q, r, a_next) in enumerate(table):
        val = (r * s_k1) ** 2 + (q * s_k2) ** 2
        if val == 0:
            continue
        if best is None or val < best[0]:
            best = (val, j, q, r, a_next)
    l1sq, j, q, r, a_next = best
    return math.sqrt(l1sq / det), l1sq, j, q, r, a_next, table


def theta_and_t(n, r, q, s_k1, s_k2):
    """theta = r*q/n  and  t = r*S_K1/(q*S_K2), the two closed-form variables."""
    theta = (r * q) / n
    t = float('inf') if q == 0 else (r * s_k1) / (q * s_k2)
    return theta, t


# ---------------------------------------------------------------------------
# test A -- exactness
# ---------------------------------------------------------------------------

def test_exactness(trials=4000, seed=20260804):
    """
    Gauss-reduced lambda_1^2 vs CF-scanned lambda_1^2, as exact integers,
    across a wide range of moduli, lambdas and (K1,K2) scale ratios.
    """
    rng = random.Random(seed)
    fails, checked = [], 0
    sizes = [16, 20, 24, 32, 48, 64, 96, 128, 192, 256, 384, 521]
    for _ in range(trials):
        bits = rng.choice(sizes)
        n = rng.randrange(2 ** (bits - 1), 2 ** bits) | 1
        lam = rng.randrange(1, n)
        eff = rng.choice([0.01, 0.05, 0.1, 0.25, 0.5, 0.9])
        k1 = max(2, int(eff * math.isqrt(n)))
        k2 = math.isqrt(n) + 1
        _, g_sq, _ = nu_hat_gauss(n, lam, k1, k2)
        _, c_sq, j, q, r, a_next, _ = nu_hat_cf(n, lam, k1, k2)
        checked += 1
        if g_sq != c_sq:
            fails.append((bits, n, lam, k1, k2, g_sq, c_sq, j, q, r))
    return checked, fails


# ---------------------------------------------------------------------------
# test B -- the minimising convergent sits at the scale sqrt(n*K2/K1)
# ---------------------------------------------------------------------------

def test_scale_localisation(trials=1500, seed=77, bits=64):
    rng = random.Random(seed)
    rows = []
    for _ in range(trials):
        n = rng.randrange(2 ** (bits - 1), 2 ** bits) | 1
        lam = rng.randrange(1, n)
        eff = rng.choice([0.02, 0.05, 0.1, 0.25, 0.5])
        k1 = max(2, int(eff * math.isqrt(n)))
        k2 = math.isqrt(n) + 1
        s_k1, _, s_k2, _ = scales(n, k1, k2)
        _, _, j, q, r, a_next, table = nu_hat_cf(n, lam, k1, k2)
        target = math.sqrt(n * s_k1 / s_k2)
        # index of the convergent whose denominator is closest to target in log
        best_j, best_d = None, None
        for jj, (qq, _, _) in enumerate(table):
            if qq == 0:
                continue
            d = abs(math.log(qq) - math.log(target))
            if best_d is None or d < best_d:
                best_j, best_d = jj, d
        rows.append((j, best_j, len(table)))
    exact = sum(1 for j, bj, _ in rows if j == bj)
    within1 = sum(1 for j, bj, _ in rows if abs(j - bj) <= 1)
    within2 = sum(1 for j, bj, _ in rows if abs(j - bj) <= 2)
    return len(rows), exact, within1, within2


# ---------------------------------------------------------------------------
# test C -- nu_hat^2 = min_j theta_j*(t_j + 1/t_j), and the sqrt(2/a) floor
# ---------------------------------------------------------------------------

def test_closed_form(trials=1200, seed=1234, bits=64):
    rng = random.Random(seed)
    max_rel_err = 0.0
    floor_ok = 0
    rows = []
    for _ in range(trials):
        n = rng.randrange(2 ** (bits - 1), 2 ** bits) | 1
        lam = rng.randrange(1, n)
        eff = rng.choice([0.02, 0.05, 0.1, 0.25, 0.5])
        k1 = max(2, int(eff * math.isqrt(n)))
        k2 = math.isqrt(n) + 1
        s_k1, _, s_k2, _ = scales(n, k1, k2)
        nu, l1sq, j, q, r, a_next, table = nu_hat_cf(n, lam, k1, k2)
        # closed form via (theta, t)
        best = None
        for (qq, rr, _) in table:
            if qq == 0 and rr == 0:
                continue
            det = n * s_k1 * s_k2
            if qq == 0:          # seed vector (n*S_K1, 0)
                val = (rr * s_k1) ** 2 / det
            elif rr == 0:        # terminal convergent, vector (0, n*S_K2)
                val = (qq * s_k2) ** 2 / det
            else:
                th, t = theta_and_t(n, rr, qq, s_k1, s_k2)
                val = th * (t + 1.0 / t)
            if best is None or val < best:
                best = val
        pred = math.sqrt(best)
        rel = abs(pred - nu) / nu
        max_rel_err = max(max_rel_err, rel)
        # floor check: nu_hat >= sqrt(2*theta_min) over convergents with q>0
        th_min = min((rr * qq) / n for (qq, rr, _) in table if qq > 0 and rr > 0)
        if nu >= math.sqrt(2 * th_min) - 1e-12:
            floor_ok += 1
        rows.append((nu, a_next, th_min))
    return trials, max_rel_err, floor_ok, rows


# ---------------------------------------------------------------------------
# test D -- why the June scale-free CF invariants failed
# ---------------------------------------------------------------------------

def max_partial_quotient(a, b):
    """max_a, the invariant falsified by Exp S (2026-06-29)."""
    best = 0
    while b:
        q, r = divmod(a, b)
        best = max(best, q)
        a, b = b, r
    return best


def spearman(xs, ys):
    def rank(v):
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
    rx, ry = rank(xs), rank(ys)
    mx, my = sum(rx) / len(rx), sum(ry) / len(ry)
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    dx = math.sqrt(sum((a - mx) ** 2 for a in rx))
    dy = math.sqrt(sum((b - my) ** 2 for b in ry))
    return num / (dx * dy) if dx and dy else 0.0


def test_scale_dependence(trials=400, seed=99, bits=20):
    """
    (i) corr(max_a, nu_hat) vs corr(a_{j*+1}, nu_hat) on one sample;
    (ii) for FIXED (n, lam), sweep K1 and show nu_hat and j* both move --
         the property no scale-free invariant can track.
    """
    rng = random.Random(seed)
    nus, maxas, ajs = [], [], []
    for _ in range(trials):
        n = rng.randrange(2 ** (bits - 1), 2 ** bits) | 1
        lam = rng.randrange(1, n)
        k1 = max(2, int(0.1 * math.isqrt(n)))
        k2 = math.isqrt(n) + 1
        nu, _, j, q, r, a_next, _ = nu_hat_cf(n, lam, k1, k2)
        nus.append(nu)
        maxas.append(max_partial_quotient(lam, n))
        ajs.append(a_next if a_next else 1)
    r_maxa = spearman(maxas, nus)
    r_aj = spearman(ajs, nus)

    # (ii) fixed curve, sweep the scale
    n = rng.randrange(2 ** 63, 2 ** 64) | 1
    lam = rng.randrange(1, n)
    sweep = []
    for eff in (0.01, 0.02, 0.05, 0.1, 0.2, 0.4, 0.8):
        k1 = max(2, int(eff * math.isqrt(n)))
        k2 = math.isqrt(n) + 1
        nu, _, j, q, r, a_next, _ = nu_hat_cf(n, lam, k1, k2)
        sweep.append((eff, nu, j, a_next))
    return r_maxa, r_aj, (n, lam, sweep), max_partial_quotient(lam, n)


# ---------------------------------------------------------------------------
# secp256k1 placement, in CF terms
# ---------------------------------------------------------------------------

SECP_P = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F
SECP_N = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
SECP_LAM = 0x5363AD4CC05C30E0A5261C028812645A122E22EA20816678DF02967C1B23BD72


def secp_report():
    assert (SECP_LAM * SECP_LAM + SECP_LAM + 1) % SECP_N == 0, "lambda check"
    out = []
    for eff in (0.05, 0.10, 0.25, 0.50):
        k1 = max(2, int(eff * math.isqrt(SECP_N)))
        k2 = math.isqrt(SECP_N) + 1
        s_k1, _, s_k2, _ = scales(SECP_N, k1, k2)
        nu_g, g_sq, _ = nu_hat_gauss(SECP_N, SECP_LAM, k1, k2)
        nu_c, c_sq, j, q, r, a_next, table = nu_hat_cf(SECP_N, SECP_LAM, k1, k2)
        th, t = theta_and_t(SECP_N, r, q, s_k1, s_k2)
        out.append((eff, nu_g, nu_c, g_sq == c_sq, j, len(table),
                    q.bit_length(), a_next, th, t))
    return out, max_partial_quotient(SECP_LAM, SECP_N)


# ---------------------------------------------------------------------------

def reachability(n, lam, k2_bound=None, effs=None):
    """
    Which convergents can the attack scale ever select?

    With the Phase 2 convention K2 = isqrt(n)+1 and K1 = eff*isqrt(n), the
    selecting scale is sqrt(n*S_K1/S_K2) ~ sqrt(n/eff), which for eff in (0,1]
    is >= sqrt(n).  So only convergents with q >~ sqrt(n) -- roughly the top
    half of the CF -- are ever reachable, and every partial quotient below that
    scale is invisible to nu_hat no matter how K1 is chosen.  That is a
    sharper form of the "max_a is scale-free" objection: max_a is often
    attained at an unreachable index.
    """
    if k2_bound is None:
        k2_bound = math.isqrt(n) + 1
    if effs is None:
        effs = [10 ** (-e / 4.0) for e in range(0, 25)]   # 1.0 down to ~1e-6
    table = convergents(lam, n)
    # index and value of the largest partial quotient, and its convergent size
    a_list = [(idx, a) for idx, (q, r, a) in enumerate(table)]
    idx_max, a_max = max(a_list, key=lambda z: z[1])
    q_at_max = table[idx_max][0]
    rows = []
    for eff in effs:
        k1 = max(2, int(eff * math.isqrt(n)))
        nu, _, j, q, r, a_next, _ = nu_hat_cf(n, lam, k1, k2_bound)
        rows.append((eff, nu, j, q.bit_length(), a_next))
    reachable_j = sorted({r[2] for r in rows})
    return rows, idx_max, a_max, q_at_max.bit_length(), len(table), reachable_j


def main():
    print("=" * 78)
    print("GLV-HNP Thread 23 -- closed form for nu_hat via continued fractions")
    print("=" * 78)

    print("\n[A] EXACTNESS: Lagrange-Gauss lambda_1^2  ==  CF-scan lambda_1^2")
    print("    (exact integer comparison; sizes 16..521 bits, eff 0.01..0.9)")
    checked, fails = test_exactness()
    print(f"    checked = {checked}   mismatches = {len(fails)}")
    if fails:
        for f in fails[:5]:
            print(f"      MISMATCH bits={f[0]} n={f[1]} lam={f[2]} "
                  f"gauss={f[5]} cf={f[6]}")
        print("    => VERDICT: the convergent-denominator claim is FALSE.")
    else:
        print("    => VERDICT: lambda_1(L2) is attained at a convergent")
        print("       denominator of lam/n in every case.  Closed form CONFIRMED.")

    print("\n[B] SCALE LOCALISATION: is the minimising j the convergent")
    print("    nearest q ~ sqrt(n*S_K1/S_K2) = sqrt(n*K2/K1)?  (64-bit)")
    tot, exact, w1, w2 = test_scale_localisation()
    print(f"    trials={tot}  j* == nearest-in-log: {exact} ({100*exact/tot:.1f}%)"
          f"  |dj|<=1: {w1} ({100*w1/tot:.1f}%)  |dj|<=2: {w2} ({100*w2/tot:.1f}%)")

    print("\n[C] CLOSED FORM  nu_hat^2 = min_j theta_j*(t_j + 1/t_j),")
    print("    theta_j = r_j*q_j/n,  t_j = r_j*S_K1/(q_j*S_K2)")
    tr, mre, floor_ok, rows = test_closed_form()
    print(f"    trials={tr}  max relative error vs Gauss nu_hat = {mre:.3e}")
    print(f"    floor  nu_hat >= sqrt(2*min_j theta_j)  holds: {floor_ok}/{tr}")
    # decile of nu_hat against a_{j*+1}
    rows_sorted = sorted(rows, key=lambda z: z[0])
    print("    nu_hat vs the partial quotient a_{j*+1} that follows the")
    print("    minimising convergent (quintiles):")
    q5 = len(rows_sorted) // 5
    print(f"      {'quintile':>9} {'nu_hat mid':>11} {'mean a_(j*+1)':>14}"
          f" {'median a':>9}")
    for i in range(5):
        chunk = rows_sorted[i * q5:(i + 1) * q5]
        if not chunk:
            continue
        mid = chunk[len(chunk) // 2][0]
        aa = sorted(x[1] if x[1] else 1 for x in chunk)
        print(f"      {i+1:>9} {mid:>11.3f} {sum(aa)/len(aa):>14.2f}"
              f" {aa[len(aa)//2]:>9}")

    print("\n[D] WHY THE JUNE SCALE-FREE CF INVARIANTS FAILED")
    r_maxa, r_aj, (n0, lam0, sweep), maxa0 = test_scale_dependence()
    print(f"    20-bit sample, eff=0.10, 400 draws:")
    print(f"      spearman(max_a,        nu_hat) = {r_maxa:+.3f}   <- June's invariant")
    print(f"      spearman(a_(j*+1),     nu_hat) = {r_aj:+.3f}   <- scale-selected one")
    print(f"    FIXED 64-bit (n, lam), max_a = {maxa0} is constant; sweep the scale:")
    print(f"      {'eff':>6} {'nu_hat':>8} {'j*':>4} {'a_(j*+1)':>9}")
    for eff, nu, j, a_next in sweep:
        print(f"      {eff:>6.2f} {nu:>8.4f} {j:>4} {a_next if a_next else 0:>9}")
    nu_min = min(s[1] for s in sweep)
    nu_max = max(s[1] for s in sweep)
    print(f"    one curve, nu_hat ranges over [{nu_min:.4f}, {nu_max:.4f}] as K1 moves.")
    print("    => a scale-free function of (n, lam) cannot predict nu_hat.")

    print("\n[E] secp256k1 in CF terms")
    rep, maxa = secp_report()
    print(f"    max_a over the whole CF of lam/n = {maxa}")
    print(f"    {'eff':>5} {'nu(Gauss)':>10} {'nu(CF)':>9} {'exact':>6} {'j*':>4}"
          f" {'/len':>5} {'q* bits':>8} {'a_(j*+1)':>9} {'theta':>7} {'t':>7}")
    for eff, ng, nc, eq, j, ln, qb, a_next, th, t in rep:
        print(f"    {eff:>5.2f} {ng:>10.4f} {nc:>9.4f} {str(eq):>6} {j:>4}"
              f" {ln:>5} {qb:>8} {a_next if a_next else 0:>9} {th:>7.4f} {t:>7.3f}")

    print("\n[F] REACHABILITY: which convergents can any choice of K1 select?")
    rows, idx_max, a_max, qbits_max, ln, reach = reachability(SECP_N, SECP_LAM)
    print(f"    secp256k1: CF of lam/n has {ln} convergents; max_a = {a_max}"
          f" at index {idx_max} (q is {qbits_max} bits)")
    print(f"    n is {SECP_N.bit_length()} bits, so sqrt(n) ~ {SECP_N.bit_length()//2} bits")
    print(f"      {'eff':>10} {'nu_hat':>8} {'j*':>4} {'q* bits':>8} {'a_(j*+1)':>9}")
    for eff, nu, j, qb, a_next in rows[::3]:
        print(f"      {eff:>10.2e} {nu:>8.4f} {j:>4} {qb:>8} "
              f"{a_next if a_next else 0:>9}")
    print(f"    reachable j over eff in [1e-6, 1]: {min(reach)}..{max(reach)}"
          f"  of 0..{ln-1}")
    if not (min(reach) <= idx_max <= max(reach)):
        print(f"    => the max_a index {idx_max} is UNREACHABLE at every eff:"
              f" secp256k1's largest")
        print(f"       partial quotient ({a_max}) sits below the sqrt(n) scale"
              f" floor and can never")
        print(f"       be selected.  max_a is not merely scale-free, it is"
              f" measuring an invisible")
        print(f"       part of the CF.")
    else:
        print(f"    => the max_a index {idx_max} IS reachable (at some eff).")
    nu_lo = min(r[1] for r in rows)
    print(f"    min nu_hat over the whole reachable scale range = {nu_lo:.4f}"
          f"  (easy tail is < 0.45)")

    print("\n" + "=" * 78)


if __name__ == "__main__":
    main()
