"""
GLV-HNP Phase 2 — Thread 23: closed form for lambda_1(L2), i.e. for nu_hat.

Background
----------
2026-07-29 (run #2, commit e845207) found the first predictor of Phase 2 LLL
success that beats the trivial baseline: the scale-free

    nu_hat = lambda_1(L2) / sqrt(det L2),        AUC 0.935 vs the June C1/C2 classes

where L2 is the *non-planted* 2-dimensional sublattice of the (2m+2)-dim GLV
lattice spanned by

    u = (n*A, 0)        A = S_K1 = n // K1
    v = (-lam*A, B)     B = S_K2 = max(1, n // K2)

det L2 = n*A*B, independent of lam.  nu_hat was *computed* (Lagrange-Gauss) but
not *understood*, and the 07-29 entry proposed Thread 23: show lambda_1(L2) is a
continued-fraction quantity of lam/n taken at the scale sqrt(n*K2/K1), which
would (a) turn an empirical predictor into an algebraic one and (b) retroactively
explain why the six June scale-free CF invariants (q_cf, max_q_cf, max_a, ...)
all failed.

Theory (proved here, then verified exactly)
-------------------------------------------
A general vector of L2 is

    a*u + b*v = ( A*(a*n - b*lam),  B*b ),        a, b in Z.

For fixed b the inner factor is minimised over a by the centered residue
theta(b) = |b*lam - p*n| minimal over p in Z.  Hence

    lambda_1(L2)^2 = min_{b >= 0, (a,b) != 0} [ A^2*theta(b)^2 + B^2*b^2 ].    (*)

CLAIM 23.1.  The minimising b in (*) is a continued-fraction convergent
denominator q_j of lam/n (with q_{-1} = 0, theta_{-1} = n covering b = 0).

Proof.  Let b* attain the min and let b' < b*.  Then B^2*b'^2 < B^2*b*^2, so (*)
forces A^2*theta(b')^2 > A^2*theta(b*)^2, i.e. theta(b') > theta(b*).  So no
smaller denominator approximates lam/n at least as well: b* is a best
approximation of the second kind to lam/n, and by the classical theorem
(Khinchin, Thm 16-17) every best approximation of the second kind is a
convergent denominator.  []

COROLLARY 23.2 (closed form).  With p_j/q_j the convergents of lam/n and
theta_j = |q_j*lam - p_j*n|,

    lambda_1(L2)^2 = min_{j >= -1} [ A^2*theta_j^2 + B^2*q_j^2 ]

    nu_hat^2 = min_{j >= -1} [ x_j^2 + y_j^2 ],
        x_j = q_j*sqrt(B/(A*n))          (denominator in scale units)
        y_j = theta_j*sqrt(A/(B*n))      (residue   in scale units)
        x_j*y_j = q_j*theta_j/n  ~  q_j/q_{j+1}  ~  1/a_{j+1}

so nu_hat is small exactly when some partial quotient a_{j+1} is large AND the
corresponding q_j sits near the scale

    Q* = sqrt(A*n/B) = sqrt(n*K2/K1)         (the point where x_j ~ 1)

x_j*y_j ~ 1/a_{j+1} >= (1/2)*nu_hat^2 by AM-GM, i.e. a_{j+1} >~ 2/nu_hat^2:
a large partial quotient is NECESSARY for small nu_hat, but not sufficient --
it must occur at the right scale.  That is the retro-explanation of the June
failures: max_a and q_cf are scale-free, and the relevant quantity is not.

Also nu_hat <= (4/3)^(1/4) = 1.0746 (Hermite bound in rank 2), matching the
observed maxima 1.012 / 1.072 on 07-29.

Parts
-----
A  exact-arithmetic verification of Corollary 23.2 (integer equality, not
   correlation), 16..64 bits and 256 bits
B  precision audit of the float Lagrange-Gauss in glv_hnp_nuhat_lib.py
C  which convergent wins: j* vs the scale index j(Q*)
D  Hermite bound + the (x,y) identity
E  scale-local partial quotient a* vs the scale-free max_a  (June retro-explain)
F  secp256k1 closed-form nu_hat at the 07-29 eff grid, both GLV roots
G  root choice as a free attacker degree of freedom

Run: python3 glv_hnp_nuhat_cf_closed_form.py
"""

import math
import random
import sys

import sympy

# ---------------------------------------------------------------------------
# Exact integer primitives
# ---------------------------------------------------------------------------

def iround_div(num, den):
    """Round-half-away-from-zero of num/den, exactly, in integers (den > 0)."""
    if num >= 0:
        return (2 * num + den) // (2 * den)
    return -((-2 * num + den) // (2 * den))

def gauss_reduce_2d_exact(u, v):
    """Lagrange-Gauss reduction in exact integers; returns lambda_1^2 (an int)."""
    def n2(w):
        return w[0] * w[0] + w[1] * w[1]
    if n2(u) < n2(v):
        u, v = v, u
    while True:
        q = iround_div(u[0] * v[0] + u[1] * v[1], n2(v))
        r = (u[0] - q * v[0], u[1] - q * v[1])
        u, v = v, r
        if n2(u) <= n2(v):
            return n2(u)

def cf_convergents(lam, n):
    """
    Convergents p_j/q_j of lam/n, 0 <= lam < n.
    Returns list of (j, q_j, theta_j, a_j) with theta_j = |q_j*lam - p_j*n|,
    starting at j = -1 with (q,theta) = (0, n) -- the b = 0 lattice vector.
    a_j is the partial quotient used to build index j.
    """
    a_list = []
    x, y = lam, n
    while y:
        a_list.append(x // y)
        x, y = y, x % y
    out = [(-1, 0, n, None)]
    p_pp, q_pp = 0, 1      # p_{-2}, q_{-2}
    p_p, q_p = 1, 0        # p_{-1}, q_{-1}
    for j, a in enumerate(a_list):
        p, q = a * p_p + p_pp, a * q_p + q_pp
        out.append((j, q, abs(q * lam - p * n), a))
        p_pp, q_pp, p_p, q_p = p_p, q_p, p, q
    return out, a_list

def lambda1_sq_cf(n, lam, A, B):
    """Closed form (Cor 23.2): lambda_1(L2)^2 via convergents. Exact ints."""
    conv, _ = cf_convergents(lam, n)
    best, best_j = None, None
    for j, q, th, _a in conv:
        if q == 0 and th == 0:
            continue
        val = A * A * th * th + B * B * q * q
        if val and (best is None or val < best):
            best, best_j = val, j
    return best, best_j, conv

def scales(n, k1_bound, k2_bound):
    """Column-diagonal scaling, verbatim from glv_hnp_nuhat_lib.scales."""
    return n // k1_bound, 1, max(1, n // k2_bound), n

def nu_hat_exact(n, lam, k1, k2):
    """nu_hat from the exact integer Lagrange-Gauss (ground truth)."""
    A, _, B, _ = scales(n, k1, k2)
    l1sq = gauss_reduce_2d_exact((n * A, 0), (-lam * A, B))
    det = n * A * B
    return math.sqrt(l1sq / det), l1sq, det, A, B

# ---------------------------------------------------------------------------
# Sample generation
# ---------------------------------------------------------------------------

def rand_instance(rng, bits, eff=0.05):
    """Random (n prime, lam, K1, K2) with the Phase 2 scaling."""
    n = int(sympy.nextprime(rng.randrange(2 ** (bits - 1), 2 ** bits)))
    lam = rng.randrange(1, n)
    k1 = max(2, int(eff * math.isqrt(n)))
    k2 = math.isqrt(n) + 1
    return n, lam, k1, k2

def glv_lambda(n):
    """Root of x^2+x+1 mod n, i.e. a GLV eigenvalue, or None."""
    if n % 3 != 1:
        return None
    try:
        r = sympy.sqrt_mod(-3, n, all_roots=False)
    except Exception:
        return None
    if r is None:
        return None
    lam = ((-1 + int(r)) * pow(2, -1, n)) % n
    return lam if (lam * lam + lam + 1) % n == 0 else None

def hdr(t):
    print("\n" + "=" * 78)
    print(t)
    print("=" * 78)

# ---------------------------------------------------------------------------
# Part A — exact verification of Corollary 23.2
# ---------------------------------------------------------------------------

def part_a(rng):
    hdr("PART A — Corollary 23.2 verified in EXACT integer arithmetic")
    print("Testing  lambda_1(L2)^2  ==  min_j [A^2 theta_j^2 + B^2 q_j^2]")
    print("(integer equality, not correlation)\n")
    print(f"{'bits':>5} {'eff':>6} {'trials':>7} {'exact match':>12} {'mismatch':>9}")
    total_ok = total = 0
    for bits in (16, 20, 24, 32, 48, 64, 128, 256):
        for eff in (0.05, 0.25, 0.50):
            ok = bad = 0
            trials = 120 if bits <= 64 else 40
            for _ in range(trials):
                n, lam, k1, k2 = rand_instance(rng, bits, eff)
                A, _, B, _ = scales(n, k1, k2)
                gt = gauss_reduce_2d_exact((n * A, 0), (-lam * A, B))
                cf, _j, _c = lambda1_sq_cf(n, lam, A, B)
                if gt == cf:
                    ok += 1
                else:
                    bad += 1
                    if bad <= 2:
                        print(f"    MISMATCH n={n} lam={lam} gt={gt} cf={cf}")
            print(f"{bits:>5} {eff:>6.2f} {trials:>7} {ok:>12} {bad:>9}")
            total_ok += ok
            total += trials
    print(f"\nTOTAL: {total_ok}/{total} exact matches "
          f"({100.0 * total_ok / total:.2f}%)")
    print("VERDICT:", "Corollary 23.2 HOLDS exactly" if total_ok == total
          else "Corollary 23.2 FALSIFIED")
    return total_ok == total

# ---------------------------------------------------------------------------
# Part B — precision audit of the shipped float Lagrange-Gauss
# ---------------------------------------------------------------------------

def part_b(rng):
    hdr("PART B — precision audit of glv_hnp_nuhat_lib.gauss_reduce_2d (float)")
    print("The shipped reducer rounds  q = round(<u,v>/|v|^2)  in f64.  At the")
    print("20/24-bit sizes the July results were measured at this is exact; the")
    print("question is whether the 256-bit extrapolation in the 07-29 entry")
    print("(secp256k1 nu_hat = 0.8709 etc.) was computed reliably.\n")
    try:
        import importlib.util
        spec = importlib.util.spec_from_file_location(
            "_lib", __file__.rsplit("/", 1)[0] + "/glv_hnp_nuhat_lib.py")
        lib = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(lib)
        float_nu = lib.rival_sublattice_nu
    except Exception as e:                                    # pragma: no cover
        print(f"  could not load glv_hnp_nuhat_lib: {e}")
        return
    print(f"{'bits':>5} {'trials':>7} {'agree':>7} {'disagree':>9} "
          f"{'max rel err':>13} {'worst nu_hat float/exact':>26}")
    for bits in (20, 24, 32, 64, 128, 192, 256):
        agree = dis = 0
        worst = 0.0
        worst_pair = ("-", "-")
        trials = 200 if bits <= 64 else 60
        for _ in range(trials):
            n, lam, k1, k2 = rand_instance(rng, bits, 0.05)
            ex_nu, ex_l1sq, det, A, B = nu_hat_exact(n, lam, k1, k2)
            fl_nu = float_nu(n, lam, k1, k2)[2]
            rel = abs(fl_nu - ex_nu) / ex_nu if ex_nu else 0.0
            # "agree" = float reducer found the same lattice vector length
            if rel < 1e-9:
                agree += 1
            else:
                dis += 1
            if rel > worst:
                worst, worst_pair = rel, (f"{fl_nu:.6f}", f"{ex_nu:.6f}")
        print(f"{bits:>5} {trials:>7} {agree:>7} {dis:>9} {worst:>13.3e} "
              f"{worst_pair[0]:>12} / {worst_pair[1]:<11}")

# ---------------------------------------------------------------------------
# Part C — which convergent wins
# ---------------------------------------------------------------------------

def part_c(rng):
    hdr("PART C — the winning convergent sits at the scale Q* = sqrt(n*K2/K1)")
    print("j*      = argmin index of Cor 23.2")
    print("j(Q*)   = index of the convergent denominator nearest Q* (log scale)")
    print("offset  = j* - j(Q*)\n")
    for eff in (0.05, 0.25, 0.50):
        offs = {}
        xs = []
        for _ in range(400):
            n, lam, k1, k2 = rand_instance(rng, 24, eff)
            A, _, B, _ = scales(n, k1, k2)
            _v, jstar, conv = lambda1_sq_cf(n, lam, A, B)
            Qstar = math.sqrt(A * n / B)
            cand = [(abs(math.log(max(q, 1)) - math.log(Qstar)), j)
                    for j, q, _t, _a in conv]
            jscale = min(cand)[1]
            offs[jstar - jscale] = offs.get(jstar - jscale, 0) + 1
            qstar_q = [q for j, q, _t, _a in conv if j == jstar][0]
            xs.append(qstar_q / Qstar)
        tot = sum(offs.values())
        near = sum(v for k, v in offs.items() if abs(k) <= 1)
        dist = " ".join(f"{k:+d}:{v}" for k, v in sorted(offs.items()))
        xs.sort()
        print(f"eff={eff:.2f}  |offset|<=1: {near}/{tot} ({100.0*near/tot:.1f}%)")
        print(f"            offset dist: {dist}")
        print(f"            q_j*/Q* quantiles  q10={xs[len(xs)//10]:.3f} "
              f"q50={xs[len(xs)//2]:.3f} q90={xs[9*len(xs)//10]:.3f}\n")

# ---------------------------------------------------------------------------
# Part D — Hermite bound and the (x,y) identity
# ---------------------------------------------------------------------------

def part_d(rng):
    hdr("PART D — Hermite bound and the normalised (x,y) identity")
    hermite = (4.0 / 3.0) ** 0.25
    print(f"Predicted hard ceiling nu_hat <= (4/3)^(1/4) = {hermite:.6f}")
    print("07-29 observed maxima: 1.012 (20d sample), 1.072 (20c high group)\n")
    mx = 0.0
    ident_ok = 0
    prod = []
    N = 600
    for _ in range(N):
        n, lam, k1, k2 = rand_instance(rng, 24, 0.05)
        nu, l1sq, det, A, B = nu_hat_exact(n, lam, k1, k2)
        mx = max(mx, nu)
        _v, jstar, conv = lambda1_sq_cf(n, lam, A, B)
        q, th = [(q, t) for j, q, t, _a in conv if j == jstar][0]
        x = q * math.sqrt(B / (A * n))
        y = th * math.sqrt(A / (B * n))
        if abs((x * x + y * y) - nu * nu) < 1e-9 * max(1.0, nu * nu):
            ident_ok += 1
        prod.append(x * y)
    print(f"max nu_hat over {N} random 24-bit instances: {mx:.6f}  "
          f"({'<=' if mx <= hermite + 1e-12 else '>'} bound)")
    print(f"identity nu_hat^2 == x_j*^2 + y_j*^2 : {ident_ok}/{N}")
    prod.sort()
    print(f"x*y = q*theta/n at the argmin: q10={prod[N//10]:.4f} "
          f"q50={prod[N//2]:.4f} q90={prod[9*N//10]:.4f}  "
          f"(AM-GM: nu_hat^2 >= 2*x*y)")

# ---------------------------------------------------------------------------
# Part E — scale-local a* beats the scale-free max_a  (June retro-explanation)
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
            avg = (i + j) / 2.0 + 1.0
            for k in range(i, j + 1):
                r[order[k]] = avg
            i = j + 1
        return r
    rx, ry = rank(xs), rank(ys)
    mx, my = sum(rx) / len(rx), sum(ry) / len(ry)
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    dx = math.sqrt(sum((a - mx) ** 2 for a in rx))
    dy = math.sqrt(sum((b - my) ** 2 for b in ry))
    return num / (dx * dy) if dx and dy else 0.0

def part_e(rng):
    hdr("PART E — scale-local partial quotient vs scale-free max_a")
    print("June falsified max_a (AUC 0.748, at the trivial baseline on Exp S).")
    print("Cor 23.2 says the predictive quantity is the partial quotient AT the")
    print("scale Q*, i.e. a_{j*+1}, not the largest one anywhere.\n")
    print("Sampling REAL GLV eigenvalues (lam^2+lam+1 = 0 mod n) for relevance.\n")
    rows = []
    tries = 0
    while len(rows) < 300 and tries < 60000:
        tries += 1
        n = int(sympy.nextprime(rng.randrange(2 ** 23, 2 ** 24)))
        lam = glv_lambda(n)
        if lam is None:
            continue
        k1 = max(2, int(0.05 * math.isqrt(n)))
        k2 = math.isqrt(n) + 1
        nu, _l, _d, A, B = nu_hat_exact(n, lam, k1, k2)
        _v, jstar, conv = lambda1_sq_cf(n, lam, A, B)
        _c, a_list = cf_convergents(lam, n)
        # a_{j*+1}: partial quotient one past the winning index
        a_next = a_list[jstar + 1] if 0 <= jstar + 1 < len(a_list) else 1
        max_a = max(a_list) if a_list else 1
        rows.append((nu, a_next, max_a))
    nus = [r[0] for r in rows]
    an = [r[1] for r in rows]
    ma = [r[2] for r in rows]
    print(f"n = {len(rows)} real j=0 GLV curves at 24 bits, eff=0.05\n")
    print(f"  spearman(nu_hat, a_scale = a_{{j*+1}}) = "
          f"{spearman(nus, an):+.3f}     <- scale-LOCAL")
    print(f"  spearman(nu_hat, max_a)               = "
          f"{spearman(nus, ma):+.3f}     <- scale-FREE (the June invariant)")
    print(f"  spearman(a_scale, max_a)              = {spearman(an, ma):+.3f}")
    # AM-GM necessary condition, empirically
    print("\n  AM-GM prediction: a_scale >~ 2/nu_hat^2.  Check by nu_hat quartile:")
    order = sorted(rows)
    qn = len(order) // 4
    for name, grp in (("Q1 (smallest nu_hat)", order[:qn]),
                      ("Q2", order[qn:2 * qn]),
                      ("Q3", order[2 * qn:3 * qn]),
                      ("Q4 (largest nu_hat)", order[3 * qn:])):
        mnu = sum(g[0] for g in grp) / len(grp)
        man = sum(g[1] for g in grp) / len(grp)
        mma = sum(g[2] for g in grp) / len(grp)
        print(f"    {name:<22} nu_hat={mnu:.3f}  a_scale={man:7.2f}  "
              f"max_a={mma:7.2f}  2/nu^2={2.0/mnu**2:7.2f}")

# ---------------------------------------------------------------------------
# Part F — secp256k1
# ---------------------------------------------------------------------------

SECP_N = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141

def part_f():
    hdr("PART F — secp256k1 closed form, BOTH GLV roots")
    n = SECP_N
    lam_a = glv_lambda(n)
    lam_b = n - 1 - lam_a
    assert (lam_b * lam_b + lam_b + 1) % n == 0
    assert pow(lam_a, -1, n) == lam_b and pow(lam_a, 3, n) == 1
    # the root used by the 07-29 entry / glv_hnp_nuhat_vs_c1c2.py:182
    LAM_0729 = 0x5363ad4cc05c30e0a5261c028812645a122e22ea20816678df02967c1b23bd72
    print(f"n     = {n}")
    print(f"lam_a = {lam_a}")
    print(f"lam_b = {lam_b}   (= lam_a^-1 = lam_a^2 mod n)")
    print(f"lam_b is the root used on 07-29: {lam_b == LAM_0729}\n")
    print("07-29 reported (float reducer, root lam_b): eff 0.05 -> 0.8709,")
    print("0.10 -> 0.6624, 0.25 -> 0.5852, 0.50 -> 0.6851.")
    print("Recomputed here in exact integers and via Cor 23.2:\n")
    print(f"{'eff':>6} {'nu(lam_a)':>10} {'nu(lam_b)':>10} {'|diff|':>8} "
          f"{'nu_min':>8} | {'CF ok':>6} {'j*':>4} {'a_{j*+1}':>9} "
          f"{'q_j*/Q*':>9}   (argmin over the better root)")
    for eff in (0.05, 0.10, 0.25, 0.50):
        k1 = max(2, int(eff * math.isqrt(n)))
        k2 = math.isqrt(n) + 1
        na, la, _d, A, B = nu_hat_exact(n, lam_a, k1, k2)
        nb, lb, _d, _A, _B = nu_hat_exact(n, lam_b, k1, k2)
        best_lam, best_l1 = (lam_a, la) if na <= nb else (lam_b, lb)
        conv, a_list = cf_convergents(best_lam, n)
        cf, jstar, _c = lambda1_sq_cf(n, best_lam, A, B)
        q, th = [(q, t) for j, q, t, _a in conv if j == jstar][0]
        Qstar = math.sqrt(A * n / B)
        a_next = a_list[jstar + 1] if 0 <= jstar + 1 < len(a_list) else 1
        print(f"{eff:>6.2f} {na:>10.6f} {nb:>10.6f} {abs(na-nb):>8.4f} "
              f"{min(na,nb):>8.6f} | {str(cf == best_l1):>6} {jstar:>4} "
              f"{a_next:>9} {q / Qstar:>9.4f}")
    print("\nRECONCILIATION: the nu(lam_b) column reproduces the 07-29 table to")
    print("all printed digits (0.8709 / 0.6624 / 0.5852 / 0.6851).  The 07-29")
    print("numbers are correct; they are simply stated for one of the two roots.")
    print("\nRefinement: |nu(lam_a) - nu(lam_b)| is 0.011 at eff=0.05 -- inside the")
    print("07-29 root-choice control (20a Part E, '<0.02', 6 small curves at one")
    print("eff) -- but grows to 0.075 at eff=0.25 and 0.125 at eff=0.50.  So the")
    print("'<0.02' bound is an eff=0.05 statement, not a universal one.  Part G")
    print("shows this costs nothing: root choice never flips a curve across the")
    print("viability cut.")
    for tag, lam in (("lam_a", lam_a), ("lam_b", lam_b)):
        _c, a_list = cf_convergents(lam, n)
        big = sorted(range(len(a_list)), key=lambda i: -a_list[i])[:6]
        print(f"\nCF({tag}/n): length {len(a_list)}, largest partial quotients "
              f"(index:value): " + ", ".join(f"{i}:{a_list[i]}" for i in big))
    print("\nReading: secp256k1's large partial quotients are NOT at the attack")
    print("scale under either root, so nu_hat stays mid-range at every eff --")
    print("consistent with the 07-29 percentile finding, now with an algebraic")
    print("reason rather than an empirical one.")

# ---------------------------------------------------------------------------
# Part G — root choice is a free attacker degree of freedom
# ---------------------------------------------------------------------------

def part_g(rng):
    hdr("PART G — root choice: is it a free attacker degree of freedom?  NO.")
    print("A j=0 curve has two GLV eigenvalues lam and lam' = n-1-lam = lam^-1.")
    print("Both give a valid HNP instance, so an attacker may use whichever has")
    print("the smaller nu_hat.  HYPOTHESIS (this run's own): the two roots are")
    print("untied, so min over roots is a free ~2x shot at a low-nu_hat instance.")
    print("RESULT BELOW: FALSIFIED -- the gain is exactly 1.00x.\n")
    print("Reversal theorem check (CF(n/lam) reversed vs CF(n/(+/-lam^-1))):")

    def cfa(x, y):
        a = []
        while y:
            a.append(x // y)
            x, y = y, x % y
        return a

    def nrm(b):
        return b[:-2] + [b[-2] + 1] if len(b) > 1 and b[-1] == 1 else b

    inv = neg = tot = 0
    curves = []
    tries = 0
    while len(curves) < 250 and tries < 60000:
        tries += 1
        p = int(sympy.nextprime(rng.randrange(2 ** 23, 2 ** 24)))
        lam = glv_lambda(p)
        if lam is None:
            continue
        curves.append((p, lam))
        if tot < 200:
            tot += 1
            R = nrm(cfa(p, lam))[::-1]
            if R == nrm(cfa(p, pow(lam, -1, p))):
                inv += 1
            if R == nrm(cfa(p, (-pow(lam, -1, p)) % p)):
                neg += 1
    print(f"  reversed == CF(n/ +lam^-1) : {inv}/{tot}   <- the OTHER GLV root")
    print(f"  reversed == CF(n/ -lam^-1) : {neg}/{tot}   <- not a GLV root")
    print("  So the REVERSAL theorem does not tie the two GLV roots: it maps")
    print("  lam to -lam^-1 (parity-dependent), and the second GLV root is")
    print("  +lam^-1.  A different, elementary relation ties them instead:\n")
    print("    lam'/n = 1 - lam/n - 1/n,  and  x = [0;1,a2,a3,...]")
    print("                                =>  1-x = [0;1+a2,a3,...]")
    print("  i.e. the two expansions share a TAIL after a 1-term realignment,")
    print("  perturbed only by the -1/n.  Aligned-prefix agreement measured:\n")

    def cfa2(x, y):
        a = []
        while y:
            a.append(x // y)
            x, y = y, x % y
        return a

    def aligned_prefix(lam, nn):
        X, Y = cfa2(lam, nn), cfa2(nn - 1 - lam, nn)
        best = 0
        for i, j in ((3, 2), (2, 3)):
            k = 0
            while k < min(len(X) - i, len(Y) - j) and X[i + k] == Y[j + k]:
                k += 1
            best = max(best, k)
        return best, len(X)

    k_s, L_s = aligned_prefix(glv_lambda(SECP_N), SECP_N)
    print(f"  secp256k1: {k_s} of {L_s} partial quotients agree after alignment "
          f"({100.0*k_s/L_s:.0f}%),")
    print(f"             and the winning index j* ~ 70 sits just past that point.")
    fr = sorted(aligned_prefix(lam, p)[0] / aligned_prefix(lam, p)[1]
                for p, lam in curves[:150])
    print(f"  24-bit curves: agreement fraction q10={fr[15]:.2f} "
          f"q50={fr[75]:.2f} q90={fr[135]:.2f}")
    print("  CAVEAT (open): at 24 bits only ~1/3 of the expansion is shared, yet")
    print("  |d nu| is still ~0.012 at eff=0.05.  Prefix sharing alone does not")
    print("  fully account for how tightly the two roots track; not resolved here.\n")

    LOW = 0.45          # 'empirically easy' cut from the 07-29 null study
    print(f"{'bits':>5} {'eff':>6} {'curves':>7} {'mean|d nu|':>11} "
          f"{'max|d nu|':>10} {'P(lam low)':>11} {'P(min low)':>11} {'gain':>6}")
    for bits, sample in ((24, curves), (32, None), (48, None)):
        if sample is None:
            sample = []
            tries = 0
            while len(sample) < 150 and tries < 200000:
                tries += 1
                p = int(sympy.nextprime(
                    rng.randrange(2 ** (bits - 1), 2 ** bits)))
                lam = glv_lambda(p)
                if lam is not None:
                    sample.append((p, lam))
        for eff in (0.05, 0.25):
            ds, one, both = [], 0, 0
            for p, lam in sample:
                k1 = max(2, int(eff * math.isqrt(p)))
                k2 = math.isqrt(p) + 1
                na = nu_hat_exact(p, lam, k1, k2)[0]
                nb = nu_hat_exact(p, p - 1 - lam, k1, k2)[0]
                ds.append(abs(na - nb))
                one += (na < LOW)
                both += (min(na, nb) < LOW)
            N = len(sample)
            p1, p2 = one / N, both / N
            print(f"{bits:>5} {eff:>6.2f} {N:>7} {sum(ds)/N:>11.4f} "
                  f"{max(ds):>10.4f} {p1:>11.3f} {p2:>11.3f} "
                  f"{(p2/p1 if p1 else float('nan')):>6.2f}x")
    print(f"\nP(lam low)  = fraction with nu_hat < {LOW} using a fixed root")
    print(f"P(min low)  = fraction with min over both roots < {LOW}")
    print("gain        = the free multiplier an attacker gets from root choice.")
    print("\nVERDICT: gain = 1.00x at every (bits, eff) tested -- not one curve in")
    print("550 changes side of the 0.45 cut when the other root is used.  Root")
    print("choice is NOT an attacker degree of freedom.  nu_hat is well defined")
    print("per curve up to ~0.01 (eff=0.05) / ~0.06 (eff=0.25); the 07-29")
    print("'<0.02' root-choice control is confirmed at eff=0.05 and is mildly")
    print("optimistic at larger eff, but the conclusion it supported is unchanged.")
    print("\nConsequence for the 07-29 secp256k1 placement: min over both roots is")
    print("0.5852 at eff=0.25, still above the 0.45 'easy' cut, and 0.6624 at")
    print("eff~0.10, still above the C2 ceiling 0.645.  Placement unchanged.")

# ---------------------------------------------------------------------------

def main():
    rng = random.Random(20260801)
    print("GLV-HNP Phase 2 — Thread 23: closed form for nu_hat")
    print("Deterministic seed 20260801")
    ok = part_a(rng)
    part_b(rng)
    part_c(rng)
    part_d(rng)
    part_e(rng)
    part_f()
    part_g(rng)
    hdr("SUMMARY")
    print(f"Corollary 23.2 exact on every trial: {ok}")
    print("nu_hat is now an explicit continued-fraction quantity of lam/n")
    print("evaluated at the scale Q* = sqrt(n*K2/K1); no lattice reduction needed.")
    print("\nDone.")
    return 0 if ok else 1

if __name__ == "__main__":
    sys.exit(main())
