"""
GLV-HNP Phase 2 -- Thread 23: closed form for lambda_1(L2), and why the
scale-free continued-fraction invariants failed.

Background (Thread 20, 2026-07-29)
----------------------------------
The (2m+2)-dimensional GLV-HNP lattice contains m disjoint copies of a
NON-planted 2-dimensional sublattice

    L2 = < u = (n*S1, 0),  v = (-lam*S1, S2) >,     S1 = n//K1,  S2 = max(1, n//K2)

on each coordinate pair (i, m+1+i).  det(L2) = n*S1*S2 does not depend on lam,
so the scale-free

    nu_hat = lambda_1(L2) / sqrt(det L2)

isolates the lam-dependence of the geometry.  Thread 20 measured nu_hat as a
strong empirical predictor of Phase 2 attack success (AUC 0.935 against the
June C1/C2 classes; spearman(nu_hat, p_hat) = -0.583 on a fixed curve with
120 synthetic lam), but nu_hat was *computed* (one Lagrange-Gauss reduction),
not *understood*.

This script derives and verifies the closed form.

Derivation
----------
A general element of L2 is

    a*v + b*u = ( (b*n - a*lam)*S1,  a*S2 ),      a, b in Z,

with squared norm  S1^2 (a*lam - b*n)^2 + S2^2 a^2.  For a != 0 the optimal b
makes |a*lam - b*n| = E(a) := the centered residue of a*lam mod n.  For a = 0
the shortest nonzero element is u itself, of norm n*S1.  Hence exactly

    lambda_1(L2)^2 = min{ (n*S1)^2,  min_{a>=1} [ S1^2 E(a)^2 + S2^2 a^2 ] }.   (*)

Claim: the minimising a in (*) is a continued-fraction convergent denominator
q_j of lam/n.  Proof: let a* be the minimiser.  If some 1 <= a' < a* had
E(a') < E(a*), then a' would give a strictly smaller value of both terms, so a*
would not be minimal.  Hence E(a') >= E(a*) for every a' < a*, i.e. a* is a
best approximation of the second kind to lam/n, and those are precisely the
convergent denominators.  Therefore

    lambda_1(L2)^2 = min{ (n*S1)^2,  min_j [ S1^2 e_j^2 + S2^2 q_j^2 ] },       (**)
    e_j = |q_j*lam - p_j*n|,

a closed form over the O(log n) convergents of lam/n -- no reduction needed.

Scale-free form.  Put rho = S1/S2 (~ K2/K1) and

    e_hat_j = e_j*sqrt(rho/n),   q_hat_j = q_j/sqrt(rho*n),
    theta_j = e_hat_j * q_hat_j = e_j*q_j/n  in (0,1).

Then nu_hat^2 = min_j (e_hat_j^2 + q_hat_j^2) >= 2*min_j theta_j, with equality
iff e_hat_j = q_hat_j, i.e. iff q_j = rho*e_j.  Since e_j ~ n/q_{j+1}, that
balance point is

    q_j * q_{j+1} ~ rho * n = n*S1/S2 ~ n*K2/K1,      so    q_j ~ sqrt(n*K2/K1).

THIS IS THE POINT.  theta_j ~ q_j/q_{j+1} ~ 1/a_{j+1}, so nu_hat is small --
the attack is easy -- exactly when lam/n has a LARGE PARTIAL QUOTIENT AT THE
SCALE q ~ sqrt(n*K2/K1).  A large partial quotient anywhere else is worthless.
That is why every scale-free CF invariant tried in June (q_cf, max_q_cf,
max_a) sat at the trivial baseline: they are blind to where the large partial
quotient sits, and the relevant location depends on the K1/K2 split, not on
the curve.

Experiments
-----------
A  Exactness of (**) against exact-integer Lagrange-Gauss, 20/24/64/128/256
   bits and the real secp256k1.  Falsifier as posed on 2026-07-29: correlation
   of predicted vs computed nu_hat < 0.9.  Prediction: bit-exact, r = 1.
B  Does the winning convergent sit at the predicted scale sqrt(n*K2/K1)?
C  Local partial quotient a_{j*+1} vs the global max_a, both against nu_hat.
D  Causal replication of Thread 20 Arm 5 (one FIXED curve, synthetic lam):
   spearman of p_hat against a_local, max_a and nu_hat.
E  secp256k1: the explicit convergent that determines nu_hat, per eff.

Run: python3 glv_hnp_nuhat_closed_form.py
"""

import importlib.util
import math
import random
import sys
import time

import sympy

_HERE = __file__.rsplit("/", 1)[0]

def _load(name, fname):
    spec = importlib.util.spec_from_file_location(name, _HERE + "/" + fname)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

_t20a = _load("_t20a", "glv_hnp_phase2_lambda_threshold.py")
_t20b = _load("_t20b", "glv_hnp_phase2_mu_response.py")

eisenstein_decompose = _t20a.eisenstein_decompose
j0_traces = _t20a.j0_traces
glv_eigenvalues = _t20a.glv_eigenvalues
identify_twist = _t20a.identify_twist
find_generator = _t20a.find_generator
mu_of = _t20a.mu_of
scales = _t20a.scales
run_experiment = _t20a.run_experiment
spearman = _t20b.spearman

# ---------------------------------------------------------------------------
# Exact integer primitives (no floats anywhere on the critical path)
# ---------------------------------------------------------------------------

def _n2(w):
    return w[0] * w[0] + w[1] * w[1]

def round_div(a, b):
    """Nearest integer to a/b for b > 0, exact for arbitrary-size ints."""
    return (2 * a + b) // (2 * b)

def gauss_reduce_exact(u, v):
    """Lagrange-Gauss reduction in exact integer arithmetic.

    Returns lambda_1(L)^2 as an exact int.  The float version in
    glv_hnp_phase2_lambda_threshold.py:216 divides via float and takes
    math.sqrt of a quantity that is ~2^1012 for a 256-bit n -- within an order
    of magnitude of the f64 ceiling.  This one has no such ceiling.
    """
    if _n2(u) < _n2(v):
        u, v = v, u
    while True:
        nv = _n2(v)
        if nv == 0:
            return _n2(u)
        q = round_div(u[0] * v[0] + u[1] * v[1], nv)
        r = (u[0] - q * v[0], u[1] - q * v[1])
        u, v = v, r
        if _n2(u) <= _n2(v):
            return _n2(u)

def cf_data(a, b):
    """Continued fraction of a/b.

    Returns (partial_quotients, convergents) where convergents is a list of
    (q_j, e_j) with q_j the denominator of the j-th convergent and
    e_j = |q_j*a - p_j*b|.  The final entry has q = b/gcd(a,b), e = 0.
    """
    quots = []
    convs = []
    x, y = a, b
    h1, h2 = 1, 0      # p_{j-1}, p_{j-2}
    k1, k2 = 0, 1      # q_{j-1}, q_{j-2}
    while y:
        c = x // y
        quots.append(c)
        h = c * h1 + h2
        k = c * k1 + k2
        h2, h1 = h1, h
        k2, k1 = k1, k
        convs.append((k, abs(k * a - h * b)))
        x, y = y, x - c * y
    return quots, convs

def lambda1_sq_closed_form(n, lam, s1, s2):
    """Formula (**).  Returns (lambda_1^2, j_star, quots, convs).

    j_star is the index of the winning convergent, or -1 if the a = 0 vector
    u = (n*S1, 0) is itself shortest.
    """
    quots, convs = cf_data(lam % n, n)
    best = (n * s1) ** 2
    j_star = -1
    for j, (q, e) in enumerate(convs):
        val = s1 * s1 * e * e + s2 * s2 * q * q
        if val < best:
            best, j_star = val, j
    return best, j_star, quots, convs

def brute_lambda1_sq(n, lam, s1, s2, a_max):
    """Reference: exhaustive min over 0 <= a <= a_max.  Small n only."""
    best = (n * s1) ** 2
    arg = 0
    for a in range(1, a_max + 1):
        r = (a * lam) % n
        e = min(r, n - r)
        val = s1 * s1 * e * e + s2 * s2 * a * a
        if val < best:
            best, arg = val, a
    return best, arg

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def k1k2(n, eff):
    return max(2, int(eff * math.isqrt(n))), math.isqrt(n) + 1

def nu_hat_from_sq(l1sq, n, s1, s2):
    """sqrt(lambda_1^2 / (n*s1*s2)) computed without overflowing f64."""
    det = n * s1 * s2
    # integer ratio first, then float -- both l1sq and det can exceed 1e308
    num, den = l1sq, det
    g = math.gcd(num, den)
    num, den = num // g, den // g
    shift = max(num.bit_length(), den.bit_length()) - 900
    if shift > 0:
        num >>= shift
        den >>= shift
    return math.sqrt(num / den)

def isqrt_scale(n, s1, s2):
    """The predicted balance scale q ~ sqrt(n*S1/S2)."""
    return math.isqrt(n * s1 // s2)

def local_a(quots, j_star):
    """Partial quotient immediately after the winning convergent."""
    if j_star < 0 or j_star + 1 >= len(quots):
        return 0
    return quots[j_star + 1]

def analyse(n, lam, k1, k2):
    s1, _, s2, _ = scales(n, k1, k2)
    l1sq_cf, j_star, quots, convs = lambda1_sq_closed_form(n, lam, s1, s2)
    l1sq_lg = gauss_reduce_exact((n * s1, 0), (-lam * s1, s2))
    q_j = convs[j_star][0] if j_star >= 0 else 0
    e_j = convs[j_star][1] if j_star >= 0 else 0
    return {
        'n': n, 'lam': lam, 's1': s1, 's2': s2,
        'l1sq_cf': l1sq_cf, 'l1sq_lg': l1sq_lg, 'exact': l1sq_cf == l1sq_lg,
        'nu_cf': nu_hat_from_sq(l1sq_cf, n, s1, s2),
        'nu_lg': nu_hat_from_sq(l1sq_lg, n, s1, s2),
        'j_star': j_star, 'q_star': q_j, 'e_star': e_j,
        'theta': (e_j * q_j) / n if j_star >= 0 else float('nan'),
        'a_local': local_a(quots, j_star),
        'max_a': max(quots) if quots else 0,
        'n_convs': len(convs), 'quots': quots, 'convs': convs,
        'q_scale': isqrt_scale(n, s1, s2),
    }

# ---------------------------------------------------------------------------
# Experiment A -- exactness of the closed form
# ---------------------------------------------------------------------------

SECP_N = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
SECP_LAM = 0x5363AD4CC05C30E0A5261C028812645A122E22EA20816678DF02967C1B23BD72

def exp_a():
    print("\n" + "=" * 78)
    print("Exp A -- closed form (**) vs exact-integer Lagrange-Gauss")
    print("=" * 78)
    rng = random.Random(20260729)
    rows = []
    for bits in (20, 24, 64, 128, 256):
        n = int(sympy.nextprime(rng.getrandbits(bits) | (1 << (bits - 1))))
        ok = mismatch = 0
        nus_cf, nus_lg = [], []
        for _ in range(400):
            lam = rng.randint(2, n - 2)
            for eff in (0.02, 0.05, 0.1, 0.25, 0.5):
                k1, k2 = k1k2(n, eff)
                r = analyse(n, lam, k1, k2)
                nus_cf.append(r['nu_cf'])
                nus_lg.append(r['nu_lg'])
                if r['exact']:
                    ok += 1
                else:
                    mismatch += 1
                    if mismatch <= 3:
                        print(f"    MISMATCH n={n} lam={lam} eff={eff} "
                              f"cf={r['l1sq_cf']} lg={r['l1sq_lg']}")
        r_sp = spearman(nus_cf, nus_lg)
        rows.append((bits, ok, mismatch, r_sp))
        print(f"  {bits:>4}-bit  n={n}")
        print(f"           exact matches {ok}/{ok+mismatch}   "
              f"spearman(nu_cf, nu_lg) = {r_sp:.6f}")

    print("\n  secp256k1 (real n, real GLV lambda):")
    assert (SECP_LAM * SECP_LAM + SECP_LAM + 1) % SECP_N == 0, "lambda check"
    print("    lambda^2 + lambda + 1 = 0 mod n   [verified]")
    for eff in (0.02, 0.05, 0.1, 0.25, 0.5):
        k1, k2 = k1k2(SECP_N, eff)
        r = analyse(SECP_N, SECP_LAM, k1, k2)
        print(f"    eff={eff:<5} nu_cf={r['nu_cf']:.6f}  nu_lg={r['nu_lg']:.6f}  "
              f"exact={r['exact']}")

    print("\n  Brute-force cross-check (small n, exhaustive a <= 20000):")
    bad = 0
    for _ in range(200):
        n = int(sympy.nextprime(rng.randint(10 ** 4, 10 ** 5)))
        lam = rng.randint(2, n - 2)
        eff = rng.choice((0.02, 0.05, 0.1, 0.25))
        k1, k2 = k1k2(n, eff)
        s1, _, s2, _ = scales(n, k1, k2)
        cf, _, _, _ = lambda1_sq_closed_form(n, lam, s1, s2)
        bf, _ = brute_lambda1_sq(n, lam, s1, s2, min(20000, n))
        if cf != bf:
            bad += 1
            if bad <= 3:
                print(f"    BRUTE MISMATCH n={n} lam={lam} cf={cf} bf={bf}")
    print(f"    {200-bad}/200 agree with exhaustive search")
    return rows

# ---------------------------------------------------------------------------
# Experiment B -- is the winning convergent at the predicted scale?
# ---------------------------------------------------------------------------

def exp_b():
    print("\n" + "=" * 78)
    print("Exp B -- does the winning convergent sit at q ~ sqrt(n*K2/K1)?")
    print("=" * 78)
    rng = random.Random(4242)
    for bits in (24, 64, 256):
        n = int(sympy.nextprime(rng.getrandbits(bits) | (1 << (bits - 1))))
        for eff in (0.05, 0.25):
            k1, k2 = k1k2(n, eff)
            logs, off, tot = [], {}, 0
            for _ in range(500):
                lam = rng.randint(2, n - 2)
                r = analyse(n, lam, k1, k2)
                if r['j_star'] < 0:
                    continue
                tot += 1
                logs.append(math.log(r['q_star'] / r['q_scale'], 2))
                # index offset from the convergent nearest the predicted scale
                qs = [q for q, _ in r['convs']]
                j_sc = min(range(len(qs)),
                           key=lambda j: abs(math.log(qs[j] / r['q_scale'])))
                d = r['j_star'] - j_sc
                off[d] = off.get(d, 0) + 1
            mean = sum(logs) / len(logs)
            var = sum((x - mean) ** 2 for x in logs) / len(logs)
            within1 = sum(1 for x in logs if abs(x) <= 1) / len(logs)
            within2 = sum(1 for x in logs if abs(x) <= 2) / len(logs)
            print(f"  {bits:>3}-bit eff={eff:<5} n={tot}  "
                  f"log2(q*/q_scale): mean={mean:+.3f} sd={math.sqrt(var):.3f}  "
                  f"|.|<=1: {within1:.1%}  |.|<=2: {within2:.1%}")
            near = sum(v for d, v in off.items() if abs(d) <= 1) / tot
            print(f"                    index offset j*-j_scale: "
                  + " ".join(f"{d:+d}:{off[d]}" for d in sorted(off))
                  + f"   |offset|<=1: {near:.1%}")

# ---------------------------------------------------------------------------
# Experiment C -- local partial quotient vs global max_a
# ---------------------------------------------------------------------------

def exp_c():
    print("\n" + "=" * 78)
    print("Exp C -- a_local = a_{j*+1} vs global max_a, both against nu_hat")
    print("=" * 78)
    rng = random.Random(90210)
    for bits in (24, 64, 256):
        n = int(sympy.nextprime(rng.getrandbits(bits) | (1 << (bits - 1))))
        for eff in (0.05, 0.25):
            k1, k2 = k1k2(n, eff)
            nus, als, mxs, ths = [], [], [], []
            for _ in range(600):
                lam = rng.randint(2, n - 2)
                r = analyse(n, lam, k1, k2)
                if r['j_star'] < 0:
                    continue
                nus.append(r['nu_cf'])
                als.append(r['a_local'])
                mxs.append(r['max_a'])
                ths.append(r['theta'])
            print(f"  {bits:>3}-bit eff={eff:<5} n={len(nus)}")
            print(f"      spearman(nu_hat, a_local) = {spearman(nus, als):+.3f}   "
                  f"<- scale-aware")
            print(f"      spearman(nu_hat, max_a)   = {spearman(nus, mxs):+.3f}   "
                  f"<- scale-free (June invariant)")
            print(f"      spearman(nu_hat, theta*)  = {spearman(nus, ths):+.3f}")
            # bound check nu^2 >= 2 theta
            viol = sum(1 for nu, th in zip(nus, ths) if nu * nu < 2 * th - 1e-9)
            ratio = [nu * nu / (2 * th) for nu, th in zip(nus, ths) if th > 0]
            ratio.sort()
            print(f"      bound nu^2 >= 2*theta*: violations {viol}/{len(nus)}; "
                  f"ratio nu^2/(2 theta*) median={ratio[len(ratio)//2]:.3f} "
                  f"max={ratio[-1]:.3f}")

def exp_c_decay():
    """Why max_a dies as n grows, measured directly.

    a_local is the partial quotient at index j*+1, where j* is scale-locked.
    max_a is the largest over all ~O(log n) partial quotients.  The chance the
    global argmax coincides with the scale-locked index falls like 1/len(CF),
    so max_a's correlation with nu_hat must decay while a_local's does not.
    This is a noise-free measurement (no attack, no seeds).
    """
    print("\n" + "=" * 78)
    print("Exp C-decay -- max_a dies with n, a_local does not (nu_hat level)")
    print("=" * 78)
    rng = random.Random(777001)
    print(f"  {'bits':>5} {'|CF|':>6} {'sp(nu,a_local)':>15} {'sp(nu,max_a)':>13} "
          f"{'P(argmax=j*+1)':>15}")
    for bits in (16, 20, 24, 32, 48, 64, 96, 128, 192, 256, 384, 521):
        n = int(sympy.nextprime(rng.getrandbits(bits) | (1 << (bits - 1))))
        k1, k2 = k1k2(n, 0.05)
        nus, als, mxs, cls, hits, tot = [], [], [], [], 0, 0
        for _ in range(1500):
            lam = rng.randint(2, n - 2)
            r = analyse(n, lam, k1, k2)
            if r['j_star'] < 0:
                continue
            nus.append(r['nu_cf'])
            als.append(r['a_local'])
            mxs.append(r['max_a'])
            cls.append(len(r['quots']))
            tot += 1
            if r['quots'].index(max(r['quots'])) == r['j_star'] + 1:
                hits += 1
        print(f"  {bits:>5} {sum(cls)/len(cls):>6.1f} "
              f"{spearman(nus, als):>+15.3f} {spearman(nus, mxs):>+13.3f} "
              f"{hits/tot:>15.1%}")

def exp_c_threshold():
    """What a_local does the empirical C2 ceiling nu_hat <= 0.645 correspond to?"""
    print("\n" + "=" * 78)
    print("Exp C2 -- translating the empirical C2 ceiling nu_hat <= 0.645")
    print("=" * 78)
    rng = random.Random(1234567)
    n = int(sympy.nextprime(rng.getrandbits(20) | (1 << 19)))
    k1, k2 = k1k2(n, 0.05)
    buckets = {}
    for _ in range(20000):
        lam = rng.randint(2, n - 2)
        r = analyse(n, lam, k1, k2)
        if r['j_star'] < 0:
            continue
        a = min(r['a_local'], 10)
        buckets.setdefault(a, []).append(r['nu_cf'])
    print(f"  n={n} (20-bit), eff=0.05, 20000 synthetic lambda")
    print(f"  {'a_local':>8} {'count':>6} {'mean nu':>9} {'P(nu<=0.645)':>13}")
    for a in sorted(buckets):
        v = buckets[a]
        lab = f">={a}" if a == 10 else str(a)
        print(f"  {lab:>8} {len(v):>6} {sum(v)/len(v):>9.3f} "
              f"{sum(1 for x in v if x <= 0.645)/len(v):>13.1%}")

# ---------------------------------------------------------------------------
# Experiment D -- causal replication with the real attack
# ---------------------------------------------------------------------------

def first_j0_curve(lo):
    p = int(sympy.nextprime(lo))
    while True:
        if p % 3 == 1:
            eis = eisenstein_decompose(p)
            if eis is not None:
                for t in j0_traces(*eis):
                    n = p + 1 - t
                    if n < 2 or n % 3 != 1 or not sympy.isprime(n):
                        continue
                    if glv_eigenvalues(n) is None:
                        continue
                    b = identify_twist(p, n)
                    if b is None:
                        continue
                    G = find_generator(p, b, n)
                    if G is None:
                        continue
                    return p, b, n, G
        p = int(sympy.nextprime(p))

def spearman_perm_p(xs, ys, iters, rng):
    obs = spearman(xs, ys)
    ys2 = list(ys)
    hits = 0
    for _ in range(iters):
        rng.shuffle(ys2)
        if abs(spearman(xs, ys2)) >= abs(obs) - 1e-15:
            hits += 1
    return obs, (hits + 1) / (iters + 1)

def exp_d(n_lam=400, seeds_n=48, m=9, eff=0.05, bits_lo=2 ** 19, label="Exp D"):
    print("\n" + "=" * 78)
    print(f"{label} -- causal design (ONE fixed curve, synthetic lambda) with "
          f"CF invariants")
    print("=" * 78)
    p, b, n, G = first_j0_curve(bits_lo)
    k1, k2 = k1k2(n, eff)
    print(f"  fixed curve p={p} b={b} n={n}  K1={k1} K2={k2} m={m} "
          f"seeds={seeds_n} lambdas={n_lam}")
    print(f"  predicted balance scale sqrt(n*S1/S2) = "
          f"{isqrt_scale(n, *[scales(n,k1,k2)[i] for i in (0,2)])}")
    rng = random.Random(31337)
    rows = []
    t0 = time.time()
    for _ in range(n_lam):
        lam = rng.randint(2, n - 2)
        r = analyse(n, lam, k1, k2)
        wins = 0
        for s in range(seeds_n):
            seed = 1000 * m + s
            d = random.Random(seed + 7777).randint(1, n - 1)
            wins += run_experiment(p, n, lam, G, m, d, k1, k2, seed)
        r['p_hat'] = wins / seeds_n
        rows.append(r)
    print(f"  {len(rows)} lambdas x {seeds_n} seeds in {time.time()-t0:.1f}s")

    ys = [r['p_hat'] for r in rows]
    prng = random.Random(555)
    out = {}
    for name, key in (('nu_hat ', 'nu_cf'), ('a_local', 'a_local'),
                      ('max_a  ', 'max_a'), ('theta* ', 'theta'),
                      ('mu     ', None)):
        if key is None:
            xs = [mu_of(r['lam'], n) for r in rows]
        else:
            xs = [r[key] for r in rows]
        r_sp, pv = spearman_perm_p(xs, ys, 5000, prng)
        out[name.strip()] = r_sp
        star = '***' if pv < 0.001 else '**' if pv < 0.01 else '*' if pv < 0.05 else ''
        print(f"    spearman({name}, p_hat) = {r_sp:+.3f}  perm p={pv:.4f} {star}")

    print("\n  p_hat by a_local (the scale-aware partial quotient):")
    bk = {}
    for r in rows:
        bk.setdefault(min(r['a_local'], 5), []).append(r['p_hat'])
    print(f"    {'a_local':>8} {'count':>6} {'mean p_hat':>11} {'mean nu':>9}")
    for a in sorted(bk):
        v = bk[a]
        nus = [r['nu_cf'] for r in rows if min(r['a_local'], 5) == a]
        lab = f">={a}" if a == 5 else str(a)
        print(f"    {lab:>8} {len(v):>6} {sum(v)/len(v):>11.3f} "
              f"{sum(nus)/len(nus):>9.3f}")
    return n, out

def exp_d_scaling():
    """The sharp prediction: a_local's power is stable in n, max_a's decays.

    a_local is the partial quotient at the scale sqrt(n*K2/K1); max_a is the
    largest one anywhere.  As n grows the CF of lam/n gets longer, so the
    chance that the global maximum happens to sit at the relevant scale falls
    like 1/len(CF).  a_local is scale-locked and should not decay.
    """
    print("\n" + "=" * 78)
    print("Exp D-scaling -- a_local vs max_a as n grows (attack-measured)")
    print("=" * 78)
    res = []
    for bits, n_lam, seeds_n in ((20, 400, 48), (24, 300, 48), (28, 250, 32),
                                 (32, 200, 24)):
        n, out = exp_d(n_lam=n_lam, seeds_n=seeds_n, bits_lo=2 ** (bits - 1),
                       label=f"Exp D[{bits}-bit]")
        res.append((bits, n, out))
    print("\n" + "-" * 78)
    print("  SUMMARY  spearman(., p_hat) vs curve size")
    print(f"  {'bits':>5} {'n':>12} {'nu_hat':>9} {'a_local':>9} {'max_a':>9} "
          f"{'mu':>9}")
    for bits, n, out in res:
        print(f"  {bits:>5} {n:>12} {out['nu_hat']:>+9.3f} "
              f"{out['a_local']:>+9.3f} {out['max_a']:>+9.3f} {out['mu']:>+9.3f}")
    return res

# ---------------------------------------------------------------------------
# Experiment E -- secp256k1, explicitly
# ---------------------------------------------------------------------------

def exp_e():
    print("\n" + "=" * 78)
    print("Exp E -- secp256k1: the convergent that determines nu_hat")
    print("=" * 78)
    quots, convs = cf_data(SECP_LAM, SECP_N)
    print(f"  lam/n has {len(quots)} partial quotients, max_a = {max(quots)} "
          f"at index {quots.index(max(quots))}")
    print(f"  {'eff':>6} {'K1':>22} {'j*':>4} {'log2 q*':>8} {'log2 q_scale':>13} "
          f"{'a_local':>8} {'theta*':>8} {'nu_hat':>8}")
    for eff in (0.02, 0.05, 0.0993, 0.25, 0.5):
        k1, k2 = k1k2(SECP_N, eff)
        r = analyse(SECP_N, SECP_LAM, k1, k2)
        print(f"  {eff:>6} {k1:>22} {r['j_star']:>4} "
              f"{math.log2(r['q_star']):>8.2f} {math.log2(r['q_scale']):>13.2f} "
              f"{r['a_local']:>8} {r['theta']:>8.4f} {r['nu_cf']:>8.4f}")
    print("\n  The winning index moves with eff: nu_hat is a property of the")
    print("  (n, lam, K1/K2) triple, not of the curve.  Any curve-level")
    print("  invariant is therefore structurally incapable of predicting it.")

# ---------------------------------------------------------------------------

def main():
    print("=" * 78)
    print("GLV-HNP Thread 23 -- closed form for lambda_1(L2)")
    print("=" * 78)
    t0 = time.time()
    exp_a()
    exp_b()
    exp_c()
    exp_c_decay()
    exp_c_threshold()
    exp_e()
    exp_d_scaling()
    print(f"\nTotal wall time {time.time()-t0:.1f}s")
    print("Done.")

if __name__ == '__main__':
    sys.exit(main())
