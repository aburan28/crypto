"""
GLV-HNP Phase 2, Thread 23: is the K1 wall an artefact of the embedding, or
the Gaussian-heuristic (unique-SVP) density bound?

Background (RESEARCH_AUTOLAB_LOG.md, 2026-07-29, Thread 20):
  T4  found a "K1 wall": the Phase-2 attack recovers d for K1 <= 4-6 and
      fails for K1 >= 8-12, on both 12-bit historical curves.
  T4b found that at K1 = 8 more signatures do not rescue it
      (m = 8/12/16/24/32 -> 0,0,1,0,1 of 5), and concluded
      "the K1 wall is genuine".
  T5  found that the planted vector is NEVER lambda_1: the trivial vector
      n*S_D*e_m (norm n) is always 2-3x shorter.  The proposed next step
      (Thread 23) was to reformulate the lattice so that the target IS
      lambda_1 -- by projecting along e_m or by an explicit CVP -- with the
      falsifier "if the K1 wall moves outward, the reformulation is a real
      improvement; if it stays, the wall is information-theoretic".

This script answers that.

RESULT (2026-08-02 run; see glv_hnp_phase23_gh_law_output.txt):
  H23 is FALSIFIED, and so is the Thread-23 reformulation proposal.
  * U1/U3/U5: LLL succeeds far past gamma = 1 (5/5 at gamma = 1.72 on
    12-bit/2557, K1=12, eff=0.235 -- above the supposed 0.1757 ceiling).
    Parameter-free "gamma < 1" scores 0.825 on 40 17-bit cells vs a 0.575
    majority baseline, but is systematically over-conservative.  The reason
    is structural: the Phase-2 lattice is NOT GH-random.  It is m identical
    copies of one 2-D block (see glv_hnp_phase23_coset_exact.py), so its GS
    profile is nothing like a random lattice's and lambda_1^GH is meaningless
    for it.  (P2) is dead on arrival for the same reason.
  * U2: at K1=8 the two 12-bit curves behave completely differently at every
    m from 8 to 64 -- 2557 gives 5/5 throughout, 2677 gives 0-2/5 throughout.
    Same n-size, same K2, same eff.  eff is therefore NOT the controlling
    variable; the curve (i.e. lam) is.
  * U4: the Thread-23 reformulation R2 -- delete the d-column outright, which
    removes T5's trivial vector n*S_D*e_m and drops dim from 2m+2 to 2m+1 --
    gives BIT-IDENTICAL success rates in all 10 cells tested (5/5<->5/5,
    2/5<->2/5, 0/5<->0/5).  Predicted gain 1.5%, measured gain 0%.
    Thread 23's proposed reformulation is CLOSED as a dead end.
  Continued in glv_hnp_phase23_coset_exact.py (exact boundary) and
  glv_hnp_phase23_block_predictor.py (the predictor that does work).

The falsified hypothesis, kept for the record:

  H23  Recovery succeeds iff the planted vector is shorter than the
       Gaussian-heuristic first minimum of the Phase-2 lattice:

           gamma(m, n, eff) = E||v_planted|| / lambda_1^GH(L)  <  1

       with  eff = K1*K2/n  the bias strength,  K2 = isqrt(n)+1, and

           E||v_planted||^2 = m*(K1*S_K1)^2/3 + (n*S_D)^2/3
                              + m*(K2*S_K2)^2/3 + S_KANNAN^2
           det L            = (n*S_K1)^m * S_D * S_K2^m * S_KANNAN
           lambda_1^GH      = sqrt(dim/(2*pi*e)) * (det L)^(1/dim),
                              dim = 2m+2.

       Substituting S_K1 = n/K1, S_K2 = n/K2, S_D = 1, S_KANNAN = n gives
       the closed form

           gamma = C(m) * n^(1/(2m+2)) * eff^(m/(2m+2)),
           C(m)  = sqrt( 2*pi*e * (2m+4) / (3*(2m+2)) ),

       hence a data-dependent critical bias

           eff_crit(m, n) = C(m)^(-(2m+2)/m) * n^(-1/m)
                          -> 3/(2*pi*e) = 0.17573...   as m -> infinity.

  Two consequences, both sharp and both testable:
   (P1) eff_crit is INCREASING in m and saturates at 0.1757.  So T4b's
        "more data does not rescue K1=8" should be a statement about m<=32
        only: at 12-bit n, K1=8 is eff=0.157 < 0.1757, so a LARGE ENOUGH m
        must work.  Predicted crossover m ~ 45.
   (P2) Any eff > 0.1757 is unreachable at ANY m.  K1=12 (eff=0.244) must
        fail for every m, including m = 64.

Experiments:
  U1  gamma vs. measured success on the T4 K1-grid (replication + prediction)
  U2  P1: m-sweep at K1=8 (eff=0.157) out to m=64, and at K1=6
  U3  P2: m-sweep at K1=12/16 (eff>0.1757) out to m=64 -- must stay dead
  U4  Thread-23 reformulation R2: drop the d-column entirely (quotient out
      the trivial vector n*e_m), recover d by back-substitution.  Does the
      wall move?  gamma predicts a ~1.5% gain, i.e. no.
  U5  fresh 17-bit curve sweep: gamma as a classifier vs eff and lam*.

Run: python3 glv_hnp_phase23_gh_law.py
"""

import math
import random
import sympy
from fpylll import IntegerMatrix, LLL, BKZ

TWO_PI_E = 2.0 * math.pi * math.e

# ---------------------------------------------------------------------------
# Minimal EC arithmetic (verbatim from glv_hnp_phase2_lambda_threshold.py)
# ---------------------------------------------------------------------------

def modinv(a, m):
    return pow(a, -1, m)

def ec_add(P, Q, p):
    if P is None: return Q
    if Q is None: return P
    x1, y1 = P; x2, y2 = Q
    if x1 == x2:
        if (y1 + y2) % p == 0: return None
        s = 3 * x1 * x1 * modinv(2 * y1, p) % p
    else:
        s = (y2 - y1) * modinv(x2 - x1, p) % p
    x3 = (s * s - x1 - x2) % p
    y3 = (s * (x1 - x3) - y1) % p
    return (x3, y3)

def ec_mul(P, k, p):
    if k == 0: return None
    R, Q = None, P
    while k > 0:
        if k & 1: R = ec_add(R, Q, p)
        Q = ec_add(Q, Q, p)
        k >>= 1
    return R

def tonelli_shanks(n, p):
    n %= p
    if n == 0: return 0
    if pow(n, (p - 1) // 2, p) != 1: return None
    if p % 4 == 3: return pow(n, (p + 1) // 4, p)
    q, s = p - 1, 0
    while q % 2 == 0: q //= 2; s += 1
    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1: z += 1
    m, c, t, r = s, pow(z, q, p), pow(n, q, p), pow(n, (q + 1) // 2, p)
    while True:
        if t == 0: return 0
        if t == 1: return r
        i, tmp = 0, t
        while tmp != 1: tmp = tmp * tmp % p; i += 1
        b = pow(c, 1 << (m - i - 1), p)
        m, c, t, r = i, b * b % p, t * b * b % p, r * b % p

def find_generator(p, b, n, seed=12345):
    rng = random.Random(seed)
    for _ in range(10000):
        x = rng.randint(0, p - 1)
        rhs = (pow(x, 3, p) + b) % p
        y = tonelli_shanks(rhs, p)
        if y is not None and y != 0:
            P = (x, y)
            if ec_mul(P, n, p) is None:
                return P
    return None

# ---------------------------------------------------------------------------
# CM theory for j=0 (verbatim)
# ---------------------------------------------------------------------------

def eisenstein_decompose(p):
    for a in range(1, 2 * math.isqrt(p // 3) + 3):
        disc = 4 * p - 3 * a * a
        if disc < 0: break
        s = math.isqrt(disc)
        if s * s != disc: continue
        for num in (a + s, a - s):
            if num % 2 == 0:
                b = num // 2
                if b >= 0 and a * a - a * b + b * b == p:
                    return (a, b)
    return None

def j0_traces(a, b):
    return [a + b, -(a + b), a - 2 * b, -(a - 2 * b), 2 * a - b, -(2 * a - b)]

def glv_roots(n):
    """Both roots of x^2+x+1 = 0 mod n (n prime, n = 1 mod 3)."""
    if n % 3 != 1: return None
    for g in range(2, 200):
        c = pow(g, (n - 1) // 3, n)
        if c != 1 and (c * c + c + 1) % n == 0:
            return (min(c, (n - 1 - c) % n), max(c, (n - 1 - c) % n))
    return None

def lam_star(lam, n):
    return min(lam % n, (n - lam) % n) / n

def build_curve(p, n, seed=12345):
    rng = random.Random(seed)
    for _ in range(400):
        b = rng.randint(1, p - 1)
        x = rng.randint(0, p - 1)
        rhs = (pow(x, 3, p) + b) % p
        y = tonelli_shanks(rhs, p)
        if y is None or y == 0: continue
        if ec_mul((x, y), n, p) is None:
            G = find_generator(p, b, n, seed=seed)
            if G is not None:
                return (p, b, n, None, G)
    return None

def search_curves(lo, hi, want, seed=1):
    """Collect j=0 GLV curves with p in [lo,hi), n prime, n = 1 mod 3."""
    out = []
    p = int(sympy.nextprime(lo))
    while p < hi and len(out) < want:
        if p % 3 == 1:
            eis = eisenstein_decompose(p)
            if eis is not None:
                for t in j0_traces(*eis):
                    n_c = p + 1 - t
                    if n_c < 2 or n_c % 3 != 1 or not sympy.isprime(n_c):
                        continue
                    roots = glv_roots(n_c)
                    if roots is None: continue
                    cur = build_curve(p, n_c)
                    if cur is None: continue
                    out.append((p, cur[1], n_c, roots[0], cur[4]))
                    break
        p = int(sympy.nextprime(p))
    return out

# ---------------------------------------------------------------------------
# Phase-2 lattice (verbatim construction from glv_hnp_phase2_20bit.py:263)
# ---------------------------------------------------------------------------

def scales(n, k1_bound, k2_bound):
    return (n // k1_bound, 1, max(1, n // k2_bound), n)

def gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed=42):
    rng = random.Random(seed)
    sigs, attempts = [], 0
    while len(sigs) < m and attempts < 400000:
        attempts += 1
        k1 = rng.randint(0, k1_bound - 1)
        k2 = rng.randint(0, k2_bound - 1)
        k_full = (k1 + lam * k2) % n
        if k_full == 0: continue
        R = ec_mul(G, k_full, p)
        if R is None: continue
        r = R[0] % n
        if r == 0: continue
        h = rng.randint(0, n - 1)
        s = modinv(k_full, n) * (h + d_secret * r) % n
        if s == 0: continue
        s_inv = modinv(s, n)
        A, B = h * s_inv % n, r * s_inv % n
        assert (A + B * d_secret) % n == k_full
        sigs.append({'A': A, 'B': B, 'k1': k1, 'k2': k2})
    return sigs

def build_glv_lattice(sigs, n, lam, k1_bound, k2_bound):
    """Original Phase-2 Kannan embedding.  dim = 2m+2."""
    m = len(sigs)
    dim = 2 * m + 2
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    M = [[0] * dim for _ in range(dim)]
    for i in range(m):
        M[i][i] = n * S_K1
    for i in range(m):
        M[m][i] = sigs[i]['B'] * S_K1
    M[m][m] = S_D
    for i in range(m):
        M[m + 1 + i][i] = -lam * S_K1
        M[m + 1 + i][m + 1 + i] = S_K2
    for i in range(m):
        M[2 * m + 1][i] = sigs[i]['A'] * S_K1
    M[2 * m + 1][dim - 1] = S_KAN
    return M, S_K1, S_D, S_K2, S_KAN

def recover_d(M_reduced, m, n, S_KAN, d_secret):
    dim = 2 * m + 2
    for row in M_reduced:
        last = row[dim - 1]
        if abs(last) != S_KAN: continue
        sign = 1 if last > 0 else -1
        d_cand = (sign * row[m]) % n
        if d_cand and d_cand == d_secret:
            return d_cand
    return None

# ---------------------------------------------------------------------------
# Thread-23 reformulation R2: quotient out the trivial vector n*S_D*e_m by
# deleting the d-column outright.  d is recovered by back-substitution from
# any single signature once (k1_i, k2_i) are known.  dim = 2m+1 (was 2m+2).
# ---------------------------------------------------------------------------

def build_glv_lattice_R2(sigs, n, lam, k1_bound, k2_bound):
    """
    Columns: [k1_0..k1_{m-1} | k2_0..k2_{m-1} | kannan].
    Rows   : n*S_K1*e_i ; (B_i*S_K1) ; (-lam*S_K1*e_i , S_K2*e_i) ; (A_i*S_K1 , 0 , S_KAN)
    Planted vector = (k1_i*S_K1 , k2_i*S_K2 , S_KAN).
    """
    m = len(sigs)
    dim = 2 * m + 1
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    rows = []
    for i in range(m):
        r = [0] * dim; r[i] = n * S_K1; rows.append(r)
    r = [0] * dim
    for i in range(m): r[i] = sigs[i]['B'] * S_K1
    rows.append(r)
    for i in range(m):
        r = [0] * dim; r[i] = -lam * S_K1; r[m + i] = S_K2; rows.append(r)
    r = [0] * dim
    for i in range(m): r[i] = sigs[i]['A'] * S_K1
    r[dim - 1] = S_KAN
    rows.append(r)
    return rows, S_K1, S_K2, S_KAN

def recover_d_R2(M_reduced, sigs, n, lam, S_K1, S_K2, S_KAN, d_secret):
    m = len(sigs)
    dim = 2 * m + 1
    for row in M_reduced:
        last = row[dim - 1]
        if abs(last) != S_KAN: continue
        sign = 1 if last > 0 else -1
        for i in range(m):
            if row[i] % S_K1 or row[m + i] % S_K2:
                continue
            k1 = sign * row[i] // S_K1
            k2 = sign * row[m + i] // S_K2
            B = sigs[i]['B']
            if math.gcd(B, n) != 1: continue
            d_cand = (k1 - sigs[i]['A'] + lam * k2) * modinv(B, n) % n
            if d_cand and d_cand == d_secret:
                return d_cand
    return None

# ---------------------------------------------------------------------------
# Gaussian-heuristic predictor
# ---------------------------------------------------------------------------

def gamma_orig(m, n, k1_bound, k2_bound):
    """E||v_planted|| / lambda_1^GH for the original dim-(2m+2) lattice."""
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    dim = 2 * m + 2
    v2 = (m * (k1_bound * S_K1) ** 2 / 3.0 + (n * S_D) ** 2 / 3.0
          + m * (k2_bound * S_K2) ** 2 / 3.0 + S_KAN ** 2)
    log_det = m * math.log(n * S_K1) + math.log(S_D) + m * math.log(S_K2) + math.log(S_KAN)
    log_gh = 0.5 * math.log(dim / TWO_PI_E) + log_det / dim
    return math.exp(0.5 * math.log(v2) - log_gh)

def gamma_R2(m, n, k1_bound, k2_bound):
    """Same for the reformulated dim-(2m+1) lattice.
    det = S_K1^m * n^(m-1) * S_K2^m * S_KAN."""
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    dim = 2 * m + 1
    v2 = (m * (k1_bound * S_K1) ** 2 / 3.0
          + m * (k2_bound * S_K2) ** 2 / 3.0 + S_KAN ** 2)
    log_det = (m * math.log(S_K1) + (m - 1) * math.log(n)
               + m * math.log(S_K2) + math.log(S_KAN))
    log_gh = 0.5 * math.log(dim / TWO_PI_E) + log_det / dim
    return math.exp(0.5 * math.log(v2) - log_gh)

def eff_crit(m, n):
    """Closed form: solve gamma_orig = 1 for eff = K1*K2/n (S_K1=n/K1 exact)."""
    C = math.sqrt(TWO_PI_E * (2 * m + 4) / (3.0 * (2 * m + 2)))
    return C ** (-(2 * m + 2) / m) * n ** (-1.0 / m)

EFF_CRIT_INF = 3.0 / TWO_PI_E   # 0.175731...

# ---------------------------------------------------------------------------
# Experiment driver
# ---------------------------------------------------------------------------

def trial(curve, m, k1_bound, seed, variant='orig', use_bkz=False, beta=20):
    p, b, n, lam, G = curve
    k2 = math.isqrt(n) + 1
    rng = random.Random(seed ^ 0x5eed)
    d_secret = rng.randrange(1, n)
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2, seed)
    if len(sigs) < m: return None
    if variant == 'orig':
        M, S_K1, S_D, S_K2, S_KAN = build_glv_lattice(sigs, n, lam, k1_bound, k2)
        dim = 2 * m + 2
    else:
        M, S_K1, S_K2, S_KAN = build_glv_lattice_R2(sigs, n, lam, k1_bound, k2)
        dim = 2 * m + 1
    A = IntegerMatrix.from_matrix(M)
    if use_bkz:
        BKZ.reduction(A, BKZ.Param(beta))
    else:
        LLL.reduction(A)
    red = [[A[i][j] for j in range(A.ncols)] for i in range(A.nrows)]
    if variant == 'orig':
        return recover_d(red, m, n, S_KAN, d_secret) is not None
    return recover_d_R2(red, sigs, n, lam, S_K1, S_K2, S_KAN, d_secret) is not None

def rate(curve, m, k1_bound, seeds, variant='orig', use_bkz=False, beta=20):
    w = t = 0
    for s in seeds:
        r = trial(curve, m, k1_bound, s, variant, use_bkz, beta)
        if r is None: continue
        t += 1; w += 1 if r else 0
    return w, t

# ===========================================================================

if __name__ == '__main__':

    print("=" * 78)
    print("Thread 23 -- the K1 wall is the Gaussian-heuristic bound")
    print("=" * 78)
    print(f"asymptotic eff_crit = 3/(2*pi*e) = {EFF_CRIT_INF:.5f}")

    SEEDS = [42, 1234, 9999, 555, 31337]

    # Historical Phase-2 curves (RESEARCH_AUTOLAB_LOG.md 2026-06-15 / 2026-07-26)
    HIST = [
        ("8-bit/199",   211,  2, 199,  106),
        ("12-bit/2557", 2557, 2, 2659, 1755),
        ("12-bit/2677", 2677, 2, 2647, 185),
    ]
    hist = []
    for label, p, b, n, lam in HIST:
        G = find_generator(p, b, n)
        assert G is not None and (lam * lam + lam + 1) % n == 0, label
        hist.append((label, (p, b, n, lam, G)))

    # ---------------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("U1: gamma vs measured success on the T4 K1-grid (m=12, 5 seeds)")
    print("-" * 78)
    print("Thread 20 T4 measured (m as logged): 2557 -> 5/5 up to K1=8, dies at 16;")
    print("2677 -> 5/5 up to K1=4, dies at K1=8.  H23 predicts success iff gamma<1.")
    u1 = []
    for label, curve in hist[1:]:
        p, b, n, lam, G = curve
        k2 = math.isqrt(n) + 1
        print(f"\n{label}: n={n} lam*={lam_star(lam,n):.3f} K2={k2} "
              f"eff_crit(m=12)={eff_crit(12,n):.4f}")
        print(f"{'K1':>4} {'eff':>7} {'gamma':>7} {'pred':>5} {'LLL':>6} {'BKZ20':>6}")
        for k1 in (2, 3, 4, 6, 8, 12, 16, 24):
            eff = k1 * k2 / n
            g = gamma_orig(12, n, k1, k2)
            w, t = rate(curve, 12, k1, SEEDS)
            bk = "-" if w == t else "%d/%d" % rate(curve, 12, k1, SEEDS, use_bkz=True)
            print(f"{k1:>4} {eff:>7.3f} {g:>7.3f} {'OK' if g < 1 else 'no':>5} "
                  f"{str(w)+'/'+str(t):>6} {bk:>6}")
            u1.append({'label': label, 'k1': k1, 'eff': eff, 'gamma': g,
                       'ok': w == t, 'any': w > 0})

    # ---------------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("U2 (P1): does MORE DATA cross the K1=8 wall?  eff_crit rises with m.")
    print("-" * 78)
    print("T4b tested m<=32 at K1=8 and concluded 'the K1 wall is genuine'.")
    print("H23: eff_crit(m,n) is increasing in m, so a large enough m must work")
    print("whenever eff < 0.17573.  Predicted crossover printed as 'pred'.")
    for label, curve in hist[1:]:
        p, b, n, lam, G = curve
        k2 = math.isqrt(n) + 1
        for k1 in (6, 8):
            eff = k1 * k2 / n
            print(f"\n{label}  K1={k1}  eff={eff:.4f}  "
                  f"({'reachable' if eff < EFF_CRIT_INF else 'UNREACHABLE'}, "
                  f"asymptote {EFF_CRIT_INF:.4f})")
            print(f"{'m':>4} {'dim':>5} {'eff_crit':>9} {'gamma':>7} {'pred':>5} {'LLL':>6}")
            for m in (8, 12, 16, 24, 32, 48, 64):
                ec = eff_crit(m, n)
                g = gamma_orig(m, n, k1, k2)
                w, t = rate(curve, m, k1, SEEDS)
                print(f"{m:>4} {2*m+2:>5} {ec:>9.4f} {g:>7.3f} "
                      f"{'OK' if g < 1 else 'no':>5} {str(w)+'/'+str(t):>6}")

    # ---------------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("U3 (P2): eff > 0.17573 must be dead at EVERY m (the true ceiling)")
    print("-" * 78)
    label, curve = hist[1]
    p, b, n, lam, G = curve
    k2 = math.isqrt(n) + 1
    for k1 in (12, 16):
        eff = k1 * k2 / n
        print(f"\n{label}  K1={k1}  eff={eff:.4f} > {EFF_CRIT_INF:.4f}")
        print(f"{'m':>4} {'gamma':>7} {'pred':>5} {'LLL':>6}")
        for m in (12, 24, 48, 64):
            g = gamma_orig(m, n, k1, k2)
            w, t = rate(curve, m, k1, SEEDS)
            print(f"{m:>4} {g:>7.3f} {'OK' if g < 1 else 'no':>5} {str(w)+'/'+str(t):>6}")

    # ---------------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("U4: Thread-23 reformulation R2 (d-column deleted, dim 2m+1).")
    print("-" * 78)
    print("The trivial vector n*S_D*e_m of T5 no longer exists in R2.  gamma_R2")
    print("predicts only a ~1.5% gain, i.e. the wall should NOT move.")
    for label, curve in hist[1:]:
        p, b, n, lam, G = curve
        k2 = math.isqrt(n) + 1
        print(f"\n{label} (m=12)")
        print(f"{'K1':>4} {'eff':>7} {'g_orig':>7} {'g_R2':>7} {'ratio':>6} "
              f"{'orig':>6} {'R2':>6}")
        for k1 in (2, 4, 6, 8, 12):
            eff = k1 * k2 / n
            go, gr = gamma_orig(12, n, k1, k2), gamma_R2(12, n, k1, k2)
            wo, to = rate(curve, 12, k1, SEEDS, variant='orig')
            wr, tr = rate(curve, 12, k1, SEEDS, variant='R2')
            print(f"{k1:>4} {eff:>7.3f} {go:>7.3f} {gr:>7.3f} {gr/go:>6.3f} "
                  f"{str(wo)+'/'+str(to):>6} {str(wr)+'/'+str(tr):>6}")

    # ---------------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("U5: fresh 17-bit curves -- gamma as a classifier vs eff and lam*")
    print("-" * 78)
    curves = search_curves(1 << 16, 1 << 17, 10)
    print(f"{len(curves)} curves collected")
    rows = []
    print(f"{'p':>7} {'n':>7} {'lam*':>6} {'K1':>4} {'eff':>7} {'gamma':>7} "
          f"{'pred':>5} {'LLL':>6}")
    for (p, b, n, lam, G) in curves:
        k2 = math.isqrt(n) + 1
        for eff_t in (0.05, 0.10, 0.15, 0.25):
            k1 = max(2, int(eff_t * n / k2))
            eff = k1 * k2 / n
            g = gamma_orig(12, n, k1, k2)
            w, t = rate((p, b, n, lam, G), 12, k1, SEEDS)
            rows.append({'gamma': g, 'eff': eff, 'lam_star': lam_star(lam, n),
                         'ok': w == t})
            print(f"{p:>7} {n:>7} {lam_star(lam,n):>6.3f} {k1:>4} {eff:>7.3f} "
                  f"{g:>7.3f} {'OK' if g<1 else 'no':>5} {str(w)+'/'+str(t):>6}")

    def best_acc(rows, key, greater_is_success):
        vals = sorted({r[key] for r in rows})
        best = (0.0, None)
        for i in range(len(vals) + 1):
            thr = vals[0] - 1 if i == 0 else vals[i - 1]
            c = sum(1 for r in rows
                    if (r[key] > thr if greater_is_success else r[key] <= thr) == r['ok'])
            if c / len(rows) > best[0]:
                best = (c / len(rows), thr)
        return best

    if rows:
        base = max(sum(1 for r in rows if r['ok']), sum(1 for r in rows if not r['ok'])) / len(rows)
        print(f"\nclassifier accuracy over {len(rows)} (curve,K1) cells:")
        print(f"  majority baseline           : {base:.3f}")
        for key, gis in (('gamma', False), ('eff', False), ('lam_star', True)):
            a, thr = best_acc(rows, key, gis)
            print(f"  best threshold on {key:<10}: {a:.3f}  (thr={thr:.4f})")
        fixed = sum(1 for r in rows if (r['gamma'] < 1.0) == r['ok']) / len(rows)
        print(f"  PARAMETER-FREE gamma < 1    : {fixed:.3f}")

    print("\n" + "=" * 78)
    print("done")
    print("=" * 78)
