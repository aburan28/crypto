"""
GLV-HNP Phase 2 -- Thread 23: is the K1 wall the Gaussian-heuristic uSVP
boundary, and can any reformulation of the lattice move it?

Background
----------
2026-07-29 (Thread 20) established two things about the Phase-2 lattice:

  (i)  lam* = min(lam, n-lam)/n is NOT a viability threshold (falsified);
       nu_hat = lambda_1(L2)/sqrt(det L2) is a real residual separator
       (AUC 0.935 against the June C1/C2 classes).
  (ii) The planted vector is NEVER lambda_1: the trivial vector n*S_D*e_m
       is in the lattice and is shorter for every m >= 1.  Recovery is a
       BDD/coset condition, not an SVP condition.

The proposed next step (log line 6089) was:

    "Thread 23 -- reformulate the Phase-2 lattice so the target is lambda_1.
     Falsifier: if sv/pv rises above 1 after the reformulation AND the K1
     wall in T4 moves outward on the lam*=0.07 curve (currently K1~4-6),
     the reformulation is a real improvement; if the wall stays at K1~4-6,
     then the wall is information-theoretic and Phase 2 is at its ceiling."

This script executes that falsifier, and answers a question the falsifier's
second branch got wrong: the wall is NOT information-theoretic.  It is the
Gaussian-heuristic wall, which sits a CONSTANT factor 2*pi*e/3 = 5.69x
(2.51 bits) inside the information-theoretic boundary.

The closed form (derived in EXP A, verified numerically against the exact
floored scales)
----------------------------------------------------------------------
With S_K1 = n//K1, S_D = 1, S_K2 = n//K2, S_KANNAN = n, dim = 2m+2:

    det L   = (n*S_K1)^m * S_D * S_K2^m * S_KANNAN  =  n^(3m+1) / (K1*K2)^m
    E||v||^2 = n^2 * (2m/3 + 4/3)
    GH(L)   = sqrt(dim/(2*pi*e)) * (det L)^(1/dim)

Writing eff = K1*K2/n, the ratio gamma = E||v_planted|| / GH(L) is

    gamma^2 = 2*pi*e * (m+2)/(3*(m+1)) * n^(1/(m+1)) * eff^(m/(m+1))

so gamma = 1 (the uSVP boundary) is at

    eff*(m, n) = [ 3(m+1) / (2*pi*e*(m+2)) ]^((m+1)/m) * n^(-1/m)
    eff*(inf)  = 3 / (2*pi*e) = 0.17565            <-- hard asymptotic ceiling

Every scale S_* cancels out of gamma at leading order, which is why no
re-scaling reformulation can move the wall.

Experiments
-----------
A  closed form vs exact-floored numeric gamma; the eff*(m,n) table.
B  retro-fit: eff*, gamma, nu_hat against the 2026-07-29 T4 K1 grid.
C  bit-length invariance: empirical eff wall at 12/16/20/24 bits vs eff*.
D  m-invariance: GH says the wall creeps to 0.1756; information theory says
   m >= log n / log(1/eff) always suffices.  These disagree sharply.
E  the two reformulations the falsifier asked for:
     V1  S_D rescale (makes the planted vector lambda_1 -- sv/pv -> 1)
     V2  drop the Kannan row, solve BDD directly by Babai nearest-plane
   measured as sv/pv AND as wall position on the lam*=0.07 curve.
F  what DOES move the wall: BKZ beta sweep.

Run: python3 glv_hnp_phase2_gh_boundary.py
"""

import math
import random
import sympy

from fpylll import IntegerMatrix, LLL, BKZ, GSO

TWO_PI_E = 2.0 * math.pi * math.e

# ---------------------------------------------------------------------------
# Minimal EC arithmetic + CM curve construction
# (copied verbatim from glv_hnp_phase2_lambda_threshold.py so that the
#  comparison against the 2026-07-29 numbers is exact)
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
        i, t2 = 0, t
        while t2 != 1:
            t2 = t2 * t2 % p; i += 1
        b = pow(c, 1 << (m - i - 1), p)
        m, c = i, b * b % p
        t, r = t * c % p, r * b % p

def find_point(p, b):
    for x in range(p):
        rhs = (pow(x, 3, p) + b) % p
        y = tonelli_shanks(rhs, p)
        if y is not None and y != 0:
            return (x, y)
    return None

def find_generator(p, b, n, seed=12345):
    rng = random.Random(seed)
    for _ in range(600):
        x = rng.randint(0, p - 1)
        rhs = (pow(x, 3, p) + b) % p
        y = tonelli_shanks(rhs, p)
        if y is None or y == 0: continue
        P = (x, y)
        if ec_mul(P, n, p) is None and P is not None:
            return P
    return None

def eisenstein_decompose(p):
    for a in range(1, 2 * math.isqrt(p // 3) + 3):
        disc = 4 * p - 3 * a * a
        if disc < 0: break
        s = math.isqrt(disc)
        if s * s != disc: continue
        for num in [a + s, a - s]:
            if num % 2 == 0:
                b = num // 2
                if b >= 0 and a * a - a * b + b * b == p:
                    return (a, b)
    return None

def j0_traces(a, b):
    return [2*a - b, -2*a + b, -(a + b), a + b, 2*b - a, a - 2*b]

def glv_roots(n):
    """Both roots of x^2+x+1 = 0 mod n, or None."""
    sq = tonelli_shanks((n - 3) % n, n)
    if sq is None: return None
    inv2 = modinv(2, n)
    r1 = (n - 1 + sq) * inv2 % n
    r2 = (n - 1 + (n - sq)) * inv2 % n
    if (r1 * r1 + r1 + 1) % n != 0: return None
    if (r2 * r2 + r2 + 1) % n != 0: return None
    return (min(r1, r2), max(r1, r2))

def lam_star(lam, n):
    return min(lam, n - lam) / n

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

def search_curves_at_bits(bits, want, seed=12345):
    """Collect `want` j=0 GLV curves with p ~ 2^bits, n prime, n = 1 mod 3."""
    out = []
    p = int(sympy.nextprime(2 ** (bits - 1) + 2 ** (bits - 2)))
    while len(out) < want and p < 2 ** bits:
        if p % 3 == 1:
            eis = eisenstein_decompose(p)
            if eis is not None:
                for t in j0_traces(*eis):
                    n_cand = p + 1 - t
                    if n_cand < 8 or n_cand % 3 != 1: continue
                    if not sympy.isprime(n_cand): continue
                    roots = glv_roots(n_cand)
                    if roots is None: continue
                    cur = build_curve(p, n_cand, seed=seed)
                    if cur is None: continue
                    out.append((p, cur[1], n_cand, roots[0], cur[4]))
                    break
        p = int(sympy.nextprime(p))
    return out

# ---------------------------------------------------------------------------
# Phase-2 lattice (verbatim), plus the two Thread-23 variants
# ---------------------------------------------------------------------------

def scales(n, k1_bound, k2_bound, sd=1):
    return (n // k1_bound, sd, max(1, n // k2_bound), n)

def gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed=42):
    rng = random.Random(seed)
    sigs = []
    attempts = 0
    while len(sigs) < m and attempts < 200000:
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
        A = h * s_inv % n
        B = r * s_inv % n
        assert (A + B * d_secret) % n == k_full
        sigs.append({'A': A, 'B': B, 'k1': k1, 'k2': k2})
    return sigs

def build_glv_lattice(sigs, n, lam, k1_bound, k2_bound, sd=1):
    """V0 / V1.  sd=1 reproduces the 2026-06-15 construction exactly."""
    m = len(sigs)
    dim = 2 * m + 2
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound, sd)
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
    return M

def planted_vector(sigs, d_secret, n, k1_bound, k2_bound, sd=1):
    m = len(sigs)
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound, sd)
    v = [0] * (2 * m + 2)
    for i in range(m):
        v[i] = sigs[i]['k1'] * S_K1
        v[m + 1 + i] = sigs[i]['k2'] * S_K2
    v[m] = d_secret * S_D
    v[2 * m + 1] = S_KAN
    return v

def build_bdd_lattice(sigs, n, lam, k1_bound, k2_bound):
    """V2: the same lattice with the Kannan row removed -- dim 2m+1.
    The instance becomes an explicit BDD/CVP with target t[i] = -A_i*S_K1."""
    m = len(sigs)
    dim = 2 * m + 1
    S_K1, S_D, S_K2, _ = scales(n, k1_bound, k2_bound)
    M = [[0] * dim for _ in range(dim)]
    for i in range(m):
        M[i][i] = n * S_K1
    for i in range(m):
        M[m][i] = sigs[i]['B'] * S_K1
    M[m][m] = S_D
    for i in range(m):
        M[m + 1 + i][i] = -lam * S_K1
        M[m + 1 + i][m + 1 + i] = S_K2
    t = [0] * dim
    for i in range(m):
        t[i] = -sigs[i]['A'] * S_K1
    return M, t

def norm(v):
    return math.sqrt(sum(float(x) * float(x) for x in v))

def recover_d(M_reduced, m, n, S_KANNAN, d_secret, sd=1):
    dim = 2 * m + 2
    for row in M_reduced:
        last = row[dim - 1]
        if abs(last) != S_KANNAN: continue
        sign = 1 if last > 0 else -1
        val = sign * row[m]
        if val % sd: continue
        d_cand = (val // sd) % n
        if d_cand == 0: continue
        if d_cand == d_secret:
            return d_cand
    return None

# ---------------------------------------------------------------------------
# Predictors
# ---------------------------------------------------------------------------

def log_det_lattice(m, n, k1_bound, k2_bound, sd=1):
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound, sd)
    return (m * math.log(n * S_K1) + math.log(S_D)
            + m * math.log(S_K2) + math.log(S_KAN))

def gh_norm(m, n, k1_bound, k2_bound, sd=1):
    """Gaussian heuristic lambda_1 of the (2m+2)-dim Phase-2 lattice."""
    dim = 2 * m + 2
    ld = log_det_lattice(m, n, k1_bound, k2_bound, sd)
    return math.sqrt(dim / TWO_PI_E) * math.exp(ld / dim)

def planted_norm_expected(m, n, k1_bound, k2_bound, sd=1):
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound, sd)
    return math.sqrt(m * (k1_bound * S_K1) ** 2 / 3.0
                     + (n * S_D) ** 2 / 3.0
                     + m * (k2_bound * S_K2) ** 2 / 3.0
                     + S_KAN ** 2)

def gamma_exact(m, n, k1_bound, k2_bound, sd=1):
    """gamma = E||v_planted|| / GH(L), using the exact floored scales."""
    return planted_norm_expected(m, n, k1_bound, k2_bound, sd) / \
           gh_norm(m, n, k1_bound, k2_bound, sd)

def gamma_closed(m, n, eff):
    """Closed form: gamma^2 = 2pi e (m+2)/(3(m+1)) n^(1/(m+1)) eff^(m/(m+1))."""
    return math.sqrt(TWO_PI_E * (m + 2) / (3.0 * (m + 1))
                     * n ** (1.0 / (m + 1)) * eff ** (m / (m + 1.0)))

def eff_star(m, n):
    """eff at which gamma = 1 (closed form)."""
    return (3.0 * (m + 1) / (TWO_PI_E * (m + 2))) ** ((m + 1.0) / m) * n ** (-1.0 / m)

def eff_star_numeric(m, n, k2_bound, sd=1):
    """eff at gamma=1 using the exact floored scales (bisection on K1)."""
    lo, hi = 1.0, float(n)
    for _ in range(80):
        mid = (lo + hi) / 2
        k1 = max(1, int(round(mid)))
        g = gamma_exact(m, n, k1, k2_bound, sd)
        if g < 1.0: lo = mid
        else: hi = mid
    return lo * k2_bound / n

def m_infotheoretic(n, eff):
    """Counting bound: m*log n >= log n + m*log(K1*K2)."""
    if eff >= 1.0: return float('inf')
    return math.log(n) / math.log(1.0 / eff)

def gauss_reduce_2d(u, v):
    def nrm2(w): return w[0] * w[0] + w[1] * w[1]
    def dot(w, z): return w[0] * z[0] + w[1] * z[1]
    if nrm2(u) > nrm2(v): u, v = v, u
    while True:
        num, den = dot(v, u), nrm2(u)
        if den == 0: break
        q = (2 * num + den) // (2 * den) if num >= 0 else -((-2 * num + den) // (2 * den))
        v = (v[0] - q * u[0], v[1] - q * u[1])
        if nrm2(v) >= nrm2(u): break
        u, v = v, u
    return u

def nu_hat(n, lam, k1_bound, k2_bound):
    """Thread 20b residual separator: lambda_1(L2)/sqrt(det L2)."""
    S_K1, _, S_K2, _ = scales(n, k1_bound, k2_bound)
    w = gauss_reduce_2d((n * S_K1, 0), (-(lam % n) * S_K1, S_K2))
    return math.sqrt(w[0] ** 2 + w[1] ** 2) / math.sqrt(float(n) * S_K1 * S_K2)

# ---------------------------------------------------------------------------
# Trials
# ---------------------------------------------------------------------------

def trial_v0(curve, m, d_secret, k1_bound, seed, sd=1, use_bkz=False, beta=20):
    """Returns (ok, planted_norm, shortest_norm, trivial_norm)."""
    p, b, n, lam, G = curve
    k2b = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2b, seed)
    if len(sigs) < m: return None
    M = build_glv_lattice(sigs, n, lam, k1_bound, k2b, sd)
    dim = 2 * m + 2
    A = IntegerMatrix.from_matrix(M)
    if use_bkz: BKZ.reduction(A, BKZ.Param(beta))
    else: LLL.reduction(A)
    red = [[A[i][j] for j in range(dim)] for i in range(dim)]
    _, S_D, _, S_KAN = scales(n, k1_bound, k2b, sd)
    ok = recover_d(red, m, n, S_KAN, d_secret, sd) is not None
    pn = norm(planted_vector(sigs, d_secret, n, k1_bound, k2b, sd))
    sn = min(norm(r) for r in red)
    return (ok, pn, sn, float(n) * S_D)

def trial_v2_bdd(curve, m, d_secret, k1_bound, seed):
    """Babai nearest-plane on the Kannan-free lattice.  Returns (ok, dist, gh)."""
    p, b, n, lam, G = curve
    k2b = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2b, seed)
    if len(sigs) < m: return None
    M, t = build_bdd_lattice(sigs, n, lam, k1_bound, k2b)
    dim = 2 * m + 1
    A = IntegerMatrix.from_matrix(M)
    LLL.reduction(A)
    Mg = GSO.Mat(A); Mg.update_gso()
    coeffs = Mg.babai([int(x) for x in t])
    w = [sum(coeffs[r] * A[r][c] for r in range(dim)) for c in range(dim)]
    e = [w[c] - t[c] for c in range(dim)]
    S_K1, S_D, S_K2, _ = scales(n, k1_bound, k2b)
    d_cand = (w[m] // S_D) % n
    return (d_cand == d_secret, norm(e), dim)

def rate_v0(curve, m, k1_bound, seeds, sd=1, use_bkz=False, beta=20):
    wins, ratios = 0, []
    for seed in seeds:
        d = random.Random(seed + 7777).randint(1, curve[2] - 1)
        r = trial_v0(curve, m, d, k1_bound, seed, sd, use_bkz, beta)
        if r is None: continue
        ok, pn, sn, tn = r
        wins += bool(ok)
        ratios.append(sn / pn)
    return wins, len(seeds), (sum(ratios) / len(ratios) if ratios else float('nan'))

def rate_v2(curve, m, k1_bound, seeds):
    wins = 0
    for seed in seeds:
        d = random.Random(seed + 7777).randint(1, curve[2] - 1)
        r = trial_v2_bdd(curve, m, d, k1_bound, seed)
        if r is None: continue
        wins += bool(r[0])
    return wins, len(seeds)


SEEDS = [42, 1234, 9999, 555, 31337]

# Historical Phase-2 anchor curves (RESEARCH_AUTOLAB_LOG 2026-06-15 / 07-26 / 07-29)
HIST = [
    ("12-bit/2557", 2557, 2, 2659, 1755),   # lam* = 0.340
    ("12-bit/2677", 2677, 2, 2647, 185),    # lam* = 0.070  (the T4 hard curve)
]

def hist_curves():
    out = []
    for label, p, b, n, lam in HIST:
        G = find_generator(p, b, n)
        assert G is not None, f"no generator for {label}"
        assert (lam * lam + lam + 1) % n == 0, f"bad lam for {label}"
        out.append((label, (p, b, n, lam, G)))
    return out


print("=" * 78)
print("Thread 23 -- the Phase-2 wall is the Gaussian-heuristic uSVP boundary")
print("=" * 78)

CURVES = hist_curves()

# ===========================================================================
# EXP A -- closed form vs exact floored scales; the eff* table
# ===========================================================================
print("\n" + "-" * 78)
print("EXP A: gamma = E||v_planted|| / GH(L),  and the boundary eff*(m,n)")
print("-" * 78)
print("gamma^2 = 2*pi*e * (m+2)/(3(m+1)) * n^(1/(m+1)) * eff^(m/(m+1))")
print(f"eff*(inf) = 3/(2*pi*e) = {3.0/TWO_PI_E:.5f}   <-- asymptotic ceiling\n")

print(f"{'n':>10} {'m':>4} {'K1':>5} {'eff':>8} {'gamma_exact':>12} "
      f"{'gamma_closed':>13} {'rel.err':>9}")
for n in (2647, 131101, 16777259):
    k2b = math.isqrt(n) + 1
    for m in (8, 12, 20):
        for k1 in (2, 4, 8, 16):
            eff = k1 * k2b / n
            ge = gamma_exact(m, n, k1, k2b)
            gc = gamma_closed(m, n, eff)
            print(f"{n:>10} {m:>4} {k1:>5} {eff:>8.4f} {ge:>12.4f} "
                  f"{gc:>13.4f} {abs(ge-gc)/ge:>8.2%}")

print(f"\neff*(m, n) -- closed form (cf) vs exact-floored bisection (num):")
print(f"{'n':>10} " + " ".join(f"{'m='+str(m):>16}" for m in (6, 8, 12, 16, 24, 32)))
for n in (2647, 131101, 16777259):
    k2b = math.isqrt(n) + 1
    cells = []
    for m in (6, 8, 12, 16, 24, 32):
        cells.append(f"{eff_star(m, n):.4f}/{eff_star_numeric(m, n, k2b):.4f}")
    print(f"{n:>10} " + " ".join(f"{c:>16}" for c in cells))

print("\nInformation-theoretic comparison (counting bound m >= log n / log(1/eff)):")
print(f"{'eff':>8} {'m_info(n=2647)':>16} {'m_info(n=2^24)':>16} {'gamma<1 at m=12?':>18}")
for eff in (0.05, 0.10, 0.1756, 0.30, 0.60, 0.90):
    mi1 = m_infotheoretic(2647, eff)
    mi2 = m_infotheoretic(2 ** 24, eff)
    g = gamma_closed(12, 2647, eff)
    print(f"{eff:>8.4f} {mi1:>16.1f} {mi2:>16.1f} "
          f"{('YES' if g < 1 else 'no') + f' (g={g:.2f})':>18}")
print("\n=> information theory permits ANY eff < 1 given enough signatures;")
print("   the GH boundary caps eff at 3/(2*pi*e) = 0.1756 no matter how large m is.")
print(f"   gap = 2*pi*e/3 = {TWO_PI_E/3.0:.3f}x = {math.log2(TWO_PI_E/3.0):.3f} bits.")


# ===========================================================================
# EXP B -- retro-fit against the 2026-07-29 T4 K1 grid
# ===========================================================================
print("\n" + "-" * 78)
print("EXP B: the 2026-07-29 T4 grid, annotated with eff, gamma, nu_hat")
print("-" * 78)
print("T4 held m=12, K2=isqrt(n)+1, and swept K1.  Prediction: the wall sits")
print("where gamma crosses 1, i.e. eff ~ eff*(12, n), modulo the nu_hat residual.\n")

M_T4 = 12
K1_GRID = [2, 3, 4, 6, 8, 12, 16, 24]
b_rows = []
for label, curve in CURVES:
    p, b, n, lam, G = curve
    k2b = math.isqrt(n) + 1
    es = eff_star(M_T4, n)
    print(f"{label}  n={n}  lam*={lam_star(lam, n):.4f}  K2={k2b}  "
          f"eff*(m=12)={es:.4f}")
    print(f"  {'K1':>4} {'eff':>8} {'gamma':>8} {'nu_hat':>8} {'LLL':>6} {'sv/pv':>7}")
    for k1 in K1_GRID:
        eff = k1 * k2b / n
        g = gamma_exact(M_T4, n, k1, k2b)
        nh = nu_hat(n, lam, k1, k2b)
        w, t, r = rate_v0(curve, M_T4, k1, SEEDS)
        b_rows.append((label, k1, eff, g, nh, w, t))
        print(f"  {k1:>4} {eff:>8.4f} {g:>8.3f} {nh:>8.4f} "
              f"{str(w)+'/'+str(t):>6} {r:>7.3f}")
    print()

print("gamma at the empirical wall (last K1 with >=3/5, first K1 with <=1/5):")
for label, _ in CURVES:
    rows = [r for r in b_rows if r[0] == label]
    last_ok = max([r for r in rows if r[5] >= 3], key=lambda r: r[1], default=None)
    first_bad = min([r for r in rows if r[5] <= 1], key=lambda r: r[1], default=None)
    lo = f"K1={last_ok[1]} gamma={last_ok[3]:.3f} eff={last_ok[2]:.4f}" if last_ok else "-"
    hi = f"K1={first_bad[1]} gamma={first_bad[3]:.3f} eff={first_bad[2]:.4f}" if first_bad else "-"
    print(f"  {label:<14} last-pass [{lo}]   first-fail [{hi}]")


# ===========================================================================
# EXP C -- bit-length invariance of the wall
# ===========================================================================
print("\n" + "-" * 78)
print("EXP C: bit-length invariance -- empirical eff wall at 12/16/20/24 bits")
print("-" * 78)
print("If the wall were a fixed K1 (or a fixed lam* band) it would not track")
print("eff*(m,n).  eff*(12,n) FALLS with n like n^(-1/12), so the wall must too.\n")

M_C = 12
EFF_GRID = [0.02, 0.035, 0.05, 0.07, 0.10, 0.14, 0.20, 0.28]
print(f"{'bits':>5} {'n':>10} {'eff*(12,n)':>11}  " +
      " ".join(f"{e:>6.3f}" for e in EFF_GRID))
c_rows = []
for bits in (12, 16, 20, 24):
    curves = search_curves_at_bits(bits, 3)
    for (p, b, n, lam, G) in curves:
        curve = (p, b, n, lam, G)
        k2b = math.isqrt(n) + 1
        cells = []
        for eff in EFF_GRID:
            k1 = max(1, int(round(eff * n / k2b)))
            w, t, _ = rate_v0(curve, M_C, k1, SEEDS)
            cells.append(w)
            c_rows.append((bits, n, eff, k1, gamma_exact(M_C, n, k1, k2b), w, t))
        print(f"{bits:>5} {n:>10} {eff_star(M_C, n):>11.4f}  " +
              " ".join(f"{c:>5}/5" for c in cells))

print("\nempirical wall (highest eff with mean wins >= 3/5, pooled over the")
print("3 curves at each bit length) vs eff*(12,n):")
for bits in (12, 16, 20, 24):
    rows = [r for r in c_rows if r[0] == bits]
    ns = sorted({r[1] for r in rows})
    pooled = {}
    for r in rows:
        pooled.setdefault(r[2], []).append(r[5])
    ok_effs = [e for e, ws in sorted(pooled.items())
               if sum(ws) / len(ws) >= 3.0]
    emp = max(ok_effs) if ok_effs else float('nan')
    pred = sum(eff_star(M_C, n) for n in ns) / len(ns)
    print(f"  {bits:>2} bits: empirical wall eff <= {emp:.3f}   "
          f"predicted eff* = {pred:.3f}   ratio {emp/pred:.2f}")


# ===========================================================================
# EXP D -- m-invariance: GH vs the counting bound
# ===========================================================================
print("\n" + "-" * 78)
print("EXP D: does more data move the wall?  (GH says barely; counting says yes)")
print("-" * 78)
label_d, curve_d = CURVES[1]      # 12-bit/2677, the T4 hard curve
p, b, n, lam, G = curve_d
k2b = math.isqrt(n) + 1
K1_D = 8                          # eff ~ 0.157, the 2026-07-29 T4b setting
eff_d = K1_D * k2b / n
print(f"curve {label_d}  n={n}  K1={K1_D}  K2={k2b}  eff={eff_d:.4f}")
print(f"counting bound: m >= {m_infotheoretic(n, eff_d):.1f} signatures suffice "
      f"information-theoretically\n")
print(f"{'m':>4} {'dim':>5} {'eff*(m,n)':>10} {'gamma':>8} {'LLL':>6} {'BKZ(30)':>9}")
for m in (6, 8, 12, 16, 24, 32, 40):
    g = gamma_exact(m, n, K1_D, k2b)
    w, t, _ = rate_v0(curve_d, m, K1_D, SEEDS)
    wb, tb, _ = rate_v0(curve_d, m, K1_D, SEEDS, use_bkz=True, beta=30)
    print(f"{m:>4} {2*m+2:>5} {eff_star(m, n):>10.4f} {g:>8.3f} "
          f"{str(w)+'/'+str(t):>6} {str(wb)+'/'+str(tb):>9}")


# ===========================================================================
# EXP E -- the two reformulations the Thread-23 falsifier asked for
# ===========================================================================
print("\n" + "-" * 78)
print("EXP E: reformulations.  Falsifier: does sv/pv rise above 1 AND does")
print("       the K1 wall move outward on the lam*=0.07 curve?")
print("-" * 78)

print("\nE1 -- V1: rescale S_D so the planted vector becomes lambda_1.")
print("     ||n*S_D*e_m|| = n*S_D  vs  ||v_planted||^2 = n^2(2m/3+1) + d^2 S_D^2,")
print("     so S_D > sqrt((2m/3+1)/(1-(d/n)^2)) makes the planted vector shortest.")
print("     Cost: it also inflates ||v_planted|| relative to GH.\n")
M_E = 12
print(f"{'curve':<14} {'K1':>4} " +
      " ".join(f"{'S_D='+str(s):>12}" for s in (1, 2, 4, 8, 16, 64)))
for label, curve in CURVES:
    p, b, n, lam, G = curve
    k2b = math.isqrt(n) + 1
    for k1 in (4, 8):
        cells = []
        for sd in (1, 2, 4, 8, 16, 64):
            w, t, r = rate_v0(curve, M_E, k1, SEEDS, sd=sd)
            cells.append(f"{w}/{t} sv/pv{r:.2f}")
        print(f"{label:<14} {k1:>4} " + " ".join(f"{c:>12}" for c in cells))

print("\nE2 -- V2: drop the Kannan row, solve BDD directly (Babai nearest-plane)")
print("     on the (2m+1)-dim lattice.  No trivial-vector competition at all.\n")
print(f"{'curve':<14} {'K1':>4} {'eff':>8} {'gamma':>8} {'V0 Kannan+LLL':>15} "
      f"{'V2 Babai BDD':>14}")
for label, curve in CURVES:
    p, b, n, lam, G = curve
    k2b = math.isqrt(n) + 1
    for k1 in (2, 3, 4, 6, 8, 12):
        eff = k1 * k2b / n
        g = gamma_exact(M_E, n, k1, k2b)
        w0, t0, _ = rate_v0(curve, M_E, k1, SEEDS)
        w2, t2 = rate_v2(curve, M_E, k1, SEEDS)
        print(f"{label:<14} {k1:>4} {eff:>8.4f} {g:>8.3f} "
              f"{str(w0)+'/'+str(t0):>15} {str(w2)+'/'+str(t2):>14}")


# ===========================================================================
# EXP F -- what DOES move the wall: stronger reduction
# ===========================================================================
print("\n" + "-" * 78)
print("EXP F: BKZ beta sweep -- the only knob that moves the GH wall")
print("-" * 78)
label_f, curve_f = CURVES[1]
p, b, n, lam, G = curve_f
k2b = math.isqrt(n) + 1
print(f"curve {label_f}  n={n}  m={M_E}  dim={2*M_E+2}\n")
print(f"{'K1':>4} {'eff':>8} {'gamma':>8} " +
      " ".join(f"{'b='+str(x):>7}" for x in (2, 20, 30, 40)))
for k1 in (4, 6, 8, 12, 16):
    eff = k1 * k2b / n
    g = gamma_exact(M_E, n, k1, k2b)
    cells = []
    for beta in (2, 20, 30, 40):
        if beta == 2:
            w, t, _ = rate_v0(curve_f, M_E, k1, SEEDS)
        else:
            w, t, _ = rate_v0(curve_f, M_E, k1, SEEDS, use_bkz=True, beta=beta)
        cells.append(f"{w}/{t}")
    print(f"{k1:>4} {eff:>8.4f} {g:>8.3f} " + " ".join(f"{c:>7}" for c in cells))


# ===========================================================================
# EXP G -- two-factor model: gamma sets the scale, nu_hat the curve offset
# ===========================================================================
print("\n" + "-" * 78)
print("EXP G: is (gamma, nu_hat) a better predictor than either alone?")
print("-" * 78)
print("Pooled cells: EXP B (2 anchor curves x 8 K1) + EXP C (12 curves x 8 eff),")
print("each cell = 5 seeds.  A cell is 'easy' if wins >= 3/5.")
print("Score s = gamma + c*nu_hat; predict easy iff s < threshold.")
print("Reported: best achievable accuracy over all thresholds, and AUC.\n")

cells = []   # (gamma, nu_hat, easy)
for label, curve in CURVES:
    p, b, n, lam, G = curve
    k2b = math.isqrt(n) + 1
    for k1 in K1_GRID:
        w, t, _ = rate_v0(curve, M_T4, k1, SEEDS)
        cells.append((gamma_exact(M_T4, n, k1, k2b), nu_hat(n, lam, k1, k2b), w >= 3))
for bits in (12, 16, 20, 24):
    for (p, b, n, lam, G) in search_curves_at_bits(bits, 3):
        curve = (p, b, n, lam, G)
        k2b = math.isqrt(n) + 1
        for eff in EFF_GRID:
            k1 = max(1, int(round(eff * n / k2b)))
            w, t, _ = rate_v0(curve, M_C, k1, SEEDS)
            cells.append((gamma_exact(M_C, n, k1, k2b), nu_hat(n, lam, k1, k2b), w >= 3))

n_easy = sum(1 for c in cells if c[2])
print(f"{len(cells)} cells, {n_easy} easy / {len(cells)-n_easy} hard; "
      f"majority baseline = {max(n_easy, len(cells)-n_easy)/len(cells):.3f}")

def auc(scores, labels):
    """AUC for 'lower score => easy'."""
    pos = [s for s, l in zip(scores, labels) if l]
    neg = [s for s, l in zip(scores, labels) if not l]
    if not pos or not neg: return float('nan')
    wins = sum((1.0 if a < b else 0.5 if a == b else 0.0)
               for a in pos for b in neg)
    return wins / (len(pos) * len(neg))

def best_acc(scores, labels):
    cand = sorted(set(scores))
    best = 0.0
    for th in cand:
        acc = sum(1 for s, l in zip(scores, labels) if (s < th) == l) / len(labels)
        best = max(best, acc)
    return best

labels = [c[2] for c in cells]
print(f"\n{'predictor':<28} {'AUC':>7} {'best acc':>9}")
g_only = [c[0] for c in cells]
nu_only = [c[1] for c in cells]
print(f"{'gamma alone':<28} {auc(g_only, labels):>7.3f} {best_acc(g_only, labels):>9.3f}")
print(f"{'nu_hat alone':<28} {auc(nu_only, labels):>7.3f} {best_acc(nu_only, labels):>9.3f}")
print()
best_c, best_a = None, -1.0
for c in [0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 3.0]:
    sc = [g + c * nh for g, nh, _ in cells]
    a, ac = auc(sc, labels), best_acc(sc, labels)
    print(f"{'gamma + ' + f'{c:.2f}' + '*nu_hat':<28} {a:>7.3f} {ac:>9.3f}")
    if a > best_a: best_a, best_c = a, c
print(f"\nbest combined weight c = {best_c:.2f} (AUC {best_a:.3f})")

print("\n" + "=" * 78)
print("Done.")
print("=" * 78)
