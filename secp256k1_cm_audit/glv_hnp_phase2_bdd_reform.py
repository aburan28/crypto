"""
GLV-HNP Phase 2, Thread 23: reformulate the Phase-2 lattice so that the
planted vector is actually the target of the search.

Background (RESEARCH_AUTOLAB_LOG.md, 2026-07-29 run, T5):
  In the Phase-2 Kannan lattice (build_glv_lattice, dim 2m+2) the shortest
  vector is ALWAYS the trivial vector n*S_D*e_m -- 100% of its energy in the
  d column, |sv[m]|/n = 1.0000 exactly -- and it is 2-3x shorter than the
  planted vector on every curve tested, success and failure alike.  The
  planted vector is therefore never lambda_1, and recovery is a BDD/coset
  condition, not an SVP condition.  That run proposed Thread 23: quotient out
  the trivial direction and solve BDD directly.

  Falsifier as stated on 2026-07-29:
    "if sv/pv rises above 1 after the reformulation and the K1 wall in T4
     moves outward on the lam*=0.07 curve (currently K1 ~ 4-6), the
     reformulation is a real improvement; if the wall stays at K1 ~ 4-6, then
     the wall is information-theoretic and Phase 2 is at its ceiling."

THE REFORMULATION
  The trivial vector exists because d gets its own lattice COLUMN (scale S_D),
  and d is uniform in [0,n), so the planted vector carries an inherently
  O(n) coordinate while n*S_D*e_m carries exactly n.  No choice of S_D helps:
  both scale linearly in S_D.

  Fix: eliminate d algebraically instead of giving it a column.  With
  c_i = B_i * B_0^{-1} mod n and w_i = k1_i + lam*k2_i = A_i + B_i*d (mod n),

        w_i - c_i * w_0 = A_i - c_i * A_0  =: C_i   (mod n),   i = 1..m-1

  which is d-free.  The unknowns are just (k1, k2) in Z^{2m}, all small.  The
  solution set of the homogeneous system is a rank-2m lattice L0 with the
  explicit triangular basis below (no HNF needed), and the planted (k1,k2)
  lies in the coset t0 + L0 for the explicit particular solution
  t0 = (0, C_1, .., C_{m-1}, 0_m).  Recovery is then a genuine CVP:

        find v in L0 minimising || (v - T) ||,   T = -t0
        error e = v - T = (k1_i * S_K1, k2_i * S_K2)
        d = (k1_0 + lam*k2_0 - A_0) * B_0^{-1}  mod n

  Dimension drops 2m+2 -> 2m, the trivial vector is gone by construction, and
  d is read off the coefficient side rather than being a coordinate.

  Two further variants are tested:
    * CENTERING.  k1_i in [0,K1) and k2_i in [0,K2) are one-sided, so the
      planted error has a large mean.  Shifting the target by the box centre
      makes the error entries uniform on [-n/2, n/2] instead of [0, n),
      cutting E[e_j^2] from n^2/3 to n^2/12: a factor 2 in ||e||, which by the
      Gaussian heuristic should buy a factor ~4 in the bias budget K1*K2/n.
    * EXACT CVP (fpylll CVP.closest_vector) vs Babai nearest-plane.  This is
      what separates "our reduction is too weak" from "the planted vector is
      genuinely not the closest vector".  If exact CVP returns something
      strictly closer than the planted error, the instance is
      information-theoretically lost and no better reduction can help.

EXPERIMENTS
  E1  construction check: coset membership + recovery on a known-good curve.
  E2  gap diagnostic: ||e_found||/||e_planted|| and lambda_1(L0)/(2||e_pl||),
      i.e. the T5 "sv/pv" statistic recomputed for the new lattice.
  E3  K1 wall sweep on both historical 12-bit curves, 4 methods x 5 seeds.
  E4  eff = K1*K2/n scaling on fresh 17-bit curves: does the wall sit where
      the Gaussian heuristic says it should?

Run: python3 glv_hnp_phase2_bdd_reform.py
"""

import math
import random
from fractions import Fraction

import sympy
from fpylll import CVP, IntegerMatrix, LLL

# ---------------------------------------------------------------------------
# Minimal EC arithmetic (short Weierstrass, a=0) -- verbatim from
# glv_hnp_phase2_lambda_threshold.py so comparisons stay exact.
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
        if y is None or y == 0: continue
        P = (x, y)
        if ec_mul(P, n, p) is None and P is not None:
            return P
    return None

def eisenstein_decompose(p):
    """p = a^2 - a*b + b^2 for p = 1 mod 3 (Cornacchia on x^2+3y^2=4p)."""
    if p % 3 != 1: return None
    r = tonelli_shanks((p - 3) % p, p)
    if r is None: return None
    # Cornacchia for 4p = x^2 + 3 y^2
    a, b = 2 * p, None
    # sqrt(-3) mod p -> lift to a solution of x^2 = -3 mod 4p via x odd
    x0 = r * modinv(1, p) % p
    # standard: run the Euclidean algorithm on (2p, x) with x^2 = -3 mod 4p
    for cand in (x0, p - x0):
        x = cand if cand % 2 == 1 else cand + p
        u, v = 2 * p, x
        lim = math.isqrt(4 * p)
        while v > lim:
            u, v = v, u % v
        rem = 4 * p - v * v
        if rem % 3 != 0: continue
        w2 = rem // 3
        w = math.isqrt(w2)
        if w * w != w2: continue
        # 4p = v^2 + 3 w^2 ; a = (v+w)/2, b = w  (a^2-ab+b^2 = p)
        if (v + w) % 2 != 0: continue
        a = (v + w) // 2
        b = w
        if a * a - a * b + b * b == p:
            return (a, b)
    return None

def j0_traces(a, b):
    return [2*a - b, -2*a + b, -(a + b), a + b, 2*b - a, a - 2*b]

def glv_roots(n):
    sq = tonelli_shanks((n - 3) % n, n)
    if sq is None: return None, None
    inv2 = modinv(2, n)
    r1 = (n - 1 + sq) * inv2 % n
    r2 = (n - 1 + (n - sq)) * inv2 % n
    if (r1 * r1 + r1 + 1) % n != 0: r1, r2 = r2, r1
    if (r1 * r1 + r1 + 1) % n != 0: return None, None
    return r1, r2

def lam_star(lam, n):
    lam %= n
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

# ---------------------------------------------------------------------------
# Scales and signatures (identical to the 2026-07-29 run)
# ---------------------------------------------------------------------------

def scales(n, k1_bound, k2_bound):
    return (n // k1_bound, 1, max(1, n // k2_bound), n)

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

# ---------------------------------------------------------------------------
# BASELINE: the 2026-07-29 Kannan lattice (dim 2m+2)
# ---------------------------------------------------------------------------

def build_kannan_lattice(sigs, n, lam, k1_bound, k2_bound):
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
    return M

def kannan_recover(sigs, n, lam, k1_bound, k2_bound, d_secret, use_bkz=False, beta=20):
    from fpylll import BKZ
    m = len(sigs)
    dim = 2 * m + 2
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    A = IntegerMatrix.from_matrix(build_kannan_lattice(sigs, n, lam, k1_bound, k2_bound))
    if use_bkz:
        BKZ.reduction(A, BKZ.Param(beta))
    else:
        LLL.reduction(A)
    for i in range(dim):
        last = A[i][dim - 1]
        if abs(last) != S_KAN: continue
        sign = 1 if last > 0 else -1
        d_cand = (sign * A[i][m]) % n
        if d_cand == d_secret:
            return True
    return False

# ---------------------------------------------------------------------------
# REFORMULATION: d-free CVP lattice (dim 2m)
# ---------------------------------------------------------------------------

def build_bdd_instance(sigs, n, lam, k1_bound, k2_bound):
    """
    Returns (basis, t0, S_K1, S_K2) where
      basis: 2m x 2m integer basis of the SCALED homogeneous solution lattice
             L0 = {(u,v) in Z^m x Z^m :
                     u_i + lam*v_i - c_i*(u_0 + lam*v_0) = 0 mod n, i=1..m-1}
      t0:    the scaled particular solution (0, C_1..C_{m-1}, 0_m)
    Coordinates: 0..m-1 are the u (=k1) columns scaled by S_K1,
                 m..2m-1 are the v (=k2) columns scaled by S_K2.
    """
    m = len(sigs)
    S_K1, _S_D, S_K2, _S_KAN = scales(n, k1_bound, k2_bound)
    lam %= n
    B0inv = modinv(sigs[0]['B'] % n, n)
    c = [sigs[i]['B'] * B0inv % n for i in range(m)]         # c[0] == 1
    C = [(sigs[i]['A'] - c[i] * sigs[0]['A']) % n for i in range(m)]

    dim = 2 * m
    basis = []
    # r_A : u_0 = 1, u_i = c_i
    row = [0] * dim
    row[0] = S_K1
    for i in range(1, m):
        row[i] = c[i] * S_K1
    basis.append(row)
    # r_B : v_0 = 1, u_i = c_i * lam
    row = [0] * dim
    row[m] = S_K2
    for i in range(1, m):
        row[i] = (c[i] * lam % n) * S_K1
    basis.append(row)
    # r_C,i : v_i = 1, u_i = -lam            (i = 1..m-1)
    for i in range(1, m):
        row = [0] * dim
        row[m + i] = S_K2
        row[i] = ((-lam) % n) * S_K1
        basis.append(row)
    # r_D,i : u_i = n                        (i = 1..m-1)
    for i in range(1, m):
        row = [0] * dim
        row[i] = n * S_K1
        basis.append(row)
    assert len(basis) == dim

    t0 = [0] * dim
    for i in range(1, m):
        t0[i] = C[i] * S_K1
    return basis, t0, S_K1, S_K2

def planted_error(sigs, k1_bound, k2_bound, n):
    """The scaled error vector the attack is supposed to find."""
    m = len(sigs)
    S_K1, _, S_K2, _ = scales(n, k1_bound, k2_bound)
    e = [0] * (2 * m)
    for i in range(m):
        e[i] = sigs[i]['k1'] * S_K1
        e[m + i] = sigs[i]['k2'] * S_K2
    return e

def in_lattice_coset(basis, t0, e):
    """Check e - t0 lies in the lattice spanned by basis (exact, via solve)."""
    import sympy as sp
    Bm = sp.Matrix(basis).T                      # columns = basis vectors
    rhs = sp.Matrix([e[i] - t0[i] for i in range(len(e))])
    sol = Bm.solve(rhs)
    return all(x.q == 1 for x in sol)

def center_vector(m, k1_bound, k2_bound, n):
    S_K1, _, S_K2, _ = scales(n, k1_bound, k2_bound)
    ctr = [0] * (2 * m)
    for i in range(m):
        ctr[i] = (k1_bound // 2) * S_K1
        ctr[m + i] = (k2_bound // 2) * S_K2
    return ctr

def norm(v):
    return math.sqrt(sum(float(x) * float(x) for x in v))

# --- Babai nearest plane, exact rational GSO ------------------------------

def gso_exact(basis):
    n = len(basis)
    bstar, mu = [], [[Fraction(0)] * n for _ in range(n)]
    B = []
    for i in range(n):
        v = [Fraction(x) for x in basis[i]]
        for j in range(i):
            if B[j] == 0:
                mu[i][j] = Fraction(0); continue
            num = sum(Fraction(basis[i][k]) * bstar[j][k] for k in range(len(v)))
            mu[i][j] = num / B[j]
            for k in range(len(v)):
                v[k] -= mu[i][j] * bstar[j][k]
        bstar.append(v)
        B.append(sum(x * x for x in v))
    return bstar, mu, B

def babai_nearest_plane(basis, target):
    """Closest-ish vector to `target`; basis should already be LLL-reduced."""
    n = len(basis)
    bstar, mu, B = gso_exact(basis)
    w = [Fraction(x) for x in target]
    coeffs = [0] * n
    for i in range(n - 1, -1, -1):
        if B[i] == 0:
            continue
        cnum = sum(w[k] * bstar[i][k] for k in range(len(w)))
        c = cnum / B[i]
        ci = int((2 * c.numerator + c.denominator) // (2 * c.denominator))  # round-half-up
        coeffs[i] = ci
        for k in range(len(w)):
            w[k] -= ci * Fraction(basis[i][k])
    v = [sum(coeffs[i] * basis[i][k] for i in range(n)) for k in range(len(target))]
    return v

def lll_reduce_rows(basis):
    A = IntegerMatrix.from_matrix([[int(x) for x in r] for r in basis])
    LLL.reduction(A)
    return [[A[i][j] for j in range(A.ncols)] for i in range(A.nrows)], A

# --- the four attack methods ----------------------------------------------

def bdd_attack(sigs, n, lam, k1_bound, k2_bound, d_secret,
               centered=True, exact=False, want_diag=False):
    """
    Solve the d-free CVP and try to recover d.
    Returns (success, diag) with diag = dict of gap statistics (if want_diag).
    """
    m = len(sigs)
    basis, t0, S_K1, S_K2 = build_bdd_instance(sigs, n, lam, k1_bound, k2_bound)
    red, A = lll_reduce_rows(basis)

    shift = center_vector(m, k1_bound, k2_bound, n) if centered else [0] * (2 * m)
    # want lattice v close to T ; then e = v - T + shift
    T = [shift[i] - t0[i] for i in range(2 * m)]

    if exact:
        v = list(CVP.closest_vector(A, tuple(int(x) for x in T)))
    else:
        v = babai_nearest_plane(red, T)

    e = [v[i] - T[i] + shift[i] for i in range(2 * m)]

    # recover d from coordinate 0 of the error
    ok = False
    if e[0] % S_K1 == 0 and e[m] % S_K2 == 0:
        k1_0 = e[0] // S_K1
        k2_0 = e[m] // S_K2
        B0inv = modinv(sigs[0]['B'] % n, n)
        d_cand = (k1_0 + lam * k2_0 - sigs[0]['A']) * B0inv % n
        ok = (d_cand == d_secret)

    diag = None
    if want_diag:
        e_pl = planted_error(sigs, k1_bound, k2_bound, n)
        e_pl_c = [e_pl[i] - shift[i] for i in range(2 * m)]     # centred planted error
        e_found_c = [e[i] - shift[i] for i in range(2 * m)]
        lam1 = norm(red[0])
        diag = {
            'n_pl': norm(e_pl_c),
            'n_found': norm(e_found_c),
            'ratio': norm(e_found_c) / norm(e_pl_c) if norm(e_pl_c) else float('inf'),
            'lam1': lam1,
            'lam1_over_2pl': lam1 / (2 * norm(e_pl_c)) if norm(e_pl_c) else float('inf'),
        }
    return ok, diag

METHODS = {
    'kannan-LLL':   lambda s, n, l, k1, k2, d: kannan_recover(s, n, l, k1, k2, d),
    'bdd-babai':    lambda s, n, l, k1, k2, d: bdd_attack(s, n, l, k1, k2, d, centered=False)[0],
    'bdd-babai-C':  lambda s, n, l, k1, k2, d: bdd_attack(s, n, l, k1, k2, d, centered=True)[0],
    'bdd-cvp-C':    lambda s, n, l, k1, k2, d: bdd_attack(s, n, l, k1, k2, d, centered=True, exact=True)[0],
}

SEEDS = [42, 1234, 9999, 555, 31337]

def success_rate(curve, m, k1_bound, method, seeds=SEEDS):
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    fn = METHODS[method]
    ok = 0
    for sd in seeds:
        d_secret = random.Random(sd ^ 0xABCDEF).randrange(1, n)
        sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed=sd)
        if len(sigs) < m: continue
        try:
            if fn(sigs, n, lam, k1_bound, k2_bound, d_secret):
                ok += 1
        except Exception as exc:                       # noqa: BLE001
            print(f"      [{method} k1={k1_bound} seed={sd}] EXC {type(exc).__name__}: {exc}")
    return ok

# ---------------------------------------------------------------------------

HIST = [
    # label,             p,    b, n,    lam,  K1, m
    ("12-bit/2557",      2557, 2, 2659, 1755, 8,  8),
    ("12-bit/2677 FAIL", 2677, 2, 2647, 185,  8,  10),
]

def load_hist():
    out = []
    for label, p, b, n, lam, k1, m in HIST:
        G = find_generator(p, b, n)
        assert G is not None, f"no generator for {label}"
        assert (lam * lam + lam + 1) % n == 0, f"bad lam for {label}"
        out.append((label, (p, b, n, lam, G), k1, m))
    return out

def search_curves(lo, hi, want=6, seed=7):
    """j=0 GLV curves with p in [lo,hi), n prime, n = 1 mod 3."""
    found = []
    p = sympy.nextprime(lo)
    while p < hi and len(found) < want:
        if p % 3 == 1:
            eis = eisenstein_decompose(p)
            if eis is not None:
                for t in j0_traces(*eis):
                    n_c = p + 1 - t
                    if n_c < 8 or not sympy.isprime(n_c) or n_c % 3 != 1:
                        continue
                    r1, r2 = glv_roots(n_c)
                    if r1 is None: continue
                    lam = min(r1, r2)
                    cur = build_curve(p, n_c, seed=seed)
                    if cur is None: continue
                    found.append((f"{p}/{n_c}", (p, cur[1], n_c, lam, cur[4])))
                    break
        p = sympy.nextprime(p)
    return found


def main():
    print("=" * 78)
    print("GLV-HNP Phase 2, Thread 23: d-free BDD reformulation")
    print("=" * 78)

    hist = load_hist()

    # ---------------- E1: construction check ----------------
    print("\n--- E1  construction check (coset membership + recovery) ---")
    for label, curve, k1_default, m in hist:
        p, b, n, lam, G = curve
        k2_bound = math.isqrt(n) + 1
        k1_bound = 2
        d_secret = random.Random(42 ^ 0xABCDEF).randrange(1, n)
        sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed=42)
        basis, t0, S_K1, S_K2 = build_bdd_instance(sigs, n, lam, k1_bound, k2_bound)
        e_pl = planted_error(sigs, k1_bound, k2_bound, n)
        member = in_lattice_coset(basis, t0, e_pl)
        okB, _ = bdd_attack(sigs, n, lam, k1_bound, k2_bound, d_secret, centered=True)
        print(f"  {label:20s} m={m:2d} dim(Kannan)={2*m+2:3d} dim(BDD)={2*m:3d} "
              f"planted in t0+L0: {member}   recover(K1=2): {okB}")

    # ---------------- E2: gap diagnostic ----------------
    print("\n--- E2  gap diagnostic on the new lattice ---")
    print("  ratio = ||e_found||/||e_planted||   (<1 => planted is NOT closest:")
    print("                                       information-theoretically lost)")
    print(f"  {'curve':20s} {'K1':>3s} {'||e_pl||':>11s} {'||e_cvp||':>11s} "
          f"{'ratio':>7s} {'lam1':>11s} {'lam1/2||e||':>11s} {'rec':>4s}")
    e2_rows = []
    for label, curve, _k1, m in hist:
        p, b, n, lam, G = curve
        k2_bound = math.isqrt(n) + 1
        for k1_bound in (2, 4, 8, 16, 32):
            d_secret = random.Random(42 ^ 0xABCDEF).randrange(1, n)
            sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed=42)
            if len(sigs) < m: continue
            ok, dg = bdd_attack(sigs, n, lam, k1_bound, k2_bound, d_secret,
                                centered=True, exact=True, want_diag=True)
            print(f"  {label:20s} {k1_bound:3d} {dg['n_pl']:11.3e} {dg['n_found']:11.3e} "
                  f"{dg['ratio']:7.4f} {dg['lam1']:11.3e} {dg['lam1_over_2pl']:11.4f} "
                  f"{'Y' if ok else 'n':>4s}")
            e2_rows.append((label, k1_bound, dg, ok))

    # ---------------- E3: K1 wall sweep ----------------
    print("\n--- E3  K1 wall sweep, 4 methods x 5 seeds (successes out of 5) ---")
    K1_GRID = [2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64]
    for label, curve, _k1, m in hist:
        p, b, n, lam, G = curve
        print(f"\n  {label}   n={n} lam={lam} lam*={lam_star(lam, n):.4f} m={m}")
        hdr = "    " + f"{'method':14s}" + "".join(f"{k:>5d}" for k in K1_GRID)
        print(hdr)
        for meth in METHODS:
            cells = []
            for k1_bound in K1_GRID:
                cells.append(success_rate(curve, m, k1_bound, meth))
            print("    " + f"{meth:14s}" + "".join(f"{c:>5d}" for c in cells))

    # ---------------- E4: eff scaling on fresh 17-bit curves ----------------
    print("\n--- E4  eff = K1*K2/n scaling, fresh 17-bit curves, m=12 ---")
    curves17 = search_curves(2**16, 2**17, want=4, seed=7)
    print(f"  found {len(curves17)} curves: " +
          ", ".join(f"{lb}(lam*={lam_star(c[3], c[2]):.3f})" for lb, c in curves17))
    m = 12
    EFFS = [0.05, 0.10, 0.20, 0.40, 0.80]
    print(f"    {'curve':16s} {'lam*':>6s} " +
          " ".join(f"{'eff=' + str(e):>16s}" for e in EFFS))
    print(f"    {'':16s} {'':>6s} " +
          " ".join(f"{'kan/bab/babC/cvp':>16s}" for _ in EFFS))
    for lb, curve in curves17:
        p, b, n, lam, G = curve
        k2_bound = math.isqrt(n) + 1
        cells = []
        for eff in EFFS:
            k1_bound = max(2, int(eff * n / k2_bound))
            res = [success_rate(curve, m, k1_bound, meth, seeds=SEEDS[:3])
                   for meth in METHODS]
            cells.append("/".join(str(r) for r in res))
        print(f"    {lb:16s} {lam_star(lam, n):6.3f} " +
              " ".join(f"{c:>16s}" for c in cells))
    print("    (each cell: kannan-LLL / bdd-babai / bdd-babai-C / bdd-cvp-C, out of 3)")

    # ---------------- E5: Gaussian-heuristic decoding prediction ----------------
    print("\n--- E5  Gaussian-heuristic prediction of the information-theoretic wall ---")
    print("  det(L0) = S_K1^m * S_K2^m * n^(m-1);  centred ||e|| ~ n*sqrt(2m/12)")
    print("  N_pred  = vol(ball_2m(||e||)) / det(L0)  = expected # lattice points")
    print("            at least as close to the target as the planted vector.")
    print("  N_pred << 1  =>  planted vector is the unique CVP solution (recoverable)")
    print("  N_pred >> 1  =>  information-theoretically lost, no reduction can help")
    print(f"  {'curve':20s} {'K1':>4s} {'m':>3s} {'N_pred':>12s}")
    for label, curve, _k1, m in hist:
        p, b, n, lam, G = curve
        k2_bound = math.isqrt(n) + 1
        for k1_bound in (2, 4, 8, 16, 32):
            S_K1, _, S_K2, _ = scales(n, k1_bound, k2_bound)
            logdet = m * math.log(S_K1) + m * math.log(S_K2) + (m - 1) * math.log(n)
            R = n * math.sqrt(2 * m / 12)
            d = 2 * m
            logvol = (d / 2) * math.log(math.pi) - math.lgamma(d / 2 + 1) + d * math.log(R)
            print(f"  {label:20s} {k1_bound:4d} {m:3d} {math.exp(logvol - logdet):12.3e}")

    print("\ndone.")


if __name__ == "__main__":
    main()
