"""
GLV-HNP Phase 2 — Thread 23: does eliminating d from the lattice move the K1 wall?

Background (RESEARCH_AUTOLAB_LOG.md, 2026-07-29 entry, T5)
----------------------------------------------------------
The Phase-2 lattice (glv_hnp_phase2_lambda_threshold.py:254, build_glv_lattice)
has dim 2m+2 with a dedicated coordinate for the *large* unknown d, scaled by
S_D = 1 (scales(), line 227).  Because row m carries S_D on the diagonal and
B_i*S_K1 on the k1-columns, the combination

    n * row_m  -  sum_i B_i * (modular row i)   =   (0, ..., 0, n*S_D, 0, ..., 0)

is a lattice vector of norm exactly n*S_D = n, while

    ||v_planted|| ~ n * sqrt(2m/3 + 4/3)   (~3n for m = 12).

So the trivial vector n*S_D*e_m is ~3x shorter than the planted vector on every
instance.  Thread 20/T5 measured this: |sv[m]|/n = 1.0000 on every curve tested,
Kannan coordinate 0, i.e. LLL's shortest output is always this trivial vector and
never the planted one.  The 2026-07-29 entry concluded recovery is a BDD/coset
condition rather than a unique-SVP condition, and proposed:

  "Thread 23 — reformulate the Phase-2 lattice so the target is lambda_1.
   Falsifier: if sv/pv rises above 1 after the reformulation AND the K1 wall in
   T4 moves outward on the lam*=0.07 curve (currently K1 ~ 4-6), the
   reformulation is a real improvement; if the wall stays at K1 ~ 4-6, then the
   wall is ... Phase 2 is at its ceiling."

This script builds the reformulated lattice and runs that falsifier.

The reformulation (d-eliminated lattice, dim 2m+1)
--------------------------------------------------
d is a *full-size* unknown, so giving it a lattice coordinate is a design flaw:
it adds ~n to ||v_planted|| and it manufactures the rival vector n*S_D*e_m.  The
classical HNP treatment eliminates d instead.  With pivot i = 0, from
A_i + d*B_i = k1_i + lam*k2_i (mod n):

    d    = B_0^{-1} (k1_0 + lam*k2_0 - A_0)                        (mod n)
    k1_i = h_i + g_i*k1_0 + g_i*lam*k2_0 - lam*k2_i   (i >= 1)     (mod n)
      g_i = B_i * B_0^{-1} mod n,     h_i = A_i - g_i*A_0 mod n.

Free integer parameters: k1_0, k2_0, ..., k2_{m-1}.  All of them are *small*.
Coordinates (dim D = 2m+1):
    col i        (0 <= i < m)   :  k1_i * S_K1
    col m+j      (0 <= j < m)   :  k2_j * S_K2
    col 2m                      :  Kannan * S_KANNAN
Generators:
    R_k1_0 :  k1_0 += 1  =>  k1_i += g_i          (i >= 1)
    R_k2_0 :  k2_0 += 1  =>  k1_i += g_i*lam      (i >= 1)
    R_k2_j :  k2_j += 1  =>  k1_j -= lam          (j >= 1)
    m modular rows  n*S_K1*e_i
    const row (h_i*S_K1 on i>=1, S_KANNAN on the last column)
Same scales S_K1 = n//K1, S_K2 = n//K2, S_KANNAN = n as the original, so the
comparison is like-for-like.

Determinants (both with S_D = 1, S_KANNAN = n):
    old:  dim 2m+2,  det = n^m     * S_K1^m * S_K2^m * S_KANNAN
    new:  dim 2m+1,  det = n^(m-1) * S_K1^m * S_K2^m * S_KANNAN
The new lattice drops one dimension and one factor of n; ||target|| drops from
n*sqrt(2m/3+4/3) to n*sqrt(2m/3+1).  Whether that is enough to move the K1 wall
is exactly what E3 measures.

Experiments
-----------
  E1  correctness: planted vector is in the new lattice; d recovers on the easy
      historical curve; new lattice has no n*e_m-style trivial vector.
  E2  geometry on the three historical curves: ||planted||, shortest reduced
      norm, sv/pv, Gaussian-heuristic ratio tau = ||planted|| / GH, old vs new.
  E3  THE FALSIFIER.  T4's K1 grid re-run with both lattices, same curves, same
      seeds, same m.  Does the wall move outward?
  E4  control for the 2026-07-29 diagnosis: sweep S_D in the OLD lattice.  If the
      trivial vector n*S_D*e_m were the binding constraint, raising S_D (which
      lengthens it relative to the rest) should help.
  E5  scaling: independent 17-bit curves, K1 wall for both formulations.

Run: python3 glv_hnp_phase2_thread23.py
"""

import math
import random

import sympy
from fpylll import IntegerMatrix, LLL, BKZ

# ---------------------------------------------------------------------------
# Helpers — copied verbatim from glv_hnp_phase2_lambda_threshold.py so that the
# old-lattice arm is bit-identical to the 2026-07-29 run.  (That module runs its
# whole T1-T5 suite at import time, so it cannot be imported as a library.)
# ---------------------------------------------------------------------------

def modinv(a, m):
    return pow(a, -1, m)

def ec_add(P, Q, p):
    if P is None: return Q
    if Q is None: return P
    x1, y1 = P; x2, y2 = Q
    if x1 == x2 and (y1 + y2) % p == 0: return None
    if P == Q:
        if y1 == 0: return None
        lam = 3 * x1 * x1 % p * modinv(2 * y1 % p, p) % p
    else:
        lam = (y2 - y1) % p * modinv((x2 - x1) % p, p) % p
    x3 = (lam * lam - x1 - x2) % p
    return (x3, (lam * (x1 - x3) - y1) % p)

def ec_mul(P, k, p):
    R, Q = None, P
    while k > 0:
        if k & 1: R = ec_add(R, Q, p)
        Q = ec_add(Q, Q, p); k >>= 1
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
    for _ in range(20000):
        x = rng.randint(0, p - 1)
        rhs = (pow(x, 3, p) + b) % p
        y = tonelli_shanks(rhs, p)
        if y is not None and y != 0:
            P = (x, y)
            if ec_mul(P, n, p) is None:
                return P
    return None

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
    return [2*a - b, -2*a + b, -(a + b), a + b, 2*b - a, a - 2*b]

def glv_roots(n):
    sq = tonelli_shanks((n - 3) % n, n)
    if sq is None: return None
    inv2 = modinv(2, n)
    r1 = (n - 1 + sq) * inv2 % n
    r2 = (n - 1 + (n - sq)) * inv2 % n
    if (r1 * r1 + r1 + 1) % n != 0: r1, r2 = r2, r1
    if (r1 * r1 + r1 + 1) % n != 0: return None
    return (min(r1, r2), max(r1, r2))

def lam_star(lam, n):
    return min(lam % n, n - (lam % n)) / n

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

def norm(v):
    return math.sqrt(sum(x * x for x in v))

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
# OLD lattice (verbatim, dim 2m+2) — S_D is exposed so E4 can sweep it
# ---------------------------------------------------------------------------

def build_glv_lattice(sigs, n, lam, k1_bound, k2_bound, s_d_mult=1):
    m = len(sigs)
    dim = 2 * m + 2
    S_K1, S_D, S_K2, S_KANNAN = scales(n, k1_bound, k2_bound)
    S_D *= s_d_mult
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
    M[2 * m + 1][dim - 1] = S_KANNAN
    return M

def planted_vector_old(sigs, d_secret, n, k1_bound, k2_bound, s_d_mult=1):
    m = len(sigs)
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    S_D *= s_d_mult
    v = [0] * (2 * m + 2)
    for i in range(m):
        v[i] = sigs[i]['k1'] * S_K1
        v[m + 1 + i] = sigs[i]['k2'] * S_K2
    v[m] = d_secret * S_D
    v[2 * m + 1] = S_KAN
    return v

def recover_d_old(M_reduced, m, n, S_KANNAN, d_secret, s_d_mult=1):
    dim = 2 * m + 2
    for row in M_reduced:
        last = row[dim - 1]
        if abs(last) != S_KANNAN: continue
        sign = 1 if last > 0 else -1
        val = sign * row[m]
        if val % s_d_mult: continue
        d_cand = (val // s_d_mult) % n
        if d_cand == 0: continue
        if d_cand == d_secret:
            return d_cand
    return None


# ---------------------------------------------------------------------------
# NEW lattice (Thread 23, d eliminated, dim 2m+1)
# ---------------------------------------------------------------------------

def build_glv_lattice_nod(sigs, n, lam, k1_bound, k2_bound):
    """d-eliminated Phase-2 lattice.  Pivot signature index 0.

    Columns: [k1_0..k1_{m-1}] [k2_0..k2_{m-1}] [kannan],  dim = 2m+1.
    Returns (rows, aux) with aux = (g, h, B0inv, A0) for recovery.
    """
    m = len(sigs)
    D = 2 * m + 1
    S_K1, _S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    lam %= n
    B0inv = modinv(sigs[0]['B'] % n, n)
    A0 = sigs[0]['A'] % n
    g = [0] * m
    h = [0] * m
    for i in range(1, m):
        g[i] = sigs[i]['B'] % n * B0inv % n
        h[i] = (sigs[i]['A'] - g[i] * A0) % n

    rows = []
    # R_k1_0 : k1_0 += 1  =>  k1_i += g_i
    r = [0] * D
    r[0] = S_K1
    for i in range(1, m):
        r[i] = g[i] * S_K1
    rows.append(r)
    # R_k2_0 : k2_0 += 1  =>  k1_i += g_i*lam
    r = [0] * D
    for i in range(1, m):
        r[i] = g[i] * lam % n * S_K1
    r[m] = S_K2
    rows.append(r)
    # R_k2_j (j >= 1) : k2_j += 1  =>  k1_j -= lam
    for j in range(1, m):
        r = [0] * D
        r[j] = -lam * S_K1
        r[m + j] = S_K2
        rows.append(r)
    # modular rows
    for i in range(m):
        r = [0] * D
        r[i] = n * S_K1
        rows.append(r)
    # constant / Kannan row
    r = [0] * D
    for i in range(1, m):
        r[i] = h[i] * S_K1
    r[D - 1] = S_KAN
    rows.append(r)
    return rows, (g, h, B0inv, A0)

def planted_vector_nod(sigs, n, k1_bound, k2_bound):
    m = len(sigs)
    S_K1, _S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    v = [0] * (2 * m + 1)
    for i in range(m):
        v[i] = sigs[i]['k1'] * S_K1
        v[m + i] = sigs[i]['k2'] * S_K2
    v[2 * m] = S_KAN
    return v

def recover_d_nod(M_reduced, m, n, k1_bound, k2_bound, aux, d_secret):
    """Read (k1_0, k2_0) off any row with |last| == S_KANNAN and rebuild d."""
    S_K1, _S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    _g, _h, B0inv, A0 = aux
    D = 2 * m + 1
    lam_local = recover_d_nod.lam
    for row in M_reduced:
        last = row[D - 1]
        if abs(last) != S_KAN: continue
        sign = 1 if last > 0 else -1
        if (sign * row[0]) % S_K1 or (sign * row[m]) % S_K2: continue
        k1_0 = sign * row[0] // S_K1
        k2_0 = sign * row[m] // S_K2
        d_cand = B0inv * ((k1_0 + lam_local * k2_0 - A0) % n) % n
        if d_cand == 0: continue
        if d_cand == d_secret:
            return d_cand
    return None
recover_d_nod.lam = None


# ---------------------------------------------------------------------------
# Gaussian-heuristic bookkeeping
# ---------------------------------------------------------------------------

def gh_ratio(det_log2, dim, target_norm):
    """tau = ||target|| / GH(lattice).  GH = sqrt(dim/(2*pi*e)) * det^(1/dim)."""
    gh_log2 = 0.5 * math.log2(dim / (2 * math.pi * math.e)) + det_log2 / dim
    return target_norm / (2.0 ** gh_log2)

def det_log2_old(n, m, S_K1, S_K2, S_KAN, s_d_mult=1):
    return (m * math.log2(n * S_K1) + math.log2(s_d_mult)
            + m * math.log2(S_K2) + math.log2(S_KAN))

def det_log2_new(n, m, S_K1, S_K2, S_KAN):
    return ((m - 1) * math.log2(n) + m * math.log2(S_K1)
            + m * math.log2(S_K2) + math.log2(S_KAN))

def expected_target_old(m, n, K1, K2, S_K1, S_K2, S_KAN, s_d_mult=1):
    return math.sqrt(m * (K1 * S_K1) ** 2 / 3.0
                     + (n * s_d_mult) ** 2 / 3.0
                     + m * (K2 * S_K2) ** 2 / 3.0
                     + S_KAN ** 2)

def expected_target_new(m, n, K1, K2, S_K1, S_K2, S_KAN):
    return math.sqrt(m * (K1 * S_K1) ** 2 / 3.0
                     + m * (K2 * S_K2) ** 2 / 3.0
                     + S_KAN ** 2)


# ---------------------------------------------------------------------------
# Experiment drivers
# ---------------------------------------------------------------------------

def run_old(curve, m, d_secret, K1, seed, use_bkz=False, beta=20, s_d_mult=1):
    p, b, n, lam, G = curve
    K2 = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, K1, K2, seed)
    if len(sigs) < m: return None
    M = build_glv_lattice(sigs, n, lam, K1, K2, s_d_mult=s_d_mult)
    dim = 2 * m + 2
    A = IntegerMatrix.from_matrix(M)
    if use_bkz: BKZ.reduction(A, BKZ.Param(beta))
    else: LLL.reduction(A)
    red = [[A[i][j] for j in range(dim)] for i in range(dim)]
    _, _, _, S_KAN = scales(n, K1, K2)
    ok = recover_d_old(red, m, n, S_KAN, d_secret, s_d_mult) is not None
    pn = norm(planted_vector_old(sigs, d_secret, n, K1, K2, s_d_mult))
    nz = [norm(r) for r in red if any(r)]
    return ok, pn, (min(nz) if nz else float('nan'))

def run_new(curve, m, d_secret, K1, seed, use_bkz=False, beta=20):
    p, b, n, lam, G = curve
    K2 = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, K1, K2, seed)
    if len(sigs) < m: return None
    rows, aux = build_glv_lattice_nod(sigs, n, lam, K1, K2)
    D = 2 * m + 1
    A = IntegerMatrix.from_matrix(rows)
    if use_bkz: BKZ.reduction(A, BKZ.Param(beta))
    else: LLL.reduction(A)
    red = [[A[i][j] for j in range(D)] for i in range(A.nrows)]
    recover_d_nod.lam = lam % n
    ok = recover_d_nod(red, m, n, K1, K2, aux, d_secret) is not None
    pn = norm(planted_vector_nod(sigs, n, K1, K2))
    nz = [norm(r) for r in red if any(r)]
    return ok, pn, (min(nz) if nz else float('nan'))

def rate(fn, curve, m, K1, seeds, **kw):
    wins, ratios = 0, []
    for seed in seeds:
        d = random.Random(seed + 7777).randint(1, curve[2] - 1)
        res = fn(curve, m, d, K1, seed, **kw)
        if res is None: continue
        ok, pn, sn = res
        wins += bool(ok)
        ratios.append(sn / pn if pn else float('nan'))
    return wins, len(seeds), (sum(ratios) / len(ratios) if ratios else float('nan'))


# ===========================================================================
if __name__ == "__main__":
    print("=" * 78)
    print("Thread 23 — d-eliminated Phase-2 lattice: does the K1 wall move?")
    print("=" * 78)

    SEEDS = [42, 1234, 9999, 555, 31337]

    HIST = [
        ("8-bit/199",        211,  2, 199,  106,  2,  6),
        ("12-bit/2557",      2557, 2, 2659, 1755, 8,  8),
        ("12-bit/2677 FAIL", 2677, 2, 2647, 185,  8,  10),
    ]
    hist = []
    for label, p, b, n, lam, k1, m in HIST:
        G = find_generator(p, b, n)
        assert G is not None and (lam * lam + lam + 1) % n == 0, label
        hist.append((label, (p, b, n, lam, G), k1, m))


    # ---------------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP E1: correctness of the d-eliminated construction")
    print("-" * 78)
    print("For each historical curve: (a) the planted vector lies in the new lattice")
    print("(checked by exact HNF membership), (b) the new lattice contains NO vector")
    print("of the shape c*e_j with |c| = n*S_K1 shorter than the planted vector other")
    print("than the modular rows themselves, (c) d recovers where it should.")

    def in_lattice(rows, v):
        """Exact membership test: v in rowspan_Z(rows), via HNF of [rows] vs [rows; v]."""
        from sympy import Matrix
        from sympy.matrices.normalforms import hermite_normal_form
        H1 = hermite_normal_form(Matrix(rows).T)
        H2 = hermite_normal_form(Matrix(rows + [v]).T)
        return H1 == H2

    print(f"\n{'curve':<20} {'m':>3} {'K1':>3} {'dim_old':>7} {'dim_new':>7} "
          f"{'planted in L_new':>17} {'d recovered':>12}")
    for label, curve, K1, m in hist:
        p, b, n, lam, G = curve
        K2 = math.isqrt(n) + 1
        d = random.Random(42 + 7777).randint(1, n - 1)
        sigs = gen_signatures(G, d, m, n, lam, p, K1, K2, 42)
        rows, aux = build_glv_lattice_nod(sigs, n, lam, K1, K2)
        v = planted_vector_nod(sigs, n, K1, K2)
        mem = in_lattice(rows, v)
        res = run_new(curve, m, d, K1, 42)
        print(f"{label:<20} {m:>3} {K1:>3} {2*m+2:>7} {2*m+1:>7} "
              f"{str(mem):>17} {str(bool(res and res[0])):>12}")

    # also: confirm the old lattice really does contain n*S_D*e_m and the new one
    # has no analogue (there is no d column at all).
    label, curve, K1, m = hist[1]
    p, b, n, lam, G = curve
    K2 = math.isqrt(n) + 1
    d = random.Random(42 + 7777).randint(1, n - 1)
    sigs = gen_signatures(G, d, m, n, lam, p, K1, K2, 42)
    M_old = build_glv_lattice(sigs, n, lam, K1, K2)
    triv = [0] * (2 * m + 2); triv[m] = n
    print(f"\nold lattice contains n*S_D*e_m (norm {n}): {in_lattice(M_old, triv)}")
    print(f"  ||planted_old|| = {norm(planted_vector_old(sigs, d, n, K1, K2)):.0f}"
          f"   ratio n/||planted|| = {n / norm(planted_vector_old(sigs, d, n, K1, K2)):.3f}")
    print(f"  ||planted_new|| = {norm(planted_vector_nod(sigs, n, K1, K2)):.0f}"
          f"   (d coordinate and its rival are both gone)")


    # ---------------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP E2: lattice geometry, old vs new, on the historical curves")
    print("-" * 78)
    print("tau = ||planted|| / GH(det,dim).  Unique-SVP via LLL needs roughly")
    print("tau < 1/1.02^dim; tau > 1 means the planted vector is not even expected")
    print("to be lambda_1.  sv/pv = ||shortest reduced|| / ||planted|| (measured).")

    print(f"\n{'curve':<20} {'m':>3} {'K1':>3} | {'tau_old':>8} {'sv/pv_old':>10} "
          f"{'win_old':>8} | {'tau_new':>8} {'sv/pv_new':>10} {'win_new':>8}")
    for label, curve, K1, m in hist:
        p, b, n, lam, G = curve
        K2 = math.isqrt(n) + 1
        S_K1, S_D, S_K2, S_KAN = scales(n, K1, K2)
        t_old = gh_ratio(det_log2_old(n, m, S_K1, S_K2, S_KAN), 2 * m + 2,
                         expected_target_old(m, n, K1, K2, S_K1, S_K2, S_KAN))
        t_new = gh_ratio(det_log2_new(n, m, S_K1, S_K2, S_KAN), 2 * m + 1,
                         expected_target_new(m, n, K1, K2, S_K1, S_K2, S_KAN))
        wo, to, ro = rate(run_old, curve, m, K1, SEEDS)
        wn, tn, rn = rate(run_new, curve, m, K1, SEEDS)
        print(f"{label:<20} {m:>3} {K1:>3} | {t_old:>8.3f} {ro:>10.3f} "
              f"{f'{wo}/{to}':>8} | {t_new:>8.3f} {rn:>10.3f} {f'{wn}/{tn}':>8}")


    # ---------------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP E3: THE FALSIFIER — T4's K1 grid, old lattice vs new lattice")
    print("-" * 78)
    print("Same two 12-bit curves, same K2 = isqrt(n)+1, same 5 seeds, m = 12 for")
    print("both so the two arms are directly comparable.  T4 (2026-07-29) found the")
    print("wall at K1 ~ 12-16 (lam*=0.34) and K1 ~ 4-6 (lam*=0.07) with the old")
    print("lattice.  If the new lattice moves the wall outward, Thread 23 is a real")
    print("improvement; if not, Phase 2 is at its lattice-geometric ceiling.")

    K1_GRID = [2, 3, 4, 6, 8, 12, 16, 24, 32]
    M_GRID = 12
    grid_rows = []
    for label, curve, _K1, _m in hist[1:]:
        p, b, n, lam, G = curve
        K2 = math.isqrt(n) + 1
        print(f"\ncurve {label}   n={n}  lam*={lam_star(lam, n):.4f}  K2={K2}  m={M_GRID}")
        print(f"  {'K1':>4} | {'old wins':>9} {'tau_old':>8} | {'new wins':>9} "
              f"{'tau_new':>8} | {'new BKZ20':>10}")
        for K1 in K1_GRID:
            S_K1, S_D, S_K2, S_KAN = scales(n, K1, K2)
            t_old = gh_ratio(det_log2_old(n, M_GRID, S_K1, S_K2, S_KAN), 2 * M_GRID + 2,
                             expected_target_old(M_GRID, n, K1, K2, S_K1, S_K2, S_KAN))
            t_new = gh_ratio(det_log2_new(n, M_GRID, S_K1, S_K2, S_KAN), 2 * M_GRID + 1,
                             expected_target_new(M_GRID, n, K1, K2, S_K1, S_K2, S_KAN))
            wo, to, _ = rate(run_old, curve, M_GRID, K1, SEEDS)
            wn, tn, _ = rate(run_new, curve, M_GRID, K1, SEEDS)
            wb, tb, _ = rate(run_new, curve, M_GRID, K1, SEEDS, use_bkz=True, beta=20)
            print(f"  {K1:>4} | {f'{wo}/{to}':>9} {t_old:>8.3f} | {f'{wn}/{tn}':>9} "
                  f"{t_new:>8.3f} | {f'{wb}/{tb}':>10}")
            grid_rows.append((label, K1, wo, wn, wb, t_old, t_new))


    # ---------------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP E4: control — sweep S_D in the OLD lattice")
    print("-" * 78)
    print("The 2026-07-29 diagnosis was that the trivial vector n*S_D*e_m (norm")
    print("n*S_D) crowds out the planted vector (norm ~ n*sqrt(2m/3+4/3)).  Both")
    print("scale linearly in S_D *only if d ~ n*S_D dominates the planted norm*; at")
    print("moderate S_D the planted norm grows more slowly, so if that vector were")
    print("the binding constraint, S_D in 2..8 should help.  It does not (see below),")
    print("which is the same conclusion E3 reaches by a different route.")

    label, curve, _K1, _m = hist[2]
    p, b, n, lam, G = curve
    K2 = math.isqrt(n) + 1
    print(f"\ncurve {label}  n={n}  m={M_GRID}  K2={K2}")
    print(f"  {'S_D':>4} | " + " ".join(f"K1={k:<2}".rjust(7) for k in [2, 4, 6, 8, 12]))
    for s_d in [1, 2, 4, 8, 16]:
        cells = []
        for K1 in [2, 4, 6, 8, 12]:
            w, t, _ = rate(run_old, curve, M_GRID, K1, SEEDS, s_d_mult=s_d)
            cells.append(f"{w}/{t}".rjust(7))
        print(f"  {s_d:>4} | " + " ".join(cells))


    # ---------------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP E5: scaling check on independent 17-bit curves")
    print("-" * 78)
    print("Three fresh j=0 GLV curves near 2^17.  Wall = largest K1 with >= 4/5.")

    def find_curves(lo, hi, want=3):
        out = []
        p = int(sympy.nextprime(lo))
        while p < hi and len(out) < want:
            if p % 3 == 1:
                eis = eisenstein_decompose(p)
                if eis is not None:
                    for t in j0_traces(*eis):
                        nc = p + 1 - t
                        if nc < 2 or nc % 3 != 1 or not sympy.isprime(nc): continue
                        r = glv_roots(nc)
                        if r is None: continue
                        cur = build_curve(p, nc)
                        if cur is None: continue
                        out.append((p, cur[1], nc, r[0], cur[4]))
                        break
            p = int(sympy.nextprime(p))
        return out[:want]

    def wall(fn, curve, m, seeds, **kw):
        """Largest K1 in the grid reaching >= 4/5, scanning upward and stopping
        after two consecutive grid points below threshold."""
        best, misses = 0, 0
        for K1 in [2, 3, 4, 6, 8, 12, 16, 24, 32, 48]:
            w, _t, _ = rate(fn, curve, m, K1, seeds, **kw)
            if w >= 4:
                best, misses = K1, 0
            else:
                misses += 1
                if misses >= 2:
                    break
        return best

    c17 = find_curves(2 ** 17, 2 ** 17 + 6000, want=3)
    print(f"\n{'p':>8} {'n':>8} {'lam*':>7} {'K2':>5} {'m':>3} | "
          f"{'wall_old':>9} {'wall_new':>9} {'gain':>6}")
    for cur in c17:
        p, b, n, lam, G = cur
        K2 = math.isqrt(n) + 1
        wo = wall(run_old, cur, M_GRID, SEEDS)
        wn = wall(run_new, cur, M_GRID, SEEDS)
        gain = (wn / wo) if wo else float('nan')
        print(f"{p:>8} {n:>8} {lam_star(lam, n):>7.4f} {K2:>5} {M_GRID:>3} | "
              f"{wo:>9} {wn:>9} {gain:>6.2f}")

    print("\n" + "=" * 78)
    print("done")
    print("=" * 78)
