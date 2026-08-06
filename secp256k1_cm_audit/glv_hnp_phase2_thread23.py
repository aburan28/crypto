"""
Thread 23 — reformulate the GLV-HNP Phase-2 lattice so the planted vector is lambda_1.

Motivation (RESEARCH_AUTOLAB_LOG.md 2026-07-29, T5): in the Phase-2 lattice
L_A (`glv_hnp_phase2_20bit.py:262 build_glv_lattice`) the shortest vector is
ALWAYS the trivial `n*S_D*e_m`, because n*row_d - sum_i B_i*row_i = n*S_D*e_m
and ||n*S_D*e_m|| = n*S_D while ||v_planted|| ~ n*S_D*sqrt(2m/3+4/3).
Both scale linearly in S_D, so no rescaling removes it. Recovery in L_A is
therefore a BDD/coset condition, not an SVP condition.

This script:
  L_B  = d eliminated algebraically (pivot on signature 0), dim 2m+1.
         The trivial vector is gone by construction.
  L_Bc = L_B with the unknowns centred (k1 -> k1 - K1/2 etc), which shrinks
         ||v_planted|| by ~2x.

Experiments:
  T23-1  correctness + is v_planted actually lambda_1 now?   (the log's falsifier)
  T23-2  T4 grid replay: does the K1 wall move?  L_A vs L_B vs L_Bc
  T23-3  Gaussian-heuristic threshold  eff* = (3/2 pi e) * n^(-1/m)  vs measurement
  T23-4  sharp falsifier: 12-bit/2677 at K1=8 (eff=0.157) needs m~70 by the
         formula.  Sweep m in {8,...,128} and see whether it ever recovers.

Run: python3 glv_hnp_phase2_thread23.py
"""

import math
import random
import sys

import sympy
from fpylll import IntegerMatrix, LLL, BKZ

# ---------------------------------------------------------------------------
# EC arithmetic (verbatim from glv_hnp_phase2_20bit.py so results are comparable)
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
# j=0 CM machinery
# ---------------------------------------------------------------------------

def eisenstein_decompose(p):
    for a in range(1, 2 * math.isqrt(p // 3) + 3):
        disc = 4 * p - 3 * a * a
        if disc < 0:
            break
        s = math.isqrt(disc)
        if s * s != disc:
            continue
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
    if sq is None:
        return None
    inv2 = modinv(2, n)
    r1 = (n - 1 + sq) * inv2 % n
    r2 = (n - 1 + (n - sq)) * inv2 % n
    if (r1 * r1 + r1 + 1) % n != 0:
        r1, r2 = r2, r1
    if (r1 * r1 + r1 + 1) % n != 0:
        return None
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
        if y is None or y == 0:
            continue
        if ec_mul((x, y), n, p) is None:
            G = find_generator(p, b, n, seed=seed)
            if G is not None:
                return (p, b, n, None, G)
    return None

def find_curves_at_bits(bits, count, seed_offset=0):
    """j=0 GLV curves with n prime, n = 1 mod 3, of the requested bit length."""
    out = []
    p = int(sympy.nextprime(2 ** (bits - 1) + seed_offset))
    while len(out) < count and p < 2 ** bits:
        if p % 3 == 1:
            eis = eisenstein_decompose(p)
            if eis is not None:
                for t in j0_traces(*eis):
                    n_cand = p + 1 - t
                    if n_cand < 8 or n_cand % 3 != 1 or n_cand.bit_length() != bits:
                        continue
                    if not sympy.isprime(n_cand):
                        continue
                    roots = glv_roots(n_cand)
                    if roots is None:
                        continue
                    cur = build_curve(p, n_cand)
                    if cur is None:
                        continue
                    out.append((p, cur[1], n_cand, roots[0], cur[4]))
                    break
        p = int(sympy.nextprime(p))
    return out

# ---------------------------------------------------------------------------
# Signatures
# ---------------------------------------------------------------------------

def gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed=42):
    rng = random.Random(seed)
    sigs = []
    attempts = 0
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
        A = h * s_inv % n
        B = r * s_inv % n
        assert (A + B * d_secret) % n == k_full
        sigs.append({'A': A, 'B': B, 'k1': k1, 'k2': k2})
    return sigs

def scales(n, k1_bound, k2_bound):
    return n // k1_bound, 1, max(1, n // k2_bound), n

# ---------------------------------------------------------------------------
# L_A : the existing Phase-2 lattice (baseline, dim 2m+2)
# ---------------------------------------------------------------------------

def build_LA(sigs, n, lam, k1_bound, k2_bound):
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
    return M, (S_K1, S_D, S_K2, S_KAN)

def planted_LA(sigs, d, n, k1_bound, k2_bound):
    m = len(sigs)
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    v = [0] * (2 * m + 2)
    for i in range(m):
        v[i] = sigs[i]['k1'] * S_K1
        v[m + 1 + i] = sigs[i]['k2'] * S_K2
    v[m] = d * S_D
    v[2 * m + 1] = S_KAN
    return v

def recover_LA(rows, m, n, S_KAN, d_secret):
    for row in rows:
        last = row[2 * m + 1]
        if abs(last) != S_KAN:
            continue
        sign = 1 if last > 0 else -1
        if (sign * row[m]) % n == d_secret:
            return True
    return False

# ---------------------------------------------------------------------------
# L_B : d eliminated by pivoting on signature 0 (dim 2m+1)
# ---------------------------------------------------------------------------
#
#   A_i + B_i d = k1_i + lam k2_i          (mod n),  i = 0..m-1
#   t_i = B_i / B_0,  u_i = A_i - t_i A_0                 (i = 1..m-1)
#   =>  k1_i = u_i + t_i k1_0 + t_i lam k2_0 - lam k2_i   (mod n)
#
# Free unknowns: k1_0, k2_0..k2_{m-1}.  Witnesses: k1_1..k1_{m-1}.
# Columns: [k1_1..k1_{m-1}] [k1_0] [k2_0..k2_{m-1}] [Kannan]   -> dim 2m+1.
#
# `centre=True` replaces k1_j = c1 + x_j, k2_j = c2 + y_j with c1=K1//2,
# c2=K2//2, which shrinks every planted coordinate by ~2x.

def build_LB(sigs, n, lam, k1_bound, k2_bound, centre=False, kan=None):
    m = len(sigs)
    assert m >= 2
    dim = 2 * m + 1
    S_K1, _, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    if centre:
        S_KAN = S_KAN // 2 if kan is None else kan
    elif kan is not None:
        S_KAN = kan
    c1 = k1_bound // 2 if centre else 0
    c2 = k2_bound // 2 if centre else 0

    B0inv = modinv(sigs[0]['B'], n)
    t = [0] * m
    u = [0] * m
    for i in range(1, m):
        t[i] = sigs[i]['B'] * B0inv % n
        u[i] = (sigs[i]['A'] - t[i] * sigs[0]['A']) % n
        if centre:
            u[i] = (u[i] - c1 + t[i] * c1 + t[i] * lam * c2 - lam * c2) % n

    W = m - 1                       # number of witness columns
    COL_K10 = W                     # column of k1_0
    COL_K2 = W + 1                  # first k2 column
    COL_KAN = dim - 1

    M = [[0] * dim for _ in range(dim)]
    r = 0
    for j in range(W):              # reduction rows on the witness columns
        M[r][j] = n * S_K1; r += 1
    for j in range(W):              # k1_0 row
        M[r][j] = t[j + 1] * S_K1
    M[r][COL_K10] = S_K1; r += 1
    for j in range(W):              # k2_0 row
        M[r][j] = t[j + 1] * lam % n * S_K1
    M[r][COL_K2] = S_K2; r += 1
    for i in range(1, m):           # k2_i rows  (i = 1..m-1)
        M[r][i - 1] = -lam * S_K1
        M[r][COL_K2 + i] = S_K2; r += 1
    for j in range(W):              # Kannan row
        M[r][j] = u[j + 1] * S_K1
    M[r][COL_KAN] = S_KAN; r += 1
    assert r == dim, (r, dim)
    return M, (S_K1, S_K2, S_KAN, c1, c2)

def planted_LB(sigs, n, k1_bound, k2_bound, centre=False, kan=None):
    m = len(sigs)
    dim = 2 * m + 1
    S_K1, _, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    if centre:
        S_KAN = S_KAN // 2 if kan is None else kan
    elif kan is not None:
        S_KAN = kan
    c1 = k1_bound // 2 if centre else 0
    c2 = k2_bound // 2 if centre else 0
    v = [0] * dim
    for i in range(1, m):
        v[i - 1] = (sigs[i]['k1'] - c1) * S_K1
    v[m - 1] = (sigs[0]['k1'] - c1) * S_K1
    for j in range(m):
        v[m + j] = (sigs[j]['k2'] - c2) * S_K2
    v[dim - 1] = S_KAN
    return v

def recover_LB(rows, sigs, n, lam, k1_bound, k2_bound, d_secret,
               centre=False, kan=None):
    """Read k1_0,k2_0 off any row with |last| = S_KAN and rebuild d."""
    m = len(sigs)
    dim = 2 * m + 1
    S_K1, _, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    if centre:
        S_KAN = S_KAN // 2 if kan is None else kan
    elif kan is not None:
        S_KAN = kan
    c1 = k1_bound // 2 if centre else 0
    c2 = k2_bound // 2 if centre else 0
    B0inv = modinv(sigs[0]['B'], n)
    for row in rows:
        last = row[dim - 1]
        if abs(last) != S_KAN:
            continue
        sg = 1 if last > 0 else -1
        a, b = sg * row[m - 1], sg * row[m]
        if a % S_K1 or b % S_K2:
            continue
        k10 = a // S_K1 + c1
        k20 = b // S_K2 + c2
        k_full = (k10 + lam * k20) % n
        if (k_full - sigs[0]['A']) % n == 0 and sigs[0]['B'] % n == 0:
            continue
        d_cand = (k_full - sigs[0]['A']) * B0inv % n
        if d_cand == d_secret:
            return True
    return False

# ---------------------------------------------------------------------------
# Lattice utilities
# ---------------------------------------------------------------------------

def norm(v):
    return math.sqrt(sum(x * x for x in v))

def reduce_matrix(M, beta=0):
    dim = len(M)
    A = IntegerMatrix.from_matrix(M)
    if beta:
        BKZ.reduction(A, BKZ.Param(beta))
    else:
        LLL.reduction(A)
    return [[A[i][j] for j in range(dim)] for i in range(dim)]

def log_det(M):
    """log|det| via exact integer determinant on the (triangular-ish) basis."""
    d = sympy.Matrix(M).det()
    return math.log(abs(d)) if d else float('-inf')

def gaussian_heuristic(logdet, dim):
    return math.exp(0.5 * math.log(dim / (2 * math.pi * math.e)) + logdet / dim)

def shortest_row(rows):
    best, bn = None, None
    for row in rows:
        if all(x == 0 for x in row):
            continue
        nv = norm(row)
        if bn is None or nv < bn:
            best, bn = row, nv
    return best, bn

# ---------------------------------------------------------------------------
# One trial
# ---------------------------------------------------------------------------

def trial(curve, m, k1_bound, seed, variant, beta=0, want_diag=False):
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    d_secret = (seed * 7919 + 12345) % (n - 2) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return None
    if variant == 'A':
        M, sc = build_LA(sigs, n, lam, k1_bound, k2_bound)
        pv = planted_LA(sigs, d_secret, n, k1_bound, k2_bound)
        rows = reduce_matrix(M, beta)
        ok = recover_LA(rows, m, n, sc[3], d_secret)
    else:
        centre = (variant == 'Bc')
        M, sc = build_LB(sigs, n, lam, k1_bound, k2_bound, centre=centre)
        pv = planted_LB(sigs, n, k1_bound, k2_bound, centre=centre)
        rows = reduce_matrix(M, beta)
        ok = recover_LB(rows, sigs, n, lam, k1_bound, k2_bound, d_secret,
                        centre=centre)
    if not want_diag:
        return ok
    sv, svn = shortest_row(rows)
    ld = log_det(M)
    gh = gaussian_heuristic(ld, len(M))
    return {'ok': ok, 'pv_norm': norm(pv), 'sv_norm': svn,
            'sv_over_pv': svn / norm(pv), 'gh': gh,
            'pv_over_gh': norm(pv) / gh, 'dim': len(M)}

def rate(curve, m, k1_bound, seeds, variant, beta=0):
    hits = 0
    for s in seeds:
        r = trial(curve, m, k1_bound, s, variant, beta)
        if r:
            hits += 1
    return hits

def k1_threshold(curve, m, seeds, variant, beta=0, cap=1 << 20):
    """Largest K1 with >= 50% recovery, by doubling then bisection.

    Returns None if K1=1 already fails, and -cap if the threshold exceeds
    `cap` (never silently reports 0 -- that was the bug in the first run,
    see RESEARCH_AUTOLAB_LOG.md 2026-08-06).
    """
    half = (len(seeds) + 1) // 2
    if rate(curve, m, 1, seeds, variant, beta) < half:
        return None
    lo = 1
    hi = 2
    while hi <= cap and rate(curve, m, hi, seeds, variant, beta) >= half:
        lo, hi = hi, hi * 2
    if hi > cap:
        return -cap
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if rate(curve, m, mid, seeds, variant, beta) >= half:
            lo = mid
        else:
            hi = mid
    return lo

# ---------------------------------------------------------------------------

SEEDS = [42, 1234, 9999, 555, 31337]

# The two 12-bit curves carried through the whole GLV-HNP thread
# (log 2026-06-15 .. 2026-07-29).  (label, p, b, n, lam)
C_2557 = ("12-bit/2557", 2557, 2, 2659, 1755)
C_2677 = ("12-bit/2677", 2677, 2, 2647, 185)

def load_known(spec):
    label, p, b, n, lam = spec
    G = find_generator(p, b, n)
    assert G is not None, label
    roots = glv_roots(n)
    assert lam in roots, (label, lam, roots)
    return label, (p, b, n, lam, G)


def main():
    print("=" * 78)
    print("Thread 23 — is the Phase-2 wall the Gaussian-heuristic uSVP bound?")
    print("=" * 78)

    lab1, c1 = load_known(C_2557)
    lab2, c2 = load_known(C_2677)
    for lab, c in ((lab1, c1), (lab2, c2)):
        n = c[2]
        print(f"  {lab}: p={c[0]} n={n} lam={c[3]} lam*={lam_star(c[3], n):.4f} "
              f"K2={math.isqrt(n)+1}")

    # -----------------------------------------------------------------
    print("\n" + "=" * 78)
    print("T23-1  L_B correctness, and is v_planted lambda_1 now?")
    print("=" * 78)
    print(f"{'curve':<12} {'var':<4} {'dim':>4} {'K1':>3} {'ok':>3} "
          f"{'sv/pv':>7} {'pv/GH':>7}")
    for lab, c in ((lab1, c1), (lab2, c2)):
        for K1 in (2, 4, 8):
            for var in ('A', 'B', 'Bc'):
                r = trial(c, 12, K1, 42, var, want_diag=True)
                print(f"{lab:<12} {var:<4} {r['dim']:>4} {K1:>3} "
                      f"{str(r['ok']):>5} {r['sv_over_pv']:>7.3f} "
                      f"{r['pv_over_gh']:>7.3f}")

    # -----------------------------------------------------------------
    print("\n" + "=" * 78)
    print("T23-2  T4 grid replay — does the K1 wall move?  (m=12, 5 seeds)")
    print("=" * 78)
    K1S = [2, 3, 4, 6, 8, 12, 16, 24]
    print(f"{'curve':<12} {'var':<4} " + " ".join(f"{k:>4}" for k in K1S))
    for lab, c in ((lab1, c1), (lab2, c2)):
        for var in ('A', 'B', 'Bc'):
            row = [rate(c, 12, K1, SEEDS, var) for K1 in K1S]
            print(f"{lab:<12} {var:<4} " + " ".join(f"{v:>4}" for v in row))
    print("\n(eff = K1*K2/n, K2 = isqrt(n)+1)")
    for lab, c in ((lab1, c1), (lab2, c2)):
        n = c[2]; K2 = math.isqrt(n) + 1
        print(f"{lab:<12} eff  " + " ".join(f"{K1*K2/n:>4.2f}" for K1 in K1S))

    # -----------------------------------------------------------------
    print("\n" + "=" * 78)
    print("T23-3  Gaussian-heuristic threshold   eff* = (3/2 pi e) * n^(-1/m)")
    print("=" * 78)
    GH_CONST = 3.0 / (2 * math.pi * math.e)
    print(f"  asymptotic constant 3/(2 pi e) = {GH_CONST:.4f}")

    print("  counting bound   eff_count = n^(-1/m)  (unique-solution bound)")
    print("  reported C = eff_measured / n^(-1/m); C -> 1 means the attack")
    print("  sits exactly on the counting bound.")

    print("\n(a) sweep m at fixed curve (12-bit/2677), variant Bc, 5 seeds")
    print(f"{'m':>4} {'dim':>4} {'eff*_GH':>9} {'eff_count':>10} "
          f"{'K1@50%':>7} {'eff_meas':>9} {'C':>6}")
    n2 = c2[2]; K2_2 = math.isqrt(n2) + 1
    for m in (4, 6, 8, 12, 16, 24, 32):
        base = n2 ** (-1.0 / m)
        K1s = k1_threshold(c2, m, SEEDS, 'Bc')
        eff_meas = K1s * K2_2 / n2 if K1s and K1s > 0 else float('nan')
        print(f"{m:>4} {2*m+1:>4} {GH_CONST*base:>9.4f} {base:>10.4f} "
              f"{str(K1s):>7} {eff_meas:>9.4f} {eff_meas/base:>6.3f}")
        sys.stdout.flush()

    print("\n(b) sweep n bit-length at fixed m=12, variant Bc, 5 seeds")
    print(f"{'bits':>5} {'n':>10} {'K2':>5} {'eff*_GH':>9} {'eff_count':>10} "
          f"{'K1@50%':>7} {'eff_meas':>9} {'C':>6}")
    for bits in (10, 12, 14, 16, 18, 20):
        cs = find_curves_at_bits(bits, 1)
        if not cs:
            print(f"{bits:>5}  (no curve found)")
            continue
        cur = cs[0]
        n = cur[2]; K2 = math.isqrt(n) + 1
        base = n ** (-1.0 / 12)
        K1s = k1_threshold(cur, 12, SEEDS, 'Bc')
        eff_meas = K1s * K2 / n if K1s and K1s > 0 else float('nan')
        print(f"{bits:>5} {n:>10} {K2:>5} {GH_CONST*base:>9.4f} {base:>10.4f} "
              f"{str(K1s):>7} {eff_meas:>9.4f} {eff_meas/base:>6.3f}")
        sys.stdout.flush()

    # -----------------------------------------------------------------
    print("\n" + "=" * 78)
    print("T23-4  sharp falsifier: 12-bit/2677 at K1=8 (eff=0.157)")
    print("=" * 78)
    eff_target = 8 * K2_2 / n2
    m_needed = math.log(n2) / math.log(GH_CONST / eff_target) if eff_target < GH_CONST else float('inf')
    print(f"  eff = 8*{K2_2}/{n2} = {eff_target:.4f}")
    print(f"  formula eff* = {GH_CONST:.4f} * n^(-1/m)  =>  m_needed = {m_needed:.1f}")
    print(f"  2026-07-29 T4b measured 0,0,1,0,1 of 5 at m = 8,12,16,24,32.")
    print(f"\n{'m':>5} {'dim':>5} {'eff*_pred':>10} {'A':>4} {'B':>4} {'Bc':>4} "
          f"{'Bc/BKZ20':>9}")
    for m in (8, 16, 32, 48, 64, 96, 128):
        eff_pred = GH_CONST * n2 ** (-1.0 / m)
        rA = rate(c2, m, 8, SEEDS, 'A')
        rB = rate(c2, m, 8, SEEDS, 'B')
        rC = rate(c2, m, 8, SEEDS, 'Bc')
        rK = rate(c2, m, 8, SEEDS[:3], 'Bc', beta=20) if m <= 64 else -1
        print(f"{m:>5} {2*m+1:>5} {eff_pred:>10.4f} {rA:>4} {rB:>4} {rC:>4} "
              f"{(str(rK)+'/3' if rK >= 0 else 'skip'):>9}")
        sys.stdout.flush()

    # -----------------------------------------------------------------
    print("\n" + "=" * 78)
    print("T23-5  A vs B vs Bc threshold across bit lengths (m=12, 5 seeds)")
    print("=" * 78)
    print(f"{'bits':>5} {'n':>10} {'K2':>5} {'A':>6} {'B':>6} {'Bc':>6} "
          f"{'effA':>7} {'effBc':>7} {'gain':>6}")
    for bits in (10, 12, 14, 16, 18, 20):
        cs = find_curves_at_bits(bits, 1)
        if not cs:
            continue
        cur = cs[0]
        n = cur[2]; K2 = math.isqrt(n) + 1
        r = {v: k1_threshold(cur, 12, SEEDS, v) for v in ('A', 'B', 'Bc')}
        eA = (r['A'] or 0) * K2 / n
        eC = (r['Bc'] or 0) * K2 / n
        g = (r['Bc'] / r['A']) if r['A'] else float('nan')
        print(f"{bits:>5} {n:>10} {K2:>5} {str(r['A']):>6} {str(r['B']):>6} "
              f"{str(r['Bc']):>6} {eA:>7.4f} {eC:>7.4f} {g:>5.1f}x")
        sys.stdout.flush()

    print("\ndone.")


if __name__ == '__main__':
    main()
