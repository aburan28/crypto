"""
Thread 23 — reformulate the GLV-HNP Phase-2 lattice so the planted vector is
lambda_1.

Motivation (RESEARCH_AUTOLAB_LOG.md 2026-07-29, finding T5): in the Phase-2
lattice built by `glv_hnp_phase2_20bit.py:262` (`build_glv_lattice`) the
shortest vector after LLL is ALWAYS the trivial vector `n*S_D*e_m`, never the
planted vector.  Cause: the secret `d` is unbounded in [0,n), so its column
carries scale S_D = 1, and the lattice therefore contains
    n * row_d  -  sum_i B_i * row_i   =   n * S_D * e_m       (norm n)
while the planted vector has norm ~ n*sqrt(2m/3 + 4/3).  The Kannan embedding
is therefore provably not solving SVP: it only ever works when LLL happens to
surface the planted vector in some *other* basis row.

V2 reformulation: eliminate `d` with signature 0 (standard HNP substitution),
so that EVERY remaining unknown is bounded and every column can be scaled up.

    k_i = A_i + B_i*d,  k_0 = A_0 + B_0*d
    => d = (k_0 - A_0) * B_0^{-1}  (mod n)
    => k_i = D_i + C_i * k_0       (mod n),  C_i = B_i*B_0^{-1}, D_i = A_i - C_i*A_0

With k_i = k1_i + lam*k2_i this gives, for i = 1..m-1,

    k1_i  =  D_i + C_i*k1_0 + lam*C_i*k2_0 - lam*k2_i   (mod n)

Free unknowns : k1_0 (<K1), k2_0..k2_{m-1} (<K2)          -> m+1
Derived, bounded: k1_1..k1_{m-1} (<K1)                    -> m-1
Lattice dimension 2m+1 (V1 was 2m+2), no unbounded coordinate, so the trivial
vector is n*S_K1*e_{k1_0} with norm n*S_K1 = n^2/K1 instead of n.

Falsifier (stated 2026-07-29): "if sv/pv rises above 1 after the reformulation
AND the K1 wall in T4 moves outward on the lam*=0.07 curve (currently K1~4-6),
the reformulation is a real improvement; if the wall stays at K1~4-6, then the
wall is information-theoretic and Phase 2 is at its ceiling."

Run: python3 glv_hnp_phase2_v2_lattice.py
"""

import math
import random
import sympy
from fpylll import IntegerMatrix, LLL, BKZ

# ---------------------------------------------------------------------------
# EC arithmetic / CM helpers  (verbatim from glv_hnp_phase2_20bit.py:22-160)
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

def find_generator(p, b, n):
    rng = random.Random(12345)
    for _ in range(10000):
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
        if disc < 0:
            break
        s = math.isqrt(disc)
        if s * s != disc:
            continue
        for num in [a + s, a - s]:
            if num % 2 == 0:
                b = num // 2
                if b >= 0 and a * a - a * b + b * b == p:
                    return (a, b)
    return None

def j0_traces(a, b):
    return [2*a - b, -2*a + b, -(a + b), a + b, 2*b - a, a - 2*b]

def glv_eigenvalue(n):
    neg3 = (n - 3) % n
    sq = tonelli_shanks(neg3, n)
    if sq is None:
        return None
    inv2 = modinv(2, n)
    r1 = (n - 1 + sq) * inv2 % n
    r2 = (n - 1 + (n - sq)) * inv2 % n
    if (r1 * r1 + r1 + 1) % n != 0:
        r1, r2 = r2, r1
    if (r1 * r1 + r1 + 1) % n != 0:
        return None
    return min(r1, r2)

def lam_star(lam, n):
    return min(lam, n - lam) / n

def build_curve(p, n):
    """Find twist parameter b with #E(F_p: y^2=x^3+b) = n, plus a generator."""
    for b_try in range(1, 200):
        G = find_generator(p, b_try, n)
        if G is not None:
            return (p, b_try, n, glv_eigenvalue(n), G)
    return None

# ---------------------------------------------------------------------------
# Signatures  (verbatim from glv_hnp_phase2_20bit.py:236)
# ---------------------------------------------------------------------------

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
        sigs.append({'A': A, 'B': B, 'k1': k1, 'k2': k2, 'k_full': k_full})
    return sigs

# ---------------------------------------------------------------------------
# V1 lattice — verbatim copy of build_glv_lattice (glv_hnp_phase2_20bit.py:262)
# ---------------------------------------------------------------------------

def build_v1(sigs, n, lam, k1_bound, k2_bound):
    m = len(sigs)
    dim = 2 * m + 2
    S_K1 = n // k1_bound
    S_D = 1
    S_K2 = max(1, n // k2_bound)
    S_KANNAN = n

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
    return M, dict(S_K1=S_K1, S_D=S_D, S_K2=S_K2, S_KANNAN=S_KANNAN, dim=dim, m=m)


def planted_v1(sigs, n, lam, sc):
    """The vector LLL is supposed to find, built from the true (k1_i,k2_i,d)."""
    m = sc['m']; dim = sc['dim']
    v = [0] * dim
    for i in range(m):
        v[i] = sigs[i]['k1'] * sc['S_K1']
    # d is recovered from sig 0:  d = (k_full - A)/B
    d = (sigs[0]['k_full'] - sigs[0]['A']) * modinv(sigs[0]['B'], n) % n
    v[m] = d * sc['S_D']
    for i in range(m):
        v[m + 1 + i] = sigs[i]['k2'] * sc['S_K2']
    v[dim - 1] = sc['S_KANNAN']
    return v


def recover_v1(rows, sigs, n, lam, sc, d_secret):
    m = sc['m']; dim = sc['dim']
    for r_idx, row in enumerate(rows):
        last = row[dim - 1]
        if abs(last) != sc['S_KANNAN']:
            continue
        sign = 1 if last > 0 else -1
        d_cand = (sign * row[m]) % n
        if d_cand == 0:
            continue
        if d_cand == d_secret:
            return d_cand, r_idx
    return None, None

# ---------------------------------------------------------------------------
# V2 lattice — d eliminated via signature 0.  dim = 2m+1.
# ---------------------------------------------------------------------------
#   col 0 .. m-2   : derived k1_i, i=1..m-1        scale S_K1
#   col m-1        : k1_0                          scale S_K1
#   col m          : k2_0                          scale S_K2
#   col m+1..2m-1  : k2_1..k2_{m-1}                scale S_K2
#   col 2m         : Kannan                        scale S_KANNAN
#
#   row 0 .. m-2   : n*S_K1 * e_col                (modulus rows)
#   row m-1        : k1_0 unknown
#   row m          : k2_0 unknown
#   row m+1..2m-1  : k2_j unknown, j=1..m-1
#   row 2m         : Kannan / constant row

def build_v2_plain(sigs, n, lam, k1_bound, k2_bound):
    m = len(sigs)
    assert m >= 2
    dim = 2 * m + 1
    S_K1 = n // k1_bound
    S_K2 = max(1, n // k2_bound)
    S_KANNAN = n

    B0inv = modinv(sigs[0]['B'], n)
    C = [None] * m
    D = [None] * m
    for i in range(1, m):
        C[i] = sigs[i]['B'] * B0inv % n
        D[i] = (sigs[i]['A'] - C[i] * sigs[0]['A']) % n

    M = [[0] * dim for _ in range(dim)]
    # modulus rows for the m-1 derived k1 columns
    for j in range(m - 1):
        M[j][j] = n * S_K1
    # k1_0 row
    for i in range(1, m):
        M[m - 1][i - 1] = C[i] * S_K1
    M[m - 1][m - 1] = S_K1
    # k2_0 row
    for i in range(1, m):
        M[m][i - 1] = lam * C[i] % n * S_K1
    M[m][m] = S_K2
    # k2_j rows, j = 1..m-1
    for j in range(1, m):
        M[m + j][j - 1] = -lam * S_K1
        M[m + j][m + j] = S_K2
    # Kannan row
    for i in range(1, m):
        M[2 * m][i - 1] = D[i] * S_K1
    M[2 * m][2 * m] = S_KANNAN

    return M, dict(S_K1=S_K1, S_K2=S_K2, S_KANNAN=S_KANNAN, dim=dim, m=m,
                   C=C, D=D, B0inv=B0inv)


# EXP I swaps this name to build_v3 to test a pre-reduced starting basis.
build_v2 = build_v2_plain


def planted_v2(sigs, n, lam, sc):
    m = sc['m']; dim = sc['dim']
    v = [0] * dim
    for i in range(1, m):
        v[i - 1] = sigs[i]['k1'] * sc['S_K1']
    v[m - 1] = sigs[0]['k1'] * sc['S_K1']
    v[m] = sigs[0]['k2'] * sc['S_K2']
    for j in range(1, m):
        v[m + j] = sigs[j]['k2'] * sc['S_K2']
    v[2 * m] = sc['S_KANNAN']
    return v


def recover_v2(rows, sigs, n, lam, sc, d_secret):
    """Read (k1_0,k2_0) off any row with |Kannan coord| = S_KANNAN, then
    d = (k1_0 + lam*k2_0 - A_0) * B_0^{-1} mod n."""
    m = sc['m']; dim = sc['dim']
    for r_idx, row in enumerate(rows):
        last = row[dim - 1]
        if abs(last) != sc['S_KANNAN']:
            continue
        sign = 1 if last > 0 else -1
        x1 = sign * row[m - 1]
        x2 = sign * row[m]
        if x1 % sc['S_K1'] or x2 % sc['S_K2']:
            continue
        k1_0 = x1 // sc['S_K1']
        k2_0 = x2 // sc['S_K2']
        k0 = (k1_0 + lam * k2_0) % n
        d_cand = (k0 - sigs[0]['A']) * sc['B0inv'] % n
        if d_cand == 0:
            continue
        if d_cand == d_secret:
            return d_cand, r_idx
    return None, None

# ---------------------------------------------------------------------------
# Instrumented single run
# ---------------------------------------------------------------------------

def norm(v):
    return math.sqrt(sum(float(x) * float(x) for x in v))


def run_one(version, curve, m, k1_bound, seed, d_secret=None, algo='LLL',
            beta=20, want_stats=False):
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    if d_secret is None:
        d_secret = random.Random(seed + 7777).randint(1, n - 1)
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return None

    if version == 1:
        M, sc = build_v1(sigs, n, lam, k1_bound, k2_bound)
        pv = planted_v1(sigs, n, lam, sc)
        rec = recover_v1
    else:
        M, sc = build_v2(sigs, n, lam, k1_bound, k2_bound)
        pv = planted_v2(sigs, n, lam, sc)
        rec = recover_v2

    A = IntegerMatrix.from_matrix(M)
    if algo == 'BKZ':
        BKZ.reduction(A, BKZ.Param(beta))
    else:
        LLL.reduction(A)
    dim = sc['dim']
    rows = [[A[i][j] for j in range(dim)] for i in range(dim)]

    d_cand, r_idx = rec(rows, sigs, n, lam, sc, d_secret)
    ok = d_cand is not None

    if not want_stats:
        return ok

    nrm = [norm(r) for r in rows]
    sv_i = min(range(dim), key=lambda i: nrm[i] if nrm[i] > 0 else float('inf'))
    pvn = norm(pv)
    sv = rows[sv_i]
    stats = dict(ok=ok, sv_over_pv=nrm[sv_i] / pvn, pv_norm=pvn,
                 sv_norm=nrm[sv_i], row=r_idx, dim=dim,
                 sv_kannan=abs(sv[dim - 1]) / sc['S_KANNAN'],
                 sv_is_planted=(sv == pv or sv == [-x for x in pv]),
                 det_root=None)
    # where does the shortest vector's energy sit?
    if version == 1:
        stats['sv_dcol'] = abs(sv[sc['m']]) / max(1.0, nrm[sv_i])
    else:
        stats['sv_dcol'] = 0.0
    return stats


def success_rate(version, curve, m, k1_bound, seeds, algo='LLL', beta=20):
    w = 0
    for s in seeds:
        r = run_one(version, curve, m, k1_bound, s, algo=algo, beta=beta)
        if r:
            w += 1
    return w, len(seeds)

# ---------------------------------------------------------------------------
# Gaussian-heuristic prediction
# ---------------------------------------------------------------------------

def gh_ratio(version, n, m, K1, K2):
    """||planted|| / lambda_1^GH.  <1 predicts the planted vector is findable."""
    S_K1 = n // K1
    S_K2 = max(1, n // K2)
    if version == 1:
        dim = 2 * m + 2
        log_det = m * math.log(n * S_K1) + math.log(1) + m * math.log(S_K2) + math.log(n)
        pv2 = (m / 3.0) * (K1 * S_K1) ** 2 + (n ** 2) / 3.0 \
            + (m / 3.0) * (K2 * S_K2) ** 2 + n ** 2
    else:
        dim = 2 * m + 1
        log_det = (m - 1) * math.log(n * S_K1) + math.log(S_K1) \
            + m * math.log(S_K2) + math.log(n)
        pv2 = (m / 3.0) * (K1 * S_K1) ** 2 + (m / 3.0) * (K2 * S_K2) ** 2 + n ** 2
    lam1 = math.sqrt(dim / (2 * math.pi * math.e)) * math.exp(log_det / dim)
    return math.sqrt(pv2) / lam1


# ===========================================================================
# MAIN
# ===========================================================================

print("=" * 78)
print("Thread 23 — V2 lattice (d eliminated) vs V1 (Kannan on unbounded d)")
print("=" * 78)

SEEDS = [42, 1234, 9999, 555, 31337]

HIST = [
    # label,             p,    b, n,    lam,  K1, m
    ("8-bit/199",        211,  2, 199,  106,  2,  6),
    ("12-bit/2557",      2557, 2, 2659, 1755, 8,  8),
    ("12-bit/2677 FAIL", 2677, 2, 2647, 185,  8,  10),
]
hist = []
for label, p, b, n, lam, k1, m in HIST:
    G = find_generator(p, b, n)
    assert G is not None, label
    assert (lam * lam + lam + 1) % n == 0, label
    hist.append((label, (p, b, n, lam, G), k1, m))

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("EXP A: is the planted vector lambda_1?   (sv = shortest reduced row)")
print("-" * 78)
print("sv/pv < 1 means some non-planted vector is shorter than the planted one,")
print("i.e. the Kannan embedding is NOT solving SVP.  T5 (2026-07-29) measured")
print("sv/pv in [0.34, 0.61] for V1.  V2 should push it above 1.\n")

print(f"{'curve':<18} {'K1':>3} {'m':>3} | {'V1 sv/pv':>9} {'V1 d-col':>8} "
      f"{'V1 ok':>6} | {'V2 sv/pv':>9} {'V2 row':>7} {'V2 ok':>6}")
for label, curve, k1, m in hist:
    s1 = run_one(1, curve, m, k1, 42, want_stats=True)
    s2 = run_one(2, curve, m, k1, 42, want_stats=True)
    print(f"{label:<18} {k1:>3} {m:>3} | {s1['sv_over_pv']:>9.3f} "
          f"{s1['sv_dcol']:>8.3f} {str(s1['ok']):>6} | "
          f"{s2['sv_over_pv']:>9.3f} {str(s2['row']):>7} {str(s2['ok']):>6}")

print("\nV1 'd-col' = fraction of the shortest vector's energy in the d column")
print("(1.000 == the trivial vector n*S_D*e_m).  V2 has no d column.")

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("EXP B: does the planted vector become the SHORTEST row in V2?")
print("-" * 78)
print(f"{'curve':<18} {'ver':>4} {'dim':>4} {'|pv|/n':>8} {'|sv|/n':>8} "
      f"{'sv/pv':>7} {'sv==pv':>7} {'sv Kan':>7}")
for label, curve, k1, m in hist:
    n = curve[2]
    for ver in (1, 2):
        s = run_one(ver, curve, m, k1, 42, want_stats=True)
        print(f"{label:<18} {('V'+str(ver)):>4} {s['dim']:>4} "
              f"{s['pv_norm']/n:>8.2f} {s['sv_norm']/n:>8.2f} "
              f"{s['sv_over_pv']:>7.3f} {str(s['sv_is_planted']):>7} "
              f"{s['sv_kannan']:>7.2f}")

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("EXP C: T4 replication — K1 wall on the two 12-bit curves, V1 vs V2")
print("-" * 78)
print("2026-07-29 T4 measured V1.  If V2 moves the wall outward the")
print("reformulation is a real improvement; if not, the wall is")
print("information-theoretic.\n")

K1_GRID = [2, 3, 4, 6, 8, 12, 16, 24]
for label, curve, _k1, m in hist[1:]:
    n = curve[2]; lam = curve[3]
    K2 = math.isqrt(n) + 1
    print(f"{label}   lam*={lam_star(lam,n):.3f}  n={n}  m={m}  K2={K2}")
    hdr = f"  {'K1':>4} " + " ".join(f"{k:>7}" for k in K1_GRID)
    print(hdr)
    for ver in (1, 2):
        cells = []
        for k1 in K1_GRID:
            w, t = success_rate(ver, curve, m, k1, SEEDS)
            cells.append(f"{w}/{t}")
        print(f"  {('V'+str(ver)):>4} " + " ".join(f"{c:>7}" for c in cells))
    ghs = []
    for k1 in K1_GRID:
        ghs.append(f"{gh_ratio(2, n, m, k1, K2):.2f}")
    print(f"  {'GH2':>4} " + " ".join(f"{g:>7}" for g in ghs))
    print()

# ---------------------------------------------------------------------------
print("-" * 78)
print("EXP D: eff = K1*K2/n sweep on fresh 17-bit curves, V1 vs V2")
print("-" * 78)


def find_17bit_curves(count, lo=2**16, hi=2**17):
    out = []
    p = int(sympy.nextprime(lo))
    while p < hi and len(out) < count:
        if p % 3 == 1:
            eis = eisenstein_decompose(p)
            if eis is not None:
                for t in j0_traces(*eis):
                    n_cand = p + 1 - t
                    if n_cand < 2 or not sympy.isprime(n_cand) or n_cand % 3 != 1:
                        continue
                    lam = glv_eigenvalue(n_cand)
                    if lam is None:
                        continue
                    cur = build_curve(p, n_cand)
                    if cur is None:
                        continue
                    out.append(cur)
                    break
        p = int(sympy.nextprime(p))
    return out


curves17 = find_17bit_curves(12)
print(f"found {len(curves17)} 17-bit j=0 GLV curves "
      f"(n bits {min(c[2].bit_length() for c in curves17)}"
      f"-{max(c[2].bit_length() for c in curves17)})\n")

M17 = 12
print(f"{'eff':>6} | {'V1 curves 5/5':>14} {'V1 trials':>10} | "
      f"{'V2 curves 5/5':>14} {'V2 trials':>10} | {'GH1':>5} {'GH2':>5}")
for eff in (0.05, 0.10, 0.15, 0.25):
    v1c = v2c = 0
    v1t = v2t = 0
    tot = 0
    gh1s = []; gh2s = []
    for cur in curves17:
        n = cur[2]
        K2 = math.isqrt(n) + 1
        K1 = max(2, int(eff * n / K2))
        gh1s.append(gh_ratio(1, n, M17, K1, K2))
        gh2s.append(gh_ratio(2, n, M17, K1, K2))
        w1, t = success_rate(1, cur, M17, K1, SEEDS)
        w2, _ = success_rate(2, cur, M17, K1, SEEDS)
        v1t += w1; v2t += w2; tot += t
        v1c += (w1 == t); v2c += (w2 == t)
    print(f"{eff:>6.2f} | {str(v1c)+'/'+str(len(curves17)):>14} "
          f"{str(v1t)+'/'+str(tot):>10} | "
          f"{str(v2c)+'/'+str(len(curves17)):>14} "
          f"{str(v2t)+'/'+str(tot):>10} | "
          f"{sum(gh1s)/len(gh1s):>5.2f} {sum(gh2s)/len(gh2s):>5.2f}")

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("EXP E: does BKZ change the V2 verdict at the wall?")
print("-" * 78)
label, curve, _k1, m = hist[2]     # 12-bit/2677, lam*=0.07
n = curve[2]
print(f"{label}  m={m}")
print(f"{'K1':>4} {'V2 LLL':>8} {'V2 BKZ20':>9} {'V2 BKZ40':>9} {'V1 LLL':>8}")
for k1 in [4, 6, 8, 12]:
    a, t = success_rate(2, curve, m, k1, SEEDS)
    b20, _ = success_rate(2, curve, m, k1, SEEDS, algo='BKZ', beta=20)
    b40, _ = success_rate(2, curve, m, k1, SEEDS, algo='BKZ', beta=40)
    c, _ = success_rate(1, curve, m, k1, SEEDS)
    print(f"{k1:>4} {str(a)+'/'+str(t):>8} {str(b20)+'/'+str(t):>9} "
          f"{str(b40)+'/'+str(t):>9} {str(c)+'/'+str(t):>8}")

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("EXP F: per-instance agreement V1 vs V2 (not just matching totals)")
print("-" * 78)
print("If V2 is the quotient of V1 by the trivial direction, the SAME instances")
print("should succeed.  Contingency over all (curve, seed, eff) triples.\n")

both = v1only = v2only = neither = 0
for eff in (0.05, 0.10, 0.15, 0.25):
    for cur in curves17:
        n = cur[2]
        K2 = math.isqrt(n) + 1
        K1 = max(2, int(eff * n / K2))
        for s in SEEDS:
            a = bool(run_one(1, cur, M17, K1, s))
            b = bool(run_one(2, cur, M17, K1, s))
            if a and b: both += 1
            elif a: v1only += 1
            elif b: v2only += 1
            else: neither += 1
tot = both + v1only + v2only + neither
print(f"  both succeed : {both:>4}")
print(f"  V1 only      : {v1only:>4}")
print(f"  V2 only      : {v2only:>4}")
print(f"  neither      : {neither:>4}")
print(f"  total        : {tot:>4}   agreement = "
      f"{(both+neither)/tot*100:.1f}%")

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("EXP G: what IS lambda_1 of V2?  (test the lambda-block hypothesis)")
print("-" * 78)
print("L2 = <(n*S_K1, 0), (-lam*S_K1, S_K2)> is a 2-D sublattice sitting inside")
print("BOTH V1 and V2 (modulus row i paired with the k2_i row).  Its shortest")
print("vector is the nu_hat quantity from commit e845207.  If |sv(V2)| matches")
print("lambda_1(L2), then V2's ceiling is the lambda-block, not the d column.\n")


def lagrange_gauss(u, v):
    """Exact 2-D shortest vector by Gauss reduction (integer input)."""
    def dot(a, b): return a[0] * b[0] + a[1] * b[1]
    while True:
        if dot(v, v) < dot(u, u):
            u, v = v, u
        q = round(dot(u, v) / dot(u, u))
        if q == 0:
            return u, v
        v = [v[0] - q * u[0], v[1] - q * u[1]]


print(f"{'curve':<18} {'K1':>3} | {'|sv(V2)|/n':>11} {'l1(L2)/n':>9} "
      f"{'ratio':>7} {'nu_hat':>7} {'|pv|/n':>7}")
for label, curve, k1, m in hist:
    p, b, n, lam, G = curve
    K2 = math.isqrt(n) + 1
    S_K1 = n // k1
    S_K2 = max(1, n // K2)
    u, v = lagrange_gauss([n * S_K1, 0], [-lam * S_K1, S_K2])
    l1 = math.sqrt(min(u[0] ** 2 + u[1] ** 2, v[0] ** 2 + v[1] ** 2))
    det = float(n) * S_K1 * S_K2
    s2 = run_one(2, curve, m, k1, 42, want_stats=True)
    print(f"{label:<18} {k1:>3} | {s2['sv_norm']/n:>11.3f} {l1/n:>9.3f} "
          f"{s2['sv_norm']/l1:>7.3f} {l1/math.sqrt(det):>7.3f} "
          f"{s2['pv_norm']/n:>7.3f}")

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("EXP H: is r = ||planted|| / lambda_1(L2) the recovery criterion?")
print("-" * 78)
print("L2 is the ONLY competing short family left in V2, and it encodes a real")
print("ambiguity: (k1,k2) -> (k1 - lam*t mod n, k2 + t) leaves k_full fixed and")
print("is excluded only by the bounds.  So r < 1 should mean 'planted is")
print("lambda_1' and r > 1 should mean 'a bounded-looking decoy exists'.\n")


def l1_L2(n, lam, S_K1, S_K2):
    u, v = lagrange_gauss([n * S_K1, 0], [-lam * S_K1, S_K2])
    return math.sqrt(min(u[0] ** 2 + u[1] ** 2, v[0] ** 2 + v[1] ** 2))


rows_h = []
for eff in (0.05, 0.10, 0.15, 0.25):
    for cur in curves17:
        p, b, n, lam, G = cur
        K2 = math.isqrt(n) + 1
        K1 = max(2, int(eff * n / K2))
        S_K1 = n // K1
        S_K2 = max(1, n // K2)
        l1 = l1_L2(n, lam, S_K1, S_K2)
        for s in SEEDS:
            st = run_one(2, cur, M17, K1, s, want_stats=True)
            if st is None:
                continue
            nu = l1 / math.sqrt(float(n) * S_K1 * S_K2)
            rows_h.append((st['pv_norm'] / l1, bool(st['ok']), eff, nu))


def auc_of(pairs):
    pos = [r for r, o in pairs if o]
    neg = [r for r, o in pairs if not o]
    if not pos or not neg:
        return float('nan'), len(pos), len(neg)
    wins = ties = 0
    for a in pos:
        for c in neg:
            if a < c: wins += 1
            elif a == c: ties += 1
    return (wins + 0.5 * ties) / (len(pos) * len(neg)), len(pos), len(neg)


pos = [r for r, o, _, _ in rows_h if o]
neg = [r for r, o, _, _ in rows_h if not o]
auc, _, _ = auc_of([(r, o) for r, o, _, _ in rows_h])
print(f"  n instances = {len(rows_h)}   success = {len(pos)}   fail = {len(neg)}")
print(f"  r | success : min {min(pos):.3f}  med {sorted(pos)[len(pos)//2]:.3f}"
      f"  max {max(pos):.3f}")
print(f"  r | failure : min {min(neg):.3f}  med {sorted(neg)[len(neg)//2]:.3f}"
      f"  max {max(neg):.3f}")
print(f"  AUC (lower r = success) = {auc:.3f}")

best_t, best_acc = None, 0.0
cands = sorted(set(round(r, 3) for r, _, _, _ in rows_h))
for t in cands:
    acc = (sum(1 for r in pos if r <= t) + sum(1 for r in neg if r > t)) / len(rows_h)
    if acc > best_acc:
        best_acc, best_t = acc, t
base = max(len(pos), len(neg)) / len(rows_h)
print(f"  best threshold r <= {best_t:.3f}: accuracy {best_acc*100:.1f}% "
      f"(majority baseline {base*100:.1f}%)")

print("\n  Stratified by eff (eff is the known dominant variable — an honest")
print("  test of r must hold eff fixed):")
print(f"    {'eff':>5} {'succ':>5} {'fail':>5} {'AUC(r)':>7} {'AUC(1/r)':>9} "
      f"{'AUC(nu)':>8}")
for eff in (0.05, 0.10, 0.15, 0.25):
    sub = [(r, o) for r, o, e, _ in rows_h if e == eff]
    subn = [(nu, o) for _, o, e, nu in rows_h if e == eff]
    a, np_, nn_ = auc_of(sub)
    an, _, _ = auc_of(subn)
    inv = float('nan') if a != a else 1.0 - a
    print(f"    {eff:>5.2f} {np_:>5} {nn_:>5} {a:>7.3f} {inv:>9.3f} "
          f"{an:>8.3f}")
print("\n  AUC(r) < 0.5 means the naive reading is BACKWARDS: within fixed eff,")
print("  a LARGER ||planted||/lambda_1(L2) predicts SUCCESS.  Since ||planted||")
print("  is ~constant at fixed eff, r varies through 1/lambda_1(L2), so this is")
print("  the same signal as nu_hat = lambda_1(L2)/sqrt(det L2) with the sign")
print("  reported in commit e845207 (low nu_hat -> high success).  AUC(nu) is")
print("  printed for the direct comparison.")

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("EXP I: V3 — V2 with the m-1 lambda-blocks pre-reduced (Lagrange-Gauss)")
print("-" * 78)
print("Rows (modulus j, k2_j) span L2 in columns (j-1, m+j).  Replacing that")
print("pair by its Gauss-reduced basis is unimodular, so the LATTICE is")
print("identical and any exact-SVP oracle is unaffected; only LLL's starting")
print("basis changes.  If the K1 wall moves, the wall was an LLL-budget")
print("artefact.  If not, it is information-theoretic.\n")


def build_v3(sigs, n, lam, k1_bound, k2_bound):
    M, sc = build_v2_plain(sigs, n, lam, k1_bound, k2_bound)
    m = sc['m']
    for j in range(1, m):
        r_mod, r_k2 = j - 1, m + j
        c1, c2 = j - 1, m + j
        u = [M[r_mod][c1], M[r_mod][c2]]
        v = [M[r_k2][c1], M[r_k2][c2]]
        a, b2 = lagrange_gauss(u, v)
        M[r_mod][c1], M[r_mod][c2] = a
        M[r_k2][c1], M[r_k2][c2] = b2
    return M, sc


def run_v3(curve, m, k1_bound, seed, algo='LLL', beta=20):
    global build_v2
    build_v2 = build_v3          # run_one() resolves build_v2 at call time
    try:
        return run_one(2, curve, m, k1_bound, seed, algo=algo, beta=beta)
    finally:
        build_v2 = build_v2_plain


for label, curve, _k1, m in hist[1:]:
    n = curve[2]; lam = curve[3]
    print(f"{label}   lam*={lam_star(lam,n):.3f}  m={m}")
    print(f"  {'K1':>4} " + " ".join(f"{k:>7}" for k in K1_GRID))
    for tag, fn in (("V2", lambda c, mm, k, s: run_one(2, c, mm, k, s)),
                    ("V3", lambda c, mm, k, s: run_v3(c, mm, k, s))):
        cells = []
        for k1 in K1_GRID:
            w = sum(1 for s in SEEDS if fn(curve, m, k1, s))
            cells.append(f"{w}/{len(SEEDS)}")
        print(f"  {tag:>4} " + " ".join(f"{c:>7}" for c in cells))
    print()

print("Done.")
