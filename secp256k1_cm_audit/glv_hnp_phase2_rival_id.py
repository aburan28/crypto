"""
GLV-HNP Phase 2 / Thread 23b -- identify the rival shortest vector of the
PROJECTED lattice.

P5 of `glv_hnp_phase2_projected.py` found that after the d-column is quotiented
out, LLL's b_1 has (a) Kannan coordinate exactly 0 and (b) a k1/k2 block-energy
split that is IDENTICAL across signature seeds for a given curve.  A
signature-independent rival can only come from the part of the basis that does
not involve the signatures, i.e. the m copies of the 2D lambda-block

    L2 = < (n*S_K1, 0), (-lam*S_K1, S_K2) >      (columns i and m+i)

which is exactly the sublattice whose normalised shortest vector
nu_hat = lambda_1(L2)/sqrt(det L2) was found on 2026-07-29 to separate the
C1/C2 success classes at AUC 0.935 (commit e845207).

This script tests the identification directly: is LLL's b_1 equal to the
Lagrange-reduced shortest vector of L2, placed in some coordinate pair
(i, m+i)?  If yes then lambda_1(projected lattice) = lambda_1(L2) exactly,
nu_hat is not a heuristic proxy but the exact first minimum, and the Thread-23
reformulation cannot help: it removes one trivial vector only to expose another.

Self-contained (repo idiom: each glv_hnp_*.py carries its own EC arithmetic) so
that importing it never re-runs another script's experiments.

Run: python3 glv_hnp_phase2_rival_id.py
Requires: fpylll   (pip install fpylll cysignals)
"""

import math
import random
from fpylll import IntegerMatrix, LLL

# ---------------------------------------------------------------------------
# Minimal EC arithmetic (short Weierstrass, a=0)
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
        if y is not None and y != 0 and ec_mul((x, y), n, p) is None:
            return (x, y)
    return None

def lam_star(lam, n):
    return min(lam % n, n - (lam % n)) / n

def norm(v):
    return math.sqrt(sum(x * x for x in v))

# ---------------------------------------------------------------------------
# Phase-2 scaling, signatures, projected lattice (same as Thread 23)
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

def build_lattice_proj(sigs, n, lam, k1_bound, k2_bound):
    m = len(sigs)
    dim = 2 * m + 1
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    rows = []
    for i in range(m):
        r = [0] * dim; r[i] = n * S_K1; rows.append(r)
    r = [0] * dim
    for i in range(m):
        r[i] = sigs[i]['B'] * S_K1
    rows.append(r)
    for i in range(m):
        r = [0] * dim
        r[i] = -lam * S_K1
        r[m + i] = S_K2
        rows.append(r)
    r = [0] * dim
    for i in range(m):
        r[i] = sigs[i]['A'] * S_K1
    r[dim - 1] = S_KAN
    rows.append(r)
    return rows

def planted_proj(sigs, d_secret, n, k1_bound, k2_bound):
    m = len(sigs)
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    v = [0] * (2 * m + 1)
    for i in range(m):
        v[i] = sigs[i]['k1'] * S_K1
        v[m + i] = sigs[i]['k2'] * S_K2
    v[2 * m] = S_KAN
    return v

def recover_proj(rows, sigs, n, lam, k1_bound, k2_bound, d_secret):
    m = len(sigs)
    dim = 2 * m + 1
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    for row in rows:
        if abs(row[dim - 1]) != S_KAN: continue
        sign = 1 if row[dim - 1] > 0 else -1
        w = [sign * x for x in row]
        if any(w[i] % S_K1 for i in range(m)): continue
        if any(w[m + i] % S_K2 for i in range(m)): continue
        for i in range(m):
            B = sigs[i]['B'] % n
            if math.gcd(B, n) != 1: continue
            d_cand = (w[i] // S_K1 - sigs[i]['A'] + lam * (w[m + i] // S_K2)) \
                * modinv(B, n) % n
            if d_cand and d_cand == d_secret:
                return True
    return False

def lll_rows(M):
    A = IntegerMatrix.from_matrix(M)
    LLL.reduction(A)
    return [[A[i][j] for j in range(A.ncols)] for i in range(A.nrows)]

def shortest_nonzero(rows):
    best, bv = None, None
    for r in rows:
        nr = norm(r)
        if nr > 0 and (best is None or nr < best):
            best, bv = nr, r
    return best, bv

# ---------------------------------------------------------------------------
# Exact 2D shortest vector (Lagrange/Gauss), verbatim from
# glv_hnp_phase2_lambda_threshold.py:188
# ---------------------------------------------------------------------------

def gauss_reduce_2d(u, v):
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
        q = (2 * num + den) // (2 * den) if num >= 0 else -((-2 * num + den) // (2 * den))
        v = (v[0] - q * u[0], v[1] - q * u[1])
        if nrm2(v) >= nrm2(u):
            break
        u, v = v, u
    return u

def lambda_block(n, lam, S_K1, S_K2):
    """(shortest vector w of L2, ||w||, sqrt(det L2))."""
    w = gauss_reduce_2d((n * S_K1, 0), (-(lam % n) * S_K1, S_K2))
    return w, math.sqrt(w[0] ** 2 + w[1] ** 2), math.sqrt(n * S_K1 * S_K2)

# ===========================================================================

SEEDS = [42, 1234, 9999, 555, 31337]
HIST = [
    ("8-bit/199",        211,  2, 199,  106,  2,  6),
    ("12-bit/2557",      2557, 2, 2659, 1755, 8,  8),
    ("12-bit/2677 FAIL", 2677, 2, 2647, 185,  8,  10),
]

print("=" * 78)
print("Thread 23b - is LLL b_1 of the projected lattice the lambda-block vector?")
print("=" * 78)

hist = []
for label, p, b, n, lam, k1, m in HIST:
    G = find_generator(p, b, n)
    assert G is not None, label
    assert (lam * lam + lam + 1) % n == 0, label
    hist.append((label, (p, b, n, lam, G), k1, m))

print("\nR1: exact identification of b_1.  'match' = b_1 equals +-w placed in")
print("    coordinate pair (i, m+i);  w = Lagrange-shortest vector of L2.")
print("-" * 78)
print(f"{'curve':<18} {'K1':>3} {'seed':>6} | {'||b1||':>11} {'||w||':>11} "
      f"{'match':>6} {'pair':>5} | {'nu_hat':>7} {'sv/pv':>7}")
all_match = True
for label, curve, k1, m in hist:
    p, b, n, lam, G = curve
    k2b = math.isqrt(n) + 1
    S_K1, S_D, S_K2, S_KAN = scales(n, k1, k2b)
    w, wn, sdet = lambda_block(n, lam, S_K1, S_K2)
    nu = wn / sdet
    for seed in SEEDS:
        d = random.Random(seed ^ 0xD).randrange(1, n)
        sigs = gen_signatures(G, d, m, n, lam, p, k1, k2b, seed)
        if len(sigs) < m:
            continue
        red = lll_rows(build_lattice_proj(sigs, n, lam, k1, k2b))
        pv = planted_proj(sigs, d, n, k1, k2b)
        sn, sv = shortest_nonzero(red)
        pair, match = None, False
        for i in range(m):
            cand = [0] * (2 * m + 1)
            cand[i], cand[m + i] = w[0], w[1]
            if sv == cand or sv == [-x for x in cand]:
                pair, match = i, True
                break
        all_match &= match
        print(f"{label:<18} {k1:>3} {seed:>6} | {sn:>11.4f} {wn:>11.4f} "
              f"{str(match):>6} {str(pair):>5} | {nu:>7.4f} {sn/norm(pv):>7.3f}")
print(f"\nR1 verdict: {'ALL MATCH' if all_match else 'not all match'}")

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("R2: gap ratio  g = ||v_planted|| / ||w||  vs recovery, across the K1 grid.")
print("    g < 1 would mean the planted vector is the first minimum.")
print("-" * 78)
K1_GRID = [2, 3, 4, 6, 8, 12, 16, 24]
for label, curve, _k1, m in hist[1:]:
    p, b, n, lam, G = curve
    k2b = math.isqrt(n) + 1
    print(f"\n{label}  lam*={lam_star(lam,n):.3f}  m={m}")
    print(f"  {'K1':<9}" + "  ".join(f"{k:>7}" for k in K1_GRID))
    grow, nrow, srow = [], [], []
    for k1 in K1_GRID:
        S_K1, S_D, S_K2, S_KAN = scales(n, k1, k2b)
        w, wn, sdet = lambda_block(n, lam, S_K1, S_K2)
        wins = tot = 0
        pvs = []
        for seed in SEEDS:
            d = random.Random(seed ^ 0xD).randrange(1, n)
            sigs = gen_signatures(G, d, m, n, lam, p, k1, k2b, seed)
            if len(sigs) < m:
                continue
            red = lll_rows(build_lattice_proj(sigs, n, lam, k1, k2b))
            tot += 1
            wins += bool(recover_proj(red, sigs, n, lam, k1, k2b, d))
            pvs.append(norm(planted_proj(sigs, d, n, k1, k2b)))
        mean_pv = sum(pvs) / len(pvs) if pvs else float('nan')
        grow.append(f"{mean_pv/wn:.3f}")
        nrow.append(f"{wn/sdet:.3f}")
        srow.append(f"{wins}/{tot}")
    print(f"  {'g':<9}" + "  ".join(f"{x:>7}" for x in grow))
    print(f"  {'nu_hat':<9}" + "  ".join(f"{x:>7}" for x in nrow))
    print(f"  {'recover':<9}" + "  ".join(f"{x:>7}" for x in srow))

# ---------------------------------------------------------------------------
# R3 -- the sign test.
#
# Two readings of R1/R2 make OPPOSITE predictions:
#   (i)  SVP-gap reading: recovery happens when the planted vector beats the
#        rival w, i.e. when g = ||v_planted||/||w|| is SMALL.  Since
#        ||v_planted|| ~ sqrt(m/3)*n and ||w|| = nu_hat*sqrt(n*S_K1*S_K2),
#        g ~ sqrt(m*eff/3)/nu_hat, so this predicts HIGH nu_hat = easier.
#   (ii) the 2026-07-29 nu_hat result (commit e845207, AUC 0.935): LOW nu_hat
#        = easier.
# These cannot both hold.  Sweep fresh curves at fixed eff and m and see which
# sign the data takes.
#
# R3 also tests the m-scaling corollary of (i): with eff fixed, g ~ sqrt(m),
# so (i) predicts recovery must DEGRADE as m grows -- the 2026-07-29 T4b
# observation that more signatures never rescue the K1 wall.
# ---------------------------------------------------------------------------
import sympy

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

print("\n" + "-" * 78)
print("R3a: m-scaling at fixed eff.  12-bit/2677, K1=8.  Reading (i) predicts")
print("     g ~ sqrt(m) and therefore no rescue from more signatures.")
print("-" * 78)
_, curve2677, _, _ = hist[2]
p, b, n, lam, G = curve2677
k2b = math.isqrt(n) + 1
S_K1, S_D, S_K2, S_KAN = scales(n, 8, k2b)
w, wn, sdet = lambda_block(n, lam, S_K1, S_K2)
print(f"  nu_hat = {wn/sdet:.4f} (independent of m);  eff = {8*k2b/n:.4f}")
print(f"  {'m':>4} {'g':>8} {'g/sqrt(m)':>10} {'recover':>9}")
for mm in (8, 12, 16, 24, 32):
    wins = tot = 0
    pvs = []
    for seed in SEEDS:
        d = random.Random(seed ^ 0xD).randrange(1, n)
        sigs = gen_signatures(G, d, mm, n, lam, p, 8, k2b, seed)
        if len(sigs) < mm:
            continue
        red = lll_rows(build_lattice_proj(sigs, n, lam, 8, k2b))
        tot += 1
        wins += bool(recover_proj(red, sigs, n, lam, 8, k2b, d))
        pvs.append(norm(planted_proj(sigs, d, n, 8, k2b)))
    g = (sum(pvs) / len(pvs)) / wn
    print(f"  {mm:>4} {g:>8.3f} {g/math.sqrt(mm):>10.4f} {wins:>5}/{tot:<3}")

print("\n" + "-" * 78)
print("R3b: cross-curve sign test.  Fresh 17-bit j=0 GLV curves, m=12,")
print("     eff = K1*K2/n held at ~0.15.  Which of nu_hat / g tracks recovery?")
print("-" * 78)

def search_curves(lo, hi, want=14):
    out = []
    pp = int(sympy.nextprime(lo))
    while pp < hi and len(out) < want:
        if pp % 3 == 1:
            eis = eisenstein_decompose(pp)
            if eis is not None:
                a_e, b_e = eis
                for t in j0_traces(a_e, b_e):
                    nc = pp + 1 - t
                    if nc < 2 or nc % 3 != 1 or not sympy.isprime(nc):
                        continue
                    roots = glv_roots(nc)
                    if roots is None:
                        continue
                    rng = random.Random(nc)
                    bb = None
                    for _ in range(400):
                        cand = rng.randint(1, pp - 1)
                        x = rng.randint(0, pp - 1)
                        rhs = (pow(x, 3, pp) + cand) % pp
                        y = tonelli_shanks(rhs, pp)
                        if y is None or y == 0:
                            continue
                        if ec_mul((x, y), nc, pp) is None:
                            bb = cand
                            break
                    if bb is None:
                        continue
                    Gg = find_generator(pp, bb, nc)
                    if Gg is None:
                        continue
                    out.append((pp, bb, nc, roots[0], Gg))
                    break
        pp = int(sympy.nextprime(pp))
    return out

rows17 = []
print(f"  {'p':>7} {'n':>7} {'lam*':>6} {'K1':>4} {'eff':>6} {'nu_hat':>7} "
      f"{'g':>7} {'recover':>9}")
for cur in search_curves(2 ** 16, 2 ** 17, want=14):
    pp, bb, nn, ll, Gg = cur
    kb = math.isqrt(nn) + 1
    k1 = max(2, round(0.15 * nn / kb))
    sk1, _sd, sk2, _sk = scales(nn, k1, kb)
    ww, wwn, wsdet = lambda_block(nn, ll, sk1, sk2)
    nu = wwn / wsdet
    wins = tot = 0
    pvs = []
    for seed in SEEDS:
        d = random.Random(seed ^ 0xD).randrange(1, nn)
        sigs = gen_signatures(Gg, d, 12, nn, ll, pp, k1, kb, seed)
        if len(sigs) < 12:
            continue
        red = lll_rows(build_lattice_proj(sigs, nn, ll, k1, kb))
        tot += 1
        wins += bool(recover_proj(red, sigs, nn, ll, k1, kb, d))
        pvs.append(norm(planted_proj(sigs, d, nn, k1, kb)))
    if not tot:
        continue
    g = (sum(pvs) / len(pvs)) / wwn
    rows17.append((nu, g, wins / tot))
    print(f"  {pp:>7} {nn:>7} {lam_star(ll,nn):>6.3f} {k1:>4} {k1*kb/nn:>6.3f} "
          f"{nu:>7.4f} {g:>7.3f} {wins:>5}/{tot:<3}")

def spearman(xs, ys):
    def rk(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
                j += 1
            avg = (i + j) / 2.0 + 1
            for k in range(i, j + 1):
                r[order[k]] = avg
            i = j + 1
        return r
    a, bq = rk(xs), rk(ys)
    na = len(xs)
    ma, mb = sum(a) / na, sum(bq) / na
    num = sum((a[i] - ma) * (bq[i] - mb) for i in range(na))
    da = math.sqrt(sum((a[i] - ma) ** 2 for i in range(na)))
    db = math.sqrt(sum((bq[i] - mb) ** 2 for i in range(na)))
    return num / (da * db) if da and db else float('nan')

if rows17:
    nus = [r[0] for r in rows17]
    gs = [r[1] for r in rows17]
    ps = [r[2] for r in rows17]
    print(f"\n  n_curves = {len(rows17)}")
    print(f"  spearman(nu_hat, success) = {spearman(nus, ps):+.3f}   "
          f"(2026-07-29 sign: NEGATIVE, low nu_hat = easier)")
    print(f"  spearman(g,      success) = {spearman(gs, ps):+.3f}   "
          f"(SVP-gap reading (i) predicts NEGATIVE, small g = easier)")
    print("  note g = sqrt(m*eff/3)/nu_hat up to O(1), so the two are")
    print("  near-deterministically anti-related at fixed (m, eff);")
    print("  whichever sign appears, only one reading survives.")

print("\n" + "=" * 78)
print("done")
print("=" * 78)
