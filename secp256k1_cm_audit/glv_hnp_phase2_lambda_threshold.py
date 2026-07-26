"""
GLV-HNP Phase 2 — λ/n threshold bisection.

Finds j=0 GLV curves with lam/n in 7 bins spanning [0.07, 0.50],
runs LLL attack on each, identifies the success/failure boundary.

Prior data (2026-07-26 autolab run):
  lam/n = 0.07 (n=2647):  FAIL at m=5..12 (BKZ-40 also fails)
  lam/n = 0.34 (n=523969): OK at m=9
  lam/n = 0.53 (n=199):   OK at m=4
  lam/n = 0.66 (n=2659):  OK at m=7

Run: python3 secp256k1_cm_audit/glv_hnp_phase2_lambda_threshold.py
"""

import math, random, sympy
from fpylll import IntegerMatrix, LLL

# ─── EC arithmetic ────────────────────────────────────────────────────────────

def modinv(a, m): return pow(a, -1, m)

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
        Q = ec_add(Q, Q, p); k >>= 1
    return R

def tonelli(n, p):
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

# ─── CM theory ────────────────────────────────────────────────────────────────

def eisenstein_decompose(p):
    for a in range(1, 2 * math.isqrt(p // 3) + 3):
        disc = 4 * p - 3 * a * a
        if disc < 0: break
        s = math.isqrt(disc)
        if s * s != disc: continue
        for num in [a + s, a - s]:
            if num % 2 == 0:
                b = num // 2
                if b >= 0 and a * a - a * b + b * b == p: return (a, b)
    return None

def j0_traces(a, b):
    return [2*a - b, -2*a + b, -(a + b), a + b, 2*b - a, a - 2*b]

def glv_eigenvalue(n):
    sq = tonelli((n - 3) % n, n)
    if sq is None: return None
    inv2 = modinv(2, n)
    r1 = (n - 1 + sq) * inv2 % n
    r2 = (n - 1 + (n - sq)) * inv2 % n
    if (r1 * r1 + r1 + 1) % n != 0: r1, r2 = r2, r1
    assert (r1 * r1 + r1 + 1) % n == 0, f"Bad GLV root for n={n}"
    return min(r1, r2)

def quick_twist_check(p, b, n, seed=777):
    """Fast check: does y²=x³+b over F_p have order divisible by n?
    Uses one non-zero point: n*P=O iff yes (n prime → order exactly n)."""
    rng = random.Random(seed)
    for _ in range(80):
        x = rng.randint(0, p - 1)
        rhs = (pow(x, 3, p) + b) % p
        y = tonelli(rhs, p)
        if y is None or y == 0: continue
        return ec_mul((x, y), n, p) is None
    return False

def find_generator(p, b, n, seed=12345):
    rng = random.Random(seed)
    for _ in range(30000):
        x = rng.randint(0, p - 1)
        rhs = (pow(x, 3, p) + b) % p
        y = tonelli(rhs, p)
        if y is not None and y != 0:
            P = (x, y)
            if ec_mul(P, n, p) is None: return P
    return None

# ─── Curve search ─────────────────────────────────────────────────────────────

def find_curve_for_ratio(lo, hi, min_n=128, max_n=2**15):
    """Return (p,b,n,lam,G) with lam/n in [lo,hi), n prime, n∈[min_n,max_n]."""
    p = sympy.nextprime(2**11)
    while p < 2**17:
        if p % 3 == 1:
            eis = eisenstein_decompose(p)
            if eis is not None:
                a_e, b_e = eis
                for t in j0_traces(a_e, b_e):
                    n = p + 1 - t
                    if n < min_n or n > max_n: continue
                    if not sympy.isprime(n): continue
                    if n % 3 != 1: continue
                    lam = glv_eigenvalue(n)
                    if lam is None: continue
                    if not (lo <= lam / n < hi): continue
                    # Found a candidate; find the right twist b
                    for b_try in range(1, 400):
                        if not quick_twist_check(p, b_try, n): continue
                        G = find_generator(p, b_try, n)
                        if G is not None:
                            return (p, b_try, n, lam, G)
        p = sympy.nextprime(p)
    return None

# ─── Attack ───────────────────────────────────────────────────────────────────

def gen_sigs(G, d, m, n, lam, p, b, k1b, k2b, seed):
    rng = random.Random(seed); sigs = []; att = 0
    while len(sigs) < m and att < 500000:
        att += 1
        k1 = rng.randint(0, k1b - 1); k2 = rng.randint(0, k2b - 1)
        kf = (k1 + lam * k2) % n
        if kf == 0: continue
        R = ec_mul(G, kf, p)
        if R is None: continue
        r = R[0] % n
        if r == 0: continue
        h = rng.randint(0, n - 1)
        s = modinv(kf, n) * (h + d * r) % n
        if s == 0: continue
        sinv = modinv(s, n)
        A = h * sinv % n; B = r * sinv % n
        assert (A + B * d) % n == kf
        sigs.append({'A': A, 'B': B})
    return sigs

def run_lll_once(G, d, m, n, lam, p, b, k1b, k2b, seed):
    sigs = gen_sigs(G, d, m, n, lam, p, b, k1b, k2b, seed)
    if len(sigs) < m: return False
    dim = 2 * m + 2
    S1 = max(1, n // k1b); S2 = max(1, n // k2b); SK = n
    M = [[0] * dim for _ in range(dim)]
    for i in range(m): M[i][i] = n * S1
    for i in range(m): M[m][i] = sigs[i]['B'] * S1
    M[m][m] = 1
    for i in range(m):
        M[m + 1 + i][i] = -lam * S1
        M[m + 1 + i][m + 1 + i] = S2
    for i in range(m): M[2 * m + 1][i] = sigs[i]['A'] * S1
    M[2 * m + 1][dim - 1] = SK
    A = IntegerMatrix.from_matrix(M); LLL.reduction(A)
    for i in range(dim):
        last = A[i][dim - 1]
        if abs(last) != SK: continue
        sign = 1 if last > 0 else -1
        dc = (sign * A[i][m]) % n
        if dc == d: return True
    return False

def test_curve(curve, k1b, m_max, seeds):
    p, b, n, lam, G = curve; k2b = math.isqrt(n) + 1
    eff = k1b * k2b / n
    m_thresh = math.ceil(math.log(n) / math.log(1 / eff)) if eff < 1 else 1
    for m in range(3, m_max + 1):
        wins = sum(
            run_lll_once(G, random.Random(s + 999).randint(1, n - 1),
                         m, n, lam, p, b, k1b, k2b, s)
            for s in seeds
        )
        if wins == len(seeds):
            return m, eff, m_thresh
    return None, eff, m_thresh

# ─── Main ─────────────────────────────────────────────────────────────────────

K1B   = 8
SEEDS = [42, 1234, 9999]
M_MAX = 16

BINS = [
    (0.07, 0.10, "0.07-0.10"),
    (0.10, 0.14, "0.10-0.14"),
    (0.14, 0.19, "0.14-0.19"),
    (0.19, 0.25, "0.19-0.25"),
    (0.25, 0.31, "0.25-0.31"),
    (0.31, 0.38, "0.31-0.38"),
    (0.38, 0.50, "0.38-0.50"),
]

print("=" * 65)
print("GLV-HNP Phase 2 — λ/n threshold bisection")
print(f"K1_BOUND={K1B} (fixed), K2_BOUND=√n+1, m=3..{M_MAX}, seeds={SEEDS}")
print("=" * 65)
print()

rows = []
for lo, hi, label in BINS:
    print(f"[{label}] Searching for j=0 curve with lam/n in [{lo:.2f},{hi:.2f})...",
          flush=True)
    curve = find_curve_for_ratio(lo, hi)
    if curve is None:
        print(f"  NO CURVE FOUND in search range")
        rows.append((label, None, None, None, None))
        continue
    p, b, n, lam, G = curve
    ratio = lam / n
    k2b = math.isqrt(n) + 1
    print(f"  Found: p={p}, n={n} ({n.bit_length()}b), lam={lam}, lam/n={ratio:.4f}")
    print(f"  K1={K1B}, K2={k2b}  → running LLL at m=3..{M_MAX}...", flush=True)
    m_ok, eff, m_thresh = test_curve(curve, K1B, M_MAX, SEEDS)
    status = f"OK at m={m_ok}" if m_ok else f"FAIL (m>{M_MAX})"
    print(f"  Result: {status}  (eff={eff:.4f}, m_thresh≈{m_thresh})")
    rows.append((label, ratio, m_ok, eff, m_thresh))

print()
print("=" * 65)
print("SUMMARY")
print("=" * 65)
print(f"{'bin':<13} {'lam/n':>8} {'m_ok':>8} {'eff':>8} {'m_thresh':>10}")
print("-" * 52)
for label, ratio, m_ok, eff, mt in rows:
    rs  = f"{ratio:.4f}" if ratio is not None else "N/A"
    ms  = str(m_ok) if m_ok is not None else "FAIL"
    es  = f"{eff:.4f}"  if eff   is not None else "N/A"
    mts = str(mt)       if mt    is not None else "N/A"
    print(f"{label:<13} {rs:>8} {ms:>8} {es:>8} {mts:>10}")

print()
print("Reference (prior runs):")
print("  lam/n=0.07 (n=2647):   FAIL at m=5..12 (BKZ-40 also fails)")
print("  lam/n=0.34 (n=523969): OK at m=9")
print("  lam/n=0.53 (n=199):    OK at m=4")
print("  lam/n=0.66 (n=2659):   OK at m=7")
