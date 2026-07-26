"""
GLV-HNP Phase 2 — extended check for borderline bins.

From the threshold bisection run, three bins FAIL at m<=16:
  0.07-0.10 (lam/n=0.0941, n=2137)
  0.10-0.14 (lam/n=0.1294, n=2203)
  0.31-0.38 (lam/n=0.3145, n=2251)  ← unexpected: 0.29 and 0.40 both succeed

Questions:
  Q1. Does 0.31-0.38 succeed at m=17-25? (borderline, not structural?)
  Q2. Does the failure depend on the specific curve, or is it bin-wide?

This script:
  - Tests 0.31-0.38, 0.25-0.31 (control), and 0.38-0.50 (control) at m=3..25
    with 5 seeds, using the same curves from the threshold bisection
  - Searches for a SECOND curve in the 0.31-0.38 bin and tests it
  - Reports whether the failure is curve-specific or bin-wide

Run: python3 secp256k1_cm_audit/glv_hnp_phase2_extended_check.py
"""

import math, random, sympy
from fpylll import IntegerMatrix, LLL

# ─── Reuse core functions ─────────────────────────────────────────────────────

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
    return min(r1, r2)

def quick_twist_check(p, b, n, seed=777):
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

def find_curves_for_ratio(lo, hi, min_n=128, max_n=2**15, count=3):
    """Return up to `count` (p,b,n,lam,G) curves with lam/n in [lo,hi)."""
    curves = []
    p = sympy.nextprime(2**11)
    while p < 2**17 and len(curves) < count:
        if p % 3 == 1:
            eis = eisenstein_decompose(p)
            if eis is not None:
                a_e, b_e = eis
                for t in j0_traces(a_e, b_e):
                    if len(curves) >= count: break
                    n = p + 1 - t
                    if n < min_n or n > max_n: continue
                    if not sympy.isprime(n): continue
                    if n % 3 != 1: continue
                    lam = glv_eigenvalue(n)
                    if lam is None: continue
                    if not (lo <= lam / n < hi): continue
                    for b_try in range(1, 400):
                        if not quick_twist_check(p, b_try, n): continue
                        G = find_generator(p, b_try, n)
                        if G is not None:
                            curves.append((p, b_try, n, lam, G))
                            break
        p = sympy.nextprime(p)
    return curves

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

def sweep(curve, k1b, m_range, seeds, label=""):
    p, b, n, lam, G = curve; k2b = math.isqrt(n) + 1
    ratio = lam / n
    eff = k1b * k2b / n
    print(f"\n  [{label}] p={p}, n={n}({n.bit_length()}b), lam={lam}, lam/n={ratio:.4f}")
    print(f"  K1={k1b}, K2={k2b}, eff={eff:.4f}")
    first_ok = None
    for m in m_range:
        wins = sum(
            run_lll_once(G, random.Random(s + 999).randint(1, n - 1),
                         m, n, lam, p, b, k1b, k2b, s)
            for s in seeds
        )
        marker = " ← FIRST SUCCESS" if (wins == len(seeds) and first_ok is None) else ""
        print(f"    m={m:2d}: {wins}/{len(seeds)}{marker}")
        if wins == len(seeds) and first_ok is None:
            first_ok = m
    return first_ok

# ─── Main ─────────────────────────────────────────────────────────────────────

K1B   = 8
SEEDS = [42, 1234, 9999, 7, 31415]   # 5 seeds for better statistics
M_MAX = 25

print("=" * 65)
print("GLV-HNP Phase 2 — extended check for borderline bins")
print(f"K1_BOUND={K1B}, K2_BOUND=√n+1, m=3..{M_MAX}, {len(SEEDS)} seeds")
print("=" * 65)

# ── Q1: Does 0.31-0.38 eventually succeed at m>16? ───────────────────────────
print("\nQ1: 0.31-0.38 bin — extended to m=25 (FAIL at m<=16 before)")
print("Searching for up to 3 curves in [0.31, 0.38)...")
curves_31 = find_curves_for_ratio(0.31, 0.38, count=3)
print(f"Found {len(curves_31)} curve(s).")
results_31 = {}
for i, curve in enumerate(curves_31):
    m_ok = sweep(curve, K1B, range(3, M_MAX + 1), SEEDS,
                 label=f"0.31-0.38 curve #{i+1}")
    results_31[i] = (curve[2], curve[3] / curve[2], m_ok)  # (n, ratio, m_ok)

# ── Q2: Control — 0.25-0.31 and 0.38-0.50 at extended m ─────────────────────
print("\nQ2 (control): 0.25-0.31 bin — was OK at m=12 with 3 seeds; verify with 5 seeds")
curves_25 = find_curves_for_ratio(0.25, 0.31, count=1)
results_25 = {}
for i, curve in enumerate(curves_25):
    m_ok = sweep(curve, K1B, range(3, M_MAX + 1), SEEDS,
                 label=f"0.25-0.31 curve #{i+1}")
    results_25[i] = (curve[2], curve[3] / curve[2], m_ok)

print("\nQ2 (control): 0.38-0.50 bin — was OK at m=7; verify with 5 seeds")
curves_38 = find_curves_for_ratio(0.38, 0.50, count=1)
results_38 = {}
for i, curve in enumerate(curves_38):
    m_ok = sweep(curve, K1B, range(3, M_MAX + 1), SEEDS,
                 label=f"0.38-0.50 curve #{i+1}")
    results_38[i] = (curve[2], curve[3] / curve[2], m_ok)

# ── Summary ───────────────────────────────────────────────────────────────────
print("\n" + "=" * 65)
print("SUMMARY")
print("=" * 65)
print(f"\n0.31-0.38 bin ({len(curves_31)} curves):")
for i, (n, ratio, m_ok) in results_31.items():
    status = f"OK at m={m_ok}" if m_ok else f"FAIL (m>{M_MAX})"
    print(f"  curve {i+1}: n={n}, lam/n={ratio:.4f} → {status}")

print(f"\n0.25-0.31 bin (control, {len(curves_25)} curve):")
for i, (n, ratio, m_ok) in results_25.items():
    status = f"OK at m={m_ok}" if m_ok else f"FAIL (m>{M_MAX})"
    print(f"  curve {i+1}: n={n}, lam/n={ratio:.4f} → {status}")

print(f"\n0.38-0.50 bin (control, {len(curves_38)} curve):")
for i, (n, ratio, m_ok) in results_38.items():
    status = f"OK at m={m_ok}" if m_ok else f"FAIL (m>{M_MAX})"
    print(f"  curve {i+1}: n={n}, lam/n={ratio:.4f} → {status}")

all_31_fail = all(m_ok is None for _, _, m_ok in results_31.values())
if all_31_fail:
    print(f"\nCONCLUSION: 0.31-0.38 bin fails at ALL curves, m<={M_MAX}")
    print("  → Either structural failure or requires m >> 25")
else:
    print(f"\nCONCLUSION: 0.31-0.38 bin succeeds for some curves")
    print("  → Failure in prior run was curve-specific (high variance)")
