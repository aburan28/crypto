"""
GLV-HNP Phase 2: lam/n isolation experiment (2026-06-18).

Hypothesis from 2026-06-17: the sharp K1 phase transition at 20-bit (eff ceiling
~0.05, vs ~0.15 for 8/12-bit) might be due to lam/n≈0.33 (unbalanced GLV),
not due to the larger lattice dimension per se.

Evidence against pure-dimension explanation:
 - 8-bit:  lam/n≈0.53, eff ceiling ~0.15, dim=12 at m=5
 - 12-bit: lam/n≈0.66, eff ceiling ~0.15, dim=14 at m=6
 - 20-bit: lam/n≈0.33, eff ceiling ~0.05, fails at eff=0.076 (K1=55)

Experiment design:
 A) lam/n≈0.33 (known bad): yesterday's 20-bit curve p=524347, n=523969, lam=177902
 B) lam/n≈0.50 (balanced):  new 20-bit curve with lam/n ∈ [0.45, 0.55]
 C) lam/n≈0.33 but DIFFERENT 20-bit curve: rule out curve-specific artefact

For each: run K1 scan (K1 = 36, 55, 72, 90, 109) at m=5..18, 3 seeds.
Report the eff ceiling (highest K1 achieving 3/3).

Run: python3 glv_hnp_lamn_isolation.py
"""

import math
import random
import sympy
from fpylll import IntegerMatrix, LLL

# ---------------------------------------------------------------------------
# EC arithmetic
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
    mv, c, t, r = s, pow(z, q, p), pow(n, q, p), pow(n, (q + 1) // 2, p)
    while True:
        if t == 0: return 0
        if t == 1: return r
        i, tmp = 0, t
        while tmp != 1: tmp = tmp * tmp % p; i += 1
        b = pow(c, 1 << (mv - i - 1), p)
        mv, c, t, r = i, b * b % p, t * b * b % p, r * b % p

def find_generator(p, b, n, seed=12345):
    rng = random.Random(seed)
    for _ in range(200000):
        x = rng.randint(0, p - 1)
        rhs = (pow(x, 3, p) + b) % p
        y = tonelli_shanks(rhs, p)
        if y is not None and y != 0:
            pt = (x, y)
            if ec_mul(pt, n, p) is None:
                return pt
    return None

# ---------------------------------------------------------------------------
# CM theory for j=0 curves (p≡1 mod 3)
# ---------------------------------------------------------------------------

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
                b_val = num // 2
                if b_val >= 0 and a * a - a * b_val + b_val * b_val == p:
                    return (a, b_val)
    return None

def j0_traces(a, b):
    return [2*a - b, -2*a + b, -(a + b), a + b, 2*b - a, a - 2*b]

def glv_eigenvalue(n):
    neg3 = (n - 3) % n
    sq = tonelli_shanks(neg3, n)
    if sq is None:
        return None, None
    inv2 = modinv(2, n)
    r1 = (n - 1 + sq) * inv2 % n
    r2 = (n - 1 + (n - sq)) * inv2 % n
    if (r1 * r1 + r1 + 1) % n != 0:
        r1, r2 = r2, r1
    assert (r1 * r1 + r1 + 1) % n == 0
    lam = min(r1, r2)
    return lam, n - 1 - lam

def find_20bit_curve(lam_ratio_lo, lam_ratio_hi, start_offset=0, label=""):
    """Find a 20-bit j=0 prime-order curve with lam/n in [lo, hi]."""
    p_start = 2**19 + start_offset
    p = sympy.nextprime(p_start - 1)
    p_max = 2**21
    found = []
    while p < p_max and len(found) < 5:
        if p % 3 == 1:
            eis = eisenstein_decompose(p)
            if eis is not None:
                a_e, b_e = eis
                for t in j0_traces(a_e, b_e):
                    n_cand = p + 1 - t
                    if n_cand < 2:
                        continue
                    if sympy.isprime(n_cand) and n_cand % 3 == 1 and n_cand.bit_length() <= 21:
                        lam_val, _ = glv_eigenvalue(n_cand)
                        if lam_val is None:
                            continue
                        ratio = lam_val / n_cand
                        if lam_ratio_lo <= ratio <= lam_ratio_hi:
                            for b_try in range(1, 100):
                                G = find_generator(p, b_try, n_cand)
                                if G is not None:
                                    found.append((p, b_try, n_cand, lam_val, G))
                                    print(f"  {label} candidate: p={p}, b={b_try}, n={n_cand} "
                                          f"({n_cand.bit_length()}b), lam/n={ratio:.4f}")
                                    break
        p = sympy.nextprime(p)
    return found[0] if found else None

# ---------------------------------------------------------------------------
# GLV-HNP experiment infrastructure (same as k1_diagnostic.py)
# ---------------------------------------------------------------------------

def gen_sigs(G, d_secret, m, n, lam, p, b, k1_bound, k2_bound, seed):
    rng = random.Random(seed)
    sigs = []
    attempts = 0
    while len(sigs) < m and attempts < 500000:
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
        sigs.append({
            'A': h * modinv(s, n) % n,
            'B': r * modinv(s, n) % n,
            'k1': k1, 'k2': k2,
        })
    return sigs

def build_lattice(sigs, n, lam, k1_bound, k2_bound):
    m = len(sigs)
    dim = 2 * m + 2
    S_K1 = n // k1_bound
    S_K2 = max(1, n // k2_bound)
    S_KANNAN = n
    M = [[0] * dim for _ in range(dim)]
    for i in range(m):
        M[i][i] = n * S_K1
    for i in range(m):
        M[m][i] = sigs[i]['B'] * S_K1
    M[m][m] = 1
    for i in range(m):
        M[m + 1 + i][i] = -lam * S_K1
        M[m + 1 + i][m + 1 + i] = S_K2
    for i in range(m):
        M[2*m+1][i] = sigs[i]['A'] * S_K1
    M[2*m+1][dim-1] = S_KANNAN
    return M, S_KANNAN

def try_recover(M_red, m, n, S_KANNAN, d_secret):
    dim = 2 * m + 2
    for row in M_red:
        if abs(row[dim-1]) != S_KANNAN: continue
        sign = 1 if row[dim-1] > 0 else -1
        d_cand = (sign * row[m]) % n
        if d_cand == d_secret:
            return True
    return False

def experiment(curve, m, d_secret, k1_bound, seed):
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    sigs = gen_sigs(G, d_secret, m, n, lam, p, b, k1_bound, k2_bound, seed)
    if len(sigs) < m: return False
    M, S_KANNAN = build_lattice(sigs, n, lam, k1_bound, k2_bound)
    dim = 2 * m + 2
    A = IntegerMatrix.from_matrix(M)
    LLL.reduction(A)
    reduced = [[A[i][j] for j in range(dim)] for i in range(dim)]
    return try_recover(reduced, m, n, S_KANNAN, d_secret)

def k1_sweep(label, curve, k1_values, m_range, seeds):
    """For each K1 in k1_values, sweep m and find first 3/3."""
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    print(f"\n{'='*60}")
    print(f"Curve: {label}")
    print(f"  n={n} ({n.bit_length()}b), lam/n={lam/n:.4f}, K2={k2_bound}")
    print(f"{'='*60}")

    results = {}
    for k1 in k1_values:
        eff = k1 * k2_bound / n
        m_thresh = math.ceil(math.log(n) / (-math.log(eff)))
        print(f"\n  K1={k1:3d}  eff={eff:.4f}  m_thresh={m_thresh}")
        first_full = None
        for m in m_range:
            wins = sum(
                experiment(curve, m, random.Random(s + 9001).randint(1, n - 1), k1, s)
                for s in seeds
            )
            marker = " ← m_thresh" if m == m_thresh else ""
            print(f"    m={m:2d}: {wins}/{len(seeds)}{marker}")
            if wins == len(seeds) and first_full is None:
                first_full = m
        ratio = f"{first_full/m_thresh:.2f}" if first_full else "N/A"
        results[k1] = (first_full, m_thresh, eff)
        print(f"    → first {len(seeds)}/{len(seeds)}: m={first_full}, ratio={ratio}")
    return results

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

print("=" * 60)
print("GLV-HNP lam/n isolation: fix 20-bit, vary lam/n (2026-06-18)")
print("=" * 60)

SEEDS = [42, 1234, 9999]
K1_VALUES = [36, 55, 72, 90, 109]
M_RANGE = range(5, 20)

# ---- Curve A: lam/n≈0.33 (known bad from 2026-06-17) ----------------------
print("\n[A] Using known 20-bit curve with lam/n≈0.33 (p=524347)...")
# Already verified: p=524347, b=2, n=523969, lam=177902
P_A, B_A, N_A, LAM_A = 524347, 2, 523969, 177902
G_A = find_generator(P_A, B_A, N_A)
assert G_A is not None and ec_mul(G_A, N_A, P_A) is None
CURVE_A = (P_A, B_A, N_A, LAM_A, G_A)
print(f"  OK: lam/n = {LAM_A/N_A:.4f}")

# ---- Curve B: lam/n≈0.50 (balanced, need to find) --------------------------
print("\n[B] Searching for 20-bit curve with lam/n ∈ [0.45, 0.55]...")
CURVE_B = find_20bit_curve(0.45, 0.55, label="B")
if CURVE_B is None:
    print("  FALLBACK: searching lam/n ∈ [0.40, 0.60]...")
    CURVE_B = find_20bit_curve(0.40, 0.60, label="B")

# ---- Curve C: different 20-bit j=0 curve with lam/n≈0.33 ------------------
print("\n[C] Searching for second 20-bit curve with lam/n ∈ [0.28, 0.39] "
      "(to rule out curve-specific artefact)...")
CURVE_C = find_20bit_curve(0.28, 0.39, start_offset=100000, label="C")
if CURVE_C is None:
    print("  FALLBACK: lam/n ∈ [0.25, 0.40]...")
    CURVE_C = find_20bit_curve(0.25, 0.40, start_offset=50000, label="C")

# ---- Run K1 sweeps --------------------------------------------------------
results_A = k1_sweep(f"A: lam/n={LAM_A/N_A:.4f} p={P_A}", CURVE_A, K1_VALUES, M_RANGE, SEEDS)

if CURVE_B:
    pB, bB, nB, lamB, _ = CURVE_B
    results_B = k1_sweep(f"B: lam/n={lamB/nB:.4f} p={pB}", CURVE_B, K1_VALUES, M_RANGE, SEEDS)
else:
    print("\n[B] FAILED to find balanced-lam/n curve. Cannot run Curve B sweep.")
    results_B = None

if CURVE_C:
    pC, bC, nC, lamC, _ = CURVE_C
    results_C = k1_sweep(f"C: lam/n={lamC/nC:.4f} p={pC}", CURVE_C, K1_VALUES, M_RANGE, SEEDS)
else:
    print("\n[C] FAILED to find second lam/n≈0.33 curve.")
    results_C = None

# ---- Summary ---------------------------------------------------------------
print("\n" + "=" * 60)
print("SUMMARY TABLE: eff ceiling by lam/n")
print("=" * 60)
print("'eff ceiling' = highest eff achieving 3/3 success in m≤19")
print()
print(f"{'Curve':<8} {'lam/n':>6} {'K1=36(e=0.05)':>14} {'K1=55(e=0.076)':>15} "
      f"{'K1=72(e=0.10)':>14} {'K1=90(e=0.124)':>15} {'K1=109(e=0.15)':>15}")
print("-" * 90)

for tag, res, curve in [("A", results_A, CURVE_A),
                         ("B", results_B, CURVE_B),
                         ("C", results_C, CURVE_C)]:
    if res is None or curve is None:
        print(f"  {tag:<6}  (not available)")
        continue
    _, _, n, lam, _ = curve
    ratio = lam / n
    row = f"  {tag:<6}  {ratio:.4f}  "
    for k1 in K1_VALUES:
        first_full, m_thresh, eff = res[k1]
        cell = str(first_full) if first_full else "FAIL"
        row += f"  {cell:>12}"
    print(row)

print()
print("Interpretation:")
print("  If Curve B (lam/n≈0.50) succeeds at K1=55 but Curve A (lam/n≈0.33) fails,")
print("  this confirms lam/n is the determining factor for the eff ceiling.")
print("  If both fail at K1=55, the failure is due to lattice dimension (bit size), not lam/n.")
print()
print("Done.")
