"""
GLV-HNP Phase 2 — λ/n threshold bisection study.

From 2026-07-26 log:
  λ/n=0.07 (n=2647): LLL FAILS even at m=12, BKZ(40) FAILS
  λ/n=0.34 (n=523969): LLL succeeds at m=9
  λ/n=0.53 (n=199):   LLL succeeds at m=4
  λ/n=0.66 (n=2659):  LLL succeeds at m=7

Goal: bisect the threshold between 0.07 and 0.34, holding eff≈0.05 fixed.

Key simplification: do NOT need an actual elliptic curve. The lattice attack
only requires (n, lam, A_i, B_i) where A_i + B_i*d ≡ k1_i + lam*k2_i (mod n).
We can simulate this with random (h_i, r_i) pairs:
  k_full = k1 + lam*k2 mod n;  s = modinv(k_full,n)*(h+d*r) mod n;
  A = h/s mod n;  B = r/s mod n.
The lattice structure is identical — no actual EC arithmetic needed.

Method: Eisenstein CM search finds 12-13 bit prime-order j=0 curves
(n prime, n≡1 mod 3) covering λ/n ∈ [0.01, 0.50] in steps of ~0.01.
Attack is run for each curve at m = 2..20, 3 seeds, eff≈0.05.

Run: python3 glv_hnp_phase2_lambda_threshold.py
"""

import math
import random
from fpylll import IntegerMatrix, LLL, BKZ

# ---------------------------------------------------------------------------
# Primality (Miller-Rabin, deterministic for n < 3.3e24)
# ---------------------------------------------------------------------------
def isprime(n):
    if n < 2: return False
    if n < 4: return True
    if n % 2 == 0: return False
    r, d = 0, n - 1
    while d % 2 == 0: r += 1; d //= 2
    for a in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]:
        if a >= n: continue
        x = pow(a, d, n)
        if x == 1 or x == n - 1: continue
        for _ in range(r - 1):
            x = x * x % n
            if x == n - 1: break
        else: return False
    return True

# ---------------------------------------------------------------------------
# Tonelli-Shanks sqrt mod p
# ---------------------------------------------------------------------------
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
        if t == 1: return r
        i, tmp = 0, t
        while tmp != 1: tmp = tmp * tmp % p; i += 1
        b = pow(c, 1 << (m - i - 1), p)
        m, c, t, r = i, b*b%p, t*b*b%p, r*b%p

# ---------------------------------------------------------------------------
# GLV eigenvalue: smallest root of x²+x+1≡0 mod n
# ---------------------------------------------------------------------------
def glv_lambda(n):
    if n < 7: return None
    sq = tonelli((n - 3) % n, n)
    if sq is None: return None
    inv2 = pow(2, -1, n)
    r1 = (n - 1 + sq) * inv2 % n
    r2 = (n - 1 + n - sq) * inv2 % n
    if (r1 * r1 + r1 + 1) % n != 0: r1, r2 = r2, r1
    if (r1 * r1 + r1 + 1) % n != 0: return None
    return min(r1, r2)

# ---------------------------------------------------------------------------
# Eisenstein CM: find all prime-order j=0 curves in a prime range
# ---------------------------------------------------------------------------
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

def find_cm_curves(p_lo, p_hi):
    """Return list of (ratio, p, n, lam) for 12-13 bit prime-order j=0 curves."""
    curves = []
    for p in range(p_lo | 1, p_hi, 2):
        if not isprime(p) or p % 3 != 1: continue
        eis = eisenstein_decompose(p)
        if eis is None: continue
        a_e, b_e = eis
        for t in j0_traces(a_e, b_e):
            n = p + 1 - t
            if n < 500: continue
            if not isprime(n) or n % 3 != 1: continue
            lam = glv_lambda(n)
            if lam is None: continue
            ratio = lam / n
            curves.append((ratio, p, n, lam))
    curves.sort()
    return curves

# ---------------------------------------------------------------------------
# Simulated HNP signatures — no actual EC needed.
# Gen: k1 ∈ [0,K1), k2 ∈ [1,K2); k_full = k1+lam*k2 mod n.
# s = modinv(k_full,n)*(h+d*r) mod n; A=h/s mod n; B=r/s mod n.
# Then A+B*d ≡ k_full (mod n). ✓
# ---------------------------------------------------------------------------
def gen_sigs(d, m, n, lam, K1, K2, seed=42):
    rng = random.Random(seed)
    sigs = []
    for _ in range(500000):
        if len(sigs) == m: break
        k1 = rng.randint(0, K1 - 1)
        k2 = rng.randint(1, K2 - 1)
        k_full = (k1 + lam * k2) % n
        if k_full == 0: continue
        h = rng.randint(1, n - 1)
        r = rng.randint(1, n - 1)
        s = pow(k_full, -1, n) * (h + d * r) % n
        if s == 0: continue
        si = pow(s, -1, n)
        A, B = h * si % n, r * si % n
        assert (A + B * d) % n == k_full
        sigs.append({'A': A, 'B': B, 'k1': k1, 'k2': k2})
    return sigs

# ---------------------------------------------------------------------------
# Build (2m+2)×(2m+2) lattice with column-diagonal scaling
# ---------------------------------------------------------------------------
def build_lattice(sigs, n, lam, K1, K2):
    m = len(sigs)
    dim = 2 * m + 2
    S1 = n // K1
    S2 = max(1, n // K2)
    Sk = n
    M = [[0] * dim for _ in range(dim)]
    for i in range(m): M[i][i] = n * S1
    for i in range(m): M[m][i] = sigs[i]['B'] * S1
    M[m][m] = 1
    for i in range(m):
        M[m + 1 + i][i] = -lam * S1
        M[m + 1 + i][m + 1 + i] = S2
    for i in range(m): M[2*m+1][i] = sigs[i]['A'] * S1
    M[2*m+1][dim - 1] = Sk
    return M, Sk

# ---------------------------------------------------------------------------
# Try to recover d from LLL/BKZ reduced basis
# ---------------------------------------------------------------------------
def try_recover(M_list, d, n, Sk, use_bkz=False, beta=20):
    A = IntegerMatrix.from_matrix(M_list)
    if use_bkz:
        BKZ.reduction(A, BKZ.Param(beta))
    else:
        LLL.reduction(A)
    m = (A.nrows - 2) // 2
    for i in range(A.nrows):
        if abs(A[i][A.ncols - 1]) == Sk:
            sign = 1 if A[i][A.ncols - 1] > 0 else -1
            cand = sign * A[i][m] % n
            if cand == d:
                return True
    return False

# ---------------------------------------------------------------------------
# Attack one curve across m values
# ---------------------------------------------------------------------------
def attack_curve(n, lam, seeds, m_max=18, use_bkz=False, beta=20):
    K1 = max(2, int(0.05 * math.isqrt(n)))
    K2 = math.isqrt(n) + 1
    eff = K1 * K2 / n
    d_secret = random.Random(0xDEAD).randint(2, n - 2)
    for m in range(2, m_max + 1):
        wins = 0
        for seed in seeds:
            sigs = gen_sigs(d_secret, m, n, lam, K1, K2, seed=seed)
            if len(sigs) < m: continue
            M, Sk = build_lattice(sigs, n, lam, K1, K2)
            if try_recover(M, d_secret, n, Sk, use_bkz=use_bkz, beta=beta):
                wins += 1
        if wins == len(seeds):
            return m, eff
    return None, eff

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    print("=" * 70)
    print("GLV-HNP Phase 2 — λ/n threshold bisection (no EC needed)")
    print("=" * 70)
    print()
    print("Finding 12-13 bit j=0 prime-order curves via Eisenstein CM ...")
    all_curves = find_cm_curves(2003, 8500)
    print(f"Found {len(all_curves)} curves with n prime, n≡1 mod 3.")
    print()

    # Bin by λ/n in steps of 0.01, pick smallest n in each bin
    bins = {}
    for ratio, p, n, lam in all_curves:
        bk = round(ratio * 100) / 100
        if bk not in bins or n < bins[bk][2]:
            bins[bk] = (ratio, p, n, lam)

    # Collect representatives with λ/n in [0.05, 0.50], step 0.02
    targets = [k for k in sorted(bins.keys()) if 0.03 <= k <= 0.50]
    reps = [(bins[k][0], bins[k][1], bins[k][2], bins[k][3]) for k in targets]

    SEEDS = (42, 137, 999)

    print(f"{'λ/n':>7} {'p':>7} {'n':>8} {'λ':>8} {'K1':>5} {'K2':>6} {'eff':>7} {'m*':>6} status")
    print("-" * 75)

    results = []
    for ratio, p, n, lam in reps:
        K1 = max(2, int(0.05 * math.isqrt(n)))
        K2 = math.isqrt(n) + 1
        m_star, eff = attack_curve(n, lam, SEEDS, m_max=18)
        status = f"m*={m_star}" if m_star else "FAIL"
        results.append((ratio, n, lam, m_star, eff))
        print(f"{ratio:>7.3f} {p:>7} {n:>8} {lam:>8} {K1:>5} {K2:>6} {eff:>7.4f} {str(m_star) if m_star else '-':>6} {status}")

    # Summary
    print()
    print("=" * 70)
    print("THRESHOLD ANALYSIS")
    print("=" * 70)
    print()
    print("λ/n   n        m*    eff")
    for ratio, n, lam, m_star, eff in results:
        m_str = str(m_star) if m_star else "FAIL"
        print(f"{ratio:.3f}  {n:<8} {m_str:<6} {eff:.4f}")

    # Find threshold bracket
    fails = [(r, n, eff) for r, n, lam, m, eff in results if m is None]
    passes = [(r, n, eff) for r, n, lam, m, eff in results if m is not None]

    print()
    print("FAIL region:")
    for r, n, eff in fails:
        print(f"  λ/n={r:.3f}, n={n}, eff={eff:.4f}")
    print("PASS region (sample):")
    for r, n, eff in passes[:8]:
        print(f"  λ/n={r:.3f}, n={n}, eff={eff:.4f}")

    # Threshold
    fail_ratios = sorted([r for r, n, eff in fails])
    pass_ratios = sorted([r for r, n, eff in passes])
    if fail_ratios and pass_ratios:
        hi_fail = max(fail_ratios)
        lo_pass = min(pass_ratios)
        if hi_fail > lo_pass:
            print(f"\nWARNING: non-monotone results (hi_fail={hi_fail:.3f} > lo_pass={lo_pass:.3f})")
            print("  This suggests threshold depends on more than just λ/n (n-size effects).")
        else:
            print(f"\nThreshold bracket: ({hi_fail:.3f}, {lo_pass:.3f})")
            print(f"Estimated threshold: λ/n ≈ {(hi_fail + lo_pass)/2:.3f}")

    print()
    print("Prior data points (2026-07-26 log, 12-bit curves):")
    print("  λ/n=0.07 (n=2647, eff≈0.039): FAIL (LLL and BKZ(40))")
    print("  λ/n=0.66 (n=2659, eff≈0.050): PASS at m=7")

if __name__ == "__main__":
    main()
