"""
GLV-HNP Phase 2: λ/n threshold sweep.

Find j=0 toy curves (y^2 = x^3 + b) with prime order n ≡ 1 mod 3 and
test the GLV-aware HNP attack across a range of λ/n ratios.

Questions:
  Q1. Below what λ/n does LLL fail?
  Q2. Does switching to the complementary eigenvalue λ₂ = n-1-λ₁ fix failures?

Approach:
  1. Search primes p ≡ 1 mod 3 and b values for prime-order j=0 curves.
  2. Find GLV eigenvalues λ (roots of x²+x+1 ≡ 0 mod n).
  3. Run column-scaled GLV-HNP lattice attack with BOTH eigenvalues.
  4. Report success rate vs λ/n for representative curves at each ratio bucket.

Run: python3 glv_hnp_phase2_lambda_threshold.py
"""

import math
import random
from fpylll import IntegerMatrix, LLL

# ---------------------------------------------------------------------------
# Primality (trial division — fine for n < 10000)
# ---------------------------------------------------------------------------
def isprime(n):
    if n < 2: return False
    if n == 2: return True
    if n % 2 == 0: return False
    for d in range(3, int(n**0.5) + 1, 2):
        if n % d == 0: return False
    return True

# ---------------------------------------------------------------------------
# Elliptic curve over F_p  (y^2 = x^3 + b, a=0)
# ---------------------------------------------------------------------------
def count_points(p, b):
    """Brute-force #E(F_p) for y^2 = x^3 + b."""
    count = 1  # point at infinity
    for x in range(p):
        rhs = (x**3 + b) % p
        if rhs == 0:
            count += 1
        elif pow(rhs, (p - 1) // 2, p) == 1:
            count += 2
    return count

def ec_add(P, Q, p):
    if P is None: return Q
    if Q is None: return P
    x1, y1 = P; x2, y2 = Q
    if x1 == x2:
        if (y1 + y2) % p == 0: return None
        s = 3 * x1 * x1 * pow(2 * y1, -1, p) % p
    else:
        s = (y2 - y1) * pow(x2 - x1, -1, p) % p
    x3 = (s * s - x1 - x2) % p
    y3 = (s * (x1 - x3) - y1) % p
    return (x3, y3)

def ec_mul(P, k, p):
    R, Q = None, P
    while k > 0:
        if k & 1: R = ec_add(R, Q, p)
        Q = ec_add(Q, Q, p)
        k >>= 1
    return R

def find_generator(p, n, b):
    """Find an affine point of order n on y^2 = x^3 + b over F_p."""
    for x in range(p):
        rhs = (x**3 + b) % p
        if rhs == 0 or pow(rhs, (p - 1) // 2, p) != 1:
            continue
        for y in range(p):
            if y * y % p == rhs:
                P = (x, y)
                if ec_mul(P, n, p) is None:
                    return P
    return None

# ---------------------------------------------------------------------------
# GLV eigenvalues: roots of x^2 + x + 1 ≡ 0 mod n
# ---------------------------------------------------------------------------
def glv_roots(n):
    """Return sorted list of roots of x^2+x+1 mod n. Only works if n ≡ 1 mod 3."""
    roots = []
    for x in range(1, n):
        if (x * x + x + 1) % n == 0:
            roots.append(x)
            if len(roots) == 2:
                break
    return sorted(roots)

# ---------------------------------------------------------------------------
# Curve database: search for j=0 prime-order curves with GLV eigenvalues
# ---------------------------------------------------------------------------
def build_curve_db(p_max=800):
    """
    Search primes p ≡ 1 mod 3 and b ∈ [1,p) for j=0 curves
    with prime group order n ≡ 1 mod 3.
    Return list of (p, b, n, lam1, lam2, ratio1, ratio2).
    """
    db = []
    seen_n = set()  # avoid duplicates with same (n, lam) but different p

    for p in range(7, p_max + 1):
        if not isprime(p): continue
        if p % 3 != 1: continue
        for b in range(1, p):
            n = count_points(p, b)
            if n in seen_n: continue
            if not isprime(n): continue
            if n % 3 != 1: continue
            if n < 20: continue  # too tiny for meaningful lattice
            roots = glv_roots(n)
            if len(roots) < 2: continue
            lam1, lam2 = roots[0], roots[1]
            r1 = lam1 / n
            r2 = lam2 / n
            db.append((p, b, n, lam1, lam2, r1, r2))
            seen_n.add(n)
            break  # one b per prime p

    return db

# ---------------------------------------------------------------------------
# GLV-HNP attack
# ---------------------------------------------------------------------------
def run_attack(p, n, lam, b, n_trials=5, k1_bound=2):
    """
    Run GLV-HNP lattice attack.
    Returns (success_rate, m_sigs) or (None, None) on error.
    """
    k2_bound = int(math.isqrt(n)) + 1
    G = find_generator(p, n, b)
    if G is None:
        return None, None
    if (lam * lam + lam + 1) % n != 0:
        return None, None

    eff = k1_bound * k2_bound / n
    if eff >= 1: return None, None
    m_thresh = math.ceil(math.log(n) / math.log(1.0 / eff))
    m_sigs = max(4, m_thresh + 2)

    s_k1 = n // k1_bound
    s_d = 1
    s_k2 = n // k2_bound
    s_kannan = n

    wins = 0
    for trial in range(n_trials):
        rng = random.Random(trial * 7919 + int(lam) + int(p) * 17)
        d_secret = rng.randint(1, n - 1)

        sigs = []
        attempts = 0
        while len(sigs) < m_sigs and attempts < 200000:
            attempts += 1
            k1 = rng.randint(0, k1_bound - 1)
            k2 = rng.randint(0, k2_bound - 1)
            k_full = (k1 + lam * k2) % n
            if k_full == 0: continue
            R = ec_mul(G, k_full, p)
            if R is None: continue
            r_val = R[0] % n
            if r_val == 0: continue
            h = rng.randint(0, n - 1)
            s = pow(k_full, -1, n) * (h + d_secret * r_val) % n
            if s == 0: continue
            s_inv = pow(s, -1, n)
            A = h * s_inv % n
            B = r_val * s_inv % n
            sigs.append({'A': A, 'B': B, 'k1': k1, 'k2': k2})

        if len(sigs) < m_sigs:
            return None, None

        # Build (2m+2)×(2m+2) lattice with column-diagonal scaling
        dim = 2 * m_sigs + 2
        M = [[0] * dim for _ in range(dim)]

        for i in range(m_sigs):
            M[i][i] = n * s_k1
        M[m_sigs][m_sigs] = s_d
        for i in range(m_sigs):
            M[m_sigs][i] = sigs[i]['B'] * s_k1
        for i in range(m_sigs):
            M[m_sigs + 1 + i][i] = -lam * s_k1
            M[m_sigs + 1 + i][m_sigs + 1 + i] = s_k2
        for i in range(m_sigs):
            M[2 * m_sigs + 1][i] = sigs[i]['A'] * s_k1
        M[2 * m_sigs + 1][dim - 1] = s_kannan

        A_mat = IntegerMatrix.from_matrix(M)
        LLL.reduction(A_mat)
        reduced = [[A_mat[i][j] for j in range(dim)] for i in range(dim)]

        d_found = False
        for row in reduced:
            last = row[dim - 1]
            if abs(last) != s_kannan: continue
            sign = 1 if last > 0 else -1
            d_cand = (sign * row[m_sigs]) % n
            if d_cand == 0: continue
            if d_cand == d_secret:
                d_found = True
                break
        if d_found:
            wins += 1

    return wins / n_trials, m_sigs

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    print("=" * 72)
    print("GLV-HNP Phase 2 — λ/n threshold sweep")
    print("=" * 72)
    print()

    print("Building curve database (p ≡ 1 mod 3, prime n ≡ 1 mod 3)...")
    db = build_curve_db(p_max=1000)
    print(f"Found {len(db)} distinct curves.")
    print()

    if not db:
        print("No curves found. Increase p_max.")
        exit(1)

    # Show full database sorted by ratio1
    db_sorted = sorted(db, key=lambda x: x[5])
    print(f"{'p':>5}  {'b':>4}  {'n':>6}  {'lam1':>6}  {'r1':>6}  {'lam2':>6}  {'r2':>6}")
    print("-" * 48)
    for (p, b, n, l1, l2, r1, r2) in db_sorted[:30]:
        print(f"{p:>5}  {b:>4}  {n:>6}  {l1:>6}  {r1:>6.3f}  {l2:>6}  {r2:>6.3f}")
    print()

    # Select representative curves: pick one per 0.05-wide bin of min(r1,r2)
    bucket_width = 0.05
    bins = {}
    for entry in db_sorted:
        p, b, n, l1, l2, r1, r2 = entry
        r_min = min(r1, r2)
        bin_key = round(r_min / bucket_width) * bucket_width
        bin_key = round(bin_key, 2)
        if bin_key not in bins:
            bins[bin_key] = []
        bins[bin_key].append(entry)

    # For each bin, pick the curve with n closest to 200
    target_n = 200
    selected = {}
    for bk, entries in sorted(bins.items()):
        entries.sort(key=lambda x: abs(x[2] - target_n))
        selected[bk] = entries[0]

    print("=" * 72)
    print("ATTACK SWEEP")
    print(f"{'bin':>5}  {'p':>5}  {'n':>6}  {'lam1':>5}  {'r1':>5}  "
          f"{'r2':>5}  {'m':>3}  {'lam1 success':>13}  {'lam2 success':>13}  {'rescue?':>7}")
    print("-" * 72)

    results = []
    for bk in sorted(selected.keys()):
        p, b, n, l1, l2, r1, r2 = selected[bk]
        sr1, m1 = run_attack(p, n, l1, b, n_trials=5)
        sr2, m2 = run_attack(p, n, l2, b, n_trials=5)

        m_used = m1 or m2
        s1 = f"{sr1:.2f}" if sr1 is not None else "ERR"
        s2 = f"{sr2:.2f}" if sr2 is not None else "ERR"

        rescue = ""
        if sr1 is not None and sr2 is not None:
            if sr1 < 0.6 and sr2 >= 0.6:
                rescue = "YES"
            elif sr1 < 0.6 and sr2 < 0.6:
                rescue = "NO"
            elif sr1 >= 0.6:
                rescue = "N/A"

        print(f"{bk:>5.2f}  {p:>5}  {n:>6}  {l1:>5}  {r1:>5.3f}  "
              f"{r2:>5.3f}  {m_used!s:>3}  {s1:>13}  {s2:>13}  {rescue:>7}")

        results.append({
            'bin': bk, 'p': p, 'n': n, 'b': b,
            'lam1': l1, 'lam2': l2, 'r1': r1, 'r2': r2,
            'sr1': sr1, 'sr2': sr2, 'm': m_used, 'rescue': rescue
        })

    print()
    print("=" * 72)
    print("ANALYSIS")
    print("=" * 72)

    valid = [r for r in results if r['sr1'] is not None]
    fails1 = [r for r in valid if r['sr1'] < 0.6]
    wins1  = [r for r in valid if r['sr1'] >= 0.6]

    if fails1 and wins1:
        max_fail = max(r['r1'] for r in fails1)
        min_win  = min(r['r1'] for r in wins1)
        print(f"\nλ₁ attack threshold:")
        print(f"  Highest failing r1 = {max_fail:.3f}")
        print(f"  Lowest passing r1  = {min_win:.3f}")
        print(f"  Threshold in [{max_fail:.2f}, {min_win:.2f}]")
    elif not fails1:
        print("\nλ₁ attack: no failures in tested range (all pass).")
    elif not wins1:
        print("\nλ₁ attack: all failures in tested range.")

    rescued = [r for r in results if r.get('rescue') == 'YES']
    not_rescued = [r for r in results if r.get('rescue') == 'NO']
    print(f"\nComplementary λ₂ rescue:")
    print(f"  Rescued (λ₁ fails, λ₂ succeeds): {len(rescued)}")
    print(f"  Not rescued (both fail):           {len(not_rescued)}")
    for r in rescued:
        print(f"    r1={r['r1']:.3f} → r2={r['r2']:.3f}: λ₂ success={r['sr2']:.2f}")
    for r in not_rescued:
        print(f"    r1={r['r1']:.3f} → r2={r['r2']:.3f}: λ₂ success={r['sr2']:.2f}")

    print()
    print("=" * 72)
    print("KEY NUMERICAL RESULTS")
    print("=" * 72)
    for r in sorted(results, key=lambda x: x['r1']):
        s1 = f"{r['sr1']:.2f}" if r['sr1'] is not None else "ERR"
        s2 = f"{r['sr2']:.2f}" if r['sr2'] is not None else "ERR"
        print(f"  r1={r['r1']:.3f} r2={r['r2']:.3f}  "
              f"lam1_success={s1}  lam2_success={s2}  rescue={r.get('rescue','')}")
