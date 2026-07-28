"""
GLV-HNP Phase 2: lambda/n threshold bisection.

Goal: identify where in [0.07, 0.34] the LLL attack transitions from
always-fail to sometimes-succeed.  Scans j=0 curves with varying lam/n
ratios (using the Eisenstein CM search), bins them, and runs LLL sweeps
at m_thresh .. m_thresh+8.

Known data points (from prior runs):
  lam/n=0.07 (p=2677, n=2647):  FAIL at all m, even BKZ(40)
  lam/n=0.34 (p=524347, n=523969): SUCCEED at m=9
  lam/n=0.53 (n=199):  SUCCEED at m=4
  lam/n=0.66 (n=2659): SUCCEED at m=7

Run: python3 glv_hnp_phase2_lambda_threshold.py
"""

import math
import random
import sympy
from fpylll import IntegerMatrix, LLL

# ── Minimal EC arithmetic ────────────────────────────────────────────────────

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

def tonelli(n, p):
    n %= p
    if n == 0: return 0
    if pow(n, (p - 1) // 2, p) != 1: return None
    if p % 4 == 3: return pow(n, (p + 1) // 4, p)
    q, s = p - 1, 0
    while q % 2 == 0: q //= 2; s += 1
    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1: z += 1
    mc, c, t, r = s, pow(z, q, p), pow(n, q, p), pow(n, (q + 1) // 2, p)
    while True:
        if t == 0: return 0
        if t == 1: return r
        i, tmp = 0, t
        while tmp != 1: tmp = tmp * tmp % p; i += 1
        b = pow(c, 1 << (mc - i - 1), p)
        mc, c, t, r = i, b*b%p, t*b*b%p, r*b%p

def find_generator(p, b, n):
    rng = random.Random(12345)
    for _ in range(50000):
        x = rng.randint(0, p - 1)
        rhs = (pow(x, 3, p) + b) % p
        y = tonelli(rhs, p)
        if y is not None and y != 0:
            P = (x, y)
            if ec_mul(P, n, p) is None:
                return P
    return None

# ── CM theory ────────────────────────────────────────────────────────────────

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
    sq = tonelli(neg3, n)
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

# ── Curve catalog per lam/n bin ───────────────────────────────────────────────

BINS = [
    (0.05, 0.10),
    (0.10, 0.15),
    (0.15, 0.20),
    (0.20, 0.25),
    (0.25, 0.30),
    (0.30, 0.35),
]

def build_catalog(n_bits_min=8, n_bits_max=22, max_per_bin=2):
    """
    Scan j=0 primes and collect up to max_per_bin qualifying curves per bin.
    Returns dict: bin_label -> list of (p, b, n, lam, G).
    """
    catalog = {b: [] for b in BINS}
    p = sympy.nextprime(2**(n_bits_min - 1))
    p_max = 2**n_bits_max
    print(f"Scanning primes {2**(n_bits_min-1)}-{p_max} for j=0 curves ...")

    checked = 0
    while p < p_max:
        checked += 1
        if p % 3 != 1:
            p = sympy.nextprime(p)
            continue
        eis = eisenstein_decompose(p)
        if eis is None:
            p = sympy.nextprime(p)
            continue
        a_e, b_e = eis
        traces = j0_traces(a_e, b_e)
        for t in traces:
            n_cand = p + 1 - t
            if n_cand < 8:
                continue
            if not sympy.isprime(n_cand):
                continue
            if n_cand % 3 != 1:
                continue
            lam, lam2 = glv_eigenvalue(n_cand)
            if lam is None:
                continue
            ratio = lam / n_cand
            for bn in BINS:
                lo, hi = bn
                if lo <= ratio < hi:
                    if len(catalog[bn]) < max_per_bin:
                        # Find b for this curve
                        b_try = None
                        for bv in range(1, min(p, 200)):
                            rhs = pow(1, 3, p)  # dummy
                            # Quick filter: try bv
                            found_pt = None
                            rng_tmp = random.Random(42)
                            for _ in range(500):
                                x = rng_tmp.randint(0, p - 1)
                                rhs = (pow(x, 3, p) + bv) % p
                                y = tonelli(rhs, p)
                                if y is not None and y != 0:
                                    found_pt = (x, y)
                                    break
                            if found_pt is None:
                                continue
                            if ec_mul(found_pt, n_cand, p) is None:
                                b_try = bv
                                break
                        if b_try is None:
                            continue
                        G = find_generator(p, b_try, n_cand)
                        if G is None:
                            continue
                        catalog[bn].append((p, b_try, n_cand, lam, G))
                        print(f"  Bin [{lo:.2f},{hi:.2f}): p={p}, n={n_cand} "
                              f"({n_cand.bit_length()}b), lam={lam}, lam/n={ratio:.4f}")
        # Check if all bins full
        if all(len(catalog[bn]) >= max_per_bin for bn in BINS):
            break
        p = sympy.nextprime(p)

    print(f"  Scanned {checked} primes total.")
    return catalog

# ── GLV lattice & recovery ────────────────────────────────────────────────────

def gen_signatures(G, d_secret, m, n, lam, p, b, k1_bound, k2_bound, seed=42):
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
        sigs.append({'A': A, 'B': B})
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

def recover_d(M_red, m, n, S_KANNAN, d_secret):
    dim = 2 * m + 2
    for row in M_red:
        if abs(row[dim-1]) != S_KANNAN: continue
        sign = 1 if row[dim-1] > 0 else -1
        d_cand = (sign * row[m]) % n
        if d_cand == 0: continue
        if d_cand == d_secret:
            return True
    return False

def run_lll(curve, m, k1_bound, seed=42):
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    d_secret = random.Random(seed + 7777).randint(1, n - 1)
    sigs = gen_signatures(G, d_secret, m, n, lam, p, b, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return False
    M, S_KANNAN = build_lattice(sigs, n, lam, k1_bound, k2_bound)
    A = IntegerMatrix.from_matrix(M)
    LLL.reduction(A)
    dim = 2*m+2
    reduced = [[A[i][j] for j in range(dim)] for i in range(dim)]
    return recover_d(reduced, m, n, S_KANNAN, d_secret)

def sweep(label, curve, k1_bound, m_range, seeds):
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    eff = k1_bound * k2_bound / n
    m_thresh = math.ceil(math.log(n) / math.log(1.0 / eff)) if eff < 1 else 1
    print(f"  {label}: n={n}b, lam/n={lam/n:.4f}, K1={k1_bound}, "
          f"eff={eff:.4f}, m_thresh≈{m_thresh}")
    results = {}
    first_success = None
    for m in m_range:
        wins = sum(run_lll(curve, m, k1_bound, seed) for seed in seeds)
        results[m] = (wins, len(seeds))
        marker = "**3/3**" if wins == len(seeds) else f"{wins}/{len(seeds)}"
        thresh_tag = "†" if m == m_thresh else " "
        print(f"    m={m}{thresh_tag}: {marker}")
        if wins == len(seeds) and first_success is None:
            first_success = m
    return results, first_success, m_thresh

# ── Main ──────────────────────────────────────────────────────────────────────

print("=" * 70)
print("GLV-HNP Phase 2 — lambda/n threshold bisection")
print("=" * 70)
print()

SEEDS = [42, 1234, 9999]

# Build catalog
catalog = build_catalog(n_bits_min=8, n_bits_max=22, max_per_bin=2)

print()
print("=" * 70)
print("SWEEP RESULTS")
print("=" * 70)
print()

# Known failure point (pin it first)
print("=== Known anchor: lam/n=0.07 (FAIL) ===")
G_anchor = find_generator(2677, 2, 2647)
curve_anchor = (2677, 2, 2647, 185, G_anchor)
k1_anchor = max(2, int(0.05 * math.sqrt(2647)))
k2_anchor = math.isqrt(2647) + 1
eff_anchor = k1_anchor * k2_anchor / 2647
mt_anchor = math.ceil(math.log(2647) / math.log(1.0 / eff_anchor)) if eff_anchor < 1 else 1
res_anchor, fs_anchor, mt_anchor = sweep(
    "anchor 12-bit/2647", curve_anchor, k1_anchor,
    range(mt_anchor, mt_anchor + 9), SEEDS
)

print()
summary_rows = []

for bn in BINS:
    lo, hi = bn
    curves = catalog[bn]
    if not curves:
        print(f"=== Bin [{lo:.2f},{hi:.2f}): NO CURVE FOUND ===")
        summary_rows.append((f"[{lo:.2f},{hi:.2f})", None, None, None, None))
        continue

    for idx, curve in enumerate(curves):
        p, b, n, lam, G = curve
        k1_bound = max(2, int(0.05 * math.sqrt(n)))
        k2_bound = math.isqrt(n) + 1
        eff = k1_bound * k2_bound / n
        m_thresh_pre = math.ceil(math.log(n) / math.log(1.0 / eff)) if eff < 1 else 1

        label = f"{n.bit_length()}-bit/{n} (lam/n={lam/n:.4f})"
        print(f"=== Bin [{lo:.2f},{hi:.2f}) curve {idx+1}: {label} ===")
        res, first_s, mt = sweep(
            label, curve, k1_bound,
            range(max(2, m_thresh_pre), m_thresh_pre + 9), SEEDS
        )
        summary_rows.append((f"[{lo:.2f},{hi:.2f})", lam/n, n, first_s, mt))
        print()

print()
print("=" * 70)
print("THRESHOLD SUMMARY")
print("=" * 70)
print(f"{'Bin':<14} {'lam/n':>8} {'n':>10} {'m_thresh':>9} {'first_3/3':>10}  verdict")
print("-" * 70)
# Anchor
print(f"{'anchor':14} {0.07:>8.4f} {2647:>10} {mt_anchor:>9} {'FAIL':>10}  FAIL")
for (bname, ratio, n, fs, mt) in summary_rows:
    if ratio is None:
        print(f"{bname:<14} {'?':>8} {'N/A':>10} {'?':>9} {'N/A':>10}  NO CURVE")
        continue
    verdict = f"m={fs}" if fs is not None else "FAIL"
    overhead = f"+{fs-mt}" if fs is not None else "—"
    print(f"{bname:<14} {ratio:>8.4f} {n:>10} {mt:>9} {verdict:>10}  overhead={overhead}")

print()
print("Done.")
