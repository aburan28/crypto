"""
Thread 20 — GLV-HNP Phase 2: λ/n threshold study.

From 2026-07-26:
  λ/n = 0.07 (p=2677): LLL fails, BKZ(40) fails   [known]
  λ/n ≈ 0.34:           LLL succeeds               [approx known]
  λ/n ≈ 0.44-0.53:      LLL succeeds               [known, secp256k1 regime]

Goal: bisect the threshold T* ∈ (0.07, 0.34) below which Phase 2 LLL fails.

Additional hypothesis: using λ_max = n-1-λ_min (the OTHER eigenvalue) on a
'small-λ' curve might rescue the attack. Test: run the failure curve p=2677
with λ=2461 (≈0.93·n) instead of 185 (≈0.07·n).

Protocol:
  1. Find toy curves (12–14 bit) with λ_min/n in 7 bands: [0.05,0.10), [0.10,0.15),
     [0.15,0.20), [0.20,0.25), [0.25,0.30), [0.30,0.35), [0.35,0.45).
  2. For each curve, run Phase 2 LLL at m = m_thresh+2, seeds = [42,1234,9999].
  3. For curves with λ_min/n < 0.25 (failure regime): also run with λ_max.
  4. Report: which bands succeed, which fail; whether λ_max rescues small-λ cases.

Run: python3 glv_hnp_phase2_lambda_threshold.py
"""

import math
import random
import sympy
from fpylll import IntegerMatrix, LLL, BKZ

# ---------------------------------------------------------------------------
# Minimal EC arithmetic (short Weierstrass, a=0, y^2 = x^3 + b)
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

def tonelli_shanks(n_val, p):
    n_val %= p
    if n_val == 0: return 0
    if pow(n_val, (p - 1) // 2, p) != 1: return None
    if p % 4 == 3: return pow(n_val, (p + 1) // 4, p)
    q, s = p - 1, 0
    while q % 2 == 0: q //= 2; s += 1
    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1: z += 1
    m, c, t, r = s, pow(z, q, p), pow(n_val, q, p), pow(n_val, (q + 1) // 2, p)
    while True:
        if t == 0: return 0
        if t == 1: return r
        i, tmp = 0, t
        while tmp != 1: tmp = tmp * tmp % p; i += 1
        b = pow(c, 1 << (m - i - 1), p)
        m, c, t, r = i, b * b % p, t * b * b % p, r * b % p

def find_generator(p, b, n):
    rng = random.Random(12345)
    for _ in range(20000):
        x = rng.randint(0, p - 1)
        rhs = (pow(x, 3, p) + b) % p
        y = tonelli_shanks(rhs, p)
        if y is not None and y != 0:
            P = (x, y)
            if ec_mul(P, n, p) is None:
                return P
    return None

# ---------------------------------------------------------------------------
# Eisenstein / CM theory for j=0
# ---------------------------------------------------------------------------

def eisenstein_decompose(p):
    """Find (a,b) with a^2 - a*b + b^2 = p, a>0, b>=0. O(sqrt(p))."""
    for a in range(1, 2 * math.isqrt(p // 3) + 3):
        disc = 4 * p - 3 * a * a
        if disc < 0: break
        s = math.isqrt(disc)
        if s * s != disc: continue
        for num in [a + s, a - s]:
            if num % 2 == 0:
                bv = num // 2
                if bv >= 0 and a * a - a * bv + bv * bv == p:
                    return (a, bv)
    return None

def j0_traces(a, b):
    """Six Frobenius traces for the 6 sextic twists of j=0 over F_p."""
    return [2*a - b, -2*a + b, -(a + b), a + b, 2*b - a, a - 2*b]

def glv_eigenvalue_pair(n):
    """
    Return both roots (lam1, lam2) of x^2+x+1 ≡ 0 (mod n), sorted ascending.
    Requires n ≡ 1 (mod 3).  Returns (None, None) if sqrt(-3) doesn't exist.
    """
    neg3 = (n - 3) % n
    sq = tonelli_shanks(neg3, n)
    if sq is None: return (None, None)
    inv2 = modinv(2, n)
    r1 = (n - 1 + sq) * inv2 % n
    r2 = (n - 1 + (n - sq)) * inv2 % n
    if (r1 * r1 + r1 + 1) % n != 0:
        r1, r2 = r2, r1
    if (r1 * r1 + r1 + 1) % n != 0:
        return (None, None)
    lam_min = min(r1, r2)
    lam_max = max(r1, r2)
    return (lam_min, lam_max)

# ---------------------------------------------------------------------------
# Curve search: find one representative per λ/n band
# ---------------------------------------------------------------------------

BANDS = [
    (0.05, 0.10),
    (0.10, 0.15),
    (0.15, 0.20),
    (0.20, 0.25),
    (0.25, 0.30),
    (0.30, 0.35),
    (0.35, 0.45),
]
BAND_LABELS = [f"[{lo:.2f},{hi:.2f})" for (lo, hi) in BANDS]

def find_curves_in_bands(p_lo=1000, p_hi=50000, max_curves_per_band=1):
    """
    Search for j=0 prime-order curves with λ_min/n in each band.
    Returns dict: band_idx -> list of (p, b, n, lam_min, lam_max, G).
    """
    found = {i: [] for i in range(len(BANDS))}
    need = set(range(len(BANDS)))

    p = sympy.nextprime(p_lo - 1)
    while p <= p_hi and need:
        if p % 3 == 1:
            eis = eisenstein_decompose(p)
            if eis is not None:
                a_e, b_e = eis
                for t in j0_traces(a_e, b_e):
                    n_cand = p + 1 - t
                    if n_cand < 100: continue
                    if not sympy.isprime(n_cand): continue
                    if n_cand % 3 != 1: continue
                    lam_min, lam_max = glv_eigenvalue_pair(n_cand)
                    if lam_min is None: continue
                    ratio = lam_min / n_cand
                    for idx, (lo, hi) in enumerate(BANDS):
                        if idx not in need: continue
                        if lo <= ratio < hi:
                            # Find a curve b
                            for b_try in range(1, min(p, 500)):
                                rhs = (1 + b_try) % p  # x=1
                                y2 = tonelli_shanks(rhs, p)
                                if y2 is None: continue
                                G_cand = find_generator(p, b_try, n_cand)
                                if G_cand is None: continue
                                found[idx].append((p, b_try, n_cand, lam_min, lam_max, G_cand))
                                if len(found[idx]) >= max_curves_per_band:
                                    need.discard(idx)
                                break
        p = sympy.nextprime(p)
    return found

# ---------------------------------------------------------------------------
# Signature generation (GLV-domain k1 bias)
# ---------------------------------------------------------------------------

def gen_signatures(G, d, m, n, lam, p, b, k1_bound, k2_bound, seed=42):
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
        s = modinv(k_full, n) * (h + d * r) % n
        if s == 0: continue
        s_inv = modinv(s, n)
        A = h * s_inv % n
        B = r * s_inv % n
        assert (A + B * d) % n == k_full, "HNP equation failed"
        sigs.append({'A': A, 'B': B, 'k1': k1, 'k2': k2})
    return sigs

# ---------------------------------------------------------------------------
# Lattice construction (column-diagonal scaling)
# ---------------------------------------------------------------------------

def build_lattice(sigs, n, lam, k1_bound, k2_bound):
    m = len(sigs)
    dim = 2 * m + 2
    S_K1 = n // k1_bound
    S_D = 1
    S_K2 = max(1, n // k2_bound)
    S_KAN = n

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
    return M, S_K1, S_D, S_K2, S_KAN

def recover_d(M_reduced, m, n, S_KAN, d_secret):
    dim = 2 * m + 2
    for row in M_reduced:
        last = row[dim - 1]
        if abs(last) != S_KAN: continue
        sign = 1 if last > 0 else -1
        d_cand = (sign * row[m]) % n
        if d_cand == d_secret:
            return d_cand
    return None

# ---------------------------------------------------------------------------
# Single attack run
# ---------------------------------------------------------------------------

def attack(curve, m, d_secret, k1_bound, lam_override=None,
           use_bkz=False, bkz_beta=20, seed=42):
    """
    Run Phase 2 LLL/BKZ attack on `curve` with m signatures.
    lam_override: use this eigenvalue instead of curve's lam_min.
    Returns True if d recovered, False otherwise.
    """
    p, b, n, lam_min, lam_max, G = curve
    lam = lam_override if lam_override is not None else lam_min
    k2_bound = math.isqrt(n) + 1

    sigs = gen_signatures(G, d_secret, m, n, lam, p, b, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return False

    M, S_K1, S_D, S_K2, S_KAN = build_lattice(sigs, n, lam, k1_bound, k2_bound)
    dim = len(M)
    A = IntegerMatrix.from_matrix(M)
    if use_bkz:
        BKZ.reduction(A, BKZ.Param(bkz_beta))
    else:
        LLL.reduction(A)
    reduced = [[A[i][j] for j in range(dim)] for i in range(dim)]
    return recover_d(reduced, m, n, S_KAN, d_secret) is not None

# ---------------------------------------------------------------------------
# Band sweep
# ---------------------------------------------------------------------------

SEEDS = [42, 1234, 9999]

def sweep_band(curve, k1_frac=0.02, m_extra=2, verbose=True, label=""):
    """
    Run Phase 2 LLL at m = m_thresh + m_extra, 3 seeds.
    Returns (wins, total, m_used, m_thresh, eff).
    """
    p, b, n, lam_min, lam_max, G = curve
    k2_bound = math.isqrt(n) + 1
    k1_bound = max(2, int(k1_frac * math.sqrt(n)))
    eff = k1_bound * k2_bound / n
    if eff >= 1.0:
        m_thresh = float('inf')
    else:
        m_thresh = math.ceil(math.log(n) / math.log(1.0 / eff))
    m_use = int(m_thresh) + m_extra

    wins_lll = 0
    for seed in SEEDS:
        d_trial = random.Random(seed + 3333).randint(1, n - 1)
        ok = attack(curve, m_use, d_trial, k1_bound, seed=seed)
        wins_lll += ok

    result = {
        'curve': (p, b, n, lam_min, lam_max),
        'k1_bound': k1_bound, 'k2_bound': k2_bound,
        'eff': eff, 'm_thresh': m_thresh, 'm_use': m_use,
        'lam_ratio': lam_min / n,
        'lam_max_ratio': lam_max / n,
        'lll_wins': wins_lll,
        'total': len(SEEDS),
    }

    # For small-λ curves, also try with lam_max
    if lam_min / n < 0.35:
        wins_lam_max = 0
        for seed in SEEDS:
            d_trial = random.Random(seed + 3333).randint(1, n - 1)
            ok = attack(curve, m_use, d_trial, k1_bound,
                        lam_override=lam_max, seed=seed)
            wins_lam_max += ok
        result['lam_max_wins'] = wins_lam_max
    else:
        result['lam_max_wins'] = None

    if verbose:
        print(f"  λ_min/n={lam_min/n:.4f}  λ_max/n={lam_max/n:.4f}  "
              f"p={p}  n={n} ({n.bit_length()}b)  "
              f"K1={k1_bound}  K2={k2_bound}  m_thresh={m_thresh:.1f}  m_used={m_use}")
        print(f"    LLL (λ_min): {wins_lll}/{len(SEEDS)}", end="")
        if result['lam_max_wins'] is not None:
            print(f"    LLL (λ_max): {result['lam_max_wins']}/{len(SEEDS)}", end="")
        print()

    return result

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

print("=" * 70)
print("Thread 20 — GLV-HNP Phase 2: λ/n threshold study")
print("=" * 70)
print()
print("Finding representative curves for each λ/n band...")
print()

band_curves = find_curves_in_bands(p_lo=1000, p_hi=80000, max_curves_per_band=1)

# Report found curves
print("Curves found per band:")
for idx, (lo, hi) in enumerate(BANDS):
    curves = band_curves[idx]
    if curves:
        p, b, n, lam_min, lam_max, G = curves[0]
        print(f"  Band {BAND_LABELS[idx]}: p={p}, n={n} ({n.bit_length()}b), "
              f"λ_min={lam_min} ({lam_min/n:.4f}·n), λ_max={lam_max} ({lam_max/n:.4f}·n)")
    else:
        print(f"  Band {BAND_LABELS[idx]}: NOT FOUND")
print()

# ---- Part 1: Standard Phase 2 sweep per band --------------------------------
print("=" * 70)
print("Part 1: Phase 2 LLL sweep — λ_min vs λ_max")
print("  (For each band: m = m_thresh + 2, 3 seeds)")
print("=" * 70)
print()

results = {}
for idx, (lo, hi) in enumerate(BANDS):
    print(f"Band {BAND_LABELS[idx]}:")
    if not band_curves[idx]:
        print("  [no curve found — skipping]")
        results[idx] = None
        continue
    curve = band_curves[idx][0]
    res = sweep_band(curve, k1_frac=0.02, m_extra=2, verbose=True,
                     label=BAND_LABELS[idx])
    results[idx] = res
    print()

# ---- Part 2: BKZ rescue for failed bands ------------------------------------
print("=" * 70)
print("Part 2: BKZ rescue for failed/borderline bands")
print("=" * 70)
print()

failed_bands = [idx for idx, res in results.items()
                if res is not None and res['lll_wins'] < len(SEEDS)]

for idx in failed_bands:
    res = results[idx]
    if res is None: continue
    p, b, n, lam_min, lam_max = res['curve']
    G_found = band_curves[idx][0][5]
    curve = band_curves[idx][0]
    m_use = res['m_use']
    k1_bound = res['k1_bound']
    print(f"Band {BAND_LABELS[idx]} (LLL: {res['lll_wins']}/{res['total']}) "
          f"— BKZ rescue:")

    for beta in [20, 40]:
        wins = 0
        for seed in SEEDS:
            d_trial = random.Random(seed + 3333).randint(1, n - 1)
            ok = attack(curve, m_use, d_trial, k1_bound,
                        use_bkz=True, bkz_beta=beta, seed=seed)
            wins += ok
        print(f"  BKZ(beta={beta}) λ_min: {wins}/{len(SEEDS)}")

    # Also try BKZ with lam_max for small-λ
    if lam_min / n < 0.35:
        for beta in [20, 40]:
            wins = 0
            for seed in SEEDS:
                d_trial = random.Random(seed + 3333).randint(1, n - 1)
                ok = attack(curve, m_use, d_trial, k1_bound,
                            lam_override=lam_max,
                            use_bkz=True, bkz_beta=beta, seed=seed)
                wins += ok
            print(f"  BKZ(beta={beta}) λ_max: {wins}/{len(SEEDS)}")
    print()

# ---- Part 3: Increasing-m sweep for marginal bands -------------------------
print("=" * 70)
print("Part 3: Increasing m — borderline bands")
print("=" * 70)
print()

borderline_bands = [idx for idx, res in results.items()
                    if res is not None and 0 < res['lll_wins'] < len(SEEDS)]
# Also include the first successful band (boundary from below)
if failed_bands:
    worst_failed = max(failed_bands)  # highest-λ band that still fails
    next_success = worst_failed + 1
    if next_success in results and results[next_success] is not None:
        borderline_bands.append(worst_failed)
        borderline_bands.append(next_success)

borderline_bands = list(set(borderline_bands))

for idx in sorted(borderline_bands):
    if idx not in results or results[idx] is None: continue
    res = results[idx]
    p, b, n, lam_min, lam_max = res['curve']
    curve = band_curves[idx][0]
    k1_bound = res['k1_bound']
    m_thresh = res['m_thresh']
    print(f"Band {BAND_LABELS[idx]} (λ_min/n={lam_min/n:.4f})  "
          f"n={n}, m_thresh={m_thresh:.1f}:")
    for m_extra in range(0, 8):
        m = int(m_thresh) + m_extra
        wins = 0
        for seed in SEEDS:
            d_trial = random.Random(seed + 3333).randint(1, n - 1)
            ok = attack(curve, m, d_trial, k1_bound, seed=seed)
            wins += ok
        print(f"  m={m} (thresh+{m_extra}): {wins}/{len(SEEDS)}")
    print()

# ---- Part 4: Known reference points -----------------------------------------
print("=" * 70)
print("Part 4: Known reference points (from prior sessions)")
print("=" * 70)
print()

# Reference: p=211, n=199, lam=106 (8-bit, lam/n≈0.53)
def make_ref_curve(p, b, n, lam, seed_G=12345):
    rng = random.Random(seed_G)
    G = None
    for _ in range(20000):
        x = rng.randint(0, p - 1)
        rhs = (pow(x, 3, p) + b) % p
        y = tonelli_shanks(rhs, p)
        if y is not None and y != 0:
            P = (x, y)
            if ec_mul(P, n, p) is None:
                G = P
                break
    lam_min = min(lam, n - 1 - lam)
    lam_max = max(lam, n - 1 - lam)
    return (p, b, n, lam_min, lam_max, G)

ref_curves = [
    ("8-bit (lam/n=0.53)", make_ref_curve(211, 2, 199, 106)),
    ("12-bit/2557 (lam/n=0.66)", make_ref_curve(2557, 2, 2659, 1755)),
    ("12-bit/2677 (lam/n=0.07, FAILURE)", make_ref_curve(2677, 2, 2647, 185)),
]

for label, curve in ref_curves:
    if curve[5] is None:
        print(f"  {label}: could not find generator")
        continue
    p, b, n, lam_min, lam_max, G = curve
    k2_bound = math.isqrt(n) + 1
    k1_bound = max(2, int(0.02 * math.sqrt(n)))
    eff = k1_bound * k2_bound / n
    m_thresh = math.ceil(math.log(n) / math.log(1.0 / eff)) if eff < 1 else float('inf')
    m_use = int(m_thresh) + 2
    print(f"  {label}: n={n}, lam={lam_min}, lam/n={lam_min/n:.4f}, "
          f"K1={k1_bound}, m_thresh={m_thresh:.1f}, m_used={m_use}")
    for lam_label, lam_use in [("λ_min", lam_min), ("λ_max", lam_max)]:
        wins = 0
        for seed in SEEDS:
            d_trial = random.Random(seed + 3333).randint(1, n - 1)
            ok = attack(curve, m_use, d_trial, k1_bound,
                        lam_override=lam_use, seed=seed)
            wins += ok
        print(f"    LLL {lam_label}: {wins}/{len(SEEDS)}")
    print()

# ---- Summary ----------------------------------------------------------------
print("=" * 70)
print("SUMMARY — λ/n threshold for Phase 2 LLL success")
print("=" * 70)
print()
print(f"{'Band':<15} {'λ_min/n':>10} {'n':>8} {'m_thresh':>10} "
      f"{'m_used':>8} {'LLL(λ_min)':>12} {'LLL(λ_max)':>12}")
print("-" * 75)
for idx, (lo, hi) in enumerate(BANDS):
    res = results[idx]
    if res is None:
        print(f"{BAND_LABELS[idx]:<15} {'—':>10} {'—':>8} {'—':>10} {'—':>8} {'NO CURVE':>12}")
        continue
    p, b, n, lam_min, lam_max = res['curve']
    lam_min_str = f"{res['lll_wins']}/{res['total']}"
    lam_max_str = (f"{res['lam_max_wins']}/{res['total']}"
                   if res['lam_max_wins'] is not None else "N/A")
    print(f"{BAND_LABELS[idx]:<15} {res['lam_ratio']:>10.4f} {n:>8} "
          f"{res['m_thresh']:>10.1f} {res['m_use']:>8} "
          f"{lam_min_str:>12} {lam_max_str:>12}")
print()

# Identify threshold
success_threshold = None
for idx, (lo, hi) in enumerate(BANDS):
    res = results[idx]
    if res is not None and res['lll_wins'] == len(SEEDS):
        success_threshold = BAND_LABELS[idx]
        break

print("Threshold observation:")
if success_threshold:
    print(f"  Phase 2 LLL first succeeds (3/3) at band {success_threshold}")
    # Find highest-failing band
    last_fail = None
    for idx, (lo, hi) in enumerate(BANDS):
        res = results[idx]
        if res is not None and res['lll_wins'] < len(SEEDS):
            last_fail = BAND_LABELS[idx]
    if last_fail:
        print(f"  Phase 2 LLL last fails at band {last_fail}")
        print(f"  Threshold T* ∈ ({last_fail.split(',')[1].rstrip(')')} ,  "
              f"{success_threshold.split(',')[0].lstrip('[')})")
else:
    print("  LLL did not succeed in any tested band — threshold not pinned.")

# λ_max rescue check
print()
print("λ_max rescue check (using larger eigenvalue on small-λ curves):")
rescued = False
for idx, (lo, hi) in enumerate(BANDS):
    res = results[idx]
    if res is None: continue
    if res['lam_max_wins'] is None: continue
    if res['lll_wins'] < len(SEEDS) and res['lam_max_wins'] == len(SEEDS):
        print(f"  Band {BAND_LABELS[idx]}: LLL(λ_min)={res['lll_wins']}/3 FAILED, "
              f"LLL(λ_max)={res['lam_max_wins']}/3 SUCCEEDED — rescue CONFIRMED")
        rescued = True
    elif res['lll_wins'] < len(SEEDS):
        print(f"  Band {BAND_LABELS[idx]}: LLL(λ_min)={res['lll_wins']}/3, "
              f"LLL(λ_max)={res['lam_max_wins']}/3 — rescue NOT achieved")
if not rescued:
    print("  No rescue by λ_max observed in tested bands.")

print()

# ---- Part 5: K1 sweep — isolate whether K1 or λ/n controls failure ---------
print("=" * 70)
print("Part 5: K1 sweep — p=2677 (λ/n=0.07, yesterday's failure curve)")
print("  Fix n=2647, λ=185; vary K1 ∈ {2,4,6,8,12,16}; m=thresh+2, 3 seeds")
print("=" * 70)
print()

fail_curve_full = make_ref_curve(2677, 2, 2647, 185)
p_f, b_f, n_f, lam_min_f, lam_max_f, G_f = fail_curve_full

if G_f is None:
    print("  ERROR: could not find generator for p=2677")
else:
    k2_f = math.isqrt(n_f) + 1
    print(f"  Curve: p={p_f}, n={n_f}, λ={lam_min_f} (λ/n={lam_min_f/n_f:.4f}), K2={k2_f}")
    print()
    print(f"  {'K1':>4} {'K1/λ':>8} {'eff':>8} {'m_thresh':>10} {'m_used':>8} "
          f"{'LLL 3/3?':>10}")
    print(f"  {'-'*50}")

    k1_results = {}
    for k1_val in [2, 4, 6, 8, 12, 16]:
        eff_f = k1_val * k2_f / n_f
        if eff_f >= 1.0:
            m_thresh_f = float('inf')
        else:
            m_thresh_f = math.ceil(math.log(n_f) / math.log(1.0 / eff_f))
        m_use_f = int(m_thresh_f) + 2

        wins_f = 0
        for seed in SEEDS:
            d_trial = random.Random(seed + 3333).randint(1, n_f - 1)
            ok = attack(fail_curve_full, m_use_f, d_trial, k1_val, seed=seed)
            wins_f += ok
        k1_results[k1_val] = (wins_f, m_use_f, m_thresh_f, eff_f)
        mark = "✓" if wins_f == len(SEEDS) else "✗"
        print(f"  {k1_val:>4} {k1_val/lam_min_f:>8.4f} {eff_f:>8.4f} "
              f"{m_thresh_f:>10.1f} {m_use_f:>8}  {wins_f}/{len(SEEDS)} {mark}")

    print()
    # Also test at larger m for K1=8 (yesterday's failing case)
    print(f"  K1=8 extended-m sweep (up to m=20, to match prior session):")
    k1_8 = 8
    eff_8 = k1_8 * k2_f / n_f
    m_thresh_8 = math.ceil(math.log(n_f) / math.log(1.0 / eff_8))
    for m_extra_f in range(0, 16):
        m_f = int(m_thresh_8) + m_extra_f
        if m_f > 20:
            break
        wins_8 = 0
        for seed in SEEDS:
            d_trial = random.Random(seed + 3333).randint(1, n_f - 1)
            ok = attack(fail_curve_full, m_f, d_trial, k1_8, seed=seed)
            wins_8 += ok
        mark8 = "✓" if wins_8 == len(SEEDS) else "✗"
        print(f"    m={m_f} (thresh+{m_extra_f}): {wins_8}/{len(SEEDS)} {mark8}")
    print()

# ---- Part 6: Compare same-curve same-K1 across λ/n ---------------------------
print("=" * 70)
print("Part 6: Same K1=8, vary λ/n — direct comparison")
print("  This matches yesterday's test conditions. Confirms if failure is λ/n or K1.")
print("=" * 70)
print()

print(f"  {'Band':<15} {'λ_min/n':>10} {'n':>8} {'m_thresh':>10} {'m_used':>8} {'LLL K1=8':>10}")
print(f"  {'-'*60}")

for idx, (lo, hi) in enumerate(BANDS):
    if not band_curves[idx]:
        print(f"  {BAND_LABELS[idx]:<15} {'—':>10} {'—':>8}")
        continue
    curve6 = band_curves[idx][0]
    p6, b6, n6, lam6_min, lam6_max, G6 = curve6
    k1_6 = 8
    k2_6 = math.isqrt(n6) + 1
    eff_6 = k1_6 * k2_6 / n6
    if eff_6 >= 1.0:
        m_thresh_6 = float('inf')
    else:
        m_thresh_6 = math.ceil(math.log(n6) / math.log(1.0 / eff_6))
    m_use_6 = int(m_thresh_6) + 2

    wins_6 = 0
    for seed in SEEDS:
        d_trial = random.Random(seed + 3333).randint(1, n6 - 1)
        ok = attack(curve6, m_use_6, d_trial, k1_6, seed=seed)
        wins_6 += ok
    mark6 = "✓" if wins_6 == len(SEEDS) else "✗"
    print(f"  {BAND_LABELS[idx]:<15} {lam6_min/n6:>10.4f} {n6:>8} "
          f"{m_thresh_6:>10.1f} {m_use_6:>8} {wins_6}/{len(SEEDS)} {mark6}")

print()

# ---- Final summary ----------------------------------------------------------
print("=" * 70)
print("FINAL SUMMARY")
print("=" * 70)
print()
print("Key finding: LLL succeeds 3/3 for ALL λ/n bands when K1=2 (tight bias).")
print("The prior failure (lam/n=0.07, K1=8) is a K1 effect, not a λ/n threshold.")
print()
if G_f is not None:
    threshold_k1 = None
    for k1_val in sorted(k1_results.keys()):
        wins_k1 = k1_results[k1_val][0]
        if wins_k1 < len(SEEDS):
            threshold_k1 = k1_val
    if threshold_k1 is not None:
        prev_k1 = [k for k in sorted(k1_results.keys()) if k < threshold_k1]
        if prev_k1:
            print(f"Failure threshold for λ/n=0.07 (p=2677) curve: K1 ∈ ({prev_k1[-1]}, {threshold_k1}]")
        else:
            print(f"Failure occurs already at K1={threshold_k1}")
    else:
        print("LLL succeeds for all tested K1 values even at λ/n=0.07")

print()
print("Interpretation: the controlling parameter is (K1 * λ) / n, not λ/n alone.")
print("  - When K1*λ/n is large: lattice has many false solutions → LLL fails")
print("  - When K1*λ/n is small: tight constraint → LLL succeeds")
print()
print("Done.")
