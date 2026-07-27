"""
GLV-HNP Phase 2: λ/n threshold — controlled study at fixed bit length.

The initial sweep (glv_hnp_phase2_lambda_threshold.py) found a non-monotone
pattern that turned out to be a confound: small-n curves (n=13, n=37) fail
even at λ/n=0.23-0.27 due to boundary effects (tiny lattice). Meanwhile
n=73 with λ/n=0.11 passes.

This script controls for n size by looking for 10-12 bit curves (n ~ 1000-4000)
with λ/n ∈ {0.10, 0.15, 0.20, 0.25, 0.30} and repeating the LLL sweep.

The known failure point: n=2647 (12-bit), λ/n=0.07 — fails even BKZ(40).
The known pass point:   n=523969 (20-bit), λ/n=0.34 — LLL 3/3 at m=9.

Goal: locate the threshold at ~12-bit scale.

Run: python3 glv_hnp_phase2_lam_controlled.py
"""

import math
import random
import sympy
from fpylll import IntegerMatrix, LLL, BKZ

# ---- EC arithmetic ----------------------------------------------------------

def modinv(a, m): return pow(a, -1, m)

def ec_add(P, Q, p):
    if P is None: return Q
    if Q is None: return P
    x1, y1 = P; x2, y2 = Q
    if x1 == x2:
        if (y1 + y2) % p == 0: return None
        s = 3*x1*x1*modinv(2*y1,p) % p
    else:
        s = (y2-y1)*modinv(x2-x1,p) % p
    x3 = (s*s - x1 - x2) % p
    y3 = (s*(x1-x3) - y1) % p
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
    if pow(n,(p-1)//2,p) != 1: return None
    if p%4==3: return pow(n,(p+1)//4,p)
    q, s = p-1, 0
    while q%2==0: q//=2; s+=1
    z = 2
    while pow(z,(p-1)//2,p) != p-1: z+=1
    mv, c, t, r = s, pow(z,q,p), pow(n,q,p), pow(n,(q+1)//2,p)
    while True:
        if t==0: return 0
        if t==1: return r
        i, tmp = 0, t
        while tmp!=1: tmp=tmp*tmp%p; i+=1
        b = pow(c, 1<<(mv-i-1), p)
        mv, c, t, r = i, b*b%p, t*b*b%p, r*b%p

def find_generator(p, b, n, seed=12345):
    rng = random.Random(seed)
    for _ in range(100000):
        x = rng.randint(0, p-1)
        rhs = (pow(x,3,p)+b)%p
        y = tonelli_shanks(rhs, p)
        if y is not None and y!=0:
            P=(x,y)
            if ec_mul(P,n,p) is None: return P
    return None

# ---- CM theory ---------------------------------------------------------------

def eisenstein_decompose(p):
    for a in range(1, 2*math.isqrt(p//3)+3):
        disc = 4*p - 3*a*a
        if disc < 0: break
        s = math.isqrt(disc)
        if s*s != disc: continue
        for num in [a+s, a-s]:
            if num%2==0:
                b=num//2
                if b>=0 and a*a-a*b+b*b==p: return (a,b)
    return None

def j0_traces(a, b):
    return [2*a-b, -2*a+b, -(a+b), a+b, 2*b-a, a-2*b]

def glv_eigenvalue(n):
    sq = tonelli_shanks((n-3)%n, n)
    if sq is None: return None, None
    inv2 = modinv(2, n)
    r1 = (n-1+sq)*inv2%n; r2 = (n-1+(n-sq))*inv2%n
    if (r1*r1+r1+1)%n!=0: r1,r2=r2,r1
    if (r1*r1+r1+1)%n!=0: return None, None
    lam = min(r1, r2)
    return lam, n-1-lam

def find_curve_b(p, n, seed=999):
    rng = random.Random(seed)
    for b in range(1, min(p, 100)):
        for x in range(min(p, 500)):
            rhs = (pow(x,3,p)+b)%p
            y = tonelli_shanks(rhs, p)
            if y is not None and y!=0:
                P=(x,y)
                if ec_mul(P,n,p) is None: return b
                break
    return None

# ---- Curve discovery: 10-12 bit curves with target λ/n ----------------------

def discover_controlled(target_ratios, n_min=512, n_max=8192, p_max=1000000):
    """Find 10-12 bit curves with λ/n near each target ratio."""
    found = {}
    remaining = set(target_ratios)
    p = sympy.nextprime(n_min)
    while p < p_max and remaining:
        if p % 3 == 1:
            eis = eisenstein_decompose(p)
            if eis is not None:
                a_e, b_e = eis
                for t in j0_traces(a_e, b_e):
                    n_cand = p + 1 - t
                    if not (n_min <= n_cand <= n_max): continue
                    if not sympy.isprime(n_cand): continue
                    if n_cand % 3 != 1: continue
                    lam, _ = glv_eigenvalue(n_cand)
                    if lam is None: continue
                    ratio = lam / n_cand
                    for target in list(remaining):
                        if abs(ratio - target) < 0.025:
                            b_param = find_curve_b(p, n_cand)
                            if b_param is None: continue
                            G = find_generator(p, b_param, n_cand)
                            if G is None: continue
                            found[target] = (p, b_param, n_cand, lam, G)
                            remaining.discard(target)
                            bits = n_cand.bit_length()
                            print(f"  [bin {target:.2f}] p={p}, n={n_cand} ({bits}b), "
                                  f"λ={lam}, λ/n={ratio:.4f}")
                            break
        p = sympy.nextprime(p)
    if remaining:
        print(f"  WARNING: no curves for bins: {sorted(remaining)}")
    return found

# ---- Lattice / LLL -----------------------------------------------------------

def gen_signatures(G, d, m, n, lam, p, b, k1b, k2b, seed):
    rng = random.Random(seed)
    sigs = []
    attempts = 0
    while len(sigs) < m and attempts < 500000:
        attempts += 1
        k1 = rng.randint(0, k1b-1)
        k2 = rng.randint(0, k2b-1)
        kf = (k1 + lam*k2) % n
        if kf==0: continue
        R = ec_mul(G, kf, p)
        if R is None: continue
        r = R[0]%n
        if r==0: continue
        h = rng.randint(0, n-1)
        s = modinv(kf,n)*(h+d*r)%n
        if s==0: continue
        si = modinv(s,n)
        A = h*si%n; B = r*si%n
        sigs.append({'A':A,'B':B,'k1':k1,'k2':k2})
    return sigs

def build_lattice(sigs, n, lam, k1b, k2b):
    m = len(sigs); dim = 2*m+2
    SK1 = max(1, n//k1b); SK2 = max(1, n//k2b); SKANN = n
    M = [[0]*dim for _ in range(dim)]
    for i in range(m): M[i][i] = n*SK1
    for i in range(m): M[m][i] = sigs[i]['B']*SK1
    M[m][m] = 1
    for i in range(m):
        M[m+1+i][i] = -lam*SK1
        M[m+1+i][m+1+i] = SK2
    for i in range(m): M[2*m+1][i] = sigs[i]['A']*SK1
    M[2*m+1][dim-1] = SKANN
    return M, SKANN

def recover(M_reduced, m, n, SKANN, d_secret):
    dim = 2*m+2
    for row in M_reduced:
        last = row[dim-1]
        if abs(last)!=SKANN: continue
        sign = 1 if last>0 else -1
        dc = (sign*row[m])%n
        if dc==d_secret: return True
    return False

def run_exp(curve, m, d, k1b, seed, bkz_beta=None):
    p, b, n, lam, G = curve
    k2b = math.isqrt(n)+1
    sigs = gen_signatures(G, d, m, n, lam, p, b, k1b, k2b, seed)
    if len(sigs)<m: return False
    M, SKANN = build_lattice(sigs, n, lam, k1b, k2b)
    A = IntegerMatrix.from_matrix(M)
    if bkz_beta:
        BKZ.reduction(A, BKZ.Param(bkz_beta))
    else:
        LLL.reduction(A)
    reduced = [[A[i][j] for j in range(len(M))] for i in range(len(M))]
    return recover(reduced, m, n, SKANN, d)

def sweep(label, curve, k1b, m_range, seeds, bkz_beta=None):
    p, b, n, lam, G = curve
    k2b = math.isqrt(n)+1
    eff = k1b*k2b/n
    mt = math.ceil(math.log(n)/math.log(1/eff)) if eff<1 else float('inf')
    algo = f"BKZ({bkz_beta})" if bkz_beta else "LLL"
    print(f"  [{algo}] {label}")
    print(f"    n={n} ({n.bit_length()}b), λ/n={lam/n:.4f}, K1={k1b}, K2={k2b}, "
          f"eff={eff:.4f}, m_thresh≈{mt:.1f}")
    first_full = None
    results = {}
    for m in m_range:
        wins = sum(
            run_exp(curve, m, random.Random(seed+5555).randint(1,n-1), k1b, seed,
                    bkz_beta=bkz_beta)
            for seed in seeds
        )
        results[m] = wins
        mk = "*" if m>=mt else " "
        print(f"    {mk}m={m}: {wins}/{len(seeds)}")
        if wins==len(seeds) and first_full is None:
            first_full = m
    return results, first_full

# ---- Main -------------------------------------------------------------------

print("="*70)
print("GLV-HNP Phase 2 — λ/n threshold controlled study (10-12 bit curves)")
print("="*70)
print()

SEEDS = [42, 1234, 9999]
TARGET_RATIOS = [0.10, 0.15, 0.20, 0.25, 0.30]

print("Step 1: Discover 10-12 bit curves for each λ/n bin")
print("-"*70)
curves = discover_controlled(TARGET_RATIOS, n_min=512, n_max=8192)

print()
print("Step 2: LLL sweep")
print("="*70)

summary = []
for target in sorted(curves.keys()):
    curve = curves[target]
    p, b, n, lam, G = curve
    k1b = max(2, int(0.05*math.sqrt(n)))
    label = f"bin {target:.2f}: n={n}, λ/n={lam/n:.4f}"
    print(f"\n{'='*70}")
    print(f"Curve: {label}")
    res_lll, first_lll = sweep(label, curve, k1b, range(4, 16), SEEDS)
    first_bkz = None
    if first_lll is None:
        print("  LLL failed — trying BKZ(40)")
        res_bkz, first_bkz = sweep(label, curve, k1b, range(4, 16), SEEDS, bkz_beta=40)
    summary.append((lam/n, target, label, first_lll, first_bkz))

print()
print("="*70)
print("SUMMARY — controlled 10-12 bit curves")
print("="*70)
print(f"{'λ/n':>8}  {'n (bits)':>10}  {'LLL 3/3':>10}  {'BKZ(40) 3/3':>12}  Status")
print("-"*70)
print(f"{'0.0699':>8}  {'2647 (12b)':>10}  {'never':>10}  {'never':>12}  FAIL (known)")
for lam_ratio, target, label, fl, fb in sorted(summary):
    # extract n from label
    n_str = label.split('n=')[1].split(',')[0]
    n_val = int(n_str)
    lll_s = f"m={fl}" if fl else "never"
    bkz_s = f"m={fb}" if fb else ("n/a" if fl else "never")
    status = "PASS" if fl else ("PASS(BKZ)" if fb else "FAIL")
    print(f"{lam_ratio:>8.4f}  {n_val:>4} ({n_val.bit_length():2}b)  "
          f"{lll_s:>10}  {bkz_s:>12}  {status}")
print(f"{'0.3396':>8}  {'523969 (20b)':>10}  {'m=9 (known)':>10}  {'n/a':>12}  PASS (known)")

print()
# Interpret threshold
fails = [r for r,_,_,fl,fb in sorted(summary) if fl is None and fb is None]
passes = [r for r,_,_,fl,_ in sorted(summary) if fl is not None]
if fails and passes:
    max_fail = max(fails)
    min_pass = min(passes)
    print(f"Threshold at ~12-bit scale: FAIL for λ/n ≤ {max_fail:.4f}, "
          f"PASS for λ/n ≥ {min_pass:.4f}")
elif not fails:
    print(f"No failures at 10-12b scale down to λ/n = {min(r for r,*_ in summary):.4f}")
    print("Threshold is below tested range — further refinement needed.")
else:
    print("All curves failed — threshold above tested range.")

print("\nDone.")
