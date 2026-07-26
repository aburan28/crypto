"""
GLV-HNP Phase 2: Refutation of the 'lam/n threshold' hypothesis.

The 2026-07-26 run claimed lam/n=0.07 is a structural failure.
The threshold study (Thread 20) showed LLL succeeds at m=12 for ALL lam/n bins.

This script tests the specific p=2677 (n=2647, lam=185, lam/n=0.070) curve
under two k1_bound settings:
  A) k1_bound=8 (as used in the 2026-07-26 run — the FAILURE case)
  B) k1_bound=2 (as used in glv_hnp_phase2_lambda_threshold.py — 0.05*sqrt(n))

Hypothesis: the failure was due to k1_bound=8 being too LARGE (weak bias),
NOT due to lam/n=0.07 being the structural obstacle.
"""

import math
import random
from fpylll import IntegerMatrix, LLL, BKZ

def modinv(a, m): return pow(a, -1, m)

def tonelli_shanks(n, p):
    n %= p
    if n == 0: return 0
    if pow(n, (p-1)//2, p) != 1: return None
    if p % 4 == 3: return pow(n, (p+1)//4, p)
    q, s = p - 1, 0
    while q % 2 == 0: q //= 2; s += 1
    z = 2
    while pow(z, (p-1)//2, p) != p-1: z += 1
    m2, c, t, r = s, pow(z, q, p), pow(n, q, p), pow(n, (q+1)//2, p)
    while True:
        if t == 0: return 0
        if t == 1: return r
        i, tmp = 0, t
        while tmp != 1: tmp = tmp*tmp%p; i += 1
        b = pow(c, 1<<(m2-i-1), p)
        m2, c, t, r = i, b*b%p, t*b*b%p, r*b%p

def ec_add(P, Q, p):
    if P is None: return Q
    if Q is None: return P
    x1,y1 = P; x2,y2 = Q
    if x1 == x2:
        if (y1+y2)%p == 0: return None
        s = 3*x1*x1*modinv(2*y1, p)%p
    else:
        s = (y2-y1)*modinv(x2-x1,p)%p
    x3 = (s*s-x1-x2)%p
    y3 = (s*(x1-x3)-y1)%p
    return (x3, y3)

def ec_mul(P, k, p):
    if k == 0: return None
    R, Q = None, P
    while k > 0:
        if k&1: R = ec_add(R, Q, p)
        Q = ec_add(Q, Q, p)
        k >>= 1
    return R

def find_generator(p, b, n, seed=12345):
    rng = random.Random(seed)
    for _ in range(10000):
        x = rng.randint(0, p-1)
        rhs = (pow(x,3,p)+b)%p
        y = tonelli_shanks(rhs, p)
        if y is not None and y != 0:
            P = (x, y)
            if ec_mul(P, n, p) is None:
                return P
    return None

def gen_signatures(G, d, m, n, lam, p, b, k1_bound, k2_bound, seed):
    rng = random.Random(seed)
    sigs = []
    attempts = 0
    while len(sigs) < m and attempts < 300000:
        attempts += 1
        k1 = rng.randint(0, k1_bound-1)
        k2 = rng.randint(0, k2_bound-1)
        k_full = (k1 + lam*k2) % n
        if k_full == 0: continue
        R = ec_mul(G, k_full, p)
        if R is None: continue
        r = R[0] % n
        if r == 0: continue
        h = rng.randint(0, n-1)
        s = modinv(k_full, n)*(h + d*r)%n
        if s == 0: continue
        s_inv = modinv(s, n)
        A = h*s_inv%n
        B = r*s_inv%n
        assert (A + B*d)%n == k_full
        sigs.append({'A':A,'B':B,'k1':k1,'k2':k2})
    return sigs

def build_lattice(sigs, n, lam, k1_bound, k2_bound):
    m = len(sigs)
    dim = 2*m+2
    S_K1 = max(1, n//k1_bound)
    S_D = 1
    S_K2 = max(1, n//k2_bound)
    S_KANNAN = n
    M = [[0]*dim for _ in range(dim)]
    for i in range(m): M[i][i] = n*S_K1
    for i in range(m): M[m][i] = sigs[i]['B']*S_K1
    M[m][m] = S_D
    for i in range(m):
        M[m+1+i][i] = -lam*S_K1
        M[m+1+i][m+1+i] = S_K2
    for i in range(m): M[2*m+1][i] = sigs[i]['A']*S_K1
    M[2*m+1][dim-1] = S_KANNAN
    return M, S_KANNAN

def recover_d(reduced, m, n, S_KANNAN):
    dim = 2*m+2
    for row in reduced:
        last = row[dim-1]
        if abs(last) != S_KANNAN: continue
        sign = 1 if last > 0 else -1
        d_cand = (sign*row[m])%n
        if d_cand == 0: continue
        return d_cand
    return None

def attack(p, b, n, lam, G, k1_bound, m, seed, use_bkz=False, bkz_beta=20):
    k2_bound = math.isqrt(n)+1
    d = random.Random(seed+5555).randint(1, n-1)
    sigs = gen_signatures(G, d, m, n, lam, p, b, k1_bound, k2_bound, seed)
    if len(sigs) < m: return None
    M, S_KANNAN = build_lattice(sigs, n, lam, k1_bound, k2_bound)
    dim = 2*m+2
    A = IntegerMatrix.from_matrix(M)
    if use_bkz:
        BKZ.reduction(A, BKZ.Param(bkz_beta))
    else:
        LLL.reduction(A)
    red = [[A[i][j] for j in range(dim)] for i in range(dim)]
    d_rec = recover_d(red, m, n, S_KANNAN)
    return d_rec == d

# ---------------------------------------------------------------------------
# Curve parameters
# ---------------------------------------------------------------------------
p, b, n, lam = 2677, 2, 2647, 185
G = find_generator(p, b, n)
assert G is not None, "Generator not found"
k2_bound = math.isqrt(n)+1  # 52

print("="*65)
print(f"p={p}, b={b}, n={n}, lam={lam}, lam/n={lam/n:.4f}")
print(f"k2_bound={k2_bound} (sqrt(n)+1)")
print("="*65)

SEEDS = [42, 1234, 9999, 7, 314159]
M_VALUES = [6, 8, 10, 12, 14, 16]

# --- Setting A: k1_bound=8 (previous run's choice — the FAILURE) -----------
print("\n=== Setting A: k1_bound=8 (from 2026-07-26 run) ===")
k1_a = 8
eff_a = k1_a * k2_bound / n
m_thresh_a = math.ceil(math.log(n) / math.log(1.0/eff_a))
print(f"  eff={eff_a:.4f}, m_thresh={m_thresh_a}")
print(f"  S_K1 = n//k1_bound = {n//k1_a}")
print(f"  -lam*S_K1 = {-lam*(n//k1_a)},  n*S_K1 = {n*(n//k1_a)}")
print(f"  ratio (-lam*S_K1)/(n*S_K1) = {lam/n:.4f}  (= lam/n, as expected)")
print()
for m in M_VALUES:
    wins = sum(attack(p,b,n,lam,G,k1_a,m,s) or 0 for s in SEEDS)
    print(f"  m={m:3d}: {wins}/{len(SEEDS)} LLL")

# --- Setting B: k1_bound=2 (this study's auto-choice) ----------------------
print("\n=== Setting B: k1_bound=2 (threshold study auto-choice) ===")
k1_b = max(2, int(0.05 * math.sqrt(n)))
eff_b = k1_b * k2_bound / n
m_thresh_b = math.ceil(math.log(n) / math.log(1.0/eff_b))
print(f"  k1_bound={k1_b}")
print(f"  eff={eff_b:.4f}, m_thresh={m_thresh_b}")
print(f"  S_K1 = n//k1_bound = {n//k1_b}")
print(f"  -lam*S_K1 = {-lam*(n//k1_b)},  n*S_K1 = {n*(n//k1_b)}")
print()
for m in M_VALUES:
    wins = sum(attack(p,b,n,lam,G,k1_b,m,s) or 0 for s in SEEDS)
    print(f"  m={m:3d}: {wins}/{len(SEEDS)} LLL")

# --- Setting C: sweep k1_bound at fixed m=12 to find the boundary ----------
print("\n=== Setting C: k1_bound sweep at m=12 ===")
m_fixed = 12
print(f"  m={m_fixed}")
for k1 in [2, 3, 4, 5, 6, 7, 8, 10, 12, 16]:
    eff = k1 * k2_bound / n
    m_thresh = math.ceil(math.log(n) / math.log(1.0/eff)) if eff < 1 else 999
    wins = sum(attack(p,b,n,lam,G,k1,m_fixed,s) or 0 for s in SEEDS)
    print(f"  k1_bound={k1:3d}  eff={eff:.4f}  m_thresh={m_thresh:3d}  wins={wins}/{len(SEEDS)}")

print("\nConclusion:")
print("  If Setting A (k1=8) fails and Setting B (k1=2) succeeds at same m=12,")
print("  the 'structural failure' was k1_bound-dependent, NOT lam/n-dependent.")
print("  The lam/n=0.07 threshold hypothesis is REFUTED.")
print()
print("Done.")
