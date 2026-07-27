"""
GLV-HNP Phase 2: follow-up to lambda_threshold.py findings.

The threshold sweep found ALL bins (lam/n ∈ [0.04, 0.50)) succeed with
K1=4, eff≈0.10.  The original failure (p=2677, n=2647, lam=185, K1=8)
used different K1.

This script investigates:

(A) Direct re-test of the failure curve under K1=4 (same eff as success cases).
    If it now PASSES, the failure was K1-dependent, not λ/n-dependent.

(B) Very small λ/n: extend bin 0 down to [0.02, 0.04) to see if the attack
    can still succeed with eff≈0.10.

(C) Root choice: re-run the failure curve (p=2677, n=2647) using the LARGER
    root (lam=2461, lam/n=0.930) to test whether root choice matters.

(D) K1 sweep on the failure curve: test K1 ∈ {4,6,8,12} to locate the
    K1 threshold above which LLL fails for p=2677, n=2647, lam=185.

Run: python3 glv_hnp_phase2_lambda_followup.py
"""

import math
import random
from fpylll import IntegerMatrix, LLL, BKZ

# ---------------------------------------------------------------------------
# Reuse EC helpers from lambda_threshold.py
# ---------------------------------------------------------------------------

def modinv(a, m): return pow(a, -1, m)

def ec_add(P, Q, p):
    if P is None: return Q
    if Q is None: return P
    x1, y1 = P; x2, y2 = Q
    if x1 == x2:
        if (y1 + y2) % p == 0: return None
        s = 3*x1*x1*modinv(2*y1, p) % p
    else:
        s = (y2-y1)*modinv(x2-x1, p) % p
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
    if p%4 == 3: return pow(n,(p+1)//4,p)
    q,s = p-1,0
    while q%2==0: q//=2; s+=1
    z=2
    while pow(z,(p-1)//2,p) != p-1: z+=1
    m,c,t,r = s,pow(z,q,p),pow(n,q,p),pow(n,(q+1)//2,p)
    while True:
        if t==0: return 0
        if t==1: return r
        i,tmp = 0,t
        while tmp!=1: tmp=tmp*tmp%p; i+=1
        b=pow(c,1<<(m-i-1),p)
        m,c,t,r = i,b*b%p,t*b*b%p,r*b%p

def is_prime(n):
    if n<2: return False
    if n<4: return True
    if n%2==0 or n%3==0: return False
    i=5
    while i*i<=n:
        if n%i==0 or n%(i+2)==0: return False
        i+=6
    return True

def find_generator(p, b, n, seed=12345):
    rng = random.Random(seed)
    for _ in range(50000):
        x = rng.randint(0, p-1)
        rhs = (pow(x,3,p)+b) % p
        y = tonelli_shanks(rhs, p)
        if y is not None and y != 0:
            P = (x, y)
            if ec_mul(P, n, p) is None:
                return P
    return None

def glv_eigenvalue_both(n):
    neg3 = (n-3) % n
    sq = tonelli_shanks(neg3, n)
    if sq is None: return None, None
    inv2 = modinv(2, n)
    r1 = (n-1+sq)*inv2 % n
    r2 = (n-1+(n-sq))*inv2 % n
    if (r1*r1+r1+1)%n != 0: return None, None
    return r1, r2

def eisenstein_decompose(p):
    for a in range(1, 2*math.isqrt(p//3)+3):
        disc = 4*p - 3*a*a
        if disc < 0: break
        s = math.isqrt(disc)
        if s*s != disc: continue
        for num in [a+s, a-s]:
            if num%2 == 0:
                b = num//2
                if b>=0 and a*a-a*b+b*b==p:
                    return (a, b)
    return None

def j0_traces(a, b):
    return [2*a-b, -2*a+b, -(a+b), a+b, 2*b-a, a-2*b]

def gen_signatures(G, d_secret, m, n, lam, p, b, k1_bound, k2_bound, seed=42):
    rng = random.Random(seed)
    sigs = []
    attempts = 0
    while len(sigs) < m and attempts < 500000:
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
        s = modinv(k_full, n) * (h + d_secret*r) % n
        if s == 0: continue
        s_inv = modinv(s, n)
        A = h*s_inv % n
        B = r*s_inv % n
        assert (A + B*d_secret) % n == k_full
        sigs.append({'A': A, 'B': B, 'k1': k1, 'k2': k2, 'k_full': k_full})
    return sigs

def build_glv_lattice(sigs, n, lam, k1_bound, k2_bound):
    m = len(sigs)
    dim = 2*m+2
    S_K1 = max(1, n//k1_bound)
    S_K2 = max(1, n//k2_bound)
    M = [[0]*dim for _ in range(dim)]
    for i in range(m): M[i][i] = n*S_K1
    for i in range(m): M[m][i] = sigs[i]['B']*S_K1
    M[m][m] = 1
    for i in range(m):
        M[m+1+i][i] = -lam*S_K1
        M[m+1+i][m+1+i] = S_K2
    for i in range(m): M[2*m+1][i] = sigs[i]['A']*S_K1
    M[2*m+1][dim-1] = n
    return M

def run_lll_attack(p, b, n, lam, G, m, d_secret, k1_bound, seed=42):
    k2_bound = math.isqrt(n)+1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, b, k1_bound, k2_bound, seed)
    if len(sigs) < m: return False
    M = build_glv_lattice(sigs, n, lam, k1_bound, k2_bound)
    dim = 2*m+2
    A = IntegerMatrix.from_matrix(M)
    LLL.reduction(A)
    for row in A:
        last = row[dim-1]
        if abs(last) != n: continue
        sign = 1 if last > 0 else -1
        d_cand = (sign * row[m]) % n
        if d_cand == d_secret:
            return True
    return False

def sweep(label, p, b, n, lam, G, k1_bound, m_range, seeds, stop_at_3=True):
    k2_bound = math.isqrt(n)+1
    eff = k1_bound*k2_bound/n
    print(f"\n{label}")
    print(f"  p={p}, n={n} ({n.bit_length()}b), lam={lam}, lam/n={lam/n:.4f}")
    print(f"  K1={k1_bound}, K2={k2_bound}, eff={eff:.4f}")
    d_secret = random.Random(77777).randint(1, n-1)
    first3 = None
    for m in m_range:
        wins = 0
        for seed in seeds:
            ok = run_lll_attack(p, b, n, lam, G, m, d_secret, k1_bound, seed)
            wins += ok
        marker = " ← 3/3 ✓" if wins == len(seeds) else ""
        print(f"  m={m}: {wins}/{len(seeds)}{marker}")
        if wins == len(seeds) and first3 is None:
            first3 = m
            if stop_at_3:
                break
    if first3 is None:
        print(f"  → never 3/3 in {list(m_range)}")
    return first3

# ---------------------------------------------------------------------------
# Setup: the known failure curve
# ---------------------------------------------------------------------------

SEEDS = [42, 1234, 9999]
M_RANGE = list(range(5, 13))

p_fail = 2677
b_fail = 2
n_fail = 2647
lam_fail_min = 185          # smaller root, lam/n = 0.070
lam_fail_max = 2647-1-185   # = 2461, lam/n = 0.930

G_fail = find_generator(p_fail, b_fail, n_fail)
assert G_fail is not None, "Could not find generator for failure curve"
print("=" * 65)
print("GLV-HNP Phase 2: Lambda-threshold follow-up")
print("=" * 65)
print(f"\nFailure curve: p={p_fail}, n={n_fail}, b={b_fail}")
print(f"  lam_min={lam_fail_min} (lam/n={lam_fail_min/n_fail:.4f})")
print(f"  lam_max={lam_fail_max} (lam/n={lam_fail_max/n_fail:.4f})")

# ---------------------------------------------------------------------------
# (A) Re-test failure curve with K1=4 (matching success-case eff)
# ---------------------------------------------------------------------------
print("\n" + "=" * 65)
print("Part A: Re-test failure curve with K1=4 (reduced K1)")
print("=" * 65)

sweep("  A1: p=2677, n=2647, lam=185 (small), K1=4",
      p_fail, b_fail, n_fail, lam_fail_min, G_fail,
      k1_bound=4, m_range=M_RANGE, seeds=SEEDS)

sweep("  A2: p=2677, n=2647, lam=185 (small), K1=8 (original)",
      p_fail, b_fail, n_fail, lam_fail_min, G_fail,
      k1_bound=8, m_range=M_RANGE, seeds=SEEDS)

# ---------------------------------------------------------------------------
# (C) Root choice: use lam_max on the failure curve
# ---------------------------------------------------------------------------
print("\n" + "=" * 65)
print("Part C: Root choice — use lam_max=2461 instead of lam_min=185")
print("=" * 65)

sweep("  C1: p=2677, n=2647, lam=2461 (large root), K1=4",
      p_fail, b_fail, n_fail, lam_fail_max, G_fail,
      k1_bound=4, m_range=M_RANGE, seeds=SEEDS)

sweep("  C2: p=2677, n=2647, lam=2461 (large root), K1=8",
      p_fail, b_fail, n_fail, lam_fail_max, G_fail,
      k1_bound=8, m_range=M_RANGE, seeds=SEEDS)

# ---------------------------------------------------------------------------
# (D) K1 sweep on failure curve (lam_min=185)
# ---------------------------------------------------------------------------
print("\n" + "=" * 65)
print("Part D: K1 sweep — failure curve, lam_min=185")
print("=" * 65)

print(f"\n  {'K1':<6} {'eff':<8} {'1st m (3/3)'}")
for k1 in [2, 4, 6, 8, 10, 12]:
    k2_bound = math.isqrt(n_fail)+1
    eff = k1*k2_bound/n_fail
    d_secret = random.Random(77777).randint(1, n_fail-1)
    first3 = None
    for m in M_RANGE:
        wins = sum(
            run_lll_attack(p_fail, b_fail, n_fail, lam_fail_min, G_fail,
                           m, d_secret, k1, seed)
            for seed in SEEDS
        )
        if wins == len(SEEDS):
            first3 = m
            break
    status = str(first3) if first3 else f"never (>{M_RANGE[-1]})"
    print(f"  K1={k1:<4} eff={eff:.4f}  → {status}")

# ---------------------------------------------------------------------------
# (B) Very small lam/n: extend search to [0.01, 0.04)
# ---------------------------------------------------------------------------
print("\n" + "=" * 65)
print("Part B: Very small λ/n — search [0.01, 0.04)")
print("=" * 65)

found_very_small = []
p = 2049
while p < 2**14 and len(found_very_small) < 3:
    if p % 2 == 0: p += 1; continue
    if not is_prime(p) or p % 3 != 1: p += 2; continue
    eis = eisenstein_decompose(p)
    if eis is None: p += 2; continue
    a_e, b_e = eis
    for t in j0_traces(a_e, b_e):
        n_cand = p + 1 - t
        if n_cand < 1024 or not is_prime(n_cand) or n_cand % 3 != 1:
            continue
        r1, r2 = glv_eigenvalue_both(n_cand)
        if r1 is None: continue
        lam_min = min(r1, r2)
        ratio = lam_min / n_cand
        if 0.01 <= ratio < 0.04:
            for b_try in range(1, 50):
                rng_tmp = random.Random(b_try*999+7)
                for _ in range(20):
                    x = rng_tmp.randint(0, p-1)
                    rhs = (pow(x,3,p)+b_try) % p
                    y = tonelli_shanks(rhs, p)
                    if y is not None and y != 0:
                        Q = ec_mul((x,y), n_cand, p)
                        if Q is None:
                            G_found = find_generator(p, b_try, n_cand)
                            if G_found:
                                found_very_small.append((p, b_try, n_cand, lam_min, G_found))
                            break
            if found_very_small and found_very_small[-1][2] == n_cand:
                break
    p += 2

if not found_very_small:
    print("  No curves found in [0.01, 0.04)")
else:
    for (p0, b0, n0, lam0, G0) in found_very_small:
        sweep(f"  B: p={p0}, n={n0}, lam={lam0}, K1=4",
              p0, b0, n0, lam0, G0, k1_bound=4,
              m_range=M_RANGE, seeds=SEEDS)

print("\nDone.")
