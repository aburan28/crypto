"""
GLV-HNP Phase 2 — IT threshold vs λ/n diagnosis.

The 2026-07-26 run found λ/n = 0.07 (p=2677, n=2647) always fails with ≤12 sigs.
Hypothesis A (original): failure due to small λ/n (structural).
Hypothesis B (new):      failure due to IT threshold not met (K1*K2/n ≈ 0.51 needs m≥17).

Today's sweep (K1=2, m=6) showed the known-fail curve PASSES 3/3, suggesting Hyp B.
This script does a definitive test:
  - Curve: p=2677, n=2647, λ=185, λ/n=0.07 (the "known-fail" curve)
  - K1_BOUND = n//100 = 26  (replicating original 1%-bias setup)
  - K2_BOUND = isqrt(n)+1 = 52
  - K1*K2/n = 26*52/2647 ≈ 0.51  →  IT threshold m* ≈ 17
  - Sweep m from 5 to 25 with 3 seeds; find minimum m for ≥2/3 success.

If Hyp B: failure for m<17, success for m≥17, λ/n irrelevant.
If Hyp A: failure for all m at this λ/n.
"""
import math, random
from fpylll import IntegerMatrix, LLL, BKZ

def modinv(a, m): return pow(a, -1, m)

def tonelli_shanks(n, p):
    n %= p
    if n == 0: return 0
    if pow(n, (p-1)//2, p) != 1: return None
    if p % 4 == 3: return pow(n, (p+1)//4, p)
    q, s = p-1, 0
    while q % 2 == 0: q //= 2; s += 1
    z = 2
    while pow(z, (p-1)//2, p) != p-1: z += 1
    mc, c, t, r = s, pow(z,q,p), pow(n,q,p), pow(n,(q+1)//2,p)
    while True:
        if t == 0: return 0
        if t == 1: return r
        i, tmp = 1, t*t % p
        while tmp != 1: tmp = tmp*tmp % p; i += 1
        if mc <= i: return None   # guard
        b = pow(c, 1 << (mc-i-1), p)
        mc, c, t, r = i, b*b%p, t*c*c%p, r*b%p

def ec_add(P, Q, p):
    if P is None: return Q
    if Q is None: return P
    x1,y1=P; x2,y2=Q
    if x1==x2:
        if (y1+y2)%p==0: return None
        s=3*x1*x1*modinv(2*y1,p)%p
    else:
        s=(y2-y1)*modinv(x2-x1,p)%p
    x3=(s*s-x1-x2)%p; y3=(s*(x1-x3)-y1)%p
    return (x3,y3)

def ec_mul(P,k,p):
    if k==0: return None
    R,Q=None,P
    while k>0:
        if k&1: R=ec_add(R,Q,p)
        Q=ec_add(Q,Q,p); k>>=1
    return R

def find_gen(p,b,n):
    for x in range(p):
        rhs=(x*x*x+b)%p
        sq=tonelli_shanks(rhs,p)
        if sq is None: continue
        for y in [sq,p-sq]:
            if ec_mul((x,y),n,p) is None: return (x,y)
    return None

def gen_sigs(G,d,n,lam,p,K1,K2,m,rng):
    sigs=[]; att=0
    while len(sigs)<m and att<500_000:
        att+=1
        k1=rng.randint(0,K1-1); k2=rng.randint(0,K2-1)
        kf=(k1+lam*k2)%n
        if kf==0: continue
        R=ec_mul(G,kf,p)
        if R is None: continue
        r=R[0]%n
        if r==0: continue
        h=rng.randint(0,n-1)
        s=modinv(kf,n)*(h+d*r)%n
        if s==0: continue
        si=modinv(s,n)
        A=h*si%n; B=r*si%n
        assert (A+B*d)%n==kf
        sigs.append({'A':A,'B':B,'k1':k1,'k2':k2})
    return sigs if len(sigs)==m else None

def build_lattice(sigs,n,lam,K1,K2):
    m=len(sigs); dim=2*m+2
    SK1=n//K1; SD=1; SK2=n//K2; SKAN=n
    M=[[0]*dim for _ in range(dim)]
    for i in range(m): M[i][i]=n*SK1
    for i in range(m): M[m][i]=sigs[i]['B']*SK1
    M[m][m]=SD
    for i in range(m):
        M[m+1+i][i]=-lam*SK1
        M[m+1+i][m+1+i]=SK2
    for i in range(m): M[2*m+1][i]=sigs[i]['A']*SK1
    M[2*m+1][dim-1]=SKAN
    return M,SK1,SD,SK2,SKAN

def run_one(G,d,n,lam,p,K1,K2,m,seed):
    rng=random.Random(seed)
    sigs=gen_sigs(G,d,n,lam,p,K1,K2,m,rng)
    if sigs is None: return False
    M,SK1,SD,SK2,SKAN=build_lattice(sigs,n,lam,K1,K2)
    dim=len(M)
    A=IntegerMatrix.from_matrix(M)
    LLL.reduction(A)
    for i in range(dim):
        row=[A[i][j] for j in range(dim)]
        last=row[dim-1]
        if abs(last)!=SKAN: continue
        sign=1 if last>0 else -1
        dc=(sign*row[m]*modinv(SD,n))%n
        if dc==d: return True
        # sig-verify fallback
        ok=True
        for s in sigs:
            kf=(s['A']+s['B']*dc)%n
            found=any((kf-lam*k2)%n<K1 for k2 in range(K2))
            if not found: ok=False; break
        if ok and dc!=0:
            if dc==d: return True
    return False

# ---------------------------------------------------------------------------
# Experiment
# ---------------------------------------------------------------------------
print("="*65)
print("GLV-HNP Phase 2 — IT threshold diagnosis for λ/n = 0.07 curve")
print("="*65)

# Known-fail curve
P_FIELD = 2677
B_COEFF = 7
N       = 2647
LAM     = 185

K1_BOUND = N // 100   # 1% bias
K2_BOUND = math.isqrt(N) + 1  # GLV bound

print(f"\nCurve: y^2 = x^3 + {B_COEFF} / F_{P_FIELD}")
print(f"n={N}, λ={LAM}, λ/n={LAM/N:.4f}")
print(f"K1_BOUND={K1_BOUND} (≈1%), K2_BOUND={K2_BOUND}")
k_frac = K1_BOUND * K2_BOUND / N
it_m = math.log(N) / math.log(1.0 / k_frac) if k_frac < 1 else float('inf')
print(f"K1*K2/n = {K1_BOUND*K2_BOUND}/{N} ≈ {k_frac:.4f}")
print(f"IT threshold: (K1*K2/n)^m < 1/n  =>  m ≥ {it_m:.1f}")
print()

G = find_gen(P_FIELD, B_COEFF, N)
assert G is not None, "Generator not found"

SEEDS = [42, 1234, 9999]
rng_d = random.Random(77)
D_SECRET = rng_d.randint(1, N-1)
print(f"Planted d = {D_SECRET}")
print()

print(f"{'m':>4}  {'K1*K2/n':>8}  {'wins':>5}  verdict")
print("-" * 45)

# Also test a "high λ/n" curve for comparison
# Find one with λ/n ≈ 0.44 in same n range
# Use p=2683 and search
def find_high_lam_curve(p_min, p_max):
    from math import isqrt
    def is_prime(n):
        if n<2: return False
        if n<4: return True
        if n%2==0: return False
        d,s=n-1,0
        while d%2==0: d//=2; s+=1
        for a in [2,3,5,7,11,13,17,19,23,29,31,37]:
            if a>=n: continue
            x=pow(a,d,n)
            if x==1 or x==n-1: continue
            for _ in range(s-1):
                x=x*x%n
                if x==n-1: break
            else: return False
        return True
    def curve_order(p,b):
        count=1
        for x in range(p):
            rhs=(x*x*x+b)%p
            if rhs==0: count+=1
            elif pow(rhs,(p-1)//2,p)==1: count+=2
        return count
    def glv_eig(n):
        for l in range(1,n):
            if (l*l+l+1)%n==0: return l
        return None
    for p in range(p_min, p_max+1):
        if not is_prime(p): continue
        if p%6!=1: continue
        for b in [1,2,3,5,7]:
            n=curve_order(p,b)
            if not is_prime(n): continue
            if n%3!=1: continue
            lam=glv_eig(n)
            if lam is None: continue
            lam_s=min(lam,n-1-lam)
            if lam_s/n > 0.35:  # high λ/n
                return p,b,n,lam_s
    return None,None,None,None

high = find_high_lam_curve(2600, 2900)

threshold_m = None
for m in range(5, 26):
    wins = sum(run_one(G, D_SECRET, N, LAM, P_FIELD, K1_BOUND, K2_BOUND, m, s) for s in SEEDS)
    verdict = "PASS" if wins >= 2 else "FAIL"
    print(f"{m:>4}  {k_frac:>8.4f}  {wins}/{len(SEEDS)}  {verdict}")
    if wins >= 2 and threshold_m is None:
        threshold_m = m

print()
if threshold_m:
    print(f"First PASS at m={threshold_m}  (IT threshold: {it_m:.1f})")
    print(f"Gap: m_first_pass - m_IT = {threshold_m - it_m:.1f}")
    print()
    if abs(threshold_m - it_m) <= 3:
        print("CONCLUSION: Failure was IT-threshold limited, NOT λ/n structural.")
        print("  The λ/n = 0.07 hypothesis from 2026-07-26 is REFUTED.")
        print("  With enough signatures (m ≥ IT threshold), LLL succeeds for any λ/n.")
    else:
        print("CONCLUSION: First pass is far from IT threshold. May be structural.")
else:
    print("FAIL at all m tested. Structural issue confirmed (or need even more sigs).")

# Compare to a high-λ/n curve in same n range
if high[0] is not None:
    hp, hb, hn, hlam = high
    print(f"\nComparison: high-λ/n curve p={hp}, b={hb}, n={hn}, λ={hlam}, λ/n={hlam/hn:.4f}")
    hK1 = hn // 100
    hK2 = math.isqrt(hn) + 1
    hfrac = hK1 * hK2 / hn
    hit_m = math.log(hn) / math.log(1.0 / hfrac)
    print(f"K1={hK1}, K2={hK2}, K1*K2/n≈{hfrac:.4f}, IT threshold m*≈{hit_m:.1f}")
    hG = find_gen(hp, hb, hn)
    hd = random.Random(88).randint(1, hn-1)
    print(f"{'m':>4}  {'K1*K2/n':>8}  {'wins':>5}  verdict")
    for m in range(5, 22):
        wins = sum(run_one(hG, hd, hn, hlam, hp, hK1, hK2, m, s) for s in SEEDS)
        verdict = "PASS" if wins >= 2 else "FAIL"
        print(f"{m:>4}  {hfrac:>8.4f}  {wins}/{len(SEEDS)}  {verdict}")
else:
    print("\n(No high-λ/n comparison curve found in range 2600-2900)")

print("\nDone.")
