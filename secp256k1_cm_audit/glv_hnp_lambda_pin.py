"""
GLV-HNP Phase 2: pin the lambda/n threshold in (0.096, 0.130).

From bisect.py:
  lam/n = 0.096 (n=2053): FAIL (>m=15)
  lam/n = 0.130 (n=2203): m*=13

Find a 12-bit curve with lam/n in (0.096, 0.113) to confirm direction.
Also try 5 seeds instead of 3 to reduce variance.

Run: python3 glv_hnp_lambda_pin.py
"""

import math
import random
import sympy
from fpylll import IntegerMatrix, LLL

def modinv(a, m): return pow(a, -1, m)
def ec_add(P, Q, p):
    if P is None: return Q
    if Q is None: return P
    x1, y1 = P; x2, y2 = Q
    if x1 == x2:
        if (y1 + y2) % p == 0: return None
        s = 3*x1*x1*modinv(2*y1,p)%p
    else:
        s=(y2-y1)*modinv(x2-x1,p)%p
    x3=(s*s-x1-x2)%p; y3=(s*(x1-x3)-y1)%p
    return(x3,y3)
def ec_mul(P,k,p):
    if k==0:return None
    R,Q=None,P
    while k>0:
        if k&1:R=ec_add(R,Q,p)
        Q=ec_add(Q,Q,p); k>>=1
    return R
def tonelli_shanks(n,p):
    n%=p
    if n==0:return 0
    if pow(n,(p-1)//2,p)!=1:return None
    if p%4==3:return pow(n,(p+1)//4,p)
    q,s=p-1,0
    while q%2==0:q//=2;s+=1
    z=2
    while pow(z,(p-1)//2,p)!=p-1:z+=1
    m,c,t,r=s,pow(z,q,p),pow(n,q,p),pow(n,(q+1)//2,p)
    while True:
        if t==0:return 0
        if t==1:return r
        i,tmp=0,t
        while tmp!=1:tmp=tmp*tmp%p;i+=1
        b=pow(c,1<<(m-i-1),p)
        m,c,t,r=i,b*b%p,t*b*b%p,r*b%p
def find_generator(p,b,n):
    rng=random.Random(12345)
    for _ in range(50000):
        x=rng.randint(0,p-1); rhs=(pow(x,3,p)+b)%p; y=tonelli_shanks(rhs,p)
        if y is not None and y!=0:
            P=(x,y)
            if ec_mul(P,n,p) is None:return P
    return None
def eisenstein_decompose(p):
    for a in range(1,2*math.isqrt(p//3)+3):
        disc=4*p-3*a*a
        if disc<0:break
        s=math.isqrt(disc)
        if s*s!=disc:continue
        for num in[a+s,a-s]:
            if num%2==0:
                b=num//2
                if b>=0 and a*a-a*b+b*b==p:return(a,b)
    return None
def j0_traces(a,b):return[2*a-b,-2*a+b,-(a+b),a+b,2*b-a,a-2*b]
def glv_eigenvalue(n):
    neg3=(n-3)%n; sq=tonelli_shanks(neg3,n)
    if sq is None:return None,None
    inv2=modinv(2,n)
    r1=(n-1+sq)*inv2%n; r2=(n-1+(n-sq))*inv2%n
    if(r1*r1+r1+1)%n!=0:r1,r2=r2,r1
    if(r1*r1+r1+1)%n!=0:return None,None
    lam=min(r1,r2); return lam,n-1-lam

def run_lll(p,b,n,lam,G,K1,m,d_secret,seed):
    K2=math.isqrt(n)+1
    rng=random.Random(seed); sigs=[]; att=0
    while len(sigs)<m and att<500000:
        att+=1
        k1=rng.randint(0,K1-1); k2=rng.randint(0,K2-1)
        k_full=(k1+lam*k2)%n
        if k_full==0:continue
        R=ec_mul(G,k_full,p)
        if R is None:continue
        r=R[0]%n
        if r==0:continue
        h=rng.randint(0,n-1)
        s=modinv(k_full,n)*(h+d_secret*r)%n
        if s==0:continue
        s_inv=modinv(s,n)
        A=h*s_inv%n; B=r*s_inv%n
        sigs.append({'A':A,'B':B})
    if len(sigs)<m:return False
    S_K1=n//K1; S_K2=max(1,n//K2); S_KANNAN=n
    dim=2*m+2
    M=[[0]*dim for _ in range(dim)]
    for i in range(m):M[i][i]=n*S_K1
    for i in range(m):M[m][i]=sigs[i]['B']*S_K1
    M[m][m]=1
    for i in range(m):M[m+1+i][i]=-lam*S_K1; M[m+1+i][m+1+i]=S_K2
    for i in range(m):M[2*m+1][i]=sigs[i]['A']*S_K1
    M[2*m+1][dim-1]=S_KANNAN
    A=IntegerMatrix.from_matrix(M); LLL.reduction(A)
    for row_i in range(dim):
        last=A[row_i][dim-1]
        if abs(last)!=S_KANNAN:continue
        sign=1 if last>0 else -1
        d_cand=(sign*A[row_i][m])%n
        if d_cand==d_secret:return True
    return False

def sweep(label,p,b,n,lam,G,K1,m_range,seeds):
    ratio=lam/n
    print(f"\n[{label}] lam/n={ratio:.4f} n={n}({n.bit_length()}b) K1={K1}")
    first_full=None
    for m in m_range:
        wins=sum(run_lll(p,b,n,lam,G,K1,m,
                         random.Random(seed+5555).randint(1,n-1),seed)
                 for seed in seeds)
        marker=" *** 5/5" if wins==len(seeds) and first_full is None else ""
        print(f"  m={m:2d}: {wins}/{len(seeds)}{marker}")
        if wins==len(seeds) and first_full is None:first_full=m
    return first_full,ratio

# ---- Find 12-bit curves in (0.096, 0.130) ----------------------------------

TARGET = [(0.096, 0.110, "12b_0.096-0.110"), (0.110, 0.130, "12b_0.110-0.130")]

def find_narrow_window():
    found = {label: None for (_,_,label) in TARGET}
    print("Searching for 12-bit curves in (0.096, 0.130)...")
    p = sympy.nextprime(2047)
    while p < 4097 and any(v is None for v in found.values()):
        p = sympy.nextprime(p)
        if p % 3 != 1: continue
        eis = eisenstein_decompose(p)
        if eis is None: continue
        a_e, b_e = eis
        for t in j0_traces(a_e, b_e):
            n_cand = p + 1 - t
            if n_cand < 2 or not sympy.isprime(n_cand): continue
            if n_cand.bit_length() != 12: continue
            if n_cand % 3 != 1: continue
            lam, _ = glv_eigenvalue(n_cand)
            if lam is None: continue
            ratio = lam / n_cand
            for (lo, hi, label) in TARGET:
                if lo <= ratio < hi and found[label] is None:
                    found_b = None
                    for b_try in range(1, 300):
                        rng_b = random.Random(b_try*17)
                        for _ in range(40):
                            x = rng_b.randint(0, p-1)
                            rhs_x = (pow(x,3,p)+b_try)%p; y_x = tonelli_shanks(rhs_x,p)
                            if y_x is not None and y_x != 0:
                                if ec_mul((x,y_x),n_cand,p) is None:
                                    found_b = b_try; break
                        if found_b is not None: break
                    if found_b is None: continue
                    G_cand = find_generator(p, found_b, n_cand)
                    if G_cand is None: continue
                    found[label] = (p, found_b, n_cand, lam, G_cand, ratio)
                    print(f"  {label}: p={p}, b={found_b}, n={n_cand}, lam={lam},"
                          f" ratio={ratio:.5f}")
    return found

# ---- Main ------------------------------------------------------------------
print("=" * 65)
print("GLV-HNP: pinning threshold in (0.096, 0.130) at 12-bit / K1=8")
print("=" * 65)
SEEDS = [42, 1234, 9999, 7777, 31337]
K1 = 8
M_RANGE = range(4, 18)

found = find_narrow_window()

summary = []

# Retest the two boundary curves
print("\n--- Boundary confirmations (5 seeds) ---")
G_l = find_generator(2053, 10, 2053)  # ratio 0.096 (from bisect run — b=10)
# from bisect: 12b_0.07-0.10: p=2137, b=10, n=2053, lam=197 — wait, let me recheck
# The bisect output said: 12b_0.07-0.10: p=2137, b=10, n=2053, lam=197, ratio=0.0960
# But p=2137 != n=2053? That's fine: p and n are different for EC groups.
G_fail2 = find_generator(2137, 10, 2053)
m_f2, r_f2 = sweep("lower_bound lam/n=0.096 (p=2137)", 2137, 10, 2053, 197,
                    G_fail2, K1, M_RANGE, SEEDS)
summary.append((r_f2, m_f2, "lam/n=0.096"))

G_ok = find_generator(2281, 7, 2203)  # from bisect: 12b_0.10-0.14: p=2281, b=7, n=2203, lam=285
m_ok, r_ok = sweep("upper_bound lam/n=0.130 (p=2281)", 2281, 7, 2203, 285,
                    G_ok, K1, M_RANGE, SEEDS)
summary.append((r_ok, m_ok, "lam/n=0.130"))

print("\n--- Narrow window curves ---")
for (lo, hi, label) in TARGET:
    entry = found.get(label)
    if entry is None:
        print(f"\n[SKIP] {label}: no curve found in 12-bit range")
        summary.append(((lo+hi)/2, None, label + "_notfound"))
        continue
    p, b, n, lam, G, ratio = entry
    m_star, ratio = sweep(label, p, b, n, lam, G, K1, M_RANGE, SEEDS)
    summary.append((ratio, m_star, label))

print("\n" + "=" * 65)
print("SUMMARY — pinning run (5 seeds, m up to 17)")
print("=" * 65)
print(f"{'lam/n':>10}  {'m*':>9}  label")
print("-" * 45)
for ratio, m_star, label in sorted(summary, key=lambda x: x[0]):
    m_str = str(m_star) if m_star is not None else "FAIL(>17)"
    print(f"{ratio:10.5f}  {m_str:>9}  {label}")

sorted_s = sorted(summary, key=lambda x: x[0])
thr_lo = None; thr_hi = None
for ratio, m_star, label in sorted_s:
    if m_star is None:
        thr_lo = ratio
    elif thr_lo is not None and thr_hi is None:
        thr_hi = ratio; break

print()
if thr_lo is not None and thr_hi is not None:
    print(f"** Threshold in ({thr_lo:.5f}, {thr_hi:.5f})")
elif thr_lo is not None:
    print(f"** All curves above {thr_lo:.5f} succeed")
else:
    print(f"** All curves tested succeed")
print("\nDone.")
