"""
GLV-HNP Phase 2: Focused λ/n threshold sweep on 12-14 bit curves.

Selects representative j=0 prime-order curves at each 0.05-wide λ/n bin,
verifies the b value via brute-force point count, then runs 10 trials of
the column-scaled GLV-HNP LLL attack.

Goal: determine at what λ/n the attack starts failing for ~12-bit n.
"""

import math, random
from fpylll import IntegerMatrix, LLL

# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------
def isprime(n):
    if n < 2: return False
    if n == 2: return True
    if n % 2 == 0: return False
    for d in range(3, int(n**0.5)+1, 2):
        if n % d == 0: return False
    return True

def glv_roots(n):
    if n % 3 != 1: return []
    disc = (-3) % n
    for sq in range(n):
        if sq*sq % n == disc:
            inv2 = pow(2,-1,n)
            r1 = (-1+sq)*inv2 % n
            r2 = (-1-sq)*inv2 % n
            return sorted([r for r in [r1,r2] if (r*r+r+1)%n==0])
    return []

def count_points(p, b):
    count = 1
    for x in range(p):
        rhs = (x**3+b)%p
        if rhs == 0: count += 1
        elif pow(rhs,(p-1)//2,p)==1: count += 2
    return count

def find_b_for_order(p, n, max_b=100):
    """Find b ∈ [1, max_b] such that #E(F_p, b) = n."""
    for b in range(1, max_b+1):
        if count_points(p, b) == n:
            return b
    return None

def ec_add(P, Q, p):
    if P is None: return Q
    if Q is None: return P
    x1,y1=P; x2,y2=Q
    if x1==x2:
        if (y1+y2)%p==0: return None
        s=3*x1*x1*pow(2*y1,-1,p)%p
    else:
        s=(y2-y1)*pow(x2-x1,-1,p)%p
    x3=(s*s-x1-x2)%p; y3=(s*(x1-x3)-y1)%p
    return (x3,y3)

def ec_mul(P, k, p):
    R,Q=None,P
    while k>0:
        if k&1: R=ec_add(R,Q,p)
        Q=ec_add(Q,Q,p); k>>=1
    return R

def find_generator(p, n, b):
    for x in range(p):
        rhs=(x**3+b)%p
        if pow(rhs,(p-1)//2,p)!=1 and rhs!=0: continue
        for y in range(p):
            if y*y%p==rhs:
                P=(x,y)
                if ec_mul(P,n,p) is None: return P
    return None

# ---------------------------------------------------------------------------
# GLV-HNP attack
# ---------------------------------------------------------------------------
def run_attack(p, n, lam, b, n_trials=10, k1_bound=2, verbose=False):
    k2_bound = int(math.isqrt(n))+1
    G = find_generator(p, n, b)
    if G is None: return None, None
    if (lam*lam+lam+1)%n!=0: return None, None

    eff = k1_bound*k2_bound/n
    if eff>=1: return None, None
    m_thresh = math.ceil(math.log(n)/math.log(1.0/eff))
    m_sigs = max(4, m_thresh+2)

    s_k1=n//k1_bound; s_d=1; s_k2=n//k2_bound; s_kannan=n

    wins=0
    for trial in range(n_trials):
        rng=random.Random(trial*31337+int(lam)+int(p)*13)
        d_secret=rng.randint(1,n-1)
        sigs=[]; attempts=0
        while len(sigs)<m_sigs and attempts<300000:
            attempts+=1
            k1=rng.randint(0,k1_bound-1)
            k2=rng.randint(0,k2_bound-1)
            k_full=(k1+lam*k2)%n
            if k_full==0: continue
            R=ec_mul(G,k_full,p)
            if R is None: continue
            rv=R[0]%n
            if rv==0: continue
            h=rng.randint(0,n-1)
            s=pow(k_full,-1,n)*(h+d_secret*rv)%n
            if s==0: continue
            si=pow(s,-1,n)
            A=h*si%n; B=rv*si%n
            sigs.append({'A':A,'B':B,'k1':k1,'k2':k2})
        if len(sigs)<m_sigs: return None, None

        dim=2*m_sigs+2
        M=[[0]*dim for _ in range(dim)]
        for i in range(m_sigs): M[i][i]=n*s_k1
        M[m_sigs][m_sigs]=s_d
        for i in range(m_sigs): M[m_sigs][i]=sigs[i]['B']*s_k1
        for i in range(m_sigs):
            M[m_sigs+1+i][i]=-lam*s_k1
            M[m_sigs+1+i][m_sigs+1+i]=s_k2
        for i in range(m_sigs): M[2*m_sigs+1][i]=sigs[i]['A']*s_k1
        M[2*m_sigs+1][dim-1]=s_kannan

        A_mat=IntegerMatrix.from_matrix(M)
        LLL.reduction(A_mat)
        reduced=[[A_mat[i][j] for j in range(dim)] for i in range(dim)]

        d_found=False
        for row in reduced:
            last=row[dim-1]
            if abs(last)!=s_kannan: continue
            sign=1 if last>0 else -1
            d_cand=(sign*row[m_sigs])%n
            if d_cand==0: continue
            if d_cand==d_secret: d_found=True; break
        if d_found: wins+=1

    return wins/n_trials, m_sigs

# ---------------------------------------------------------------------------
# Candidate curves from CM formula (verified range)
# ---------------------------------------------------------------------------
def find_curves_cm(n_min=2000, n_max=16000):
    curves=[]
    seen_n=set()
    t_max=int(2*math.sqrt(n_max))+1
    u_max=int(2*math.sqrt(n_max/3))+1
    for t in range(-t_max, t_max+1):
        for u in range(1, u_max+1):
            val=t*t+3*u*u
            if val%4!=0: continue
            p=val//4
            if p<7 or not isprime(p) or p%3!=1: continue
            for sign in [1,-1]:
                n=p+1-sign*t
                if n<n_min or n>n_max or not isprime(n) or n%3!=1: continue
                if n in seen_n: continue
                roots=glv_roots(n)
                if len(roots)<2: continue
                lam1,lam2=roots[0],roots[1]
                seen_n.add(n)
                curves.append((p,n,lam1,lam2,lam1/n,lam2/n))
    return sorted(curves, key=lambda x: x[4])

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
if __name__=="__main__":
    print("="*72)
    print("GLV-HNP Phase 2: Focused λ/n threshold sweep (12-14 bit curves)")
    print("="*72)
    print()

    print("Building CM curve database...")
    curves=find_curves_cm(n_min=2000, n_max=12000)
    print(f"Found {len(curves)} candidates.")

    # Select one curve per 0.05-wide bin of r1 (smaller eigenvalue ratio)
    bucket_width=0.05
    bins={}
    for entry in curves:
        p,n,l1,l2,r1,r2=entry
        bk=round(round(r1/bucket_width)*bucket_width, 2)
        if bk not in bins: bins[bk]=[]
        bins[bk].append(entry)

    # Pick curve with n closest to 4000 in each bin
    target_n=4000
    selected={}
    for bk,entries in sorted(bins.items()):
        if bk>0.5: continue  # skip r1>0.5 (those are λ₂ values)
        entries.sort(key=lambda x: abs(x[1]-target_n))
        selected[bk]=entries[0]

    print(f"\n{'Bin':>6}  {'p':>7}  {'n':>6}  {'lam1':>6}  {'r1':>7}  {'r2':>7}  "
          f"{'b':>4}  {'m':>3}  {'lam1 rate':>10}  {'lam2 rate':>10}  rescue?")
    print("-"*80)

    results=[]
    for bk in sorted(selected.keys()):
        p,n,l1,l2,r1,r2=selected[bk]

        # Find b (brute-force for p < 16000)
        b=None
        if p < 16000:
            b=find_b_for_order(p,n,max_b=200)
        if b is None:
            print(f"{bk:>6.2f}  {p:>7}  {n:>6}  {l1:>6}  {r1:>7.4f}  {r2:>7.4f}  "
                  f"{'?':>4}  {'?':>3}  {'NO_B':>10}  {'NO_B':>10}  ?")
            continue

        # Attack with lam1 (small ratio)
        sr1,m1=run_attack(p,n,l1,b,n_trials=10)
        # Attack with lam2 (large ratio, complement)
        sr2,m2=run_attack(p,n,l2,b,n_trials=10)

        m_used=m1 or m2
        s1=f"{sr1:.1f}" if sr1 is not None else "ERR"
        s2=f"{sr2:.1f}" if sr2 is not None else "ERR"

        rescue="N/A"
        if sr1 is not None and sr2 is not None:
            if sr1<0.5 and sr2>=0.5: rescue="YES"
            elif sr1<0.5 and sr2<0.5: rescue="NO"
            elif sr1>=0.5: rescue="N/A"

        print(f"{bk:>6.2f}  {p:>7}  {n:>6}  {l1:>6}  {r1:>7.4f}  {r2:>7.4f}  "
              f"{b:>4}  {m_used!s:>3}  {s1:>10}  {s2:>10}  {rescue}")

        results.append({'bin':bk,'p':p,'n':n,'b':b,'l1':l1,'l2':l2,
                        'r1':r1,'r2':r2,'sr1':sr1,'sr2':sr2,'m':m_used,'rescue':rescue})

    print()
    print("="*72)
    print("SUMMARY")
    print("="*72)
    valid=[r for r in results if r['sr1'] is not None]
    fails=[r for r in valid if r['sr1']<0.5]
    wins=[r for r in valid if r['sr1']>=0.5]
    print(f"\nλ₁ attack on {len(valid)} curves:")
    win_ratios = sorted([f"{r['r1']:.3f}" for r in wins])
    fail_ratios = sorted([f"{r['r1']:.3f}" for r in fails])
    print(f"  Passing (>=50%): {len(wins)} -- ratios: {win_ratios}")
    print(f"  Failing (<50%):  {len(fails)} -- ratios: {fail_ratios}")

    rescued=[r for r in fails if r.get('rescue')=='YES']
    not_rescued=[r for r in fails if r.get('rescue')=='NO']
    print(f"\nComplementary λ₂ rescue:")
    print(f"  Rescued (λ₁ fails, λ₂ succeeds): {len(rescued)}")
    print(f"  Not rescued:                       {len(not_rescued)}")
    if rescued:
        for r in rescued:
            print(f"    r1={r['r1']:.3f}→r2={r['r2']:.3f}: lam2_rate={r['sr2']:.1f}")
    if not_rescued:
        for r in not_rescued:
            print(f"    r1={r['r1']:.3f}→r2={r['r2']:.3f}: lam2_rate={r['sr2']:.1f}")

    if fails and wins:
        max_fail=max(r['r1'] for r in fails)
        min_win=min(r['r1'] for r in wins)
        print(f"\nThreshold estimate: λ/n in [{max_fail:.3f}, {min_win:.3f}]")
