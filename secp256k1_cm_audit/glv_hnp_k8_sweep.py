"""
GLV-HNP k1_bound=8 sweep: determine if n=2647 failure is λ/n-specific or n-size-specific.

Uses the CM curve database (find_curves_cm) to select curves at different λ/n bins,
runs the attack with k1_bound=8 (matching glv_hnp_phase2_20bit.py settings),
and checks whether failure correlates with λ/n or with n-size.
"""

import math, random
from fpylll import IntegerMatrix, LLL

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

def find_b_for_order(p, n, max_b=200):
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

def run_attack(p, n, lam, b, k1_bound, n_trials=10):
    k2_bound = int(math.isqrt(n))+1
    G = find_generator(p, n, b)
    if G is None: return None, None, None
    if (lam*lam+lam+1)%n!=0: return None, None, None

    eff = k1_bound*k2_bound/n
    if eff>=1: return None, None, None
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
        if len(sigs)<m_sigs: return None, None, None

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

    return wins/n_trials, m_sigs, eff

def find_curves_cm(n_min=2000, n_max=5000):
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

if __name__=="__main__":
    print("="*76)
    print("GLV-HNP k1_bound=8 sweep: is failure λ/n-specific or n-size-specific?")
    print("="*76)
    print()

    # Include the known failing case
    known_fail = (2677, 2647, 185, 2461, 185/2647, 2461/2647)

    print("Building CM curve database (n in [2000,5000])...")
    curves = find_curves_cm(n_min=2000, n_max=5000)
    print(f"Found {len(curves)} candidates.")

    # Select curves per 0.05-wide λ/n bin, picking ones with n closest to 3000
    bucket_width=0.05
    bins={}
    for entry in curves:
        p,n,l1,l2,r1,r2=entry
        bk=round(round(r1/bucket_width)*bucket_width, 2)
        if bk>0.5: continue
        if bk not in bins: bins[bk]=[]
        bins[bk].append(entry)

    target_n=3000
    selected=[]
    for bk in sorted(bins.keys()):
        entries=bins[bk]
        entries.sort(key=lambda x: abs(x[1]-target_n))
        # Take top 3 per bin for robustness
        for e in entries[:3]:
            selected.append((bk, e))

    print()
    print(f"{'bin':>5}  {'p':>6}  {'n':>5}  {'lam1':>6}  {'r1':>7}  {'eff':>6}  "
          f"{'m':>3}  {'k8_rate':>8}  note")
    print("-"*72)

    # Also test the known failing case
    p,n,l1,l2,r1,r2 = known_fail
    b = find_b_for_order(p, n, max_b=200)
    if b is not None:
        sr,m,eff = run_attack(p, n, l1, b, k1_bound=8, n_trials=10)
        s = f"{sr:.1f}" if sr is not None else "ERR"
        e_str = f"{eff:.3f}" if eff is not None else "?"
        print(f"{'0.07':>5}  {p:>6}  {n:>5}  {l1:>6}  {r1:>7.4f}  {e_str:>6}  "
              f"{m!s:>3}  {s:>8}  KNOWN_FAIL")

    seen_pairs=set()
    results=[]
    for bk, entry in selected:
        p,n,l1,l2,r1,r2=entry
        if (p,n) in seen_pairs: continue
        seen_pairs.add((p,n))

        if p>=5000: continue
        b = find_b_for_order(p, n, max_b=200)
        if b is None: continue

        sr,m,eff = run_attack(p, n, l1, b, k1_bound=8, n_trials=10)
        if sr is None: continue

        s = f"{sr:.1f}"
        e_str = f"{eff:.3f}"
        print(f"{bk:>5.2f}  {p:>6}  {n:>5}  {l1:>6}  {r1:>7.4f}  {e_str:>6}  "
              f"{m!s:>3}  {s:>8}")
        results.append({'bk':bk,'p':p,'n':n,'l1':l1,'r1':r1,'eff':eff,'m':m,'sr':sr})

    print()
    print("="*76)
    print("ANALYSIS")
    print("="*76)

    if results:
        fails = [r for r in results if r['sr'] < 0.5]
        wins  = [r for r in results if r['sr'] >= 0.5]
        print(f"\nWith k1_bound=8:")
        print(f"  Pass (>=50%): {len(wins)}  bins={sorted(set(r['bk'] for r in wins))}")
        print(f"  Fail (<50%):  {len(fails)}  bins={sorted(set(r['bk'] for r in fails))}")
        if fails:
            print(f"\n  Failing curves:")
            for r in sorted(fails, key=lambda x: x['r1']):
                print(f"    n={r['n']}, r1={r['r1']:.3f}, eff={r['eff']:.3f}, m={r['m']}, rate={r['sr']:.1f}")
        if wins:
            print(f"\n  Passing curves:")
            for r in sorted(wins, key=lambda x: x['r1']):
                print(f"    n={r['n']}, r1={r['r1']:.3f}, eff={r['eff']:.3f}, m={r['m']}, rate={r['sr']:.1f}")
