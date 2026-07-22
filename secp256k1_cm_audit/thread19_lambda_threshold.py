"""
Thread 19: GLV-HNP Phase 2 — λ/n threshold scan.

Uses pre-validated curves from PARI scan (find_glv_curves.gp) with
lam verified: lam^2+lam+1≡0 mod n.

Tests LLL recovery at m=9, 5 seeds per curve.

PARI raw output (p,n,lam_raw,b,ratio) — note PARI normalizes lam<n/2
which may give the wrong eigenvalue. We re-verify and fix here.
"""

import random
from fpylll import IntegerMatrix, LLL

K1B = 2  # k1 in {0, 1}

def modinv(a, m):
    return pow(a, m-2, m)

def tonelli_shanks(n, p):
    if n % p == 0: return 0
    if pow(n, (p-1)//2, p) != 1: return None
    if p % 4 == 3: return pow(n, (p+1)//4, p)
    q, s = p-1, 0
    while q % 2 == 0: q //= 2; s += 1
    z = 2
    while pow(z, (p-1)//2, p) != p-1: z += 1
    M, c, t, r = s, pow(z, q, p), pow(n, q, p), pow(n, (q+1)//2, p)
    while True:
        if t == 1: return r
        i, tt = 1, t*t % p
        while tt != 1: i += 1; tt = tt*tt % p
        b = pow(c, 1 << (M-i-1), p)
        M, c, t, r = i, b*b%p, t*b*b%p, r*b%p

def get_lam(n):
    """Return lam in (0, n/2] with lam^2+lam+1≡0 mod n, or None."""
    sq = tonelli_shanks(-3 % n, n)
    if sq is None: return None
    for r in [sq, n-sq]:
        for lam in [(-1+r)*modinv(2,n)%n, (-1-r)*modinv(2,n)%n]:
            if (lam*lam + lam + 1) % n == 0:
                return lam if lam <= n//2 else n-lam if (n-lam)*(n-lam)%(n)+(n-lam)%n+1 == 0 else lam
    return None

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

def ec_mul(P, k, p):
    R,Q=None,P
    while k>0:
        if k&1: R=ec_add(R,Q,p)
        Q=ec_add(Q,Q,p); k>>=1
    return R

def find_gen(n, p, b):
    for x in range(p):
        rhs=(x**3+b)%p
        for y in range(p):
            if y*y%p==rhs:
                pt=(x,y)
                if ec_mul(pt,n,p) is None: return pt
    return None

def build_lattice(sigs, lam, n, k2b):
    S_K1=n//K1B; S_K2=n//k2b; S_KAN=n
    m=len(sigs); dim=2*m+2
    M=[[0]*dim for _ in range(dim)]
    for i in range(m): M[i][i]=n*S_K1
    for i in range(m): M[m][i]=sigs[i]['B']*S_K1
    M[m][m]=1
    for i in range(m):
        M[m+1+i][i]=-lam*S_K1; M[m+1+i][m+1+i]=S_K2
    for i in range(m): M[2*m+1][i]=sigs[i]['A']*S_K1
    M[2*m+1][dim-1]=S_KAN
    return M

def recover_d(A_mat, m, n):
    dim=2*m+2
    for i in range(dim):
        last=A_mat[i][dim-1]
        if abs(last)==n:
            sign=1 if last>0 else -1
            d_cand=(sign*A_mat[i][m])%n
            if d_cand>0: return d_cand
    return None

def run_attack(p, n, lam, b, k2b, d_secret, m, seed):
    G=find_gen(n,p,b)
    if G is None: return False
    rng=random.Random(seed); sigs=[]
    for _ in range(500000):
        k1=rng.randint(0,K1B-1); k2=rng.randint(0,k2b-1)
        k_full=(k1+lam*k2)%n
        if k_full==0: continue
        R=ec_mul(G,k_full,p)
        if R is None: continue
        r=R[0]%n
        if r==0: continue
        h=rng.randint(0,n-1)
        s=modinv(k_full,n)*(h+d_secret*r)%n
        if s==0: continue
        si=modinv(s,n)
        sigs.append({'A':h*si%n,'B':r*si%n,'k1':k1,'k2':k2})
        if len(sigs)==m: break
    if len(sigs)<m: return False
    M=build_lattice(sigs,lam,n,k2b)
    A_mat=IntegerMatrix.from_matrix(M)
    LLL.reduction(A_mat)
    return recover_d(A_mat,m,n)==d_secret

def test_curve(p, n, lam, b, m=9, seeds=5, label=""):
    k2b=int(n**0.5)+2
    wins=0; rng=random.Random(42424242)
    for seed in range(seeds):
        d=rng.randint(1,n-1)
        wins+=run_attack(p,n,lam,b,k2b,d,m,seed)
    ratio=lam/n; tag="PASS" if wins==seeds else ("PARTIAL" if wins>0 else "FAIL")
    print(f"  lam/n={ratio:.3f}  p={p},n={n}  {wins}/{seeds}  [{tag}]  {label}")
    return ratio, wins/seeds

# --------------------------------------------------------------------------
# Pre-validated curves: lam^2+lam+1≡0 mod n confirmed below.
# PARI raw: (p,n,lam_pari,b) — we re-derive lam from scratch here.
# --------------------------------------------------------------------------
PARI_CURVES = [
    # (p, n, b)  — lam derived by get_lam(n)
    (787, 757, 2),    # lam/n ≈ 0.036
    (547, 547, 2),    # lam/n ≈ 0.073
    (331, 331, 2),    # lam/n ≈ 0.094
    (577, 613, 7),    # lam/n ≈ 0.106
    (277, 283, 5),    # lam/n ≈ 0.156
    (229, 223, 6),    # lam/n ≈ 0.175
    (379, 409, 2),    # lam/n ≈ 0.130
    (829, 823, 2),    # lam/n ≈ 0.211
    (691, 643, 3),    # lam/n ≈ 0.277
    (349, 313, 2),    # lam/n ≈ 0.316
    (877, 937, 2),    # lam/n ≈ 0.345
    (313, 349, 10),   # lam/n ≈ 0.352
    (283, 277, 3),    # lam/n ≈ 0.422
    (211, 199, 2),    # lam/n ≈ 0.533
]

print("=" * 70)
print("Thread 19: GLV-HNP Phase 2 — λ/n threshold scan")
print("LLL attack at m=9, 5 seeds, k2b=ceil(sqrt(n))")
print("=" * 70)
print()

# Verify and collect valid curves
valid = []
for p, n, b in PARI_CURVES:
    lam = get_lam(n)
    if lam is None:
        print(f"  SKIP p={p},n={n}: no GLV eigenvalue mod n")
        continue
    # Double-check
    if (lam*lam+lam+1)%n != 0:
        # Try n-lam
        lam2 = n - lam
        if (lam2*lam2+lam2+1)%n == 0:
            lam = lam2
        else:
            print(f"  SKIP p={p},n={n}: eigenvalue check failed lam={lam}")
            continue
    valid.append((p, n, lam, b, lam/n))

valid.sort(key=lambda x: x[4])
print(f"Testing {len(valid)} curves:")
print()

results = []
for p, n, lam, b, ratio in valid:
    r, rate = test_curve(p, n, lam, b, m=9, seeds=5, label=f"(b={b})")
    results.append((r, rate))

print()
print("Summary:")
for r, rate in results:
    bar = "#" * int(rate*10) + "." * (10-int(rate*10))
    print(f"  lam/n={r:.3f}  [{bar}] {rate:.0%}")

failing = [r for r,rate in results if rate < 0.4]
passing = [r for r,rate in results if rate >= 0.6]
print()
if failing and passing:
    print(f"Transition: last clear FAIL at lam/n≤{max(failing):.3f}, "
          f"first clear PASS at lam/n≥{min(passing):.3f}")
print("secp256k1 lam/n ≈ 0.330")
