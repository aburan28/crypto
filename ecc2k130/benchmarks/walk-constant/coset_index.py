"""Coset confinement of the sigma walk, ecc2k130/WALK-CONSTANT.md section 3.

The sigma walk multiplies the scalar by 1 + s^j (j = 3..10) at every step,
and classes identify s and -1, so a walk never leaves the coset of the
subgroup <1 + s^j, s, -1> of Z_l^* it starts in.  Walks started in
different cosets never meet; an index above 1 would be a penalty no
branch statistic shows.  This prints the index at every degree the
measurement ran, and at ECC2K-130.  Run: python3 coset_index.py
"""
import math, random
def is_prime(n):
    if n < 2: return False
    for p in [2,3,5,7,11,13,17,19,23,29,31,37,41,43]:
        if n % p == 0: return n == p
    d, r = n-1, 0
    while d % 2 == 0: d//=2; r+=1
    for a in [2,3,5,7,11,13,17,19,23,29,31,37,41]:
        x = pow(a, d, n)
        if x in (1, n-1): continue
        for _ in range(r-1):
            x = x*x % n
            if x == n-1: break
        else: return False
    return True
def rho(n):
    if n % 2 == 0: return 2
    while True:
        y, c, m = random.randrange(1,n), random.randrange(1,n), 128
        g, r, q = 1, 1, 1
        while g == 1:
            x = y
            for _ in range(r): y = (y*y + c) % n
            k = 0
            while k < r and g == 1:
                ys = y
                for _ in range(min(m, r-k)):
                    y = (y*y + c) % n; q = q*abs(x-y) % n
                g = math.gcd(q, n); k += m
            r *= 2
        if g == n:
            g = 1
            while g == 1:
                ys = (ys*ys + c) % n; g = math.gcd(abs(x-ys), n)
        if g != n: return g
def factor(n, out):
    for p in range(2, 100000):
        while n % p == 0: out[p] = out.get(p,0)+1; n//=p
        if p*p > n: break
    if n == 1: return
    stack=[n]
    while stack:
        m = stack.pop()
        if m == 1: continue
        if is_prime(m): out[m]=out.get(m,0)+1; continue
        d = rho(m); stack += [d, m//d]
def index(l, s, n):
    M = l-1; f = {}; factor(M, f)
    gens = [(1 + pow(s, j, l)) % l for j in range(3, 11)] + [s % l, l-1]
    assert pow(s, n, l) == 1, "s is not an n-th root of unity"
    k = 1
    for q, e in f.items():
        # largest q^t (t <= e) such that every generator is a q^t-th power
        t = 0
        while t < e and all(pow(g, M // q**(t+1), l) == 1 for g in gens): t += 1
        k *= q**t
    return f, k
def sqrt_mod(a, p):
    a %= p
    if p % 4 == 3: return pow(a, (p+1)//4, p)
    q, s = p-1, 0
    while q % 2 == 0: q//=2; s+=1
    z = 2
    while pow(z, (p-1)//2, p) != p-1: z += 1
    m, c, t, r = s, pow(z, q, p), pow(a, q, p), pow(a, (q+1)//2, p)
    while t != 1:
        i, tt = 0, t
        while tt != 1: tt = tt*tt % p; i += 1
        b = pow(c, 1 << (m-i-1), p); m, c, t, r = i, b*b % p, t*b*b % p, r*b % p
    return r

CURVES = [  # (n, l): the largest prime factor of #E(F_2^n), y^2 + xy = x^3 + 1
    (19, 130873), (23, 2095853), (31, 1439393), (37, 230603167),
    (41, 549756390943), (59, 10063074221),
    (131, 680564733841876926932320129493409985129),
]
for n, l in CURVES:
    r = sqrt_mod(-7, l)
    inv2 = pow(2, -1, l)
    for s in [(-1 + r) * inv2 % l, (-1 - r) * inv2 % l]:
        assert (s * s + s + 2) % l == 0
        if pow(s, n, l) != 1:  # the Frobenius eigenvalue has order dividing n
            continue
        f, k = index(l, s, n)
        print(f"n = {n}: l = {l}, s = {s}; index of <1 + s^j, s, -1> in Z_l^* = {k}")
