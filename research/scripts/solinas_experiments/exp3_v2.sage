"""Exp3 v2: wider b search + check both cube roots of λ."""
import time, sys

ALL = [
    ("P-192", 64,  0xffffffffdfffffff),
    ("P-192", 96,  0xffffffffffffffdfffffffff),
    ("P-192", 112, 0xfffffffffffffffffff7ffffffff),
    ("P-192", 120, 0xffffffffffffffffffffff7fffffff),
    ("P-224", 56,  0xffffffff000001),
    ("P-224", 64,  0xffffffffff000001),
    ("P-224", 70,  0x3fffffff0000000001),
    ("P-224", 80,  0xffffffff000000000001),
    ("P-224", 84,  0xffff00000000000000001),
    ("P-224", 96,  0xffffffffffffffff00000001),
    ("P-224", 126, 0x3ffffffffffffffffffffffffffc0001),
    ("P-256", 128, 0xfffe0000800000000000ffffffffffff),
    ("P-384", 96,  0xfffffffffffffffbff8000ff),
    ("P-384", 112, 0xffffffffffffffffffdff80001ff),
    ("P-384", 120, 0xfffffffffffffffffffdffc00003ff),
    ("P-521", 61,  0x1fffffffffffffff),
    ("P-521", 89,  0x1ffffffffffffffffffffff),
    ("P-521", 107, 0x7ffffffffffffffffffffffffff),
    ("P-521", 127, 0x7fffffffffffffffffffffffffffffff),
]

print(f"  {'curve':14s} {'b':>5s}  {'h':>5s}  log2(n')  endo OK")
print("-" * 60)
sys.stdout.flush()

# Wider range, larger cofactor tolerance, also try b values that are perfect cubes/squares
for fam, n, p in ALL:
    p = ZZ(p); Fp = GF(p)
    # cube root of unity
    zeta = None
    for c in range(2, 100):
        cand = Fp(c)^((p-1)//3)
        if cand != 1: zeta = cand; break
    found = None
    bmax = 500 if n <= 96 else 300
    cof_max = 16 if n <= 96 else 24
    for b in range(1, bmax+1):
        try:
            E = EllipticCurve(Fp, [0, b]); N = E.order()
        except Exception: continue
        # try to peel off small primes
        rem = N; cof = 1
        for q in prime_range(2, 200):
            while rem % q == 0: rem //= q; cof *= q
            if cof > 10^5: break
        if rem > 1 and rem.is_prime() and cof <= cof_max:
            found = (b, N, cof, rem); break
    if not found:
        print(f"  {fam}/{n:<4d}    -- no GLV curve in b≤{bmax} h≤{cof_max}")
        sys.stdout.flush(); continue
    b, N, cof, nprime = found
    # Solve λ
    R = Zmod(nprime)
    if nprime % 3 != 1:
        print(f"  {fam}/{n:<4d}  b={b:3d} h={cof:3d} log2(n')={float(log(nprime,2)):.2f}  λ undefined (n' ≢ 1 mod 3)")
        sys.stdout.flush(); continue
    s = R(-3).sqrt()
    lam = R((-1 + ZZ(s))) / R(2)
    lam2 = R(-1) - lam   # the other root: -1 - λ
    # Verify on prime-order subgroup point
    while True:
        P = E.random_point()
        Pp = cof * P
        if Pp != E(0) and Pp.order() == nprime: break
    xP, yP = Pp.xy()
    phi_P = E(zeta * xP, yP)
    ok1 = (phi_P == ZZ(lam) * Pp)
    ok2 = (phi_P == ZZ(lam2) * Pp)
    # also check other zeta
    phi2 = E(zeta**2 * xP, yP)
    ok3 = (phi2 == ZZ(lam) * Pp)
    ok4 = (phi2 == ZZ(lam2) * Pp)
    endo = ok1 or ok2 or ok3 or ok4
    print(f"  {fam}/{n:<4d}  b={b:3d} h={cof:3d} log2(n')={float(log(nprime,2)):.2f}  endo OK = {endo}  ({['','λ','-1-λ',''][ok1+2*ok2]} via {['ζ','ζ²'][ok3 or ok4]})" if endo else f"  {fam}/{n:<4d}  b={b:3d} h={cof:3d} endo CHECK FAILED")
    sys.stdout.flush()
print("DONE")
