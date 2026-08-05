"""
Thread 23e — cmin is the Eisenstein decomposition of n (closed form), and what
it says about secp256k1.

Thread 23d established (AUC 0.953, matched 8-vs-8 at 12 bits) that
    r = cmin/sqrt(n),   cmin = min_{0<|e|<=K2} |lam*e mod n|_centered
predicts the GLV-HNP Phase-2 K1 wall: low-r curves walled at K1 ~ 12, typical-r
curves at K1 ~ 6.

The low-r curves found there look like Eisenstein decompositions:
    n=3907, cmin=1, e*=63:   1^2 - 1*63 + 63^2 = 3907
    n=2551, cmin=1, e*=51:   1^2 - 1*51 + 51^2 = 2551
    n=2503, cmin=2, e*=49:   2^2 + 2*49 + 49^2 = 2503
    n=2557, cmin=3, e*=49:   3^2 + 3*49 + 49^2 = 2557

Conjecture (C-eis).  For n prime = 1 mod 3 with lam^2 + lam + 1 = 0 mod n, the
lattice L_lam = {(x,y) : x = -lam*y mod n} has covolume n and is a scaled copy
of the Eisenstein lattice Z[omega].  Its minimal vectors are the six associates
of the Eisenstein integer of norm n, so with n = a^2 - a*b + b^2,
    cmin = min |x| over the six associates (x,y) of a + b*omega
           SUBJECT TO (i) the box constraint |y| <= K2 = floor(sqrt(n)) + 1
                  and (ii) x = +-lam*y (mod n).
Both side conditions are essential.  (i) because the associate with the
smallest |x| typically has |y| ~ sqrt(n) just above K2, hence unreachable;
(ii) because n splits as pi*pi_bar in Z[omega] and only one of the two
unit-orbits lies in L_lam.
That is computable in polylog time by Cornacchia -- no lattice reduction, no
signatures, no curve arithmetic.

This matters because it makes the K1 wall a published, precomputable property
of the group order.  Part 2 evaluates it on the real curves.

Run: python3 glv_hnp_phase2_cmin_eisenstein.py
"""

import math

import sympy

_src = open(__file__.replace('cmin_eisenstein', 'projected')).read()
_head = _src.split('# ' + '=' * 75)[0]
_G = {}
exec(compile(_head, 'glv_hnp_phase2_projected.py', 'exec'), _G)

eisenstein_decompose = _G['eisenstein_decompose']
j0_traces = _G['j0_traces']
glv_roots = _G['glv_roots']
lam_star = _G['lam_star']


def cmin_bruteforce(n, lam, k2_bound):
    best, best_e = n, 0
    for e in range(1, k2_bound + 1):
        r = lam * e % n
        r = min(r, n - r)
        if r < best:
            best, best_e = r, e
    return best, best_e


def associates(a, b):
    """All 12 isometric images of a + b*omega for N(x+y*omega) = x^2-xy+y^2:
    the 6 unit multiples (x,y) -> (-y, x-y) with negation, TIMES the conjugation
    (x,y) -> (y,x).  The conjugate orbit is NOT reachable from the unit orbit
    and omitting it was the bug in the first version of this script."""
    orb = []
    for (x0, y0) in ((a, b), (b, a)):
        x, y = x0, y0
        for _ in range(3):
            orb.append((x, y))
            orb.append((-x, -y))
            x, y = -y, x - y
    return orb


def cornacchia_d3(n):
    """n prime = 1 mod 3  ->  (a,b) with a^2 - a*b + b^2 = n.  Polylog time."""
    if n % 3 != 1:
        return None
    z = _G['tonelli_shanks']((n - 3) % n, n)   # sqrt(-3) mod n
    if z is None:
        return None
    if z % 2 == 0:                              # match parity of D = -3 (odd)
        z = n - z
    u, v = 2 * n, z
    limit = math.isqrt(4 * n)
    while v > limit:
        u, v = v, u % v
    t = v
    rem = 4 * n - t * t
    if rem <= 0 or rem % 3 != 0:
        return None
    s = math.isqrt(rem // 3)
    if s * s != rem // 3:
        return None
    # 4n = t^2 + 3s^2, and with u=2a-b, v=b: b = s, a = (t+s)/2
    if (t + s) % 2 != 0:
        t = -t
        if (t + s) % 2 != 0:
            return None
    a, b = (t + s) // 2, s
    if a * a - a * b + b * b != n:
        a, b = (-t + s) // 2, s
        if a * a - a * b + b * b != n:
            return None
    return (a, b)


def cmin_closed_form(n, lam, k2_bound=None, big=False):
    """cmin from the Eisenstein decomposition of n.

    n splits in Z[omega] as pi * pi_bar.  The lattice L_lam = {(x,y): x = lam*y
    mod n} is the ideal lattice of exactly ONE of the two primes, so only one of
    the two unit-orbits (associates of pi, associates of pi_bar) actually lies
    in L_lam.  Rather than track which, filter all 12 candidates by the defining
    congruence -- that is self-verifying and costs one modmul each.

    Returns ((cmin, e_star), (a,b)), with the box constraint |e| <= K2 applied.
    """
    if k2_bound is None:
        k2_bound = math.isqrt(n) + 1
    eis = cornacchia_d3(n) if big else eisenstein_decompose(n)
    if eis is None:
        return None, None
    a, b = eis
    cands = []
    for (x, y) in associates(a, b):
        if not (0 < abs(y) <= k2_bound):
            continue
        if (x - lam * y) % n and (x + lam * y) % n:
            continue                       # not in L_lam: wrong prime
        cands.append((abs(x), abs(y)))
    if not cands:
        return None, (a, b)
    return min(cands), (a, b)


print("=" * 78)
print("Thread 23e — cmin in closed form via the Eisenstein decomposition of n")
print("=" * 78)

# ---------------------------------------------------------------------------
# Part 1 — verify C-eis by brute force over all 12-bit and 16-bit GLV orders
# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("Part 1: brute-force cmin vs closed form, all j=0 GLV group orders")
print("-" * 78)

for lo, hi, label in ((2 ** 11, 2 ** 12, "12-bit"), (2 ** 15, 2 ** 16, "16-bit")):
    tested = agree = 0
    mismatches = []
    p = int(sympy.nextprime(lo))
    while p < hi:
        if p % 3 == 1:
            eis = eisenstein_decompose(p)
            if eis is not None:
                for t in j0_traces(*eis):
                    n = p + 1 - t
                    if n < 2 or n % 3 != 1 or not sympy.isprime(n):
                        continue
                    roots = glv_roots(n)
                    if roots is None:
                        continue
                    lam = roots[0]
                    k2 = math.isqrt(n) + 1
                    bf, e_star = cmin_bruteforce(n, lam, k2)
                    cf, ab = cmin_closed_form(n, lam, k2)
                    cf = cf[0] if cf else None
                    tested += 1
                    if bf == cf:
                        agree += 1
                    else:
                        mismatches.append((n, lam, bf, cf, ab, e_star))
        p = int(sympy.nextprime(p))
    print(f"{label}: {agree}/{tested} agree")
    for mm in mismatches[:6]:
        n, lam, bf, cf, ab, e_star = mm
        print(f"   MISMATCH n={n} lam={lam} brute={bf} (e*={e_star}) "
              f"closed={cf} eis={ab}")

# ---------------------------------------------------------------------------
# Part 2 — the real curves
# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("Part 2: cmin for deployed j=0 / GLV curves")
print("-" * 78)

CURVES = [
    ("secp256k1",
     0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141),
    ("BN254 (Fr)",
     0x30644E72E131A029B85045B68181585D2833E84879B9709143E1F593F0000001),
    ("BLS12-381 (Fr)",
     0x73EDA753299D7D483339D80809A1D80553BDA402FFFE5BFEFFFFFFFF00000001),
]

print(f"{'curve':<18} {'bits':>5} {'cmin':>42} {'cmin/sqrt(n)':>14}")
for label, n in CURVES:
    if n % 3 != 1:
        print(f"{label:<18} {n.bit_length():>5}   n != 1 mod 3 (no GLV eigenvalue)")
        continue
    lam = glv_roots(n)[0]
    res, ab = cmin_closed_form(n, lam, big=True)
    if res is None:
        print(f"{label:<18} {n.bit_length():>5}   no Eisenstein decomposition found")
        continue
    cf, e_star = res
    k2 = math.isqrt(n) + 1
    ok = (lam * e_star - cf) % n == 0 or (lam * e_star + cf) % n == 0
    r = cf / math.isqrt(n)
    print(f"{label:<18} {n.bit_length():>5} {cf:>42} {r:>14.6f}")
    print(f"{'':<18} {'':>5}   n = a^2-ab+b^2 with a={ab[0]}, b={ab[1]}")
    print(f"{'':<18} {'':>5}   attained at e*={e_star}  (e*/K2 = {e_star/k2:.6f})")
    print(f"{'':<18} {'':>5}   direct check lam*e* = +-cmin (mod n): {ok}")

print("\nInterpretation: the Thread-23d mechanism needs cmin < K1, and K1 at the")
print("information-theoretic wall is ~0.1*sqrt(n).  So a curve is anomalous only")
print("if cmin/sqrt(n) << 0.1.  Values of order 0.1-1.0 mean NO anomaly: the")
print("Phase-2 wall sits at the generic information-theoretic point.")

print("\n" + "=" * 78)
print("done")
print("=" * 78)
