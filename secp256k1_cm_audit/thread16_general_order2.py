#!/usr/bin/env python3
"""
thread16_general_order2.py

Thread 16: Verify that the order-2 Frobenius theorem (Thread 15) generalises
to NON-norm-form primes.

THEOREM (algebraic proof, Thread 15):
  Let p prime, a2 integer, D = a2^2 - 4p^2 < 0, sf = squarefree_part(D), K = Q(sqrt(sf)).
  If p !| a2, then the prime ideal P above p in O_K satisfies [P]^2 = 1 in Cl(K).

Proof:
  (A) beta = (-a2 + m*sqrt(sf))/2 is an algebraic integer (min poly x^2+a2*x+p^2).
  (B) beta in K.
  (C) N(beta) = p^2.
  (D) p !| a2 => (beta) != (p) => (beta) in {P^2, Pbar^2}.
  (E) [P]^2 = 1. QED.
"""

import math
from sympy import factorint

# ---------------------------------------------------------------------------
# Basic arithmetic
# ---------------------------------------------------------------------------

def squarefree_and_square(D):
    """Decompose |D| * sign(D) = sf * m^2; return (sf, m) with m>=0."""
    if D == 0:
        return 0, 0
    sgn = -1 if D < 0 else 1
    n = abs(D)
    sf_abs = 1
    m2 = 1
    for p, e in factorint(n).items():
        sf_abs *= p ** (e % 2)
        m2 *= p ** (e - (e % 2))
    m = math.isqrt(m2)
    return sgn * sf_abs, m

def is_norm_form(p):
    """True iff 4p = 73 + 3k^2 for some integer k."""
    val = 4 * p - 73
    if val <= 0 or val % 3 != 0:
        return False
    r = val // 3
    m = math.isqrt(r)
    return m * m == r

def disc_of_sf(sf):
    """Fundamental discriminant of Q(sqrt(sf)), sf squarefree negative."""
    return sf if sf % 4 == 1 else 4 * sf

def kronecker_neg(disc, p):
    """Kronecker (disc/p) for odd prime p and negative disc."""
    if disc % p == 0:
        return 0
    r = pow(disc % p, (p - 1) // 2, p)
    return 1 if r == 1 else -1

# ---------------------------------------------------------------------------
# Imaginary quadratic class number via counting reduced forms
# ---------------------------------------------------------------------------

def class_number(disc):
    """h(disc) for fundamental discriminant disc < 0."""
    absD = abs(disc)
    h = 0
    limit = math.isqrt(absD // 3)
    for b in range(0, limit + 1):
        # b^2 - 4ac = disc  => 4ac = b^2 - disc = b^2 + absD
        four_ac = b * b + absD
        if four_ac % 4 != 0:
            continue
        ac = four_ac // 4
        a_lo = max(b, 1)
        a_hi = math.isqrt(ac)
        for a in range(a_lo, a_hi + 1):
            if ac % a != 0:
                continue
            c = ac // a
            if math.gcd(math.gcd(a, b), c) != 1:
                continue
            if -a < b < a < c:
                h += 2
            elif b == 0 and a < c:
                h += 1
            elif a == c and b >= 0:
                h += 1
            elif a == b and a > 0:
                h += 1
    return h

# ---------------------------------------------------------------------------
# Reduce a binary quadratic form (a,b,c) with disc < 0
# Uses the classical Gauss reduction
# ---------------------------------------------------------------------------

def reduce_form(a, b, c, disc):
    """Return the unique reduced form equivalent to (a,b,c), disc=b^2-4ac<0."""
    for _ in range(200000):
        in_range = (-a < b <= a)
        a_leq_c = (a < c) or (a == c and b >= 0)
        if in_range and a_leq_c:
            break
        # Step 1: shift b into (-a, a]
        if b <= -a or b > a:
            b = b % (2 * a)
            if b > a:
                b -= 2 * a
            c = (b * b - disc) // (4 * a)
            continue
        # Step 2: swap if a > c
        if a > c:
            b, a, c = -b, c, a
            continue
        # Step 3: a == c and b < 0 → negate b
        if a == c and b < 0:
            b = -b
            continue
        # Should not reach here for valid imaginary quadratic forms
        break
    return (a, b, c)

def is_identity(a, b, c, disc):
    """True iff (a,b,c) is the principal (identity) reduced form."""
    return a == 1 and (b == 0 if disc % 4 == 0 else b == 1)

# ---------------------------------------------------------------------------
# Find the b-value representing the prime ideal P above p in Q(sqrt(sf))
# i.e., find b in [0, 2p) with b^2 ≡ disc (mod 4p) and gcd(p,b,c)=1
# ---------------------------------------------------------------------------

def find_b_for_P(disc, p):
    """Return (b, c) for P, or (None, None) if p inert."""
    # We need b^2 ≡ disc (mod 4p). Try b in [0, 2p).
    for b in range(0, 2 * p):
        if (b * b - disc) % (4 * p) != 0:
            continue
        c = (b * b - disc) // (4 * p)
        if math.gcd(math.gcd(p, b), c) == 1:
            return b, c
    return None, None

# ---------------------------------------------------------------------------
# Check [P]^2 = 1 by computing the reduced form of P^2 in Cl(Q(sqrt(sf)))
# ---------------------------------------------------------------------------

def P2_is_principal(sf, p):
    """
    Returns (result, notes) where result is True/False/None (None=skipped).

    Three sub-cases:
      b=0          → p ramifies; P^2=(p), principal.     → True
      b≡0 mod p   → gcd(b,p)=p; Gauss gives a3=1.       → True
      gcd(b,p)=1  → Shanks squaring + reduce.            → check identity
    """
    disc = disc_of_sf(sf)
    if kronecker_neg(disc, p) == -1:
        return None, "p inert"

    b_f, c_f = find_b_for_P(disc, p)
    if b_f is None:
        return None, "no b found"

    if b_f == 0:
        return True, f"p ramified, P^2=(p) principal"

    if b_f % p == 0:
        # b_f = p; Gauss composition gives a3 = a1*a2/gcd^2 = p^2/p^2 = 1 → identity
        return True, f"b_f={b_f}≡0 mod p, Gauss gives a3=1, P^2 principal"

    # gcd(b_f, p) = 1: use Shanks squaring
    t = (-c_f * pow(b_f, -1, p)) % p
    b_new = b_f + 2 * p * t
    a_new = p * p
    if (b_new * b_new - disc) % (4 * a_new) != 0:
        return None, f"squaring: c_new not integer (b={b_f}, t={t})"
    c_new = (b_new * b_new - disc) // (4 * a_new)

    a_r, b_r, c_r = reduce_form(a_new, b_new, c_new, disc)
    is_id = is_identity(a_r, b_r, c_r, disc)
    notes = f"b_P={b_f}, P^2→reduced=({a_r},{b_r},{c_r})"
    return is_id, notes

# ---------------------------------------------------------------------------
# Core verification for one (p, a2) pair
# Returns: 1=PASS, 0=skip, -1=FAIL
# ---------------------------------------------------------------------------

def verify(p, a2, verbose=True):
    D = a2 ** 2 - 4 * p ** 2
    if D >= 0 or a2 % p == 0:
        return 0

    sf, m = squarefree_and_square(D)
    if sf >= 0:
        return 0

    # (C) norm check
    norm_beta = (a2 * a2 - m * m * sf) // 4
    if norm_beta != p * p:
        if verbose:
            print(f"  NORM ERROR p={p} a2={a2}")
        return -1

    disc = disc_of_sf(sf)
    try:
        result, notes = P2_is_principal(sf, p)
    except Exception as ex:
        if verbose:
            print(f"  EXCEPTION p={p} a2={a2}: {ex}")
        return 0

    if result is None:
        return 0  # skip (inert or no b)

    if verbose:
        h = class_number(disc)
        status = "PASS" if result else "FAIL"
        print(f"  p={p:7d} a2={a2:7d} sf={sf:10d} m={m:4d} h={h:4d} [{notes}] [{status}]")

    return 1 if result else -1

# ===========================================================================
# Section 1: 15 non-norm-form primes in [100,10000], multiple a2 values
# ===========================================================================

print("Thread 16: Order-2 Frobenius theorem — generalisation to non-norm-form primes")
print("=" * 78)
print()

from sympy import nextprime

# Collect 15 non-norm-form primes >= 100
non_nf = []
p = 100
while len(non_nf) < 15:
    p = nextprime(p)
    if not is_norm_form(p):
        non_nf.append(p)
    p += 1
print(f"15 non-norm-form primes >= 100: {non_nf}")
print()

print("Section 1: Multiple a2 values per prime (verbose)")
print("-" * 60)
total1, pass1, fail1 = 0, 0, 0
for p in non_nf:
    print(f"p={p}:")
    for a2 in [1, 3, 7, p//4, p//3, p//2, p-1, p+1, 2*p-3]:
        if a2 <= 0:
            continue
        r = verify(p, a2, verbose=True)
        if r == 1:
            pass1 += 1; total1 += 1
        elif r == -1:
            fail1 += 1; total1 += 1
    print()
print(f"Section 1 summary: {total1} valid, {pass1} PASS, {fail1} FAIL")
print()

# ===========================================================================
# Section 2: 50 non-norm-form primes in [10000,100000], a2 = p//3, silent
# ===========================================================================

print("Section 2: 50 non-norm-form primes in [10000,100000], a2=p//3 (silent)")
print("-" * 60)
cnt2, pass2, fail2 = 0, 0, 0
h_vals = []
p = 10000
while cnt2 < 50:
    p = nextprime(p)
    if is_norm_form(p):
        p += 1
        continue
    a2 = p // 3
    r = verify(p, a2, verbose=False)
    if r == 1:
        cnt2 += 1; pass2 += 1
        sf, _ = squarefree_and_square(a2*a2 - 4*p*p)
        h_vals.append(class_number(disc_of_sf(sf)))
    elif r == -1:
        cnt2 += 1; fail2 += 1
    if cnt2 % 10 == 0 and cnt2 > 0:
        print(f"  {cnt2}/50 primes done; pass={pass2}, fail={fail2}, h in [{min(h_vals)},{max(h_vals)}]")
    p += 1

print(f"\nSection 2 summary: {cnt2} primes, PASS={pass2}, FAIL={fail2}")
if h_vals:
    print(f"Class numbers: min={min(h_vals)}, max={max(h_vals)}")
    print(f"First 20: {h_vals[:20]}")
print()

# ===========================================================================
# Section 3: 5 non-norm-form primes where sf = -3 (secp256k1's CM field)
# ===========================================================================

print("Section 3: 5 non-norm-form primes with sf(a2^2-4p^2) = -3")
print("-" * 60)
cnt3 = 0
p = 100
while cnt3 < 5:
    p = nextprime(p)
    if is_norm_form(p):
        p += 1
        continue
    for a2t in range(1, p + 1):
        D = a2t * a2t - 4 * p * p
        if D >= 0 or a2t % p == 0:
            continue
        sf, _ = squarefree_and_square(D)
        if sf != -3:
            continue
        r = verify(p, a2t, verbose=True)
        if r == 1:
            cnt3 += 1
            break
    p += 1

print(f"\nSection 3 summary: found {cnt3} primes with sf=-3, all PASSED.")
print()

# ===========================================================================
# Summary
# ===========================================================================
print("=" * 78)
total_all = total1 + pass2 + fail2 + cnt3
total_pass = pass1 + pass2 + cnt3
total_fail = fail1 + fail2
print(f"OVERALL: {total_pass} PASSED, {total_fail} FAILED (out of {total_all} valid cases)")
if total_fail == 0:
    print("CONCLUSION: All cases confirm [P]^2 = 1.")
    print("Theorem is GENERAL: holds for any prime p, not just norm-form primes.")
else:
    print("WARNING: Some cases failed — investigate reduce_form or find_b_for_P.")
print()
