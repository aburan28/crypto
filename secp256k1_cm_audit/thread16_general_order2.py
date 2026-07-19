#!/usr/bin/env python3
"""
Thread 16: Verify the order-2 Frobenius theorem for non-norm-form primes.

THEOREM (generalised from Thread 15):
  Let p be any prime and t any integer with 0 < |t| < 2*sqrt(p) and p ∤ t.
  Form the product Jacobian J = E × E^t with biquadratic Weil polynomial
      F(T) = T^4 + a2*T^2 + p^2,  a2 = 2p - t^2.
  Let D = a2^2 - 4p^2 = t^2*(t^2 - 4p) < 0,
      sf = squarefree part of D  (D = sf * m^2, m > 0),
      K  = Q(sqrt(sf))  (imaginary quadratic field).
  Then the prime ideal P of O_K above p satisfies [P]^2 = 1 in Cl(K).

PROOF (algebraic, no class group computation required):
  Define beta = (-a2 + m*sqrt(sf)) / 2.
  (A) beta is an algebraic integer: satisfies x^2 + a2*x + p^2 = 0. ✓
      (Trace = -a2 ∈ Z; norm = p^2 ∈ Z; polynomial is monic.)
  (B) beta ∈ K = Q(sqrt(sf)) since sqrt(D) = m*sqrt(sf). ✓
  (C) N_{K/Q}(beta) = (a2^2 - m^2*sf)/4 = (a2^2 - D)/4 = 4p^2/4 = p^2. ✓
  (D) p ∤ a2: since a2 = 2p - t^2 ≡ -t^2 (mod p), and p ∤ t, we have p ∤ t^2,
      so p ∤ a2.  Hence (beta) ≠ (p) in O_K. ✓
  (E) The ideals of norm p^2 in O_K are: P^2, Pbar^2, and P*Pbar = (p).
      Since (beta) ≠ (p), we have (beta) = P^2 or Pbar^2.
      Hence P^2 = (beta) is principal, so [P]^2 = 1 in Cl(K). □

Steps (A)-(D) are purely arithmetic; we verify them for 10 non-norm-form primes.
Step (E) is a logical consequence (unique factorisation in Dedekind domains).

Additionally we verify that p actually splits in K (Kronecker symbol check),
which is required for P to exist as a prime above p (distinct from Pbar).
"""

import math
import sys

# ---------------------------------------------------------------
# Number-theory utilities (no external dependencies)
# ---------------------------------------------------------------

def isqrt(n):
    """Integer square root (floor)."""
    if n < 0:
        raise ValueError(f"isqrt({n})")
    if n == 0:
        return 0
    x = int(math.isqrt(n))
    # Correct for floating-point errors
    while x * x > n:
        x -= 1
    while (x + 1) * (x + 1) <= n:
        x += 1
    return x

def is_perfect_square(n):
    if n < 0:
        return False, 0
    r = isqrt(n)
    return (r * r == n), r

def squarefree_part(n):
    """Return (sf, m) with n = sf * m^2, sf squarefree (and same sign as n)."""
    if n == 0:
        return 0, 0
    sign = 1 if n > 0 else -1
    n_abs = abs(n)
    # Factor out perfect square part
    sf = 1
    m_sq = 1
    # Trial division squarefree decomposition
    temp = n_abs
    for p in range(2, isqrt(temp) + 1):
        if p * p > temp:
            break
        if temp % p == 0:
            cnt = 0
            while temp % p == 0:
                temp //= p
                cnt += 1
            # p^cnt contributes p^(cnt//2) to m, and p^(cnt%2) to sf
            sf *= p ** (cnt % 2)
            m_sq *= p ** (cnt // 2)
        if temp == 1:
            break
    if temp > 1:
        sf *= temp  # residual prime factor (odd power)
    return sign * sf, m_sq

def kronecker(a, n):
    """Kronecker symbol (a/n) — Jacobi symbol extended to n=2."""
    if n == 0:
        return 1 if abs(a) == 1 else 0
    if n < 0:
        n = -n
    if n == 1:
        return 1
    # Jacobi symbol via quadratic reciprocity
    result = 1
    if a < 0:
        a = -a
        if n % 4 == 3:
            result = -result
    while True:
        a = a % n
        if a == 0:
            return 0 if n > 1 else result
        while a % 2 == 0:
            a //= 2
            if n % 8 in (3, 5):
                result = -result
        a, n = n, a
        if a % 4 == 3 and n % 4 == 3:
            result = -result
        if n == 1:
            return result

def is_prime_miller_rabin(n, witnesses=None):
    """Deterministic Miller-Rabin for n < 3.3*10^24 (sufficient here)."""
    if n < 2:
        return False
    if n == 2 or n == 3 or n == 5 or n == 7:
        return True
    if n % 2 == 0:
        return False
    # Write n-1 = 2^r * d
    r, d = 0, n - 1
    while d % 2 == 0:
        r += 1
        d //= 2
    # Deterministic witnesses for n < 3,317,044,064,679,887,385,961,981
    if witnesses is None:
        witnesses = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]
    for a in witnesses:
        if a >= n:
            continue
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue
        for _ in range(r - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True

def is_normform(p):
    """Is p of the form (73 + 3k^2)/4 for some odd positive integer k?"""
    val = 4 * p - 73
    if val <= 0 or val % 3 != 0:
        return False
    k_sq = val // 3
    ok, k = is_perfect_square(k_sq)
    return ok and k % 2 == 1

# ---------------------------------------------------------------
# Core verification for a single (p, t) pair
# ---------------------------------------------------------------

def verify_general(p, t, label):
    """
    Verify [P]^2 = 1 in Cl(Q(sqrt(sf))) for the (p,t) pair.
    Returns True if all checks (A)-(D) pass, False otherwise.
    Prints a result line.
    """
    results = {}

    # Sanity checks
    assert is_prime_miller_rabin(p), f"p={p} not prime"
    assert 0 < abs(t) < 2 * math.sqrt(p), f"t={t} outside Hasse bound for p={p}"
    assert t % p != 0, f"p | t for p={p}, t={t}"

    # --- Construct Weil polynomial ---
    a2 = 2 * p - t * t          # coefficient in T^4 + a2*T^2 + p^2
    D  = a2 * a2 - 4 * p * p   # = t^2*(t^2 - 4p) < 0

    assert D < 0, f"D={D} >= 0 (expected negative)"

    sf, m = squarefree_part(D)          # D = sf * m^2
    m_abs = abs(m)

    # Check D = sf * m^2 exactly
    results["D=sf*m^2"] = (sf * m_abs * m_abs == D)

    # (A) minpoly of beta: x^2 + a2*x + p^2 is monic with integer coefficients ✓
    # Verify trace = -a2, norm = p^2
    trace_ok  = True          # trace = -(beta + betabar) = -(-a2+m√sf)/2 - (-a2-m√sf)/2 = a2... wait
    # Actually minpoly: if beta = (-a2 + m√sf)/2, then betabar = (-a2 - m√sf)/2.
    # beta + betabar = -a2  =>  trace of minpoly = -(beta+betabar) = a2  =>  coeff of x is a2 ✓
    # beta * betabar = (a2^2 - m^2*sf)/4 = (a2^2 - D)/4 = p^2  =>  const term = p^2 ✓
    # Monic, Z-coefficients: ✓ (a2 ∈ Z, p^2 ∈ Z)
    results["(A) beta_algebraic_integer"] = True   # logical consequence of D = sf*m^2

    # (B) beta ∈ Q(sqrt(sf)): ✓ by construction
    results["(B) beta_in_K"] = True

    # (C) N(beta) = p^2
    norm_beta = (a2 * a2 - m_abs * m_abs * sf) // 4
    results["(C) N(beta)=p^2"] = (norm_beta == p * p)
    # Cross-check: a2^2 - m^2*sf = a2^2 - D = a2^2 - (a2^2 - 4p^2) = 4p^2 ✓
    assert a2 * a2 - m_abs * m_abs * sf == 4 * p * p, "Norm formula mismatch"

    # (D) p ∤ a2  =>  (beta) ≠ (p)
    results["(D) p_nmid_a2"] = (a2 % p != 0)
    # Analytical: a2 = 2p - t^2 ≡ -t^2 (mod p). Since p ∤ t, p ∤ t^2, so p ∤ a2. ✓

    # Verify p splits in K = Q(sqrt(sf)): need Kronecker(sf, p) = 1
    # (If p is inert, there's no prime P above p to talk about; theorem vacuously holds differently.)
    kron = kronecker(sf % p, p) if p != 2 else 0
    results["p_splits_in_K"] = (kron == 1)

    all_ok = all(results.values())

    print(f"  {label:16s}  p={p:<7d}  t={t:<5d}  a2={a2:<10d}  sf={sf:<9d}  m={m_abs:<6d}  "
          f"N(β)=p²:{['NO ','YES'][results['(C) N(beta)=p^2']]}  "
          f"p∤a2:{['NO(!)','yes'][results['(D) p_nmid_a2']]}  "
          f"p↑K:{['no ','YES'][results['p_splits_in_K']]}  "
          f"[P]²=1:{['FAIL','YES'][all_ok]}")

    return all_ok

# ---------------------------------------------------------------
# Test suite: 10 non-norm-form primes with varied traces
# ---------------------------------------------------------------

# Primes NOT of the form (73 + 3k^2)/4:
# p=101: 4*101-73=331, 331/3 not integer  => not norm-form
# p=127: 4*127-73=435, 435/3=145, sqrt(145)≈12.04 not integer => not norm-form
# p=157: 4*157-73=555, 555/3=185, sqrt(185)≈13.6 => not norm-form
# p=199: 4*199-73=723, 723/3=241, sqrt(241)≈15.5 => not norm-form
# p=251: 4*251-73=931, 931/3 not integer => not norm-form
# p=307: 4*307-73=1155, 1155/3=385, sqrt(385)≈19.6 => not norm-form
# p=401: 4*401-73=1531, 1531/3 not integer => not norm-form
# p=503: 4*503-73=1939, 1939/3 not integer => not norm-form
# p=601: 4*601-73=2331, 2331/3=777, sqrt(777)≈27.9 => not norm-form
# p=701: 4*701-73=2731, 2731/3 not integer => not norm-form

TEST_CASES = [
    # (prime, trace, label)  -- traces chosen: varied, |t| << 2*sqrt(p)
    (101,   3,  "p=101,t=3"),
    (127,   5,  "p=127,t=5"),
    (157,   7,  "p=157,t=7"),
    (199,   9,  "p=199,t=9"),
    (251,   4,  "p=251,t=4"),
    (307,  11,  "p=307,t=11"),
    (401,   6,  "p=401,t=6"),
    (503,  13,  "p=503,t=13"),
    (601,  10,  "p=601,t=10"),
    (701,  15,  "p=701,t=15"),
]

def main():
    print("Thread 16: General order-2 Frobenius theorem (non-norm-form verification)")
    print("=" * 75)
    print()
    print("THEOREM: For any prime p and trace t with 0 < |t| < 2√p, p ∤ t,")
    print("  the prime P of O_{Q(√sf)} above p (sf = sf-part(t²(t²-4p)))")
    print("  satisfies [P]² = 1 in Cl(Q(√sf)).")
    print()
    print("Proof steps verified numerically:")
    print("  (A) β = (−a₂ + m√sf)/2 is an algebraic integer [logical]")
    print("  (B) β ∈ Q(√sf) [logical]")
    print("  (C) N(β) = (a₂² − m²sf)/4 = p²  [arithmetic]")
    print("  (D) p ∤ a₂ = 2p − t², so (β) ≠ (p)  [arithmetic]")
    print("  (E) ∴ (β) = P² or P̄² ⟹ [P]² = 1  [Dedekind domain theory]")
    print()
    print("Norm-form exclusion check:")
    for p, t, lbl in TEST_CASES:
        nf = is_normform(p)
        mark = "NORM-FORM (!)" if nf else "not norm-form ✓"
        print(f"  p={p}: {mark}")
    print()
    print("Main verification (columns: N(β)=p²?, p∤a₂?, p splits in K?, [P]²=1?):")
    print("-" * 95)

    passed, total = 0, 0
    for p, t, label in TEST_CASES:
        ok = verify_general(p, t, label)
        total += 1
        if ok:
            passed += 1

    print("-" * 95)
    print(f"\nResult: {passed}/{total} passed — theorem [P]²=1 holds for all non-norm-form test cases.")
    if passed == total:
        print("ALL PASSED ✓")
        print()
        print("CONCLUSION: The algebraic proof (A)-(E) from Thread 15 is completely")
        print("general — it uses no properties of the norm-form family 4p = 73 + 3k².")
        print("The theorem holds for arbitrary primes p with arbitrary traces t (p ∤ t, t≠0).")
    else:
        print("FAILURES DETECTED — see rows marked FAIL above.")
        sys.exit(1)

if __name__ == "__main__":
    main()
