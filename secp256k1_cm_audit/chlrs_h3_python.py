"""
chlrs_h3_python.py

Analyse Howe (H3) condition gcd(#E1,#E2) vs Rosenhain formula success.
Uses a pure-Python approach (no fpylll/PARI required) for counting points
and factoring — for small primes only.

Key question:
  For j=0 curves y^2=x^3+b over F_p (p≡1 mod 3),
  does gcd(#E1, #E2) = 1 predict whether the Rosenhain formula for
  the Howe-glued curve works?
"""

def is_prime(n):
    if n < 2: return False
    if n == 2: return True
    if n % 2 == 0: return False
    for i in range(3, int(n**0.5)+1, 2):
        if n % i == 0: return False
    return True

def primes_up_to(n):
    return [p for p in range(2, n+1) if is_prime(p)]

def pow_mod(a, e, m):
    return pow(a, e, m)

def modinv(a, p):
    # Extended Euclid
    return pow(a, p-2, p)

def kronecker(a, p):
    """Legendre symbol (a/p) for odd prime p"""
    if a % p == 0: return 0
    return 1 if pow_mod(a % p, (p-1)//2, p) == 1 else -1

def is_cubic_residue(a, p):
    return pow_mod(a % p, (p-1)//3, p) == 1

def cube_roots_mod_p(b, p):
    """Return list of x with x^3 ≡ -b (mod p), i.e. roots of x^3+b=0."""
    target = (-b) % p
    # First find one cube root of target
    # For p ≡ 1 mod 3, use Tonelli-Shanks variant for cube roots
    roots = [x for x in range(p) if pow_mod(x, 3, p) == target]
    return roots

def count_points_j0(b, p):
    """Count #E(F_p) for E: y^2 = x^3 + b."""
    count = 1  # point at infinity
    for x in range(p):
        rhs = (x**3 + b) % p
        if rhs == 0:
            count += 1
        elif pow_mod(rhs, (p-1)//2, p) == 1:
            count += 2
    return count

def mobius_image(x_val, r1, r2, s1, p):
    """T(x) = (x-r1)/(x-r2) * (s1-r2)/(s1-r1) mod p"""
    num = (x_val - r1) % p
    den = (x_val - r2) % p
    scale_num = (s1 - r2) % p
    scale_den = (s1 - r1) % p
    return (num * modinv(den, p) % p) * (scale_num * modinv(scale_den, p) % p) % p

def hyperelliptic_charpoly_frobenius_check(lambdas, p, t1):
    """
    Count #Jac(C) for C: y^2=x(x-1)(x-l1)(x-l2)(x-l3) over F_p.
    Compare to expected (p+1-t1)(p+1+t1).
    """
    l1, l2, l3 = lambdas
    # Count points on genus-2 curve over F_p: #C(F_p) = sum over x of (1 + legendre(f(x)/p))
    # Use #Jac via counting on C and applying genus-2 formula.
    # For simplicity, count #C(F_p) directly:
    count_C = 2  # two points at infinity for genus-2 in some models, but Rosenhain: 1 inf pt + finite pts
    # Actually for y^2 = f(x) with f degree 5, there's 1 point at infinity.
    count_C = 1  # point at infinity
    for xi in range(p):
        rhs = xi * ((xi - 1) % p) % p * ((xi - l1) % p) % p * ((xi - l2) % p) % p * ((xi - l3) % p) % p
        rhs %= p
        if rhs == 0:
            count_C += 1
        elif pow_mod(rhs, (p-1)//2, p) == 1:
            count_C += 2
    target_n = (p + 1 - t1) * (p + 1 + t1)
    # We can't easily get #Jac from #C directly without more arithmetic.
    # Instead use: #Jac(C) = p^2 + c*p + ... from char poly, but we need the char poly.
    # Simpler: check if #Jac = target_n using Weil's formulas for genus-2:
    # #Jac(F_p) = p^2 - (s2-s1^2/2)*p + ... this is complex.
    # Use: for genus-2 over F_p, if char poly is T^4 + aT^3 + bT^2 + aT + 1 (after normalizing),
    # then #Jac = char_poly(1).
    # We DON'T have the char poly without PARI, so let's skip direct Jac count.
    # Instead, just verify #C matches a split Jacobian expectation.
    return count_C, target_n

def main():
    print("p    b1  d   #E1  #E2  gcd  #E1*#E2  #rts1  #rts2  H3?")
    print("-" * 60)
    
    results = []
    for p in primes_up_to(100):
        if p % 3 != 1:
            continue
        # Find cube root of unity mod p (to identify omega)
        # omega^2 + omega + 1 = 0 mod p
        omega = None
        for w in range(2, p):
            if (w*w + w + 1) % p == 0:
                omega = w
                break
        if omega is None:
            continue
        
        for b1 in range(1, 8):
            # Check if x^3+b1 has 3 roots mod p (i.e., -b1 is a cube mod p)
            if (p-1) % 3 != 0:
                continue
            if not is_cubic_residue(-b1 % p, p):
                continue
            roots1 = cube_roots_mod_p(b1, p)
            if len(roots1) != 3:
                continue
            
            # Find first non-square d
            d = None
            for dd in range(2, p):
                if kronecker(dd, p) == -1:
                    d = dd
                    break
            if d is None:
                continue
            
            b2 = (b1 * d * d * d) % p
            roots2 = cube_roots_mod_p(b2, p)
            if len(roots2) != 3:
                continue
            
            n1 = count_points_j0(b1, p)
            n2 = count_points_j0(b2, p)
            t1 = p + 1 - n1
            import math
            g = math.gcd(n1, n2)
            
            print(f"  {p:<4} {b1:<3} {d:<3} {n1:<5} {n2:<5} {g:<4} {n1*n2:<8} {len(roots1):<6} {len(roots2):<6} {'H3:ok' if g==1 else 'H3:FAIL'}")
            results.append((p, b1, d, n1, n2, g, t1))
    
    print()
    print("Summary:")
    h3_ok = [(r[0], r[1], r[2], r[3], r[4], r[5]) for r in results if r[5] == 1]
    h3_fail = [(r[0], r[1], r[2], r[3], r[4], r[5]) for r in results if r[5] != 1]
    print(f"  H3 satisfied (gcd=1): {len(h3_ok)} cases")
    print(f"  H3 failed  (gcd>1):   {len(h3_fail)} cases")
    print()
    print("H3-satisfied cases (candidates for Rosenhain to work):")
    for r in h3_ok:
        print(f"    p={r[0]}, b1={r[1]}, d={r[2]}: #E1={r[3]}, #E2={r[4]}, gcd={r[5]}")

if __name__ == "__main__":
    main()
