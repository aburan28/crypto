"""
Biehl-Meyer-Mueller invalid-curve attack on a P-224-shaped curve.

Scaled-down: same equation form as NIST P-224 (y^2 = x^3 - 3x + b),
small prime p so we can brute-count orders. Attack steps are identical
to the real-curve case.

Threat model: a server holds static private key d and offers an ECDH
oracle that returns d*Q for an attacker-supplied "public key" Q,
WITHOUT checking that Q lies on the real curve.
"""

import random
from math import prod
from sympy import isprime, sqrt_mod

random.seed(1)

p = 65521                  # prime, ~2^16
a = (-3) % p               # same coefficient as P-224

# ---- EC arithmetic on y^2 = x^3 - 3x + b ----
# Note: the formulas reference a (= -3) but NOT b. This is exactly why
# sending a point on an invalid curve E_b' still computes correctly
# on E_b' when the server thinks it's working on E_b.

def add(P, Q):
    if P is None: return Q
    if Q is None: return P
    x1, y1 = P; x2, y2 = Q
    if x1 == x2:
        if (y1 + y2) % p == 0:
            return None
        m = ((3 * x1 * x1 + a) * pow(2 * y1, -1, p)) % p
    else:
        m = ((y2 - y1) * pow((x2 - x1) % p, -1, p)) % p
    x3 = (m * m - x1 - x2) % p
    y3 = (m * (x1 - x3) - y1) % p
    return (x3, y3)

def mul(k, P):
    R = None; Q = P
    while k > 0:
        if k & 1: R = add(R, Q)
        Q = add(Q, Q)
        k >>= 1
    return R

# ---- Point counting via QR-table enumeration ----
qr = bytearray(p)
for i in range(1, p):
    qr[(i * i) % p] = 1

def count_points(b):
    n = 1
    for x in range(p):
        rhs = (x*x*x + a*x + b) % p
        if rhs == 0:   n += 1
        elif qr[rhs]:  n += 2
    return n

# ---- Build target ----
def find_prime_order_b():
    for b in range(2, p):
        n = count_points(b)
        if isprime(n):
            return b, n
    raise RuntimeError("no prime-order curve found")

print("[*] Building target curve (P-224 analog)...")
b_real, n_real = find_prime_order_b()
print(f"    E: y^2 = x^3 - 3x + {b_real}  (mod {p})")
print(f"    #E = {n_real}  (prime, {n_real.bit_length()} bits)")

d = random.randrange(2, n_real)
print(f"[*] Server's static private key:  d = {d}")

def vulnerable_oracle(Q):
    """Returns d*Q. Does NOT verify Q lies on E_{b_real}."""
    return mul(d, Q)

# ---- Attacker scans invalid curves for small-subgroup factors ----
print("\n[*] Scanning invalid curves E_{b'} for small order factors...")
small_primes = [2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97]
factors = {}

def find_point_of_order(b_prime, q, n_prime):
    cof = n_prime // q
    for _ in range(200):
        x = random.randrange(p)
        rhs = (x*x*x + a*x + b_prime) % p
        if rhs != 0 and not qr[rhs]: continue
        y = sqrt_mod(rhs, p)
        if y is None: continue
        P = (x, int(y))
        Pq = mul(cof, P)
        if Pq is not None and mul(q, Pq) is None:
            return Pq
    return None

for b_prime in range(1, p):
    if b_prime == b_real: continue
    if prod(factors.keys()) >= n_real: break
    n_prime = count_points(b_prime)
    for q in small_primes:
        if q in factors: continue
        if n_prime % q != 0: continue
        Pq = find_point_of_order(b_prime, q, n_prime)
        if Pq is None: continue
        factors[q] = (b_prime, Pq, n_prime)
        print(f"    invalid b'={b_prime:5d}  |E_b'|={n_prime:6d}  "
              f"-> point of order {q}")
        if prod(factors.keys()) >= n_real: break

M = prod(factors.keys())
print(f"\n[*] Collected subgroups; combined modulus M = {M}  "
      f"({M.bit_length()} bits, key is {n_real.bit_length()} bits)")

# ---- Query oracle, do DLP in each small subgroup ----
print("\n[*] Querying oracle and solving DLP in each small subgroup...")
mods, rems = [], []
for q in sorted(factors.keys()):
    b_prime, Pq, _ = factors[q]
    R = vulnerable_oracle(Pq)        # = d * Pq, computed on E_{b'}
    # Brute-force the unique k in [0, q) with k*Pq == R
    acc = None
    for k in range(q):
        if acc == R:
            rems.append(k); mods.append(q)
            print(f"    via b'={b_prime:5d}:  d mod {q:2d} = {k}")
            break
        acc = add(acc, Pq)
    else:
        raise RuntimeError(f"DLP failed for q={q}")

# ---- CRT ----
d_rec = 0
for r, m in zip(rems, mods):
    Mi = M // m
    d_rec = (d_rec + r * Mi * pow(Mi, -1, m)) % M

print(f"\n[*] CRT result:    d ≡ {d_rec}  (mod {M})")
print(f"[*] True d mod M = {d % M}")
print(f"[*] Match: {d_rec == d % M}")
if M >= n_real:
    print(f"[*] FULL KEY RECOVERED: d = {d_rec}")
else:
    print(f"[*] Recovered {M.bit_length()} of {n_real.bit_length()} bits; "
          f"remaining {(n_real // M).bit_length()} bits via Pollard kangaroo.")
