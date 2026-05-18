"""
Biehl-Meyer-Mueller invalid-curve attack on real NIST P-224.

End-to-end demonstration of partial private-key recovery against a
vulnerable ECDH oracle (no on-curve validation). We scan invalid
curves E_{b'}, find small prime factors q in #E_{b'}, query the
oracle with points of order q, solve the DLP in each small subgroup
(BSGS), and CRT the residues.

Real-curve scaling is the same as the small demo — just with
PARI's Schoof-Elkies-Atkin for point counting instead of brute
enumeration. Each ellcard on P-224 takes ~3s.
"""

import time, math, random
from math import prod
from sympy import factorint, sqrt_mod
from cypari2 import Pari

pari = Pari()
pari.allocatemem(1 << 30)
random.seed(0xC0FFEE)

# ---- NIST P-224 parameters ----
p  = 2**224 - 2**96 + 1
a  = (-3) % p
b  = 0xB4050A850C04B3ABF54132565044B0B7D7BFD8BA270B39432355FFB4
n  = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFF16A2E0B8F03E13DD29455C5C2A3D
Gx = 0xB70E0CBD6BB4BF7F321390B94A03C1D356C21122343280D6115C1D21
Gy = 0xBD376388B5F723FB4C22DFE6CD4375A05A07476444D5819985007E34

# ---- EC arithmetic — b-free, as expected for the attack ----
def add(P, Q):
    if P is None: return Q
    if Q is None: return P
    x1, y1 = P; x2, y2 = Q
    if x1 == x2:
        if (y1 + y2) % p == 0: return None
        m = ((3 * x1 * x1 + a) * pow(2 * y1, -1, p)) % p
    else:
        m = ((y2 - y1) * pow((x2 - x1) % p, -1, p)) % p
    x3 = (m * m - x1 - x2) % p
    y3 = (m * (x1 - x3) - y1) % p
    return (x3, y3)

def neg(P):
    return None if P is None else (P[0], (-P[1]) % p)

def mul(k, P):
    R = None; Q = P
    while k > 0:
        if k & 1: R = add(R, Q)
        Q = add(Q, Q)
        k >>= 1
    return R

def card(b_prime):
    E = pari.ellinit([a, b_prime], p)
    return int(pari.ellcard(E))

def bsgs(Q, P, ord_P):
    """Solve k*P = Q for k in [0, ord_P)."""
    if Q is None: return 0
    m_step = int(math.isqrt(ord_P)) + 1
    baby = {None: 0}
    cur = P
    for j in range(1, m_step):
        baby[cur] = j
        cur = add(cur, P)
    giant = neg(mul(m_step, P))
    cur = Q
    for i in range(m_step + 1):
        if cur in baby:
            return (i * m_step + baby[cur]) % ord_P
        cur = add(cur, giant)
    raise RuntimeError("BSGS failed")

def find_point_of_order(b_prime, q, n_prime):
    cof = n_prime // q
    for _ in range(200):
        x = random.randrange(p)
        rhs = (x*x*x + a*x + b_prime) % p
        y = sqrt_mod(rhs, p)
        if y is None: continue
        P = (x, int(y))
        T = mul(cof, P)
        if T is not None and mul(q, T) is None:
            return T
    return None

# ---- Target ----
d = random.randrange(2, n)
print(f"[*] P-224 server private key:")
print(f"    d = 0x{d:056x}\n")

def vulnerable_oracle(Q):
    """ECDH responder: returns d*Q with NO check that Q lies on E_b."""
    return mul(d, Q)

# ---- Scan invalid curves ----
TARGET_BITS = 40
print(f"[*] Scanning invalid curves; target = recover {TARGET_BITS} bits of d.\n")

factors = {}
M = 1
t0 = time.time()
b_prime = 0
scanned = 0
useful = 0

while M.bit_length() < TARGET_BITS:
    b_prime += 1
    if b_prime == b: continue
    # Skip singular curves: Δ = -16(4a^3 + 27b'^2) = 0 mod p
    if (4 * pow(a, 3, p) + 27 * b_prime * b_prime) % p == 0: continue
    n_prime = card(b_prime)
    scanned += 1
    fs = factorint(n_prime, limit=10**5)
    small_qs = sorted(q for q in fs if q < 10**5 and q not in factors and q > 2)
    if not small_qs:
        continue
    progress = False
    for q in small_qs:
        Pq = find_point_of_order(b_prime, q, n_prime)
        if Pq is None: continue
        factors[q] = (b_prime, Pq, n_prime)
        M *= q
        progress = True
        print(f"    b'={b_prime:3d}  q={q:>6d} ({q.bit_length():2d}b)   "
              f"cumulative M = {M.bit_length():3d}b   "
              f"[{time.time()-t0:5.1f}s]")
        if M.bit_length() >= TARGET_BITS: break
    if progress: useful += 1

print(f"\n[*] {scanned} curves counted, {useful} useful, "
      f"{len(factors)} subgroups, M = {M.bit_length()} bits  "
      f"({time.time()-t0:.1f}s)\n")

# ---- Query oracle and solve DLPs ----
print("[*] Querying oracle and solving DLPs (BSGS) per subgroup...")
mods, rems = [], []
for q in sorted(factors.keys()):
    b_prime, Pq, _ = factors[q]
    R = vulnerable_oracle(Pq)
    k = bsgs(R, Pq, q)
    rems.append(k); mods.append(q)
    print(f"    via b'={b_prime:3d}:  d mod {q:>6d} = {k}")

# ---- CRT ----
d_rec = 0
for r, m in zip(rems, mods):
    Mi = M // m
    d_rec = (d_rec + r * Mi * pow(Mi, -1, m)) % M

print(f"\n[*] CRT result:    d ≡ 0x{d_rec:x}  (mod M, M is {M.bit_length()} bits)")
print(f"[*] True d mod M = 0x{d % M:x}")
print(f"[*] MATCH: {d_rec == d % M}")
print(f"\n[*] Bits recovered: {M.bit_length()} / {n.bit_length()}")
print(f"    To recover the full key, continue scanning until prod(q) > n.")
print(f"    SafeCurves' 2^58.4 figure is the rho cost on the LAST remaining")
print(f"    factor once cheap small-subgroup contributions are exhausted.")
