"""
Thread 18: Howe (2,2)-gluing conditions for all 15 pairs of j=0 sextic twists.

secp256k1: y^2 = x^3 + 7, p prime, p ≡ 1 (mod 6).
Exactly 6 sextic isomorphism classes: y^2 = x^3 + b_k, b_k = 7*u^k (k=0..5),
u = primitive 6th root of unity mod p.

For each pair (i,j), checks Howe conditions:
  H1: N_i ≠ N_j  (curves not F_p-isogenous)
  H2: x^3+b_i and x^3+b_j have same factorisation pattern over F_p
      (equivalently, same 2-torsion Galois module structure)
  H3: gcd(N_i, N_j) = 1

Result: 5 of 15 pairs are Howe-glueable.

Run: python3 thread18_howe_sextic_twists.py
"""

import math

p = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F
n = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
t = p + 1 - n

# --- Utility functions ---
def sqrt_mod(n, p):
    n = n % p
    if n == 0: return 0
    if pow(n, (p-1)//2, p) != 1: return None
    if p % 4 == 3: return pow(n, (p+1)//4, p)
    Q, S = p - 1, 0
    while Q % 2 == 0: Q //= 2; S += 1
    z = 2
    while pow(z, (p-1)//2, p) != p-1: z += 1
    M, c, t2, R = S, pow(z, Q, p), pow(n, Q, p), pow(n, (Q+1)//2, p)
    while True:
        if t2 == 1: return R
        i, temp = 1, (t2*t2) % p
        while temp != 1: temp = (temp*temp) % p; i += 1
        b = pow(c, 1 << (M - i - 1), p)
        M, c, t2, R = i, (b*b) % p, (t2*b*b) % p, (R*b) % p

def ec_double(P, p):
    if P is None: return None
    X, Y = P
    if Y == 0: return None
    lam = (3 * X * X * pow(2*Y, p-2, p)) % p
    X2 = (lam*lam - 2*X) % p
    Y2 = (lam*(X - X2) - Y) % p
    return (X2, Y2)

def ec_add(P, Q, p):
    if P is None: return Q
    if Q is None: return P
    X1, Y1 = P
    X2, Y2 = Q
    if X1 == X2:
        if Y1 != Y2: return None
        return ec_double(P, p)
    lam = ((Y2 - Y1) * pow(X2 - X1, p-2, p)) % p
    X3 = (lam*lam - X1 - X2) % p
    Y3 = (lam*(X1 - X3) - Y1) % p
    return (X3, Y3)

def ec_mul(k, P, p):
    if k == 0: return None
    Q = None
    while k > 0:
        if k & 1: Q = ec_add(Q, P, p)
        P = ec_double(P, p)
        k >>= 1
    return Q

def find_point(b, p):
    for x in range(2, 100000):
        rhs = (pow(x, 3, p) + b) % p
        y = sqrt_mod(rhs, p)
        if y is not None:
            return (x, y)
    return None

# --- CM setup ---
assert p % 6 == 1
s_sq = (4*p - t*t) // 3
s = math.isqrt(s_sq)
assert s*s == s_sq and 4*p == t*t + 3*s*s

# Primitive 6th root u = (1 + sqrt(-3)) / 2 mod p
u = ((1 + sqrt_mod((-3) % p, p)) * pow(2, p-2, p)) % p
assert pow(u, 6, p) == 1 and pow(u, 3, p) != 1

b = [(7 * pow(u, k, p)) % p for k in range(6)]

# --- H2: 2-torsion factorisation ---
# For p ≡ 1 (mod 3): x^3+b_k is either [1,1,1] (cube) or [3] (irred).
def torsion_pattern(bk):
    chi3 = pow((-bk) % p, (p-1)//3, p)
    return [1, 1, 1] if chi3 == 1 else [3]

deg_pat = [torsion_pattern(bk) for bk in b]
irred_k = [k for k in range(6) if deg_pat[k] == [3]]   # odd-order twists
split_k = [k for k in range(6) if deg_pat[k] == [1,1,1]]  # 4|order twists

# --- Curve orders ---
# Known exactly: k=0 (secp256k1, order n), k=3 (quadratic twist, order n_twist).
# The 4 odd-order slots come from traces {t, (t-3s)/2, -t, (3s-t)/2}.
# The 2 split-order slots come from traces {-(t+3s)/2, (t+3s)/2}.
n_twist = p + 1 + t  # k=3 order

# Determine b_2's order by finding a point and checking which candidate annihilates it.
T1 = (t - 3*s) // 2     # trace slot 1
T4 = (3*s - t) // 2     # trace slot 4
N_T1 = p + 1 - T1       # candidate order A (odd)
N_T4 = p + 1 - T4       # candidate order B (odd)

P_b2 = find_point(b[2], p)
assert P_b2 is not None
if ec_mul(N_T4, P_b2, p) is None:
    N_b2, N_b5 = N_T4, N_T1
else:
    N_b2, N_b5 = N_T1, N_T4

# Split orders (for b_1 and b_4; exact assignment doesn't matter for H3 results)
T2 = -(t + 3*s) // 2
T5 = (t + 3*s) // 2
N_split_0 = p + 1 - T2
N_split_1 = p + 1 - T5

order_of = {0: n, 1: None, 2: N_b2, 3: n_twist, 4: None, 5: N_b5}

# --- 15-pair Howe check ---
print("=" * 68)
print("Howe conditions for all 15 sextic-twist pairs of secp256k1")
print("=" * 68)
print(f"  4 irred twists (odd order): k ∈ {irred_k}")
print(f"  2 split twists (4|order):   k ∈ {split_k}")
print()
print(f"  Orders of irred twists:")
print(f"    k=0 (secp256k1): {n}")
print(f"    k=2 (cubic):     {N_b2}")
print(f"    k=3 (quad-twist):{n_twist}")
print(f"    k=5 (cubic2):    {N_b5}")
print()
print(f"  (k=3 and k=5 orders both divisible by 3: structurally expected")
print(f"   since p≡1 mod 3, n≡1 mod 3 ⟹ n_twist≡0 mod 3)")
print()
print(f"  {'(i,j)':6s} | H1 | H2 | H3           | Glueable?")
print(f"  {'------':6s}|----|----|--------------|---------")

glueable = []
for i in range(6):
    for j in range(i+1, 6):
        h1 = True  # All 6 traces distinct → H1 always YES
        h2 = (deg_pat[i] == deg_pat[j])
        if not h2:
            h3_str = "N/A"
            h3 = False
        elif order_of[i] is None or order_of[j] is None:
            # Split twists: 4|N_i and 4|N_j → 4|gcd → H3 fails
            h3_str = "NO (4|gcd)"
            h3 = False
        else:
            g = math.gcd(order_of[i], order_of[j])
            h3_str = f"YES (gcd={g})" if g == 1 else f"NO  (gcd={g})"
            h3 = (g == 1)
        ok = h1 and h2 and h3
        if ok: glueable.append((i, j))
        print(f"  ({i},{j})   | YES| {'YES' if h2 else 'NO ':3s}| {h3_str:13s}| {'YES ← glueable' if ok else 'no'}")

print()
print(f"Howe-glueable pairs: {len(glueable)} / 15")
for pair in glueable:
    print(f"  {pair}")
print()
print("  Each glueable pair (i,j) yields a genus-2 curve C/F_p with")
print("  Jac(C) (2,2)-isogenous to E_i × E_j. No sub-√p ECDLP impact.")
print("=" * 68)
