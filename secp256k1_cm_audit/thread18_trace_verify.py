"""
Thread 18: Verify correct CM trace assignment for secp256k1 sextic twists.
Uses Python's pow() and a simple point addition to confirm which trace
belongs to each of the 6 sextic twist curves.
"""

p = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F
n_secp = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
t_secp = p + 1 - n_secp
# 4p = t^2 + 3s^2
s_sq = (4*p - t_secp**2) // 3
s = int(s_sq**0.5)
# verify with a few Newton iterations
while (s+1)**2 <= s_sq:
    s += 1
while s**2 > s_sq:
    s -= 1
assert s*s == s_sq, "s not integer"
assert 4*p == t_secp**2 + 3*s**2, "CM formula check failed"

# Primitive 6th root of unity in F_p*
# g = 3 is a primitive root (same as PARI found)
g = 3
u = pow(g, (p-1)//6, p)
assert pow(u, 6, p) == 1
assert pow(u, 3, p) == p-1  # u^3 = -1
assert pow(u, 2, p) != 1

# Six b-values
b = [7 * pow(u, k, p) % p for k in range(6)]
assert b[3] == p - 7  # quadratic twist

# Six abstract CM traces
T1 = t_secp
T2 = (t_secp - 3*s) // 2
T3 = -(t_secp + 3*s) // 2
T4 = -t_secp
T5 = (3*s - t_secp) // 2
T6 = (t_secp + 3*s) // 2
T_all = [T1, T2, T3, T4, T5, T6]

# Modular square root via Tonelli-Shanks
def sqrt_mod(n, p):
    n = n % p
    if n == 0:
        return 0
    assert pow(n, (p-1)//2, p) == 1, "not a QR"
    if p % 4 == 3:
        return pow(n, (p+1)//4, p)
    # Tonelli-Shanks
    q = p - 1
    r = 0
    while q % 2 == 0:
        q //= 2
        r += 1
    z = 2
    while pow(z, (p-1)//2, p) != p-1:
        z += 1
    m = r
    c = pow(z, q, p)
    t2 = pow(n, q, p)
    rr = pow(n, (q+1)//2, p)
    while True:
        if t2 == 1:
            return rr
        i = 1
        tt = (t2 * t2) % p
        while tt != 1:
            tt = (tt * tt) % p
            i += 1
        b2 = pow(c, 1 << (m-i-1), p)
        m = i
        c = (b2 * b2) % p
        t2 = (t2 * c) % p
        rr = (rr * b2) % p

# Elliptic curve point addition for y^2 = x^3 + b (a=0)
def point_add(P1, P2, p):
    if P1 is None:
        return P2
    if P2 is None:
        return P1
    x1, y1 = P1
    x2, y2 = P2
    if x1 == x2:
        if y1 != y2 or y1 == 0:
            return None  # point at infinity
        # doubling
        lam = (3 * x1 * x1) * pow(2 * y1, p-2, p) % p
    else:
        lam = (y2 - y1) * pow(x2 - x1, p-2, p) % p
    x3 = (lam*lam - x1 - x2) % p
    y3 = (lam*(x1 - x3) - y1) % p
    return (x3, y3)

def scalar_mult(n, P, p):
    if n == 0:
        return None
    if n < 0:
        x, y = P
        P = (x, (-y) % p)
        n = -n
    result = None
    addend = P
    while n:
        if n & 1:
            result = point_add(result, addend, p)
        addend = point_add(addend, addend, p)
        n >>= 1
    return result

def find_point(bk, p):
    """Find an F_p-rational point on y^2 = x^3 + bk."""
    for x in range(2, 10000):
        rhs = (pow(x, 3, p) + bk) % p
        if pow(rhs, (p-1)//2, p) == 1:
            y = sqrt_mod(rhs, p)
            assert (y*y - pow(x,3,p) - bk) % p == 0
            return (x, y)
    raise ValueError("No point found")

def test_trace(bk, T_cand, p):
    """Test whether T_cand is the correct trace for E: y^2 = x^3 + bk."""
    N_cand = p + 1 - T_cand
    P = find_point(bk, p)
    result = scalar_mult(N_cand, P, p)
    return result is None  # True if N_cand * P = O

print("=== Thread 18: CM trace assignment verification ===")
print()
print(f"t_secp = {t_secp}")
print(f"s = {s}")
print()
print("Abstract traces and N mod 4:")
for i, Ti in enumerate(T_all, 1):
    Ni = p + 1 - Ti
    print(f"  T{i} = {Ti:+d}  N{i} mod 4 = {Ni % 4}")
print()
print("k=1 and k=4 have [1,1,1] 2-torsion -> need 4|N.")
print("Traces with 4|N: T3 and T6.")
print()

# Test k=1
print("=== k=1: Testing T3 vs T6 ===")
for Ti_cand, label in [(T3, "T3"), (T6, "T6")]:
    correct = test_trace(b[1], Ti_cand, p)
    Ni = p + 1 - Ti_cand
    print(f"  {label} = {Ti_cand:+d}  N = {Ni % (10**20)}...  N mod 4 = {Ni % 4}  ->  N*P = O? {correct}")

print()

# Test k=2
print("=== k=2: Testing T2 vs T5 ===")
for Ti_cand, label in [(T2, "T2"), (T5, "T5")]:
    correct = test_trace(b[2], Ti_cand, p)
    Ni = p + 1 - Ti_cand
    print(f"  {label} = {Ti_cand:+d}  N = ...{Ni % (10**20)}  N mod 4 = {Ni % 4}  ->  N*P = O? {correct}")

print()

# Determine correct assignment
print("=== Determining correct assignment ===")
if test_trace(b[1], T3, p):
    T_k1, T_k4 = T3, T6
    print("k=1: trace = T3 = -(t+3s)/2")
    print("k=4: trace = T6 =  (t+3s)/2")
else:
    T_k1, T_k4 = T6, T3
    print("k=1: trace = T6 =  (t+3s)/2")
    print("k=4: trace = T3 = -(t+3s)/2")

if test_trace(b[2], T2, p):
    T_k2, T_k5 = T2, T5
    print("k=2: trace = T2 = (t-3s)/2")
    print("k=5: trace = T5 = (3s-t)/2")
else:
    T_k2, T_k5 = T5, T2
    print("k=2: trace = T5 = (3s-t)/2")
    print("k=5: trace = T2 = (t-3s)/2")

print()

T_correct = [T1, T_k1, T_k2, T4, T_k4, T_k5]
N_correct = [p + 1 - T for T in T_correct]

print("=== Correct orders ===")
for k in range(6):
    print(f"  k={k}: trace = {T_correct[k]:+d}  N = {N_correct[k]}  N mod 4 = {N_correct[k] % 4}")
print()

# 2-torsion patterns (from PARI factorization result)
# k=0: [3]  k=1: [1,1,1]  k=2: [3]  k=3: [3]  k=4: [1,1,1]  k=5: [3]
deg_pat = [[3], [1,1,1], [3], [3], [1,1,1], [3]]

print("=== Recomputing all 15 Howe conditions with CORRECT orders ===")
print()
print(f"  (i,j) | H1: n_i≠n_j | H2: same 2-tor | H3: gcd=1 | Glueable?")
print(f"  ------|-------------|----------------|-----------|----------")

glueable = []
for i in range(6):
    for j in range(i+1, 6):
        ni, nj = N_correct[i], N_correct[j]
        h1 = (ni != nj)
        h2 = (deg_pat[i] == deg_pat[j])
        h3 = (__import__('math').gcd(ni, nj) == 1)
        g = h1 and h2 and h3
        label = "YES <- Howe-glueable" if g else "no"
        print(f"  ({i},{j})  | {'YES' if h1 else 'NO '}  | {'YES' if h2 else 'NO '}  | {'YES' if h3 else 'NO '}  | {label}")
        if g:
            glueable.append((i,j))

print()
print(f"Howe-glueable pairs (CORRECT assignment): {len(glueable)} / 15")
print(f"Glueable pairs: {glueable}")

# Compare to script's answer
script_glueable = [(0,2),(0,3),(0,5),(1,4),(2,3)]
print()
print(f"Script's answer (assumed T[k+1]<->b_k): {script_glueable}")
print(f"Correct answer:                          {glueable}")
print(f"Answers match: {sorted(glueable) == sorted(script_glueable)}")
