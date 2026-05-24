"""
Howe (2,2)-conditions for all 15 pairs of j=0 sextic twists of secp256k1.

For each pair (E_i, E_j) of the 6 sextic twists of y²=x³+7 over F_p,
check:
  H1: #E_i ≠ #E_j          (not F_p-isogenous)
  H2: x³+b_i, x³+b_j same  (same F_p-Galois structure of E[2])
       factorisation pattern
  H3: gcd(#E_i, #E_j) = 1  (no shared F_p-order obstruction)

Strategy:
  - CM formula 4p = t² + 3s² gives all 6 traces ±t, ±(t±3s)/2.
  - Trace matching: for each twist b_k, find a random affine point P,
    test which candidate trace T satisfies (p+1-T)*P = O.
  - H2: factor x³+b_k over F_p by testing for roots (x³+b_k ≡ 0 mod p)
    and reducibility.
"""

import math
import sys
from itertools import combinations

# secp256k1 parameters
p = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F
n_known = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
t_known = p + 1 - n_known

# ---------------------------------------------------------------------------
# Modular arithmetic helpers
# ---------------------------------------------------------------------------

def legendre(a, p):
    return pow(a % p, (p - 1) // 2, p)

def modsqrt(a, p):
    """Tonelli-Shanks square root mod p. Returns r with r²≡a (mod p), or None."""
    a %= p
    if a == 0:
        return 0
    if legendre(a, p) != 1:
        return None
    if p % 4 == 3:
        return pow(a, (p + 1) // 4, p)
    # General Tonelli-Shanks
    s, q = 0, p - 1
    while q % 2 == 0:
        s += 1
        q //= 2
    z = 2
    while legendre(z, p) != p - 1:
        z += 1
    m_ts, c, t_ts, r = s, pow(z, q, p), pow(a, q, p), pow(a, (q + 1) // 2, p)
    while True:
        if t_ts == 1:
            return r
        i, temp = 1, (t_ts * t_ts) % p
        while temp != 1:
            temp = (temp * temp) % p
            i += 1
        b = pow(c, 1 << (m_ts - i - 1), p)
        m_ts, c, t_ts, r = i, (b * b) % p, (t_ts * b * b) % p, (r * b) % p

# ---------------------------------------------------------------------------
# Weierstrass short-form EC arithmetic over F_p: y² = x³ + B
# Points as (x, y) tuples; None = point at infinity.
# ---------------------------------------------------------------------------

def ec_add(P, Q, B, p):
    if P is None:
        return Q
    if Q is None:
        return P
    px, py = P
    qx, qy = Q
    if px == qx:
        if (py + qy) % p == 0:
            return None  # P + (-P) = O
        # P == Q: doubling
        lam = (3 * px * px) * pow(2 * py, p - 2, p) % p
    else:
        lam = (qy - py) * pow(qx - px, p - 2, p) % p
    rx = (lam * lam - px - qx) % p
    ry = (lam * (px - rx) - py) % p
    return (rx, ry)

def ec_mul(k, P, B, p):
    """Scalar multiplication k*P on y²=x³+B mod p (no reduction — exact)."""
    if k == 0 or P is None:
        return None
    result = None
    addend = P
    while k:
        if k & 1:
            result = ec_add(result, addend, B, p)
        addend = ec_add(addend, addend, B, p)
        k >>= 1
    return result

def find_point(B, p):
    """Find the first affine point on y²=x³+B mod p by scanning x=1,2,..."""
    for x in range(1, 10000):
        ysq = (x * x * x + B) % p
        y = modsqrt(ysq, p)
        if y is not None and y != 0:
            return (x, y)
    return None

# ---------------------------------------------------------------------------
# 2-torsion factorisation pattern of x³+b mod p
# ---------------------------------------------------------------------------

def torsion_pattern(b, p):
    """
    Degree pattern of x³+b over F_p.
    Returns one of: [1,1,1], [1,2], [3].
    """
    # Count roots of x³+b ≡ 0 mod p, i.e., x³ ≡ -b mod p
    roots = [x for x in range(p) if pow(x, 3, p) == (-b) % p] if p < 1000 else None
    if roots is not None:
        n_roots = len(roots)
        if n_roots == 3:
            return [1, 1, 1]
        elif n_roots == 1:
            return [1, 2]
        else:
            return [3]
    # For large p: x³ ≡ -b mod p has a solution iff (-b)^{(p-1)/3} ≡ 1 mod p
    # (since gcd(3, p-1) = 3 when p ≡ 1 mod 3, which holds for secp256k1)
    neg_b = (-b) % p
    if neg_b == 0:
        return [1, 1, 1]  # b=0 → x³ = 0 triple root (degenerate, shouldn't happen)
    cube_test = pow(neg_b, (p - 1) // 3, p)
    if cube_test == 1:
        # -b is a cube → 3 roots (but are they in F_p? yes, if -b has a cube root in F_p)
        # Actually: x³ = -b has exactly 3 solutions in F_p when -b is a cube and p ≡ 1 mod 3
        # (since (p-1)/3 divides into 3 cube roots of unity times one cube root of -b)
        # BUT this requires the cube root of -b to be in F_p.  If -b is a cube in F_p*,
        # then it has exactly 3 cube roots in F_p*.
        return [1, 1, 1]
    else:
        # -b is not a cube in F_p → x³+b is irreducible over F_p iff it has no root.
        # (A cubic is irreducible iff it has no root.)
        return [3]

# ---------------------------------------------------------------------------
# Main computation
# ---------------------------------------------------------------------------

print("=" * 72)
print("Howe conditions for all 15 pairs of secp256k1 sextic twists")
print("=" * 72)
print()
print(f"p   = {hex(p)}")
print(f"t   = {t_known}")
print()

# Step 1: p ≡ 1 (mod 6)
assert p % 6 == 1, f"p ≡ {p%6} (mod 6), expected 1"
print(f"p mod 6 = {p % 6}  → exactly 6 sextic twist classes ✓")
print()

# Step 2: CM structure 4p = t² + 3s²
rem = 4 * p - t_known ** 2
assert rem % 3 == 0, "4p - t² not divisible by 3"
s_sq = rem // 3
s = math.isqrt(s_sq)
assert s * s == s_sq, "s is not an integer; CM formula broken"
print(f"4p = t² + 3s²  with  s = {s}")
print(f"Verify: {4*p == t_known**2 + 3*s**2}")
print()

# Step 3: Six traces
t_plus  = (t_known + 3 * s) // 2
t_minus = (t_known - 3 * s) // 2
assert (t_known + 3 * s) % 2 == 0, "(t+3s) must be even"
# Assign in ζ₆^k order (from CM unit analysis):
# k=0: t, k=1: (t-3s)/2, k=2: -(t+3s)/2, k=3: -t, k=4: (3s-t)/2, k=5: (t+3s)/2
T = [t_known, t_minus, -t_plus, -t_known, s - t_minus, t_plus]
# (Note: (3s-t)/2 = s - (t-3s)/2 = s - t_minus ... let me recompute:
#  (3s-t)/2 = -(t-3s)/2 = -t_minus )
T = [t_known, t_minus, -t_plus, -t_known, -t_minus, t_plus]
N = [p + 1 - tk for tk in T]

print("Six CM traces and orders:")
print(f"  k | {'trace':>66} | order (p+1-trace)")
print(f"  --|{'':->68}|-------------------")
for k in range(6):
    print(f"  {k} | {T[k]:>66} | {N[k]}")
print()
print(f"N[0] = n_known? {N[0] == n_known}")
print(f"N[3] = p+1+t (quad twist)? {N[3] == p + 1 + t_known}")
print()

# Step 4: Find primitive 6th root u mod p
# u = smallest integer with ord(u) = 6 mod p
# u^6 ≡ 1, u^3 ≡ -1 (order-2 subgroup), u^2 ≢ 1
print("Finding primitive 6th root u mod p...")
u = None
for cand in range(2, 10000):
    if pow(cand, (p - 1) // 6, p) != 1:
        continue
    if pow(cand, (p - 1) // 2, p) != 1:
        continue  # must be a square (6th power implies square)
    # Check order is exactly 6: u^6≡1, u^3≡-1 (≢1), u^2≢1
    u6 = pow(cand, 6 * ((p - 1) // 6), p)
    u3 = pow(cand, 3 * ((p - 1) // 6), p)
    u2 = pow(cand, 2 * ((p - 1) // 6), p)
    if u6 == 1 and u3 == p - 1 and u2 != 1:
        u = pow(cand, (p - 1) // 6, p)
        print(f"  primitive root candidate: g={cand}")
        break

if u is None:
    # Fallback: use any generator of (F_p*)/(F_p*)^6
    # u = g^((p-1)/6) for some g that is a primitive root mod p
    # Since we can't easily find a primitive root for 256-bit p, use a known trick:
    # Try small values as generators (primitive roots mod p)
    for g_cand in range(2, 1000):
        # Check that g_cand is a primitive root: g^((p-1)/q) ≢ 1 for all prime q | p-1
        # For secp256k1, p-1 = 2 * (p-1)/2 where (p-1)/2 is the cofactor
        # Quick check: just verify order is not a proper divisor of 6
        gu = pow(g_cand, (p - 1) // 6, p)
        if pow(gu, 3, p) == p - 1 and pow(gu, 2, p) != 1:
            u = gu
            print(f"  u = {g_cand}^{{(p-1)/6}} mod p")
            break

print(f"u = {u}")
print(f"u^6 ≡ 1 (mod p): {pow(u, 6, p) == 1}")
print(f"u^3 ≡ -1 (mod p): {pow(u, 3, p) == p - 1}")
print(f"u^2 ≢ 1 (mod p): {pow(u, 2, p) != 1}")
print()

# Step 5: b-values
b_vals = [(7 * pow(u, k, p)) % p for k in range(6)]
print("Six b-values (b_k = 7 * u^k mod p):")
for k in range(6):
    print(f"  k={k}: b = {b_vals[k]}")
print(f"b_3 = p-7? {b_vals[3] == p - 7}")
print()

# Step 6: Match b_k to trace via point-order test
print("Matching b_k to trace (scalar-multiplication test)...")
print()
trace_of_k = [None] * 6
trace_of_k[0] = T[0]  # known: secp256k1
trace_of_k[3] = T[3]  # known: quadratic twist (b_3 = -7 mod p)

for k in range(6):
    if trace_of_k[k] is not None:
        print(f"  k={k}: trace = {trace_of_k[k]}  (known a priori)")
        continue
    B = b_vals[k]
    P = find_point(B, p)
    if P is None:
        print(f"  k={k}: WARNING: no affine point found in x=1..10000")
        continue
    px, py = P
    print(f"  k={k}: testing point ({px}, ...)")
    found = False
    for ti, T_cand in enumerate(T):
        n_cand = p + 1 - T_cand
        result = ec_mul(n_cand, P, B, p)
        if result is None:
            trace_of_k[k] = T_cand
            print(f"  k={k}: trace matched T[{ti}] = {T_cand}  (order {n_cand})")
            found = True
            break
    if not found:
        print(f"  k={k}: WARNING: no trace match found!")

print()

# Step 7: H2 — 2-torsion factorisation patterns
print("2-torsion factorisation patterns (degree pattern of x³+b_k over F_p):")
print()
deg_pats = []
for k in range(6):
    pat = torsion_pattern(b_vals[k], p)
    deg_pats.append(pat)
    print(f"  k={k}: x³ + {b_vals[k]} → pattern {pat}")

print()

# Step 8: All 15 pairwise checks
print("=" * 72)
print("All 15 pairwise Howe checks")
print("=" * 72)
print()
print(f"  {'(i,j)':6} | {'H1: n_i≠n_j':12} | {'H2: same 2-tor':15} | {'H3: gcd=1':10} | Glueable?")
print(f"  {'------':6}+{'------------':13}+{'---------------':16}+{'----------':11}+----------")

glueable = []
for i, j in combinations(range(6), 2):
    if trace_of_k[i] is None or trace_of_k[j] is None:
        print(f"  ({i},{j})    | SKIP (no trace match)")
        continue
    ni = p + 1 - trace_of_k[i]
    nj = p + 1 - trace_of_k[j]
    h1 = (ni != nj)
    h2 = (deg_pats[i] == deg_pats[j])
    h3 = (math.gcd(ni, nj) == 1)
    ok = h1 and h2 and h3
    h1s = "✓" if h1 else "✗"
    h2s = "✓" if h2 else "✗"
    h3s = "✓" if h3 else "✗"
    verdict = "YES ← Howe-glueable" if ok else "no"
    print(f"  ({i},{j})    | {h1s:12} | {h2s:15} | {h3s:10} | {verdict}")
    if ok:
        glueable.append((i, j))

print()
print(f"Howe-glueable pairs: {len(glueable)} / 15")
if glueable:
    print(f"Glueable pairs: {glueable}")
print()

# Summary table: gcd values for all pairs
print("gcd(n_i, n_j) table:")
print(f"  {'':4}", end="")
for j in range(6):
    print(f"  k={j}", end="")
print()
for i in range(6):
    print(f"  k={i}", end="")
    for j in range(6):
        if trace_of_k[i] is None or trace_of_k[j] is None:
            print(f"   {'?':4}", end="")
        elif i == j:
            print(f"   {'—':4}", end="")
        else:
            ni = p + 1 - trace_of_k[i]
            nj = p + 1 - trace_of_k[j]
            g = math.gcd(ni, nj)
            if g == 1:
                print(f"   {'1':4}", end="")
            else:
                # Print small factor
                f = 2
                while f * f <= g:
                    if g % f == 0:
                        break
                    f += 1
                print(f"   {f if f*f<=g else g}…", end="")
    print()
print()

print("=" * 72)
print("Done.")
print("=" * 72)
