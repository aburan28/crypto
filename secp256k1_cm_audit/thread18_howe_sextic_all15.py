"""
Thread 18: Howe (2,2)-conditions for all 15 pairs of j=0 sextic twists.
secp256k1: y²=x³+7 over F_p (j=0, CM disc=-3).

For p ≡ 1 (mod 6) there are exactly 6 sextic-twist iso-classes of j=0
curves over F_p.  They are y²=x³+b_k where b_k = 7*u^k mod p and
u = primitive 6th root of unity in F_p*.

For each of C(6,2)=15 pairs (i,j) we check Howe (1996) conditions:
  H1: n_i ≠ n_j          (E_i ≇ E_j over F_p)
  H2: x³+b_i and x³+b_j have the same F_p-factorisation pattern
      (same Galois module structure of E[2])
  H3: gcd(n_i, n_j) = 1  (no common F_p-rational order obstruction)

Trace assignment: CM formula 4p = t²+3s² gives all 6 traces directly:
  T_0 = t,  T_1 = (t-3s)/2,  T_2 = -(t+3s)/2,
  T_3 = -t, T_4 = -(t-3s)/2, T_5 = (t+3s)/2.
No scalar multiplication needed.

2-torsion pattern for p≡1(mod 3): x³+b is [1,1,1] iff -b is a cube
in F_p*, otherwise [3].  The intermediate case [1,2] is impossible.
"""

import math
from itertools import combinations

p = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F
n_known = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
t_known = p + 1 - n_known

print("=" * 72)
print("Thread 18: Howe conditions for all 15 pairs of secp256k1 sextic twists")
print("=" * 72)
print()
print(f"p       = {hex(p)}")
print(f"t_known = {t_known}")
print()

# --- Step 1: p ≡ 1 (mod 6) ---
assert p % 6 == 1, f"Expected p≡1(mod 6); got p≡{p%6}(mod 6)"
print(f"Step 1: p mod 6 = {p % 6}  → exactly 6 sextic twist classes ✓")
print()

# --- Step 2: CM structure 4p = t² + 3s² ---
rem = 4 * p - t_known ** 2
assert rem % 3 == 0, "4p - t^2 not divisible by 3"
s_sq = rem // 3
s = math.isqrt(s_sq)
assert s * s == s_sq, "s is not a perfect square"
assert 4 * p == t_known ** 2 + 3 * s ** 2, "CM formula check failed"
print(f"Step 2: s = {s}  (from 4p = t² + 3s²)")
print(f"        4p = t² + 3s²: {4*p == t_known**2 + 3*s**2} ✓")
print()

# --- Step 3: Six traces (CM formula) ---
# The 6 sextic twists of y²=x³+b over F_p for p≡1(mod 6) have traces:
#   k=0: t_known       k=1: (t-3s)/2     k=2: -(t+3s)/2
#   k=3: -t_known      k=4: -(t-3s)/2    k=5: (t+3s)/2
t_p = (t_known + 3 * s)
t_m = (t_known - 3 * s)
assert t_p % 2 == 0 and t_m % 2 == 0, "(t±3s) must both be even"
T = [
    t_known,           # k=0
    t_m // 2,          # k=1  (t-3s)/2
    -t_p // 2,         # k=2  -(t+3s)/2
    -t_known,          # k=3
    -t_m // 2,         # k=4  -(t-3s)/2
    t_p // 2,          # k=5  (t+3s)/2
]
N = [p + 1 - tk for tk in T]

print("Step 3: Six CM traces and curve orders")
print(f"  {'k':2} | {'trace':>66} | order")
print(f"  ---+{'-'*68}+------")
for k in range(6):
    print(f"  {k:2} | {T[k]:>66} | {N[k]}")
print()
assert N[0] == n_known, "N[0] ≠ n_known"
assert N[3] == p + 1 + t_known, "N[3] ≠ p+1+t (quad twist)"
print(f"  Sanity N[0]=n_known: ✓   N[3]=p+1+t: ✓")
print()

# --- Step 4: Primitive 6th root u mod p ---
# Find smallest integer g such that u = g^((p-1)/6) has order 6 in F_p*.
# u has order 6 iff u^6=1, u^3=-1 (=-1 mod p), u^2≠1.
print("Step 4: Finding primitive 6th root u mod p...")
u = None
for g_cand in range(2, 100000):
    cand_u = pow(g_cand, (p - 1) // 6, p)
    if cand_u == 1:
        continue
    if pow(cand_u, 3, p) == p - 1 and pow(cand_u, 2, p) != 1:
        u = cand_u
        print(f"  g = {g_cand},  u = g^((p-1)/6) mod p")
        break
assert u is not None, "Failed to find primitive 6th root"
assert pow(u, 6, p) == 1,     "u^6 ≢ 1"
assert pow(u, 3, p) == p - 1, "u^3 ≢ -1"
assert pow(u, 2, p) != 1,     "u^2 = 1 (order < 6)"
print(f"  u^6 ≡ 1 (mod p): ✓   u^3 ≡ -1: ✓   u^2 ≢ 1: ✓")
print()

# --- Step 5: Six b-values ---
b = [(7 * pow(u, k, p)) % p for k in range(6)]
print("Step 5: b-values (b_k = 7·u^k mod p):")
for k in range(6):
    print(f"  k={k}: b_{k} = {b[k]}")
assert b[3] == p - 7, f"b_3 should equal p-7; got {b[3]}"
print(f"  b_3 = p-7: ✓")
print()

# --- Step 6: 2-torsion factorisation patterns ---
# For p≡1(mod 3): x³+b has pattern [1,1,1] iff -b is a cube in F_p*,
# i.e., (-b)^((p-1)/3) ≡ 1 (mod p).  Otherwise: [3] (irreducible).
# The intermediate [1,2] pattern is impossible for p≡1(mod 3).
exp3 = (p - 1) // 3
def torsion_pattern(bk):
    neg_b = (-bk) % p
    if neg_b == 0:
        return [1, 1, 1]  # degenerate (b=0), shouldn't happen
    cube_test = pow(neg_b, exp3, p)
    return [1, 1, 1] if cube_test == 1 else [3]

deg_pats = [torsion_pattern(bk) for bk in b]
print("Step 6: 2-torsion factorisation patterns (degree of x³+b_k over F_p):")
for k in range(6):
    cube_flag = "cube" if deg_pats[k] == [1, 1, 1] else "non-cube"
    print(f"  k={k}: pattern={deg_pats[k]}  (-b_{k} is {cube_flag} mod p)")
print()

# Note: pattern [3] means x³+b is irreducible → E has no F_p-rational 2-torsion
# → E has prime order (necessary for cryptographic use).
# For secp256k1 (k=0), pattern should be [3] (since n is prime → no 2-torsion).
if deg_pats[0] == [3]:
    print("  k=0 (secp256k1): pattern [3] → no F_p-rational 2-torsion ✓ (prime order)")
else:
    print(f"  WARNING: k=0 pattern = {deg_pats[0]}  (unexpected for prime-order curve)")
print()

# --- Step 7: All 15 pairwise Howe checks ---
print("=" * 72)
print("Step 7: All 15 pairwise Howe checks")
print("=" * 72)
print()
print(f"  {'(i,j)':6} | H1(n_i≠n_j) | H2(same 2-tor) | H3(gcd=1) | gcd(n_i,n_j) | Glueable?")
print(f"  ------+-------------+----------------+-----------+--------------+-----------")

glueable_pairs = []
gcd_table = {}

for i, j in combinations(range(6), 2):
    ni, nj = N[i], N[j]
    h1 = (ni != nj)
    h2 = (deg_pats[i] == deg_pats[j])
    g = math.gcd(ni, nj)
    h3 = (g == 1)
    gcd_table[(i, j)] = g
    ok = h1 and h2 and h3
    marks = (
        ("✓" if h1 else "✗"),
        ("✓" if h2 else "✗"),
        ("✓" if h3 else "✗"),
    )
    # Show small prime factor of gcd if gcd > 1
    if g == 1:
        gcd_str = "1"
    else:
        f = 2
        while f * f <= g and g % f != 0:
            f += 1
        if f * f <= g:
            gcd_str = f"{f}…"
        else:
            gcd_str = str(g)

    verdict = "YES ← glueable" if ok else "no"
    print(
        f"  ({i},{j})    | {marks[0]:11} | {marks[1]:14} | {marks[2]:9} | "
        f"{gcd_str:12} | {verdict}"
    )
    if ok:
        glueable_pairs.append((i, j))

print()
print(f"Howe-glueable pairs: {len(glueable_pairs)} / 15")
if glueable_pairs:
    print(f"Glueable pair indices (0-based k): {glueable_pairs}")
print()

# --- Step 8: Summary of H1/H2/H3 failures ---
print("=" * 72)
print("Step 8: Analysis of condition failures")
print("=" * 72)
print()

# H1 failures: pairs with same curve order
h1_failures = [(i,j) for i,j in combinations(range(6),2) if N[i] == N[j]]
print(f"H1 failures (n_i = n_j): {len(h1_failures)} pairs")
for i,j in h1_failures:
    print(f"  ({i},{j}): n_{i} = n_{j} = {N[i]}")
print()

# H2 failures: pairs with different 2-torsion pattern
h2_failures = [(i,j) for i,j in combinations(range(6),2)
               if deg_pats[i] != deg_pats[j]]
print(f"H2 failures (different 2-torsion pattern): {len(h2_failures)} pairs")
for i,j in h2_failures:
    print(f"  ({i},{j}): k={i} → {deg_pats[i]},  k={j} → {deg_pats[j]}")
print()

# H3 failures: pairs with gcd > 1
h3_failures = [(i,j) for i,j in combinations(range(6),2)
               if gcd_table[(i,j)] > 1]
print(f"H3 failures (gcd > 1): {len(h3_failures)} pairs")
for i,j in h3_failures:
    g = gcd_table[(i,j)]
    print(f"  ({i},{j}): gcd(n_{i}, n_{j}) = {g}")
    # Factorize g (small factor check)
    facs = []
    gg = g
    for q in range(2, min(10000, gg+1)):
        if gg % q == 0:
            e = 0
            while gg % q == 0:
                gg //= q
                e += 1
            facs.append((q, e))
        if gg == 1:
            break
    if gg > 1:
        facs.append((gg, 1))
    print(f"         small prime factors: {facs[:5]}{'...' if len(facs)>5 else ''}")
print()

# --- Step 9: Structural interpretation ---
print("=" * 72)
print("Step 9: Structural interpretation")
print("=" * 72)
print()

all_pats = set(str(p) for p in deg_pats)
print(f"Distinct 2-torsion patterns across 6 twists: {sorted(all_pats)}")
print()
pat_groups = {}
for k in range(6):
    key = str(deg_pats[k])
    pat_groups.setdefault(key, []).append(k)
for pat, ks in sorted(pat_groups.items()):
    print(f"  Pattern {pat}: twists k ∈ {ks}")
print()

print("Key observations:")
print(f"  1. Number of distinct 2-torsion patterns: {len(pat_groups)}")
print(f"     (H2 holds within same-pattern group; fails across groups)")
print()
print(f"  2. 6 traces: {T}")
print(f"     All distinct (H1): {len(set(T)) == 6}")
print()

# Check cube residue pattern among b_k
cubes = [(k, deg_pats[k] == [1,1,1]) for k in range(6)]
cube_ks = [k for k, is_cube in cubes if is_cube]
noncube_ks = [k for k, is_cube in cubes if not is_cube]
print(f"  3. Cubic residue of -b_k mod p:")
print(f"     Cube    (-b_k is a cube): k ∈ {cube_ks}   → x³+b_k fully splits [1,1,1]")
print(f"     Non-cube: k ∈ {noncube_ks} → x³+b_k irreducible [3]")
print()
print(f"     Note: b_k = 7·u^k, so -b_k = -7·u^k.")
print(f"     (-7·u^k) is a cube mod p iff (-7)·(u^k) is a cube,")
print(f"     iff k ≡ k_0 (mod 3) for the unique k_0 where -7·u^k₀ is a cube.")
print()

# Verify the cube / non-cube split has the expected Z/3Z structure
# Among k=0..5, exactly 2 values of k mod 3 partition the cubes:
# if -7 is not a cube, then -7*u^k is a cube iff u^k is (1/-7-th of a cube),
# i.e., iff k ≡ k₀ (mod 3) for some k₀.
# So the pattern of cubes among k=0..5 follows a 2-vs-4 or 0-vs-6 split.
print(f"  4. Cube/non-cube split: {len(cube_ks)} cubes, {len(noncube_ks)} non-cubes")
if len(cube_ks) == 2:
    print(f"     → 2 twists with [1,1,1] (E[2] all F_p-rational, E has even order)")
    print(f"     → 4 twists with [3] (E[2] trivial, E has odd (possibly prime) order)")
elif len(cube_ks) == 0:
    print(f"     → All 6 twists have [3] (no F_p-rational 2-torsion)")
elif len(cube_ks) == 6:
    print(f"     → All 6 twists have [1,1,1] (all have F_p-rational 2-torsion)")
else:
    print(f"     → Unexpected split {len(cube_ks)} cubes")
print()

print("=" * 72)
print("CONCLUSION")
print("=" * 72)
print()
if len(glueable_pairs) == 0:
    print("NO Howe-glueable pairs among the 15 pairs of sextic twists.")
    print()
    # Determine which condition blocked each pair
    failed_by = {"H1": 0, "H2": 0, "H3": 0, "H2+H3": 0, "H1+H2": 0, "H1+H3": 0}
    for i,j in combinations(range(6),2):
        ni, nj = N[i], N[j]
        h1 = (ni != nj)
        h2 = (deg_pats[i] == deg_pats[j])
        h3 = (gcd_table[(i,j)] == 1)
        if not h1 and not h2 and not h3:
            failed_by["H1+H2"] += 1
        elif not h1 and not h3:
            failed_by["H1+H3"] += 1
        elif not h2 and not h3:
            failed_by["H2+H3"] += 1
        elif not h1:
            failed_by["H1"] += 1
        elif not h2:
            failed_by["H2"] += 1
        elif not h3:
            failed_by["H3"] += 1
    print("  Failure breakdown (why pairs are not glueable):")
    for reason, cnt in sorted(failed_by.items()):
        if cnt > 0:
            print(f"    {reason}: {cnt} pair(s)")
    print()
    print("  Implication: NO (2,2)-cover of type E_i × E_j (sextic-twist pair)")
    print("  exists for secp256k1 that passes all Howe conditions simultaneously.")
    print("  This is consistent with the structural-completeness paper:")
    print("  the obstruction here is structural (torsion mismatch or gcd>1),")
    print("  not a gap in the analysis.")
else:
    print(f"{len(glueable_pairs)} Howe-glueable pair(s) found: {glueable_pairs}")
    print()
    print("  For each glueable pair (i,j), a smooth genus-2 curve C/F_p")
    print("  exists with Jac(C) →(2,2)→ E_i × E_j (by Howe 1996, Theorem 1).")
    print()
    print("  ECDLP implication: NONE directly.")
    print("  Jac(C)(F_p) has order n_i*n_j ≈ p²; Pollard-ρ on Jac(C) costs ~p,")
    print("  which is WORSE than ECDLP on E_0=secp256k1 (cost ~√p).")
    print()
    print("  The existence of a (2,2)-cover is a structural fact, not an attack.")

print()
print("Done.")
print("=" * 72)
