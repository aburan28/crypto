#!/usr/bin/env python3
"""
Thread 20: Python verification of Igusa-Clebsch invariants for secp256k1 naive cover.

Implements Clebsch A,B,C,D transvectants then converts to I2,I4,I6,I10 using the
Sage formula (clebsch_to_igusa, attributed to [LY2001]):
  I2  = -120 * A
  I4  = -720 * A^2 + 6750 * B
  I6  = 8640 * A^3 - 108000 * A*B + 202500 * C
  I10 = -62208*A^5 + 972000*A^3*B + 1620000*A^2*C
         - 3037500*A*B^2 - 6075000*B*C - 4556250*D

References:
  - igusa_clebsch_complete.gp (PARI port in this repo)
  - Sage source: sage/schemes/hyperelliptic_curves/invariants.py
  - Cardona-Quer 2005 Appendix A
  - chlrs_fp3_rosenhain.gp claim: I2=-87024 for h=(x^3+7)(x^3+189)

Goal: verify/correct the hardcoded values in chlrs_fp3_rosenhain.gp line 231.
"""

from fractions import Fraction
import math

# ----------------------------------------------------------------
# Binary form representation: h(x,y) of degree d is stored as
#   [c_0, c_1, ..., c_d]  where c_i = coefficient of x^i y^{d-i}
# For a univariate polynomial h(x) of degree ≤ d, homogenize with y.
# ----------------------------------------------------------------

def poly_to_binary(coeffs_list, deg):
    """coeffs_list[i] = coefficient of x^i; pad to length deg+1."""
    v = [Fraction(0)] * (deg + 1)
    for i, c in enumerate(coeffs_list):
        if i <= deg:
            v[i] = Fraction(c)
    return v

def poly_coeffs(poly_str_or_list):
    """Accept a coefficient list [a0, a1, ..., a6] for a0 + a1*x + ... + a6*x^6."""
    return poly_str_or_list

# ----------------------------------------------------------------
# k-th transvectant of binary forms f (degree m) and g (degree n).
# Returns coefficient vector of the result (degree m+n-2k).
#
# Formula:
#   (f,g)_k = ((m-k)!(n-k)! / (m!n!)) * sum_{i=0}^{k} (-1)^i C(k,i)
#             * (d^{k-i}/dx^{k-i} d^i/dy^i f) * (d^i/dx^i d^{k-i}/dy^{k-i} g)
#
# The partial d^s/dx^s d^t/dy^t of c_j * x^j y^{deg-j} gives
#   c_j * j!/(j-s)! * (deg-j)!/(deg-j-t)! * x^{j-s} y^{deg-j-t}   [if j>=s, deg-j>=t]
# ----------------------------------------------------------------

def transvectant(f, m, g, n, k):
    out_deg = m + n - 2 * k
    if out_deg < 0:
        return [Fraction(0)]
    out = [Fraction(0)] * (out_deg + 1)
    scale = Fraction(math.factorial(m - k) * math.factorial(n - k),
                     math.factorial(m) * math.factorial(n))

    for i in range(k + 1):
        sign = (-1) ** i
        binom = math.comb(k, i)
        s_f = k - i  # order of x-deriv on f
        t_f = i      # order of y-deriv on f
        s_g = i      # order of x-deriv on g
        t_g = k - i  # order of y-deriv on g

        for a in range(s_f, m - t_f + 1):   # a = x-exponent of f term (after deriv: a-s_f)
            if a < 0 or a > m:
                continue
            # falling factorial for x: a*(a-1)*...*(a-s_f+1)
            fall_x_f = 1
            for jj in range(a, a - s_f, -1):
                fall_x_f *= jj
            # falling factorial for y: (m-a)*(m-a-1)*...*(m-a-t_f+1)
            fall_y_f = 1
            for jj in range(m - a, m - a - t_f, -1):
                fall_y_f *= jj
            coeff_f = f[a] * fall_x_f * fall_y_f
            if coeff_f == 0:
                continue

            for b in range(s_g, n - t_g + 1):   # b = x-exponent of g term
                if b < 0 or b > n:
                    continue
                fall_x_g = 1
                for jj in range(b, b - s_g, -1):
                    fall_x_g *= jj
                fall_y_g = 1
                for jj in range(n - b, n - b - t_g, -1):
                    fall_y_g *= jj
                coeff_g = g[b] * fall_x_g * fall_y_g
                if coeff_g == 0:
                    continue

                j_out = (a - s_f) + (b - s_g)   # x-exponent in output
                if 0 <= j_out <= out_deg:
                    out[j_out] += sign * binom * coeff_f * coeff_g

    for j in range(out_deg + 1):
        out[j] *= scale
    return out


def scalar_from_form(v):
    """Extract constant from a degree-0 binary form."""
    assert len(v) == 1, f"Expected degree-0 form, got len={len(v)}: {v}"
    return v[0]


# ----------------------------------------------------------------
# Clebsch invariants
# ----------------------------------------------------------------

def clebsch_ABCD(h_coeffs):
    """
    h_coeffs = [a0, a1, a2, a3, a4, a5, a6] (coefficient of x^i)
    Returns (A, B, C, D) as Fraction objects.
    """
    hv = [Fraction(c) for c in h_coeffs]
    # i = (f,f)_4, degree 4
    iv = transvectant(hv, 6, hv, 6, 4)
    # A = (f,f)_6, degree 0
    A = scalar_from_form(transvectant(hv, 6, hv, 6, 6))
    # B = (i,i)_4, degree 0
    B = scalar_from_form(transvectant(iv, 4, iv, 4, 4))
    # Delta = (i,i)_2, degree 4
    Delta = transvectant(iv, 4, iv, 4, 2)
    # C = (i, Delta)_4, degree 0
    C = scalar_from_form(transvectant(iv, 4, Delta, 4, 4))
    # y1 = (f,i)_4, degree 2
    y1 = transvectant(hv, 6, iv, 4, 4)
    # y2 = (i,y1)_2, degree 2
    y2 = transvectant(iv, 4, y1, 2, 2)
    # y3 = (i,y2)_2, degree 2
    y3 = transvectant(iv, 4, y2, 2, 2)
    # D = (y3,y1)_2, degree 0
    D = scalar_from_form(transvectant(y3, 2, y1, 2, 2))
    return A, B, C, D


def igusa_quadruple(h_coeffs):
    """Returns (I2, I4, I6, I10) as Fraction objects."""
    A, B, C, D = clebsch_ABCD(h_coeffs)
    I2  = -120 * A
    I4  = -720 * A**2 + 6750 * B
    I6  =  8640 * A**3 - 108000 * A * B + 202500 * C
    I10 = (-62208 * A**5 + 972000 * A**3 * B + 1620000 * A**2 * C
           - 3037500 * A * B**2 - 6075000 * B * C - 4556250 * D)
    return I2, I4, I6, I10


def poldisc(h_coeffs):
    """Polynomial discriminant of h = sum h_coeffs[i]*x^i (exact integer/fraction)."""
    # Use resultant(h, h') / leading_coeff
    # For our purposes, use a simple product formula over roots would be complex.
    # Instead, implement via sylvester matrix determinant.
    from fractions import Fraction
    n = len(h_coeffs) - 1  # degree
    # trim trailing zeros
    while n > 0 and h_coeffs[n] == 0:
        n -= 1
    coeffs = [Fraction(h_coeffs[i]) for i in range(n + 1)]

    # derivative coefficients
    dcoeffs = [Fraction(i) * coeffs[i] for i in range(1, n + 1)]

    # Sylvester matrix of h (degree n) and h' (degree n-1): size (2n-1) x (2n-1)
    size = 2 * n - 1
    mat = [[Fraction(0)] * size for _ in range(size)]
    # First n-1 rows: shift of h' (degree n-1)
    for row in range(n - 1):
        for col in range(n):
            if row + col < size:
                mat[row][row + col] = dcoeffs[n - 1 - col]
    # Last n rows: shift of h (degree n)
    for row in range(n):
        for col in range(n + 1):
            if (n - 1 + row) < size and col < size:
                mat[n - 1 + row][(n - 1 + row) - (n - col)] = coeffs[n - col] if (n - 1 + row) - (n - col) >= 0 else Fraction(0)

    # Actually, let's use a direct sylvester approach more carefully
    # Sylvester matrix: rows 0..n-2 are f shifted, rows n-1..2n-2 are f' shifted
    mat = [[Fraction(0)] * size for _ in range(size)]
    # n-1 shifts of f' (degree n-1, so n coefficients) starting at column positions 0..n-2
    for row in range(n - 1):
        for j in range(n):
            mat[row][row + j] = dcoeffs[n - 1 - j]
    # n shifts of f (degree n, so n+1 coefficients) starting at column positions 0..n-1
    for row in range(n):
        for j in range(n + 1):
            mat[n - 1 + row][row + j] = coeffs[n - j]

    # Gaussian elimination to get determinant
    det = Fraction(1)
    m = [row[:] for row in mat]
    for col in range(size):
        # Find pivot
        pivot_row = None
        for row in range(col, size):
            if m[row][col] != 0:
                pivot_row = row
                break
        if pivot_row is None:
            return Fraction(0)
        if pivot_row != col:
            m[col], m[pivot_row] = m[pivot_row], m[col]
            det *= -1
        det *= m[col][col]
        inv = Fraction(1, 1) / m[col][col]
        for row in range(col + 1, size):
            if m[row][col] != 0:
                factor = m[row][col] * inv
                for c in range(col, size):
                    m[row][c] -= factor * m[col][c]

    # disc(h) = (-1)^{n(n-1)/2} / a_n * Res(h, h')
    an = coeffs[n]
    sign = (-1) ** (n * (n - 1) // 2)
    return Fraction(sign) * det / an


# ----------------------------------------------------------------
# Tests
# ----------------------------------------------------------------

print("=" * 60)
print("Thread 20: Python Igusa-Clebsch verification")
print("=" * 60)
print()

# Test 1: h = x^6 + 1  =>  coeffs = [1,0,0,0,0,0,1]
print("---- Test 1: h = x^6 + 1 ----")
h1 = [1, 0, 0, 0, 0, 0, 1]
I2, I4, I6, I10 = igusa_quadruple(h1)
print(f"  I2  = {I2}  (expected: -120)")
print(f"  I4  = {I4}")
print(f"  I6  = {I6}")
print(f"  I10 = {I10}")
print(f"  I2 == -120? {I2 == -120}")
print()

# Test 2: h = x^6 + 2x^5 + 3x^4 + 5x^3 + 7x^2 + 11x + 13
print("---- Test 2: h = x^6 + 2x^5 + 3x^4 + 5x^3 + 7x^2 + 11x + 13 ----")
h2 = [13, 11, 7, 5, 3, 2, 1]   # [a0, a1, ..., a6]
I2, I4, I6, I10 = igusa_quadruple(h2)
# Expected I2 = -120*13*1 + 20*11*2 - 8*7*3 + 3*5^2 = -1560 + 440 - 168 + 75 = -1213
print(f"  I2  = {I2}  (expected: -1213)")
print(f"  I4  = {I4}")
print(f"  I6  = {I6}")
print(f"  I10 = {I10}")
print(f"  I2 == -1213? {I2 == -1213}")
print()

# Test 3: isomorphism invariance under x -> 2x + 1
print("---- Test 3: isomorphism invariance (x -> 2x+1 on h_orig) ----")
h_orig_coeffs = [7, 5, 0, 2, 3, 0, 1]  # x^6 + 3x^4 + 2x^3 + 5x + 7
I2o, I4o, I6o, I10o = igusa_quadruple(h_orig_coeffs)
print(f"  h_orig: I2={I2o}, I4={I4o}, I6={I6o}")

# Apply x -> 2x+1: h_iso(x) = 2^6 * h_orig((x-1)/2) ... computed symbolically
# For a simpler check, just verify the ratios I4/I2^2 and I6/I2^3 are preserved.
# h_iso: substitute x -> 2x+1 in h_orig and extract coefficients
def compose_affine(h_c, u, v):
    """Return coefficients of u^6 * h((x-v)/u) as a polynomial in x, for isomorphism."""
    # h(y) = sum h_c[i] y^i, y=(x-v)/u, result = u^6 * h((x-v)/u)
    # = u^6 * sum h_c[i] * ((x-v)/u)^i = sum h_c[i] * u^(6-i) * (x-v)^i
    # Expand (x-v)^i using binomial theorem
    n = len(h_c) - 1
    result = [Fraction(0)] * (n + 1)
    for i, ci in enumerate(h_c):
        factor = Fraction(u) ** (n - i) * Fraction(ci)
        # Expand (x-v)^i = sum_j C(i,j) x^j (-v)^{i-j}
        for j in range(i + 1):
            coeff = factor * math.comb(i, j) * (Fraction(-v) ** (i - j))
            result[j] += coeff
    return result

h_iso_coeffs = compose_affine(h_orig_coeffs, 2, 1)
I2i, I4i, I6i, I10i = igusa_quadruple(h_iso_coeffs)
print(f"  h_iso:  I2={I2i}, I4={I4i}, I6={I6i}")
r4_orig = I4o / I2o**2 if I2o != 0 else None
r4_iso  = I4i / I2i**2 if I2i != 0 else None
r6_orig = I6o / I2o**3 if I2o != 0 else None
r6_iso  = I6i / I2i**3 if I2i != 0 else None
print(f"  I4/I2^2 preserved? {r4_orig == r4_iso}  ({r4_orig})")
print(f"  I6/I2^3 preserved? {r6_orig == r6_iso}  ({r6_orig})")
print()

# ----------------------------------------------------------------
# MAIN RESULT: secp256k1 naive cover  h = (x^3+7)(x^3+189)
# ----------------------------------------------------------------
print("=" * 60)
print("MAIN: secp256k1 naive cover h = (x^3+7)(x^3+189)")
print("=" * 60)
print()
# Expand: (x^3+7)(x^3+189) = x^6 + (7+189)x^3 + 7*189
# = x^6 + 196x^3 + 1323
# coeffs: a0=1323, a1=0, a2=0, a3=196, a4=0, a5=0, a6=1
h_secp = [1323, 0, 0, 196, 0, 0, 1]
print(f"h = x^6 + 196*x^3 + 1323  (coeffs a0..a6: {h_secp})")
print()

print("Computing Clebsch ABCD (exact rational arithmetic)...")
A, B, C, D = clebsch_ABCD(h_secp)
print(f"  A = {A}")
print(f"  B = {B}")
print(f"  C = {C}")
print(f"  D = {D}")
print()

I2, I4, I6, I10 = igusa_quadruple(h_secp)
print("Igusa-Clebsch invariants over Q:")
print(f"  I2  = {I2}")
print(f"  I4  = {I4}")
print(f"  I6  = {I6}")
print(f"  I10 = {I10}")
print()

# Compare with hardcoded claim in chlrs_fp3_rosenhain.gp line 231
I2_claim  = -87024
I4_claim  = 19302628212
I6_claim  = -636669361341720
I10_claim = 46374105383717408990784
print("Comparison with chlrs_fp3_rosenhain.gp line-231 claim:")
print(f"  I2  claimed={I2_claim},  computed={I2},  match={I2==I2_claim}")
print(f"  I4  claimed={I4_claim},  computed={I4},  match={I4==I4_claim}")
print(f"  I6  claimed={I6_claim},  computed={I6},  match={I6==I6_claim}")
print(f"  I10 claimed={I10_claim},  computed={I10},  match={I10==I10_claim}")
print()

# Also check I2 via explicit formula: -120*a0*a6 + 20*a1*a5 - 8*a2*a4 + 3*a3^2
a0, a1, a2, a3, a4, a5, a6 = h_secp
I2_explicit = -120*a0*a6 + 20*a1*a5 - 8*a2*a4 + 3*a3**2
print(f"I2 via explicit formula (-120*a0*a6 + 20*a1*a5 - 8*a2*a4 + 3*a3^2):")
print(f"  = -120*{a0}*{a6} + 20*{a1}*{a5} - 8*{a2}*{a4} + 3*{a3}^2")
print(f"  = {-120*a0*a6} + 0 - 0 + {3*a3**2}")
print(f"  = {I2_explicit}")
print(f"  Matches transvectant I2? {Fraction(I2_explicit) == I2}")
print()

# Normalised ratios
if I2 != 0:
    print("Normalised invariant ratios (isomorphism class markers in A_2):")
    print(f"  I4/I2^2  = {I4/I2**2}")
    print(f"  I6/I2^3  = {I6/I2**3}")
    print(f"  I10/I2^5 = {I10/I2**5}")
    print()

# ----------------------------------------------------------------
# Naive cover with different twist coefficient (b_twist via d=3)
# b_twist = d^3 * b = 27 * 7 = 189 ← this IS what we computed
# Check if d=2 gives different result
# ----------------------------------------------------------------
print("---- Naive cover with d=2: h = (x^3+7)(x^3+56) ----")
# b_twist = 2^3 * 7 = 56
h_d2 = [7*56, 0, 0, 7+56, 0, 0, 1]  # x^6 + 63x^3 + 392
print(f"h = x^6 + {7+56}x^3 + {7*56}  (d=2, b_twist={8*7})")
I2d2, I4d2, I6d2, I10d2 = igusa_quadruple(h_d2)
print(f"  I2={I2d2}  I4={I4d2}  I6={I6d2}")
if I2d2 != 0:
    print(f"  I4/I2^2 = {I4d2/I2d2**2}  (same as d=3? {I4d2/I2d2**2 == I4/I2**2})")
    print(f"  I6/I2^3 = {I6d2/I2d2**3}  (same as d=3? {I6d2/I2d2**3 == I6/I2**3})")
print()

# ----------------------------------------------------------------
# p_secp reduction check (modular arithmetic spot check)
# ----------------------------------------------------------------
p_secp = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F
print(f"Reduction mod p_secp (secp256k1 prime):")
I2_mod  = int(I2)  % p_secp
I4_mod  = int(I4)  % p_secp
I6_mod  = int(I6)  % p_secp
I10_mod = int(I10) % p_secp
print(f"  I2  mod p = {I2_mod}")
print(f"  I4  mod p = {I4_mod}")
print(f"  I6  mod p = {I6_mod}")
print(f"  I10 mod p = {I10_mod}")
print()
print("Note: I2=0 mod p would mean naive cover is F_p-isomorphic to ANOTHER")
print("model (not a smoothness obstruction since I10 is the discriminant).")
print(f"  I2 ≠ 0 mod p? {I2_mod != 0}")
print(f"  I10 ≠ 0 mod p? {I10_mod != 0}")
print()

# ----------------------------------------------------------------
# Summary
# ----------------------------------------------------------------
print("=" * 60)
print("Summary")
print("=" * 60)
print()
print("Igusa-Clebsch quadruple of naive cover y^2 = (x^3+7)(x^3+189) over Q:")
print(f"  I2  = {int(I2)}")
print(f"  I4  = {int(I4)}")
print(f"  I6  = {int(I6)}")
print(f"  I10 = {int(I10)}")
print()
print("Comparison with chlrs_fp3_rosenhain.gp hardcoded claim:")
claim_ok = (I2 == I2_claim and I4 == I4_claim and
            I6 == I6_claim and I10 == I10_claim)
if claim_ok:
    print("  All 4 invariants MATCH — claim in chlrs_fp3_rosenhain.gp is CORRECT.")
else:
    print("  MISMATCH — claim in chlrs_fp3_rosenhain.gp needs correction.")
    print(f"  Correct values: I2={int(I2)}, I4={int(I4)},")
    print(f"    I6={int(I6)}, I10={int(I10)}")
print()
print("BLOCKED remaining: Igusa invariants of the ACTUAL Howe-glued curve")
print("  (not the naive cover) require Mestre reconstruction (Sage/Magma).")
print("  Naive cover and Howe-glued curve are in DIFFERENT isogeny classes.")
print("  Naive cover Jac: char poly = (T^2-tT+p)(T^2+tT+p) -- this holds only")
print("  when y^2=(x^3+b1)(x^3+b2) has Jac ~ E1xE2 (need to verify).")
