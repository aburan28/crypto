"""
Semaev summation polynomials over F_p (toy scale) and a small multivariate
resultant, used by the relation-yield / Groebner-proxy harness.

S_m(x_1,...,x_m) vanishes iff there exist points P_i on E with x(P_i)=x_i and
P_1 + ... + P_m = O.  Recurrence:

    S_2(x1,x2) = x1 - x2
    S_3(x1,x2,x3) = (x1-x2)^2 x3^2
                    - 2((x1+x2)(x1 x2 + a) + 2b) x3
                    + ((x1 x2 - a)^2 - 4 b (x1+x2))
    S_m(x1,...,x_m) = Res_X( S_{m-1}(x1,...,x_{m-2},X), S_3(x_{m-1},x_m,X) )

Polynomials are represented as dicts {monomial_exponent_tuple: coeff} in the
ring F_p[x1,...,xk].  Only used at toy scale.
"""

from __future__ import annotations


def _zero():
    return {}


def padd(A, B, p):
    R = dict(A)
    for m, c in B.items():
        R[m] = (R.get(m, 0) + c) % p
        if R[m] == 0:
            del R[m]
    return R


def pscale(A, s, p):
    s %= p
    if s == 0:
        return {}
    return {m: (c * s) % p for m, c in A.items()}


def pmul(A, B, p, nvars):
    R = {}
    for ma, ca in A.items():
        for mb, cb in B.items():
            m = tuple(ma[i] + mb[i] for i in range(nvars))
            R[m] = (R.get(m, 0) + ca * cb) % p
            if R[m] == 0:
                del R[m]
    return R


def var(i, nvars):
    e = [0] * nvars
    e[i] = 1
    return {tuple(e): 1}


def const(c, nvars):
    if c == 0:
        return {}
    return {tuple([0] * nvars): c}


def semaev3(a, b, p, nvars=3, i1=0, i2=1, i3=2):
    """S_3 in variables x_{i1},x_{i2},x_{i3} of an nvars-variable ring."""
    x1 = var(i1, nvars)
    x2 = var(i2, nvars)
    x3 = var(i3, nvars)

    def mul(*fs):
        r = const(1, nvars)
        for f in fs:
            r = pmul(r, f, p, nvars)
        return r

    x1mx2 = padd(x1, pscale(x2, -1, p), p)
    x1px2 = padd(x1, x2, p)
    x1x2 = mul(x1, x2)

    term1 = mul(mul(x1mx2, x1mx2), mul(x3, x3))
    inner = padd(mul(x1px2, padd(x1x2, const(a, nvars), p)), const(2 * b, nvars), p)
    term2 = pscale(mul(inner, x3), -2, p)
    t3a = padd(x1x2, const(-a, nvars), p)
    term3 = padd(mul(t3a, t3a), pscale(x1px2, -4 * b, p), p)
    return padd(padd(term1, term2, p), term3, p)


def as_univariate(poly, var_index, nvars, other_values, p):
    """Specialise all variables except var_index to numeric values; return a
    dense coefficient list (low->high) in the remaining variable."""
    coeffs = {}
    for m, c in poly.items():
        val = c
        for i in range(nvars):
            if i == var_index:
                continue
            if m[i]:
                val = (val * pow(other_values[i], m[i], p)) % p
        d = m[var_index]
        coeffs[d] = (coeffs.get(d, 0) + val) % p
    if not coeffs:
        return [0]
    deg = max(coeffs)
    return [coeffs.get(i, 0) % p for i in range(deg + 1)]
