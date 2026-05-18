#!/usr/bin/env python3
"""
Castryck-Decru gluing divisor map.

Translates the bottom half of `FromProdToJac` from richelot_aux.m: solves
the polynomial system that gives the Mumford form of the image divisor of
a point (P_c, P) ∈ E_α × E_β on the Kani-glued Jacobian J(h).

Field: F_{p²} = F_p[i]/(i² + 1), with p ≡ 3 (mod 4).
Inputs are F_{p²} elements represented as [a, b] meaning a + b·i.

Method:
  1. The beta-side equation eq3 (x_2 = P_x) is LINEAR in u_0. Solve
     symbolically:  u_0 = u_0(u_1, v_1).
  2. Substitute u_0 into eq4 (y_2 = P_y), which becomes linear in v_0.
     Solve:  v_0 = v_0(u_1, v_1).
  3. Substitute both into the remaining 4 equations (alpha-side eq1, eq2,
     and the two curve coefficients eq5a/eq5b). Now 4 equations in 2
     unknowns (u_1, v_1).
  4. Sympy Groebner over GF(p) (with `i` as an extra variable + relation
     i² + 1 = 0 for F_{p²} case). The reduced system is small enough that
     Groebner runs in well under a second even for F_{p²}.

The closed-form elimination is the breakthrough that makes F_{p²} feasible
without a native GF(p²) Groebner.

I/O: JSON over stdin/stdout.
"""

import json
import sys
from itertools import product

from sympy import Poly, groebner, symbols, together, numer, expand
from sympy.polys.domains import GF


def mod_inv(x, p):
    return pow(x % p, -1, p)


def pwr_mod_u(k, u0, u1):
    """x^k mod (x² + u1 x + u0) = a·x + b, returned as (a, b) sympy exprs."""
    a, b = 0, 1
    if k == 0:
        return a, b
    a, b = 1, 0
    for _ in range(2, k + 1):
        a, b = b - a * u1, -a * u0
    return a, b


# ---------------------------------------------------------------------------
# F_{p²} arithmetic, with elements as Python tuples (a, b) for a + b·i.
# ---------------------------------------------------------------------------

class F2:
    """F_{p²} element with explicit (a, b) = a + b·i representation."""
    __slots__ = ("a", "b", "p")

    def __init__(self, a, b, p):
        self.a = a % p
        self.b = b % p
        self.p = p

    @classmethod
    def from_list(cls, pair, p):
        return cls(pair[0], pair[1], p)

    def __add__(self, o): return F2(self.a + o.a, self.b + o.b, self.p)
    def __sub__(self, o): return F2(self.a - o.a, self.b - o.b, self.p)
    def __neg__(self):    return F2(-self.a, -self.b, self.p)
    def __mul__(self, o):
        if isinstance(o, int):
            return F2(self.a * o, self.b * o, self.p)
        # (a+bi)(c+di) = (ac-bd) + (ad+bc)i
        return F2(self.a * o.a - self.b * o.b,
                  self.a * o.b + self.b * o.a, self.p)
    def __pow__(self, e):
        result = F2(1, 0, self.p)
        base = F2(self.a, self.b, self.p)
        while e > 0:
            if e & 1:
                result = result * base
            base = base * base
            e >>= 1
        return result
    def inv(self):
        # (a+bi)^{-1} = (a-bi)/(a²+b²)
        n = (self.a * self.a + self.b * self.b) % self.p
        ni = mod_inv(n, self.p)
        return F2(self.a * ni, (-self.b * ni), self.p)
    def __truediv__(self, o):
        return self * o.inv() if isinstance(o, F2) else self * mod_inv(o, self.p)
    def is_zero(self):
        return self.a == 0 and self.b == 0
    def is_in_fp(self):
        return self.b == 0
    def __repr__(self):
        return f"F2({self.a}, {self.b})"
    def to_list(self):
        return [self.a, self.b]


# ---------------------------------------------------------------------------
# Gluing constants (work over F_{p²}).
# ---------------------------------------------------------------------------

def compute_constants(alpha, beta, p):
    """alpha, beta are lists of 3 F2 elements. Returns dict of F2 constants."""
    alp1, alp2, alp3 = alpha
    bet1, bet2, bet3 = beta

    a1 = (
        (alp3 - alp2) * (alp3 - alp2) / (bet3 - bet2)
        + (alp2 - alp1) * (alp2 - alp1) / (bet2 - bet1)
        + (alp1 - alp3) * (alp1 - alp3) / (bet1 - bet3)
    )
    b1 = (
        (bet3 - bet2) * (bet3 - bet2) / (alp3 - alp2)
        + (bet2 - bet1) * (bet2 - bet1) / (alp2 - alp1)
        + (bet1 - bet3) * (bet1 - bet3) / (alp1 - alp3)
    )
    a2 = (
        alp1 * (bet3 - bet2) + alp2 * (bet1 - bet3) + alp3 * (bet2 - bet1)
    )
    b2 = (
        bet1 * (alp3 - alp2) + bet2 * (alp1 - alp3) + bet3 * (alp2 - alp1)
    )
    delta_alp = (alp1 - alp2) * (alp1 - alp2) * (alp1 - alp3) * (alp1 - alp3) * (alp2 - alp3) * (alp2 - alp3)
    delta_bet = (bet1 - bet2) * (bet1 - bet2) * (bet1 - bet3) * (bet1 - bet3) * (bet2 - bet3) * (bet2 - bet3)
    big_a = delta_bet * a1 / a2
    big_b = delta_alp * b1 / b2
    s1 = (-big_b) / big_a * a2 / a1
    s2 = (
        alp1 * (alp3 - alp2) * (alp3 - alp2) / (bet3 - bet2)
        + alp2 * (alp1 - alp3) * (alp1 - alp3) / (bet1 - bet3)
        + alp3 * (alp2 - alp1) * (alp2 - alp1) / (bet2 - bet1)
    ) / a1
    t1 = (-big_a) / big_b * b2 / b1
    t2 = (
        bet1 * (bet3 - bet2) * (bet3 - bet2) / (alp3 - alp2)
        + bet2 * (bet1 - bet3) * (bet1 - bet3) / (alp1 - alp3)
        + bet3 * (bet2 - bet1) * (bet2 - bet1) / (alp2 - alp1)
    ) / b1
    return dict(
        a1=a1, b1=b1, a2=a2, b2=b2,
        delta_alp=delta_alp, delta_bet=delta_bet,
        big_a=big_a, big_b=big_b,
        s1=s1, s2=s2, t1=t1, t2=t2,
    )


def compute_h_coeffs(alpha, beta, p):
    """Compute h(x) coefficients (low to high), as F2 elements."""
    alp1, alp2, alp3 = alpha
    bet1, bet2, bet3 = beta
    K = compute_constants(alpha, beta, p)
    big_a, big_b = K["big_a"], K["big_b"]

    def f_factor(a_coef, b_coef):
        # x² · a + b → coeffs [b, 0, a] (low → high)
        return [b_coef, F2(0, 0, p), a_coef]

    f1 = f_factor(
        big_a * (alp2 - alp1) * (alp1 - alp3),
        big_b * (bet2 - bet1) * (bet1 - bet3),
    )
    f2 = f_factor(
        big_a * (alp3 - alp2) * (alp2 - alp1),
        big_b * (bet3 - bet2) * (bet2 - bet1),
    )
    f3 = f_factor(
        big_a * (alp1 - alp3) * (alp3 - alp2),
        big_b * (bet1 - bet3) * (bet3 - bet2),
    )

    def poly_mul(p1, p2):
        r = [F2(0, 0, p)] * (len(p1) + len(p2) - 1)
        for i_, ci in enumerate(p1):
            for j_, cj in enumerate(p2):
                r[i_ + j_] = r[i_ + j_] + ci * cj
        return r

    prod = poly_mul(poly_mul(f1, f2), f3)
    return [(-c) for c in prod]  # leading minus


def f2_to_sympy(x: F2, i_sym):
    """F2 element a + b·i → sympy expression a + b * i_sym."""
    return x.a + x.b * i_sym


# ---------------------------------------------------------------------------
# Main routine
# ---------------------------------------------------------------------------

def gluing_divisor_map(p, alpha_list, beta_list, pc_list, p_pt_list):
    alpha = [F2.from_list(x, p) for x in alpha_list]
    beta = [F2.from_list(x, p) for x in beta_list]
    pc = [F2.from_list(x, p) for x in pc_list]
    p_pt = [F2.from_list(x, p) for x in p_pt_list]

    K = compute_constants(alpha, beta, p)
    s1, s2, t1, t2 = K["s1"], K["s2"], K["t1"], K["t2"]
    big_a, big_b = K["big_a"], K["big_b"]
    delta_alp, delta_bet = K["delta_alp"], K["delta_bet"]
    h_coeffs = compute_h_coeffs(alpha, beta, p)

    # Detect if everything is in F_p (b=0 throughout). If yes, drop i to
    # speed up Groebner significantly.
    in_fp = all(
        x.is_in_fp() for x in alpha + beta + pc + p_pt
    ) and all(c.is_in_fp() for c in [s1, s2, t1, t2, big_a, big_b, delta_alp, delta_bet]) \
      and all(c.is_in_fp() for c in h_coeffs)

    if in_fp:
        sols = _solve_fp(p, alpha, beta, pc, p_pt, K, h_coeffs)
    else:
        sols = _solve_fp2(p, alpha, beta, pc, p_pt, K, h_coeffs)
    return {"solutions": sols, "h_coeffs": [[c.a, c.b] for c in h_coeffs]}


def _build_reduced_system(p, alpha, beta, pc, p_pt, K, h_coeffs, use_i):
    """Build the 4-equation, 2-unknown system after eliminating u_0, v_0 in
    closed form. Returns (polys, u0_expr, v0_expr, vars) where polys is a
    list of 4 (+1 if use_i) sympy expressions in (u_1, v_1[, i_sym]).

    `use_i = True` adds `i_sym` as an extra variable with i² + 1 = 0.
    """
    from sympy import together as _together, numer as _numer, solve as _solve
    s1, s2 = K["s1"], K["s2"]
    t1, t2 = K["t1"], K["t2"]
    big_a, big_b = K["big_a"], K["big_b"]
    delta_alp, delta_bet = K["delta_alp"], K["delta_bet"]

    if use_i:
        i_sym = symbols("i")
        def F(x): return f2_to_sympy(x, i_sym)
        s1_, s2_ = F(s1), F(s2)
        t1_, t2_ = F(t1), F(t2)
        big_a_, big_b_ = F(big_a), F(big_b)
        delta_alp_, delta_bet_ = F(delta_alp), F(delta_bet)
        alp_sum = F(alpha[0]) + F(alpha[1]) + F(alpha[2])
        bet_sum = F(beta[0]) + F(beta[1]) + F(beta[2])
        pc_x_, pc_y_ = F(pc[0]), F(pc[1])
        p_x_, p_y_ = F(p_pt[0]), F(p_pt[1])
        h_sym = [F(c) for c in h_coeffs]
    else:
        i_sym = None
        s1_, s2_ = s1.a, s2.a
        t1_, t2_ = t1.a, t2.a
        big_a_, big_b_ = big_a.a, big_b.a
        delta_alp_, delta_bet_ = delta_alp.a, delta_bet.a
        alp_sum = alpha[0].a + alpha[1].a + alpha[2].a
        bet_sum = beta[0].a + beta[1].a + beta[2].a
        pc_x_, pc_y_ = pc[0].a, pc[1].a
        p_x_, p_y_ = p_pt[0].a, p_pt[1].a
        h_sym = [c.a for c in h_coeffs]

    u0, u1, v0, v1 = symbols("u0 u1 v0 v1")

    u0_t = 1 / u0
    u1_t = u1 / u0
    v0_t = (u1 * v0 - u0 * v1) / u0 ** 2
    v1_t = (u1 ** 2 * v0 - u0 * v0 - u0 * u1 * v1) / u0 ** 2

    # Use symbolic inverses to keep field consistent
    lamb1 = -delta_bet_ * v1_t / (s1_ * u1_t * big_a_ ** 3)
    lamb2 = -delta_alp_ * v1 / (t1_ * u1 * big_b_ ** 3)

    x1 = lamb1 ** 2 + alp_sum - s1_ * (u1_t ** 2 - 2 * u0_t) - 2 * s2_
    y1 = -lamb1 * (x1 - s2_ + (u0_t * v1_t - u1_t * v0_t) * s1_ / v1_t)
    x2 = lamb2 ** 2 + bet_sum - t1_ * (u1 ** 2 - 2 * u0) - 2 * t2_
    y2 = -lamb2 * (x2 - t2_ + (u0 * v1 - u1 * v0) * t1_ / v1)

    eq1 = _numer(_together(x1 - pc_x_))
    eq2 = _numer(_together(y1 - pc_y_))
    eq3 = _numer(_together(x2 - p_x_))
    eq4 = _numer(_together(y2 - p_y_))

    v_sq_a = 2 * v0 * v1 - v1 ** 2 * u1
    v_sq_b = v0 ** 2 - v1 ** 2 * u0
    h_a, h_b = 0, 0
    for k, hk in enumerate(h_sym):
        ak, bk = pwr_mod_u(k, u0, u1)
        h_a += hk * ak
        h_b += hk * bk
    eq5a = expand(v_sq_a - h_a)
    eq5b = expand(v_sq_b - h_b)

    # Closed-form elimination: u_0 from eq3 (linear in u_0)
    u0_sols = _solve(eq3, u0)
    if not u0_sols:
        return None
    u0_expr = u0_sols[0]
    # v_0 from eq4 after substituting u_0
    eq4_sub = _numer(_together(eq4.subs(u0, u0_expr)))
    v0_sols = _solve(eq4_sub, v0)
    if not v0_sols:
        return None
    v0_expr = v0_sols[0]

    def sub_uv(e):
        return _numer(_together(e.subs([(v0, v0_expr), (u0, u0_expr)])))

    reduced = [sub_uv(eq1), sub_uv(eq2), sub_uv(eq5a), sub_uv(eq5b)]
    if use_i:
        reduced.append(i_sym ** 2 + 1)

    return {
        "polys": reduced,
        "u0_expr": u0_expr,
        "v0_expr": v0_expr,
        "vars": (u0, u1, v0, v1, i_sym) if use_i else (u0, u1, v0, v1),
    }


def _solve_fp(p, alpha, beta, pc, p_pt, K, h_coeffs):
    """F_p subcase: full Groebner (slower but covers trivial-solution case
    that the closed-form-elimination approach misses due to u_1=0
    singularity in the eliminating expression)."""
    s1, s2, t1, t2 = K["s1"].a, K["s2"].a, K["t1"].a, K["t2"].a
    big_a, big_b = K["big_a"].a, K["big_b"].a
    delta_alp, delta_bet = K["delta_alp"].a, K["delta_bet"].a
    alp = [a.a for a in alpha]
    bet = [b.a for b in beta]
    pc_x, pc_y = pc[0].a, pc[1].a
    p_x, p_y = p_pt[0].a, p_pt[1].a
    h_int = [c.a for c in h_coeffs]

    u0, u1, v0, v1 = symbols("u0 u1 v0 v1")

    u0_t = 1 / u0
    u1_t = u1 / u0
    v0_t = (u1 * v0 - u0 * v1) / u0 ** 2
    v1_t = (u1 ** 2 * v0 - u0 * v0 - u0 * u1 * v1) / u0 ** 2

    inv_a3 = mod_inv(pow(big_a, 3, p), p)
    inv_b3 = mod_inv(pow(big_b, 3, p), p)
    inv_s1 = mod_inv(s1, p)
    inv_t1 = mod_inv(t1, p)

    lamb1 = -delta_bet * inv_a3 * v1_t * inv_s1 / u1_t
    lamb2 = -delta_alp * inv_b3 * v1 * inv_t1 / u1

    x1 = lamb1 ** 2 + alp[0] + alp[1] + alp[2] - s1 * (u1_t ** 2 - 2 * u0_t) - 2 * s2
    y1 = -lamb1 * (x1 - s2 + (u0_t * v1_t - u1_t * v0_t) * s1 / v1_t)
    x2 = lamb2 ** 2 + bet[0] + bet[1] + bet[2] - t1 * (u1 ** 2 - 2 * u0) - 2 * t2
    y2 = -lamb2 * (x2 - t2 + (u0 * v1 - u1 * v0) * t1 / v1)

    eq1 = numer(together(x1 - pc_x))
    eq2 = numer(together(y1 - pc_y))
    eq3 = numer(together(x2 - p_x))
    eq4 = numer(together(y2 - p_y))

    v_sq_a = 2 * v0 * v1 - v1 ** 2 * u1
    v_sq_b = v0 ** 2 - v1 ** 2 * u0
    h_a, h_b = 0, 0
    for k, hk in enumerate(h_int):
        ak, bk = pwr_mod_u(k, u0, u1)
        h_a += hk * ak
        h_b += hk * bk
    eq5a = expand(v_sq_a - h_a)
    eq5b = expand(v_sq_b - h_b)

    polys = [Poly(e, u0, u1, v0, v1, domain=GF(p)) for e in [eq1, eq2, eq3, eq4, eq5a, eq5b]]
    g = groebner(polys, u0, u1, v0, v1, order="lex", domain=GF(p))
    basis = list(g)

    return _extract_solutions_fp(basis, [u0, u1, v0, v1], p)


def _extract_solutions_fp(basis, vars_, p):
    """Extract all solutions over GF(p) from a lex Groebner basis."""
    u0, u1, v0, v1 = vars_
    # Last poly should be univariate in v1
    last = basis[-1]
    last_poly = Poly(last, v1, domain=GF(p))
    if last_poly.free_symbols != {v1}:
        return [{"error": "last basis not univariate in v1", "last": str(last_poly)}]
    roots = list(last_poly.ground_roots().keys())
    sols = []
    for v1_val in roots:
        v1_int = int(v1_val) % p
        partial = {v1: v1_int}
        for target in [v0, u1, u0]:
            found = False
            for b in basis:
                sub = Poly(b, u0, u1, v0, v1, domain=GF(p)).subs(partial)
                try:
                    pol = Poly(sub, target, domain=GF(p))
                except Exception:
                    continue
                if pol.degree() >= 1 and pol.free_symbols == {target}:
                    rs = list(pol.ground_roots().keys())
                    if rs:
                        partial[target] = int(rs[0]) % p
                        found = True
                        break
            if not found:
                partial[target] = None
                break
        sols.append({
            "u0": [partial[u0], 0] if partial.get(u0) is not None else None,
            "u1": [partial[u1], 0] if partial.get(u1) is not None else None,
            "v0": [partial[v0], 0] if partial.get(v0) is not None else None,
            "v1": [partial[v1], 0] if partial.get(v1) is not None else None,
        })
    return sols


def _poly_to_fp2_coeffs(poly_expr, var, i_sym, p):
    """Parse a sympy expression (univariate in `var` after reducing mod i²+1)
    into a list of (a, b) int tuples, where coeffs[k] = (a, b) represents
    the F_{p²} coefficient a + b·i of var^k.
    """
    from sympy import Poly as SP, expand as _expand
    p_obj = SP(_expand(poly_expr), var, i_sym, domain=GF(p))
    max_var_deg = p_obj.degree(var)
    # Tuples (a, b) for each power of var
    coeffs = [(0, 0) for _ in range(max_var_deg + 1)]
    for monom, coef in p_obj.terms():
        var_deg, i_deg = monom
        c = int(coef) % p
        i_red = i_deg % 4
        a, b = coeffs[var_deg]
        if i_red == 0:
            a = (a + c) % p
        elif i_red == 1:
            b = (b + c) % p
        elif i_red == 2:
            a = (a - c) % p
        else:
            b = (b - c) % p
        coeffs[var_deg] = (a, b)
    return coeffs


def _fp2_roots_brute(coeffs, p):
    """Brute-force search for all (a, b) ∈ F_{p²} where Σ coeffs[k] · x^k = 0.
    Uses numpy vectorization across all p² candidates simultaneously.
    For p=431, finds all roots in ~50ms vs. ~30s for pure Python.
    """
    import numpy as np
    n = len(coeffs)
    if n == 0:
        return []
    # Build grid of (a, b) for all of F_{p²}
    a_grid, b_grid = np.meshgrid(np.arange(p, dtype=np.int64),
                                  np.arange(p, dtype=np.int64),
                                  indexing='ij')
    res_a = np.full((p, p), coeffs[n - 1][0], dtype=np.int64)
    res_b = np.full((p, p), coeffs[n - 1][1], dtype=np.int64)
    for k in range(n - 2, -1, -1):
        # res = res · x + coeffs[k]
        new_a = (res_a * a_grid - res_b * b_grid) % p
        new_b = (res_a * b_grid + res_b * a_grid) % p
        res_a = (new_a + coeffs[k][0]) % p
        res_b = (new_b + coeffs[k][1]) % p
    zero_mask = (res_a == 0) & (res_b == 0)
    ia, ib = np.where(zero_mask)
    return list(zip(ia.tolist(), ib.tolist()))


def _eval_fp2_poly(coeffs, x: F2, p):
    """Wrapper kept for backward compat — uses fast tuple form internally."""
    n = len(coeffs)
    if n == 0:
        return F2(0, 0, p)
    if isinstance(coeffs[0], F2):
        # Old style: convert
        tup = [(c.a, c.b) for c in coeffs]
    else:
        tup = coeffs
    res_a, res_b = tup[n - 1]
    for k in range(n - 2, -1, -1):
        new_a = (res_a * x.a - res_b * x.b) % p
        new_b = (res_a * x.b + res_b * x.a) % p
        res_a = (new_a + tup[k][0]) % p
        res_b = (new_b + tup[k][1]) % p
    return F2(res_a, res_b, p)


def _solve_fp2(p, alpha, beta, pc, p_pt, K, h_coeffs):
    """F_{p²} general case: closed-form elimination of u_0, v_0, then
    Groebner with `i` as extra variable (i²+1=0) over GF(p). Root
    extraction uses fast F_{p²} Horner evaluation after parsing coefficients
    once (vs. sympy-per-call which is ~1000× slower)."""
    built = _build_reduced_system(p, alpha, beta, pc, p_pt, K, h_coeffs, use_i=True)
    if built is None:
        return [{"error": "no closed-form for u_0 or v_0"}]
    u0, u1, v0, v1, i_sym = built["vars"]
    u0_expr = built["u0_expr"]
    v0_expr = built["v0_expr"]

    polys = [Poly(expand(e), i_sym, u1, v1, domain=GF(p)) for e in built["polys"]]
    g = groebner(polys, i_sym, u1, v1, order="lex", domain=GF(p))
    basis = list(g)

    # Find the basis poly univariate in v1 over GF(p²).
    candidate = None
    for b in reversed(basis):
        free = b.free_symbols
        if v1 in free and free.issubset({v1, i_sym}):
            candidate = b
            break
    if candidate is None:
        return [{"error": "no univariate-in-v1 basis poly found"}]
    candidate_expr = candidate.as_expr() if hasattr(candidate, "as_expr") else candidate

    # Parse to F_{p²} coefficients once, then evaluate fast via tuple form.
    v1_coeffs_tup = _poly_to_fp2_coeffs(candidate_expr, v1, i_sym, p)
    import time as _time
    _t0 = _time.time()
    v1_root_tuples = _fp2_roots_brute(v1_coeffs_tup, p)
    print(f"# v1 brute force ({len(v1_root_tuples)} roots): {_time.time()-_t0:.2f}s",
          file=sys.stderr)
    if not v1_root_tuples:
        return [{"error": "no v1 roots in F_{p²}"}]
    v1_roots = [F2(a, b, p) for (a, b) in v1_root_tuples]

    solutions = []
    for v1_val in v1_roots:
        # Find a basis poly that's univariate in u_1 after substituting v_1.
        u1_val = None
        for b in basis:
            b_expr = b.as_expr() if hasattr(b, "as_expr") else b
            free = b.free_symbols if hasattr(b, "free_symbols") else b_expr.free_symbols
            if u1 not in free:
                continue
            b_sub = b_expr.subs(v1, f2_to_sympy(v1_val, i_sym))
            try:
                u1_coeffs_tup = _poly_to_fp2_coeffs(b_sub, u1, i_sym, p)
            except Exception:
                continue
            if len(u1_coeffs_tup) < 2:
                continue
            _t0 = _time.time()
            u1_root_tuples = _fp2_roots_brute(u1_coeffs_tup, p)
            print(f"# u1 brute force ({len(u1_root_tuples)} roots): {_time.time()-_t0:.2f}s",
                  file=sys.stderr)
            if u1_root_tuples:
                u1_val = F2(u1_root_tuples[0][0], u1_root_tuples[0][1], p)
                break
        if u1_val is None:
            continue
        # u_0, v_0 from closed-form
        u0_sym = u0_expr.subs([(u1, f2_to_sympy(u1_val, i_sym)),
                               (v1, f2_to_sympy(v1_val, i_sym))])
        v0_sym = v0_expr.subs([(u1, f2_to_sympy(u1_val, i_sym)),
                               (v1, f2_to_sympy(v1_val, i_sym))])
        try:
            u0_val = _eval_in_fp2(u0_sym, {i_sym: F2(0, 1, p)}, p)
            v0_val = _eval_in_fp2(v0_sym, {i_sym: F2(0, 1, p)}, p)
        except Exception:
            continue
        solutions.append({
            "u0": u0_val.to_list(),
            "u1": u1_val.to_list(),
            "v0": v0_val.to_list(),
            "v1": v1_val.to_list(),
        })
    return solutions if solutions else [{"error": "no F_{p²} solutions extracted"}]


def _solve_fp2_OLD_FULL_SYSTEM(p, alpha, beta, pc, p_pt, K, h_coeffs):
    """[deprecated] Full 6-eq, 4-unknown Groebner over GF(p) with i."""
    i_sym = symbols("i")
    u0, u1, v0, v1 = symbols("u0 u1 v0 v1")

    # Cast all F2 constants to sympy expressions in i_sym.
    def F(x: F2):
        return f2_to_sympy(x, i_sym)

    s1, s2 = F(K["s1"]), F(K["s2"])
    t1, t2 = F(K["t1"]), F(K["t2"])
    big_a, big_b = F(K["big_a"]), F(K["big_b"])
    delta_alp, delta_bet = F(K["delta_alp"]), F(K["delta_bet"])

    u0_t = 1 / u0
    u1_t = u1 / u0
    v0_t = (u1 * v0 - u0 * v1) / u0 ** 2
    v1_t = (u1 ** 2 * v0 - u0 * v0 - u0 * u1 * v1) / u0 ** 2

    # Reduce big_a, big_b inverses pre-emptively. In F_{p²} we don't have
    # `pow(big_a, 3, p)` directly — compute big_a³ symbolically and let
    # Groebner handle the rest. We can't easily "invert mod i²+1" inside
    # symbolic expressions without normalizing, so we'll just multiply
    # through later when clearing denominators.

    lamb1_num = -delta_bet * v1_t
    lamb1_den = s1 * u1_t * big_a ** 3
    lamb2_num = -delta_alp * v1
    lamb2_den = t1 * u1 * big_b ** 3

    lamb1 = lamb1_num / lamb1_den
    lamb2 = lamb2_num / lamb2_den

    alp_sum = F(alpha[0]) + F(alpha[1]) + F(alpha[2])
    bet_sum = F(beta[0]) + F(beta[1]) + F(beta[2])
    x1 = lamb1 ** 2 + alp_sum - s1 * (u1_t ** 2 - 2 * u0_t) - 2 * s2
    y1 = -lamb1 * (x1 - s2 + (u0_t * v1_t - u1_t * v0_t) * s1 / v1_t)
    x2 = lamb2 ** 2 + bet_sum - t1 * (u1 ** 2 - 2 * u0) - 2 * t2
    y2 = -lamb2 * (x2 - t2 + (u0 * v1 - u1 * v0) * t1 / v1)

    pc_x = F(pc[0]); pc_y = F(pc[1])
    p_x = F(p_pt[0]); p_y = F(p_pt[1])
    eq1 = numer(together(x1 - pc_x))
    eq2 = numer(together(y1 - pc_y))
    eq3 = numer(together(x2 - p_x))
    eq4 = numer(together(y2 - p_y))

    v_sq_a = 2 * v0 * v1 - v1 ** 2 * u1
    v_sq_b = v0 ** 2 - v1 ** 2 * u0
    h_a_expr, h_b_expr = 0, 0
    for k, hk in enumerate(h_coeffs):
        ak, bk = pwr_mod_u(k, u0, u1)
        h_a_expr += F(hk) * ak
        h_b_expr += F(hk) * bk
    eq5a = expand(v_sq_a - h_a_expr)
    eq5b = expand(v_sq_b - h_b_expr)

    # Add i² + 1 = 0 as an extra relation.
    polys_expr = [eq1, eq2, eq3, eq4, eq5a, eq5b, i_sym ** 2 + 1]
    polys = []
    for e in polys_expr:
        polys.append(Poly(expand(e), i_sym, u0, u1, v0, v1, domain=GF(p)))

    # Lex order with i first → eliminate i early.
    g = groebner(polys, i_sym, u0, u1, v0, v1, order="lex", domain=GF(p))
    basis = list(g)

    solutions = _extract_solutions_fp2(basis, i_sym, [u0, u1, v0, v1], p)
    return {"solutions": solutions, "h_coeffs": [[c.a, c.b] for c in h_coeffs]}


def _extract_solutions_fp2(basis, i_sym, vars_, p):
    """Extract solutions over F_{p²} from a lex Groebner basis whose
    polynomials may involve `i_sym` (with i² + 1 = 0 in the system).

    Strategy:
      1. For each variable from smallest (v1) up, find a basis polynomial
         that's univariate in that variable after substituting the larger
         vars (with i_sym treated as an F_{p²} element).
      2. Find roots of that polynomial in F_{p²} by brute force.
      3. Back-substitute.

    Brute force is fine for our toy parameters (p ≈ 431, |F_{p²}| ≈ 186k).
    """
    u0, u1, v0, v1 = vars_

    def is_only_in(b, allowed):
        """Does basis poly b only involve symbols in `allowed`?"""
        return all(s in allowed or s == i_sym for s in b.free_symbols)

    def find_univariate(target, partial):
        """Find a basis poly that's univariate in `target` after substituting
        `partial` (a dict {var: F2}). Returns the poly's evaluator (callable
        F2 → F2) or None."""
        for b in basis:
            sub = b
            try:
                # Substitute the already-known F_{p²} values as polynomials in i_sym.
                subs_map = {sym: f2_to_sympy(val, i_sym) for sym, val in partial.items()}
                if subs_map:
                    sub = expand(b.as_expr().xreplace(subs_map))
            except Exception:
                continue
            # Check: after substitution, depends only on `target` and i_sym.
            free = sub.free_symbols if hasattr(sub, "free_symbols") else set()
            if target not in free:
                continue
            allowed = {target, i_sym}
            if not free.issubset(allowed):
                continue
            # Return an evaluator: given F2 value `target_val`, returns the
            # poly evaluated in F_{p²}.
            sub_expr = sub
            def eval_at(target_val, _sub=sub_expr):
                return _eval_in_fp2(_sub, {i_sym: F2(0, 1, p), target: target_val}, p)
            return eval_at
        return None

    def f2_roots(eval_fn):
        """Brute-force find all F_{p²} roots of a poly whose evaluator is `eval_fn`."""
        roots = []
        for a in range(p):
            for b in range(p):
                if eval_fn(F2(a, b, p)).is_zero():
                    roots.append(F2(a, b, p))
        return roots

    # Find roots of last basis polynomial (in v1)
    eval_v1 = find_univariate(v1, {})
    if eval_v1 is None:
        return [{"error": "no univariate-in-v1 basis poly found"}]
    v1_roots = f2_roots(eval_v1)
    if not v1_roots:
        return [{"error": "no v1 roots in F_{p²}"}]

    solutions = []
    for v1_val in v1_roots:
        partial = {v1: v1_val}
        ok = True
        for target in [v0, u1, u0]:
            eval_t = find_univariate(target, partial)
            if eval_t is None:
                ok = False
                break
            rs = f2_roots(eval_t)
            if not rs:
                ok = False
                break
            partial[target] = rs[0]  # pick first root
        if not ok:
            continue
        solutions.append({
            "u0": partial[u0].to_list(),
            "u1": partial[u1].to_list(),
            "v0": partial[v0].to_list(),
            "v1": partial[v1].to_list(),
        })
    return solutions if solutions else [{"error": "no full solutions extracted"}]


def _eval_in_fp2(sympy_poly_expr, subs_map, p):
    """Evaluate a sympy expression at F2-valued points. Substitutions: dict
    of {sympy_symbol: F2}. Returns an F2 element."""
    # Convert sympy expression to nested arithmetic in F2.
    # We assume the expression is a polynomial in the substituted variables.
    expr = sympy_poly_expr.as_expr() if hasattr(sympy_poly_expr, "as_expr") else sympy_poly_expr
    expr = expand(expr)
    # Collect terms
    from sympy import Poly as SP, Symbol
    syms = list(subs_map.keys())
    pol = SP(expr, *syms, domain=GF(p))
    result = F2(0, 0, p)
    for monom, coef in pol.terms():
        term = F2(int(coef), 0, p)
        for s, e in zip(syms, monom):
            for _ in range(e):
                term = term * subs_map[s]
        result = result + term
    return result


def main():
    data = json.load(sys.stdin)
    res = gluing_divisor_map(
        p=data["p"],
        alpha_list=data["alpha"],
        beta_list=data["beta"],
        pc_list=data["pc"],
        p_pt_list=data["p_pt"],
    )
    print(json.dumps(res, default=str))


if __name__ == "__main__":
    main()
