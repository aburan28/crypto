#!/usr/bin/env python3
"""msolve bridge for the F_{p²} gluing divisor map.

msolve (https://msolve.lip6.fr) is a fast Groebner basis solver implemented
in C with native finite-field support. For our reduced gluing system over
F_p[i]/(i²+1), encoded as polynomials in (i, u_1, v_1) over GF(p), msolve
solves in ~0.1 seconds where sympy hangs indefinitely.

I/O is the same JSON contract as gluing_divisor_map.py but the inner solve
uses an msolve subprocess instead of sympy Groebner.
"""

import json
import os
import re
import subprocess
import sys
import tempfile
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gluing_divisor_map import (
    F2,
    compute_constants,
    compute_h_coeffs,
    _build_reduced_system,
    _eval_in_fp2,
    _solve_fp,
    f2_to_sympy,
    mod_inv,
)

import numpy as np
from sympy import Poly, expand, symbols
from sympy.polys.domains import GF


def call_msolve(polys_str, varnames, p, msolve_bin="msolve"):
    """Run msolve. polys_str is comma-newline-separated string of polynomials
    in the variables (separator '^' for powers). Returns the parsed output."""
    with tempfile.NamedTemporaryFile("w", suffix=".ms", delete=False) as f:
        f.write(", ".join(varnames) + "\n")
        f.write(f"{p}\n")
        f.write(polys_str)
        in_path = f.name
    out_path = in_path + ".out"
    try:
        subprocess.run(
            [msolve_bin, "-f", in_path, "-o", out_path],
            check=True, capture_output=True, timeout=120,
        )
        with open(out_path) as f:
            raw = f.read()
        return _parse_msolve(raw, varnames)
    finally:
        try: os.unlink(in_path)
        except: pass
        try: os.unlink(out_path)
        except: pass


def _parse_msolve(raw, varnames):
    """Parse msolve's rational-parametrization output.

    Output shape (single-string excerpt):
      [0, [p, nvars, _, [var, var, ..., A], [linform_coefs],
           [1, [[deg_w, w_coefs], [deg_d, d_coefs], [[[deg_q, q_coefs]], ...]]]]]:
    """
    raw = raw.strip().rstrip(":")
    # Strip outer [0, [...]]
    obj = _parse_lisp(raw)
    # obj = [0, [p, nvars, _, vars, linform, type_data]]
    inner = obj[1]
    p = inner[0]
    n_msolve_vars = inner[1]
    msolve_varnames = inner[3]   # includes the introduced parameter "A"
    linform = inner[4]
    type_data = inner[5]
    type_tag = type_data[0]
    if type_tag == -1:
        return {"empty": True}
    param = type_data[1]
    w_deg, w_coefs = param[0]
    d_deg, d_coefs = param[1]
    var_params = param[2]
    var_results = {}
    # Each var_params[k] is [[deg, coefs]] for var k. There are len(varnames)
    # entries (one per original variable, the last "A" doesn't get its own).
    for vk, var in enumerate(varnames):
        if vk >= len(var_params):
            break
        deg_q, q_coefs = var_params[vk][0]
        var_results[var] = (deg_q, q_coefs)
    return {
        "p": p,
        "linform": linform,
        "w": (w_deg, w_coefs),
        "denom": (d_deg, d_coefs),
        "vars": var_results,
    }


def _parse_lisp(s):
    """Quick-and-dirty parser for msolve's Python-list-style output."""
    # Use Python eval after light cleanup. msolve outputs commas, brackets,
    # ints, and 'string' literals — all valid Python.
    return eval(s, {"__builtins__": {}}, {})


def f2_poly_eval(coefs, x_a, x_b, p):
    """Horner: evaluate Σ coefs[k]·x^k at x = (x_a + x_b·i) ∈ F_{p²}.
    coefs is a list of F_p ints. Returns (a, b)."""
    n = len(coefs)
    res_a, res_b = coefs[n - 1] % p, 0
    for k in range(n - 2, -1, -1):
        new_a = (res_a * x_a - res_b * x_b) % p
        new_b = (res_a * x_b + res_b * x_a) % p
        res_a = (new_a + coefs[k]) % p
        res_b = new_b
    return res_a, res_b


def fp2_roots_numpy(coefs, p):
    """Numpy brute force: find all (a, b) ∈ F_{p²} with Σ coefs[k]·(a+bi)^k = 0."""
    n = len(coefs)
    if n == 0:
        return []
    a_grid, b_grid = np.meshgrid(np.arange(p, dtype=np.int64),
                                  np.arange(p, dtype=np.int64),
                                  indexing='ij')
    res_a = np.full((p, p), coefs[-1] % p, dtype=np.int64)
    res_b = np.zeros((p, p), dtype=np.int64)
    for k in range(n - 2, -1, -1):
        new_a = (res_a * a_grid - res_b * b_grid) % p
        new_b = (res_a * b_grid + res_b * a_grid) % p
        res_a = (new_a + coefs[k]) % p
        res_b = new_b
    zero_mask = (res_a == 0) & (res_b == 0)
    ia, ib = np.where(zero_mask)
    return list(zip(ia.tolist(), ib.tolist()))


def _fp2_pow(x_a, x_b, n, p):
    """(x_a + x_b·i)^n in F_{p²}, returns (a, b). n can be negative."""
    if n == 0:
        return (1, 0)
    if n < 0:
        # Inverse first
        nrm = (x_a * x_a + x_b * x_b) % p
        ni = mod_inv(nrm, p)
        x_a, x_b = (x_a * ni) % p, (-x_b * ni) % p
        n = -n
    r_a, r_b = 1, 0
    b_a, b_b = x_a, x_b
    while n > 0:
        if n & 1:
            r_a, r_b = ((r_a * b_a - r_b * b_b) % p,
                       (r_a * b_b + r_b * b_a) % p)
        b_a, b_b = ((b_a * b_a - b_b * b_b) % p,
                   (2 * b_a * b_b) % p)
        n >>= 1
    return (r_a, r_b)


def _eval_sympy_poly_in_fp2(poly_expr, var_vals, i_sym, p):
    """Evaluate a sympy polynomial in (var1, var2, ..., i_sym) at numeric
    F_{p²} values given in `var_vals` (dict {sympy_var: (a, b)}) and i_sym
    = (0, 1). Returns (a, b)."""
    from sympy import Poly as SP
    from sympy.polys.polyerrors import GeneratorsNeeded
    sym_list = list(var_vals.keys()) + [i_sym]
    try:
        pol = SP(poly_expr, *sym_list)
    except Exception as ex:
        raise ValueError(f"_eval_sympy_poly_in_fp2 failed to convert to Poly: {ex}")
    result_a, result_b = 0, 0
    for monom, coef in pol.terms():
        # Reduce coefficient: sympy gives a Rational; reduce as a/b mod p.
        from sympy import Rational
        if hasattr(coef, 'p') and hasattr(coef, 'q'):
            num = int(coef.p) % p
            den = int(coef.q) % p
            c_a = (num * mod_inv(den, p)) % p
        else:
            c_a = int(coef) % p
        c_b = 0
        # Multiply by each var^exp
        for s, e in zip(sym_list, monom):
            if e == 0:
                continue
            if s is i_sym:
                # i^e
                r = e % 4
                if r == 0:
                    pass  # ×1
                elif r == 1:
                    c_a, c_b = (-c_b) % p, c_a
                elif r == 2:
                    c_a, c_b = (-c_a) % p, (-c_b) % p
                else:
                    c_a, c_b = c_b, (-c_a) % p
            else:
                val_a, val_b = var_vals[s]
                pow_a, pow_b = _fp2_pow(val_a, val_b, e, p)
                c_a, c_b = ((c_a * pow_a - c_b * pow_b) % p,
                           (c_a * pow_b + c_b * pow_a) % p)
        result_a = (result_a + c_a) % p
        result_b = (result_b + c_b) % p
    return (result_a, result_b)


def _eval_rational_in_fp2(expr, var_vals, i_sym, p):
    """Evaluate a sympy rational expression in (vars + i_sym) at F_{p²} values.
    Splits into num/den polynomial parts, evaluates each, then divides in F_{p²}."""
    from sympy import together, fraction, expand
    expr2 = together(expand(expr))
    num, den = fraction(expr2)
    num_a, num_b = _eval_sympy_poly_in_fp2(num, var_vals, i_sym, p)
    den_a, den_b = _eval_sympy_poly_in_fp2(den, var_vals, i_sym, p)
    if den_a == 0 and den_b == 0:
        raise ValueError("denominator zero")
    nrm = (den_a * den_a + den_b * den_b) % p
    ni = mod_inv(nrm, p)
    den_inv_a = (den_a * ni) % p
    den_inv_b = (-den_b * ni) % p
    return ((num_a * den_inv_a - num_b * den_inv_b) % p,
            (num_a * den_inv_b + num_b * den_inv_a) % p)


def solve_fp2_via_msolve(p, alpha, beta, pc, p_pt):
    """Solve the gluing divisor map over F_{p²} using msolve."""
    K = compute_constants(alpha, beta, p)
    h_coeffs = compute_h_coeffs(alpha, beta, p)
    built = _build_reduced_system(p, alpha, beta, pc, p_pt, K, h_coeffs, use_i=True)
    if built is None:
        return [{"error": "no closed-form for u_0 or v_0"}]
    u0, u1, v0, v1, i_sym = built["vars"]
    u0_expr = built["u0_expr"]
    v0_expr = built["v0_expr"]

    # Build msolve input
    poly_strs = []
    for e in built["polys"]:
        p_obj = Poly(expand(e), i_sym, u1, v1, domain=GF(p))
        s = str(p_obj.as_expr()).replace('**', '^')
        poly_strs.append(s)
    polys_str = ",\n".join(poly_strs)

    parsed = call_msolve(polys_str, ["i", "u1", "v1"], p)
    if parsed.get("empty"):
        return [{"error": "msolve: empty solution set"}]
    if "w" not in parsed:
        return [{"error": f"msolve: unexpected output: {parsed}"}]

    # Find roots of w(A) in F_{p²}.
    w_deg, w_coefs = parsed["w"]
    a_roots = fp2_roots_numpy(w_coefs, p)
    if not a_roots:
        return [{"error": "w(A) has no F_{p²} roots"}]

    # For each root A_0 (in F_{p²}), evaluate i, u_1, v_1.
    d_deg, d_coefs = parsed["denom"]

    # msolve's rational parametrization for our system has denominator = 1
    # (parsed["denom"] = (0, [1])), so each variable x_k = q_k(A) directly,
    # no division. (When the denominator is non-trivial, the formula is
    # x_k = q_k(A) / w'(A).)
    use_division = not (d_deg == 0 and d_coefs == [1])
    if use_division:
        w_prime_coefs = [(k * w_coefs[k]) % p for k in range(1, w_deg + 1)]

    solutions = []
    for (a0, b0) in a_roots:
        def ev(coefs):
            return f2_poly_eval(coefs, a0, b0, p)
        if use_division:
            wp_a, wp_b = ev(w_prime_coefs)
            if wp_a == 0 and wp_b == 0:
                continue
            nrm = (wp_a * wp_a + wp_b * wp_b) % p
            nrm_inv = mod_inv(nrm, p)
            wp_inv_a = (wp_a * nrm_inv) % p
            wp_inv_b = (-wp_b * nrm_inv) % p
            def transform(qa, qb):
                return ((qa * wp_inv_a - qb * wp_inv_b) % p,
                        (qa * wp_inv_b + qb * wp_inv_a) % p)
        else:
            def transform(qa, qb): return (qa, qb)
        ia, ib = transform(*ev(parsed["vars"]["i"][1]))
        u1a, u1b = transform(*ev(parsed["vars"]["u1"][1]))
        v1a, v1b = transform(*ev(parsed["vars"]["v1"][1]))
        # Verify i² + 1 = 0
        i_sq_a = (ia * ia - ib * ib) % p
        i_sq_b = (2 * ia * ib) % p
        if not (i_sq_a == p - 1 and i_sq_b == 0):
            continue  # not a valid i value in F_{p²}
        # Skip the degenerate u_1=0 case (closed-form has 1/u_1²).
        if u1a == 0 and u1b == 0:
            continue
        # Compute u_0, v_0 via closed-form
        u0_sym = u0_expr.subs([
            (u1, u1a + u1b * i_sym),
            (v1, v1a + v1b * i_sym),
        ])
        v0_sym = v0_expr.subs([
            (u1, u1a + u1b * i_sym),
            (v1, v1a + v1b * i_sym),
        ])
        # Evaluate u0_expr, v0_expr at (u_1 = u1a + u1b·i, v_1 = v1a + v1b·i).
        # These are rational in (u_1, v_1, i_sym); evaluate numerically in F_{p²}.
        try:
            u0_a, u0_b = _eval_rational_in_fp2(
                u0_expr, {u1: (u1a, u1b), v1: (v1a, v1b)}, i_sym, p,
            )
            v0_a, v0_b = _eval_rational_in_fp2(
                v0_expr, {u1: (u1a, u1b), v1: (v1a, v1b)}, i_sym, p,
            )
        except Exception as ex:
            sys.stderr.write(f"  u0/v0 eval failed: {ex}\n")
            continue
        solutions.append({
            "u0": [u0_a, u0_b],
            "u1": [u1a, u1b],
            "v0": [v0_a, v0_b],
            "v1": [v1a, v1b],
        })
    if not solutions:
        return [{"error": "no valid F_{p²} solutions extracted (all had i² ≠ -1)"}]
    return solutions


def main():
    data = json.load(sys.stdin)
    p = data["p"]
    alpha = [F2.from_list(x, p) for x in data["alpha"]]
    beta = [F2.from_list(x, p) for x in data["beta"]]
    pc = [F2.from_list(x, p) for x in data["pc"]]
    p_pt = [F2.from_list(x, p) for x in data["p_pt"]]

    in_fp = all(x.is_in_fp() for x in alpha + beta + pc + p_pt)
    if in_fp:
        # Fall through to the sympy-Groebner F_p path which finds the
        # trivial solution our existing test relies on.
        K = compute_constants(alpha, beta, p)
        h_coeffs = compute_h_coeffs(alpha, beta, p)
        sols = _solve_fp(p, alpha, beta, pc, p_pt, K, h_coeffs)
    else:
        sols = solve_fp2_via_msolve(p, alpha, beta, pc, p_pt)
    h_coeffs = compute_h_coeffs(alpha, beta, p)
    print(json.dumps({
        "solutions": sols,
        "h_coeffs": [[c.a, c.b] for c in h_coeffs],
    }, default=str))


if __name__ == "__main__":
    main()
