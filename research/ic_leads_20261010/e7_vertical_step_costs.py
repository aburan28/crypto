#!/usr/bin/env python3
"""E7 (conductor gap, angle V9): cost table for one vertical isogeny step per (n, ell).

Pure arithmetic, no Sage.  For every a=0 Koblitz rung n and every prime ell | f = cond(Z[pi]):
  d        = degree of the field where E[ell] kernel points live  (= ord of pi/2 mod ell, with pi acting
             as a scalar on E[ell] when ell is inert, or on the eigenlines when split)
  kernel route   : build F_{q^d}, find a point of order ell, Velu/sqrt-elu over F_{q^d}, descend j.
                   F_q-mults ~ (ell-1)/2 kernel multiples x M(d) for Velu, or ~ sqrt(ell) x M(d) for sqrt-elu,
                   where M(d) ~ d^1.585 F_q-mults is one F_{q^d} multiplication (Karatsuba) and the
                   point-finding cost is log2(#E(F_{q^d})) x M(d) doublings.
  modular route  : Phi_ell has ~ ell^2 coefficients; evaluating Phi_ell(X, j) costs ~ ell^2 F_q-mults, root
                   finding of a degree ell+1 polynomial over F_q costs ~ ell * n F_q-mults (gcd with X^q - X).
  transfer       : evaluating the isogeny on (G, Q): ~ 2 * (ell-1)/2 F_q-mults for Velu with rational kernel,
                   or sqrt-elu ~ 2 * sqrt(ell) * log(ell).
All counts are order-of-magnitude models in F_q-multiplications (log2 reported); they are what
"reachable" means operationally and are compared to rho's 2^(0.5 log2 r - 0.5 log2(4n)).
Output: e7_vertical_step_costs.json and a Markdown table on stdout.
"""
import json, math, sys
from sympy import factorint, legendre_symbol, n_order, isprime

LADDER = [23, 37, 41, 53, 61, 71, 73, 83, 97, 131]

def lucas(n):
    V = [2, -1]; U = [0, 1]
    for k in range(1, n):
        V.append(-V[k] - 2 * V[k - 1]); U.append(-U[k] - 2 * U[k - 1])
    return V[n], abs(U[n])

def chi(ell):
    return 1 if ell == 2 else int(legendre_symbol((-7) % ell, ell))

def M(d):  # one F_{q^d} multiplication in F_q multiplications
    return max(1.0, d ** 1.585)

def l2(x):
    return round(math.log2(x), 2) if x > 0 else None

def row(n, ell, t, f):
    q = 2 ** n
    order = q + 1 - t
    r = order // 4
    rho_bits = 0.5 * math.log2(math.pi * r / (4 * n))
    c = chi(ell)
    # pi acts on E[ell]; with t^2 - 4q = -7 f^2 and ell | f, pi = t/2 mod ell is a scalar on E[ell].
    a = (t * pow(2, -1, ell)) % ell
    d = int(n_order(a, ell)) if a % ell else None
    curves_below = ell - c
    out = {"n": n, "ell": int(ell), "ell_bits": l2(ell), "chi_minus7": c, "pi_mod_ell": int(a),
           "kernel_field_degree": d, "curves_one_level_down": int(curves_below),
           "log2_r": l2(r), "rho_signed_frobenius_bits": round(rho_bits, 2)}
    if d is None:
        return out
    # kernel route
    point_find = d * math.log2(order) * M(d) * 4          # scalar mult by cofactor, ~4 F_{q^d} mults per bit
    velu = ((ell - 1) / 2) * M(d) * 12                      # ~12 F_{q^d} mults per kernel multiple (x, y, Velu sums)
    sqrtelu = math.sqrt(ell) * math.log2(ell + 1) * M(d) * 12
    descend_j = d * d * n                                   # minimal polynomial of j over F_q by linear algebra
    out["kernel_route_velu_log2_Fq_mults"] = l2(point_find + velu + descend_j)
    out["kernel_route_sqrtelu_log2_Fq_mults"] = l2(point_find + sqrtelu + descend_j)
    out["kernel_route_feasible_d_le_1000"] = d <= 1000
    # modular route
    phi_eval = ell ** 2
    root_find = (ell + 1) * n * 2
    out["modular_route_log2_Fq_mults"] = l2(phi_eval + root_find)
    out["modular_route_feasible_ell_le_1e4"] = ell <= 10_000
    out["modular_polynomial_coefficients_log2"] = l2(ell ** 2)
    # transfer (the isogeny evaluated on two points, kernel rational over F_{q^d})
    transfer_velu = 2 * ((ell - 1) / 2) * M(d) * 6
    transfer_sqrtelu = 2 * math.sqrt(ell) * math.log2(ell + 1) * M(d) * 6
    out["transfer_two_points_velu_log2"] = l2(transfer_velu)
    out["transfer_two_points_sqrtelu_log2"] = l2(transfer_sqrtelu)
    best = min(x for x in [out["kernel_route_velu_log2_Fq_mults"], out["kernel_route_sqrtelu_log2_Fq_mults"],
                           out["modular_route_log2_Fq_mults"]] if x is not None)
    out["cheapest_descent_log2"] = best
    out["descent_cheaper_than_rho"] = best < rho_bits
    return out

def main():
    rows = []
    for n in LADDER:
        t, f = lucas(n)
        assert t * t - 4 * 2 ** n == -7 * f * f
        for ell in sorted(factorint(f)):
            rows.append(row(n, ell, t, f))
    json.dump(rows, open("e7_vertical_step_costs.json", "w"), indent=1)
    print("| n | ell | bits | chi | d (kernel field) | Velu kernel route | sqrt-elu kernel route | modular route | transfer (2 pts) | cheapest | rho (A=2n) | descent < rho |")
    print("|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|:--|")
    for o in rows:
        if o["kernel_field_degree"] is None:
            continue
        print(f"| {o['n']} | {o['ell']} | {o['ell_bits']} | {o['chi_minus7']:+d} | {o['kernel_field_degree']} | "
              f"2^{o['kernel_route_velu_log2_Fq_mults']} | 2^{o['kernel_route_sqrtelu_log2_Fq_mults']} | "
              f"2^{o['modular_route_log2_Fq_mults']} | 2^{o['transfer_two_points_sqrtelu_log2']} | "
              f"2^{o['cheapest_descent_log2']} | 2^{o['rho_signed_frobenius_bits']} | {'yes' if o['descent_cheaper_than_rho'] else 'no'} |")

if __name__ == "__main__":
    main()
