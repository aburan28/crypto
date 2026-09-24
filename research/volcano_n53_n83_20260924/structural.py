"""Exact small-model and recurrence checks for K_a: y²+xy=x³+a*x²+1.

This is a structural audit. It neither enumerates isogenies at n=53/83
nor measures relation collection, regularity, or discrete logarithms.
"""
import json
from math import isqrt
from pathlib import Path


# These are irreducible over F_2; exhaustive results also verify the trace.
REDUCTION_POLYNOMIALS = {1: 0b11, 3: 0b1011, 5: 0b100101, 7: 0b10000011}


def gf_mul(x: int, y: int, degree: int, modulus: int) -> int:
    acc = 0
    while y:
        if y & 1:
            acc ^= x
        y >>= 1
        x <<= 1
        if x & (1 << degree):
            x ^= modulus
    return acc


def exhaustive_order(degree: int, a: int) -> int:
    q, modulus = 1 << degree, REDUCTION_POLYNOMIALS[degree]
    # Independent equation solving: enumerate every possible affine y.
    # For each x, compare y²+xy with x³+a*x²+1.
    square = [gf_mul(z, z, degree, modulus) for z in range(q)]
    count = 1
    for x in range(q):
        x2 = square[x]
        rhs = gf_mul(x2, x, degree, modulus) ^ (a * x2) ^ 1
        for y in range(q):
            if square[y] ^ gf_mul(x, y, degree, modulus) == rhs:
                count += 1
    return count


def recurrence(degree: int, a: int) -> tuple[int, int]:
    trace_one = -1 if a == 0 else 1
    t_prev, t_cur = 2, trace_one
    u_prev, u_cur = 0, 1
    if degree == 1:
        return t_cur, u_cur
    for _ in range(2, degree + 1):
        t_prev, t_cur = t_cur, trace_one * t_cur - 2 * t_prev
        u_prev, u_cur = u_cur, trace_one * u_cur - 2 * u_prev
    return t_cur, u_cur


def factors(n: int) -> dict[int, int]:
    result = {}
    divisor = 2
    while divisor <= isqrt(n):
        if n % divisor == 0:
            exponent = 0
            while n % divisor == 0:
                n //= divisor
                exponent += 1
            result[divisor] = exponent
        divisor += 1 if divisor == 2 else 2
    if n > 1:
        result[n] = 1
    return result


def main() -> None:
    controls = []
    for n in REDUCTION_POLYNOMIALS:
        for a in (0, 1):
            t, u = recurrence(n, a)
            exhaustive = exhaustive_order(n, a)
            assert exhaustive == (1 << n) + 1 - t, (n, a, t, exhaustive)
            assert t * t - 4 * (1 << n) == -7 * u * u
            controls.append({'n': n, 'a': a, 'count': exhaustive})
    rows = []
    for n in (53, 83):
        for a in (0, 1):
            t, u = recurrence(n, a)
            q = 1 << n
            f_pi = abs(u)
            split = factors(f_pi)
            assert t * t - 4 * q == -7 * f_pi * f_pi
            assert -7 % 4 == 1
            prod = 1
            for prime, exponent in split.items():
                assert factors(prime) == {prime: 1}
                prod *= prime ** exponent
            assert prod == f_pi
            order = q + 1 - t
            local = {}
            for ell in split:
                symbol = pow((-7) % ell, (ell - 1) // 2, ell)
                symbol = -1 if symbol == ell - 1 else symbol
                assert symbol in (-1, 1)
                # Every cyclic ell subgroup is Frobenius-stable when ell | f_pi.
                # Pointwise-rational ell-torsion would force ell | #E.
                local[str(ell)] = {'legendre_minus7': symbol,
                                   'horizontal_kernels_on_maximal_order': 1 + symbol,
                                   'descending_kernels_on_maximal_order': ell - symbol,
                                   'ell_divides_rational_group_order': order % ell == 0}
                assert order % ell != 0
            row = {'n': n, 'a': a, 'q': q, 'trace': t,
                   'curve_order': order, 'frobenius_discriminant': -7 * f_pi * f_pi,
                   'fundamental_discriminant': -7,
                   'frobenius_conductor': f_pi,
                   'frobenius_conductor_factorization': {str(k): v for k, v in split.items()},
                   'local_kernel_prediction': local,
                   'base_endomorphism_conductor': 1,
                   'candidate_nonmaximal_conductors': 'positive divisors of f_pi',
                   'relation_yield': None, 'groebner_regularity': None,
                   'isogeny_enumeration_complete': False}
            if n == 53 and a == 0:
                known = 21044858204113
                assert order % known == 0
                row['repo_n53_subgroup_order'] = known
                row['repo_n53_cofactor'] = order // known
            rows.append(row)
    out = {'status': 'exact recurrence with exhaustive small-field controls; no target-size isogeny sweep',
           'model': 'y^2 + xy = x^3 + a*x^2 + 1, a in {0,1}',
           'reference_small_field_counts': controls, 'results': rows,
           'prime_field_volcano_backend_applies': False}
    path = Path(__file__).with_name('results.json')
    path.write_text(json.dumps(out, indent=2) + '\n')
    print('Exhaustive binary-field control models:', len(controls))
    for row in rows:
        print(f"K_{row['a']}, n={row['n']}: t={row['trace']}, #E={row['curve_order']}, "
              f"f_pi={row['frobenius_conductor']} = {row['frobenius_conductor_factorization']}")


if __name__ == '__main__':
    main()
