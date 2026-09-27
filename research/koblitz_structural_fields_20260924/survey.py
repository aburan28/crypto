"""Exact Koblitz order arithmetic at field exponents 37, 53 and 83."""
import argparse
import json
from math import isqrt, prod
from pathlib import Path

EXPONENTS = (37, 53, 83)


def factor_integer(value):
    factors = []
    divisor = 2
    while divisor * divisor <= value:
        exponent = 0
        while value % divisor == 0:
            value //= divisor
            exponent += 1
        if exponent:
            factors.append((divisor, exponent))
        divisor = 3 if divisor == 2 else divisor + 2
    if value > 1:
        factors.append((value, 1))
    return factors


def is_prime(value):
    if value < 2:
        return False
    return all(value % d for d in range(2, isqrt(value) + 1))


def row(n, a):
    mu = 2*a-1
    q = 1 << n
    A, B = 1, 0
    for _ in range(n):
        A, B = -2*B, A+mu*B
    trace = 2*A+mu*B
    previous, current = 2, mu
    for _ in range(2, n+1):
        previous, current = current, mu*current-2*previous
    conductor = abs(B)
    factors = factor_integer(conductor)
    assert current == trace
    assert A*A+mu*A*B+2*B*B == q
    assert trace*trace-4*q == -7*conductor*conductor
    assert trace*trace <= 4*q and trace % 2 == 1
    assert conductor % 2 == 1
    assert all(is_prime(p) for p, _ in factors)
    assert prod(p**e for p, e in factors) == conductor
    return dict(field_exponent=n, field_size=q, a2=a, a6=1, mu=mu,
                model_id=f'KOBLITZ-F2N{n}-A{a}',
                trace=trace, point_count=q+1-trace, tau_power_A=A, tau_power_B=B,
                endomorphism_field_discriminant=-7,
                source_endomorphism_order_conductor=1,
                frobenius_order_conductor=conductor,
                frobenius_order_discriminant=-7*conductor*conductor,
                factorization=[dict(prime=p, exponent=e) for p, e in factors],
                rational_prime_degree_descents=[dict(
                    isogeny_degree=p, depth=e,
                    first_target_order_conductor=p,
                    first_target_order_discriminant=-7*p*p,
                    first_target_order_id=f'ORD-DKminus7-F{p}',
                    map_constructed=False) for p, e in factors],
                finite_field_basis=None, neighbor_j_invariants=None,
                solver_measurements=None, degree_of_regularity=None,
                end_to_end_operations=None, speedup=None)


def survey():
    rows = [row(n, a) for n in EXPONENTS for a in (0, 1)]
    for n in EXPONENTS:
        a, b = [r for r in rows if r['field_exponent'] == n]
        assert a['trace'] == -b['trace']
        assert a['frobenius_order_conductor'] == b['frobenius_order_conductor']
        assert a['point_count']+b['point_count'] == 2*((1 << n)+1)
    return dict(kind='exact structural survey; not solver evidence', rows=rows)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path('results.json'))
    args = parser.parse_args()
    result = survey()
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    for r in result['rows']:
        print(r['model_id'], r['factorization'])
