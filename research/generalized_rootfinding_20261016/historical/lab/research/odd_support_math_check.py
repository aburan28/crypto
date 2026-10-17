"""Independent tiny oracle for odd norm support fibers.

This uses an ordinary power recurrence and explicit domain multisets, rather
than the compact Frobenius construction. It is correctness evidence only.
"""

import argparse
import hashlib
from itertools import combinations_with_replacement
import json
from pathlib import Path
import random
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from rootfinder.fields import PrimeField, ExtensionField
from rootfinder import polynomial as P
from lab.research.odd_algebra_check import Arithmetic, formulas


def scale(poly, scalar, field):
    return P.trim(tuple(field.mul(value, scalar) for value in poly), field)


def domain_poly(domain, field):
    result = (1,)
    for element in domain:
        result = P.mul(result, (field.neg(element), 1), field)
    return result


def remainder_power_recurrence(domain_polynomial, e2, e3, s, field):
    powers = [((1,), (), ()), ((), (1,), ()), ((), (), (1,))]
    result = [(), (), ()]
    for exponent, coefficient in enumerate(domain_polynomial):
        if exponent >= 3:
            previous = powers
            next_power = tuple(P.add(P.sub(scale(previous[-1][j], s, field),
                                                 P.mul(e2, previous[-2][j], field), field),
                                    P.mul(e3, previous[-3][j], field), field)
                               for j in range(3))
            powers = powers[1:] + [next_power]
        current = powers[exponent] if exponent < 3 else powers[-1]
        for coordinate in range(3):
            result[coordinate] = P.add(result[coordinate],
                                       scale(current[coordinate], coefficient, field), field)
    return tuple(result)


def discriminant(e2, e3, s, field):
    a = Arithmetic(field)
    result = scale(P.mul(e2, e2, field), a.square(s), field)
    result = P.sub(result, scale(P.mul(P.mul(e2, e2, field), e2, field),
                                 a.scalar(4), field), field)
    result = P.sub(result, scale(e3, a.mul(a.scalar(4), a.cube(s)), field), field)
    result = P.sub(result, scale(P.mul(e3, e3, field), a.scalar(27), field), field)
    return P.add(result, scale(P.mul(e2, e3, field),
                               a.mul(a.scalar(18), s), field), field)


def multiset_oracle(domain, field):
    by_sum = {}
    for roots in combinations_with_replacement(domain, 3):
        polynomial = domain_poly(roots, field)
        source_sum = field.neg(polynomial[2])
        kind = "squarefree" if len(set(roots)) == 3 else "repeated"
        by_sum.setdefault(source_sum, {"squarefree": set(), "repeated": set()})[kind].add(polynomial)
    return by_sum


def run():
    rng = random.Random(2026101201)
    results = []
    for field in (PrimeField(3), PrimeField(5), ExtensionField(3, 2), ExtensionField(5, 2)):
        arithmetic = Arithmetic(field)
        counter = {"field_order": field.order, "characteristic": field.characteristic,
                   "degree": field.degree, "modulus": field.modulus,
                   "curve_cases": 0, "sum_fibers": 0, "parameter_evaluations": 0,
                   "squarefree_parameters": 0, "repeated_parameters": 0,
                   "gcd_degrees": [], "degree_bound_checks": 0}
        while counter["curve_cases"] < 8:
            a2, a4, z, y, offset = (rng.randrange(field.order) for _ in range(5))
            a6 = arithmetic.sub(arithmetic.square(y),
                                arithmetic.add(arithmetic.cube(z),
                                               arithmetic.mul(a2, arithmetic.square(z)),
                                               arithmetic.mul(a4, z)))
            curve = (a6, a4, a2, 1)
            derivative = (a4, arithmetic.mul(arithmetic.scalar(2), a2), arithmetic.scalar(3))
            if P.gcd_poly(curve, derivative, field) != (1,):
                continue
            counter["curve_cases"] += 1
            domain = tuple(sorted(field.add(offset, scalar) for scalar in range(field.characteristic)))
            domain_polynomial = domain_poly(domain, field)
            cube_domain = P.mul(P.mul(domain_polynomial, domain_polynomial, field),
                                domain_polynomial, field)
            independent = multiset_oracle(domain, field)
            sums = tuple(sorted(field.add(arithmetic.mul(arithmetic.scalar(3), offset), scalar)
                                for scalar in range(field.characteristic)))
            assert set(independent) == set(sums)
            for source_sum in sums:
                e2, e3, _, _, _ = formulas(field, a2, a4, z, y, source_sum, 0, 0)
                remainders = remainder_power_recurrence(domain_polynomial, e2, e3,
                                                       source_sum, field)
                assert any(remainders)
                for coordinate, remainder in enumerate(remainders):
                    if remainder:
                        assert len(remainder) - 1 <= 2 * (len(domain) - coordinate)
                    counter["degree_bound_checks"] += 1
                delta = discriminant(e2, e3, source_sum, field)
                assert len(delta) == 13
                assert delta[-1] == field.neg(field.inv(arithmetic.scalar(16)))
                common = ()
                for remainder in remainders:
                    common = P.gcd_poly(common, remainder, field)
                assert common and len(common) - 1 <= 2 * len(domain)
                counter["gcd_degrees"].append(len(common) - 1)
                frobenius = [0] * (field.order + 1)
                frobenius[1], frobenius[-1] = field.neg(1), 1
                restricted = P.gcd_poly(common, tuple(frobenius), field)
                squarefree, repeated, oracle_sf, oracle_rep = set(), set(), set(), set()
                exception_roots = []
                for parameter in range(field.order):
                    h = (field.neg(P.evaluate_poly(e3, parameter, field)),
                         P.evaluate_poly(e2, parameter, field), field.neg(source_sum), 1)
                    if not P.evaluate_poly(restricted, parameter, field):
                        squarefree.add(parameter)
                    is_repeated = P.evaluate_poly(delta, parameter, field) == 0
                    if is_repeated:
                        exception_roots.append(parameter)
                        if not P.divmod_poly(cube_domain, h, field)[1]:
                            repeated.add(parameter)
                    if h in independent[source_sum]["squarefree"]:
                        oracle_sf.add(parameter)
                    if h in independent[source_sum]["repeated"]:
                        oracle_rep.add(parameter)
                    counter["parameter_evaluations"] += 1
                assert squarefree == oracle_sf
                assert repeated == oracle_rep
                assert squarefree.isdisjoint(repeated)
                assert len(exception_roots) <= min(12, field.order)
                counter["sum_fibers"] += 1
                counter["squarefree_parameters"] += len(squarefree)
                counter["repeated_parameters"] += len(repeated)
        results.append(counter)
    return {"status": "PASS", "seed": 2026101201,
            "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            "scope": "Independent explicit multiset oracle for monic support specialization, degree bounds, gcd roots, and degree12 repeated branch; no speed or novelty claim.",
            "fields": results}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = run()
    encoded = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(encoded)
    print(encoded, end="")


if __name__ == "__main__":
    main()
