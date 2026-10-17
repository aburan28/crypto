"""Degree-three residual-divisor pairs for 16 P-192 Type I classes."""

import hashlib
import json
import time
from pathlib import Path

from sage.all import PolynomialRing, matrix, prod, vector
import type_i_divisor_transfer as dt


HERE = Path(__file__).resolve().parent
half_path = HERE / "type_i_half_fiber_probe.json"
plane_path = HERE / "type_i_plane_quartic.json"
half = json.loads(half_path.read_text())
plane = json.loads(plane_path.read_text())
p, L = dt.p, dt.L
anchors = [tuple(L(int(x)) for x in item["coordinates_Fp_decimal"])
           for item in plane["known_points_verified"] if item["source"] == "base"][:3]
assert len(anchors) == 3 and all(c[2] != 0 for c in anchors)
monomials = [(3, 0, 0), (0, 3, 0), (0, 0, 3), (2, 1, 0), (2, 0, 1),
             (1, 2, 0), (0, 2, 1), (1, 0, 2), (0, 1, 2), (1, 1, 1)]
X = PolynomialRing(L, "x")
x = X.gen()
Y = PolynomialRing(X, "y")
y = Y.gen()
quartic = sum(Y(L(int(c)) * x ** e[0]) * y ** e[1]
              for key, c in plane["quartic_coefficients"].items()
              for e in [tuple(map(int, key.split(",")))])
assert quartic.degree() == 4


def decode(encoded):
    return tuple(L([int(c) for c in coordinate]) for coordinate in encoded)


def plane_coordinate(point):
    values = dt.check_model(point)
    assert all(value != 0 for value in values)
    return dt.Vandermonde.solve_right(vector(L, [1 / value for value in values]))


def encode(element):
    return [str(coefficient) for coefficient in element.polynomial().list()]


def encode_polynomial(poly):
    return [encode(coefficient) for coefficient in poly.list()]


def residual(seed):
    points = [tuple(a ** (p ** i) for a in decode(item))
              for i in range(3) for item in seed]
    coordinates = [plane_coordinate(point) for point in points] + anchors
    assert len(coordinates) == 9 and all(c[2] != 0 for c in coordinates)
    evaluations = matrix(L, [[prod(c[k] ** e[k] for k in range(3))
                              for e in monomials] for c in coordinates])
    assert evaluations.rank() == 9
    coefficients = list(evaluations.right_kernel().basis()[0])
    first = next(c for c in coefficients if c != 0)
    coefficients = [c / first for c in coefficients]
    cubic = sum(Y(coefficients[i] * x ** e[0]) * y ** e[1]
                for i, e in enumerate(monomials))
    assert cubic.degree() == 3
    assert all(sum(coefficients[i] * prod(c[k] ** e[k] for k in range(3))
                   for i, e in enumerate(monomials)) == 0 for c in coordinates)
    resultant = quartic.resultant(cubic)
    known = prod((x - c[0] / c[2] for c in coordinates), X.one())
    quotient, remainder = resultant.quo_rem(known)
    assert remainder == 0 and resultant.degree() == 12
    assert known.degree() == 9 and known.gcd(known.derivative()).degree() == 0
    h = quotient.monic()
    assert h.degree() == 3 and h.gcd(h.derivative()).degree() == 0
    assert h.gcd(known).degree() == 0
    subresultant = next(item for item in quartic.subresultants(cubic)
                        if item.degree() == 1)
    linear, constant = subresultant[1], subresultant[0]
    assert linear.gcd(h).degree() == 0
    y_value = (-constant * linear.inverse_mod(h)) % h
    assert y_value.degree() < 3
    assert quartic(y_value) % h == 0 and cubic(y_value) % h == 0
    return {
        "cubic_coefficients": [encode(c) for c in coefficients],
        "cubic_coefficients_in_Fp": all(c ** p == c for c in coefficients),
        "resultant_degree": int(resultant.degree()),
        "known_x_squarefree": True,
        "residual_x_degree": int(h.degree()),
        "residual_x_squarefree_disjoint": True,
        "residual_x_polynomial": encode_polynomial(h),
        "residual_y_polynomial": encode_polynomial(y_value),
        "residual_coefficients_in_Fp": all(c ** p == c for c in h.list() + y_value.list()),
        "subresultant_y_coordinate_verified": True,
    }


start = time.monotonic_ns()
classes = []
for item in half["classes"]:
    assert item["unit_divisor_verified"]
    positive_input = residual(item["target_fiber"])
    negative_input = residual(item["half_fiber"])
    classes.append({"base_index": item["base_index"],
                    "positive_input_residual": positive_input,
                    "negative_input_residual": negative_input,
                    "reduced_pair_positive_from": "negative_input_residual",
                    "reduced_pair_negative_from": "positive_input_residual"})
    print(f"degree-three residual pair {item['base_index'] + 1}/16", flush=True)

out = {
    "status": "verified_degree_three_residual_pairs_for_16_type_i_classes",
    "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    "half_fiber_receipt_sha256": hashlib.sha256(half_path.read_bytes()).hexdigest(),
    "plane_receipt_sha256": hashlib.sha256(plane_path.read_bytes()).hexdigest(),
    "field_prime": str(p),
    "anchor_base_indices": [0, 1, 2],
    "cubic_monomial_exponents": [list(e) for e in monomials],
    "classes": classes,
    "degree_three_pairs": len(classes),
    "Fp_defined_pairs": sum(all(side["cubic_coefficients_in_Fp"] and
                                side["residual_coefficients_in_Fp"]
                                for side in (c["positive_input_residual"],
                                             c["negative_input_residual"]))
                            for c in classes),
    "elapsed_ns": time.monotonic_ns() - start,
    "class_identity": "If A+T+R_A and B+T+R_B are the two cubic sections, [A-B]=[R_B-R_A].",
    "boundary": "Exact degree-three effective divisor pairs, represented by triangular x/y polynomials; a fixed-origin normal form and generic Jacobian group law remain separate.",
}
(HERE / "type_i_cubic_reduction.json").write_text(json.dumps(out, indent=2, sort_keys=True) + "\n")
print(json.dumps({k: v for k, v in out.items() if k != "classes"}, sort_keys=True))
