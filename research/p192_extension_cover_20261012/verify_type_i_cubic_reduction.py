"""Independent replay of all 16 degree-three residual-divisor pairs."""

import hashlib
import json
from pathlib import Path

from sage.all import PolynomialRing, matrix, prod, vector
import type_i_divisor_transfer as dt


HERE = Path(__file__).resolve().parent
receipt_path = HERE / "type_i_cubic_reduction.json"
half_path = HERE / "type_i_half_fiber_probe.json"
plane_path = HERE / "type_i_plane_quartic.json"
receipt = json.loads(receipt_path.read_text())
half = json.loads(half_path.read_text())
plane = json.loads(plane_path.read_text())
assert receipt["source_sha256"] == hashlib.sha256(
    (HERE / "type_i_cubic_reduction.py").read_bytes()
).hexdigest()
assert receipt["half_fiber_receipt_sha256"] == hashlib.sha256(half_path.read_bytes()).hexdigest()
assert receipt["plane_receipt_sha256"] == hashlib.sha256(plane_path.read_bytes()).hexdigest()
assert receipt["field_prime"] == str(dt.p)
assert receipt["anchor_base_indices"] == [0, 1, 2]
assert len(receipt["classes"]) == len(half["classes"]) == 16

L, p = dt.L, dt.p
X = PolynomialRing(L, "x")
x = X.gen()
Y = PolynomialRing(X, "y")
y = Y.gen()
quartic_terms = [(tuple(map(int, key.split(","))), L(int(c)))
                 for key, c in plane["quartic_coefficients"].items()]
f = sum(Y(c * x ** e[0]) * y ** e[1] for e, c in quartic_terms)
assert f.degree() == 4
anchors = [tuple(L(int(a)) for a in item["coordinates_Fp_decimal"])
           for item in plane["known_points_verified"] if item["source"] == "base"][:3]
exponents = [tuple(e) for e in receipt["cubic_monomial_exponents"]]
assert len(exponents) == 10 and all(sum(e) == 3 for e in exponents)


def field_element(encoded):
    return L([int(c) for c in encoded])


def point(encoded):
    return tuple(field_element(a) for a in encoded)


def polynomial(encoded):
    return X([field_element(a) for a in encoded])


def quartic_value(u, v, w):
    return sum(c * u ** e[0] * v ** e[1] * w ** e[2]
               for e, c in quartic_terms)


def check_side(seed, reported):
    support = [tuple(a ** (p ** power) for a in point(encoded))
               for power in range(3) for encoded in seed]
    coordinates = []
    for source in support:
        values = dt.check_model(source)
        assert all(value != 0 for value in values)
        u, v, w = dt.Vandermonde.solve_right(vector(L, [1 / value for value in values]))
        assert w != 0 and quartic_value(u, v, w) == 0
        coordinates.append((u, v, w))
    coordinates.extend(anchors)
    assert len(coordinates) == 9 and all(c[2] != 0 for c in coordinates)

    coefficients = [field_element(c) for c in reported["cubic_coefficients"]]
    assert len(coefficients) == 10 and any(c != 0 for c in coefficients)
    assert all(c ** p == c for c in coefficients)
    rows = [[prod(c[k] ** e[k] for k in range(3)) for e in exponents]
            for c in coordinates]
    evaluation = matrix(L, rows)
    assert evaluation.rank() == 9
    assert all(sum(a * b for a, b in zip(row, coefficients)) == 0 for row in rows)

    cubic = sum(Y(coefficients[i] * x ** e[0]) * y ** e[1]
                for i, e in enumerate(exponents))
    assert cubic.degree() == 3
    resultant = f.resultant(cubic)
    known = prod((x - c[0] / c[2] for c in coordinates), X.one())
    assert resultant.degree() == 12 and known.degree() == 9
    assert known.gcd(known.derivative()).degree() == 0
    residual_x = polynomial(reported["residual_x_polynomial"])
    residual_y = polynomial(reported["residual_y_polynomial"])
    assert residual_x.is_monic() and residual_x.degree() == 3
    assert residual_y.degree() < 3
    assert residual_x.gcd(residual_x.derivative()).degree() == 0
    assert residual_x.gcd(known).degree() == 0
    quotient, remainder = resultant.quo_rem(known)
    assert remainder == 0 and quotient.monic() == residual_x
    subresultants = f.subresultants(cubic)
    linear = next(s for s in subresultants if s.degree() == 1)
    assert linear[1].gcd(residual_x).degree() == 0
    assert (linear[1] * residual_y + linear[0]) % residual_x == 0
    assert f(residual_y) % residual_x == 0
    assert cubic(residual_y) % residual_x == 0
    assert all(c ** p == c for c in residual_x.list() + residual_y.list())
    assert reported["cubic_coefficients_in_Fp"]
    assert reported["residual_coefficients_in_Fp"]
    return True


for index, (item, source) in enumerate(zip(receipt["classes"], half["classes"])):
    assert item["base_index"] == source["base_index"] == index
    assert source["unit_divisor_verified"]
    assert check_side(source["target_fiber"], item["positive_input_residual"])
    assert check_side(source["half_fiber"], item["negative_input_residual"])
    assert item["reduced_pair_positive_from"] == "negative_input_residual"
    assert item["reduced_pair_negative_from"] == "positive_input_residual"
    print(f"replayed residual pair {index + 1}/16", flush=True)

out = {
    "status": "verified_degree_three_residual_pairs_for_16_type_i_classes",
    "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    "producer_receipt_sha256": hashlib.sha256(receipt_path.read_bytes()).hexdigest(),
    "classes_replayed": 16,
    "Fp_defined_residual_divisors_replayed": 32,
    "section_identity": "For cubic sections A+T+R_A and B+T+R_B, [A-B]=[R_B-R_A].",
    "boundary": "Degree-three positive and negative effective divisors for each class; no fixed-origin normal form or general Jacobian group law.",
}
(HERE / "type_i_cubic_reduction_verification.json").write_text(
    json.dumps(out, indent=2, sort_keys=True) + "\n"
)
print(json.dumps(out, indent=2, sort_keys=True))
