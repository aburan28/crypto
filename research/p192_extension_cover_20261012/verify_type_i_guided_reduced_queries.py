"""Independent replay of frozen guided scalars and degree-three query divisors."""

import hashlib
import json
from pathlib import Path

from sage.all import EllipticCurve, PolynomialRing, matrix, prod, vector
import type_i_divisor_transfer as dt
import verify_type_i_quotient_base as vb


HERE = Path(__file__).resolve().parent
receipt_path = HERE / "type_i_guided_reduced_queries.json"
query_path = HERE / "type_i_geom_character_collector.json"
query_verify_path = HERE / "type_i_geom_character_verification.json"
plane_path = HERE / "type_i_plane_quartic.json"
receipt = json.loads(receipt_path.read_text())
query = json.loads(query_path.read_text())
plane = json.loads(plane_path.read_text())
assert receipt["source_sha256"] == hashlib.sha256(
    (HERE / "type_i_guided_reduced_queries.py").read_bytes()
).hexdigest()
assert receipt["query_receipt_sha256"] == hashlib.sha256(query_path.read_bytes()).hexdigest()
assert receipt["query_verification_sha256"] == hashlib.sha256(
    query_verify_path.read_bytes()
).hexdigest()
assert receipt["plane_receipt_sha256"] == hashlib.sha256(plane_path.read_bytes()).hexdigest()
assert len(receipt["results"]) == len(query["accepted_scalars_decimal"]) == 64

p, r, L = dt.p, dt.r, dt.L
E = EllipticCurve(L, [L(dt.gm.a), L(dt.gm.b)])
G = E(L(dt.gm.gx), L(dt.gm.gy))
O = E(0)
assert r * G == O
anchors = [tuple(L(int(v)) for v in item["coordinates_Fp_decimal"])
           for item in plane["known_points_verified"] if item["source"] == "base"][:3]
monomials = [(3, 0, 0), (0, 3, 0), (0, 0, 3), (2, 1, 0), (2, 0, 1),
             (1, 2, 0), (0, 2, 1), (1, 0, 2), (0, 1, 2), (1, 1, 1)]
X = PolynomialRing(L, "x")
x = X.gen()
Y = PolynomialRing(X, "y")
y = Y.gen()
quartic = sum(Y(L(int(c)) * x ** e[0]) * y ** e[1]
              for key, c in plane["quartic_coefficients"].items()
              for e in [tuple(map(int, key.split(",")))])


def decode(encoded):
    return L([int(c) for c in encoded])


def poly(encoded):
    return X([decode(c) for c in encoded])


def check_side(fiber, side):
    coordinates = []
    for power in range(3):
        for source in fiber:
            point = tuple(value ** (p ** power) for value in source)
            ys = dt.check_model(point)
            assert all(value != 0 for value in ys)
            coordinate = dt.Vandermonde.solve_right(vector(L, [1 / value for value in ys]))
            assert coordinate[2] != 0
            coordinates.append(tuple(coordinate))
    coordinates.extend(anchors)
    assert len(coordinates) == 9
    rows = [[prod(c[k] ** e[k] for k in range(3)) for e in monomials]
            for c in coordinates]
    assert matrix(L, rows).rank() == 9
    coefficients = [decode(c) for c in side["cubic_coefficients"]]
    assert len(coefficients) == 10 and all(c ** p == c for c in coefficients)
    assert all(sum(a * b for a, b in zip(row, coefficients)) == 0 for row in rows)
    cubic = sum(Y(coefficients[i] * x ** e[0]) * y ** e[1]
                for i, e in enumerate(monomials))
    h = poly(side["residual_x_polynomial"])
    y_value = poly(side["residual_y_polynomial"])
    known = prod((x - c[0] / c[2] for c in coordinates), X.one())
    resultant = quartic.resultant(cubic)
    quotient, remainder = resultant.quo_rem(known)
    assert remainder == 0 and resultant.degree() == 12
    assert quotient.monic() == h and h.degree() == 3 and h.is_monic()
    assert h.gcd(h.derivative()).degree() == 0
    assert known.gcd(known.derivative()).degree() == 0 and h.gcd(known).degree() == 0
    linear = next(s for s in quartic.subresultants(cubic) if s.degree() == 1)
    assert linear[1].gcd(h).degree() == 0
    assert (linear[1] * y_value + linear[0]) % h == 0
    assert quartic(y_value) % h == 0 and cubic(y_value) % h == 0
    assert all(c ** p == c for c in h.list() + y_value.list())


base = [item[0] for item in vb.base]
assert len(base) == 16
for index, item in enumerate(receipt["results"]):
    scalar = int(query["accepted_scalars_decimal"][index])
    assert item["query_index"] == index and int(item["scalar"]) == scalar
    q = -2 * scalar * G
    h = -scalar * G
    assert 2 * h == q and item["target_quotient_point"] == str(q)
    check_side(dt.fiber(q), item["positive_input_residual"])
    check_side(dt.fiber(h), item["negative_input_residual"])
    matches = []
    for i in range(16):
        for j in range(i, 16):
            for k in range(j, 16):
                if base[i] + base[j] + base[k] == q:
                    matches.append([i, j, k])
    assert matches == item["three_sum_matches"]
    print(f"replayed guided reduced query {index + 1}/64", flush=True)

assert receipt["degree_three_Fp_pairs"] == 64
assert receipt["three_sum_matches"] == query["exact_three_sum_matches"] == 0
out = {
    "status": "verified_guided_degree_three_query_reduction",
    "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    "producer_receipt_sha256": hashlib.sha256(receipt_path.read_bytes()).hexdigest(),
    "accepted_scalars_replayed": 64,
    "degree_three_Fp_pairs_replayed": 64,
    "unordered_triples_checked_per_query": 816,
    "three_sum_matches": 0,
    "boundary": "Frozen conditional known-scalar query law and degree-three reduction; no useful relation rows from the fixed base.",
}
(HERE / "type_i_guided_reduced_verification.json").write_text(
    json.dumps(out, indent=2, sort_keys=True) + "\n"
)
print(json.dumps(out, indent=2, sort_keys=True))
