"""Replay the frozen guided queries as prime-field degree-three divisor pairs."""

import hashlib
import json
import time
from pathlib import Path

from sage.all import EllipticCurve, GF, PolynomialRing, matrix, prod, vector
import type_i_divisor_transfer as dt
import verify_type_i_quotient_base as vb


HERE = Path(__file__).resolve().parent
query_path = HERE / "type_i_geom_character_collector.json"
query_verify_path = HERE / "type_i_geom_character_verification.json"
plane_path = HERE / "type_i_plane_quartic.json"
query = json.loads(query_path.read_text())
query_verify = json.loads(query_verify_path.read_text())
plane = json.loads(plane_path.read_text())
assert query_verify["status"] == "verified_type_i_f3_character_query_probe"
assert query_verify["producer_sha256"] == hashlib.sha256(
    (HERE / "type_i_geom_character_collector.py").read_bytes()
).hexdigest()
assert len(query["accepted_scalars_decimal"]) == query["accepted_queries"] == 64

p, r, L = dt.p, dt.r, dt.L
E = EllipticCurve(L, [L(dt.gm.a), L(dt.gm.b)])
G = E(L(dt.gm.gx), L(dt.gm.gy))
O = E(0)
assert r * G == O and r.is_prime()
assert len(vb.base) == 16
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


def encode(element):
    return [str(c) for c in element.polynomial().list()]


def encode_poly(poly):
    return [encode(c) for c in poly.list()]


def residual(fiber):
    points = [tuple(value ** (p ** power) for value in point)
              for power in range(3) for point in fiber]
    coordinates = []
    for point in points:
        ys = dt.check_model(point)
        assert all(value != 0 for value in ys)
        coordinates.append(dt.Vandermonde.solve_right(vector(L, [1 / value for value in ys])))
    coordinates.extend(anchors)
    assert len(coordinates) == 9 and all(c[2] != 0 for c in coordinates)
    rows = [[prod(c[k] ** e[k] for k in range(3)) for e in monomials]
            for c in coordinates]
    evaluations = matrix(L, rows)
    assert evaluations.rank() == 9
    coefficients = list(evaluations.right_kernel().basis()[0])
    first = next(c for c in coefficients if c != 0)
    coefficients = [c / first for c in coefficients]
    cubic = sum(Y(coefficients[i] * x ** e[0]) * y ** e[1]
                for i, e in enumerate(monomials))
    known = prod((x - c[0] / c[2] for c in coordinates), X.one())
    resultant = quartic.resultant(cubic)
    quotient, remainder = resultant.quo_rem(known)
    assert resultant.degree() == 12 and remainder == 0
    assert known.degree() == 9 and known.gcd(known.derivative()).degree() == 0
    h = quotient.monic()
    assert h.degree() == 3 and h.gcd(h.derivative()).degree() == 0
    assert h.gcd(known).degree() == 0
    linear = next(s for s in quartic.subresultants(cubic) if s.degree() == 1)
    assert linear[1].gcd(h).degree() == 0
    y_value = (-linear[0] * linear[1].inverse_mod(h)) % h
    assert y_value.degree() < 3
    assert quartic(y_value) % h == 0 and cubic(y_value) % h == 0
    assert all(c ** p == c for c in coefficients + h.list() + y_value.list())
    return {"cubic_coefficients": [encode(c) for c in coefficients],
            "residual_x_polynomial": encode_poly(h),
            "residual_y_polynomial": encode_poly(y_value)}


pairs = {}
for i in range(16):
    for j in range(i, 16):
        pairs.setdefault(str(vb.base[i][0] + vb.base[j][0]), []).append((i, j))
assert sum(map(len, pairs.values())) == 136

start = time.monotonic_ns()
results = []
for index, decimal in enumerate(query["accepted_scalars_decimal"]):
    scalar = int(decimal)
    q = -2 * scalar * G
    h = -scalar * G
    assert 2 * h == q and r * q == O
    positive = residual(dt.fiber(q))
    negative = residual(dt.fiber(h))
    matches = []
    for k in range(16):
        for i, j in pairs.get(str(q - vb.base[k][0]), []):
            assert all(vb.base[i][m] + vb.base[j][m] + vb.base[k][m] == q
                       for m in range(3))
            matches.append((i, j, k))
    results.append({"query_index": index, "scalar": decimal,
                    "target_quotient_point": str(q),
                    "positive_input_residual": positive,
                    "negative_input_residual": negative,
                    "three_sum_matches": matches})
    print(f"guided reduced query {index + 1}/64", flush=True)

assert sum(len(item["three_sum_matches"]) for item in results) == query["exact_three_sum_matches"]
out = {
    "status": "verified_guided_degree_three_query_reduction",
    "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    "query_receipt_sha256": hashlib.sha256(query_path.read_bytes()).hexdigest(),
    "query_verification_sha256": hashlib.sha256(query_verify_path.read_bytes()).hexdigest(),
    "plane_receipt_sha256": hashlib.sha256(plane_path.read_bytes()).hexdigest(),
    "accepted_queries": len(results),
    "attempted_scalars": query["attempted_scalars"],
    "degree_three_Fp_pairs": len(results),
    "three_sum_matches": sum(len(item["three_sum_matches"]) for item in results),
    "quotient_incidence_rank": 0,
    "results": results,
    "elapsed_ns": time.monotonic_ns() - start,
    "class_equation": "query target = a*N with N=(-2G,-2G,-2G); half h=-aG, q=-2aG, and [D_q]=[R_h-R_q]",
    "boundary": "Frozen known-scalar guided law and exact degree-three pairs for its accepted queries; fixed 16-class three-sum support and zero rank remain unchanged.",
}
(HERE / "type_i_guided_reduced_queries.json").write_text(json.dumps(out, indent=2, sort_keys=True) + "\n")
print(json.dumps({k: v for k, v in out.items() if k != "results"}, sort_keys=True))
