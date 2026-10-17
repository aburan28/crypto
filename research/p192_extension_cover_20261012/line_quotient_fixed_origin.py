"""Reduce 141 exact P-192 line-derived quotient classes to degree three.

The two-cubic identity is G_cut=B+T+U and F_cut=A+U+D, whence
[D-T]=[B-A]. A and B are the frozen quotient and half-quotient input
residual cubics; T is one fixed rational degree-three origin divisor.
"""

import hashlib
import json
import time
from pathlib import Path

from sage.all import GF, PolynomialRing, matrix

HERE = Path(__file__).resolve().parent
sha = lambda path: hashlib.sha256(Path(path).read_bytes()).hexdigest()
source_path = HERE / "type_i_line_quotient_cubic_probe.json"
plane_path = HERE / "type_i_plane_quartic.json"
source = json.loads(source_path.read_text())
plane = json.loads(plane_path.read_text())
assert source["status"] == "bounded_line_quotient_cubic_constructor_probe"
assert source["successes"] == 141 and source["failures"] == 0
p = int(source["field_prime"])
K = GF(p)
R = PolynomialRing(K, "x")
x = R.gen()
S = PolynomialRing(R, "y")
y = S.gen()
exponents = [(3, 0), (0, 3), (0, 0), (2, 1), (2, 0),
             (1, 2), (0, 2), (1, 0), (0, 1), (1, 1)]
f = sum(S(K(int(c)) * x ** a) * y ** b
        for label, c in plane["quartic_coefficients"].items()
        for a, b, _ in [tuple(map(int, label.split(",")))])
assert f.degree() == 4


def poly(encoded):
    assert all(len(c) == 1 for c in encoded)
    return R([K(int(c[0])) for c in encoded])


def encode_poly(value):
    return [[str(int(c))] for c in value.list()]


def section(coefficients):
    return sum(S(coefficients[i] * x ** a) * y ** b
               for i, (a, b) in enumerate(exponents))


def evaluation_rows(h, y_value):
    assert h.is_monic() and y_value.degree() < h.degree()
    columns = [(x ** a * y_value ** b) % h for a, b in exponents]
    return [[c[j] for c in columns] for j in range(h.degree())]


def normalized(coefficients):
    first = next(c for c in coefficients if c != 0)
    return [c / first for c in coefficients]


def divide_section(cubic, known, residual_degree):
    resultant = f.resultant(cubic)
    if resultant == 0 or resultant.degree() != 12:
        return None
    quotient, remainder = resultant.quo_rem(known)
    if remainder != 0 or quotient.degree() != residual_degree:
        return None
    h = quotient.monic()
    if h.gcd(h.derivative()).degree() != 0 or h.gcd(known).degree() != 0:
        return None
    linear = next((item for item in f.subresultants(cubic) if item.degree() == 1), None)
    if linear is None or linear[1].gcd(h).degree() != 0:
        return None
    y_value = (-linear[0] * linear[1].inverse_mod(h)) % h
    assert f(y_value) % h == 0 and cubic(y_value) % h == 0
    return h, y_value


anchors = [item["coordinates_Fp_decimal"]
           for item in plane["known_points_verified"] if item["source"] == "base"][:3]
assert len(anchors) == 3
origin_coordinates = [(K(int(c[0])) / K(int(c[2])),
                       K(int(c[1])) / K(int(c[2]))) for c in anchors]
assert len(set(origin_coordinates)) == 3
origin_h = R.prod(x - xx for xx, _ in origin_coordinates)
origin_y = R.lagrange_polynomial(origin_coordinates)
assert origin_h.degree() == 3 and f(origin_y) % origin_h == 0

trials = ([1, 1, 1, 1], [1, 2, 3, 5], [1, 3, 7, 11],
          [2, 5, 11, 17], [1, 0, 0, 0], [0, 1, 0, 0])
start = time.monotonic_ns()
representatives = []
for item in source["rows"]:
    index = item["index"]
    assert item["status"] == "verified_prime_field_cubic_pair"
    A, B = item["positive_input_residual"], item["negative_input_residual"]
    h_a, y_a = poly(A["residual_x_polynomial"]), poly(A["residual_y_polynomial"])
    h_b, y_b = poly(B["residual_x_polynomial"]), poly(B["residual_y_polynomial"])
    assert all(h.degree() == 3 and f(v) % h == 0
               for h, v in ((h_a, y_a), (h_b, y_b)))
    assert h_b.gcd(origin_h).degree() == 0
    denominator_rows = matrix(K, evaluation_rows(h_b, y_b) +
                              evaluation_rows(origin_h, origin_y))
    assert denominator_rows.rank() == 6
    basis = [list(v) for v in denominator_rows.right_kernel().basis()]
    selected = None
    for trial in trials:
        coefficients = normalized([sum(K(trial[j]) * basis[j][i] for j in range(4))
                                   for i in range(10)])
        g = section(coefficients)
        intermediate = divide_section(g, h_b * origin_h, 6)
        if intermediate is None:
            continue
        h_u, y_u = intermediate
        if h_u.gcd(h_a).degree() != 0:
            continue
        numerator_rows = matrix(K, evaluation_rows(h_a, y_a) +
                                evaluation_rows(h_u, y_u))
        if numerator_rows.rank() != 9:
            continue
        numerator_coefficients = normalized(list(numerator_rows.right_kernel().basis()[0]))
        final = divide_section(section(numerator_coefficients), h_a * h_u, 3)
        if final is None:
            continue
        h_d, y_d = final
        selected = {
            "index": index, "cell": item["cell"], "quotient": item["quotient"],
            "denominator_trial": trial,
            "denominator_cubic_coefficients": [str(int(c)) for c in coefficients],
            "intermediate_x_polynomial": encode_poly(h_u),
            "intermediate_y_polynomial": encode_poly(y_u),
            "numerator_cubic_coefficients": [str(int(c)) for c in numerator_coefficients],
            "representative_x_polynomial": encode_poly(h_d),
            "representative_y_polynomial": encode_poly(y_d),
            "status": "verified_fixed_origin_degree_three",
        }
        break
    if selected is None:
        selected = {"index": index, "cell": item["cell"], "quotient": item["quotient"],
                    "status": "bounded_fixed_origin_trial_exhausted"}
    representatives.append(selected)
    if (index + 1) % 16 == 0:
        print(f"LINE_FIXED_ORIGIN_PROGRESS {index + 1}/141", flush=True)

successes = sum(row["status"] == "verified_fixed_origin_degree_three"
                for row in representatives)
out = {
    "status": "bounded_line_quotient_fixed_origin_reduction",
    "field_prime": str(p), "subgroup_order": source["subgroup_order"],
    "input_classes": 141, "reduced_classes": successes,
    "unreduced_classes": 141 - successes,
    "representatives": representatives,
    "origin_base_indices": [0, 1, 2],
    "origin_x_polynomial": encode_poly(origin_h),
    "origin_y_polynomial": encode_poly(origin_y),
    "class_identity": "G_cut=B+T+U and F_cut=A+U+D imply [D-T]=[B-A]",
    "elapsed_ns_unisolated": time.monotonic_ns() - start,
    "source_sha256": sha(__file__),
    "input_receipt_sha256": sha(source_path),
    "plane_receipt_sha256": sha(plane_path),
    "runtime_sha256": sha(HERE / "type_i_line_fixed_origin_runtime_info.json"),
    "boundary": "Fixed-origin degree-three divisors for successful line quotient classes; independent native section replay and quotient map identity follow separately.",
}
(HERE / "type_i_line_fixed_origin_reduction.json").write_text(
    json.dumps(out, indent=2, sort_keys=True) + "\n")
print("LINE_FIXED_ORIGIN_PASS", json.dumps({
    "reduced_classes": successes, "unreduced_classes": 141 - successes}, sort_keys=True))
