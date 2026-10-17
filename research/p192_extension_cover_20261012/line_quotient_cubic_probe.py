"""Construct prime-field cubic-section input divisors for 141 quotient points.

Every point Q is an exact r-torsion quotient class. For each Q and its unique
r-torsion half H, take both inverse-cover fiber points and three p-Frobenius
images. Record the residual cubic cut over Fp when all checks pass; retain
per-class failures rather than classifying them as absent mathematics.
"""

import hashlib
import json
import time
from pathlib import Path

from sage.all import EllipticCurve, Integer, PolynomialRing, matrix, prod, vector

import type_i_divisor_transfer as dt

HERE = Path(__file__).resolve().parent
sha = lambda path: hashlib.sha256(Path(path).read_bytes()).hexdigest()
quotient_path = HERE / "type_i_line_guided_collector.json"
eligibility_path = HERE / "type_i_line_fiber_eligibility.json"
plane_path = HERE / "type_i_plane_quartic.json"
quotient = json.loads(quotient_path.read_text())
eligibility = json.loads(eligibility_path.read_text())
plane = json.loads(plane_path.read_text())
assert quotient["base_classes"] == eligibility["input_classes"] == 141

p, r, L = dt.p, dt.r, dt.L
E = EllipticCurve(L, [L(dt.gm.a), L(dt.gm.b)])
O = E(0)
inv2 = Integer(pow(2, -1, int(r)))
anchors = [tuple(L(int(v)) for v in item["coordinates_Fp_decimal"])
           for item in plane["known_points_verified"] if item["source"] == "base"][:3]
assert len(anchors) == 3
# Match the frozen cubic-section coordinate order in the earlier proof.
monomials = [(3, 0, 0), (0, 3, 0), (0, 0, 3), (2, 1, 0), (2, 0, 1),
             (1, 2, 0), (0, 2, 1), (1, 0, 2), (0, 1, 2), (1, 1, 1)]
R = PolynomialRing(L, "x")
x = R.gen()
S = PolynomialRing(R, "y")
y = S.gen()
quartic = sum(S(L(int(c)) * x ** a) * y ** b
              for label, c in plane["quartic_coefficients"].items()
              for a, b, _ in [tuple(map(int, label.split(",")))])


def encode_poly(poly):
    return [[str(int(c.polynomial()[0]))] for c in poly.list()]


def residual(P):
    fiber = dt.fiber(P)
    points = [tuple(value ** (p ** power) for value in point)
              for power in range(3) for point in fiber]
    coordinates = []
    for point in points:
        ys = dt.check_model(point)
        assert all(value != 0 for value in ys)
        coordinates.append(dt.Vandermonde.solve_right(
            vector(L, [1 / value for value in ys])))
    coordinates.extend(anchors)
    assert len(coordinates) == 9 and all(c[2] != 0 for c in coordinates)
    rows = [[prod(c[k] ** e[k] for k in range(3)) for e in monomials]
            for c in coordinates]
    evaluation = matrix(L, rows)
    assert evaluation.rank() == 9
    coefficients = list(evaluation.right_kernel().basis()[0])
    first = next(c for c in coefficients if c != 0)
    coefficients = [c / first for c in coefficients]
    cubic = sum(S(coefficients[i] * x ** e[0]) * y ** e[1]
                for i, e in enumerate(monomials))
    known = prod((x - c[0] / c[2] for c in coordinates), R.one())
    resultant = quartic.resultant(cubic)
    quotient_poly, remainder = resultant.quo_rem(known)
    assert resultant.degree() == 12 and remainder == 0
    assert known.degree() == 9 and known.gcd(known.derivative()).degree() == 0
    h = quotient_poly.monic()
    assert h.degree() == 3 and h.gcd(h.derivative()).degree() == 0
    assert h.gcd(known).degree() == 0
    linear = next(s for s in quartic.subresultants(cubic) if s.degree() == 1)
    assert linear[1].gcd(h).degree() == 0
    y_value = (-linear[0] * linear[1].inverse_mod(h)) % h
    assert y_value.degree() < 3
    assert quartic(y_value) % h == 0 and cubic(y_value) % h == 0
    assert all(c ** p == c for c in coefficients + h.list() + y_value.list())
    return {
        "cubic_coefficients": [[str(int(c.polynomial()[0]))] for c in coefficients],
        "residual_x_polynomial": encode_poly(h),
        "residual_y_polynomial": encode_poly(y_value),
    }


start = time.monotonic_ns()
rows = []
for index, ((qx, qy), cell) in enumerate(zip(
        quotient["base_points_fp"], eligibility["rows"])):
    assert cell["index"] == index
    Q = E(L(Integer(qx)), L(Integer(qy)))
    assert Q != O and r * Q == O
    H = inv2 * Q
    assert 2 * H == Q
    row = {"index": index, "cell": cell["cell"], "quotient": [qx, qy]}
    try:
        row["positive_input_residual"] = residual(Q)
        row["negative_input_residual"] = residual(H)
        row["status"] = "verified_prime_field_cubic_pair"
    except Exception as error:
        row["status"] = "bounded_constructor_failure"
        row["error_type"] = type(error).__name__
        row["error"] = str(error)[:300]
    rows.append(row)
    if (index + 1) % 16 == 0:
        print(f"LINE_QUOTIENT_CUBIC_PROGRESS {index + 1}/141", flush=True)

successes = {cell: sum(row["status"] == "verified_prime_field_cubic_pair" and
                        row["cell"] == cell for row in rows)
             for cell in eligibility["cells"]}
out = {
    "status": "bounded_line_quotient_cubic_constructor_probe",
    "field_prime": str(p), "subgroup_order": str(r),
    "classes": 141, "successes_by_cell": successes,
    "successes": sum(successes.values()), "failures": 141 - sum(successes.values()),
    "rows": rows, "elapsed_ns_unisolated": time.monotonic_ns() - start,
    "source_sha256": sha(__file__),
    "quotient_receipt_sha256": sha(quotient_path),
    "eligibility_receipt_sha256": sha(eligibility_path),
    "plane_receipt_sha256": sha(plane_path),
    "runtime_sha256": sha(HERE / "type_i_line_cubic_runtime_info.json"),
    "boundary": "Exact prime-field cubic input pairs where verified; subsequent fixed-origin degree-three reduction and independent quotient-class replay remain separate.",
}
(HERE / "type_i_line_quotient_cubic_probe.json").write_text(
    json.dumps(out, indent=2, sort_keys=True) + "\n")
print("LINE_QUOTIENT_CUBIC_PASS", json.dumps({
    "successes": out["successes"], "failures": out["failures"],
    "successes_by_cell": successes}, sort_keys=True))
