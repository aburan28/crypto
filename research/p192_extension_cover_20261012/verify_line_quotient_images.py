"""Independently replay all three quotient images of 141 fiber differences."""

import hashlib
import json
from pathlib import Path

from sage.all import EllipticCurve, Integer

import type_i_divisor_transfer as dt

HERE = Path(__file__).resolve().parent
sha = lambda path: hashlib.sha256(Path(path).read_bytes()).hexdigest()
quotient_path = HERE / "type_i_line_guided_collector.json"
cubic_path = HERE / "type_i_line_quotient_cubic_probe.json"
reduction_path = HERE / "type_i_line_fixed_origin_reduction.json"
quotient = json.loads(quotient_path.read_text())
cubic = json.loads(cubic_path.read_text())
reduction = json.loads(reduction_path.read_text())
assert len(quotient["base_points_fp"]) == len(cubic["rows"]) == len(reduction["representatives"]) == 141
p, r, L = dt.p, dt.r, dt.L
E = EllipticCurve(L, [L(dt.gm.a), L(dt.gm.b)])
O = E(0)
inv2 = Integer(pow(2, -1, int(r)))


def images(point):
    X = point[0]
    Y = dt.check_model(point)
    result = []
    for j in range(3):
        power = p ** j
        rho = dt.inf ** power
        assert X != rho
        a = dt.c0 ** power + dt.d0 ** power / (X - rho)
        b = dt.kappa ** power * Y[j] / (X - rho) ** 2
        result.append(E(a, b))
    return result


def fiber_sum(P):
    fiber = dt.fiber(P)
    assert len(fiber) == 2
    total = [O, O, O]
    for i in range(3):
        for point in fiber:
            conjugate = tuple(c ** (p ** i) for c in point)
            value = images(conjugate)
            assert value[i] == P
            total = [total[j] + value[j] for j in range(3)]
    return total


checked = 0
for index, pair in enumerate(quotient["base_points_fp"]):
    Q = E(L(Integer(pair[0])), L(Integer(pair[1])))
    assert Q != O and r * Q == O
    H = inv2 * Q
    assert 2 * H == Q
    assert cubic["rows"][index]["quotient"] == list(pair)
    assert reduction["representatives"][index]["quotient"] == list(pair)
    q_sum, h_sum = fiber_sum(Q), fiber_sum(H)
    assert all(q_sum[j] - h_sum[j] == Q for j in range(3))
    checked += 1
    if checked % 16 == 0:
        print(f"LINE_QUOTIENT_IMAGE_PROGRESS {checked}/141", flush=True)

out = {
    "status": "LINE_QUOTIENT_IMAGES_INDEPENDENT_PASS",
    "classes": checked, "components_per_class": 3,
    "identity": "Each conjugate fiber pair maps to Q or H in its matching quotient; pushforward(F_Q-F_H)=(Q,Q,Q), H=2^-1 Q on E[r]",
    "source_sha256": sha(__file__),
    "quotient_receipt_sha256": sha(quotient_path),
    "cubic_receipt_sha256": sha(cubic_path),
    "reduction_receipt_sha256": sha(reduction_path),
    "runtime_sha256": sha(HERE / "type_i_line_quotient_image_runtime_info.json"),
}
(HERE / "type_i_line_quotient_image_verification.json").write_text(
    json.dumps(out, indent=2, sort_keys=True) + "\n")
print("LINE_QUOTIENT_IMAGES_INDEPENDENT_PASS", checked)
