"""Independently replay the three quotient pushforwards for a new-point slice."""

import argparse
import hashlib
import json
import time
from pathlib import Path

from sage.all import EllipticCurve, Integer

import type_i_divisor_transfer as dt

HERE = Path(__file__).resolve().parent
parser = argparse.ArgumentParser()
parser.add_argument("--start", type=int, required=True)
parser.add_argument("--count", type=int, required=True)
parser.add_argument("--stride", type=int, default=1)
parser.add_argument("--out", type=Path, required=True)
args = parser.parse_args()
input_path = HERE / "type_i_line_full_new_quotient_points.json"
points = json.loads(input_path.read_text())
assert 0 <= args.start < len(points)
assert args.stride > 0 and args.count > 0
assert args.start + (args.count - 1) * args.stride < len(points)
p, r, L = dt.p, dt.r, dt.L
E = EllipticCurve(L, [L(dt.gm.a), L(dt.gm.b)])
O = E(0)
inv2 = Integer(pow(2, -1, int(r)))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def images(point):
    x = point[0]
    ys = dt.check_model(point)
    out = []
    for j in range(3):
        power = p ** j
        rho = dt.inf ** power
        assert x != rho
        a = dt.c0 ** power + dt.d0 ** power / (x - rho)
        b = dt.kappa ** power * ys[j] / (x - rho) ** 2
        out.append(E(a, b))
    return out


def fiber_sum(P):
    fiber = dt.fiber(P)
    assert len(fiber) == 2
    total = [O, O, O]
    for power in range(3):
        for point in fiber:
            conjugate = tuple(value ** (p ** power) for value in point)
            image = images(conjugate)
            assert image[power] == P
            total = [total[j] + image[j] for j in range(3)]
    return total


start_ns = time.monotonic_ns()
for offset in range(args.count):
    index = args.start + offset * args.stride
    pair = points[index]
    Q = E(L(Integer(pair[0])), L(Integer(pair[1])))
    assert Q != O and r * Q == O
    H = inv2 * Q
    assert 2 * H == Q
    q_sum, h_sum = fiber_sum(Q), fiber_sum(H)
    assert all(q_sum[j] - h_sum[j] == Q for j in range(3))
    if (offset + 1) % 16 == 0:
        print("NEW_QUOTIENT_IMAGE_PROGRESS", index + 1,
              len(points), flush=True)

out = {
    "status": "verified_three_quotient_images_for_new_slice",
    "start_index": args.start,
    "stride": args.stride,
    "classes": args.count,
    "components_per_class": 3,
    "identity": "pushforward(F_Q-F_H)=(Q,Q,Q), H=2^-1 Q on E[r]",
    "elapsed_ns_unisolated": time.monotonic_ns() - start_ns,
    "source_sha256": sha(Path(__file__)),
    "input_sha256": sha(input_path),
    "runtime_info_sha256": sha(HERE / "type_i_line_new64_runtime_info.json"),
}
args.out.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n")
print("NEW_QUOTIENT_IMAGES_PASS", json.dumps({
    "start_index": args.start, "stride": args.stride,
    "classes": args.count,
    "elapsed_s": round(out["elapsed_ns_unisolated"] / 1e9, 3),
}, sort_keys=True), flush=True)
