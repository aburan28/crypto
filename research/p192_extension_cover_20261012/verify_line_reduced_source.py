"""Independent polynomial replay of the fixed-origin source divisor rule."""

import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
sha = lambda path: hashlib.sha256(Path(path).read_bytes()).hexdigest()
proof = json.loads((HERE / "type_i_line_reduced_source.json").read_text())
plane = json.loads((HERE / "type_i_plane_quartic.json").read_text())
p = int(plane["p"])
assert proof["status"] == "verified_fixed_origin_degree_three_source_representatives"
assert proof["source_sha256"] == sha(HERE / "line_reduced_source.py")
paths = {"base": HERE / "type_i_line_base_points.tsv",
         "new": HERE / "type_i_line_new_points.tsv",
         "pairs": HERE / "type_i_line_new_point_pairs.tsv",
         "plane": HERE / "type_i_plane_quartic.json",
         "anchor": HERE / "type_i_line_anchor_plane.json"}
assert proof["input_sha256"] == {name: sha(path) for name, path in paths.items()}
anchor = json.loads(paths["anchor"].read_text())
assert anchor["status"] == "verified_quotient_base_anchor_plane_coordinate"
assert anchor["source_sha256"] == sha(HERE / "line_anchor_plane.py")
assert anchor["runtime_sha256"] == sha(HERE / "type_i_line_anchor_runtime_info.json")
assert anchor["base_receipt_sha256"] == sha(HERE / "type_i_quotient_base_probe.json")
assert anchor["plane_receipt_sha256"] == sha(HERE / "type_i_plane_quartic.json")
base = [tuple(map(int, line.split()))[1:] for line in paths["base"].read_text().splitlines()]
new = [tuple(map(int, line.split())) for line in paths["new"].read_text().splitlines()]
points = base + new
assert len(points) == proof["source_point_count"] == 9925
assert len(set(points)) == len(points)
A = tuple(proof["anchor_A"]["coordinates"])
B = tuple(proof["fixed_B"]["coordinates"])
C = tuple(proof["fixed_C"]["coordinates"])
assert A == tuple(anchor["plane_affine"])
assert proof["anchor_A"]["quotient_base_pilot_index"] == anchor["anchor_pilot_index"]
assert proof["anchor_A"]["cover_x"] == anchor["anchor_x"]
assert A not in points
assert base[proof["fixed_B"]["base_index_zero"]] == B
assert base[proof["fixed_C"]["base_index_zero"]] == C


def multiply(u, v):
    result = [0] * (len(u) + len(v) - 1)
    for i, x in enumerate(u):
        for j, y in enumerate(v):
            result[i + j] = (result[i + j] + x * y) % p
    return result


def power(u, exponent):
    value = [1]
    for _ in range(exponent):
        value = multiply(value, u)
    return value


m = int(proof["line_slope"])
b = int(proof["line_intercept"])
coeff = [0] * 5
for label, raw in plane["quartic_coefficients"].items():
    xpow, ypow, zpow = map(int, label.split(","))
    term = multiply(power([0, 1], xpow), power([b, m], ypow))
    for k, value in enumerate(term):
        coeff[k] = (coeff[k] + int(raw) * value) % p


def divide(u, root):
    quotient = [0] * (len(u) - 1)
    quotient[-1] = u[-1]
    for k in range(len(quotient) - 2, -1, -1):
        quotient[k] = (u[k + 1] + root * quotient[k + 1]) % p
    assert (u[0] + root * quotient[0]) % p == 0
    return quotient


assert B[0] != C[0]
assert B[1] == (m * B[0] + b) % p
assert C[1] == (m * C[0] + b) % p
assert A[1] != (m * A[0] + b) % p
assert sum(int(value) * pow(A[0], int(label.split(",")[0]), p) *
           pow(A[1], int(label.split(",")[1]), p)
           for label, value in plane["quartic_coefficients"].items()) % p == 0
residual = divide(divide(coeff, B[0]), C[0])
assert list(map(int, proof["residual_coefficients_ascending"])) == residual
disc = (residual[1] ** 2 - 4 * residual[0] * residual[2]) % p
assert pow(disc, (p - 1) // 2, p) == p - 1
for x, _ in (B, C):
    assert (residual[0] + residual[1] * x + residual[2] * x * x) % p != 0
on_line = [i for i, (x, y) in enumerate(points) if (y - m * x - b) % p == 0]
assert on_line == proof["line_point_indices_zero"] == [
    proof["fixed_B"]["base_index_zero"], proof["fixed_C"]["base_index_zero"]]
print("LINE_REDUCED_SOURCE_INDEPENDENT_PASS", len(points))
