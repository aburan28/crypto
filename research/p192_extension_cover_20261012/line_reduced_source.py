"""Certify fixed-origin degree-three source representatives for line points.

The r-primary quotient classes in the projection receipt are h times these
source classes. This screen does not reduce that h-fold class to degree three.
"""

import hashlib
import json
from math import comb
from pathlib import Path

HERE = Path(__file__).resolve().parent
sha = lambda path: hashlib.sha256(Path(path).read_bytes()).hexdigest()
plane = json.loads((HERE / "type_i_plane_quartic.json").read_text())
anchor_path = HERE / "type_i_line_anchor_plane.json"
anchor = json.loads(anchor_path.read_text())
p = int(plane["p"])
assert plane["smooth_plane_quartic_verified"] is True
assert anchor["status"] == "verified_quotient_base_anchor_plane_coordinate"
assert int(anchor["field_prime"]) == p
base_path = HERE / "type_i_line_base_points.tsv"
new_path = HERE / "type_i_line_new_points.tsv"
pair_path = HERE / "type_i_line_new_point_pairs.tsv"
base = [(u, v) for _, u, v in
        (tuple(map(int, line.split())) for line in base_path.read_text().splitlines())]
new = [tuple(map(int, line.split())) for line in new_path.read_text().splitlines()]
split_pairs = {tuple(map(int, line.split())) for line in pair_path.read_text().splitlines()}
points = base + new
assert len(base) == 141 and len(new) == 9784 and len(set(points)) == 9925


def quartic_on_line(m, b):
    coeff = [0] * 5
    for label, value in plane["quartic_coefficients"].items():
        xpow, ypow, zpow = map(int, label.split(","))
        assert xpow + ypow + zpow == 4
        for k in range(ypow + 1):
            coeff[xpow + k] += int(value) * comb(ypow, k) * m ** k * b ** (ypow - k)
    return [c % p for c in coeff]


def divide_linear(coeff, root):
    out = [0] * (len(coeff) - 1)
    out[-1] = coeff[-1]
    for i in range(len(out) - 2, -1, -1):
        out[i] = (coeff[i + 1] + root * out[i + 1]) % p
    assert (coeff[0] + root * out[0]) % p == 0
    return out


def evaluate(coeff, x):
    result = 0
    for c in reversed(coeff):
        result = (result * x + c) % p
    return result


chosen = None
for i in range(141):
    for j in range(i + 1, 141):
        if (i + 1, j + 1) in split_pairs:
            continue
        B, C = base[i], base[j]
        m = (C[1] - B[1]) * pow(C[0] - B[0], -1, p) % p
        b = (B[1] - m * B[0]) % p
        residual = divide_linear(divide_linear(quartic_on_line(m, b), B[0]), C[0])
        assert len(residual) == 3 and residual[2] != 0
        disc = (residual[1] ** 2 - 4 * residual[2] * residual[0]) % p
        assert pow(disc, (p - 1) // 2, p) == p - 1
        assert evaluate(residual, B[0]) != 0 and evaluate(residual, C[0]) != 0
        chosen = (i, j, B, C, m, b, residual)
        break
    if chosen is not None:
        break
assert chosen is not None
i, j, B, C, m, b, residual = chosen
A = tuple(anchor["plane_affine"])
assert A not in points
assert (A[1] - m * A[0] - b) % p != 0
on_line = [k for k, P in enumerate(points)
           if (P[1] - m * P[0] - b) % p == 0]
assert on_line == [i, j]

# A smooth plane quartic has canonical line sections. Because this line
# meets B and C simply and has no further rational points, P+B+C contains
# no canonical line section for every selected P, including P=B and P=C.
out = {
    "status": "verified_fixed_origin_degree_three_source_representatives",
    "field_prime": str(p), "source_point_count": len(points),
    "origin_divisor": ["A", "B", "C"],
    "anchor_A": {"quotient_base_pilot_index": anchor["anchor_pilot_index"],
                 "cover_x": anchor["anchor_x"], "coordinates": list(A)},
    "fixed_B": {"base_index_zero": i, "coordinates": list(B)},
    "fixed_C": {"base_index_zero": j, "coordinates": list(C)},
    "line_slope": str(m), "line_intercept": str(b),
    "residual_coefficients_ascending": list(map(str, residual)),
    "residual_discriminant_character": -1,
    "line_point_indices_zero": on_line,
    "representative_rule": "D_P=P+B+C relative to D_0=A+B+C; [D_P-D_0]=[P-A]",
    "nonspecial_rule": "No selected P beyond B,C lies on BC; BC meets B,C simply and has irreducible quadratic residual. Thus no length-three D_P lies in a canonical line section.",
    "projection_boundary": "A is the exact quotient-base cover anchor. The r-torsion quotient key is the cofactor h applied to [P-A]. D_P is a reduced source representative, not a degree-three representative of h[P-A].",
    "source_sha256": sha(__file__),
    "input_sha256": {"base": sha(base_path), "new": sha(new_path),
                     "pairs": sha(pair_path),
                     "plane": sha(HERE / "type_i_plane_quartic.json"),
                     "anchor": sha(anchor_path)},
}
(HERE / "type_i_line_reduced_source.json").write_text(
    json.dumps(out, indent=2, sort_keys=True) + "\n")
print("LINE_REDUCED_SOURCE_PASS", json.dumps({
    "source_point_count": len(points), "fixed_B": i, "fixed_C": j,
    "line_point_indices_zero": on_line}, sort_keys=True))
