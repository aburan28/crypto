#!/usr/bin/env python3
"""Independent prime-field replay of the P-192 line-section receipt."""

import collections
import hashlib
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
sha = lambda path: hashlib.sha256(Path(path).read_bytes()).hexdigest()
support = json.loads((HERE / "type_i_divisor_support.json").read_bytes())
plane = json.loads((HERE / "type_i_plane_quartic.json").read_bytes())
p = int(support["field_prime"])
assert p == int(plane["p"])
assert support["status"] == "exact_fixed_origin_divisor_support_audit"
assert support["squarefree_cubics"] == 160
assert plane["smooth_plane_quartic_verified"] is True

base_path = HERE / "type_i_line_base_points.tsv"
new_path = HERE / "type_i_line_new_points.tsv"
pair_path = HERE / "type_i_line_new_point_pairs.tsv"
base = [tuple(map(int, line.split())) for line in base_path.read_text().splitlines()]
new = [tuple(map(int, line.split())) for line in new_path.read_text().splitlines()]
pairids = [tuple(map(int, line.split())) for line in pair_path.read_text().splitlines()]
assert len(base) == 141 and len(new) == len(pairids) == 9784
assert all(len(row) == 3 for row in base)
assert all(len(row) == 2 for row in new + pairids)
assert all(0 <= x < p and 0 <= y < p for _, x, y in base)
assert all(0 <= x < p and 0 <= y < p for x, y in new)

terms = []
for key, value in plane["quartic_coefficients"].items():
    a, b, c = map(int, key.split(","))
    assert a + b + c == 4
    terms.append((a, b, int(value) % p))


def curve_at(x, y):
    return sum(c * pow(x, a, p) * pow(y, b, p) for a, b, c in terms) % p


def poly_at(coefficients, x):
    return sum(int(c) * pow(x, i, p) for i, c in enumerate(coefficients)) % p


def trim(poly):
    while poly and poly[-1] % p == 0:
        poly.pop()
    return poly


def line_polynomial(m, b):
    result = [0] * 5
    for xa, ya, coeff in terms:
        for j in range(ya + 1):
            result[xa + j] = (result[xa + j] + coeff * math.comb(ya, j)
                              * pow(m, j, p) * pow(b, ya - j, p)) % p
    return trim(result)


def divide_by_pair(poly, x1, x2):
    divisor = [(x1 * x2) % p, (-x1 - x2) % p, 1]
    work = poly[:]
    quotient = [0] * max(0, len(work) - 2)
    while len(work) >= 3:
        scale = work[-1] % p
        offset = len(work) - 3
        quotient[offset] = scale
        for j in range(3):
            work[offset + j] = (work[offset + j] - scale * divisor[j]) % p
        trim(work)
    assert not work, "line does not contain both selected base points"
    return trim(quotient)


for class_index, x, y in base:
    assert 1 <= class_index <= 160
    record = support["records"][class_index - 1]
    assert poly_at(record["x_coefficients"], x) == 0
    assert poly_at(record["y_coefficients"], x) == y
    assert curve_at(x, y) == 0

base_xy = [(x, y) for _, x, y in base]
assert len(set(base_xy)) == len({x for x, _ in base_xy}) == 141
assert all(curve_at(x, y) == 0 for x, y in new)
assert len(set(new)) == 9784 and not set(base_xy).intersection(new)
new_by_pair = collections.defaultdict(list)
for pair, point in zip(pairids, new):
    assert len(pair) == 2 and 1 <= pair[0] < pair[1] <= 141
    new_by_pair[pair].append(point)
assert len(new_by_pair) == 4892
assert all(len(values) == 2 and values[0] != values[1]
           for values in new_by_pair.values())

line_keys = set()
split = nonsplit = 0
for i in range(141):
    x1, y1 = base_xy[i]
    for j in range(i + 1, 141):
        x2, y2 = base_xy[j]
        assert x1 != x2
        m = ((y2 - y1) * pow((x2 - x1) % p, -1, p)) % p
        b = (y1 - m * x1) % p
        key = (m, b)
        assert key not in line_keys, "three selected points are collinear"
        line_keys.add(key)
        residual = divide_by_pair(line_polynomial(m, b), x1, x2)
        assert len(residual) == 3 and residual[2] != 0
        disc = (residual[1] ** 2 - 4 * residual[2] * residual[0]) % p
        assert disc != 0
        points = new_by_pair.get((i + 1, j + 1), [])
        if pow(disc, (p - 1) // 2, p) == 1:
            split += 1
            assert len(points) == 2
            for x, y in points:
                assert poly_at(residual, x) == 0
                assert y == (m * x + b) % p
        else:
            nonsplit += 1
            assert pow(disc, (p - 1) // 2, p) == p - 1
            assert not points

assert len(line_keys) == 9870 and (nonsplit, split) == (4978, 4892)
log = (HERE / "type_i_line_support_screen.log").read_text().splitlines()
assert len(log) == 2 and "***" not in "\n".join(log)
assert log[0] == ("LINE_RESIDUAL_GP_PASS points=141 pairs=9870 "
                  "degree_counts=[0, 0, 9870, 0] repeated_quadratic_roots=0 "
                  "nonsplit_split=[4978, 4892] new_rational_points=9784 "
                  "all_points_distinct=9925")
assert log[1] == f"FIRST_NEW_RATIONAL_POINT=[{new[0][0]}, {new[0][1]}]"

out = {
    "status": "LINE_SUPPORT_INDEPENDENT_PASS",
    "field_prime": str(p), "base_points": 141, "pair_lines": 9870,
    "nonsplit_lines": nonsplit, "split_lines": split,
    "new_rational_points": len(new), "all_affine_points_distinct": 9925,
    "source_sha256": sha(__file__),
    "input_sha256": {name: sha(path) for name, path in
                     (("support", HERE / "type_i_divisor_support.json"),
                      ("plane", HERE / "type_i_plane_quartic.json"),
                      ("base_points", base_path), ("new_points", new_path),
                      ("new_point_pairs", pair_path),
                      ("gp_source", HERE / "line_support_screen.gp"),
                      ("gp_log", HERE / "type_i_line_support_screen.log"))},
}
(HERE / "type_i_line_support_verification.json").write_text(
    json.dumps(out, indent=2, sort_keys=True) + "\n")
print(json.dumps(out, sort_keys=True))
