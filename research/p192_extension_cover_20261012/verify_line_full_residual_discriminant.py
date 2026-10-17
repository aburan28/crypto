"""Reverse-order, block-inversion replay of the full quartic residual screen."""

import hashlib
import json
from math import comb
from pathlib import Path

HERE = Path(__file__).resolve().parent
sha = lambda path: hashlib.sha256(Path(path).read_bytes()).hexdigest()
receipt_path = HERE / "type_i_line_full_residual_discriminant.json"
receipt = json.loads(receipt_path.read_text())
assert receipt["source_sha256"] == sha(HERE / "line_full_residual_discriminant.py")
assert all(digest == sha(HERE / name)
           for name, digest in receipt["input_sha256"].items())
plane = json.loads((HERE / "type_i_plane_quartic.json").read_text())
p = int(plane["p"])
old = [tuple(map(int, row.split()[1:])) for row in
       (HERE / "type_i_line_base_points.tsv").read_text().splitlines()]
new = [tuple(map(int, row.split())) for row in
       (HERE / "type_i_line_new_points.tsv").read_text().splitlines()]
points = old + new
n = len(points)
assert n == 9925 and len(set(points)) == n
terms = [(tuple(map(int, key.split(","))), int(c) % p)
         for key, c in plane["quartic_coefficients"].items()]
assert all(sum(exponents) == 4 for exponents, _ in terms)


def direction_tensor(x, y):
    xp = [pow(x, k, p) for k in range(5)]
    yp = [pow(y, k, p) for k in range(5)]
    tensor = [[0] * (k + 1) for k in range(5)]
    for k in range(5):
        for j in range(k + 1):
            ix = k - j
            tensor[k][j] = sum(
                c * comb(a, ix) * comb(b, j) * xp[a - ix] * yp[b - j]
                for (a, b, _), c in terms if ix <= a and j <= b
            ) % p
    assert tensor[0][0] == 0
    return tensor


def value(coeffs, m):
    total = 0
    for c in reversed(coeffs):
        total = (total * m + c) % p
    return total


def direct_line_polynomial(m, b):
    coeffs = [0] * 5
    for (a, e, _), c in terms:
        for j in range(e + 1):
            coeffs[a + j] = (coeffs[a + j] + c * comb(e, j)
                             * pow(m, j, p) * pow(b, e - j, p)) % p
    return coeffs


# An independent quartic substitution validates the Taylor tensor on every
# pair in the preceding 11,040-pair second-generation control selection.
control_points = old + new[:64]
controls = 0
for i in range(141, 205):
    x, y = control_points[i]
    tensor = direction_tensor(x, y)
    for j in range(i):
        u, v = control_points[j]
        d = (u - x) % p
        assert d
        m = (v - y) * pow(d, -1, p) % p
        affine = direct_line_polynomial(m, (y - m * x) % p)
        shifted = [sum(affine[k] * comb(k, h) * pow(x, k - h, p)
                       for k in range(h, 5)) % p for h in range(5)]
        assert shifted == [value(tensor[k], m) for k in range(5)]
        controls += 1
assert controls == 11040

counts = {"pairs": 0, "vertical_pairs": 0, "degree_below_four": 0,
          "residual_double_root": 0, "tangent_at_first": 0,
          "tangent_at_second": 0}
reversed_points = list(reversed(points))
for i, (x, y) in enumerate(reversed_points):
    tensor = direction_tensor(x, y)
    for start in range(i + 1, n, 1024):
        rows, prefixes, vertical = [], [], []
        product = 1
        for j in range(start, min(n, start + 1024)):
            u, v = reversed_points[j]
            d = (u - x) % p
            if d:
                rows.append((j, d, (v - y) % p))
                product = product * d % p
                prefixes.append(product)
            else:
                vertical.append((j, (v - y) % p))
        if rows:
            inverse_product = pow(product, -1, p)
            inverses = [0] * len(rows)
            for k in range(len(rows) - 1, -1, -1):
                inverses[k] = inverse_product * (prefixes[k - 1] if k else 1) % p
                inverse_product = inverse_product * rows[k][1] % p
            assert inverse_product == 1
            for (j, d, dy), inverse_d in zip(rows, inverses):
                m = dy * inverse_d % p
                c2 = value(tensor[2], m)
                c3 = value(tensor[3], m)
                A = value(tensor[4], m)
                B = (c3 + d * A) % p
                C = (c2 + d * B) % p
                counts["pairs"] += 1
                if A == 0:
                    counts["degree_below_four"] += 1
                elif (B * B - 4 * A * C) % p == 0:
                    counts["residual_double_root"] += 1
                if C == 0:
                    counts["tangent_at_second"] += 1
                if (A * d * d + B * d + C) % p == 0:
                    counts["tangent_at_first"] += 1
                if (i + j) % 1000003 == 0:
                    assert (value(tensor[1], m) + d * C) % p == 0
        for _, d in vertical:
            assert d
            A = tensor[4][4]
            B = (tensor[3][3] + d * A) % p
            C = (tensor[2][2] + d * B) % p
            counts["pairs"] += 1
            counts["vertical_pairs"] += 1
            if A == 0:
                counts["degree_below_four"] += 1
            elif (B * B - 4 * A * C) % p == 0:
                counts["residual_double_root"] += 1
            if C == 0:
                counts["tangent_at_second"] += 1
            if (A * d * d + B * d + C) % p == 0:
                counts["tangent_at_first"] += 1
    if (i + 1) % 1000 == 0:
        print("FULL_RESIDUAL_INDEPENDENT_PROGRESS", i + 1, n,
              counts["pairs"], flush=True)

assert counts == receipt["counts"]
assert all(counts[k] == 0 for k in counts if k != "pairs")
assert counts["pairs"] == n * (n - 1) // 2
assert receipt["exception_stream_sha256"] == hashlib.sha256(b"").hexdigest()
assert receipt["exception_witnesses_first_1000"] == []
out = {
    "status": "FULL_RESIDUAL_DISCRIMINANT_INDEPENDENT_PASS",
    "point_count": n,
    "unordered_pairs_visited": counts["pairs"],
    "direct_substitution_controls": controls,
    "counts": counts,
    "source_sha256": sha(__file__),
    "receipt_sha256": sha(receipt_path),
    "input_sha256": {name: sha(HERE / name) for name in
                     ("type_i_plane_quartic.json",
                      "type_i_line_base_points.tsv",
                      "type_i_line_new_points.tsv")},
}
(HERE / "type_i_line_full_residual_discriminant_verification.json").write_text(
    json.dumps(out, indent=2, sort_keys=True) + "\n")
print("FULL_RESIDUAL_DISCRIMINANT_INDEPENDENT_PASS", counts,
      flush=True)
