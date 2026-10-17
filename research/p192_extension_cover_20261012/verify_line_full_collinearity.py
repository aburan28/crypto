"""Independent reverse-order, block-inversion replay of all point lines."""

import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
sha = lambda path: hashlib.sha256(Path(path).read_bytes()).hexdigest()
source = HERE / "line_full_collinearity.py"
receipt_path = HERE / "type_i_line_full_collinearity.json"
receipt = json.loads(receipt_path.read_text())
plane = json.loads((HERE / "type_i_plane_quartic.json").read_text())
p = int(plane["p"])
assert receipt["source_sha256"] == sha(source)
assert all(digest == sha(HERE / name)
           for name, digest in receipt["input_sha256"].items())
assert receipt["field_prime"] == str(p)

old = [tuple(map(int, row.split()[1:])) for row in
       (HERE / "type_i_line_base_points.tsv").read_text().splitlines()]
new = [tuple(map(int, row.split())) for row in
       (HERE / "type_i_line_new_points.tsv").read_text().splitlines()]
points = old + new
n = len(points)
assert len(old) == 141 and len(new) == 9784 and n == 9925
assert len(set(points)) == n
terms = [(tuple(map(int, exponents.split(","))), int(c) % p)
         for exponents, c in plane["quartic_coefficients"].items()]
assert all(sum(exponents) == 4 for exponents, _ in terms)
for x, y in points:
    assert 0 <= x < p and 0 <= y < p
    assert sum(c * pow(x, a, p) * pow(y, b, p)
               for (a, b, _), c in terms) % p == 0

expected = {}
for row in receipt["line_records"]:
    slope, intercept = map(int, row["line"])
    members = set(row["point_indices_zero_based"])
    assert 3 <= len(members) <= 4
    assert (slope, intercept) not in expected
    expected[slope, intercept] = members

reversed_points = list(reversed(points))
found = {}
pair_count = 0
for i, (x0, y0) in enumerate(reversed_points):
    orig_i = n - 1 - i
    first = {}
    vertical_first = None
    for offset in range(i + 1, n, 1024):
        block = []
        prefixes = []
        product = 1
        for j in range(offset, min(offset + 1024, n)):
            xj, yj = reversed_points[j]
            orig_j = n - 1 - j
            pair_count += 1
            dx = (xj - x0) % p
            if not dx:
                if vertical_first is None:
                    vertical_first = orig_j
                else:
                    found.setdefault((-1, x0), set()).update(
                        (orig_i, vertical_first, orig_j))
                continue
            block.append((orig_j, dx, (yj - y0) % p))
            product = product * dx % p
            prefixes.append(product)
        if not block:
            continue
        inverses = [0] * len(block)
        inverse_product = pow(product, -1, p)
        for k in range(len(block) - 1, -1, -1):
            inverses[k] = inverse_product * (prefixes[k - 1] if k else 1) % p
            inverse_product = inverse_product * block[k][1] % p
        assert inverse_product == 1
        for (orig_j, _, dy), inverse_dx in zip(block, inverses):
            slope = dy * inverse_dx % p
            if slope in first:
                intercept = (y0 - slope * x0) % p
                found.setdefault((slope, intercept), set()).update(
                    (orig_i, first[slope], orig_j))
            else:
                first[slope] = orig_j
    if (i + 1) % 500 == 0:
        print("FULL_COLLINEARITY_INDEPENDENT_PROGRESS", i + 1, n,
              pair_count, len(found), flush=True)

assert pair_count == n * (n - 1) // 2 == receipt["unordered_pairs_visited"]
assert found == expected
pairs = [tuple(map(int, row.split())) for row in
         (HERE / "type_i_line_new_point_pairs.tsv").read_text().splitlines()]
assert len(pairs) == 9784
source_quartets = set()
for k in range(4892):
    a, b = pairs[2 * k]
    assert pairs[2 * k + 1] == (a, b)
    source_quartets.add(frozenset((a - 1, b - 1,
                                   len(old) + 2 * k,
                                   len(old) + 2 * k + 1)))
assert len(source_quartets) == 4892
assert {frozenset(members) for members in found.values()} == source_quartets
out = {
    "status": "FULL_COLLINEARITY_INDEPENDENT_PASS",
    "point_count": n,
    "unordered_pairs_visited": pair_count,
    "lines_with_at_least_three_selected_points": len(found),
    "source_sha256": sha(__file__),
    "receipt_sha256": sha(receipt_path),
    "input_sha256": {name: sha(HERE / name) for name in
                     ("type_i_plane_quartic.json",
                      "type_i_line_base_points.tsv",
                      "type_i_line_new_points.tsv",
                      "type_i_line_new_point_pairs.tsv")},
}
(HERE / "type_i_line_full_collinearity_verification.json").write_text(
    json.dumps(out, indent=2, sort_keys=True) + "\n")
print("FULL_COLLINEARITY_INDEPENDENT_PASS", pair_count, len(found),
      flush=True)
