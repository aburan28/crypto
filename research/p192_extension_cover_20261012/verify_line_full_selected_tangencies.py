"""Independent homogeneous-gradient replay of selected-point tangencies."""

import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
sha = lambda path: hashlib.sha256(Path(path).read_bytes()).hexdigest()
receipt_path = HERE / "type_i_line_full_selected_tangencies.json"
receipt = json.loads(receipt_path.read_text())
assert receipt["source_sha256"] == sha(HERE / "line_full_selected_tangencies.py")
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

incidences = {}
for i, (x, y) in enumerate(points):
    fx = fy = fz = 0
    for (a, b, c), value in terms:
        if a:
            fx += a * value * pow(x, a - 1, p) * pow(y, b, p)
        if b:
            fy += b * value * pow(x, a, p) * pow(y, b - 1, p)
        if c:
            fz += c * value * pow(x, a, p) * pow(y, b, p)
    fx, fy, fz = fx % p, fy % p, fz % p
    assert fx or fy
    assert (fx * x + fy * y + fz) % p == 0
    if fy:
        inv = pow(fy, -1, p)
        slope = -fx * inv % p
        intercept = -fz * inv % p
        assert (y - slope * x - intercept) % p == 0
        for j, (u, v) in enumerate(points):
            if j != i and (v - slope * u - intercept) % p == 0:
                pair = (min(i, j), max(i, j))
                incidences.setdefault(pair, [False, False])[i == pair[1]] = True
    else:
        assert fx
        abscissa = -fz * pow(fx, -1, p) % p
        assert x == abscissa
        for j, (u, _) in enumerate(points):
            if j != i and u == abscissa:
                pair = (min(i, j), max(i, j))
                incidences.setdefault(pair, [False, False])[i == pair[1]] = True
    if (i + 1) % 1000 == 0:
        print("FULL_SELECTED_TANGENCY_INDEPENDENT_PROGRESS", i + 1,
              n, len(incidences), flush=True)

expected = {}
for row in receipt["witnesses"]:
    pair = tuple(row["indices_zero_based"])
    assert pair[0] < pair[1]
    expected[pair] = [row["tangent_at_first"], row["tangent_at_second"]]
assert incidences == expected
assert receipt["unordered_pairs_visited"] == n * (n - 1) // 2
assert receipt["tangent_pair_count"] == len(expected)
out = {
    "status": "FULL_SELECTED_TANGENCY_INDEPENDENT_PASS",
    "point_count": n,
    "tangent_pair_count": len(incidences),
    "source_sha256": sha(__file__),
    "receipt_sha256": sha(receipt_path),
}
(HERE / "type_i_line_full_selected_tangencies_verification.json").write_text(
    json.dumps(out, indent=2, sort_keys=True) + "\n")
print("FULL_SELECTED_TANGENCY_INDEPENDENT_PASS", n,
      len(incidences), flush=True)
