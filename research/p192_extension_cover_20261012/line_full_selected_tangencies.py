"""Check tangent incidences between all frozen rational quartic points."""

import hashlib
import json
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sha = lambda path: hashlib.sha256(Path(path).read_bytes()).hexdigest()
inputs = ("type_i_plane_quartic.json", "type_i_line_base_points.tsv",
          "type_i_line_new_points.tsv")
plane = json.loads((HERE / inputs[0]).read_text())
p = int(plane["p"])
old = [tuple(map(int, row.split()[1:])) for row in
       (HERE / inputs[1]).read_text().splitlines()]
new = [tuple(map(int, row.split())) for row in
       (HERE / inputs[2]).read_text().splitlines()]
points = old + new
assert len(old) == 141 and len(new) == 9784
assert len(set(points)) == 9925
terms = [(tuple(map(int, key.split(","))), int(c) % p)
         for key, c in plane["quartic_coefficients"].items()]
assert all(sum(exponents) == 4 for exponents, _ in terms)

gradients = []
for x, y in points:
    assert sum(c * pow(x, a, p) * pow(y, b, p)
               for (a, b, _), c in terms) % p == 0
    fx = sum(a * c * pow(x, a - 1, p) * pow(y, b, p)
             for (a, b, _), c in terms if a) % p
    fy = sum(b * c * pow(x, a, p) * pow(y, b - 1, p)
             for (a, b, _), c in terms if b) % p
    assert fx or fy
    gradients.append((fx, fy))

start = time.monotonic_ns()
witnesses = []
n = len(points)
for i in range(n):
    x0, y0 = points[i]
    fx0, fy0 = gradients[i]
    for j in range(i + 1, n):
        x1, y1 = points[j]
        dx, dy = x1 - x0, y1 - y0
        tangent_at_i = (fx0 * dx + fy0 * dy) % p == 0
        fx1, fy1 = gradients[j]
        tangent_at_j = (fx1 * dx + fy1 * dy) % p == 0
        if tangent_at_i or tangent_at_j:
            witnesses.append({"indices_zero_based": [i, j],
                              "tangent_at_first": tangent_at_i,
                              "tangent_at_second": tangent_at_j})
    if (i + 1) % 1000 == 0:
        print("FULL_SELECTED_TANGENCY_PROGRESS", i + 1, n,
              len(witnesses), flush=True)

out = {
    "status": "FULL_SELECTED_POINT_TANGENCY_CENSUS",
    "point_count": n,
    "unordered_pairs_visited": n * (n - 1) // 2,
    "tangent_pair_count": len(witnesses),
    "witnesses": witnesses,
    "source_sha256": sha(__file__),
    "input_sha256": {name: sha(HERE / name) for name in inputs},
    "elapsed_ns_unisolated": time.monotonic_ns() - start,
    "boundary": "Checks tangency at selected endpoints only; ramification at residual points and nontransverse lines outside the selected set require separate tests.",
}
(HERE / "type_i_line_full_selected_tangencies.json").write_text(
    json.dumps(out, indent=2, sort_keys=True) + "\n")
print("FULL_SELECTED_TANGENCY_PASS", out["unordered_pairs_visited"],
      len(witnesses), flush=True)
