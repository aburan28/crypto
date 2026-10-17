"""Freeze the exact plane coordinate of the quotient-base cover anchor."""

import hashlib
import json
import random
from pathlib import Path

from sage.all import Integer, prod, vector

import type_i_divisor_transfer as dt

HERE = Path(__file__).resolve().parent
sha = lambda path: hashlib.sha256(Path(path).read_bytes()).hexdigest()
base_path = HERE / "type_i_quotient_base_probe.json"
plane_path = HERE / "type_i_plane_quartic.json"
base = json.loads(base_path.read_text())
plane = json.loads(plane_path.read_text())
p, L = dt.p, dt.L
rng = random.Random(base["pilot_seed"])
x = [Integer(rng.randrange(int(p)))
     for _ in range(base["anchor_pilot_index"] + 1)][-1]
X = L(x)
U = (X - dt.zeta) * (X - dt.beta)
y = U ** ((dt.q + 1) // 2)
N = prod(X - root ** (p ** j)
         for root in (dt.zeta, dt.beta) for j in range(3))
if y * y ** p * y ** (p * p) == -N:
    y = -y
assert y * y ** p * y ** (p * p) == N
Y = [y, y ** p, y ** (p * p)]
cover = (X, *dt.Vandermonde.solve_right(vector(L, Y)))
dt.check_model(cover)
coords = dt.Vandermonde.solve_right(vector(L, [1 / value for value in Y]))
assert coords[2] != 0
u, v = coords[0] / coords[2], coords[1] / coords[2]
assert u ** p == u and v ** p == v
up = int(u.polynomial()[0])
vp = int(v.polynomial()[0])
value = 0
for label, coefficient in plane["quartic_coefficients"].items():
    a, b, c = map(int, label.split(","))
    value += int(coefficient) * pow(up, a, int(p)) * pow(vp, b, int(p))
assert value % int(p) == 0
out = {
    "status": "verified_quotient_base_anchor_plane_coordinate",
    "field_prime": str(p), "anchor_pilot_index": base["anchor_pilot_index"],
    "anchor_x": str(x), "plane_affine": [up, vp],
    "source_sha256": sha(__file__), "base_receipt_sha256": sha(base_path),
    "plane_receipt_sha256": sha(plane_path),
    "runtime_sha256": sha(HERE / "type_i_line_anchor_runtime_info.json"),
}
(HERE / "type_i_line_anchor_plane.json").write_text(
    json.dumps(out, indent=2, sort_keys=True) + "\n")
print("LINE_ANCHOR_PLANE_PASS", json.dumps(out, sort_keys=True))
