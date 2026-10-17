"""Independent full replay of the P-192 line-point quotient projection.

This verifier does not import the producer's plane-point inverse or screen.
It reconstructs the cover point from two quartic ratios and the norm, then
recomputes the frozen first-quotient keys, old-base controls, and pair sums.
"""

import hashlib
import json
import random
from pathlib import Path

from sage.all import EllipticCurve, Integer, PolynomialRing, prod, vector

import type_i_divisor_transfer as dt

HERE = Path(__file__).resolve().parent
sha = lambda path: hashlib.sha256(Path(path).read_bytes()).hexdigest()
receipt_path = HERE / "type_i_line_projection.json"
receipt = json.loads(receipt_path.read_bytes())
base_receipt = json.loads((HERE / "type_i_quotient_base_probe.json").read_bytes())
plane_receipt = json.loads((HERE / "type_i_plane_quartic.json").read_bytes())
fixture = json.loads((HERE / "type_i_base_points_legacy_fixture.json").read_bytes())
assert receipt["status"] == "verified_line_support_first_quotient_projection_screen"
assert receipt["source_sha256"] == sha(HERE / "line_projection_screen.py")
assert receipt["inverse_source_sha256"] == sha(HERE / "plane_point_inverse.py")
assert receipt["gp_source_sha256"] == sha(HERE / "line_support_screen.gp")
assert receipt["gp_log_sha256"] == sha(HERE / "type_i_line_support_screen.log")
assert receipt["runtime_sha256"] == sha(HERE / "type_i_line_sage_runtime_info.json")
for name, path in (("base_points", HERE / "type_i_line_base_points.tsv"),
                   ("new_points", HERE / "type_i_line_new_points.tsv"),
                   ("new_point_pairs", HERE / "type_i_line_new_point_pairs.tsv")):
    assert receipt["input_sha256"][name] == sha(path)

p, L, F, r = dt.p, dt.L, dt.gm.F, dt.r
E = EllipticCurve(L, [L(dt.gm.a), L(dt.gm.b)])
O = E(0)
h = Integer(base_receipt["candidate_cofactor"])
R = PolynomialRing(L, "t")
t = R.gen()
f = [prod(t - root ** (p ** i) for root in dt.quartic_roots)
     for i in range(3)]
conjugates = [dt.zeta ** (p ** i) for i in range(3)]


def inverse(u, v):
    A = [L(u) + conjugates[i] * L(v) + conjugates[i] ** 2
         for i in range(3)]
    assert all(a != 0 for a in A)
    g1 = A[0] ** 2 * f[0] - A[1] ** 2 * f[1]
    g2 = A[0] ** 2 * f[0] - A[2] ** 2 * f[2]
    gcd = g1.gcd(g2).monic()
    assert gcd.degree() == 1
    X = -gcd[0]
    N = prod(X - root ** (p ** i)
             for root in (dt.zeta, dt.beta) for i in range(3))
    assert N != 0
    scale = A[0] * f[0](X) / (A[1] * A[2] * N)
    Y = [1 / (scale * a) for a in A]
    assert all(Y[i] ** 2 == f[i](X) for i in range(3))
    assert prod(Y) == N
    coordinates = dt.Vandermonde.solve_right(vector(L, Y))
    point = (X, *coordinates)
    dt.check_model(point)
    assert all(c ** p == c for c in point)
    back = dt.Vandermonde.solve_right(vector(L, [1 / y for y in Y]))
    assert back[0] / back[2] == L(u) and back[1] / back[2] == L(v)
    return point


def quotient(point, index):
    X = point[0]
    y = dt.check_model(point)[index]
    power = p ** index
    rho = dt.inf ** power
    c = dt.c0 ** power
    d = dt.d0 ** power
    kappa = dt.kappa ** power
    assert X != rho
    return E(c + d / (X - rho), kappa * y / (X - rho) ** 2)


def key(Q):
    if Q == O:
        return "identity"
    assert Q[0] ** p == Q[0] and Q[1] ** p == Q[1]
    return f"{int(Q[0].polynomial()[0])}:{int(Q[1].polynomial()[0])}"


seed = random.Random(base_receipt["pilot_seed"])
anchor_x = [Integer(seed.randrange(int(p)))
            for _ in range(base_receipt["anchor_pilot_index"] + 1)][-1]
X = L(anchor_x)
U = (X - dt.zeta) * (X - dt.beta)
y = U ** ((dt.q + 1) // 2)
N = prod(X - root ** (p ** i)
         for root in (dt.zeta, dt.beta) for i in range(3))
if y * y ** p * y ** (p * p) == -N:
    y = -y
assert y * y ** p * y ** (p * p) == N
coords = dt.Vandermonde.solve_right(vector(L, [y, y ** p, y ** (p * p)]))
anchor = (X, *coords)
dt.check_model(anchor)
anchor_images = [quotient(anchor, i) for i in range(3)]
fixture_points = [E(L(int(row["x"])), L(int(row["y"])))
                  for row in fixture["points"]]
assert len(fixture_points) == 16
old_keys = {key(q) for q in fixture_points} | {key(-q) for q in fixture_points}

controls = [row for row in plane_receipt["known_points_verified"]
            if row["source"] == "base"]
assert len(controls) == 16
for row in controls:
    i = row["index"]
    c = [F(int(a)) for a in row["coordinates_Fp_decimal"]]
    point = inverse(c[0] / c[2], c[1] / c[2])
    assert point[0] == L(F(int(base_receipt["base_selection"][i]["x"])))
    images = [quotient(point, j) for j in range(3)]
    projected = [h * (images[j] - anchor_images[j]) for j in range(3)]
    assert projected[0] == projected[1] == projected[2] == fixture_points[i]

base_rows = [tuple(map(int, line.split())) for line in
             (HERE / "type_i_line_base_points.tsv").read_text().splitlines()]
new_rows = [tuple(map(int, line.split())) for line in
            (HERE / "type_i_line_new_points.tsv").read_text().splitlines()]
assert len(base_rows) == 141 and len(new_rows) == 9784
base_points = []
for _, u, v in base_rows:
    point = inverse(F(u), F(v))
    images = [quotient(point, j) for j in range(3)]
    projected = [h * (images[j] - anchor_images[j]) for j in range(3)]
    assert projected[0] == projected[1] == projected[2]
    assert projected[0] != O and r * projected[0] == O
    base_points.append(projected[0])
base_keys = [key(q) for q in base_points]
base_folded = {min(key(q), key(-q)) for q in base_points}
assert len(set(base_keys)) == len(base_folded) == 141
assert not set(base_keys).intersection(old_keys)
pair_keys = [key(base_points[i] + base_points[j])
             for i in range(141) for j in range(i + 1, 141)]
assert len(pair_keys) == len(set(pair_keys)) == 9870

new_keys = []
seen = set()
seen_folded = set()
old_matches = base_matches = zeros = 0
for index, (u, v) in enumerate(new_rows):
    point = inverse(F(u), F(v))
    Q = h * (quotient(point, 0) - anchor_images[0])
    assert r * Q == O
    current = key(Q)
    folded = min(current, key(-Q))
    if index < 64:
        vector3 = [h * (quotient(point, j) - anchor_images[j])
                   for j in range(3)]
        assert vector3[0] == vector3[1] == vector3[2] == Q
    zeros += Q == O
    old_matches += current in old_keys
    base_matches += folded in base_folded
    assert current not in seen and folded not in seen_folded
    seen.add(current)
    seen_folded.add(folded)
    new_keys.append(current)
    if (index + 1) % 1000 == 0:
        print(f"LINE_PROJECTION_REPLAY_PROGRESS {index + 1}/{len(new_rows)}", flush=True)

assert zeros == old_matches == base_matches == 0
assert receipt["base_projected_distinct"] == len(set(base_keys)) == 141
assert receipt["distinct_base_pair_sums"] == len(set(pair_keys)) == 9870
assert receipt["new_projected_distinct"] == len(seen) == 9784
assert receipt["new_projected_sign_distinct"] == len(seen_folded) == 9784
assert receipt["new_zero"] == receipt["old16_signed_matches"] == 0
assert receipt["base141_signed_matches"] == 0
assert receipt["new_exact_repeats"] == receipt["new_signed_repeats"] == 0
assert receipt["base_key_sha256"] == hashlib.sha256(
    "\n".join(base_keys).encode()).hexdigest()
assert receipt["new_key_sha256"] == hashlib.sha256(
    "\n".join(new_keys).encode()).hexdigest()

out = {"status": "LINE_PROJECTION_INDEPENDENT_PASS",
       "base_points": 141, "new_points": 9784,
       "old_controls": 16, "three_component_new_controls": 64,
       "base_pair_sums": 9870,
       "receipt_sha256": sha(receipt_path),
       "source_sha256": sha(__file__)}
(HERE / "type_i_line_projection_verification.json").write_text(
    json.dumps(out, indent=2, sort_keys=True) + "\n")
print(json.dumps(out, sort_keys=True), flush=True)
