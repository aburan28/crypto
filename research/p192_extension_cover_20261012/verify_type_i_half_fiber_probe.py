"""Independent point and quotient replay of the unit-coefficient inverses."""

import hashlib
import json
from pathlib import Path

from sage.all import EllipticCurve, PolynomialRing, prod, vector
import type_i_divisor_transfer as dt
import verify_type_i_quotient_base as vb


HERE = Path(__file__).resolve().parent
receipt_path = HERE / "type_i_half_fiber_probe.json"
receipt = json.loads(receipt_path.read_text())
assert receipt["source_sha256"] == hashlib.sha256(
    (HERE / "type_i_half_fiber_probe.py").read_bytes()
).hexdigest()
assert receipt["base_receipt_sha256"] == hashlib.sha256(
    (HERE / "type_i_quotient_base_probe.json").read_bytes()
).hexdigest()
plane_path = HERE / "type_i_plane_quartic.json"
assert receipt["plane_receipt_sha256"] == hashlib.sha256(plane_path.read_bytes()).hexdigest()

p, r, L, F = dt.p, dt.r, dt.L, dt.gm.F
E = EllipticCurve(L, [L(dt.gm.a), L(dt.gm.b)])
O = E(0)
assert r.is_prime() and r % 2 == 1
half = (r + 1) // 2
plane = json.loads(plane_path.read_text())
R = PolynomialRing(F, 3, names=("Z0", "Z1", "Z2"))
Z = R.gens()
quartic = sum(F(int(c)) * prod(Z[i] ** e[i] for i in range(3))
              for key, c in plane["quartic_coefficients"].items()
              for e in [tuple(map(int, key.split(",")))])


def decode(encoded):
    return tuple(L([int(c) for c in coordinate]) for coordinate in encoded)


def quotient(point, index):
    x = point[0]
    y = dt.check_model(point)[index]
    rho = dt.inf ** (p ** index)
    c = dt.c0 ** (p ** index)
    d = dt.d0 ** (p ** index)
    scale = dt.kappa ** (p ** index)
    assert x != rho
    return E(c + d / (x - rho), scale * y / (x - rho) ** 2)


def check_plane(point):
    ys = dt.check_model(point)
    assert all(y != 0 for y in ys)
    coordinates = dt.Vandermonde.solve_right(vector(L, [1 / y for y in ys]))
    assert quartic(*coordinates) == 0


assert len(receipt["classes"]) == len(vb.base) == 16
for index, item in enumerate(receipt["classes"]):
    assert item["base_index"] == index
    assert item["half_fiber_split"] and item["target_fiber_split"]
    assert item["unit_divisor_verified"]
    q = vb.base[index][0]
    assert vb.base[index] == (q, q, q) and q != O and r * q == O
    h = half * q
    assert 2 * h == q
    negative_seed = [decode(point) for point in item["half_fiber"]]
    positive_seed = [decode(point) for point in item["target_fiber"]]
    assert len(negative_seed) == len(positive_seed) == 2
    assert all(quotient(point, 0) == h for point in negative_seed)
    assert all(quotient(point, 0) == q for point in positive_seed)
    totals = [O, O, O]
    for shift in range(3):
        positives = [tuple(a ** (p ** shift) for a in point) for point in positive_seed]
        negatives = [tuple(a ** (p ** shift) for a in point) for point in negative_seed]
        for point in positives + negatives:
            check_plane(point)
        row = [sum((quotient(point, k) for point in positives), O) -
               sum((quotient(point, k) for point in negatives), O)
               for k in range(3)]
        assert row == [q if k == shift else O for k in range(3)]
        totals = [totals[k] + row[k] for k in range(3)]
    assert tuple(totals) == vb.base[index]
    assert item["positive_points"] == item["negative_points"] == 6

out = {"status": "verified_type_i_half_fiber_unit_divisors",
       "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
       "producer_receipt_sha256": hashlib.sha256(receipt_path.read_bytes()).hexdigest(),
       "classes_replayed": 16,
       "unit_coefficient_divisors_verified": 16,
       "all_support_points_on_smooth_plane_quartic": True,
       "all_quotient_vectors_match": True,
       "reduced_degree_three_representatives": None}
(HERE / "type_i_half_fiber_verification.json").write_text(
    json.dumps(out, indent=2, sort_keys=True) + "\n"
)
print(json.dumps(out, indent=2, sort_keys=True))
