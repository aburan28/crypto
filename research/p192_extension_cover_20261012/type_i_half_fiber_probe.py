"""Test unit-coefficient inverse divisors using elliptic halves of quotient classes."""

import hashlib
import json
from pathlib import Path

from sage.all import EllipticCurve, prod, vector
import type_i_divisor_transfer as dt
import verify_type_i_quotient_base as vb

HERE = Path(__file__).resolve().parent
p, r, L, F = dt.p, dt.r, dt.L, dt.gm.F
E = EllipticCurve(L, [L(dt.gm.a), L(dt.gm.b)])
O = E(0)
half = (r + 1) // 2
plane = json.loads((HERE / "type_i_plane_quartic.json").read_text())
from sage.all import PolynomialRing
R = PolynomialRing(F, 3, names=("Z0", "Z1", "Z2"))
Z = R.gens()
quartic = sum(F(int(c)) * prod(Z[i] ** e[i] for i in range(3))
              for key, c in plane["quartic_coefficients"].items()
              for e in [tuple(map(int, key.split(",")))])


def frobenius(point, power):
    return tuple(a ** (p ** power) for a in point)


def quotient(point, index):
    x = point[0]
    y = dt.check_model(point)[index]
    rho, c, d, kappa = [a ** (p ** index)
                        for a in (dt.inf, dt.c0, dt.d0, dt.kappa)]
    return E(c + d / (x - rho), kappa * y / (x - rho) ** 2)


def plane_check(point):
    ys = dt.check_model(point)
    assert all(y != 0 for y in ys)
    coordinate = dt.Vandermonde.solve_right(vector(L, [1 / y for y in ys]))
    assert quartic(*coordinate) == 0


def encode(point):
    return [[str(coefficient) for coefficient in coordinate.polynomial().list()]
            for coordinate in point]


records = []
for index, target in enumerate(vb.base):
    q = target[0]
    assert target == (q, q, q) and r * q == O
    h = half * q
    assert 2 * h == q
    record = {"base_index": index, "half_fiber_split": False,
              "target_fiber_split": False, "unit_divisor_verified": False}
    try:
        positive_seed = dt.fiber(q)
        record["target_fiber_split"] = len(positive_seed) == 2
        negative_seed = dt.fiber(h)
        record["half_fiber_split"] = len(negative_seed) == 2
        positive = [frobenius(point, i) for i in range(3)
                    for point in positive_seed]
        negative = [frobenius(point, i) for i in range(3)
                    for point in negative_seed]
        for point in positive + negative:
            plane_check(point)
        image = tuple(sum((quotient(point, i) for point in positive), O) -
                      sum((quotient(point, i) for point in negative), O)
                      for i in range(3))
        assert image == target
        record["unit_divisor_verified"] = True
        record["positive_points"] = len(positive)
        record["negative_points"] = len(negative)
        record["half_fiber"] = [encode(point) for point in negative_seed]
        record["target_fiber"] = [encode(point) for point in positive_seed]
    except Exception as error:
        record["error_type"] = type(error).__name__
        record["error"] = str(error).splitlines()[0]
    records.append(record)
    print(json.dumps(record, sort_keys=True), flush=True)

out = {"status": "bounded_type_i_half_fiber_probe",
       "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
       "base_receipt_sha256": hashlib.sha256((HERE / "type_i_quotient_base_probe.json").read_bytes()).hexdigest(),
       "plane_receipt_sha256": hashlib.sha256((HERE / "type_i_plane_quartic.json").read_bytes()).hexdigest(),
       "classes": records,
       "unit_coefficient_divisors_verified": sum(x["unit_divisor_verified"] for x in records),
       "reduced_degree_three_representatives": None,
       "boundary": "Unit-coefficient degree-six positive and negative formal divisors when half fibers split; degree-three Jacobian reduction remains open."}
(HERE / "type_i_half_fiber_probe.json").write_text(json.dumps(out, indent=2, sort_keys=True) + "\n")
print(json.dumps({k: v for k, v in out.items() if k != "classes"}, sort_keys=True))
