#!/usr/bin/env sage -python
"""Validate the selected nondegenerate n=31 GGMP cell as curve points."""

from __future__ import annotations

import itertools
import json
from pathlib import Path

from sage.all import EllipticCurve, GF, PolynomialRing


ROOT = Path(__file__).resolve().parent
RUN = ROOT / "stage-4-selected-ggmp-20260909"
CELL = RUN / "n31-l5-m3-ggmp-a0-f0"


def point(curve, field, record):
    return curve(field.from_integer(int(record["x"])), field.from_integer(int(record["y"])))


def main() -> None:
    report = json.loads((RUN / "result.json").read_text())
    row = report["instances"][0]
    manifest = row["manifest"]
    assert row["cell"] == {
        "n": 31,
        "ell": 5,
        "m": 3,
        "basis": "ggmp",
        "curve_a": 0,
        "factor_index": 0,
    }
    assert manifest["factor_base_geometry"]["distinct_curve_points"] == 63
    assert manifest["factor_base_predicate"]["factor_bitmask"] == 37
    assert manifest["factor_base_predicate"]["linearised_exponents"] == [0, 2, 5]
    assert not manifest["factor_base_predicate"]["uses_discrete_log_labels"]
    assert not manifest["factor_base_predicate"]["enumerates_target_subgroup"]

    polynomial = PolynomialRing(GF(2), "z")
    z = polynomial.gen()
    field = GF(2**31, "a", modulus=z**31 + z**3 + 1)
    curve = EllipticCurve(field, [1, 0, 0, 0, 1])
    identity = curve(0)
    basis = [field.from_integer(int(value)) for value in manifest["factor_base_basis_bitmasks"]]
    factor_domain = [field(0)]
    for value in basis:
        factor_domain += [old + value for old in factor_domain]
    lifts = {int(x.to_integer()): curve.lift_x(x, all=True) for x in factor_domain}
    assert sum(len(points) for points in lifts.values()) == 63
    target = point(curve, field, manifest["target"])
    planted = [point(curve, field, value) for value in manifest["planted_points"]]
    assert sum(planted, identity) == target

    cms = next(value for value in row["solvers"] if value["solver"] == "cryptominisat")
    assert cms["status"] == "sat" and cms["source_model_valid"] is True
    values = {
        abs(int(token)): int(token) > 0
        for line in (CELL / "cryptominisat.stdout").read_text().splitlines()
        if line.startswith("v ")
        for token in line.split()[1:]
        if token != "0"
    }
    max_variable = manifest["exports"]["cryptominisat_xor_dimacs"]["variables"]
    assert all(index in values for index in range(1, max_variable + 1))
    xs = []
    ell = 5
    for block in range(3):
        x = field(0)
        for offset, basis_value in enumerate(basis):
            if values[block * ell + offset + 1]:
                x += basis_value
        xs.append(x)
    choices = [lifts[int(x.to_integer())] for x in xs]
    witness = next(
        (candidate for candidate in itertools.product(*choices) if sum(candidate, identity) == target),
        None,
    )
    assert witness is not None

    direct_points = [candidate for points in lifts.values() for candidate in points]
    direct_witness = next(
        (
            candidate
            for candidate in itertools.product(direct_points, repeat=3)
            if sum(candidate, identity) == target
        ),
        None,
    )
    assert direct_witness is not None

    discovery = json.loads((ROOT / "stage-4-discovery-summary.json").read_text())
    selected = discovery["selected"]
    assert (selected["curve_a"], selected["factor_index"], selected["distinct_curve_points"]) == (0, 0, 63)
    output = {
        "schema": "koblitz_ggmp_selected_cell_validation.v1",
        "selection_reproduced": True,
        "factor_domain_size": len(factor_domain),
        "factor_curve_points": sum(len(points) for points in lifts.values()),
        "target": manifest["target"],
        "planted_point_sum_valid": True,
        "cryptominisat": {
            "source_model_valid": True,
            "x_coordinates": [int(x.to_integer()) for x in xs],
            "rational_lifts_exist": [bool(choice) for choice in choices],
            "signed_point_sum_valid": True,
            "witness": [[int(value[0].to_integer()), int(value[1].to_integer())] for value in witness],
        },
        "direct_mitm_point_decomposition_reproduced": True,
        "native_sat": manifest["native_sat"]["result"],
        "wdsat": next(value for value in row["solvers"] if value["solver"] == "wdsat")["status"],
        "claim_boundary": "one nondegenerate public-parameter GGMP cell; no scaling or SOTA claim",
    }
    (ROOT / "stage-4-selected-validation.json").write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
