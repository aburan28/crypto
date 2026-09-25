#!/usr/bin/env python3
"""Independent Fq replay of a frozen degree-7 bridge producer receipt.

This verifier does not import Sage or the producer. It checks the saved
bit-polynomial arithmetic, source α, both 263 quotients, cubic quotient,
target labels, orientation and all Fq full-point equalities. The extension
exceptional-point assertions remain producer/Sage evidence.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(REPO / "research/ecc2k130_relations"))
sys.path.insert(0, str(REPO / "research/ecc2k130_oriented_transport_20260924"))
from relations import Koblitz, challenge_curve  # noqa: E402
from oriented_velu import BinaryVeluMap  # noqa: E402
from cubic_velu import CubicVeluMap  # noqa: E402


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def unhex(point):
    return None if point is None else tuple(int(x, 16) for x in point)


def alpha(curve, point):
    return curve.add(point, curve.mul(curve.frobenius(point), 2))


def matrix_action(A, point):
    return [(A[0][0] * point[0] + A[0][1] * point[1]) % 263,
            (A[1][0] * point[0] + A[1][1] * point[1]) % 263]


def normalize(v):
    return [1, v[1] * pow(v[0], -1, 263) % 263] if v[0] else [0, 1]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    if args.out.exists():
        raise SystemExit("refusing to overwrite replay receipt")
    spec = json.loads((HERE / "FROZEN.json").read_text())
    source = json.loads(args.input.read_text())
    assert source["status"] == "PASS" and source["phase"] == "complete"
    assert source["freeze_sha256"] == sha(HERE / "FROZEN.json")
    assert source["input_sha256"] == spec["input_sha256"]
    assert source["implementation_sha256"] == spec["implementation_sha256"]
    for rel, expected in spec["input_sha256"].items():
        assert sha(REPO / rel) == expected, rel
    for name, expected in spec["implementation_sha256"].items():
        assert sha(HERE / name) == expected, name

    saved = spec["saved_kernel"]
    m0, m1 = saved["twist_frobenius_matrix_columns_mod263"]
    M = [[m0[0], m1[0]], [m0[1], m1[1]]]
    A = [[(int(i == j) - 2 * M[i][j]) % 263 for j in range(2)]
         for i in range(2)]
    assert A == saved["source_norm7_matrix_rows_mod263"]
    assert normalize(matrix_action(A, saved["source_line"])) == [1, 74]
    assert normalize(matrix_action(A, [1, 74])) == saved["source_line"]
    v = saved["other_orbit_line"]
    hits = []
    for k in range(131):
        if normalize(v) == [1, 74]:
            hits.append(k)
        v = matrix_action(M, v)
    assert hits == [23]
    assert source["matrix"]["conjugacy_exponents"] == hits

    F, E, P, Q, r = challenge_curve()
    torsion = json.loads((REPO / "research/ecc2k130_direction_review_20260924/"
                          "twist_torsion_results.json").read_text())
    row = next(row for row in torsion["runs"] if row["seed"] == saved["torsion_seed"])
    G, H = [tuple(int(x, 16) for x in point) for point in row["basis_on_twist"]]
    twist = Koblitz(F, a=1, b=1)
    G1 = twist.add(G, twist.neg(twist.mul(twist.frobenius(G), 2)))
    assert G1 == unhex(source["derived_generator"])
    assert twist.mul(G1, pow(155, -1, 263)) == twist.add(G, twist.mul(H, 74))
    phi0 = BinaryVeluMap.from_generator(E, twist, G, 263)
    phi1 = BinaryVeluMap.from_generator(E, twist, G1, 263)
    assert hex(phi0.codomain.b) == source["leaf0_b"]
    assert hex(phi1.codomain.b) == source["leaf1_b"]
    coeffs = tuple(int(x, 16) for x in source["kernel_polynomial_low_to_high"][:3])
    cubic = CubicVeluMap(phi0.codomain, coeffs)
    assert hex(cubic.codomain.b) == source["cubic_codomain_b"]
    assert source["division_factor_degrees"] == spec["gate"]["predicted_7_division_factor_degrees"]
    assert source["dual_composition_controls"] == spec["gate"]["predicted_dual_composition_controls"]
    s = int(source["codomain_model_isomorphism_s"], 16)
    assert F.sqr(s) ^ s == cubic.codomain.a ^ phi1.codomain.a
    assert cubic.codomain.b == phi1.codomain.b

    labels = spec["public_target_labels"]
    assert len(source["cases"]) == 4 + labels["count"]
    controls = {"P": P, "Q": Q, "P+Q": E.add(P, Q), "2P": E.mul(P, 2)}
    checked = []
    for i, case in enumerate(source["cases"]):
        if i < 4:
            label = spec["gate"]["full_point_controls"][i]
            expected = controls[label]
            assert case["u"] is None and case["v"] is None
        else:
            j = i - 4
            label = f"target-{j}"
            prefix = labels["prefix"]
            u = int.from_bytes(hashlib.sha256(f"{prefix}|{j}|u".encode()).digest(), "big") % r
            v = 1 + int.from_bytes(hashlib.sha256(f"{prefix}|{j}|v".encode()).digest(), "big") % (r - 1)
            assert (case["u"], case["v"]) == (str(u), str(v))
            expected = E.add(E.mul(P, u), E.mul(Q, v))
        assert case["label"] == label and unhex(case["source"]) == expected
        image_alpha = alpha(E, expected)
        assert unhex(case["alpha_source"]) == image_alpha
        x0 = phi0(expected)
        direct = phi1(image_alpha)
        assert x0 == unhex(case["leaf0"])
        assert direct == unhex(case["direct_leaf1"])
        cubic_image = cubic(x0)
        assert cubic_image == unhex(case["cubic_leaf1"])
        target_image = (cubic_image[0], cubic_image[1] ^ F.mul(s, cubic_image[0]))
        assert target_image == direct == unhex(case["sage_oriented_leaf1"])
        checked.append(label)

    assert cubic(None) is None
    try:
        cubic((0, 0))
    except ValueError:
        pass
    else:
        raise AssertionError("off-curve input accepted")
    assert source["off_curve_rejected"] is True
    assert source["infinity_mapped_to_infinity"] is True
    for name in ("forward_kernel", "reverse_kernel"):
        exceptional = source[name]
        assert exceptional["root_count"] == 3
        assert exceptional["point_count"] == 6
        assert exceptional["exact_order_seven"] is True
        assert exceptional["q_frobenius_eigenvalue"] == 5
        assert exceptional["eigenvalue_checks"] == 6
        assert exceptional["all_mapped_to_infinity"] is True
    cost = source["cost"]
    assert cost["shared_263_torsion_reused_for_direct_second_leaf"] is True
    assert cost["direct_second_leaf_rederived_from_saved_basis"] is True
    assert cost["incremental_mul_equivalent_ratio"] is None
    assert cost["cold_mul_equivalent_ratio"] is None
    assert cost["promotion"] is False
    replay = {"schema": "ecc2k130-leaf7-bridge-replay-v1", "status": "PASS",
              "source_receipt_sha256": sha(args.input),
              "freeze_sha256": source["freeze_sha256"],
              "independent_Fq_case_labels": checked,
              "exceptional_extension_scope": "producer assertions recorded; no independent extension replay",
              "cost_promotion": False}
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(replay, indent=2, sort_keys=True) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
