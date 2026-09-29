#!/usr/bin/env python3
"""Hash and arithmetic preflight only. This file never imports Sage or runs a map."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
PRIME = 263


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def matmul(a: list[list[int]], b: list[list[int]]) -> list[list[int]]:
    return [[sum(a[i][k] * b[k][j] for k in range(2)) % PRIME
             for j in range(2)] for i in range(2)]


def apply(a: list[list[int]], v: tuple[int, int]) -> tuple[int, int]:
    x = (a[0][0] * v[0] + a[0][1] * v[1]) % PRIME
    y = (a[1][0] * v[0] + a[1][1] * v[1]) % PRIME
    if x:
        return (1, y * pow(x, -1, PRIME) % PRIME)
    assert y
    return (0, 1)


def orbit(a: list[list[int]], seed: tuple[int, int]) -> tuple[tuple[int, int], ...]:
    values = []
    point = seed
    while point not in values:
        values.append(point)
        point = apply(a, point)
    assert point == seed
    return tuple(values)


def check_arithmetic(spec: dict) -> dict:
    instance = spec["instance"]
    assert spec["schema"] == "ecc2k130-leaf7-metered-v1"
    assert spec["status"] == "HELD_PROTOCOL_ONLY"
    assert spec["release_main_head"] is None
    assert spec["release_main_head_role"] == "pre_release_main_ancestor_not_checkout"
    assert spec["prerequisite"]["no_accepted_certificate_yet"] is True
    assert instance["field_degree"] == 131
    assert instance["field_modulus_hex"] == "0x800000000000000000000000000002007"
    assert instance["tau_polynomial"] == "T^2 + T + 2"
    a, b = 0, 1
    for _ in range(1, 131):
        a, b = -2 * b, a - b
    assert [str(a), str(b)] == instance["tau_131_coefficients"]
    assert b % PRIME == 0 and a % PRIME == PRIME - 1
    r = int(instance["subgroup_order"])
    assert int(instance["full_group_order"]) == 4 * r
    assert int(instance["trace"]) == (1 << 131) + 1 - 4 * r
    assert int(instance["frobenius_order_conductor"]) == b
    assert instance["source_endomorphism_conductor"] == 1
    assert instance["source_endomorphism_discriminant"] == -7
    assert instance["leaf_endomorphism_conductor"] == PRIME
    assert instance["leaf_endomorphism_discriminant"] == -7 * PRIME * PRIME

    panel = spec["kernel_panel"]
    c0, c1 = panel["twist_frobenius_matrix_columns_mod263"]
    m = [[c0[0], c1[0]], [c0[1], c1[1]]]
    norm7 = [[(int(i == j) - 2 * m[i][j]) % PRIME for j in range(2)]
             for i in range(2)]
    wrong_sign = [[(int(i == j) + 2 * m[i][j]) % PRIME for j in range(2)]
                  for i in range(2)]
    assert norm7 == panel["norm7_matrix_rows_mod263"]
    assert matmul(norm7, norm7) == [[PRIME - 7, 0], [0, PRIME - 7]]
    assert matmul(norm7, m) == matmul(m, norm7)
    first = tuple(panel["first_line"])
    other = tuple(panel["other_orbit_line"])
    image = tuple(panel["norm7_image_line"])
    assert apply(norm7, first) == image
    assert apply(wrong_sign, first) != image

    first_orbit = orbit(m, first)
    other_orbit = orbit(m, other)
    assert len(first_orbit) == len(other_orbit) == 131
    assert set(first_orbit).isdisjoint(other_orbit)
    hits = [i for i, line in enumerate(other_orbit) if line == image]
    assert hits == [panel["other_orbit_conjugacy_exponent"]] == [23]
    assert all(apply(norm7, line) in other_orbit for line in first_orbit)
    assert all(apply(norm7, line) in first_orbit for line in other_orbit)
    lines = {(1, j) for j in range(PRIME)} | {(0, 1)}
    assert len(lines) == PRIME + 1
    remaining = set(lines)
    cycles = []
    while remaining:
        cycle = orbit(m, min(remaining))
        assert set(cycle) <= remaining
        remaining.difference_update(cycle)
        cycles.append(len(cycle))
    assert sorted(cycles) == panel["all_projective_line_cycle_lengths"]
    assert panel["descending_orbit_lengths"] == [131, 131]
    assert panel["horizontal_line_count"] == 2
    assert panel["descending_line_count"] == 262
    assert (pow(-7, (PRIME - 1) // 2, PRIME) == 1)

    point_panel = spec["point_panel"]
    labels = point_panel["controls"] + [
        f"target-{i}" for i in range(point_panel["target_count"])]
    assert labels == point_panel["labels"]
    assert point_panel["target_prefix"] == "leaf-seven-bridge-v1"
    coefficients = []
    for i in range(8):
        prefix = point_panel["target_prefix"]
        u = int.from_bytes(hashlib.sha256(f"{prefix}|{i}|u".encode()).digest(), "big") % r
        v = 1 + int.from_bytes(hashlib.sha256(f"{prefix}|{i}|v".encode()).digest(), "big") % (r - 1)
        assert 0 <= u < r and 1 <= v < r
        coefficients.append([str(u), str(v)])
    assert len({tuple(pair) for pair in coefficients}) == 8
    projection = hashlib.sha256(
        (json.dumps(coefficients, separators=(",", ":")) + "\n").encode()
    ).hexdigest()
    assert projection == point_panel["coefficient_pairs_sha256"]
    cost = spec["cost"]
    assert cost["unit"] == "Fq multiplication equivalent"
    assert cost["calibration_batches"] == 11
    assert cost["calibration_operands_per_batch"] == 10000
    assert cost["square_weight_label"] == "leaf-seven-cost-v1"
    assert len(cost["calibration_operand_pair_sha256"]) == 11
    for batch, expected in enumerate(cost["calibration_operand_pair_sha256"]):
        stream = hashlib.sha256()
        for i in range(cost["calibration_operands_per_batch"]):
            for suffix in ("a", "b"):
                label = f"{cost['square_weight_label']}|{batch}|{i}|{suffix}"
                h = hashlib.sha256(label.encode()).digest()
                value = 1 + int.from_bytes(h, "big") % ((1 << 131) - 1)
                stream.update(value.to_bytes(17, "big"))
        assert stream.hexdigest() == expected, f"calibration batch {batch} changed"
    assert cost["max_calibration_relative_mad"] == 0.05
    assert "first_leaf_panel_transport" in cost["bridge_phases"]
    assert "first_leaf_panel_transport" not in cost["common_phases"]
    producer = spec["metering_producer"]
    assert producer["status"] == "NOT_IMPLEMENTED"
    assert producer["sha256"] is None
    assert producer["introduction_commit"] is None
    assert all(path in spec["input_sha256"] for path in cost["field_primitive_paths"])
    assert producer["required_before_release"] is True
    assert not (REPO / producer["path"]).exists()
    assert cost["incremental_bridge_over_direct_at_most"] == 0.9
    assert cost["cold_bridge_over_direct_below"] == 1.0
    assert cost["inversion_internals_in_mul_and_sq"] is True
    assert cost["extension_arithmetic_reduced_to_Fq"] is True
    assert cost["no_unmetered_phases"] is True
    assert cost["ecdlp_or_pdp_claims"] is False
    assert spec["caps"] == {
        "child_wall_seconds": 300,
        "child_peak_rss_bytes": 2147483648,
    }
    return {
        "cycle_lengths": sorted(cycles),
        "norm7_image_line": list(image),
        "unique_conjugacy_exponent": hits[0],
        "coefficient_pairs_sha256": projection,
        "input_panel_size": len(labels),
        "calibration_batches_checked": 11,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--manifest-only", action="store_true",
        help="check arithmetic and freeze fields without unavailable source blobs")
    args = parser.parse_args()
    spec = json.loads((HERE / "FROZEN.json").read_text())
    projection = check_arithmetic(spec)
    if not args.manifest_only:
        for relative, expected in spec["input_sha256"].items():
            actual = digest(REPO / relative)
            assert actual == expected, f"{relative}: expected {expected}, got {actual}"
    print(json.dumps({
        "schema": "ecc2k130-leaf7-metered-preflight-v1",
        "status": "PASS",
        "source_hashes_checked": not args.manifest_only,
        "measurement_run": False,
        **projection,
    }, sort_keys=True))


if __name__ == "__main__":
    main()
