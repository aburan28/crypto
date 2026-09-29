#!/usr/bin/env python3
"""Independent full or partial group-law replay of both native-cyclic L384 arms."""
from __future__ import annotations

import argparse
from collections import defaultdict
import hashlib
import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
ORBIT = ROOT / "autolab_orbit_extract_20260924"
ROTATION = ROOT / "autolab_n53_rank_rotation_20260925"
RHO_STUDY = ROOT / "autolab_matched_point_rho_n53_20260925"
BASE_HASH = "2de8ec46916999dfe5e1ac95e68ea15a3f001eec27e013fa96d9525719b4f684"
POINT_KEYS = ("factor_base_point_coordinates", "factor_base_point_labels", "factor_base_representatives")


def load(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def arrays_hash(header: dict) -> str:
    data = json.dumps({key: header[key] for key in POINT_KEYS},
                      sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode() + b"\n"
    return hashlib.sha256(data).hexdigest()


def inputs():
    points = [tuple(json.loads(line)) for line in (HERE / "points_L384.jsonl").read_text().splitlines()]
    scalars = [int(line) for line in (HERE / "validator_scalars_L384.txt").read_text().splitlines()]
    assert len(points) == len(scalars) == len(set(points)) == len(set(scalars)) == 384
    return points, scalars


def replay_training(panel: Path, arm: str, rank, verifier, policy):
    root = panel / arm / "training"
    raw_file = root / "producer.stdout.jsonl"
    if not raw_file.is_file():
        return None, {"status": "NO_TRAINING_RAW"}
    receipt_path = root / "resource_receipt.json"
    if not receipt_path.is_file() or json.loads(receipt_path.read_text())["returncode"] != 0:
        return None, {"status": "CENSORED_TRAINING_CHILD", "raw_sha256": sha(raw_file)}
    raw = json.loads(raw_file.read_text())
    assert raw["factor_base_input_hash"] is None and raw["factor_base_input_path"] is None
    generated = raw["compact_orbit_base_header"]
    assert generated["selection_mode"] == "ascending_x_v1"
    assert generated["scanned_x"] == 465
    assert arrays_hash(generated) == BASE_HASH
    assert (generated["orbit_columns"], generated["factor_base_points"]) == (220, 23320)
    header_path = root / "base_header.jsonl"
    header = json.loads(header_path.read_text()) if header_path.is_file() else generated
    assert arrays_hash(header) == BASE_HASH
    if header_path.is_file():
        assert header["base_hash"] == BASE_HASH
        assert header["field_x_values_scanned"] == 465
        assert all(header[key] == value for key, value in generated.items())
    curve = verifier.Curve(header)
    assert curve.on_curve(curve.generator)
    assert curve.scalar(curve.generator, curve.order) is None
    by_point, reps, lam = rank.verify_orbit_labels(curve, header)
    assert len(by_point) == 23320 and len(reps) == 220
    by_x = defaultdict(list)
    for point in by_point:
        by_x[point[0]].append(point)
    assert len(by_x) == 11660 and all(len(pair) == 2 for pair in by_x.values())
    schedule = [int(line) for line in (HERE / f"training_{arm}_scalars.txt").read_text().splitlines()]
    assert (root / "target_scalars.txt").read_bytes() == (HERE / f"training_{arm}_scalars.txt").read_bytes()
    point_path = ROTATION / "target_points.jsonl" if arm == "A" else HERE / "training_B_points.jsonl"
    expected_points = [tuple(json.loads(line)) for line in point_path.read_text().splitlines()]
    assert len(schedule) == len(expected_points) == 512
    batch = raw["compact_orbit_batch"]
    assert batch["regular_scan_policy"] == "target_cyclic_v1"
    assert batch["targets_requested"] == 512
    queries, relations, failures = (
        batch["query_observations"], batch["relations"], batch["failed_target_scalars"]
    )
    assert len(queries) == 512 and len(relations) + len(failures) == 512
    matrix = rank.Echelon(220, curve.order)
    relation_index = failure_index = 0
    gains = []
    full_at = None
    solution = None
    post_rank = 0
    for index, (scalar, target, query) in enumerate(zip(schedule, expected_points, queries), 1):
        assert query["scalar"] == scalar
        assert curve.scalar(curve.generator, scalar) == target
        assert query["regular_scan_start"] == policy.expected_start(target, "target_cyclic_v1", batch["regular_states"])
        assert 0 <= query["regular_keys_visited"] <= batch["regular_states"]
        assert 0 <= query["group_lift_attempts"] <= query["indexed_partner_hits"] <= query["partner_roots"] <= 2 * query["s3_calls"]
        if not query["hit"]:
            assert failures[failure_index] == scalar
            failure_index += 1
            continue
        relation = relations[relation_index]
        relation_index += 1
        assert relation["scalar"] == scalar
        lift = verifier.check_witness(curve, by_x, target,
                                      relation["x_codes"], relation["pinned_intermediates"])
        row = {}
        for point in map(tuple, lift):
            column, coefficient = by_point[point]
            row[column] = (row.get(column, 0) + coefficient) % curve.order
        if solution is not None:
            assert sum(value * solution[column] for column, value in row.items()) % curve.order == scalar
            post_rank += 1
        if matrix.insert(row, scalar):
            gains.append(index)
            if len(matrix.pivots) == 220 and full_at is None:
                full_at = index
                solution = matrix.solution()
    assert relation_index == len(relations) and failure_index == len(failures)
    assert matrix.solution() == solution
    if solution is not None:
        for rep, log in zip(reps, solution):
            assert curve.scalar(curve.generator, log) == rep
    operational_path = root / "operational_solution.json"
    if operational_path.is_file():
        operational = json.loads(operational_path.read_text())
        assert operational["rank"] == len(matrix.pivots)
        assert operational["relations_used"] == len(relations)
        assert operational["failed_target_scalars"] == failures
        assert operational["first_full_rank_at"] == full_at
        assert operational["rank_gain_indices"] == gains
        assert operational["post_rank_predictions"] == post_rank
        assert operational["base_log_solution"] == solution
    report = {
        "status": "REPLAYED_FULL_RANK" if solution is not None else "REPLAYED_RANK_DEFICIENT",
        "raw_sha256": sha(raw_file), "native_base_arrays_sha256": BASE_HASH,
        "header_sha256": sha(header_path) if header_path.is_file() else None,
        "orbit_labels_replayed": len(by_point), "relations_replayed": len(relations),
        "failed_training_targets": failures, "rank": len(matrix.pivots),
        "first_full_rank_target": full_at, "frobenius_eigenvalue": lam,
        "training_s3_calls": sum(item["s3_calls"] for item in queries),
    }
    return (curve, by_point, by_x, solution, verifier), report


def replay_point(panel: Path, arm: str, context, points, labels, policy):
    raw_file = panel / arm / "ic/ic.stdout.jsonl"
    if not raw_file.is_file():
        return {"status": "NO_POINT_RAW"}
    receipt = panel / arm / "ic/resource_receipt.json"
    if not receipt.is_file() or json.loads(receipt.read_text())["returncode"] != 0:
        return {"status": "CENSORED_POINT_CHILD", "raw_sha256": sha(raw_file)}
    curve, by_point, by_x, logs, verifier = context
    assert logs is not None
    raw = json.loads(raw_file.read_text())
    assert raw["factor_base_input_hash"] == BASE_HASH
    batch = raw["compact_orbit_point_batch"]
    assert batch["regular_scan_policy"] == "target_cyclic_v1"
    assert batch["targets_requested"] == 384
    queries, relations, failures = (
        batch["query_observations"], batch["relations"], batch["failed_target_points"]
    )
    assert len(queries) == 384 and len(relations) + len(failures) == 384
    result = []
    relation_index = failure_index = 0
    for index, (target, expected, query) in enumerate(zip(points, labels, queries)):
        assert tuple(query["target_point"]) == target
        assert curve.scalar(curve.generator, expected) == target
        assert query["regular_scan_start"] == policy.expected_start(target, "target_cyclic_v1", batch["regular_states"])
        assert 0 <= query["group_lift_attempts"] <= query["indexed_partner_hits"] <= query["partner_roots"] <= 2 * query["s3_calls"]
        if not query["hit"]:
            assert tuple(failures[failure_index]) == target
            failure_index += 1
            result.append(None)
            continue
        relation = relations[relation_index]
        relation_index += 1
        assert tuple(relation["target_point"]) == target
        lift = verifier.check_witness(curve, by_x, target,
                                      relation["x_codes"], relation["pinned_intermediates"])
        value = sum(by_point[tuple(point)][1] * logs[by_point[tuple(point)][0]]
                    for point in lift) % curve.order
        assert curve.scalar(curve.generator, value) == target and value == expected, index
        result.append(value)
    assert relation_index == len(relations) and failure_index == len(failures)
    operational_file = panel / arm / "ic/operational_recovery.json"
    if operational_file.is_file():
        operational = json.loads(operational_file.read_text())
        assert operational["recovered_scalars"] == result
        assert operational["count"] == len(relations)
    return {
        "status": "REPLAYED_ALL_384_LOGS" if not failures else "REPLAYED_PARTIAL_POINT_LOGS",
        "raw_sha256": sha(raw_file), "relations_replayed": len(relations),
        "failed_point_targets": failures,
        "recovered_scalars_sha256": hashlib.sha256(json.dumps(result, separators=(",", ":")).encode()).hexdigest(),
        "s3_calls": sum(item["s3_calls"] for item in queries),
    }


def replay_rho(panel: Path, points, labels):
    raw_file = panel / "rho/rho.stdout.jsonl"
    if not raw_file.is_file():
        return {"status": "NO_RHO_RAW"}
    receipt = panel / "rho/resource_receipt.json"
    if not receipt.is_file() or json.loads(receipt.read_text())["returncode"] != 0:
        return {"status": "CENSORED_RHO_CHILD", "raw_sha256": sha(raw_file)}
    math = load(RHO_STUDY / "analyze.py", "native_l384_rho_independent_math")
    rows = [json.loads(line) for line in raw_file.read_text().splitlines() if line]
    fixtures = [row for row in rows if row["kind"] == "rho_ks_batch_fixture"]
    summaries = [row for row in rows if row["kind"] == "rho_ks_batch_summary"]
    assert len(fixtures) == 384 and len(summaries) == 1
    summary = summaries[0]
    assert summary["all_verified"] and summary["target_source"] == "explicit_public_points"
    assert summary["n"] == 53 and summary["a"] == 0
    assert summary["quotient_mode"] == "signed_frobenius"
    assert summary["batch_seed"] == 531929
    assert summary["corpus"] == "n53-native-cyclic-L384-20260929-v1"
    assert summary["fixtures"] == 384 and summary["automorphism_size"] == 106
    assert summary["dp_bits"] == 4 and summary["precompute_walks"] == 0
    assert summary["precompute_steps"] == summary["precompute_table_entries"] == 0
    assert math.on_curve(math.GENERATOR) and math.scalar(math.GENERATOR, math.ORDER) is None
    for index, (row, target, expected) in enumerate(zip(fixtures, points, labels)):
        assert row["fixture_index"] == index and row["published_fixture_scalar"] is None
        assert row["target_source"] == "explicit_public_points"
        assert tuple(row["published_q"]) == target
        value = row["recovered_fixture_scalar"]
        assert 0 < value < math.ORDER
        assert math.scalar(math.GENERATOR, value) == target and value == expected
    return {"status": "REPLAYED_ALL_384_RHO_LOGS", "raw_sha256": sha(raw_file),
            "rho_scalars_replayed": 384, "rho_walk_steps": summary["total_walk_steps"],
            "rho_table_entries": summary["table_entries"], "rho_charges": summary["charges"]}


def run(panel: Path):
    rank = load(ORBIT / "cold_batch_rank.py", "native_l384_independent_rank")
    verifier = rank.load_verifier()
    policy = load(ROTATION / "audit.py", "native_l384_independent_policy")
    points, labels = inputs()
    arms = {}
    for arm in ("B", "A"):
        context, training = replay_training(panel, arm, rank, verifier, policy)
        point = (replay_point(panel, arm, context, points, labels, policy)
                 if context is not None and context[3] is not None else {"status": "NO_FULL_RANK_POINT_REPLAY"})
        arms[arm] = {"training": training, "point": point}
    rho = replay_rho(panel, points, labels)
    complete = (rho["status"] == "REPLAYED_ALL_384_RHO_LOGS"
                and all(arms[arm]["training"]["status"] == "REPLAYED_FULL_RANK"
                        and arms[arm]["point"]["status"] == "REPLAYED_ALL_384_LOGS"
                        for arm in ("A", "B")))
    return {
        "classification": "INDEPENDENT_FULL_REPLAY" if complete else "INDEPENDENT_PARTIAL_REPLAY",
        "points_sha256": sha(HERE / "points_L384.jsonl"),
        "validator_scalars_sha256": sha(HERE / "validator_scalars_L384.txt"),
        "native_base_arrays_sha256": BASE_HASH,
        "arms": arms, "rho": rho,
        "common_operation_unit": None,
        "n131_transfer": None,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--panel", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    report = run(args.panel)
    args.out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"classification": report["classification"],
                      "B_training_replayed": report["arms"]["B"]["training"].get("relations_replayed", 0),
                      "B_point_replayed": report["arms"]["B"]["point"].get("relations_replayed", 0)},
                     sort_keys=True))


if __name__ == "__main__":
    main()
