#!/usr/bin/env python3
"""Independently recover a compact point batch or replay matched rho scalars."""
from __future__ import annotations

import argparse
from collections import defaultdict
import hashlib
import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
ORBIT = REPO / "research/sat_factor_base_review_20260908/autolab_orbit_extract_20260924"
RHO = REPO / "research/sat_factor_base_review_20260908/autolab_matched_point_rho_n53_20260925"


def load(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def target_info(block: int, count: int):
    manifest = json.loads((HERE / "TARGET_MANIFEST.json").read_text())
    row = manifest["blocks"][block]
    points_file = HERE / row[f"L{count}_file"]
    assert sha(points_file) == row[f"L{count}_sha256"]
    points = [tuple(json.loads(line)) for line in points_file.read_text().splitlines()]
    assert len(points) == count == len(set(points))
    labels = row["validator_scalars"][:count]
    assert len(labels) == count
    return points_file, points, labels


def verify_ic(args, points_file, points, labels):
    rank = load(ORBIT / "cold_batch_rank.py", "cold_batch_rank_shared_log")
    verifier = rank.load_verifier()
    header = json.loads((args.training / "base_header.jsonl").read_text())
    train_report = json.loads((args.training / "validation.json").read_text())
    assert train_report["rank"] == train_report["columns"] == 220
    assert train_report["factor_base_log_solution_verified"]
    assert train_report["all_orbit_labels_independently_verified"]
    assert train_report["all_relations_independently_group_verified"]
    assert train_report["producer_stdout_sha256"] == sha(args.training / "producer.stdout.jsonl")
    assert header["base_hash"] == rank.BASE_HASH
    curve = verifier.Curve(header)
    by_point, reps, _ = rank.verify_orbit_labels(curve, header)
    logs = train_report["factor_base_log_solution"]
    assert len(reps) == len(logs) == 220
    by_x = defaultdict(list)
    for point in by_point:
        by_x[point[0]].append(point)
    raw_file = args.raw / "ic.stdout.jsonl"
    observed = json.loads(raw_file.read_text())
    batch = observed["compact_orbit_point_batch"]
    assert observed["factor_base_input_hash"] == rank.BASE_HASH
    assert batch["targets_requested"] == args.count
    assert batch["targets_extracted"] == args.count
    assert not batch["failed_target_points"]
    assert len(batch["relations"]) == args.count
    assert batch["sat_verification_included"] is False
    recovered = []
    for index, (relation, target, expected) in enumerate(zip(batch["relations"], points, labels)):
        assert tuple(relation["target_point"]) == target, index
        assert curve.on_curve(target) and curve.scalar(target, curve.order) is None
        lift = verifier.check_witness(curve, by_x, target,
                                      relation["x_codes"], relation["pinned_intermediates"])
        value = sum(by_point[tuple(point)][1] * logs[by_point[tuple(point)][0]]
                    for point in lift) % curve.order
        assert curve.scalar(curve.generator, value) == target, index
        assert value == expected, index
        recovered.append(value)
    observations = batch["query_observations"]
    assert len(observations) == args.count
    for index, (item, target) in enumerate(zip(observations, points)):
        assert tuple(item["target_point"]) == target, index
        assert item["hit"], index
        assert item["s3_calls"] >= 0
        assert 0 <= item["group_lift_attempts"] <= item["indexed_partner_hits"] <= item["partner_roots"]
        assert item["partner_roots"] <= 2 * item["s3_calls"]
    return {
        "arm": "ic", "n": 53, "count": args.count, "block": args.block,
        "raw_sha256": sha(raw_file), "points_sha256": sha(points_file),
        "training_validation_sha256": sha(args.training / "validation.json"),
        "base_hash": rank.BASE_HASH, "recovered_scalars": recovered,
        "orbit_labels_checked": len(by_point),
        "all_relations_group_replayed": True,
        "index_entries": batch["index_entries"],
        "regular_states": batch["regular_states"],
        "query_ms_sum": batch["query_ms_sum"],
        "batch_loop_wall_ms": batch["batch_loop_wall_ms"],
        "query_observations": observations,
        "partner_trials_sum": sum(r["trials"] for r in batch["relations"]),
    }


def verify_rho(args, points_file, points, labels):
    math = load(RHO / "analyze.py", "matched_rho_independent_math")
    raw_file = args.raw / "rho.stdout.jsonl"
    rows = [json.loads(line) for line in raw_file.read_text().splitlines() if line]
    fixtures = [item for item in rows if item["kind"] == "rho_ks_batch_fixture"]
    summaries = [item for item in rows if item["kind"] == "rho_ks_batch_summary"]
    assert len(fixtures) == args.count and len(summaries) == 1
    summary = summaries[0]
    assert summary["all_verified"] and summary["target_source"] == "explicit_public_points"
    assert summary["n"] == 53 and summary["a"] == 0
    assert summary["quotient_mode"] == "signed_frobenius"
    assert summary["batch_seed"] == 531320 + args.block
    assert summary["corpus"] == f"n53-shared-log-20260925-b{args.block}-L{args.count}"
    assert summary["fixtures"] == args.count
    assert summary["automorphism_size"] == 106
    assert summary["dp_bits"] == 4 and summary["precompute_walks"] == 0
    assert summary["precompute_steps"] == 0 and summary["precompute_table_entries"] == 0
    assert all(isinstance(value, int) and value >= 0 for value in summary["charges"].values())
    assert math.on_curve(math.GENERATOR) and math.scalar(math.GENERATOR, math.ORDER) is None
    recovered = []
    for index, (row, target, expected) in enumerate(zip(fixtures, points, labels)):
        assert row["fixture_index"] == index
        assert row["published_fixture_scalar"] is None
        assert row["target_source"] == "explicit_public_points"
        assert tuple(row["published_q"]) == target
        assert math.on_curve(target) and math.scalar(target, math.ORDER) is None
        value = row["recovered_fixture_scalar"]
        assert 0 < value < math.ORDER
        assert math.scalar(math.GENERATOR, value) == target
        assert value == expected
        recovered.append(value)
    return {
        "arm": "rho", "n": 53, "count": args.count, "block": args.block,
        "raw_sha256": sha(raw_file), "points_sha256": sha(points_file),
        "recovered_scalars": recovered,
        "walk_steps": summary["total_walk_steps"],
        "table_entries": summary["table_entries"],
        "cross_target_solves": summary["cross_target_solves"],
        "charges": summary["charges"],
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--arm", choices=("ic", "rho"), required=True)
    parser.add_argument("--training", type=Path, required=True)
    parser.add_argument("--block", type=int, choices=(0, 1, 2), required=True)
    parser.add_argument("--count", type=int, choices=(32, 128), required=True)
    parser.add_argument("--raw", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    points_file, points, labels = target_info(args.block, args.count)
    report = (verify_ic if args.arm == "ic" else verify_rho)(args, points_file, points, labels)
    report["verdict"] = "PASS"
    args.out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"arm": args.arm, "block": args.block, "count": args.count,
                      "verdict": "PASS"}, sort_keys=True))


if __name__ == "__main__":
    main()
