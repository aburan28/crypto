#!/usr/bin/env python3
"""Independent full group-law replay of the n53 combined L384 control."""
from __future__ import annotations

import argparse
from collections import defaultdict
import hashlib
import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ORBIT = HERE.parent / "autolab_orbit_extract_20260924"
SHARED = HERE.parent / "autolab_shared_log_n53_20260925"
RHO = HERE.parent / "autolab_matched_point_rho_n53_20260925"
BASE_HASH = "d859319015ea405fd18aee41b51396ce4edcab64ef66265d8edcdeb5e040eb71"
POINTS_HASH = "d5185187014a12516aeef29306b65d4d20864293bb2c60667a9684fa8e51ac97"
MANIFEST_HASH = "f1843670a169d65645bf83886aeffee2e60e25627faa0002763eb1e7c13d8362"
SEED = 531384
CORPUS = "n53-combined-L384-20260925-v1"


def load(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def inputs():
    manifest_file = SHARED / "TARGET_MANIFEST.json"
    assert sha(manifest_file) == MANIFEST_HASH
    manifest = json.loads(manifest_file.read_text())
    parts = []
    labels = []
    for block in manifest["blocks"]:
        point_file = SHARED / block["L128_file"]
        assert sha(point_file) == block["L128_sha256"]
        parts.append(point_file.read_bytes())
        labels += block["validator_scalars"]
    point_file = HERE / "points_L384.jsonl"
    assert point_file.read_bytes() == b"".join(parts)
    assert sha(point_file) == POINTS_HASH
    points = [tuple(json.loads(line)) for line in point_file.read_text().splitlines()]
    assert len(points) == len(set(points)) == len(labels) == len(set(labels)) == 384
    return points, labels


def training_replay(panel: Path):
    rank = load(ORBIT / "cold_batch_rank.py", "audit_independent_rank")
    verifier = rank.load_verifier()
    training = panel / "training"
    header_file = training / "base_header.jsonl"
    header = json.loads(header_file.read_text())
    assert header["base_hash"] == BASE_HASH
    curve = verifier.Curve(header)
    by_point, reps, lam = rank.verify_orbit_labels(curve, header)
    assert len(by_point) == 23320 and len(reps) == 220
    by_x = defaultdict(list)
    for point in by_point:
        by_x[point[0]].append(point)
    scalars = rank.target_schedule(512, curve.order)
    assert scalars == [int(line) for line in (training / "target_scalars.txt").read_text().splitlines()]
    observed = json.loads((training / "producer.stdout.jsonl").read_text())
    assert observed["factor_base_input_hash"] == BASE_HASH
    batch = observed["compact_orbit_batch"]
    assert batch["targets_requested"] == batch["targets_extracted"] == 512
    assert len(batch["relations"]) == len(batch["query_observations"]) == 512
    assert not batch["failed_target_scalars"]
    matrix = rank.Echelon(220, curve.order)
    gains = []
    first_full_rank = None
    solved = None
    post_rank = 0
    for index, (relation, scalar, query) in enumerate(zip(batch["relations"], scalars, batch["query_observations"]), 1):
        assert relation["scalar"] == query["scalar"] == scalar
        assert query["hit"] and 0 <= query["partner_roots"] <= 2 * query["s3_calls"]
        assert 0 <= query["group_lift_attempts"] <= query["indexed_partner_hits"] <= query["partner_roots"]
        target = curve.scalar(curve.generator, scalar)
        assert target is not None
        lift = verifier.check_witness(curve, by_x, target, relation["x_codes"], relation["pinned_intermediates"])
        row = {}
        for point in map(tuple, lift):
            column, coefficient = by_point[point]
            row[column] = (row.get(column, 0) + coefficient) % curve.order
        if solved is not None:
            assert sum(value * solved[column] for column, value in row.items()) % curve.order == scalar
            post_rank += 1
        if matrix.insert(row, scalar):
            gains.append(index)
            if len(matrix.pivots) == 220 and first_full_rank is None:
                first_full_rank = index
                solved = matrix.solution()
    assert matrix.solution() == solved and solved is not None
    for rep, log in zip(reps, solved):
        assert curve.scalar(curve.generator, log) == rep
    operational_file = training / "operational_solution.json"
    operational = json.loads(operational_file.read_text())
    assert operational["rank"] == 220
    assert operational["relations_used"] == 512
    assert operational["first_full_rank_at"] == first_full_rank
    assert operational["rank_gain_indices"] == gains
    assert operational["post_rank_predictions"] == post_rank
    assert operational["base_log_solution"] == solved
    return rank, verifier, curve, by_point, by_x, solved, {
        "training_relations_replayed": 512,
        "orbit_labels_replayed": len(by_point),
        "rank": len(matrix.pivots),
        "first_full_rank_at": first_full_rank,
        "post_rank_predictions_replayed": post_rank,
        "training_s3_calls": sum(item["s3_calls"] for item in batch["query_observations"]),
        "training_index_entries": batch["index_entries"],
        "training_regular_states": batch["regular_states"],
        "operational_solution_sha256": sha(operational_file),
    }


def ic_replay(panel: Path, points, labels, verifier, curve, by_point, by_x, logs):
    raw_file = panel / "ic/ic.stdout.jsonl"
    observed = json.loads(raw_file.read_text())
    assert observed["factor_base_input_hash"] == BASE_HASH
    batch = observed["compact_orbit_point_batch"]
    assert batch["targets_requested"] == batch["targets_extracted"] == 384
    assert not batch["failed_target_points"]
    assert len(batch["relations"]) == len(batch["query_observations"]) == 384
    recovered = []
    for index, (relation, query, target, expected) in enumerate(
        zip(batch["relations"], batch["query_observations"], points, labels)
    ):
        assert tuple(relation["target_point"]) == tuple(query["target_point"]) == target, index
        assert query["hit"] and 0 <= query["partner_roots"] <= 2 * query["s3_calls"], index
        assert 0 <= query["group_lift_attempts"] <= query["indexed_partner_hits"] <= query["partner_roots"], index
        assert curve.on_curve(target) and curve.scalar(target, curve.order) is None
        lift = verifier.check_witness(curve, by_x, target, relation["x_codes"], relation["pinned_intermediates"])
        value = sum(by_point[tuple(point)][1] * logs[by_point[tuple(point)][0]] for point in lift) % curve.order
        assert curve.scalar(curve.generator, value) == target and value == expected, index
        recovered.append(value)
    operational_file = panel / "ic/operational_recovery.json"
    operational = json.loads(operational_file.read_text())
    assert operational["count"] == 384 and operational["recovered_scalars"] == recovered
    return recovered, {
        "ic_relations_replayed": 384,
        "ic_s3_calls": sum(item["s3_calls"] for item in batch["query_observations"]),
        "ic_partner_roots": sum(item["partner_roots"] for item in batch["query_observations"]),
        "ic_index_entries": batch["index_entries"],
        "ic_regular_states": batch["regular_states"],
        "operational_recovery_sha256": sha(operational_file),
    }


def rho_replay(panel: Path, points, labels):
    math = load(RHO / "analyze.py", "audit_rho_independent_math")
    raw_file = panel / "rho/rho.stdout.jsonl"
    rows = [json.loads(line) for line in raw_file.read_text().splitlines() if line]
    fixtures = [row for row in rows if row["kind"] == "rho_ks_batch_fixture"]
    summaries = [row for row in rows if row["kind"] == "rho_ks_batch_summary"]
    assert len(fixtures) == 384 and len(summaries) == 1
    summary = summaries[0]
    assert summary["all_verified"] and summary["target_source"] == "explicit_public_points"
    assert summary["n"] == 53 and summary["a"] == 0
    assert summary["quotient_mode"] == "signed_frobenius"
    assert summary["batch_seed"] == SEED and summary["corpus"] == CORPUS
    assert summary["fixtures"] == 384 and summary["automorphism_size"] == 106
    assert summary["dp_bits"] == 4 and summary["precompute_walks"] == 0
    assert summary["precompute_steps"] == summary["precompute_table_entries"] == 0
    assert all(isinstance(value, int) and value >= 0 for value in summary["charges"].values())
    assert math.on_curve(math.GENERATOR) and math.scalar(math.GENERATOR, math.ORDER) is None
    recovered = []
    for index, (row, target, expected) in enumerate(zip(fixtures, points, labels)):
        assert row["fixture_index"] == index and row["published_fixture_scalar"] is None
        assert row["target_source"] == "explicit_public_points"
        assert tuple(row["published_q"]) == target
        assert math.on_curve(target) and math.scalar(target, math.ORDER) is None
        value = row["recovered_fixture_scalar"]
        assert 0 < value < math.ORDER
        assert math.scalar(math.GENERATOR, value) == target and value == expected, index
        recovered.append(value)
    return recovered, {
        "rho_scalars_replayed": 384,
        "rho_walk_steps": summary["total_walk_steps"],
        "rho_table_entries": summary["table_entries"],
        "rho_cross_target_solves": summary["cross_target_solves"],
        "rho_charges": summary["charges"],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--panel", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    points, labels = inputs()
    rank, verifier, curve, by_point, by_x, logs, training = training_replay(args.panel)
    ic_logs, ic = ic_replay(args.panel, points, labels, verifier, curve, by_point, by_x, logs)
    rho_logs, rho = rho_replay(args.panel, points, labels)
    assert ic_logs == rho_logs == labels
    report = {
        "classification": "PASS_ALL_384_SAME_Q_LOGS_AND_512_TRAINING_RELATIONS",
        "points_sha256": POINTS_HASH,
        "validator_manifest_sha256": MANIFEST_HASH,
        "training": training, "ic": ic, "rho": rho,
        "all_operational_and_rho_logs_match_withheld_labels": True,
        "input_count": 384, "unique_points": len(set(points)),
    }
    args.out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"classification": report["classification"], "input_count": 384}, sort_keys=True))


if __name__ == "__main__":
    main()
