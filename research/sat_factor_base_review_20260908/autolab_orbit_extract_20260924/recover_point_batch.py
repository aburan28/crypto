#!/usr/bin/env python3
"""Replay point-only compact relations and recover synthetic target scalars."""

import argparse
from collections import defaultdict
import gzip
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess

HERE = Path(__file__).resolve().parent
BASE = HERE / "independent_replay_20260924_codex/base_header.jsonl.gz"
RANK = HERE / "cold_batch_rank.py"


def load_rank():
    spec = importlib.util.spec_from_file_location("cold_batch_rank", RANK)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify(training, point_batch):
    rank = load_rank()
    verifier = rank.load_verifier()
    with gzip.open(BASE, "rt") as stream:
        header = json.loads(next(stream))
        assert not stream.read()
    curve = verifier.Curve(header)
    by_point, reps, eigenvalue = rank.verify_orbit_labels(curve, header)
    by_x = defaultdict(list)
    for point in by_point:
        by_x[point[0]].append(point)
    train_manifest = json.loads((training / "manifest.json").read_text())
    train_report = json.loads((training / "validation.json").read_text())
    assert train_report["rank"] == train_report["columns"] == len(reps) == 220
    assert train_report["factor_base_log_solution_verified"]
    assert train_report["all_relations_independently_group_verified"]
    assert train_report["producer_stdout_sha256"] == sha(training / "producer.stdout.jsonl")
    assert train_manifest["base_gzip_sha256"] == sha(BASE)
    logs = train_report["factor_base_log_solution"]
    for rep, log in zip(reps, logs):
        assert curve.scalar(curve.generator, log) == rep
    manifest = json.loads((point_batch / "manifest.json").read_text())
    assert manifest["base_hash"] == header["base_hash"]
    # This archived point run used the exact source merged by PR #735. Later
    # telemetry-only edits must not rewrite its source identity. CI fetches the
    # immutable merge commit, and this check binds that blob to the manifest.
    archived_source = subprocess.check_output(
        ["git", "show", "d92439080c2d5c3a858f71abfc8e749e273125bb:examples/koblitz_s5_sat_instance.rs"],
        cwd=rank.REPO,
    )
    assert manifest["source_sha256"] == hashlib.sha256(archived_source).hexdigest()
    assert manifest["point_file_sha256"] == sha(point_batch / "target_points.jsonl")
    assert manifest["training_manifest_sha256"] == sha(training / "manifest.json")
    assert manifest["training_validation_sha256"] == sha(training / "validation.json")
    points = [tuple(json.loads(line)) for line in (point_batch / "target_points.jsonl").read_text().splitlines()]
    assert len(points) == len(set(points)) == len(manifest["validation_scalars"])
    receipt = json.loads((point_batch / "resource_receipt.json").read_text())
    assert receipt["producer_stdout_sha256"] == sha(point_batch / "producer.stdout.jsonl")
    producer_stderr = point_batch / "producer.stderr.txt"
    stderr_digest = (sha(producer_stderr) if producer_stderr.exists()
                     else hashlib.sha256(b"").hexdigest())
    assert receipt["producer_stderr_sha256"] == stderr_digest
    assert receipt["targets_requested"] == len(points)
    assert len(receipt["negative_controls"]) == 3
    for control in receipt["negative_controls"]:
        assert control["returncode"] != 0
        stderr = point_batch / (control["case"] + ".stderr.txt")
        assert control["stderr_sha256"] == sha(stderr)
        assert control["expected_rejection"] in stderr.read_text()
        assert json.loads((point_batch / (control["case"] + ".jsonl")).read_text()) == control["point"]
    observation = json.loads((point_batch / "producer.stdout.jsonl").read_text())
    batch = observation["compact_orbit_point_batch"]
    assert observation["factor_base_input_hash"] == header["base_hash"]
    assert observation["models_examined"] == 0
    assert batch["targets_requested"] == len(points)
    assert batch["targets_extracted"] + len(batch["failed_target_points"]) == len(points)
    assert batch["sat_verification_included"] is False
    assert batch["sat_verification_wall_ms"] is None
    assert len(batch["relations"]) == len(points), "the fixed target batch has an extraction failure"
    assert [tuple(row["target_point"]) for row in batch["relations"]] == points
    rows = []
    for relation, target, validator in zip(batch["relations"], points, manifest["validation_scalars"]):
        assert curve.on_curve(target)
        assert curve.scalar(target, curve.order) is None
        lift = verifier.check_witness(curve, by_x, target,
                                      relation["x_codes"], relation["pinned_intermediates"])
        recovered = sum(by_point[tuple(point)][1] * logs[by_point[tuple(point)][0]]
                        for point in lift) % curve.order
        assert curve.scalar(curve.generator, recovered) == target
        assert recovered == validator, "held-out validator scalar disagrees"
        rows.append({
            "target_point": list(target), "recovered_scalar": recovered,
            "validator_scalar": validator, "query_ms": relation["query_ms"],
            "independent_group_lift": lift,
        })
    report = {
        "schema_version": "1.0", "classification": "PUBLIC_SYNTHETIC_POINT_INPUT_LOG_RECOVERY",
        "scope": "Eight held-out point-only inputs recovered using previously solved n53 factor-base logs; cold IC setup is charged, but no matched rho or n131 transfer claim",
        "targets_requested": len(points), "targets_recovered": len(rows),
        "independent_orbit_labels_verified": len(by_point),
        "training_full_rank_at_extracted": train_report["full_rank_at_extracted"],
        "training_process_wall_ms": train_report["process_wall_ms"],
        "point_batch_reported_stage_ms": batch["charged_total_ms"],
        "point_batch_query_ms_sum": batch["query_ms_sum"],
        "point_batch_stdout_sha256": sha(point_batch / "producer.stdout.jsonl"),
        "point_input_sha256": sha(point_batch / "target_points.jsonl"),
        "training_validation_sha256": sha(training / "validation.json"),
        "frobenius_eigenvalue": eigenvalue,
        "rows": rows,
    }
    (point_batch / "validation.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps({key: report[key] for key in (
        "classification", "targets_recovered", "training_full_rank_at_extracted",
        "training_process_wall_ms", "point_batch_query_ms_sum")}, sort_keys=True))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--training", type=Path, required=True)
    parser.add_argument("--point-batch", type=Path, required=True)
    args = parser.parse_args()
    verify(args.training, args.point_batch)
