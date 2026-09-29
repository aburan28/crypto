#!/usr/bin/env python3
"""Independently replay frozen rho point logs and charge the cold IC archive."""

import argparse
import hashlib
import json
from pathlib import Path
import statistics

HERE = Path(__file__).resolve().parent
PROTOCOL = HERE / "protocol.json"
N = 53
ORDER = 21044858204113
GENERATOR = (198217578752339, 7929897206038174)
POLY = (1 << N) | (1 << 6) | (1 << 2) | (1 << 1) | 1
MASK = (1 << N) - 1


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def mul(left, right):
    assert 0 <= left <= MASK and 0 <= right <= MASK
    result = 0
    while right:
        if right & 1:
            result ^= left
        right >>= 1
        left <<= 1
        if left >> N:
            left ^= POLY
    return result


def square(value):
    return mul(value, value)


def power(base, exponent):
    answer = 1
    while exponent:
        if exponent & 1:
            answer = mul(answer, base)
        base = square(base)
        exponent >>= 1
    return answer


def inverse(value):
    assert value
    answer = power(value, (1 << N) - 2)
    assert mul(value, answer) == 1
    return answer


def on_curve(point):
    x, y = point
    return square(y) ^ mul(x, y) == mul(square(x), x) ^ 1


def add(left, right):
    if left is None:
        return right
    if right is None:
        return left
    x1, y1 = left
    x2, y2 = right
    if x1 == x2:
        if (y1 ^ y2) == x1:
            return None
        if x1 == 0:
            return None
        slope = x1 ^ mul(y1, inverse(x1))
        x3 = square(slope) ^ slope
        y3 = square(x1) ^ mul(slope ^ 1, x3)
        return x3, y3
    slope = mul(y1 ^ y2, inverse(x1 ^ x2))
    x3 = square(slope) ^ slope ^ x1 ^ x2
    y3 = mul(slope, x1 ^ x3) ^ x3 ^ y1
    return x3, y3


def scalar(point, number):
    result = None
    while number:
        if number & 1:
            result = add(result, point)
        point = add(point, point)
        number >>= 1
    return result


def load_run(seed, count, protocol, labels):
    out = HERE / "runs" / f"seed_{seed}_L{count}"
    manifest = json.loads((out / "manifest.json").read_text())
    receipt = json.loads((out / "resource_receipt.json").read_text())
    assert manifest["protocol_sha256"] == digest(PROTOCOL)
    assert manifest["source_sha256"] == protocol["algorithm"]["source_sha256"]
    assert manifest["point_file_sha256"] == protocol["input"][f"point_file_{count}_sha256"]
    assert manifest["seed"] == receipt["seed"] == seed
    assert manifest["targets"] == receipt["targets"] == count
    assert receipt["complete_producer_output"] and not receipt["timed_out"]
    raw_file = out / "rho.stdout.jsonl"
    assert digest(raw_file) == receipt["stdout_sha256"]
    assert digest(out / "rho.stderr.txt") == receipt["stderr_sha256"]
    rows = [json.loads(line) for line in raw_file.read_text().splitlines()]
    fixtures = [row for row in rows if row["kind"] == "rho_ks_batch_fixture"]
    summaries = [row for row in rows if row["kind"] == "rho_ks_batch_summary"]
    assert len(fixtures) == count and len(summaries) == 1
    summary = summaries[0]
    assert summary["n"] == N and summary["a"] == 0
    assert summary["quotient_mode"] == "signed_frobenius"
    assert summary["automorphism_size"] == 106
    assert summary["target_source"] == "explicit_public_points"
    assert summary["precompute_walks"] == 0
    assert summary["dp_bits"] == 4
    assert summary["all_verified"]
    points = [tuple(json.loads(line)) for line in
              (HERE / f"target_points_{count}.jsonl").read_text().splitlines()]
    validations = []
    for index, (row, point) in enumerate(zip(fixtures, points)):
        assert row["fixture_index"] == index
        assert row["published_fixture_scalar"] is None
        assert row["target_source"] == "explicit_public_points"
        assert tuple(row["published_q"]) == point
        assert all(0 <= value <= MASK for value in point)
        assert on_curve(point)
        assert scalar(point, ORDER) is None
        recovered = row["recovered_fixture_scalar"]
        assert 0 < recovered < ORDER
        assert scalar(GENERATOR, recovered) == point
        assert recovered == labels[index]
        validations.append({
            "index": index, "point": point, "recovered_scalar": recovered,
            "validator_scalar": labels[index], "independent_scalar_check": True,
        })
    return {
        "seed": seed, "targets": count,
        "manifest_sha256": digest(out / "manifest.json"),
        "resource_receipt_sha256": digest(out / "resource_receipt.json"),
        "stdout_sha256": receipt["stdout_sha256"],
        "wall_ms": receipt["process_wall_ms"],
        "user_cpu_ms": receipt["child_user_cpu_ms"],
        "system_cpu_ms": receipt["child_system_cpu_ms"],
        "peak_rss_raw": receipt["child_peak_rss_raw"],
        "rss_unit": receipt["rss_unit"],
        "load_average_before": receipt["load_average_before"],
        "load_average_after": receipt["load_average_after"],
        "charges": summary["charges"],
        "walk_steps": summary["total_walk_steps"],
        "table_entries": summary["table_entries"],
        "cross_target_solves": summary["cross_target_solves"],
        "recovered": validations,
    }


def run(args):
    protocol = json.loads(PROTOCOL.read_text())
    label_manifest = json.loads(args.labels_manifest.read_text())
    assert label_manifest["point_file_sha256"] == protocol["input"]["point_file_8_sha256"]
    labels = label_manifest["validation_scalars"]
    assert len(labels) == 8
    assert on_curve(GENERATOR) and scalar(GENERATOR, ORDER) is None
    blocks = [load_run(seed, count, protocol, labels)
              for seed, count in protocol["run_order"]]
    training_manifest = json.loads(args.ic_training_manifest.read_text())
    training = json.loads(args.ic_training_validation.read_text())
    point_receipt = json.loads(args.ic_point_resource.read_text())
    point_validation = json.loads(args.ic_point_validation.read_text())
    assert label_manifest["training_manifest_sha256"] == digest(args.ic_training_manifest)
    assert label_manifest["training_validation_sha256"] == digest(args.ic_training_validation)
    assert training_manifest["targets"] == 512
    assert training_manifest["factor_base_hash"] == label_manifest["base_hash"]
    assert point_validation["point_input_sha256"] == protocol["input"]["point_file_8_sha256"]
    assert point_validation["point_batch_stdout_sha256"] == point_receipt["producer_stdout_sha256"]
    assert len(point_validation["rows"]) == 8
    assert all(row["recovered_scalar"] == row["validator_scalar"] == labels[index]
               for index, row in enumerate(point_validation["rows"]))
    assert training["rank"] == 220
    assert training["all_relations_independently_group_verified"]
    assert training["full_rank_at_extracted"] == 461
    assert point_receipt["targets_extracted"] == 8
    assert point_receipt["targets_requested"] == 8
    ic_training_wall = training["process_wall_ms"]
    ic_point_wall = point_receipt["positive_process_wall_ms"]
    ic_total = ic_training_wall + ic_point_wall
    batch = [b for b in blocks if b["targets"] == 8]
    singles = [b for b in blocks if b["targets"] == 1]
    single_dir = HERE / "ic_point_1_raw"
    single_protocol_path = HERE / "ic_single_addendum_protocol.json"
    single_protocol = json.loads(single_protocol_path.read_text())
    single_manifest = json.loads((single_dir / "manifest.json").read_text())
    single_receipt = json.loads((single_dir / "resource_receipt.json").read_text())
    single_validation = json.loads((single_dir / "validation.json").read_text())
    assert single_protocol["point_sha256"] == protocol["input"]["point_file_1_sha256"]
    assert single_protocol["training_manifest_sha256"] == digest(args.ic_training_manifest)
    assert single_protocol["training_validation_sha256"] == digest(args.ic_training_validation)
    assert single_manifest["point_file_sha256"] == single_protocol["point_sha256"]
    assert single_manifest["source_sha256"] == single_protocol["point_source_sha256"]
    assert single_manifest["exe_sha256"] == single_protocol["point_executable_sha256"]
    assert digest(single_dir / "target_points.jsonl") == single_protocol["point_sha256"]
    assert digest(single_dir / "producer.stdout.jsonl") == single_receipt["producer_stdout_sha256"]
    assert single_receipt["targets_requested"] == single_receipt["targets_extracted"] == 1
    assert len(single_receipt["negative_controls"]) == 3
    assert len(single_validation["rows"]) == 1
    single_point = tuple(json.loads((single_dir / "target_points.jsonl").read_text()))
    assert single_point == tuple(single_protocol["point_coordinates"])
    assert single_point == tuple(batch[0]["recovered"][0]["point"])
    recovered_single = single_validation["rows"][0]["recovered_scalar"]
    assert recovered_single == single_validation["rows"][0]["validator_scalar"] == labels[0]
    assert scalar(GENERATOR, recovered_single) == single_point
    ic_single_point_wall = single_receipt["positive_process_wall_ms"]
    ic_single_total = ic_training_wall + ic_single_point_wall
    decision = {
        "L8": "cold matched-point no-crossover for this fixed eight-point stream"
        if all(b["wall_ms"] < ic_total for b in batch)
        else "inconclusive",
        "L1": "post-panel addendum: cold matched-point no-crossover for this fixed first point"
        if all(b["wall_ms"] < ic_single_total for b in singles)
        else "post-panel addendum: inconclusive",
        "condition": "Every recovered scalar independently verified; all three rho cold walls below the corresponding fully charged two-process IC wall for L1 and L8",
        "do_not_extrapolate": "No n131 verdict or crossover beyond the first and eight fixed public synthetic targets",
        "protocol_correction": "The original preregistered L1 rule used an eight-point IC numerator and was excluded. A separately preregistered post-panel one-point IC measurement supplies the valid L1 numerator.",
    }
    result = {
        "schema_version": "1.0",
        "protocol_sha256": digest(PROTOCOL),
        "rho_source_sha256": protocol["algorithm"]["source_sha256"],
        "target_points_8_sha256": protocol["input"]["point_file_8_sha256"],
        "target_points_1_sha256": protocol["input"]["point_file_1_sha256"],
        "label_manifest_sha256": digest(args.labels_manifest),
        "ic_training_manifest_sha256": digest(args.ic_training_manifest),
        "ic_training_validation_sha256": digest(args.ic_training_validation),
        "ic_training_source_sha256": training_manifest["producer_source_sha256"],
        "ic_training_executable_sha256": training_manifest["producer_executable_sha256"],
        "ic_point_source_sha256": label_manifest["source_sha256"],
        "ic_point_executable_sha256": label_manifest["exe_sha256"],
        "ic_point_validation_sha256": digest(args.ic_point_validation),
        "ic_point_resource_sha256": digest(args.ic_point_resource),
        "ic_training_wall_ms": ic_training_wall,
        "ic_point_wall_ms": ic_point_wall,
        "ic_inclusive_two_process_wall_ms": ic_total,
        "ic_training_rank": training["rank"],
        "ic_first_full_rank_relation": training["full_rank_at_extracted"],
        "ic_training_extra_holdout_relations": 51,
        "ic_single_addendum_protocol_sha256": digest(single_protocol_path),
        "ic_single_manifest_sha256": digest(single_dir / "manifest.json"),
        "ic_single_resource_receipt_sha256": digest(single_dir / "resource_receipt.json"),
        "ic_single_validation_sha256": digest(single_dir / "validation.json"),
        "ic_single_stdout_sha256": digest(single_dir / "producer.stdout.jsonl"),
        "ic_single_point_wall_ms": ic_single_point_wall,
        "ic_single_inclusive_two_process_wall_ms": ic_single_total,
        "ic_single_recovered_scalar": recovered_single,
        "rho_L8_wall_ms": [b["wall_ms"] for b in batch],
        "rho_L8_cpu_ms": [b["user_cpu_ms"] + b["system_cpu_ms"] for b in batch],
        "rho_L8_wall_median_ms": statistics.median(b["wall_ms"] for b in batch),
        "rho_L8_walk_steps": [b["walk_steps"] for b in batch],
        "rho_L1_wall_ms": [b["wall_ms"] for b in singles],
        "rho_L1_cpu_ms": [b["user_cpu_ms"] + b["system_cpu_ms"] for b in singles],
        "ic_over_rho_L8_wall_ratios": [ic_total / b["wall_ms"] for b in batch],
        "ic_over_rho_L1_wall_ratios": [ic_single_total / b["wall_ms"] for b in singles],
        "blocks": blocks,
        "decision": decision,
    }
    args.out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({key: result[key] for key in
                      ("rho_L8_wall_ms", "rho_L8_cpu_ms", "rho_L8_walk_steps",
                       "rho_L1_wall_ms", "ic_inclusive_two_process_wall_ms",
                       "ic_single_inclusive_two_process_wall_ms",
                       "ic_over_rho_L1_wall_ratios",
                       "ic_over_rho_L8_wall_ratios", "decision")}, sort_keys=True))


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    source = Path("research/sat_factor_base_review_20260908/autolab_orbit_extract_20260924")
    p.add_argument("--labels-manifest", type=Path,
                   default=source / "point_recovery_20260924_codex/manifest.json")
    p.add_argument("--ic-training-manifest", type=Path,
                   default=source / "cold_rank_20260924_codex/manifest.json")
    p.add_argument("--ic-training-validation", type=Path,
                   default=source / "cold_rank_20260924_codex/validation.json")
    p.add_argument("--ic-point-resource", type=Path,
                   default=source / "point_recovery_20260924_codex/resource_receipt.json")
    p.add_argument("--ic-point-validation", type=Path,
                   default=source / "point_recovery_20260924_codex/validation.json")
    p.add_argument("--out", type=Path, default=HERE / "analysis.json")
    run(p.parse_args())
