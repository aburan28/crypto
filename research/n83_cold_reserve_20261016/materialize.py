#!/usr/bin/env python3
"""Bind replayed four-point N83 reserve rows to exact IC1 identities."""

import argparse
import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
DYNAMIC_MATERIALIZE = ROOT / "research/n83_root11_dynamic_20261016/materialize.py"
RUNNER = HERE / "run.py"
SOURCE_COMMIT = "d67d73defd9ce1b092e1680790611fb87124b6a6"
SOURCE_SHA256 = "d094345474605075566d646ddecef8f694cdd644b0f9190a6b1ef1551fd06888"
RUN_ID = 38890911296
ARTIFACT_ID = 12045705683
ARCHIVE_SHA256 = "bcf6c55708a7e8026d760205810c47de780d97b2904d2b4a87310457b5e4b9c3"
ANALYSIS_SHA256 = "28253a17d18b6951843eca6451782e590b6d7d07c9b7a54ec96e9850bffea45d"

spec = importlib.util.spec_from_file_location("n83_dynamic_materialize", DYNAMIC_MATERIALIZE)
dynamic = importlib.util.module_from_spec(spec)
spec.loader.exec_module(dynamic)
panel = dynamic.panel
panel.SOURCE_COMMIT = SOURCE_COMMIT
panel.SOURCE_SHA256 = SOURCE_SHA256
rspec = importlib.util.spec_from_file_location("n83_reserve_runner", RUNNER)
runner = importlib.util.module_from_spec(rspec)
rspec.loader.exec_module(runner)


def write_once(path, value):
    blob = json.dumps(value, indent=2, sort_keys=True) + "\n"
    if path.exists() and path.read_text() != blob:
        raise ValueError(f"sealed file differs: {path}")
    path.write_text(blob)


def workload(index, point, answer, source, input_path):
    record = {
        "schema_version": 1, "curve_id": panel.CURVE_ID,
        "target_count": 1, "public_targets": [point],
        "input_law": {
            "domain": "n83-cold-factor-base-one-plus-four-20261016-v1",
            "fixture_index": index,
            "scalar_rule": "1+(little_endian(SHA256(domain||u32_le(index))) mod (r-1))",
            "seed_sha256": answer["seed_sha256"],
            "generator_sha256": source["hashes"][
                "research/n83_cold_factor_base_20261016/generate.py"],
            "answer_sidecar_sha256": panel.file_sha(
                ROOT / "research/n83_cold_factor_base_20261016/reserve_answers.jsonl"),
        },
        "target_sha256": panel.file_sha(input_path),
        "cold_or_warm": "one_public_point_per_fresh_process_empty_base_index_and_logs",
    }
    return panel.digest(record)[:12], record


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("results", type=Path)
    parser.add_argument("--archive", required=True, type=Path)
    args = parser.parse_args()
    results = args.results.resolve()
    if panel.file_sha(args.archive) != ARCHIVE_SHA256:
        raise ValueError("exact raw artifact ZIP hash differs")
    analysis = panel.read(results / "analysis.json")
    source = panel.read(results / "source.json")
    host = panel.read(results / "host.json")
    versioned_analysis = HERE / f"analysis_{RUN_ID}.json"
    if (panel.file_sha(results / "analysis.json") != ANALYSIS_SHA256
            or panel.file_sha(versioned_analysis) != ANALYSIS_SHA256
            or panel.read(versioned_analysis) != analysis
            or analysis["status"] != "BOUNDED_UNKNOWN"
            or analysis["evidence_status"] != "PASS_INDEPENDENT_REPLAY"
            or analysis["selection_admitted"] is not True
            or analysis["selected_pipeline"] is not None
            or analysis["selected_candidate_id"] is not None
            or host["git_sha"] != SOURCE_COMMIT
            or source["source_sha256"] != SOURCE_SHA256
            or source["hashes"]["examples/koblitz_orbit_dlp_fast.rs"] != SOURCE_SHA256
            or source["hashes"]["Cargo.lock"] != panel.LOCK_SHA256
            or source["schedule"] != runner.schedule()
            or source["runner_sha256"] != panel.file_sha(RUNNER)):
        raise ValueError("exact source, schedule or replay status differs")
    target_path = ROOT / "research/n83_cold_factor_base_20261016/reserve_targets.jsonl"
    answer_path = ROOT / "research/n83_cold_factor_base_20261016/reserve_answers.jsonl"
    points = [json.loads(line) for line in target_path.read_text().splitlines()]
    answers = [json.loads(line) for line in answer_path.read_text().splitlines()]
    if len(points) != 4 or len(answers) != 4:
        raise ValueError("four frozen reserve points and answers required")
    field, curve = panel.curve_record()
    manifests, workloads, rows = {}, {}, []
    for point_row in analysis["point_rows"]:
        index = point_row["point_index"]
        point, answer = points[index - 1], answers[index - 1]
        if (point_row["public_point"] != point
                or answer["point"] != point or answer["fixture_index"] != index):
            raise ValueError(f"point {index}: frozen identity differs")
        workload_id, record = workload(index, point, answer, source,
                                       results / "inputs" / f"point_{index}.jsonl")
        workloads[workload_id] = {"workload_id": workload_id, **record}
        for item in point_row["arms"]:
            arm = item["arm"]
            name = arm.split("_", 1)[1]
            pipeline = "k1024root22" if name.startswith("k1024root22_") else name
            k, layout = {
                "k1024root22": (1024, "root22_dynamic"),
                "k1152root22": (1152, "root22_dynamic"),
                "k1152root11": (1152, "root11_dynamic"),
                "k1536root22": (1536, "root22_dynamic"),
            }[pipeline]
            if item["status"] != "PASS_INDEPENDENT_REPLAY" or item["resources"]["status"] != "PASS":
                raise ValueError(f"{arm}: a complete exact replay is required")
            folder = results / "runs" / arm
            summary = panel.read(folder / "stdout.jsonl")
            manifest = dynamic.candidate(k, layout, summary,
                                         folder / "base_dump.jsonl", field, curve)
            if (item["base_hash"] != summary["base_hash"]
                    or item["actual_factor_base_points"] != 166 * k):
                raise ValueError(f"{arm}: replayed base count or hash differs")
            candidate_id = manifest["candidate_id"]
            if candidate_id in manifests and manifests[candidate_id] != manifest:
                raise ValueError(f"{arm}: candidate identity collision")
            manifests[candidate_id] = manifest
            repetition = 2 if arm.endswith("_b") else 1
            rows.append({
                "arm": arm, "candidate_id": candidate_id,
                "workload_id": workload_id,
                "run_id": candidate_id + "W" + workload_id + "R" + str(repetition),
                "status": item["status"], "hardware_tier": analysis["hardware_tier"],
                "source_commit": SOURCE_COMMIT,
                "executed_binary_sha256": source["binary_sha256"],
                "rho_binary_sha256": source["rho_binary_sha256"],
                "host_receipt_sha256": panel.file_sha(results / "host.json"),
                "source_receipt_sha256": panel.file_sha(results / "source.json"),
                "correctness_receipt_sha256": ANALYSIS_SHA256,
                "raw_artifact_run_id": RUN_ID, "raw_artifact_id": ARTIFACT_ID,
                "cold_process_ms": item["cold_process_ms_exploratory"],
                "one_target_online_ms": item["one_target_online_ms"],
                "online_phases_ms": item["online_phases_ms"],
                "setup_ms": item["setup_ms"], "rank_stage_ms": item["rank_stage_ms"],
                "targets_total_ms": item["targets_total_ms"],
                "rank_attempts": item["rank_attempts"],
                "rank_relations": item["rank_relations"],
                "target_probes": item["target_probes"],
                "actual_usable_points": item["actual_factor_base_points"],
                "effective_columns": item["effective_columns"],
                "peak_rss_bytes": item["resources"]["peak_rss_bytes"],
                "outer_wall_seconds": item["full_child_wall_seconds"],
            })
    output = {
        "schema": "n83.cold-reserve-ic1/v1", "source_commit": SOURCE_COMMIT,
        "raw_artifact_sha256": ARCHIVE_SHA256,
        "analysis_sha256": ANALYSIS_SHA256,
        "workloads": [workloads[key] for key in sorted(workloads)],
        "candidates": sorted(manifests), "rows": rows,
        "selection_status": analysis["status"],
        "selected_candidate_id": None,
        "aa_controls": [row["aa_control"] for row in analysis["point_rows"]],
        "paired_intervals": analysis["paired_intervals"],
        "same_point_rho": [row["rho"] for row in analysis["point_rows"]],
    }
    for candidate_id, manifest in manifests.items():
        write_once(HERE / f"candidate_{candidate_id}.json", manifest)
    write_once(HERE / f"ic1_rows_{RUN_ID}.json", output)
    print(json.dumps({"candidate_count": len(manifests),
                      "workload_count": len(workloads), "row_count": len(rows)},
                     sort_keys=True))


if __name__ == "__main__":
    main()
