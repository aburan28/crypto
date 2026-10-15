#!/usr/bin/env python3
"""Seal the independently replayed N83 root11 Bloom comparison."""

import argparse
import gzip
import hashlib
import importlib.util
import io
import json
import subprocess
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
PANEL = REPO / "research/n83_factor_base_100_results_20261014/materialize.py"
SOURCE_COMMIT = "a9b44ee1454beb2bbe8614b681714df07fd11f3c"
SOURCE_SHA = "814afa03ae1c386e62e518be3c7a55fb108b09c41b05bef64ca1371787466c39"
LOCK_SHA = "ac5cd7a3459be776180afe1eb2975c76b1a88f84bd1e99aa20475927ec260a90"
POINT_SHA = "35e34151dcc1ad6048939560adf28fc7c1b1017833fe005a656c3f4e2bc07de8"
BASE_HASH = "a4ba33098d2776520b98a326585a2301f8ccd5bc005df7c1a70fbe8496c61d8b"
WORKLOAD_PATH = REPO / "research/n83_factor_base_refine_results_20261014/workload.json.gz"
WORKLOAD_SHA = "cda8bdbbc908d8fbfd5ccef394d9afee0ccfbdeef17f54705828648487501ac5"
ARTIFACT_SHA = "098c1ea1383c3ddaf78aaf1ac28bc32c80cb2caa5c69ccfb155bb060f5d2bda2"
ANALYSIS_SHA = "2976f932e8eb2fd2d285b1f40f065c43dd6573fb6879845bd915f9b7ee675a6a"
RSS_CAP = 14 * 1024**3
ARMS = ("bloom_off", "bloom_on")

spec = importlib.util.spec_from_file_location("n83_bloom_panel_materializer", PANEL)
panel = importlib.util.module_from_spec(spec)
spec.loader.exec_module(panel)
panel.SOURCE_COMMIT = SOURCE_COMMIT
panel.SOURCE_SHA256 = SOURCE_SHA


def git_file(path):
    return subprocess.check_output(["git", "show", SOURCE_COMMIT + ":" + path], cwd=REPO)


def write_gzip_once(path, value):
    raw = (json.dumps(value, indent=2, sort_keys=True, ensure_ascii=False) + "\n").encode()
    buffer = io.BytesIO()
    with gzip.GzipFile(filename="", mode="wb", fileobj=buffer, compresslevel=9,
                       mtime=0) as stream:
        stream.write(raw)
    blob = buffer.getvalue()
    if path.exists() and path.read_bytes() != blob:
        raise ValueError(f"sealed record differs: {path}")
    if not path.exists():
        path.write_bytes(blob)
    return {"sha256": hashlib.sha256(blob).hexdigest(),
            "json_sha256": hashlib.sha256(raw).hexdigest(),
            "compressed_bytes": len(blob), "json_bytes": len(raw)}


def candidate(summary, dump_path, bloom, field, curve):
    record = panel.candidate(3300, "root11", summary, panel.read(dump_path),
                             dump_path, field, curve)
    if (summary["base_hash"] != BASE_HASH
            or summary["root_table_entries"] != 903852075
            or summary["root_bloom_enabled"] is not bloom
            or summary["root_bloom_bytes"] != (1 << 29 if bloom else 0)):
        raise ValueError("Bloom candidate base, index or prefilter differs")
    record["implementation"]["algorithm_flags"]["root_bloom"] = bloom
    record["point_decomposition"]["root_bloom_prefilter"] = (
        {"bits": 1 << 32, "hashes": 3,
         "rule": "negative_skips_exact_index_lookup; positive_uses_exact_lookup"}
        if bloom else "none")
    body = {key: value for key, value in record.items() if key != "candidate_id"}
    record["candidate_id"] = (
        "IC1N83Ckb1a1fb547800PDP4rootRCguidedLAgaussTDdirectISO0h"
        + panel.digest(body)[:12])
    return record


def seal(root, artifact_zip):
    if (panel.file_sha(artifact_zip) != ARTIFACT_SHA
            or panel.file_sha(root / "analysis.json") != ANALYSIS_SHA
            or panel.file_sha(WORKLOAD_PATH) != WORKLOAD_SHA
            or hashlib.sha256(git_file("research/n83_factor_base_refine_20261014/targets.jsonl")).hexdigest() != POINT_SHA):
        raise ValueError("raw, analysis, workload or target digest differs")
    analysis = panel.read(root / "analysis.json")
    host = panel.read(root / "host.json")
    source = panel.read(root / "source.json")
    if (host["git_sha"] != SOURCE_COMMIT
            or source["schema"] != "n83.root11-bloom-source/v1"
            or source["hashes"]["examples/koblitz_orbit_dlp_fast.rs"] != SOURCE_SHA
            or source["hashes"]["Cargo.lock"] != LOCK_SHA
            or source["hashes"]["research/n83_factor_base_refine_20261014/targets.jsonl"] != POINT_SHA
            or source["screen_runner_sha256"] != panel.file_sha(REPO / "research/n83_root11_bloom_20261015/run.py")
            or analysis["schema"] != "n83.root11-bloom-analysis/v1"
            or analysis["target_count"] != 1000
            or analysis["hardware_tier"] != "cloud_instance_exploratory"
            or set(row["arm"] for row in analysis["arms"]) != set(ARMS)):
        raise ValueError("executed source, inputs or arm set differs")
    with gzip.open(WORKLOAD_PATH) as stream:
        workload = json.load(stream)
    if workload["workload_id"] != "0d7b05d9fca6":
        raise ValueError("shared 1000-point workload identity differs")
    field, curve = panel.curve_record()
    verified = {row["arm"]: row for row in analysis["arms"]}
    rows = []
    for arm in ARMS:
        result = verified[arm]
        folder = root / "runs" / arm
        usage = panel.read(folder / "usage.json")
        row = {"arm": arm, "candidate_id": None, "run_id": None,
               "workload_id": workload["workload_id"], "status": result["status"],
               "source_commit": SOURCE_COMMIT,
               "executed_binary_sha256": source["binary_sha256"],
               "hardware_tier": "cloud_instance_exploratory",
               "workload_role": "secondary_reused_1000_point_layout_screen",
               "controlled_speedup": None,
               "correctness_receipt_sha256": ANALYSIS_SHA,
               "host_receipt_sha256": panel.file_sha(root / "host.json"),
               "exit_status": usage["exit_status"],
               "full_process_wall_seconds": usage["wall_seconds"],
               "peak_rss_bytes": usage["peak_rss_bytes"],
               "rss_guard_bytes": RSS_CAP,
               "stop_reason": usage["stop_reason"]}
        if result["status"] == "PASS_INDEPENDENT_REPLAY":
            summary = panel.read(folder / "stdout.jsonl")
            dump_path = folder / "base_dump.jsonl"
            manifest = candidate(summary, dump_path, arm == "bloom_on", field, curve)
            panel.write_once(HERE / f"candidate_{arm}.json", manifest)
            targets = [json.loads(line) for line in
                       (folder / "targets_result.jsonl").read_text().splitlines() if line]
            if len(targets) != 1000:
                raise ValueError(f"{arm}: target rows incomplete")
            row.update({"candidate_id": manifest["candidate_id"],
                        "run_id": manifest["candidate_id"] + "W" + workload["workload_id"] + "R1",
                        "actual_factor_base_points": result["actual_factor_base_points"],
                        "effective_columns": result["effective_columns"],
                        "root_table_entries": result["root_table_entries"],
                        "root_table_allocated_bytes": result["root_table_allocated_bytes"],
                        "root_state_allocated_bytes": result["root_state_allocated_bytes"],
                        "root_bloom_bytes": result["root_bloom_bytes"],
                        "setup_ms": result["setup_ms"],
                        "ic_online_mean_ms_exploratory": result["mean_online_ms_exploratory"],
                        "ic_online_ms_per_target_exploratory": [float(t["target_ms"]) for t in targets],
                        "ic_online_phases_ms_per_target": [
                            {key: target[key] for key in
                             ("target_query_ms", "target_pdp_ms", "target_relation_check_ms",
                              "target_descent_ms", "target_recovery_check_ms")}
                            for target in targets],
                        "target_probes": [target["probes"] for target in targets]})
        rows.append(row)
    measurement_index = write_gzip_once(HERE / "measurements.json.gz", rows)
    difference = analysis["paired_difference"]
    complete = (analysis["status"] == "COMPLETE"
                and analysis["same_base_and_key_count"] is True
                and all(row["status"] == "PASS_INDEPENDENT_REPLAY"
                        and row["peak_rss_bytes"] <= RSS_CAP for row in rows))
    if (difference["confirmed_faster"] is not True
            or difference["bootstrap_95pct_ms"][1] >= 0
            or analysis["selected_arm_exploratory"] != "bloom_on"):
        raise ValueError("frozen Bloom selection statistic differs")
    selection = {"schema": "n83.root11-bloom-selection/v1",
                 "status": "BLOOM_CONFIRMED_IN_REUSED_PANEL" if complete else "INCOMPLETE",
                 "selected_arm": "bloom_on" if complete else None,
                 "candidate_id": rows[1]["candidate_id"] if complete else None,
                 "paired_difference": difference,
                 "artifact_sha256": ARTIFACT_SHA, "analysis_sha256": ANALYSIS_SHA,
                 "source_commit": SOURCE_COMMIT,
                 "workload_id": workload["workload_id"]}
    panel.write_once(HERE / "selection.json", selection)
    panel.write_once(HERE / "records_index.json", {
        "schema": "n83.root11-bloom-compressed-records/v1",
        "workload_path": str(WORKLOAD_PATH.relative_to(REPO)),
        "workload_sha256": WORKLOAD_SHA,
        "measurements": measurement_index,
        "artifact_sha256": ARTIFACT_SHA,
        "analysis_sha256": ANALYSIS_SHA})
    print(json.dumps(selection, sort_keys=True))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--artifact-zip", type=Path, required=True)
    parser.add_argument("artifact_results", type=Path)
    args = parser.parse_args()
    seal(args.artifact_results.resolve(), args.artifact_zip.resolve())


if __name__ == "__main__":
    main()
