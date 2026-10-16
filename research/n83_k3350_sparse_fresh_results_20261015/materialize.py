#!/usr/bin/env python3
"""Seal replayed fresh N83 sparse K3300/K3350 candidates and runs."""

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
SOURCE_COMMIT = "01a93f468067bc5d7599c089a0e0cc53293b3b4b"
SOURCE_SHA256 = "c9824acebda0360416c3c507594bb07c17ca798c32cbf8f84bba95ca2c44d5b1"
LOCK_SHA256 = "ac5cd7a3459be776180afe1eb2975c76b1a88f84bd1e99aa20475927ec260a90"
POINT_SHA256 = "fc8be4c2c946bdbcb78a6845c8e2be3d1a51728c2ed4d812f6546ea5b000055e"
FIXTURE_SHA256 = "7121dcaf478061b87f54f28a0dc8e5af28e2637fd34a6837dd8f6434b6e29e8b"
EXPECTED_BASE_HASHES = {
    "k3350_root11": "dc2e0c5185b88b2b0c74ede0bbc7b889d756f2f21558e9b932850d97a1bdd07b",
    "k3300_root11": "a4ba33098d2776520b98a326585a2301f8ccd5bc005df7c1a70fbe8496c61d8b",
}
RSS_CAP_BYTES = 14 * 1024 ** 3
ARMS = ((3300, "root11"), (3350, "root11"))

spec = importlib.util.spec_from_file_location("n83_panel_materializer", PANEL)
panel = importlib.util.module_from_spec(spec)
spec.loader.exec_module(panel)
panel.SOURCE_COMMIT = SOURCE_COMMIT
panel.SOURCE_SHA256 = SOURCE_SHA256


def write_gzip_once(path, value):
    raw = (json.dumps(value, indent=2, sort_keys=True,
                      ensure_ascii=False) + "\n").encode("utf-8")
    buffer = io.BytesIO()
    with gzip.GzipFile(filename="", mode="wb", fileobj=buffer,
                       compresslevel=9, mtime=0) as stream:
        stream.write(raw)
    compressed = buffer.getvalue()
    if path.exists():
        if path.read_bytes() != compressed:
            raise ValueError(f"existing sealed record differs: {path}")
    else:
        path.write_bytes(compressed)
    return {"sha256": hashlib.sha256(compressed).hexdigest(),
            "json_sha256": hashlib.sha256(raw).hexdigest(),
            "compressed_bytes": len(compressed), "json_bytes": len(raw)}


def git_file(path):
    return subprocess.check_output(["git", "show", SOURCE_COMMIT + ":" + path], cwd=REPO)


def frozen_points():
    generator_path = "research/n83_k3350_sparse_fresh_20261015/generate.py"
    exclusion_path = "research/n83_k3350_sparse_fresh_20261015/excluded_bloom_targets.jsonl"
    for relative in (generator_path, exclusion_path):
        if (REPO / relative).read_bytes() != git_file(relative):
            raise ValueError(f"executed fixture source differs: {relative}")
    point_bytes = (REPO / "research/n83_k3350_sparse_fresh_20261015/targets.jsonl").read_bytes()
    answer_bytes = (REPO / "research/n83_k3350_sparse_fresh_20261015/answers.jsonl").read_bytes()
    if (hashlib.sha256(point_bytes).hexdigest() != POINT_SHA256
            or hashlib.sha256(answer_bytes).hexdigest() != FIXTURE_SHA256):
        raise ValueError("committed fresh points or answer sidecar changed")
    points = [json.loads(line) for line in point_bytes.splitlines() if line]
    answers = [json.loads(line) for line in answer_bytes.splitlines() if line]
    if (len(points) != 20000 or len(answers) != 20000
            or len({tuple(point) for point in points}) != 20000
            or any(answer["fixture_index"] != index or answer["point"] != point
                   for index, (point, answer) in enumerate(zip(points, answers)))):
        raise ValueError("fresh point count, uniqueness or answer pairing differs")
    return points


def sparse_candidate(k, variant, summary, dump_path, field, curve):
    record = panel.candidate(k, variant, summary, panel.read(dump_path),
                             dump_path, field, curve)
    arm = f"k{k}_{variant}"
    if (summary["sparse_canonical_enabled"] is not True
            or summary["root_bloom_enabled"] is not False
            or summary["root_bloom_bytes"] != 0
            or summary["base_hash"] != EXPECTED_BASE_HASHES[arm]):
        raise ValueError(f"{arm}: sparse selector or base differs")
    record["implementation"]["algorithm_flags"]["sparse_canonical"] = True
    record["point_decomposition"]["canonical_root_selector"] = "sparse_linear_map_exact"
    body = {key: value for key, value in record.items() if key != "candidate_id"}
    record["candidate_id"] = (
        f"IC1N83Ckb1a1fb{166*k}PDP4rootRCguidedLAgaussTDdirectISO0h"
        + panel.digest(body)[:12])
    return record


def seal(root, artifact_zip, artifact_sha256, analysis_sha256):
    if (panel.file_sha(artifact_zip) != artifact_sha256
            or panel.file_sha(root / "analysis.json") != analysis_sha256):
        raise ValueError("raw artifact or independent analysis digest mismatch")
    analysis = panel.read(root / "analysis.json")
    host = panel.read(root / "host.json")
    source = panel.read(root / "source.json")
    hashes = source["hashes"]
    if (host["git_sha"] != SOURCE_COMMIT
            or source["schema"] != "n83.k3350-sparse-fresh-source/v1"
            or hashes["examples/koblitz_orbit_dlp_fast.rs"] != SOURCE_SHA256
            or hashes["Cargo.lock"] != LOCK_SHA256
            or hashes["research/n83_k3350_sparse_fresh_20261015/targets.jsonl"] != POINT_SHA256
            or hashes["research/n83_k3350_sparse_fresh_20261015/answers.jsonl"] != FIXTURE_SHA256
            or source["panel_runner_sha256"] != panel.file_sha(
                REPO / "research/n83_k3350_sparse_fresh_20261015/run.py")
            or analysis["schema"] != "n83.k3350-sparse-fresh-analysis/v1"
            or analysis["hardware_tier"] != "cloud_instance_exploratory"
            or analysis["workload_role"] != "secondary_fresh_20000_point_sparse_factor_base_confirmation"
            or analysis["target_count"] != 20000
            or len(analysis["arms"]) != len(ARMS)):
        raise ValueError("source, input or allocation identity mismatch")
    points = frozen_points()
    workload = {
        "schema_version": 1, "curve_id": panel.CURVE_ID,
        "target_count": 20000, "public_targets": points,
        "input_law": {"domain": "n83-k3350-sparse-fresh-20000-20261015-v1",
                      "scalar_rule": "1+(little_endian(SHA256(domain||u32_le(index))) mod (r-1))",
                      "answer_sidecar_sha256": FIXTURE_SHA256},
        "target_sha256": POINT_SHA256,
        "cold_or_warm": "one_fresh_process_per_candidate_with_empty_index_and_logs; 20000_target_fresh_confirmation_after_shared_preparation",
    }
    workload_id = panel.digest(workload)[:12]
    workload_file = write_gzip_once(HERE / "workload.json.gz",
                                    {"workload_id": workload_id, **workload})
    field, curve = panel.curve_record()
    verified_by_arm = {row["arm"]: row for row in analysis["arms"]}
    if set(verified_by_arm) != {f"k{k}_{variant}" for k, variant in ARMS}:
        raise ValueError("analysis arms differ from frozen protocol")
    rows = []
    for k, variant in ARMS:
        arm = f"k{k}_{variant}"
        verified = verified_by_arm[arm]
        row = {"arm": arm, "candidate_id": None, "workload_id": workload_id,
               "run_id": None, "status": verified["status"],
               "hardware_tier": "cloud_instance_exploratory",
               "workload_role": "secondary_fresh_20000_point_sparse_factor_base_confirmation",
               "selection_status": analysis["status"],
               "controlled_speedup": None, "source_commit": SOURCE_COMMIT,
               "executed_binary_sha256": source["binary_sha256"],
               "host_receipt_sha256": panel.file_sha(root / "host.json"),
               "source_receipt_sha256": panel.file_sha(root / "source.json"),
               "correctness_receipt_sha256": analysis_sha256}
        arm_root = root / "runs" / arm
        usage_path = arm_root / "usage.json"
        if usage_path.exists():
            usage = panel.read(usage_path)
            row.update({"exit_status": usage["exit_status"],
                        "full_process_wall_seconds": usage["wall_seconds"],
                        "peak_rss_bytes": usage["peak_rss_bytes"],
                        "sampled_peak_rss_bytes": usage["sampled_peak_rss_bytes"],
                        "rss_guard_bytes": RSS_CAP_BYTES,
                        "stop_reason": usage["stop_reason"],
                        "major_faults": usage["major_faults"],
                        "involuntary_context_switches": usage["involuntary_context_switches"],
                        "resource_envelope_status": (
                            "OBSERVED_PEAK_ABOVE_14_GIB_GUARD"
                            if usage["peak_rss_bytes"] > RSS_CAP_BYTES else
                            "OBSERVED_PEAK_WITHIN_14_GIB_GUARD")})
        if verified["status"] == "PASS_INDEPENDENT_REPLAY":
            summary = panel.read(arm_root / "stdout.jsonl")
            targets = [json.loads(line) for line in
                       (arm_root / "targets_result.jsonl").read_text().splitlines() if line]
            if len(targets) != 20000:
                raise ValueError(f"{arm}: expected 20000 target rows")
            dump_path = arm_root / "base_dump.jsonl"
            manifest = sparse_candidate(k, variant, summary, dump_path, field, curve)
            candidate_id = manifest["candidate_id"]
            panel.write_once(HERE / f"candidate_K{k}_{variant}.json", manifest)
            row.update({
                "candidate_id": candidate_id,
                "run_id": candidate_id + "W" + workload_id + "R1",
                "actual_factor_base_points": verified["actual_factor_base_points"],
                "effective_columns": k, "final_rank": summary["rank"],
                "root_table_allocated_bytes": verified["root_table_allocated_bytes"],
                "root_state_allocated_bytes": verified["root_state_allocated_bytes"],
                "regular_states": summary["regular_states"],
                "root_table_entries": verified["root_table_entries"],
                "target_probes": [target["probes"] for target in targets],
                "ic_online_ms_per_target_exploratory":
                    [float(target["target_ms"]) for target in targets],
                "ic_online_mean_ms_exploratory": verified["mean_online_ms_exploratory"],
                "ic_online_phases_ms_per_target": [
                    {key: target[key] for key in
                     ("target_query_ms", "target_pdp_ms", "target_relation_check_ms",
                      "target_descent_ms", "target_recovery_check_ms")}
                    for target in targets],
                "setup_ms": verified["setup_ms"],
                "rank_ms": summary["timing_ms"]["rank_stage"],
                "final_la_ms": summary["timing_ms"]["linear_algebra"],
            })
        rows.append(row)
    measurements_file = write_gzip_once(HERE / "measurements.json.gz", rows)
    eligible = [row for row in rows if row["status"] == "PASS_INDEPENDENT_REPLAY"
                and row.get("peak_rss_bytes", RSS_CAP_BYTES + 1) <= RSS_CAP_BYTES]
    eligible_by_arm = {row["arm"]: row for row in eligible}
    comparisons = analysis["paired_comparisons"]
    expected_comparisons = {"k3350_root11"}
    baseline = eligible_by_arm.get("k3300_root11")
    if baseline is not None:
        if set(comparisons) != expected_comparisons:
            raise ValueError("paired contrast set changed")
    elif comparisons:
        raise ValueError("paired contrasts require the independently replayed baseline")
    else:
        comparisons = {arm: None for arm in sorted(expected_comparisons)}
    confirmed = [eligible_by_arm[arm] for arm, result in comparisons.items()
                 if (result is not None and result["confirmed_faster_than_k3300"]
                     and arm in eligible_by_arm)]
    selected = (min(confirmed, key=lambda row: row["ic_online_mean_ms_exploratory"])
                if confirmed else baseline)
    if analysis["selected_arm_confirmed"] != (selected["arm"] if selected else None):
        raise ValueError("frozen fresh factor-base selection changed")
    complete = (len(eligible) == len(ARMS)
                and analysis["status"] == "COMPLETE"
                and analysis["baseline_base_hash_matches"] is True)
    for row in eligible:
        if row["arm"] in EXPECTED_BASE_HASHES and row["arm"] in verified_by_arm:
            if verified_by_arm[row["arm"]]["base_hash"] != EXPECTED_BASE_HASHES[row["arm"]]:
                raise ValueError(f"frozen {row['arm']} base hash changed")
    selection = {
        "schema": "n83.k3350-sparse-fresh-selection/v1",
        "status": "FRESH_FACTOR_BASE_SELECTED" if complete else "INCOMPLETE",
        "selected_arm": selected["arm"] if selected else None,
        "candidate_id": selected["candidate_id"] if selected else None,
        "candidate_manifest_sha256": None, "source_base_hash_blake3": None,
        "mean_online_ms_exploratory": (selected["ic_online_mean_ms_exploratory"]
                                        if selected else None),
        "paired_comparisons": comparisons,
        "baseline_base_hash_matches": analysis["baseline_base_hash_matches"],
        "analysis_sha256": analysis_sha256,
        "artifact_sha256": artifact_sha256,
        "source_commit": SOURCE_COMMIT,
        "source_sha256": SOURCE_SHA256,
        "workload_id": workload_id,
    }
    if selected is not None:
        manifest_path = HERE / f"candidate_{selected['arm'].replace('k', 'K', 1)}.json"
        manifest = panel.read(manifest_path)
        selection["candidate_manifest_sha256"] = panel.file_sha(manifest_path)
        selection["source_base_hash_blake3"] = manifest["factor_base"]["source_base_hash_blake3"]
    panel.write_once(HERE / "selection.json", selection)
    panel.write_once(HERE / "records_index.json", {
        "schema": "n83.k3350-sparse-fresh-compressed-records/v1",
        "workload_id": workload_id,
        "measurement_rows": len(rows),
        "workload": workload_file,
        "measurements": measurements_file,
        "analysis_sha256": analysis_sha256,
        "artifact_sha256": artifact_sha256,
    })
    print(json.dumps({"workload_id": workload_id,
                      "candidate_ids": [row["candidate_id"] for row in rows],
                      "statuses": [row["status"] for row in rows],
                      "selected_arm": selection["selected_arm"],
                      "selection_status": selection["status"]}, sort_keys=True))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--artifact-zip", required=True, type=Path)
    parser.add_argument("--artifact-sha256", required=True)
    parser.add_argument("--analysis-sha256", required=True)
    parser.add_argument("artifact_results", type=Path)
    args = parser.parse_args()
    for name, value in (("artifact", args.artifact_sha256),
                        ("analysis", args.analysis_sha256)):
        if len(value) != 64 or any(ch not in "0123456789abcdef" for ch in value):
            parser.error(f"{name} SHA-256 must be 64 lowercase hex digits")
    seal(args.artifact_results.resolve(), args.artifact_zip.resolve(),
         args.artifact_sha256, args.analysis_sha256)


if __name__ == "__main__":
    main()
