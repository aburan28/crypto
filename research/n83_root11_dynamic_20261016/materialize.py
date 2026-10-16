#!/usr/bin/env python3
"""Bind the replayed dynamic one-root screen to exact IC1 identities."""

import argparse
import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
PANEL = REPO / "research/n83_factor_base_100_results_20261014/materialize.py"
SOURCE_COMMIT = "3890b32d4a5269fec0be6c65fe1fc2fce3cf8cff"
SOURCE_SHA256 = "d094345474605075566d646ddecef8f694cdd644b0f9190a6b1ef1551fd06888"
ARCHIVE_SHA256 = "ed24dcf04f751af6a7d832c21bb8e26295bd2cc4635eacfaec994751b3b82ff6"
ANALYSIS_SHA256 = "91aa67e82993dc0ebf001768b0d91be77d79ba610c713d685dc1f43805c4e0a3"
ARMS = (("k1271fixed_a", 1271, "root11_fixed"),
        ("k1271dynamic_a", 1271, "root11_dynamic"),
        ("k1271root22", 1271, "root22_dynamic"),
        ("k1024dynamic", 1024, "root11_dynamic"),
        ("k1536dynamic", 1536, "root11_dynamic"),
        ("k1271dynamic_b", 1271, "root11_dynamic"),
        ("k1271fixed_b", 1271, "root11_fixed"))

spec = importlib.util.spec_from_file_location("n83_panel_identity", PANEL)
panel = importlib.util.module_from_spec(spec)
spec.loader.exec_module(panel)
panel.SOURCE_COMMIT = SOURCE_COMMIT
panel.SOURCE_SHA256 = SOURCE_SHA256


def write_once(path, value):
    blob = json.dumps(value, sort_keys=True, indent=2) + "\n"
    if path.exists() and path.read_text() != blob:
        raise ValueError(f"sealed file differs: {path.name}")
    path.write_text(blob)


def candidate(k, layout, summary, dump_path, field, curve):
    root11 = layout.startswith("root11")
    dynamic = layout == "root11_dynamic"
    state_count = 83 * k * (k - 1) // 2 + 42 * k
    request = ((4 * state_count + 2) // 3 if layout != "root11_fixed"
               else 1 << 29)
    slots = max(16, 1 << (2 * request - 1).bit_length())
    if (summary["sparse_canonical_enabled"] is not True
            or summary["root_bloom_enabled"] is not False
            or summary["root_table_request_entries"] != request
            or summary["root_table_allocated_bytes"] != 5 * slots
            or summary["indexed5_load75_enabled"] is not
            (layout != "root11_fixed")
            or summary["indexed5_fixed30_enabled"] is not
            (layout == "root11_fixed")
            or summary["root_only11_enabled"] is not root11):
        raise ValueError(f"K{k} {layout}: executed exact-index policy differs")
    # The shared builder understands fixed one-root and dynamic root22.
    # Normalize only the two capacity flags to construct the same base and
    # field record, then replace every affected algorithm field below.
    template = dict(summary)
    if dynamic:
        template["indexed5_load75_enabled"] = False
        template["indexed5_fixed30_enabled"] = True
    record = panel.candidate(k, "root11" if root11 else "root22",
                             template, panel.read(dump_path), dump_path,
                             field, curve)
    record["implementation"]["algorithm_flags"]["sparse_canonical"] = True
    record["point_decomposition"]["canonical_root_selector"] = "sparse_linear_map_exact"
    record["point_decomposition"]["root_slot"].update({
        "request_entries": request, "allocated_bytes": 5 * slots})
    if dynamic:
        record["implementation"]["algorithm_flags"].update({
            "indexed5_load75_request": True,
            "indexed5_fixed30_request": False})
        record["point_decomposition"]["root_slot"]["request_rule"] = (
            "ceil_four_times_retained_states_over_three_then_round_twice_request_to_power_of_two")
        record["point_decomposition"]["variant"] += "_dynamic_load75"
    body = {key: value for key, value in record.items() if key != "candidate_id"}
    record["candidate_id"] = (
        f"IC1N83Ckb1a1fb{166*k}PDP4rootRCguidedLAgaussTDdirectISO0h"
        + panel.digest(body)[:12])
    return record


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("results", type=Path)
    parser.add_argument("--archive", type=Path, required=True)
    args = parser.parse_args()
    root = args.results.resolve()
    if (panel.file_sha(args.archive) != ARCHIVE_SHA256
            or panel.file_sha(HERE / "analysis_38876385086.json") != ANALYSIS_SHA256):
        raise ValueError("archived raw artifact or independent replay digest differs")
    analysis = panel.read(HERE / "analysis_38876385086.json")
    source = panel.read(root / "source.json")
    host = panel.read(root / "host.json")
    if (panel.read(root / "analysis.json") != analysis
            or analysis["status"] != "SCREEN_COMPLETE"
            or analysis["evidence_status"] != "PASS_INDEPENDENT_REPLAY"
            or host["git_sha"] != SOURCE_COMMIT
            or source["source_sha256"] != SOURCE_SHA256
            or source["hashes"]["examples/koblitz_orbit_dlp_fast.rs"] != SOURCE_SHA256
            or source["hashes"]["Cargo.lock"] != panel.LOCK_SHA256
            or source["arm_order"] != [arm for arm, _, _ in ARMS]):
        raise ValueError("executed source, arm order, or replay status differs")
    target_path = REPO / "research/n83_cold_factor_base_20261016/targets.jsonl"
    answer_path = REPO / "research/n83_cold_factor_base_20261016/answers.jsonl"
    point = json.loads(target_path.read_text())
    answer = json.loads(answer_path.read_text())
    if answer["fixture_index"] != 0 or answer["point"] != point:
        raise ValueError("primary public point and answer sidecar differ")
    workload = {
        "schema_version": 1, "curve_id": panel.CURVE_ID,
        "target_count": 1, "public_targets": [point],
        "input_law": {
            "domain": "n83-cold-factor-base-one-plus-four-20261016-v1",
            "fixture_index": 0,
            "scalar_rule": "1+(little_endian(SHA256(domain||u32_le(index))) mod (r-1))",
            "seed_sha256": answer["seed_sha256"],
            "generator_sha256": source["hashes"][
                "research/n83_cold_factor_base_20261016/generate.py"],
            "answer_sidecar_sha256": panel.file_sha(answer_path)},
        "target_sha256": panel.file_sha(target_path),
        "cold_or_warm": "one_public_point_per_fresh_process_empty_base_index_and_logs",
    }
    workload_id = panel.digest(workload)[:12]
    if workload_id != "eebf98ee1031":
        raise ValueError("primary one-point workload identity differs")
    field, curve = panel.curve_record()
    by_arm = {row["arm"]: row for row in analysis["arms"]}
    manifests, rows = {}, []
    for arm, k, layout in ARMS:
        item = by_arm[arm]
        if item["status"] != "PASS_INDEPENDENT_REPLAY" or item["resources"]["status"] != "PASS":
            raise ValueError(f"{arm}: full independent replay required")
        arm_root = root / "runs" / arm
        summary = panel.read(arm_root / "stdout.jsonl")
        record = candidate(k, layout, summary, arm_root / "base_dump.jsonl",
                           field, curve)
        if (item["base_hash"] != summary["base_hash"]
                or item["actual_factor_base_points"] != 166 * k):
            raise ValueError(f"{arm}: replayed base identity differs")
        candidate_id = record["candidate_id"]
        old = manifests.get(candidate_id)
        if old is not None and old != record:
            raise ValueError(f"{arm}: candidate hash collision")
        manifests[candidate_id] = record
        repetition = 2 if arm.endswith("_b") else 1
        rows.append({
            "arm": arm, "candidate_id": candidate_id,
            "workload_id": workload_id,
            "run_id": candidate_id + "W" + workload_id + "R" + str(repetition),
            "status": item["status"], "hardware_tier": analysis["hardware_tier"],
            "source_commit": SOURCE_COMMIT,
            "executed_binary_sha256": source["binary_sha256"],
            "rho_binary_sha256": source["rho_binary_sha256"],
            "host_receipt_sha256": panel.file_sha(root / "host.json"),
            "source_receipt_sha256": panel.file_sha(root / "source.json"),
            "correctness_receipt_sha256": ANALYSIS_SHA256,
            "raw_artifact_run_id": 38876385086,
            "raw_artifact_id": 12036957402,
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
            "raw_evidence": item,
        })
    output = {
        "schema": "n83.root11-dynamic-ic1/v1",
        "workload": {"workload_id": workload_id, **workload},
        "selection_status": analysis["status"],
        "evidence_status": analysis["evidence_status"],
        "selected_candidate_id": None,
        "fixed_aa_control": analysis["fixed_aa_control"],
        "dynamic_aa_control": analysis["dynamic_aa_control"],
        "same_point_rho": analysis["rho"],
        "rows": rows, "candidates": sorted(manifests),
        "source_commit": SOURCE_COMMIT,
        "archive_sha256": ARCHIVE_SHA256,
    }
    for candidate_id, record in manifests.items():
        write_once(HERE / f"candidate_{candidate_id}.json", record)
    write_once(HERE / "ic1_rows_38876385086.json", output)
    print(json.dumps({"candidate_count": len(manifests),
                      "row_count": len(rows),
                      "workload_id": workload_id}, sort_keys=True))


if __name__ == "__main__":
    main()
