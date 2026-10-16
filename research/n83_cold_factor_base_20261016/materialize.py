#!/usr/bin/env python3
"""Bind independently replayed cold-screen rows to exact IC1 identities."""

import argparse
import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
PANEL = REPO / "research/n83_factor_base_100_results_20261014/materialize.py"
SOURCE_COMMIT = "c5e87681d2a6e081743a0a9913965ba92227f008"
SOURCE_SHA256 = "f78382ada0ee4eabfa6820d805d0e51598dbed82154b43ea309bcd489305e2aa"
ARCHIVE_SHA256 = "57b4e0391ca86dc2af40d0b4a00ac642d41069e563b91c34a041022b5506b6fa"
ANALYSIS_SHA256 = "609b86d76e6f50054711a3ab452f55f319e8b17706a3d12eed76fafc6cf60e06"
ARMS = (("k2543root22_a", 2543, "root22"),
        ("k256root22", 256, "root22"),
        ("k512root22", 512, "root22"),
        ("k768root22", 768, "root22"),
        ("k1024root22", 1024, "root22"),
        ("k1536root22", 1536, "root22"),
        ("k2048root22", 2048, "root22"),
        ("k3450root11wide3", 3450, "root11"),
        ("k2543root22_b", 2543, "root22"))

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


def manifest(k, variant, summary, dump_path, field, curve):
    record = panel.candidate(k, variant, summary, panel.read(dump_path),
                             dump_path, field, curve)
    if summary["sparse_canonical_enabled"] is not True:
        raise ValueError("sparse canonicalization differs")
    record["implementation"]["algorithm_flags"]["sparse_canonical"] = True
    record["point_decomposition"]["canonical_root_selector"] = "sparse_linear_map_exact"
    if variant == "root11":
        if (summary["root_bloom_enabled"] is not True
                or summary["root_bloom_profile"] != "wide3"
                or summary["root_bloom_hashes"] != 3
                or summary["root_bloom_bits_per_request"] != 16
                or summary["root_bloom_bytes"] != 1 << 30):
            raise ValueError("root11 wide3 Bloom profile differs")
        record["implementation"]["algorithm_flags"].update(
            {"root_bloom": True, "root_bloom_profile": "wide3"})
        record["point_decomposition"]["root_bloom_filter"] = {
            "profile": "wide3", "hashes": 3, "bits_per_request": 16,
            "allocated_bytes": 1 << 30,
            "policy": "membership_prefilter_then_exact_root_lookup"}
    elif (summary["root_bloom_enabled"] is not False
          or summary["root_bloom_profile"] != "none"
          or summary["root_table_request_entries"] !=
          (4 * (83 * k * (k - 1) // 2 + 42 * k) + 2) // 3):
        raise ValueError("root22 exact index request differs")
    body = {key: value for key, value in record.items() if key != "candidate_id"}
    record["candidate_id"] = (
        f"IC1N83Ckb1a1fb{166*k}PDP4rootRCguidedLAgaussTDdirectISO0h"
        + panel.digest(body)[:12])
    return record


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("results", type=Path,
                        help="extracted artifact results/ with raw runs and base dumps")
    parser.add_argument("--archive", required=True, type=Path,
                        help="compressed copy of the complete original artifact")
    args = parser.parse_args()
    root = args.results.resolve()
    if (panel.file_sha(args.archive) != ARCHIVE_SHA256
            or panel.file_sha(HERE / "analysis_replay_38863274637.json")
            != ANALYSIS_SHA256):
        raise ValueError("archived artifact or independent replay digest differs")
    analysis = panel.read(HERE / "analysis_replay_38863274637.json")
    source = panel.read(root / "source.json")
    host = panel.read(root / "host.json")
    if (analysis["status"] != "BOUNDED_UNKNOWN"
            or analysis["evidence_status"] != "PASS_INDEPENDENT_REPLAY"
            or host["git_sha"] != SOURCE_COMMIT
            or source["hashes"]["examples/koblitz_orbit_dlp_fast.rs"] != SOURCE_SHA256
            or source["hashes"]["Cargo.lock"] != panel.LOCK_SHA256
            or source["hashes"]["research/n83_cold_factor_base_20261016/targets.jsonl"]
            != panel.file_sha(HERE / "targets.jsonl")
            or source["arm_order"] != [name for name, _, _ in ARMS]):
        raise ValueError("executed source, target, or replay admission differs")
    if panel.read(root / "analysis.json") != analysis:
        raise ValueError("raw artifact analysis differs from checked replay")
    point = json.loads((HERE / "targets.jsonl").read_text())
    answer = json.loads((HERE / "answers.jsonl").read_text())
    if answer["fixture_index"] != 0 or answer["point"] != point:
        raise ValueError("primary target and answer sidecar differ")
    workload = {
        "schema_version": 1,
        "curve_id": panel.CURVE_ID,
        "target_count": 1,
        "public_targets": [point],
        "input_law": {
            "domain": "n83-cold-factor-base-one-plus-four-20261016-v1",
            "fixture_index": 0,
            "scalar_rule": "1+(little_endian(SHA256(domain||u32_le(index))) mod (r-1))",
            "seed_sha256": answer["seed_sha256"],
            "generator_sha256": source["hashes"][
                "research/n83_cold_factor_base_20261016/generate.py"],
            "answer_sidecar_sha256": panel.file_sha(HERE / "answers.jsonl")},
        "target_sha256": panel.file_sha(HERE / "targets.jsonl"),
        "cold_or_warm": "one_public_point_per_fresh_process_empty_base_index_and_logs",
    }
    workload_id = panel.digest(workload)[:12]
    field, curve = panel.curve_record()
    by_arm = {row["arm"]: row for row in analysis["arms"]}
    manifests = {}
    rows = []
    for arm, k, variant in ARMS:
        row = by_arm[arm]
        candidate_id = None
        if row["status"] == "PASS_INDEPENDENT_REPLAY":
            run_dir = root / "runs" / arm
            summary = panel.read(run_dir / "stdout.jsonl")
            candidate = manifest(k, variant, summary,
                                 run_dir / "base_dump.jsonl", field, curve)
            candidate_id = candidate["candidate_id"]
            if (row["base_hash"] != summary["base_hash"]
                    or row["actual_factor_base_points"] != 166 * k):
                raise ValueError(f"{arm}: replayed base identity differs")
            old = manifests.get(candidate_id)
            if old is not None and old != candidate:
                raise ValueError(f"{arm}: candidate hash collision")
            manifests[candidate_id] = candidate
        run_number = 2 if arm == "k2543root22_b" else 1
        rows.append({
            "arm": arm, "candidate_id": candidate_id,
            "workload_id": workload_id,
            "run_id": (candidate_id + "W" + workload_id + "R" + str(run_number)
                       if candidate_id else None),
            "status": row["status"],
            "hardware_tier": analysis["hardware_tier"],
            "source_commit": SOURCE_COMMIT,
            "executed_binary_sha256": source["binary_sha256"],
            "rho_binary_sha256": source["rho_binary_sha256"],
            "host_receipt_sha256": panel.file_sha(root / "host.json"),
            "source_receipt_sha256": panel.file_sha(root / "source.json"),
            "correctness_receipt_sha256": ANALYSIS_SHA256,
            "raw_artifact_run_id": 38863274637,
            "raw_artifact_id": 12032697449,
            "cold_process_ms": row.get("cold_process_ms_exploratory"),
            "one_target_online_ms": row.get("one_target_online_ms"),
            "online_phases_ms": row.get("online_phases_ms"),
            "setup_ms": row.get("setup_ms"),
            "rank_stage_ms": row.get("rank_stage_ms"),
            "targets_total_ms": row.get("targets_total_ms"),
            "rank_attempts": row.get("rank_attempts"),
            "rank_relations": row.get("rank_relations"),
            "target_probes": row.get("target_probes"),
            "actual_usable_points": row.get("actual_factor_base_points"),
            "effective_columns": row.get("effective_columns"),
            "peak_rss_bytes": row.get("resources", {}).get("peak_rss_bytes"),
            "outer_wall_seconds": row.get("full_child_wall_seconds",
                                          row.get("usage", {}).get("wall_seconds")),
            "resource_status": row.get("resources", {}).get("status"),
            "raw_evidence": row,
        })
    output = {
        "schema": "n83.cold-factor-base-ic1/v1",
        "workload": {"workload_id": workload_id, **workload},
        "selection_status": analysis["status"],
        "evidence_status": analysis["evidence_status"],
        "selected_candidate_id": None,
        "aa_control": analysis["aa_control"],
        "same_point_rho": analysis["rho"],
        "rows": rows,
        "candidates": sorted(manifests),
        "source_commit": SOURCE_COMMIT,
        "archive_sha256": ARCHIVE_SHA256,
    }
    for candidate_id, record in manifests.items():
        path = HERE / f"candidate_{candidate_id}.json"
        write_once(path, record)
    write_once(HERE / "ic1_rows_38863274637.json", output)
    print(json.dumps({"workload_id": workload_id,
                      "candidate_count": len(manifests),
                      "row_count": len(rows)}, sort_keys=True))


if __name__ == "__main__":
    main()
