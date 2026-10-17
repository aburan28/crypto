#!/usr/bin/env python3
"""Materialize exact IC1 identities for the independently replayed N83 grid."""

import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
HOLDOUT = ROOT / "research/n83_e2e_holdout_20261016"
ARCHIVE = Path("/Volumes/SSD990-2/n83-archive-20261016/raw/n83-e2e-grid-38940108434.zip")
RESULTS = Path("/Volumes/SSD990-2/n83-archive-20261016/extracted/n83-e2e-grid-38940108434/results")
RUN_ID = 38940108434
ARTIFACT_ID = 12067378823
SOURCE_COMMIT = "551b9dfc632dce7a5ab8218b1b1dff6ae4b4d785"
ARCHIVE_SHA256 = "dde28f98f344f960c6d7d367e503583796dabfe379ca810d76c50dc24fdcbcc2"
ANALYSIS_SHA256 = "040d4890dcda06b1fd1f1c9b9e1c3db5ae19d1cad10d03293e787484dc960da6"
LOCAL_REPLAY_SHA256 = "a04c3d1ade70d0cccc63d1a650cc373b2cb92d0980b8a1d2ddd0e2dbf3abfa9d"


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    value = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(value)
    return value


runner = module("n83_e2e_grid_identity_runner", HERE / "run.py")
dynamic = module("n83_e2e_grid_identity_dynamic",
                 ROOT / "research/n83_root11_dynamic_20261016/materialize.py")
panel = dynamic.panel
panel.SOURCE_COMMIT = SOURCE_COMMIT
panel.SOURCE_SHA256 = runner.SOURCE_SHA256


def main():
    if (panel.file_sha(ARCHIVE) != ARCHIVE_SHA256
            or panel.file_sha(RESULTS / "analysis.json") != ANALYSIS_SHA256
            or panel.file_sha(HOLDOUT / "screen_analysis.json") != ANALYSIS_SHA256
            or panel.file_sha(RESULTS.parent / "local_replay.stdout.txt") !=
            LOCAL_REPLAY_SHA256):
        raise ValueError("exact provider artifact or independent replay differs")
    analysis = panel.read(RESULTS / "analysis.json")
    source = panel.read(RESULTS / "source.json")
    host = panel.read(RESULTS / "host.json")
    if (analysis["status"] != "SCREEN_COMPLETE"
            or analysis["evidence_status"] != "PASS_INDEPENDENT_REPLAY"
            or analysis["aa_control"]["status"] != "PASS_10_PERCENT"
            or host["git_sha"] != SOURCE_COMMIT
            or source["source_sha256"] != runner.SOURCE_SHA256
            or source["hashes"]["Cargo.lock"] != panel.LOCK_SHA256
            or source["arm_order"] != [arm for arm, _ in runner.ARMS]):
        raise ValueError("executed grid source, order or replay status differs")
    point = json.loads(runner.cold.TARGETS.read_text())
    answer = json.loads(runner.cold.ANSWERS.read_text())
    if answer["fixture_index"] != 0 or answer["point"] != point:
        raise ValueError("frozen exposed public point differs")
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
            "answer_sidecar_sha256": panel.file_sha(runner.cold.ANSWERS),
        },
        "target_sha256": panel.file_sha(runner.cold.TARGETS),
        "cold_or_warm": "one_public_point_per_fresh_process_empty_base_index_and_logs",
    }
    workload_id = panel.digest(workload)[:12]
    if workload_id != "eebf98ee1031":
        raise ValueError("one-point workload identity differs")
    field, curve = panel.curve_record()
    rows, manifests = [], {}
    for item in analysis["arms"]:
        arm, k = item["arm"], item["k"]
        if (item["status"] != "PASS_INDEPENDENT_REPLAY"
                or item["resources"]["status"] != "PASS"
                or item["actual_factor_base_points"] != 166 * k):
            raise ValueError(f"{arm}: complete replay required")
        folder = RESULTS / "runs" / arm
        summary = panel.read(folder / "stdout.jsonl")
        manifest = dynamic.candidate(k, "root11_dynamic", summary,
                                     folder / "base_dump.jsonl", field, curve)
        candidate_id = manifest["candidate_id"]
        if candidate_id in manifests and manifests[candidate_id] != manifest:
            raise ValueError(f"{arm}: candidate identity collision")
        manifests[candidate_id] = manifest
        rows.append({
            "arm": arm, "candidate_id": candidate_id,
            "workload_id": workload_id,
            "run_id": candidate_id + "W" + workload_id +
                      ("R2" if arm.endswith("_b") else "R1"),
            "status": item["status"],
            "hardware_tier": analysis["hardware_tier"],
            "source_commit": SOURCE_COMMIT,
            "executed_binary_sha256": source["binary_sha256"],
            "rho_binary_sha256": source["rho_binary_sha256"],
            "host_receipt_sha256": panel.file_sha(RESULTS / "host.json"),
            "source_receipt_sha256": panel.file_sha(RESULTS / "source.json"),
            "correctness_receipt_sha256": ANALYSIS_SHA256,
            "raw_artifact_run_id": RUN_ID,
            "raw_artifact_id": ARTIFACT_ID,
            "raw_artifact_sha256": ARCHIVE_SHA256,
            "cold_process_ms": item["cold_process_ms_exploratory"],
            "one_target_online_ms": item["one_target_online_ms"],
            "online_phases_ms": item["online_phases_ms"],
            "setup_ms": item["setup_ms"],
            "rank_stage_ms": item["rank_stage_ms"],
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
        "schema": "n83.e2e-factor-base-ic1/v1",
        "source_commit": SOURCE_COMMIT,
        "raw_artifact_sha256": ARCHIVE_SHA256,
        "analysis_sha256": ANALYSIS_SHA256,
        "workload": {"workload_id": workload_id, **workload},
        "candidates": sorted(manifests), "rows": rows,
        "selection_status": analysis["status"],
        "selected_candidate_id": None,
        "aa_control": analysis["aa_control"],
        "same_point_rho": analysis["rho"],
    }
    for candidate_id, manifest in manifests.items():
        dynamic.write_once(HERE / f"candidate_{candidate_id}.json", manifest)
    dynamic.write_once(HERE / f"ic1_rows_{RUN_ID}.json", output)
    print(json.dumps({"candidate_count": len(manifests),
                      "row_count": len(rows), "workload_id": workload_id},
                     sort_keys=True))


if __name__ == "__main__":
    main()
