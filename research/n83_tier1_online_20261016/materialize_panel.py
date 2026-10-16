#!/usr/bin/env python3
"""Bind the replayed five-point Tier 1 panel to exact one-target IC1 rows."""

import argparse
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
SOURCE_COMMIT = "4238a20bacba174b634cd9efdf88d048e61a77ff"
RUN_ID = 38921627515
SELECTOR_ARTIFACT_ID = 12058065194
SELECTOR_ZIP_SHA256 = "fa3f471254ed965ad3f89904dbf3f7e7dfd8f920893217053412da5211c9a030"
FULL_ARTIFACT_ID = 12059186751
FULL_ZIP_SHA256 = "73d79f6bd94671203284572aaeae81baf59048fd2c478be77977bee07b8b6e10"
ANALYSIS_SHA256 = "627b10ed0eb69fc2f9571f2e64c6291f78d5de3b5a6e4f2214a321015e7535a5"
LOCAL_REPLAY_SHA256 = "d4963ef08f16c47cd804c493136794c6005af190c3fff416d483750f1a48fc57"


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read(path):
    return json.loads(path.read_text())


def identity(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"),
                                     ensure_ascii=False).encode()).hexdigest()[:12]


def write_once(path, value):
    blob = json.dumps(value, sort_keys=True, indent=2) + "\n"
    if path.exists() and path.read_text() != blob:
        raise ValueError(f"sealed record differs: {path}")
    path.write_text(blob)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("results", type=Path)
    parser.add_argument("--archive", required=True, type=Path)
    parser.add_argument("--selector-archive", required=True, type=Path)
    parser.add_argument("--replay-log", required=True, type=Path)
    args = parser.parse_args()
    results = args.results.resolve()
    if (sha256(args.archive) != FULL_ZIP_SHA256
            or sha256(args.selector_archive) != SELECTOR_ZIP_SHA256
            or sha256(args.replay_log) != LOCAL_REPLAY_SHA256
            or sha256(results / "analysis.json") != ANALYSIS_SHA256
            or sha256(HERE / f"analysis_panel_{RUN_ID}.json") != ANALYSIS_SHA256):
        raise ValueError("raw ZIP, independent replay, or analysis identity differs")
    analysis = read(results / "analysis.json")
    host = read(results / "host.json")
    source = read(results / "source.json")
    selection = read(HERE / "selection.json")
    seal = read(HERE / "github_artifact_seal.json")
    candidate = read(HERE / f"candidate_{selection['candidate_id']}.json")
    if (analysis["status"] != "TIER1_ONLINE_SPEEDUP_PASS"
            or analysis["evidence_status"] != "PASS_INDEPENDENT_REPLAY"
            or analysis["execution_order_status"] != "PASS"
            or analysis["hardware_tier"] != "tier1_github_hosted_runner"
            or analysis["candidate_id"] != selection["candidate_id"]
            or analysis["pair_count"] != 5
            or analysis["admitted_pair_count"] != 5
            or len(analysis["point_rows"]) != 5
            or len(analysis["cross_point_base_identity"]) != 1
            or analysis["paired_log_ratio_interval"]["lower_95_log_ratio"] <= 0
            or host["git_sha"] != SOURCE_COMMIT
            or source["candidate_id"] != selection["candidate_id"]
            or host["instance_identity"]["runId"] != RUN_ID
            or host["cpu_affinity"] != seal["selector_cpus"]
            or seal["artifact_id"] != SELECTOR_ARTIFACT_ID
            or seal["artifact_digest"] != "sha256:" + SELECTOR_ZIP_SHA256
            or seal["head_sha"] != SOURCE_COMMIT
            or seal["selection_sha256"] != sha256(HERE / "selection.json")
            or seal["candidate_manifest_sha256"] != sha256(
                HERE / f"candidate_{selection['candidate_id']}.json")
            or candidate["candidate_id"] != selection["candidate_id"]
            or selection["screen_commit"] != SOURCE_COMMIT
            or selection["k"] != 3072
            or candidate["factor_base"]["actual_usable_points"] != 509952
            or candidate["factor_base"]["effective_columns"] != 3072
            or candidate["factor_base"]["source_base_hash_blake3"] !=
            analysis["cross_point_base_identity"][0][0]):
        raise ValueError("selected candidate, host, seal, or five-point gate differs")
    targets = [json.loads(line) for line in
               (HERE / "targets.jsonl").read_text().splitlines()]
    answers = [json.loads(line) for line in
               (HERE / "answers.jsonl").read_text().splitlines()]
    if len(targets) != 5 or len(answers) != 5:
        raise ValueError("five frozen target fixtures required")
    workloads, rows = [], []
    for index, (point, answer, measured) in enumerate(
            zip(targets, answers, analysis["point_rows"]), start=1):
        ic, rho = measured["ic"], measured["rho"]
        point_input = results / "inputs" / f"point_{index}.jsonl"
        if (answer["fixture_index"] != index - 1
                or answer["point"] != point
                or measured["point_index"] != index
                or measured["public_point"] != point
                or json.loads(point_input.read_text()) != point
                or source["point_input_sha256"][str(index)] != sha256(point_input)
                or not measured["pair_admitted"]
                or ic["status"] != "PASS_INDEPENDENT_REPLAY"
                or rho["status"] != "PASS_INDEPENDENT_REPLAY"
                or measured["calibration_aa"]["status"] != "PASS_10_PERCENT"
                or ic["resources"]["status"] != "PASS"
                or rho["resources"]["status"] != "PASS"
                or measured["online_ratio"] <= 1
                or ic["base_hash"] != analysis["cross_point_base_identity"][0][0]
                or ic["root_table_entries"] !=
                analysis["cross_point_base_identity"][0][1]
                or ic["actual_factor_base_points"] !=
                candidate["factor_base"]["actual_usable_points"]
                or ic["effective_columns"] !=
                candidate["factor_base"]["effective_columns"]):
            raise ValueError(f"point {index}: frozen input or paired replay differs")
        workload = {
            "schema_version": 1,
            "curve_id": candidate["curve"]["curve_id"],
            "target_count": 1,
            "public_targets": [point],
            "input_law": {
                "domain": "n83-tier1-online-five-fresh-20261016-v1",
                "fixture_index": index - 1,
                "scalar_rule": "1+(little_endian(SHA256(domain||u32_le(index))) mod (r-1))",
                "seed_sha256": answer["seed_sha256"],
                "generator_sha256": sha256(HERE / "generate.py"),
                "answer_sidecar_sha256": sha256(HERE / "answers.jsonl"),
            },
            "target_sha256": sha256(point_input),
            "cold_or_warm": "one_public_point_per_fresh_process_empty_base_index_and_logs",
        }
        workload_id = identity(workload)
        workloads.append({"workload_id": workload_id, **workload})
        rows.append({
            "candidate_id": selection["candidate_id"],
            "workload_id": workload_id,
            "run_id": selection["candidate_id"] + "W" + workload_id + "R1",
            "point_index": index,
            "public_point": point,
            "status": "PASS_INDEPENDENT_REPLAY",
            "hardware_tier": analysis["hardware_tier"],
            "source_commit": SOURCE_COMMIT,
            "executed_binary_sha256": source["binary_sha256"],
            "rho_binary_sha256": source["rho_binary_sha256"],
            "host_receipt_sha256": sha256(results / "host.json"),
            "source_receipt_sha256": sha256(results / "source.json"),
            "correctness_receipt_sha256": ANALYSIS_SHA256,
            "raw_artifact_run_id": RUN_ID,
            "raw_artifact_id": FULL_ARTIFACT_ID,
            "selector_artifact_id": SELECTOR_ARTIFACT_ID,
            "one_target_online_ms": ic["one_target_online_ms"],
            "rho_online_ms": rho["online_ms"],
            "online_speedup": measured["online_ratio"],
            "online_phases_ms": ic["online_phases_ms"],
            "cold_process_ms_supplementary": ic["cold_process_ms_exploratory"],
            "actual_usable_points": ic["actual_factor_base_points"],
            "effective_columns": ic["effective_columns"],
            "root_table_entries": ic["root_table_entries"],
            "peak_rss_bytes": ic["resources"]["peak_rss_bytes"],
            "calibration_aa": measured["calibration_aa"],
        })
    write_once(HERE / f"ic1_rows_{RUN_ID}.json", {
        "schema": "n83.tier1-online-ic1/v1",
        "candidate_id": selection["candidate_id"],
        "source_commit": SOURCE_COMMIT,
        "full_archive_sha256": FULL_ZIP_SHA256,
        "selector_archive_sha256": SELECTOR_ZIP_SHA256,
        "analysis_sha256": ANALYSIS_SHA256,
        "local_replay_log_sha256": LOCAL_REPLAY_SHA256,
        "selection_status": selection["screen_status"],
        "panel_status": analysis["status"],
        "hardware_tier": analysis["hardware_tier"],
        "paired_log_ratio_interval": analysis["paired_log_ratio_interval"],
        "median_online_ratio": analysis["median_online_ratio"],
        "workloads": workloads,
        "rows": rows,
    })
    print(json.dumps({"candidate_id": selection["candidate_id"],
                      "row_count": len(rows),
                      "workload_ids": [item["workload_id"] for item in workloads]},
                     sort_keys=True))


if __name__ == "__main__":
    main()
