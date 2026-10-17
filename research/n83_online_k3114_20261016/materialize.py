#!/usr/bin/env python3
"""Bind the replayed K3114 N83 screen to exact IC1 candidates and rows."""

import argparse
import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RUN_ID = 38916776487
SOURCE_COMMIT = "9abf6afe8bc52ad4f136a19ff1603d60c344f30d"
SOURCE_SHA256 = "3ee6bc11399b382e72789b152043e5a67abc5b783c4e998ad8c534d60baebc03"
ARTIFACT_ID = 12054666847
ARCHIVE_SHA256 = "08f09b113c7256d829a606f2e4c630ccd77e09605060002e654a615ccd32b4bc"
ANALYSIS_SHA256 = "7a5f69618b000dac4fdf2a94626ce0639620868cf0960c00bcab6dbbda9113de"


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    value = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(value)
    return value


dynamic = module("n83_k3114_candidate_helper",
                 ROOT / "research/n83_root11_dynamic_20261016/materialize.py")
panel = dynamic.panel
panel.SOURCE_COMMIT = SOURCE_COMMIT
panel.SOURCE_SHA256 = SOURCE_SHA256
runner = module("n83_k3114_runner", HERE / "run.py")


def write_once(path, value):
    blob = json.dumps(value, sort_keys=True, indent=2) + "\n"
    if path.exists() and path.read_text() != blob:
        raise ValueError(f"sealed record differs: {path}")
    path.write_text(blob)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("results", type=Path)
    parser.add_argument("--archive", required=True, type=Path)
    args = parser.parse_args()
    root = args.results.resolve()
    if (ARTIFACT_ID is None or ARCHIVE_SHA256 is None
            or ANALYSIS_SHA256 is None):
        raise ValueError("exact artifact and analysis identities are not filled")
    versioned = HERE / f"analysis_{RUN_ID}.json"
    if (panel.file_sha(args.archive) != ARCHIVE_SHA256
            or panel.file_sha(root / "analysis.json") != ANALYSIS_SHA256
            or panel.file_sha(versioned) != ANALYSIS_SHA256):
        raise ValueError("raw ZIP or independent analysis hash differs")
    analysis = panel.read(root / "analysis.json")
    source = panel.read(root / "source.json")
    host = panel.read(root / "host.json")
    if (analysis["evidence_status"] != "PASS_INDEPENDENT_REPLAY"
            or analysis["execution_order_status"] != "PASS"
            or host["git_sha"] != SOURCE_COMMIT
            or source["source_sha256"] != SOURCE_SHA256
            or source["hashes"]["examples/koblitz_orbit_dlp_fast.rs"] != SOURCE_SHA256
            or source["runner_sha256"] != panel.file_sha(HERE / "run.py")
            or source["arm_order"] != [name for name, _, _, _ in runner.ARMS]):
        raise ValueError("executed source, order or independent replay differs")
    workload = panel.read(ROOT / "research/n83_root11_dynamic_20261016"
                          / "ic1_rows_38876385086.json")["workload"]
    workload_id = workload["workload_id"]
    if workload_id != "eebf98ee1031":
        raise ValueError("same exposed one-point workload identity differs")
    field, curve = panel.curve_record()
    by_arm = {row["arm"]: row for row in analysis["arms"]}
    manifests = {}
    by_k = {}
    rows = []
    for name, k, layout, _ in runner.ARMS:
        item = by_arm[name]
        if (item["status"] != "PASS_INDEPENDENT_REPLAY"
                or item["resources"]["status"] != "PASS"):
            rows.append({"arm": name, "candidate_id": None,
                         "workload_id": workload_id, "run_id": None,
                         "status": item["status"], "raw_evidence": item})
            continue
        folder = root / "runs" / name
        summary = panel.read(folder / "stdout.jsonl")
        record = dynamic.candidate(k, layout, summary,
                                   folder / "base_dump.jsonl", field, curve)
        if (item["base_hash"] != summary["base_hash"]
                or item["actual_factor_base_points"] != 166 * k):
            raise ValueError(f"{name}: replayed base identity differs")
        candidate_id = record["candidate_id"]
        if k in by_k and by_k[k] != candidate_id:
            raise ValueError(f"K{k}: A/A candidate identities differ")
        by_k[k] = candidate_id
        manifests[candidate_id] = record
        repetition = "ab".index(name[-1]) + 1
        rows.append({
            "arm": name, "candidate_id": candidate_id,
            "workload_id": workload_id,
            "run_id": candidate_id + "W" + workload_id + "R" + str(repetition),
            "status": item["status"], "hardware_tier": analysis["hardware_tier"],
            "source_commit": SOURCE_COMMIT,
            "executed_binary_sha256": source["binary_sha256"],
            "rho_binary_sha256": source["rho_binary_sha256"],
            "host_receipt_sha256": panel.file_sha(root / "host.json"),
            "source_receipt_sha256": panel.file_sha(root / "source.json"),
            "correctness_receipt_sha256": ANALYSIS_SHA256,
            "raw_artifact_run_id": RUN_ID,
            "raw_artifact_id": ARTIFACT_ID,
            "one_target_online_ms": item["one_target_online_ms"],
            "online_phases_ms": item["online_phases_ms"],
            "cold_process_ms_exploratory": item["cold_process_ms_exploratory"],
            "actual_usable_points": item["actual_factor_base_points"],
            "effective_columns": item["effective_columns"],
            "peak_rss_bytes": item["resources"]["peak_rss_bytes"],
        })
    winner = analysis["selected_pipeline"]
    if analysis["status"] == "SCREEN_SELECTED":
        selected_k = {"k3072root11": 3072, "k3114root11": 3114}.get(winner)
        if selected_k is None or selected_k not in by_k:
            raise ValueError("selected exact candidate is absent")
        selected_candidate_id = by_k[selected_k]
    else:
        if winner is not None:
            raise ValueError("unadmitted winner was recorded")
        selected_candidate_id = None
    for candidate_id, record in manifests.items():
        write_once(HERE / f"candidate_{candidate_id}.json", record)
    write_once(HERE / f"ic1_rows_{RUN_ID}.json", {
        "schema": "n83.online-k3114-ic1/v1",
        "workload": workload, "source_commit": SOURCE_COMMIT,
        "archive_sha256": ARCHIVE_SHA256,
        "analysis_sha256": ANALYSIS_SHA256,
        "candidates": sorted(manifests), "rows": rows,
        "selection_status": analysis["status"],
        "selected_candidate_id": selected_candidate_id,
        "online_aa_controls": analysis["online_aa_controls"],
        "same_point_rho": analysis["rho"],
    })
    print(json.dumps({"candidate_count": len(manifests),
                      "row_count": len(rows),
                      "selected_candidate_id": selected_candidate_id},
                     sort_keys=True))


if __name__ == "__main__":
    main()
