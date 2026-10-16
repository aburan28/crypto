#!/usr/bin/env python3
"""Bind replayed K3300 Tier 1 receipts to five one-target IC1 rows."""

import argparse
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
STUDY = ROOT / "research/n83_tier1_refine_reserve_20261015"
PANEL = ROOT / "research/n83_factor_base_refine_results_20261014"
CURVE_ID = "EC1N83Ckb1a1haacb69ee9f37"
SOURCE_COMMIT = "d0e4aa13074c9fea910592c6e26b60b5e97fc439"
SOURCE_SHA256 = "5e9fd826d22816b429a82e658c517d9387916aabe7660f708e4c29901cce647e"
LOCK_SHA256 = "ac5cd7a3459be776180afe1eb2975c76b1a88f84bd1e99aa20475927ec260a90"
TARGETS_SHA256 = "2ccaac2baf367598e0d4e724b500ba66c5b4a711f0be8207565ce3acd6b1d0b4"
ANSWERS_SHA256 = "5b170d9bc4d34fdd21c032990ceee3a62af335c9af7fcaae2596335b8699dd8c"
SELECTION_SHA256 = "a9f7042ae415ae0d6aebf70fd82afce8baa0464db7288a9390c03487682b8c60"
PANEL_SELECTION_SHA256 = "307f07dfe9c2b222b99da44723eb28a7a5aaeb3f41e5ab348ea609c3b4104bef"
CANDIDATE_MANIFEST_SHA256 = "bf102583b9352a6a692a66551c763f274ab34ff91956ce1e465501e75ce39785"
CANDIDATE_ID = "IC1N83Ckb1a1fb547800PDP4rootRCguidedLAgaussTDdirectISO0h0f8ab8bfcff7"
RSS_CAP_BYTES = 14 * 1024**3


def file_sha(path):
    sha = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            sha.update(block)
    return sha.hexdigest()


def read(path):
    return json.loads(path.read_text())


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"),
                      ensure_ascii=False).encode("utf-8")


def digest(value):
    return hashlib.sha256(canonical(value)).hexdigest()


def write_once(path, value):
    if path.exists():
        if read(path) != value:
            raise ValueError(f"existing sealed record differs: {path}")
        return
    with path.open("x") as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write("\n")


def one_target_workload(point, index):
    record = {
        "schema_version": 1, "curve_id": CURVE_ID,
        "target_count": 1, "public_targets": [point],
        "input_law": {
            "domain": "n83-tier1-refine-reserve5-20261014-v1",
            "fixture_index": index,
            "scalar_rule": "1+(little_endian(SHA256(domain||u32_le(index))) mod (r-1))",
            "five_point_fixture_sha256": TARGETS_SHA256,
        },
        "cold_or_warm": "one_fresh_process_with_empty_index_and_logs; one_target_online_after_reusable_preparation",
    }
    return {"workload_id": digest(record)[:12], **record}


def source_checks(root):
    host, source, analysis = (read(root / name) for name in
                              ("host.json", "source.json", "analysis.json"))
    selection_path = STUDY / "selection.json"
    panel_path = PANEL / "selected_candidate.json"
    manifest_path = PANEL / "candidate_K3300_root11.json"
    target_path = ROOT / "research/n83_factor_base_refine_20261014/reserve_targets.jsonl"
    answer_path = ROOT / "research/n83_factor_base_refine_20261014/reserve_answers.jsonl"
    frozen = {
        selection_path: SELECTION_SHA256,
        panel_path: PANEL_SELECTION_SHA256,
        manifest_path: CANDIDATE_MANIFEST_SHA256,
        target_path: TARGETS_SHA256,
        answer_path: ANSWERS_SHA256,
    }
    if any(file_sha(path) != expected for path, expected in frozen.items()):
        raise ValueError("sealed selection, manifest or target fixture changed")
    selection = read(selection_path)
    panel = read(panel_path)
    manifest = read(manifest_path)
    hashes = source["hashes"]
    if (host["git_sha"] != SOURCE_COMMIT
            or source["schema"] != "n83.tier1-refine-reserve-source/v1"
            or source["selection_sha256"] != SELECTION_SHA256
            or source["selection"] != selection
            or source["runner_sha256"] != file_sha(STUDY / "run.py")
            or hashes["examples/koblitz_orbit_dlp_fast.rs"] != SOURCE_SHA256
            or hashes["Cargo.lock"] != LOCK_SHA256
            or hashes[str(target_path.relative_to(ROOT))] != TARGETS_SHA256
            or hashes[str(answer_path.relative_to(ROOT))] != ANSWERS_SHA256
            or panel["status"] != "PANEL_SELECTED"
            or panel["candidate_id"] != CANDIDATE_ID
            or manifest["candidate_id"] != CANDIDATE_ID
            or manifest["curve"]["curve_id"] != CURVE_ID
            or manifest["implementation"]["source_sha256"]["ic"] != SOURCE_SHA256
            or manifest["factor_base"]["source_base_hash_blake3"] !=
            panel["source_base_hash_blake3"]
            or selection["candidate_id"] != CANDIDATE_ID
            or selection["k"] != 3300 or selection["variant"] != "root11"
            or selection["rss_cap_bytes"] != RSS_CAP_BYTES
            or analysis["schema"] != "n83.tier1-refine-reserve-analysis/v1"
            or analysis["candidate_id"] != CANDIDATE_ID
            or analysis["target_count"] != 5):
        raise ValueError("executed source, selected candidate or analysis identity differs")
    points = [json.loads(line) for line in target_path.read_text().splitlines() if line]
    answers = [json.loads(line) for line in answer_path.read_text().splitlines() if line]
    if len(points) != 5 or len(answers) != 5:
        raise ValueError("fresh target fixture size changed")
    if any(answer["point"] != point or answer["fixture_index"] != index
           for index, (point, answer) in enumerate(zip(points, answers))):
        raise ValueError("target and answer fixture pairing changed")
    return host, source, analysis, points


def pair_row(root, source, analysis, point, index, workload,
             artifact_sha256, analysis_sha256):
    pair = analysis["pairs"][index]
    if pair["pair"] != index + 1:
        raise ValueError("pair order changed")
    pair_root = root / "pairs" / f"p{index + 1}"
    row = {
        "candidate_id": CANDIDATE_ID,
        "workload_id": workload["workload_id"],
        "run_id": CANDIDATE_ID + "W" + workload["workload_id"] + "R1",
        "pair_index": index + 1,
        "status": pair["status"],
        "hardware_tier": "tier1_cloud_instance",
        "source_commit": SOURCE_COMMIT,
        "source_sha256": SOURCE_SHA256,
        "ic_binary_sha256": source["ic_binary_sha256"],
        "rho_binary_sha256": source["rho_binary_sha256"],
        "selection_sha256": SELECTION_SHA256,
        "analysis_sha256": analysis_sha256,
        "artifact_sha256": artifact_sha256,
        "host_receipt_sha256": file_sha(root / "host.json"),
        "source_receipt_sha256": file_sha(root / "source.json"),
        "public_target": point,
        "execution_order": pair.get("order"),
        "rss_guard_bytes": RSS_CAP_BYTES,
        "ic_online_ms": None, "rho_online_ms": None,
        "rho_over_ic_online": None, "ic_online_phases_ms": None,
        "ic_setup_ms": None, "ic_peak_rss_bytes": None,
        "rho_peak_rss_bytes": None,
        "noise": pair.get("noise"),
    }
    for name in ("ic", "rho"):
        usage = pair_root / name / "usage.json"
        if usage.exists():
            receipt = read(usage)
            row[name + "_process_status"] = receipt["status"]
            row[name + "_process_wall_seconds"] = receipt["wall_seconds"]
            row[name + "_peak_rss_bytes"] = receipt["peak_rss_bytes"]
            row[name + "_stop_reason"] = receipt["stop_reason"]
    if pair["status"] == "PASS_TIER1_PAIR":
        target = read(pair_root / "ic" / "targets_result.jsonl")
        summary = read(pair_root / "ic" / "stdout.jsonl")
        rho = read(pair_root / "rho" / "receipt.json")
        if (target["target"] != point or rho["public_target_q"] != point
                or summary["base_hash"] !=
                read(PANEL / "selected_candidate.json")["source_base_hash_blake3"]
                or float(target["target_ms"]) != pair["ic"]["online_ms"]
                or float(rho["online_ms"]) != pair["rho"]["online_ms"]):
            raise ValueError(f"pair {index + 1}: raw online or base identity differs")
        row.update({
            "ic_online_ms": pair["ic"]["online_ms"],
            "rho_online_ms": pair["rho"]["online_ms"],
            "rho_over_ic_online": pair["ratio_rho_over_ic"],
            "ic_online_phases_ms": {key: target[key] for key in
                ("target_query_ms", "target_pdp_ms", "target_relation_check_ms",
                 "target_descent_ms", "target_recovery_check_ms")},
            "ic_setup_ms": pair["ic"]["setup_ms"],
            "ic_setup_scope": "producer_base_and_index_setup_total; rank_stage_recorded_separately",
            "ic_producer_timing_ms": summary["timing_ms"],
            "ic_rank_stage_ms": summary["timing_ms"]["rank_stage"],
            "ic_final_la_ms": summary["timing_ms"]["linear_algebra"],
            "ic_rank_attempts": summary["rank_attempts"],
            "ic_verified_rank_relations": summary["rank_relations"],
            "ic_rank_failures": summary["rank_failures"],
            "ic_rank_rows_without_gain": summary["rank_rows_without_gain"],
            "ic_rank_probes_mean": summary["rank_probes_mean"],
            "ic_base_scanned_x": summary["base_scanned_x"],
            "ic_actual_factor_base_points": pair["ic"]["actual_factor_base_points"],
            "ic_effective_columns": summary["orbit_columns"],
            "ic_final_rank": summary["rank"],
            "ic_root_table_entries": pair["ic"]["root_table_entries"],
            "ic_base_hash": pair["ic"]["base_hash"],
            "independent_ic_scalar_replay": "PASS",
            "independent_rho_scalar_replay": "PASS",
            "independent_four_point_relation_replay": "PASS",
        })
    return row


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
    root = args.artifact_results.resolve()
    if (file_sha(args.artifact_zip) != args.artifact_sha256
            or file_sha(root / "analysis.json") != args.analysis_sha256):
        raise ValueError("raw artifact or independently replayed analysis digest mismatch")
    _, source, analysis, points = source_checks(root)
    if len(analysis["pairs"]) != 5:
        raise ValueError("five paired result rows expected")
    workloads = []
    rows = []
    for index, point in enumerate(points):
        workload = one_target_workload(point, index)
        workloads.append(workload)
        rows.append(pair_row(root, source, analysis, point, index, workload,
                             args.artifact_sha256, args.analysis_sha256))
    write_once(HERE / "workloads.json", workloads)
    write_once(HERE / "measurements.json", rows)
    print(json.dumps({"status": analysis["status"],
                      "workload_ids": [row["workload_id"] for row in rows],
                      "pair_statuses": [row["status"] for row in rows]},
                     sort_keys=True))


if __name__ == "__main__":
    main()
