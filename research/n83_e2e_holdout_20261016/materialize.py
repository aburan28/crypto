#!/usr/bin/env python3
"""Bind a replayed N83 five-point holdout to exact IC1 candidate identities."""

import argparse
import hashlib
import importlib.util
import json
import statistics
import subprocess
import sys
import zipfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    value = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(value)
    return value


runner = module("n83_e2e_holdout_identity_run", HERE / "run.py")
dynamic = module("n83_e2e_holdout_identity_dynamic",
                 ROOT / "research/n83_root11_dynamic_20261016/materialize.py")
panel = dynamic.panel
panel.SOURCE_COMMIT = runner.SCREEN_HEAD_SHA
panel.SOURCE_SHA256 = runner.screen.SOURCE_SHA256


def file_sha(path):
    value = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def workload(index, point, answer, source, input_path):
    record = {
        "schema_version": 1, "curve_id": panel.CURVE_ID,
        "target_count": 1, "public_targets": [point],
        "input_law": {
            "domain": "n83-e2e-cold-holdout-five-20261016-v1",
            "fixture_index": index,
            "scalar_rule": "1+(little_endian(SHA256(domain||u32_le(index))) mod (r-1))",
            "seed_sha256": answer["seed_sha256"],
            "generator_sha256": file_sha(HERE / "generate.py"),
            "answer_sidecar_sha256": runner.ANSWER_SHA256,
        },
        "target_sha256": file_sha(input_path),
        "cold_or_warm": "one_public_point_per_fresh_process_empty_base_index_and_logs",
    }
    return panel.digest(record)[:12], record


def archive_analysis(archive, provider_digest):
    digest = provider_digest.removeprefix("sha256:")
    if len(digest) != 64 or any(c not in "0123456789abcdef" for c in digest):
        raise ValueError("provider digest must be lowercase SHA-256")
    if file_sha(archive) != digest:
        raise ValueError("holdout artifact differs from provider digest")
    with zipfile.ZipFile(archive) as saved:
        if saved.testzip() is not None:
            raise ValueError("holdout artifact ZIP integrity failed")
        for name in saved.namelist():
            if Path(name).is_absolute() or ".." in Path(name).parts:
                raise ValueError("unsafe artifact path")
        return saved.read("results/analysis.json"), digest


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("results", type=Path, help="exact extracted artifact results")
    parser.add_argument("--artifact-zip", type=Path, required=True)
    parser.add_argument("--provider-sha256", required=True)
    parser.add_argument("--run-id", type=int, required=True)
    parser.add_argument("--head-sha", required=True)
    parser.add_argument("--reuse-local-replay", action="store_true",
                        help="use a completed verifier log bound to this exact analysis")
    args = parser.parse_args()
    results = args.results.resolve()
    cloud_analysis, archive_sha = archive_analysis(
        args.artifact_zip.resolve(), args.provider_sha256)
    if (results / "analysis.json").read_bytes() != cloud_analysis:
        raise ValueError("extracted holdout analysis differs from exact ZIP")
    if args.reuse_local_replay:
        status = (results.parent / "local_replay.exit_status.txt").read_text().strip()
        stdout = (results.parent / "local_replay.stdout.txt").read_text()
        stderr = (results.parent / "local_replay.stderr.txt").read_text()
        if status != "0" or stderr or json.loads(stdout) != json.loads(cloud_analysis):
            raise ValueError("completed local replay log differs from cloud analysis")
    else:
        replay = subprocess.run([sys.executable, str(HERE / "verify.py"), str(results)],
                                cwd=ROOT, capture_output=True, text=True, check=False)
        (results.parent / "local_replay.stdout.txt").write_text(replay.stdout)
        (results.parent / "local_replay.stderr.txt").write_text(replay.stderr)
        (results.parent / "local_replay.exit_status.txt").write_text(
            str(replay.returncode) + "\n")
        if replay.returncode != 0:
            raise ValueError("local independent holdout replay failed")
    if (results / "analysis.json").read_bytes() != cloud_analysis:
        raise ValueError("local independent holdout replay differs from cloud")
    analysis = json.loads(cloud_analysis)
    source = panel.read(results / "source.json")
    host = panel.read(results / "host.json")
    challengers, selection = runner.finalists()
    expected = ("COLD_WINNER_CONFIRMED_WITHIN_GRID", "BOUNDED_UNKNOWN")
    if (analysis["schema"] != "n83.e2e-holdout-analysis/v1"
            or analysis["status"] not in expected
            or analysis["evidence_status"] != "PASS_INDEPENDENT_REPLAY"
            or (analysis["selected_pipeline"] is not None
                and analysis["selection_admitted"] is not True)
            or host["git_sha"] != args.head_sha
            or host["instance_id"] != f"gh-{args.run_id}-1-five-point-cold-holdout"
            or source["source_sha256"] != runner.screen.SOURCE_SHA256
            or source["hashes"]["Cargo.lock"] != panel.LOCK_SHA256
            or source["challenger_ks"] != challengers
            or source["selection_sha256"] != file_sha(runner.SELECTION)
            or source["screen_archive_sha256"] !=
            selection["screen_archive_sha256"]):
        raise ValueError("holdout source, selection or admission differs")
    points = [json.loads(line) for line in runner.TARGETS.read_text().splitlines()]
    answers = [json.loads(line) for line in runner.ANSWERS.read_text().splitlines()]
    if len(points) != 5 or len(answers) != 5:
        raise ValueError("five frozen points and answers required")
    field, curve = panel.curve_record()
    manifests, workloads, rows = {}, {}, []
    for point_row in analysis["point_rows"]:
        point_index = point_row["point_index"]
        index = point_index - 1
        point, answer = points[index], answers[index]
        rho = point_row["rho"]
        if (point_row["public_point"] != point or answer["point"] != point
                or answer["fixture_index"] != index
                or rho["status"] != "PASS_INDEPENDENT_REPLAY"
                or rho["resources"]["status"] != "PASS"):
            raise ValueError(f"point {point_index}: frozen fixture differs")
        workload_id, record = workload(index, point, answer, source,
                                       results / "inputs" / f"point_{point_index}.jsonl")
        workloads[workload_id] = {"workload_id": workload_id, **record}
        for item in point_row["arms"]:
            arm = item["arm"]
            k = item["k"]
            folder = results / "runs" / arm
            summary = panel.read(folder / "stdout.jsonl")
            if (item["status"] != "PASS_INDEPENDENT_REPLAY"
                    or item["resources"]["status"] != "PASS"
                    or item["actual_factor_base_points"] != 166 * k
                    or item["base_hash"] != summary["base_hash"]):
                raise ValueError(f"{arm}: complete replay and base identity required")
            manifest = dynamic.candidate(k, "root11_dynamic", summary,
                                         folder / "base_dump.jsonl", field, curve)
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
                "source_commit": runner.SCREEN_HEAD_SHA,
                "run_head_sha": args.head_sha,
                "executed_binary_sha256": source["binary_sha256"],
                "rho_binary_sha256": source["rho_binary_sha256"],
                "host_receipt_sha256": file_sha(results / "host.json"),
                "source_receipt_sha256": file_sha(results / "source.json"),
                "correctness_receipt_sha256": hashlib.sha256(cloud_analysis).hexdigest(),
                "raw_artifact_run_id": args.run_id,
                "raw_artifact_sha256": archive_sha,
                "cold_process_ms": item["cold_process_ms_exploratory"],
                "one_target_online_ms": item["one_target_online_ms"],
                "paired_rho_online_ms": rho["online_ms"],
                "rho_online_over_ic_online_ratio":
                    rho["online_ms"] / item["one_target_online_ms"],
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
    selected_id = None
    if analysis["selected_pipeline"] is not None:
        chosen_arm = "p1_" + analysis["selected_pipeline"]
        if analysis["selected_pipeline"] == "k1152":
            chosen_arm += "_a"
        chosen = [row["candidate_id"] for row in rows
                  if row["arm"] == chosen_arm]
        if len(chosen) != 1:
            raise ValueError("selected pipeline lacks exact candidate")
        selected_id = chosen[0]
    pipeline_summaries = {}
    rho_by_point = {row["point_index"]: row["rho"]["online_ms"]
                    for row in analysis["point_rows"]}
    for pipeline in ("k1152", *(f"k{k}" for k in challengers)):
        per_point = []
        for point_index in range(1, 6):
            prefix = f"p{point_index}_{pipeline}"
            matched = [row for row in rows if row["arm"] == prefix
                       or row["arm"] in (prefix + "_a", prefix + "_b")]
            expected_count = 2 if pipeline == "k1152" else 1
            if len(matched) != expected_count:
                raise ValueError(f"{prefix}: incomplete pipeline repetitions")
            if (len({row["candidate_id"] for row in matched}) != 1
                    or len({row["actual_usable_points"] for row in matched}) != 1
                    or len({row["effective_columns"] for row in matched}) != 1):
                raise ValueError(f"{prefix}: identity or base size differs")
            per_point.append({
                "point_index": point_index,
                "candidate_id": matched[0]["candidate_id"],
                "actual_usable_points": matched[0]["actual_usable_points"],
                "effective_columns": matched[0]["effective_columns"],
                "peak_rss_bytes": max(row["peak_rss_bytes"] for row in matched),
                "cold_process_ms": statistics.mean(
                    row["cold_process_ms"] for row in matched),
                "setup_ms": statistics.mean(row["setup_ms"] for row in matched),
                "rank_stage_ms": statistics.mean(
                    row["rank_stage_ms"] for row in matched),
                "targets_total_ms": statistics.mean(
                    row["targets_total_ms"] for row in matched),
                "one_target_online_ms": statistics.mean(
                    row["one_target_online_ms"] for row in matched),
                "rho_online_ms": rho_by_point[point_index],
            })
        pipeline_summaries[pipeline] = {
            "candidate_id": per_point[0]["candidate_id"],
            "actual_usable_points": per_point[0]["actual_usable_points"],
            "effective_columns": per_point[0]["effective_columns"],
            "max_peak_rss_bytes": max(point["peak_rss_bytes"]
                                      for point in per_point),
            "per_point": per_point,
            "mean_cold_process_ms": statistics.mean(
                point["cold_process_ms"] for point in per_point),
            "mean_setup_ms": statistics.mean(point["setup_ms"] for point in per_point),
            "mean_rank_stage_ms": statistics.mean(
                point["rank_stage_ms"] for point in per_point),
            "mean_targets_total_ms": statistics.mean(
                point["targets_total_ms"] for point in per_point),
            "mean_one_target_online_ms": statistics.mean(
                point["one_target_online_ms"] for point in per_point),
            "mean_rho_online_ms": statistics.mean(
                point["rho_online_ms"] for point in per_point),
        }
        pipeline_summaries[pipeline]["mean_cold_over_rho_online"] = (
            pipeline_summaries[pipeline]["mean_cold_process_ms"] /
            pipeline_summaries[pipeline]["mean_rho_online_ms"])
        if analysis["selection_admitted"]:
            frozen_mean = analysis["mean_cold_process_ms_by_pipeline"][pipeline]
            if abs(pipeline_summaries[pipeline]["mean_cold_process_ms"]
                   - frozen_mean) > 1e-6:
                raise ValueError(f"{pipeline}: materialized mean differs from replay")
    output = {
        "schema": "n83.e2e-holdout-ic1/v1",
        "source_commit": runner.SCREEN_HEAD_SHA,
        "run_head_sha": args.head_sha,
        "raw_artifact_sha256": archive_sha,
        "analysis_sha256": hashlib.sha256(cloud_analysis).hexdigest(),
        "workloads": [workloads[key] for key in sorted(workloads)],
        "candidates": sorted(manifests), "rows": rows,
        "selection_status": analysis["status"],
        "selected_candidate_id": selected_id,
        "evidence_status": analysis["evidence_status"],
        "selection_admitted": analysis["selection_admitted"],
        "hardware_tier": analysis["hardware_tier"],
        "pipeline_summaries": pipeline_summaries,
        "aa_controls": [row["aa_control"] for row in analysis["point_rows"]],
        "paired_intervals": analysis["paired_intervals"],
        "same_point_rho": [row["rho"] for row in analysis["point_rows"]],
    }
    for candidate_id, manifest in manifests.items():
        dynamic.write_once(HERE / f"candidate_{candidate_id}.json", manifest)
    dynamic.write_once(HERE / f"ic1_rows_{args.run_id}.json", output)
    print(json.dumps({"candidate_count": len(manifests),
                      "workload_count": len(workloads), "row_count": len(rows),
                      "selected_candidate_id": selected_id}, sort_keys=True))


if __name__ == "__main__":
    main()
