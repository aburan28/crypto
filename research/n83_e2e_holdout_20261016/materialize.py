#!/usr/bin/env python3
"""Bind a replayed N83 five-point holdout to exact IC1 candidate identities."""

import argparse
import hashlib
import importlib.util
import json
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
    args = parser.parse_args()
    results = args.results.resolve()
    cloud_analysis, archive_sha = archive_analysis(
        args.artifact_zip.resolve(), args.provider_sha256)
    if (results / "analysis.json").read_bytes() != cloud_analysis:
        raise ValueError("extracted holdout analysis differs from exact ZIP")
    replay = subprocess.run([sys.executable, str(HERE / "verify.py"), str(results)],
                            cwd=ROOT, capture_output=True, text=True, check=False)
    (results.parent / "local_replay.stdout.txt").write_text(replay.stdout)
    (results.parent / "local_replay.stderr.txt").write_text(replay.stderr)
    if replay.returncode != 0 or (results / "analysis.json").read_bytes() != cloud_analysis:
        raise ValueError("local independent holdout replay differs from cloud")
    analysis = json.loads(cloud_analysis)
    source = panel.read(results / "source.json")
    host = panel.read(results / "host.json")
    challengers, selection = runner.finalists()
    expected = ("COLD_WINNER_CONFIRMED_WITHIN_GRID", "BOUNDED_UNKNOWN")
    if (analysis["schema"] != "n83.e2e-holdout-analysis/v1"
            or analysis["status"] not in expected
            or analysis["evidence_status"] != "PASS_INDEPENDENT_REPLAY"
            or analysis["selection_admitted"] is not True
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
        if (point_row["public_point"] != point or answer["point"] != point
                or answer["fixture_index"] != index):
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
