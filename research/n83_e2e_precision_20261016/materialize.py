#!/usr/bin/env python3
"""Bind the replayed N83 precision panel to exact IC1 candidate identities."""

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
EXPECTED_CANDIDATES = {
    1104: "IC1N83Ckb1a1fb183264PDP4rootRCguidedLAgaussTDdirectISO0h6f30a77abf9e",
    1120: "IC1N83Ckb1a1fb185920PDP4rootRCguidedLAgaussTDdirectISO0hc922b499a234",
}


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    value = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(value)
    return value


run = module("n83_precision_identity_run", HERE / "run.py")
# Importing the earlier materializer also pins its canonical curve and
# candidate-record construction to the source used by this panel.
holdout = module("n83_precision_identity_holdout",
                 ROOT / "research/n83_e2e_holdout_20261016/materialize.py")

panel = holdout.panel
dynamic = holdout.dynamic


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def archive_analysis(archive, provider_sha):
    expected = provider_sha.removeprefix("sha256:")
    if (len(expected) != 64 or
            any(character not in "0123456789abcdef" for character in expected)):
        raise ValueError("provider digest must be a lowercase SHA-256")
    if sha256(archive) != expected:
        raise ValueError("precision ZIP differs from provider digest")
    with zipfile.ZipFile(archive) as saved:
        if saved.testzip() is not None:
            raise ValueError("precision ZIP integrity failed")
        names = saved.namelist()
        if len(names) != len(set(names)):
            raise ValueError("duplicate precision ZIP members")
        for name in names:
            if Path(name).is_absolute() or ".." in Path(name).parts:
                raise ValueError("unsafe precision ZIP member")
        return saved.read("results/analysis.json"), expected


def bind_extraction(archive, results):
    count = 0
    with zipfile.ZipFile(archive) as saved:
        for info in saved.infolist():
            if info.is_dir() or not info.filename.startswith("results/"):
                continue
            relative = Path(info.filename).relative_to("results")
            local = results / relative
            if not local.is_file():
                raise ValueError(f"missing extracted ZIP member: {relative}")
            digest = hashlib.sha256()
            with saved.open(info) as stream:
                for block in iter(lambda: stream.read(1024 * 1024), b""):
                    digest.update(block)
            if digest.hexdigest() != sha256(local):
                raise ValueError(f"extracted ZIP member differs: {relative}")
            count += 1
    if count < 20:
        raise ValueError("precision ZIP lacks complete result members")
    return count


def workload(index, point, answer, input_path):
    record = {
        "schema_version": 1,
        "curve_id": panel.CURVE_ID,
        "target_count": 1,
        "public_targets": [point],
        "input_law": {
            "domain": run.generator.DOMAIN.decode(),
            "fixture_index": index,
            "scalar_rule": "1+(little_endian(SHA256(domain||u32_le(index))) mod (r-1))",
            "seed_sha256": answer["seed_sha256"],
            "generator_sha256": sha256(HERE / "generate.py"),
            "answer_sidecar_sha256": sha256(run.ANSWERS),
        },
        "target_sha256": sha256(input_path),
        "cold_or_warm": "one_public_point_per_fresh_process_empty_base_index_and_logs",
    }
    return panel.digest(record)[:12], record


def replay(results, cloud_analysis, reuse):
    if reuse:
        status = (results.parent / "local_replay.exit_status.txt").read_text().strip()
        stdout = (results.parent / "local_replay.stdout.txt").read_text()
        stderr = (results.parent / "local_replay.stderr.txt").read_text()
        if (status != "0" or stderr or
                json.loads(stdout) != json.loads(cloud_analysis)):
            raise ValueError("local replay log differs from cloud analysis")
    else:
        process = subprocess.run(
            [sys.executable, str(HERE / "verify.py"), str(results)],
            cwd=ROOT, capture_output=True, text=True, check=False)
        (results.parent / "local_replay.stdout.txt").write_text(process.stdout)
        (results.parent / "local_replay.stderr.txt").write_text(process.stderr)
        (results.parent / "local_replay.exit_status.txt").write_text(
            f"{process.returncode}\n")
        if process.returncode != 0:
            raise ValueError("independent precision replay did not admit selection")
    if (results / "analysis.json").read_bytes() != cloud_analysis:
        raise ValueError("local precision replay differs from provider ZIP")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("results", type=Path)
    parser.add_argument("--artifact-zip", type=Path, required=True)
    parser.add_argument("--provider-sha256", required=True)
    parser.add_argument("--run-id", type=int, required=True)
    parser.add_argument("--head-sha", required=True)
    parser.add_argument("--reuse-local-replay", action="store_true")
    args = parser.parse_args()
    results = args.results.resolve()
    cloud_analysis, archive_sha = archive_analysis(
        args.artifact_zip.resolve(), args.provider_sha256)
    if (results / "analysis.json").read_bytes() != cloud_analysis:
        raise ValueError("extracted analysis differs from provider ZIP")
    raw_member_count = bind_extraction(args.artifact_zip.resolve(), results)
    replay(results, cloud_analysis, args.reuse_local_replay)
    analysis = json.loads(cloud_analysis)
    source = panel.read(results / "source.json")
    host = panel.read(results / "host.json")
    hashes, fixture = run.fixture()
    if (analysis["schema"] != "n83.e2e-precision-analysis/v1"
            or analysis["status"] not in
            ("COLD_WINNER_CONFIRMED_WITHIN_PAIR", "BOUNDED_UNKNOWN")
            or analysis["evidence_status"] != "PASS_INDEPENDENT_REPLAY"
            or analysis["selection_admitted"] is not True
            or analysis["point_count"] != 16
            or analysis["ic_arm_count"] != 64
            or analysis["rho_arm_count"] != 16
            or len(analysis["point_rows"]) != 16
            or host["git_sha"] != args.head_sha
            or host["instance_id"] !=
            f"gh-{args.run_id}-1-sixteen-point-cold-precision"
            or source["hashes"] != hashes
            or source["source_sha256"] != run.screen.SOURCE_SHA256
            or source["runner_sha256"] != sha256(HERE / "run.py")
            or source["targets_sha256"] != fixture["targets_sha256"]
            or source["answers_sha256"] != fixture["answers_sha256"]
            or source["schedule"] != run.schedule()):
        raise ValueError("precision source, host, fixture or admission differs")
    points = [json.loads(line) for line in run.TARGETS.read_text().splitlines()]
    answers = [json.loads(line) for line in run.ANSWERS.read_text().splitlines()]
    field, curve = panel.curve_record()
    manifests, workloads, rows = {}, {}, []
    for position, point_row in enumerate(analysis["point_rows"], start=1):
        point, answer = points[position - 1], answers[position - 1]
        rho = point_row["rho"]
        if (point_row["point_index"] != position
                or point_row["public_point"] != point
                or answer["point"] != point
                or answer["fixture_index"] != position - 1
                or rho["status"] != "PASS_INDEPENDENT_REPLAY"
                or rho["resources"]["status"] != "PASS"):
            raise ValueError(f"point {position}: fixture or rho differs")
        workload_id, record = workload(
            position - 1, point, answer,
            results / "inputs" / f"point_{position}.jsonl")
        workloads[workload_id] = {"workload_id": workload_id, **record}
        if len(point_row["arms"]) != 4:
            raise ValueError(f"point {position}: four IC arms required")
        for item in point_row["arms"]:
            arm, k = item["arm"], item["k"]
            folder = results / "runs" / arm
            summary = panel.read(folder / "stdout.jsonl")
            if (arm not in {f"p{position}_k{value}_{suffix}"
                            for value in run.K_VALUES for suffix in ("a", "b")}
                    or item["status"] != "PASS_INDEPENDENT_REPLAY"
                    or item["resources"]["status"] != "PASS"
                    or item["actual_factor_base_points"] != 166 * k
                    or item["base_hash"] != summary["base_hash"]):
                raise ValueError(f"{arm}: replay or base differs")
            manifest = dynamic.candidate(k, "root11_dynamic", summary,
                                         folder / "base_dump.jsonl", field, curve)
            candidate_id = manifest["candidate_id"]
            if candidate_id != EXPECTED_CANDIDATES[k]:
                raise ValueError(f"{arm}: candidate differs from five-point holdout")
            if candidate_id in manifests and manifests[candidate_id] != manifest:
                raise ValueError(f"{arm}: candidate identity collision")
            manifests[candidate_id] = manifest
            repetition = 1 if arm.endswith("_a") else 2
            rows.append({
                "arm": arm, "candidate_id": candidate_id,
                "workload_id": workload_id,
                "run_id": f"{candidate_id}W{workload_id}R{repetition}",
                "status": item["status"],
                "hardware_tier": analysis["hardware_tier"],
                "source_commit": holdout.runner.SCREEN_HEAD_SHA,
                "run_head_sha": args.head_sha,
                "executed_binary_sha256": source["binary_sha256"],
                "rho_binary_sha256": source["rho_binary_sha256"],
                "host_receipt_sha256": sha256(results / "host.json"),
                "source_receipt_sha256": sha256(results / "source.json"),
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
    if len(manifests) != 2 or len(workloads) != 16 or len(rows) != 64:
        raise ValueError("precision candidate, workload or run counts differ")
    by_pipeline = {}
    for k in run.K_VALUES:
        pipeline = f"k{k}"
        per_point = []
        for position in range(1, 17):
            matched = [row for row in rows
                       if row["arm"] in (f"p{position}_{pipeline}_a",
                                         f"p{position}_{pipeline}_b")]
            if len(matched) != 2 or matched[0]["candidate_id"] != matched[1]["candidate_id"]:
                raise ValueError(f"point {position}: repeated candidate differs")
            per_point.append({
                "point_index": position,
                "candidate_id": matched[0]["candidate_id"],
                "cold_process_ms": statistics.mean(row["cold_process_ms"] for row in matched),
                "setup_ms": statistics.mean(row["setup_ms"] for row in matched),
                "rank_stage_ms": statistics.mean(row["rank_stage_ms"] for row in matched),
                "targets_total_ms": statistics.mean(row["targets_total_ms"] for row in matched),
                "one_target_online_ms": statistics.mean(
                    row["one_target_online_ms"] for row in matched),
                "rho_online_ms": matched[0]["paired_rho_online_ms"],
                "peak_rss_bytes": max(row["peak_rss_bytes"] for row in matched),
            })
        means = {name: statistics.mean(point[name] for point in per_point)
                 for name in ("cold_process_ms", "setup_ms", "rank_stage_ms",
                              "targets_total_ms", "one_target_online_ms",
                              "rho_online_ms")}
        by_pipeline[pipeline] = {
            "candidate_id": per_point[0]["candidate_id"],
            "actual_usable_points": 166 * k,
            "effective_columns": next(row["effective_columns"] for row in rows
                                      if row["candidate_id"] == per_point[0]["candidate_id"]),
            "max_peak_rss_bytes": max(point["peak_rss_bytes"] for point in per_point),
            "per_point": per_point, "mean_ms": means,
            "mean_cold_over_rho_online": means["cold_process_ms"] / means["rho_online_ms"],
        }
    chosen = analysis["selected_pipeline"]
    selected_id = by_pipeline[chosen]["candidate_id"] if chosen else None
    output = {
        "schema": "n83.e2e-precision-ic1/v1",
        "source_commit": holdout.runner.SCREEN_HEAD_SHA,
        "run_head_sha": args.head_sha,
        "raw_artifact_run_id": args.run_id,
        "raw_artifact_sha256": archive_sha,
        "bound_result_members": raw_member_count,
        "analysis_sha256": hashlib.sha256(cloud_analysis).hexdigest(),
        "workloads": [workloads[key] for key in sorted(workloads)],
        "candidates": sorted(manifests), "rows": rows,
        "selection_status": analysis["status"],
        "selected_candidate_id": selected_id,
        "evidence_status": analysis["evidence_status"],
        "selection_admitted": analysis["selection_admitted"],
        "hardware_tier": analysis["hardware_tier"],
        "pipeline_summaries": by_pipeline,
        "aa_controls": [point["aa_control"] for point in analysis["point_rows"]],
        "paired_interval": analysis["paired_interval"],
        "same_point_rho": [point["rho"] for point in analysis["point_rows"]],
    }
    for candidate_id, manifest in manifests.items():
        dynamic.write_once(HERE / f"candidate_{candidate_id}.json", manifest)
    dynamic.write_once(HERE / f"ic1_rows_{args.run_id}.json", output)
    print(json.dumps({"candidate_count": len(manifests),
                      "workload_count": len(workloads), "row_count": len(rows),
                      "selected_candidate_id": selected_id}, sort_keys=True))


if __name__ == "__main__":
    main()
