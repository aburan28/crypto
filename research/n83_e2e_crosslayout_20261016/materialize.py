#!/usr/bin/env python3
"""Bind a replayed N83 cross-layout screen to exact IC1 records."""

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
EXPECTED_ROOT11 = {
    1104: "IC1N83Ckb1a1fb183264PDP4rootRCguidedLAgaussTDdirectISO0h6f30a77abf9e",
    1120: "IC1N83Ckb1a1fb185920PDP4rootRCguidedLAgaussTDdirectISO0hc922b499a234",
}


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    value = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(value)
    return value


runner = module("n83_crosslayout_identity_run", HERE / "run.py")
prior = module("n83_crosslayout_prior_identity",
               ROOT / "research/n83_e2e_holdout_20261016/materialize.py")
panel = prior.panel
dynamic = prior.dynamic


def sha256(path):
    value = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def archive_analysis(archive, provider_sha):
    expected = provider_sha.removeprefix("sha256:")
    if (len(expected) != 64 or
            any(character not in "0123456789abcdef" for character in expected)):
        raise ValueError("provider digest must be lowercase SHA-256")
    if sha256(archive) != expected:
        raise ValueError("cross-layout ZIP differs from provider digest")
    with zipfile.ZipFile(archive) as saved:
        if saved.testzip() is not None:
            raise ValueError("cross-layout ZIP integrity failed")
        names = saved.namelist()
        if len(names) != len(set(names)):
            raise ValueError("duplicate ZIP members")
        for name in names:
            if Path(name).is_absolute() or ".." in Path(name).parts:
                raise ValueError("unsafe ZIP member")
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
        raise ValueError("cross-layout ZIP lacks complete result members")
    return count


def replay(results, cloud_analysis, reuse):
    if reuse:
        status = (results.parent / "local_replay.exit_status.txt").read_text().strip()
        stdout = (results.parent / "local_replay.stdout.txt").read_text()
        stderr = (results.parent / "local_replay.stderr.txt").read_text()
        if status != "0" or stderr or json.loads(stdout) != json.loads(cloud_analysis):
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
            raise ValueError("local cross-layout replay did not pass")
    if (results / "analysis.json").read_bytes() != cloud_analysis:
        raise ValueError("local cross-layout replay differs from provider ZIP")


def workload(source):
    point = json.loads(runner.cold.TARGETS.read_text())
    answer = json.loads(runner.cold.ANSWERS.read_text())
    if answer["fixture_index"] != 0 or answer["point"] != point:
        raise ValueError("exposed target or answer sidecar differs")
    record = {
        "schema_version": 1, "curve_id": panel.CURVE_ID,
        "target_count": 1, "public_targets": [point],
        "input_law": {
            "domain": "n83-cold-factor-base-one-plus-four-20261016-v1",
            "fixture_index": 0,
            "scalar_rule": "1+(little_endian(SHA256(domain||u32_le(index))) mod (r-1))",
            "seed_sha256": answer["seed_sha256"],
            "generator_sha256": source["hashes"][
                "research/n83_cold_factor_base_20261016/generate.py"],
            "answer_sidecar_sha256": sha256(runner.cold.ANSWERS),
        },
        "target_sha256": sha256(runner.cold.TARGETS),
        "cold_or_warm": "one_public_point_per_fresh_process_empty_base_index_and_logs",
    }
    identity = panel.digest(record)[:12]
    if identity != "eebf98ee1031":
        raise ValueError("exposed one-point workload identity differs")
    return {"workload_id": identity, **record}


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
    if (analysis["schema"] != "n83.e2e-crosslayout-analysis/v1"
            or analysis["status"] != "SCREEN_COMPLETE"
            or analysis["evidence_status"] != "PASS_INDEPENDENT_REPLAY"
            or analysis["selection_admitted"] is not True
            or analysis["selected_candidate_id"] is not None
            or len(analysis["arms"]) != len(runner.ARMS)
            or host["git_sha"] != args.head_sha
            or host["instance_id"] != f"gh-{args.run_id}-1-exposed-point-crosslayout"
            or source["hashes"] != runner.frozen_inputs()
            or source["source_sha256"] != runner.SOURCE_SHA256
            or source["runner_sha256"] != sha256(HERE / "run.py")
            or source["target_sha256"] != runner.TARGET_SHA256
            or source["arm_order"] != [name for name, _, _ in runner.ARMS]):
        raise ValueError("cross-layout source, host, target or admission differs")
    workload_record = workload(source)
    workload_id = workload_record["workload_id"]
    field, curve = panel.curve_record()
    by_arm = {item["arm"]: item for item in analysis["arms"]}
    manifests, rows = {}, []
    for arm, k, layout in runner.ARMS:
        item = by_arm[arm]
        folder = results / "runs" / arm
        summary = panel.read(folder / "stdout.jsonl")
        if (item["status"] != "PASS_INDEPENDENT_REPLAY"
                or item["resources"]["status"] != "PASS"
                or item["base_hash"] != summary["base_hash"]
                or item["actual_factor_base_points"] != 166 * k):
            raise ValueError(f"{arm}: replay or factor-base identity differs")
        manifest = dynamic.candidate(k, layout, summary,
                                     folder / "base_dump.jsonl", field, curve)
        candidate_id = manifest["candidate_id"]
        if layout == "root11_dynamic" and candidate_id != EXPECTED_ROOT11[k]:
            raise ValueError(f"{arm}: root11 candidate differs from holdout")
        if candidate_id in manifests and manifests[candidate_id] != manifest:
            raise ValueError(f"{arm}: candidate identity collision")
        manifests[candidate_id] = manifest
        repetition = 2 if arm.endswith("_b") else 1
        rows.append({
            "arm": arm, "candidate_id": candidate_id,
            "workload_id": workload_id,
            "run_id": f"{candidate_id}W{workload_id}R{repetition}",
            "status": item["status"],
            "hardware_tier": analysis["hardware_tier"],
            "source_commit": prior.runner.SCREEN_HEAD_SHA,
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
            "paired_rho_online_ms": analysis["rho"]["online_ms"],
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
            "raw_evidence": item,
        })
    if len(manifests) != 4 or len(rows) != 6:
        raise ValueError("four cross-layout candidates and six arms required")
    by_pipeline = {}
    for pipeline, measured in analysis["descriptive_cold_process_ms_by_pipeline"].items():
        candidate_rows = [row for row in rows
                          if pipeline == f"k{next(item['k'] for item in analysis['arms'] if item['arm'] == row['arm'])}"
                          + next(item["layout"] for item in analysis["arms"]
                                 if item["arm"] == row["arm"])]
        if not candidate_rows:
            raise ValueError(f"{pipeline}: no exact run row")
        if abs(statistics.mean(row["cold_process_ms"] for row in candidate_rows)
               - measured) > 1e-6:
            raise ValueError(f"{pipeline}: replayed mean differs")
        by_pipeline[pipeline] = candidate_rows[0]["candidate_id"]
    output = {
        "schema": "n83.e2e-crosslayout-ic1/v1",
        "source_commit": prior.runner.SCREEN_HEAD_SHA,
        "run_head_sha": args.head_sha,
        "raw_artifact_run_id": args.run_id,
        "raw_artifact_sha256": archive_sha,
        "bound_result_members": raw_member_count,
        "analysis_sha256": hashlib.sha256(cloud_analysis).hexdigest(),
        "workload": workload_record,
        "candidates": sorted(manifests), "rows": rows,
        "pipeline_candidate_ids": by_pipeline,
        "descriptive_order": analysis["descriptive_order"],
        "selection_status": analysis["status"],
        "selected_candidate_id": None,
        "evidence_status": analysis["evidence_status"],
        "hardware_tier": analysis["hardware_tier"],
        "aa_controls": analysis["aa_controls"],
        "same_point_rho": analysis["rho"],
    }
    for candidate_id, manifest in manifests.items():
        dynamic.write_once(HERE / f"candidate_{candidate_id}.json", manifest)
    dynamic.write_once(HERE / f"ic1_rows_{args.run_id}.json", output)
    print(json.dumps({"candidate_count": len(manifests),
                      "row_count": len(rows), "workload_id": workload_id},
                     sort_keys=True))


if __name__ == "__main__":
    main()
