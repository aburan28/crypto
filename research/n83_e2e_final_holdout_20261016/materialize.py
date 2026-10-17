#!/usr/bin/env python3
"""Materialize provider-bound N83 finalist IC1 manifests and per-run rows."""

import argparse
import hashlib
import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
DOMAIN = "n83-e2e-final-holdout-sixteen-20261016-v1"
SCALAR_RULE = "1+(little_endian(SHA256(domain||u32_le(index))) mod (r-1))"


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    value = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(value)
    return value


runner = module("n83_final_materialize_runner", HERE / "run.py")
selector = module("n83_final_materialize_selector", HERE / "select.py")
finegrid = module("n83_final_materialize_finegrid",
                  ROOT / "research/n83_e2e_finegrid_20261016/materialize.py")


def read(path):
    return json.loads(path.read_text())


def sha256(path):
    value = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def write_once(path, value):
    blob = json.dumps(value, indent=2, sort_keys=True) + "\n"
    if path.exists() and path.read_text() != blob:
        raise ValueError(f"sealed materialized file differs: {path}")
    path.write_text(blob)


def workload(index, encoded, answer, point_file, field_curve_id):
    record = {
        "schema_version": 1, "curve_id": field_curve_id,
        "target_count": 1, "public_targets": [encoded],
        "input_law": {
            "domain": DOMAIN, "fixture_index": index,
            "scalar_rule": SCALAR_RULE,
            "seed_sha256": answer["seed_sha256"],
            "generator_sha256": sha256(HERE / "generate.py"),
            "answer_sidecar_sha256": sha256(HERE / "answers.jsonl"),
        },
        "target_sha256": sha256(point_file),
        "cold_or_warm": "one_public_point_per_fresh_process_empty_base_index_and_logs",
    }
    return {"workload_id": finegrid.panel.digest(record)[:12], **record}


def row(item, candidate_id, workload_id, repetition, point_index, half,
        artifact_sha, run_id, source, host, analysis_path, control,
        tier_status):
    passed = (item["status"] == "PASS_INDEPENDENT_REPLAY"
              and item.get("resources", {}).get("status") == "PASS")
    online = item.get("one_target_online_ms") if passed else None
    return {
        "candidate_id": candidate_id, "workload_id": workload_id,
        "run_id": f"{candidate_id}W{workload_id}R{repetition}",
        "point_index": point_index, "stratum_half": half,
        "control_repeat": control,
        "status": item["status"],
        "hardware_tier": "cloud_instance_tier1" if passed and
                           tier_status == "PASS_TIER1" else
                         "cloud_instance_exploratory",
        "source_snapshot_sha256": source["source_sha256"],
        "run_head_sha": host["git_sha"],
        "executed_binary_sha256": source["binary_sha256"],
        "rho_binary_sha256": source["rho_binary_sha256"],
        "host_receipt_sha256": sha256(host["path"]),
        "source_receipt_sha256": sha256(source["path"]),
        "correctness_receipt_sha256": sha256(analysis_path),
        "raw_artifact_run_id": run_id,
        "raw_artifact_sha256": artifact_sha,
        "cold_process_ms": item.get("cold_process_ms_exploratory") if passed else None,
        "one_target_online_ms": online,
        "paired_rho_online_ms": None,
        "online_phases_ms": item.get("online_phases_ms") if passed else None,
        "setup_ms": item.get("setup_ms") if passed else None,
        "rank_stage_ms": item.get("rank_stage_ms") if passed else None,
        "targets_total_ms": item.get("targets_total_ms") if passed else None,
        "rank_attempts": item.get("rank_attempts"),
        "rank_relations": item.get("rank_relations"),
        "target_probes": item.get("target_probes"),
        "actual_usable_points": item.get("actual_factor_base_points"),
        "effective_columns": item.get("effective_columns"),
        "peak_rss_bytes": item.get("resources", {}).get("peak_rss_bytes"),
        "outer_wall_seconds": item.get("full_child_wall_seconds"),
        "raw_evidence": item,
    }


def materialize(validated_path, selection_path, half_dirs):
    validated, selected = read(validated_path), read(selection_path)
    if (validated["schema"] != "n83.e2e-final-validated/v1"
            or validated["source_sha256"] != runner.SOURCE_SHA256
            or validated["targets_sha256"] != runner.TARGETS_SHA256
            or validated["answers_sha256"] != runner.ANSWERS_SHA256
            or selected != json.loads(json.dumps(selector.select(validated)))):
        raise ValueError("validated rows or selection differ from frozen rule")
    if (len(half_dirs) != 2 or len(validated["host_strata"]) != 2
            or len(validated["points"]) != 16):
        raise ValueError("two host strata and sixteen point rows required")
    target_lines = (HERE / "targets.jsonl").read_text().splitlines(keepends=True)
    encoded = [json.loads(line) for line in target_lines]
    answers = [json.loads(line) for line in (HERE / "answers.jsonl").read_text().splitlines()]
    if len(encoded) != 16 or len(answers) != 16:
        raise ValueError("sealed finalist fixture count differs")
    field, curve = finegrid.panel.curve_record()
    by_id = {arm[3]: arm for arm in runner.ARMS}
    by_name = {arm[0]: arm for arm in runner.ARMS}
    manifests, rows, rho_rows, workloads = {}, [], [], []
    half_analysis = {}
    for half, directory in enumerate(half_dirs):
        analysis_path = directory / "analysis.json"
        analysis = read(analysis_path)
        host = read(directory / "host.json")
        source = read(directory / "source.json")
        stratum = validated["host_strata"][half]
        if (analysis["schema"] != "n83.e2e-finalist-half-analysis/v1"
                or analysis["half"] != half
                or stratum["half"] != half
                or stratum["analysis_sha256"] != sha256(analysis_path)
                or host["instance_id"] != stratum["instance_id"]
                or source["source_sha256"] != runner.SOURCE_SHA256):
            raise ValueError(f"half {half}: local raw receipts differ")
        half_analysis[half] = analysis
        host["path"], source["path"] = (directory / "host.json",
                                         directory / "source.json")
        for point in analysis["points"]:
            index = point["point_index"]
            if index not in range(half * 8, half * 8 + 8):
                raise ValueError("point assigned to wrong host stratum")
            raw_point = directory / "points" / f"p{index:02d}"
            if (raw_point / "target.jsonl").read_text() != target_lines[index]:
                raise ValueError("per-point public input differs")
            workload_record = workload(index, encoded[index], answers[index],
                                       raw_point / "target.jsonl",
                                       finegrid.panel.CURVE_ID)
            workloads.append(workload_record)
            workload_id = workload_record["workload_id"]
            for arm in point["arms"]:
                candidate_id = arm["candidate_id"]
                if candidate_id not in by_id:
                    raise ValueError("unexpected finalist candidate ID")
                rho = point["rho_by_candidate"][candidate_id]
                calibration = point["calibration_by_candidate"][candidate_id]
                rho_rows.append({"point_index": index,
                                 "paired_candidate_id": candidate_id,
                                 "workload_id": workload_id,
                                 "status": rho["status"],
                                 "online_ms": rho.get("online_ms"),
                                 "calibration": calibration,
                                 "host_instance_id": host["instance_id"],
                                 "raw_artifact_sha256": stratum["raw_artifact_sha256"],
                                 "raw_evidence": rho})
                name, k, layout, _ = by_id[candidate_id]
                item = arm["raw_replay"]
                if item["status"] == "PASS_INDEPENDENT_REPLAY":
                    folder = raw_point / "runs" / name
                    manifest = finegrid.dynamic.candidate(
                        k, layout, read(folder / "stdout.jsonl"),
                        folder / "base_dump.jsonl", field, curve)
                    if manifest["candidate_id"] != candidate_id:
                        raise ValueError("raw factor-base manifest identity differs")
                    if candidate_id in manifests and manifests[candidate_id] != manifest:
                        raise ValueError("candidate manifest differs between points")
                    manifests[candidate_id] = manifest
                run_row = row(item, candidate_id, workload_id, 1, index, half,
                              stratum["raw_artifact_sha256"], validated["run_id"],
                              source, host, analysis_path, False,
                              stratum["tier_status"])
                run_row["paired_rho_online_ms"] = rho.get("online_ms")
                run_row["pair_calibration"] = calibration
                rows.append(run_row)
        first_index = half * 8
        first_workload = workloads[first_index]["workload_id"]
        for name, item in analysis["control_arms"].items():
            if name not in runner.CONTROL_KEYS:
                raise ValueError("unexpected control arm")
            candidate_id = by_name[name][3]
            run_row = row(item, candidate_id, first_workload, 2,
                          first_index, half, stratum["raw_artifact_sha256"],
                          validated["run_id"], source, host, analysis_path, True,
                          stratum["tier_status"])
            run_row["paired_rho_online_ms"] = None
            rows.append(run_row)
    if (len(rows) != 132 or len(rho_rows) != 128 or len(workloads) != 16
            or len({item["workload_id"] for item in workloads}) != 16):
        raise ValueError("finalist run or workload count differs")
    if validated["evidence_status"] == "PASS_INDEPENDENT_REPLAY" and set(manifests) != set(selector.IDS):
        raise ValueError("admitted evidence lacks an exact candidate manifest")
    output = {
        "schema": "n83.e2e-finalist-ic1/v1",
        "run_id": validated["run_id"],
        "run_attempt": validated["run_attempt"],
        "run_head_sha": validated["head_sha"],
        "validated_rows_sha256": sha256(validated_path),
        "selection_sha256": sha256(selection_path),
        "host_strata": validated["host_strata"],
        "workloads": workloads,
        "candidates": sorted(manifests),
        "rows": rows, "rho_rows": rho_rows,
        "evidence_status": validated["evidence_status"],
        "online_selection": selected["online"],
        "cold_selection": selected["cold"],
    }
    for candidate_id, manifest in manifests.items():
        write_once(HERE / f"candidate_{candidate_id}.json", manifest)
    write_once(HERE / f"ic1_rows_{validated['run_id']}.json", output)
    return {"candidate_count": len(manifests), "row_count": len(rows),
            "rho_count": len(rho_rows)}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("validated_rows", type=Path)
    parser.add_argument("selection", type=Path)
    parser.add_argument("half0_results", type=Path)
    parser.add_argument("half1_results", type=Path)
    args = parser.parse_args()
    result = materialize(args.validated_rows.resolve(), args.selection.resolve(),
                         (args.half0_results.resolve(), args.half1_results.resolve()))
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
