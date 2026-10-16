#!/usr/bin/env python3
"""Run two fresh-process IC/rho pairs on one unseen N83 public point."""

import argparse
import importlib.util
import json
import os
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
CORE = REPO / "research/n83_final_tier1_20261015/run.py"
CONFIRM = REPO / "research/n83_edge_confirmation_20261015/run.py"
RESULT = REPO / "research/n83_edge_confirmation_results_20261015/selection.json"
TARGETS = HERE / "targets.jsonl"
ANSWERS = HERE / "answers.jsonl"
TARGET_SHA256 = "dfea95ee4e423dd627adbd116828e9f478779c51476f0fa172cf4c160e645d12"
ANSWER_SHA256 = "ef4d939a30fbf26809f97acc0607fcd4fe04661e658e8b353c1ec77f5075d335"
ANALYSIS_SHA256 = "81a945b5711cf599996aa127aa3fbec3128ca8f523a8bad0b2fc8a4244e7359a"
ADMISSION_SHA256 = "803df5758dad2e7d4769554a6d6eb741c61ce7b8c8b1887c843f1616fcc9e1d8"
SOURCE_SHA256 = "f2745227d95e36ef45c73c5104a23ada6cebb0794405368581d2be17044845dc"
RHO_SHA256 = "fb14072f4fa5d08e18d8b79130b6d58a220adae0f76369cc9be949691c83eaa0"
PIPELINES = ("root4k3425std", "root11k3450wide3")

spec = importlib.util.spec_from_file_location("n83_diagnostic_core", CORE)
core = importlib.util.module_from_spec(spec)
spec.loader.exec_module(core)
cspec = importlib.util.spec_from_file_location("n83_diagnostic_confirmation", CONFIRM)
confirm = importlib.util.module_from_spec(cspec)
cspec.loader.exec_module(confirm)
parent = core.parent


def frozen():
    hashes = parent.check_inputs()
    if (hashes["examples/koblitz_orbit_dlp_fast.rs"] != SOURCE_SHA256
            or hashes["examples/n83_rho_optimization_eval.rs"] != RHO_SHA256
            or parent.digest(TARGETS) != TARGET_SHA256
            or parent.digest(ANSWERS) != ANSWER_SHA256
            or parent.digest(confirm.ADMISSION) != ADMISSION_SHA256):
        raise ValueError("one-target fixture, source or admission differs")
    admission, _ = confirm.frozen()
    result = json.loads(RESULT.read_text())
    if (result["schema"] != "n83.edge-confirmation-selection/v1"
            or result["status"] != "BOUNDED_UNKNOWN"
            or result["selected_pipeline"] is not None
            or result["analysis_sha256"] != ANALYSIS_SHA256
            or result["admission_sha256"] != ADMISSION_SHA256
            or result["balanced_comparison"]["descriptive_lowest_pipeline"] != PIPELINES[0]
            or tuple(row["pipeline"] for row in admission["candidate_order"]) != PIPELINES):
        raise ValueError("source finalist comparison differs")
    points = [json.loads(line) for line in TARGETS.read_text().splitlines() if line]
    answers = [json.loads(line) for line in ANSWERS.read_text().splitlines() if line]
    if len(points) != 1 or len(answers) != 1:
        raise ValueError("one-target fixture is incomplete")
    candidates = []
    for row in admission["candidate_order"]:
        candidates.append({
            "pipeline": row["pipeline"], "candidate_id": row["candidate_id"],
            "candidate_manifest_sha256": row["candidate_manifest_sha256"],
            "k": row["k"], "variant": row["variant"],
            "bloom_profile": row["bloom_profile"],
            "source_base_hash_blake3": row["base_hash"],
            "root_table_entries": row["root_table_entries"],
            "actual_usable_points": row["actual_usable_points"],
        })
    hashes[str(TARGETS.relative_to(REPO))] = TARGET_SHA256
    hashes[str(ANSWERS.relative_to(REPO))] = ANSWER_SHA256
    hashes[str(confirm.ADMISSION.relative_to(REPO))] = ADMISSION_SHA256
    hashes[str(RESULT.relative_to(REPO))] = parent.digest(RESULT)
    return hashes, points[0], answers[0], candidates


def run(binary, rho_binary, output):
    hashes, point, _, candidates = frozen()
    if output.exists():
        raise ValueError("output directory already exists")
    output.mkdir(parents=True)
    group = parent.cgroup_root()
    host = parent.host_receipt(group)
    host["schema"] = "n83.unresolved-tier1-host/v1"
    host["region_or_zone"] = (os.environ.get("RUNNER_REGION") or
                              os.environ.get("AZURE_REGION"))
    host["region_or_zone_status"] = (
        "EXPOSED" if host["region_or_zone"] else "NOT_EXPOSED_BY_RUNNER")
    execution_cpus, topology = core.select_execution_cpus(host["visible_affinity"] or [])
    host["execution_cpus"] = execution_cpus
    host["execution_topology"] = topology
    host["clock_ticks_per_second"] = os.sysconf("SC_CLK_TCK")
    host["affinity_set_status"] = "NOT_ATTEMPTED"
    if len(execution_cpus) == 4 and hasattr(os, "sched_setaffinity"):
        try:
            os.sched_setaffinity(0, set(execution_cpus))
            host["affinity_set_status"] = (
                "PASS" if set(os.sched_getaffinity(0)) == set(execution_cpus)
                else "MISMATCH")
        except OSError as error:
            host["affinity_set_status"] = f"FAILED:{error}"
    parent.write_json(output / "host.json", host)
    parent.write_json(output / "source.json", {
        "schema": "n83.unresolved-tier1-source/v1", "hashes": hashes,
        "candidate_order": candidates,
        "confirmation_selection_sha256": parent.digest(RESULT),
        "ic_binary_sha256": parent.digest(binary),
        "rho_binary_sha256": parent.digest(rho_binary),
        "runner_sha256": parent.digest(Path(__file__)),
    })
    memory_max = (host["cgroup_before"]["memory.max"] or "").strip()
    memory_too_small = memory_max.isdigit() and int(memory_max) < parent.RSS_CAP_BYTES
    if (host["mem_total_kib"] is None
            or host["mem_total_kib"] * 1024 < parent.RSS_CAP_BYTES + 256 * 1024**2
            or host["disk_free_bytes"] < 2 * 1024**3
            or len(execution_cpus) < 4 or memory_too_small
            or host["affinity_set_status"] != "PASS"):
        parent.write_json(output / "preflight.json", {
            "status": "REJECTED_HOST_CAPACITY", "host": host})
        return 2
    parent.write_json(output / "preflight.json", {
        "status": "PASS_CAPACITY_PREFLIGHT", "rss_cap_bytes": parent.RSS_CAP_BYTES})
    inputs = output / "inputs"
    inputs.mkdir()
    point_file = inputs / "point_1.jsonl"
    point_file.write_text(json.dumps(point, separators=(",", ":")) + "\n")
    for index, selected in enumerate(candidates):
        pair = output / "pairs" / selected["pipeline"]
        pair.mkdir(parents=True)
        parent.write_json(pair / "calibration_before.json", core.calibrate())
        parent.write_json(pair / "system_before.json", core.system_counters())
        order = ("ic", "rho") if index == 0 else ("rho", "ic")
        parent.write_json(pair / "order.json", {
            "candidate_id": selected["candidate_id"],
            "pipeline": selected["pipeline"], "order": order,
            "point_sha256": parent.digest(point_file), "started_utc": parent.now()})
        for arm in order:
            folder = pair / arm
            if arm == "ic":
                command = [str(binary), f"construct:83:1:{selected['k']}",
                           str(point_file), "29", str(folder / "targets_result.jsonl")]
                environment = core.ic_environment(folder, selected)
                cap = 1500
            else:
                command = [str(rho_binary), str(point_file),
                           str(2026101401 + index), str(parent.RHO_CONFIG),
                           str(folder / "receipt.json")]
                environment = {"PATH": os.environ.get("PATH", "/usr/bin:/bin")}
                cap = 120
            print(f"{selected['pipeline']} {arm} starting at {parent.now()}", flush=True)
            usage = parent.one_run(folder, command, environment, cap, group)
            print(f"{selected['pipeline']} {arm}: {usage['status']} after "
                  f"{usage['wall_seconds']:.3f}s", flush=True)
        parent.write_json(pair / "system_after.json", core.system_counters())
        parent.write_json(pair / "calibration_after.json", core.calibrate())
    parent.write_json(output / "allocation_end.json", {
        "recorded_utc": parent.now(),
        "instance_id": host["instance_id"],
        "cgroup_after": parent.cgroup_snapshot(group),
        "visible_affinity": sorted(os.sched_getaffinity(0)),
    })
    return 0


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--preflight", action="store_true")
    parser.add_argument("--binary", type=Path)
    parser.add_argument("--rho-binary", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.preflight:
        hashes, point, _, candidates = frozen()
        print(json.dumps({"hashes": hashes, "public_point": point,
                          "candidate_order": candidates}, sort_keys=True))
        return 0
    if not args.binary or not args.rho_binary or not args.output:
        parser.error("--binary, --rho-binary and --output are required")
    return run(args.binary.resolve(), args.rho_binary.resolve(), args.output.resolve())


if __name__ == "__main__":
    raise SystemExit(main())
