#!/usr/bin/env python3
"""Five alternating same-point N83 IC/rho pairs on one cloud allocation."""

import argparse
import hashlib
import importlib.util
import json
import os
import resource
import statistics
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
PARENT = ROOT / "research/n83_root4_20261014/run.py"
FINAL_RUNNER = ROOT / "research/n83_final_tournament_20261015/run.py"
TARGETS = HERE / "targets.jsonl"
ANSWERS = HERE / "answers.jsonl"
SELECTION = HERE / "selection.json"
TARGETS_SHA256 = "34c530851eb135b7e348ac976caaae37c5883c52b52962ad496028c058b1dcda"
ANSWERS_SHA256 = "8705b909479f4ee1eb8a19fe327c0622cd890e0d8443e08bdbdf04a4c6b99548"
SOURCE_SHA256 = "f2745227d95e36ef45c73c5104a23ada6cebb0794405368581d2be17044845dc"
RSS_CAP_BYTES = 29 * 1024**3 // 2

spec = importlib.util.spec_from_file_location("n83_final_parent", PARENT)
parent = importlib.util.module_from_spec(spec)
spec.loader.exec_module(parent)
parent.EXPECTED[parent.SOURCE] = SOURCE_SHA256
parent.RSS_CAP_BYTES = RSS_CAP_BYTES
fspec = importlib.util.spec_from_file_location("n83_final_tournament_runner", FINAL_RUNNER)
final_runner = importlib.util.module_from_spec(fspec)
fspec.loader.exec_module(final_runner)


def frozen_inputs():
    hashes = parent.check_inputs()
    if hashes["examples/koblitz_orbit_dlp_fast.rs"] != SOURCE_SHA256:
        raise ValueError("producer source changed")
    if parent.digest(TARGETS) != TARGETS_SHA256 or parent.digest(ANSWERS) != ANSWERS_SHA256:
        raise ValueError("fresh target or answer fixture changed")
    hashes[str(TARGETS.relative_to(ROOT))] = TARGETS_SHA256
    hashes[str(ANSWERS.relative_to(ROOT))] = ANSWERS_SHA256
    hashes[str(SELECTION.relative_to(ROOT))] = parent.digest(SELECTION)
    hashes[str(final_runner.ADMISSION.relative_to(ROOT))] = parent.digest(final_runner.ADMISSION)
    targets = [json.loads(line) for line in TARGETS.read_text().splitlines()]
    if len(targets) != 5 or len({tuple(p) for p in targets}) != 5:
        raise ValueError("fresh target count or uniqueness changed")
    return hashes, targets


def selected_config():
    data = json.loads(SELECTION.read_text())
    admission, _ = final_runner.frozen()
    candidates = admission["candidate_order"]
    selected = [row for row in candidates if row["pipeline"] == data["pipeline"]]
    hashes = (data["tournament_artifact_sha256"],
              data["tournament_analysis_sha256"])
    if (data["schema"] != "n83.final-tier1-selection/v1"
            or data["source_sha256"] != SOURCE_SHA256
            or len(selected) != 1
            or data["variant"] not in ("root4", "root11", "root21", "root22")
            or data["bloom_profile"] not in ("off", "standard", "wide3")
            or not isinstance(data["k"], int) or data["k"] < 1
            or data["rss_cap_bytes"] != RSS_CAP_BYTES
            or data["ic_child_cap_seconds"] != 1500
            or data["rho_child_cap_seconds"] != 120
            or not isinstance(data["tournament_artifact_id"], int)
            or data["tournament_artifact_id"] <= 0
            or data["admission_sha256"] != parent.digest(final_runner.ADMISSION)
            or any(len(value) != 64 or any(c not in "0123456789abcdef" for c in value)
                   for value in hashes)
            or not data["candidate_id"].startswith("IC1N83")):
        raise ValueError("selection manifest incomplete or inconsistent")
    row = selected[0]
    if (data["candidate_id"] != row["candidate_id"]
            or data["candidate_manifest_sha256"] != row["candidate_manifest_sha256"]
            or data["k"] != row["k"] or data["variant"] != row["variant"]
            or data["bloom_profile"] != row["bloom_profile"]
            or data["source_base_hash_blake3"] != row["base_hash"]
            or data["root_table_entries"] != row["root_table_entries"]
            or row["actual_usable_points"] != 166 * data["k"]):
        raise ValueError("selection differs from sealed panel candidate")
    return data


def calibrate():
    durations = []
    digests = []
    for repetition in range(3):
        state = (b"n83-tier1-calibration-v1" + repetition.to_bytes(4, "little"))
        before_cpu = resource.getrusage(resource.RUSAGE_SELF)
        started = time.monotonic_ns()
        for _ in range(1_000_000):
            state = hashlib.sha256(state).digest()
        durations.append((time.monotonic_ns() - started) / 1_000_000)
        after_cpu = resource.getrusage(resource.RUSAGE_SELF)
        digests.append({"sha256": state.hex(),
                        "user_cpu_seconds": after_cpu.ru_utime - before_cpu.ru_utime,
                        "system_cpu_seconds": after_cpu.ru_stime - before_cpu.ru_stime})
    return {"sample_ms": durations, "median_ms": statistics.median(durations),
            "checksums_and_cpu": digests, "recorded_utc": parent.now()}


def system_counters():
    visible = sorted(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else []
    return {"recorded_utc": parent.now(),
            "proc_stat": parent.read_text(Path("/proc/stat")),
            "proc_interrupts": parent.read_text(Path("/proc/interrupts")),
            "cpu_frequency_khz": {
                str(cpu): parent.read_text(Path(
                    f"/sys/devices/system/cpu/cpu{cpu}/cpufreq/scaling_cur_freq"))
                for cpu in visible}}


def select_execution_cpus(visible):
    topology = {}
    for cpu in visible:
        base = Path(f"/sys/devices/system/cpu/cpu{cpu}")
        node = next((entry.name for entry in base.glob("node*") if entry.name[4:].isdigit()), None)
        topology[cpu] = {
            "core_id": parent.read_text(base / "topology/core_id"),
            "package_id": parent.read_text(base / "topology/physical_package_id"),
            "thread_siblings_list": parent.read_text(base / "topology/thread_siblings_list"),
            "numa_node": node,
            "scaling_governor": parent.read_text(base / "cpufreq/scaling_governor"),
        }
    chosen = []
    seen_cores = set()
    for cpu in visible:
        key = (topology[cpu]["package_id"], topology[cpu]["core_id"])
        if key not in seen_cores:
            chosen.append(cpu)
            seen_cores.add(key)
        if len(chosen) == 4:
            break
    for cpu in visible:
        if len(chosen) == 4:
            break
        if cpu not in chosen:
            chosen.append(cpu)
    return chosen, {str(cpu): topology[cpu] for cpu in chosen}


def ic_environment(folder, selected):
    return final_runner.environment(folder, selected)


def run(binary, rho_binary, output):
    hashes, targets = frozen_inputs()
    selected = selected_config()
    if output.exists():
        raise ValueError("output directory already exists")
    output.mkdir(parents=True)
    group = parent.cgroup_root()
    host = parent.host_receipt(group)
    host["schema"] = "n83.final-tier1-host/v1"
    host["region_or_zone"] = (os.environ.get("RUNNER_REGION") or
                              os.environ.get("AZURE_REGION"))
    host["region_or_zone_status"] = (
        "EXPOSED" if host["region_or_zone"] else "NOT_EXPOSED_BY_RUNNER")
    execution_cpus, topology = select_execution_cpus(host["visible_affinity"] or [])
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
        "schema": "n83.final-tier1-source/v1", "hashes": hashes,
        "selection_sha256": parent.digest(SELECTION),
        "selection": selected,
        "ic_binary_sha256": parent.digest(binary),
        "rho_binary_sha256": parent.digest(rho_binary),
        "runner_sha256": parent.digest(Path(__file__)),
    })
    memory_max = (host["cgroup_before"]["memory.max"] or "").strip()
    memory_limit_too_small = (memory_max.isdigit()
                              and int(memory_max) < parent.RSS_CAP_BYTES)
    if (host["mem_total_kib"] is None or
            host["mem_total_kib"] * 1024 < RSS_CAP_BYTES + 256 * 1024**2
            or host["disk_free_bytes"] < 2 * 1024**3
            or len(execution_cpus) < 4
            or memory_limit_too_small
            or host["affinity_set_status"] != "PASS"):
        parent.write_json(output / "preflight.json", {
            "status": "REJECTED_HOST_CAPACITY", "host": host})
        return 2
    parent.write_json(output / "preflight.json", {
        "status": "PASS_CAPACITY_PREFLIGHT", "rss_cap_bytes": parent.RSS_CAP_BYTES})
    inputs = output / "inputs"
    inputs.mkdir()
    for index, target in enumerate(targets):
        point_file = inputs / f"point_{index+1}.jsonl"
        point_file.write_text(json.dumps(target, separators=(",", ":")) + "\n")
        pair = output / "pairs" / f"p{index+1}"
        pair.mkdir(parents=True)
        parent.write_json(pair / "calibration_before.json", calibrate())
        parent.write_json(pair / "system_before.json", system_counters())
        order = ("ic", "rho") if index % 2 == 0 else ("rho", "ic")
        parent.write_json(pair / "order.json", {
            "target_index": index, "order": order, "point_sha256": parent.digest(point_file),
            "started_utc": parent.now()})
        for arm in order:
            folder = pair / arm
            if arm == "ic":
                command = [str(binary), f"construct:83:1:{selected['k']}",
                           str(point_file), "29", str(folder / "targets_result.jsonl")]
                environment = ic_environment(folder, selected)
                cap = selected["ic_child_cap_seconds"]
            else:
                command = [str(rho_binary), str(point_file),
                           str(2026101401 + index), str(parent.RHO_CONFIG),
                           str(folder / "receipt.json")]
                environment = {"PATH": os.environ.get("PATH", "/usr/bin:/bin")}
                cap = selected["rho_child_cap_seconds"]
            print(f"pair {index+1} {arm} starting at {parent.now()}", flush=True)
            usage = parent.one_run(folder, command, environment, cap, group)
            print(f"pair {index+1} {arm}: {usage['status']} after "
                  f"{usage['wall_seconds']:.3f}s", flush=True)
        parent.write_json(pair / "system_after.json", system_counters())
        parent.write_json(pair / "calibration_after.json", calibrate())
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
        hashes, targets = frozen_inputs()
        print(json.dumps({"hashes": hashes, "target_count": len(targets),
                          "selection": selected_config(),
                          "host": parent.host_receipt(parent.cgroup_root())},
                         sort_keys=True))
        return 0
    if not args.binary or not args.rho_binary or not args.output:
        parser.error("--binary, --rho-binary and --output are required")
    return run(args.binary.resolve(), args.rho_binary.resolve(), args.output.resolve())


if __name__ == "__main__":
    raise SystemExit(main())
