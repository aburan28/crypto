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
TARGETS = ROOT / "research/n83_factor_base_refine_20261014/reserve_targets.jsonl"
ANSWERS = ROOT / "research/n83_factor_base_refine_20261014/reserve_answers.jsonl"
SELECTION = HERE / "selection.json"
PANEL_SELECTION = ROOT / "research/n83_factor_base_refine_results_20261014/selected_candidate.json"
TARGETS_SHA256 = "2ccaac2baf367598e0d4e724b500ba66c5b4a711f0be8207565ce3acd6b1d0b4"
ANSWERS_SHA256 = "5b170d9bc4d34fdd21c032990ceee3a62af335c9af7fcaae2596335b8699dd8c"
SOURCE_SHA256 = "5e9fd826d22816b429a82e658c517d9387916aabe7660f708e4c29901cce647e"

spec = importlib.util.spec_from_file_location("n83_final_parent", PARENT)
parent = importlib.util.module_from_spec(spec)
spec.loader.exec_module(parent)


def frozen_inputs():
    hashes = parent.check_inputs()
    if hashes["examples/koblitz_orbit_dlp_fast.rs"] != SOURCE_SHA256:
        raise ValueError("producer source changed")
    if parent.digest(TARGETS) != TARGETS_SHA256 or parent.digest(ANSWERS) != ANSWERS_SHA256:
        raise ValueError("fresh target or answer fixture changed")
    hashes[str(TARGETS.relative_to(ROOT))] = TARGETS_SHA256
    hashes[str(ANSWERS.relative_to(ROOT))] = ANSWERS_SHA256
    targets = [json.loads(line) for line in TARGETS.read_text().splitlines()]
    if len(targets) != 5 or len({tuple(p) for p in targets}) != 5:
        raise ValueError("fresh target count or uniqueness changed")
    return hashes, targets


def selected_config():
    data = json.loads(SELECTION.read_text())
    if (data["schema"] != "n83.tier1-refine-reserve-selection/v1"
            or data["source_sha256"] != SOURCE_SHA256
            or data["variant"] not in ("root11", "root21", "root22")
            or not isinstance(data["k"], int) or data["k"] < 1
            or data["rss_cap_bytes"] != parent.RSS_CAP_BYTES
            or not isinstance(data["ic_child_cap_seconds"], int)
            or data["ic_child_cap_seconds"] <= 0
            or data["ic_child_cap_seconds"] > 1500
            or len(data["selection_analysis_sha256"]) != 64
            or any(c not in "0123456789abcdef"
                   for c in data["selection_analysis_sha256"])
            or not data["candidate_id"].startswith("IC1N83")):
        raise ValueError("selection manifest incomplete or inconsistent")
    panel = json.loads(PANEL_SELECTION.read_text())
    manifest_path = PANEL_SELECTION.parent / f"candidate_K{data['k']}_{data['variant']}.json"
    manifest = json.loads(manifest_path.read_text())
    if (panel["schema"] != "n83.factor-base-refine-selection/v1"
            or panel["status"] != "PANEL_SELECTED"
            or panel["selected_arm"] != f"k{data['k']}_{data['variant']}"
            or parent.digest(PANEL_SELECTION) != data["panel_selection_sha256"]
            or panel["analysis_sha256"] != data["selection_analysis_sha256"]
            or panel["source_sha256"] != data["source_sha256"]
            or panel["candidate_id"] != data["candidate_id"]
            or panel["candidate_manifest_sha256"] != data["candidate_manifest_sha256"]
            or panel["source_base_hash_blake3"] != data["source_base_hash_blake3"]
            or parent.digest(manifest_path) != data["candidate_manifest_sha256"]
            or manifest["candidate_id"] != data["candidate_id"]
            or manifest["implementation"]["source_sha256"]["ic"] != SOURCE_SHA256
            or manifest["factor_base"]["actual_usable_points"] != 166 * data["k"]
            or manifest["factor_base"]["effective_columns"] != data["k"]
            or manifest["factor_base"]["source_base_hash_blake3"] !=
            data["source_base_hash_blake3"]):
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
    return {"recorded_utc": parent.now(),
            "proc_stat": parent.read_text(Path("/proc/stat")),
            "proc_interrupts": parent.read_text(Path("/proc/interrupts"))}


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


def ic_environment(folder, variant):
    return {
        "PATH": os.environ.get("PATH", "/usr/bin:/bin"),
        "KIC_PACKED_ROOT_TABLE": "1",
        "KIC_PACKED_ROOT_TABLE_DENSE": "0",
        "KIC_PACKED_ROOT_TABLE_LOAD75": "0",
        "KIC_STATE28": "0", "KIC_STATE26": "0",
        "KIC_ROOT_ONLY22": "1" if variant == "root22" else "0",
        "KIC_ROOT_ONLY21": "1" if variant == "root21" else "0",
        "KIC_ROOT_ONLY11": "1" if variant == "root11" else "0",
        "KIC_INDEXED_ROOT5": "1", "KIC_INDEXED_ROOT4": "0",
        "KIC_ROOT_BLOOM": "0",
        "KIC_SYMMETRIC_S3": "1",
        "KIC_INDEXED5_LOAD75": "1" if variant == "root22" else "0",
        "KIC_INDEXED5_FIXED30": "0" if variant == "root22" else "1",
        "KIC_INDEXED4_FIXED31": "0",
        "KIC_SPARSE_CANONICAL": "0", "KIC_RANK_THREADS": "4",
        "KIC_TARGET_THREADS": "1", "KIC_DUMP_BASE": str(folder / "base_dump.jsonl"),
    }


def run(binary, rho_binary, output):
    hashes, targets = frozen_inputs()
    selected = selected_config()
    if output.exists():
        raise ValueError("output directory already exists")
    output.mkdir(parents=True)
    group = parent.cgroup_root()
    host = parent.host_receipt(group)
    host["schema"] = "n83.tier1-refine-reserve-host/v1"
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
        "schema": "n83.tier1-refine-reserve-source/v1", "hashes": hashes,
        "selection_sha256": parent.digest(SELECTION),
        "selection": selected,
        "ic_binary_sha256": parent.digest(binary),
        "rho_binary_sha256": parent.digest(rho_binary),
        "runner_sha256": parent.digest(Path(__file__)),
    })
    memory_max = (host["cgroup_before"]["memory.max"] or "").strip()
    memory_limit_too_small = (memory_max.isdigit()
                              and int(memory_max) < parent.RSS_CAP_BYTES)
    if (host["mem_total_kib"] is None or host["mem_total_kib"] < 15_000_000
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
                environment = ic_environment(folder, selected["variant"])
                cap = selected["ic_child_cap_seconds"]
            else:
                command = [str(rho_binary), str(point_file),
                           str(2026101401 + index), str(parent.RHO_CONFIG),
                           str(folder / "receipt.json")]
                environment = {"PATH": os.environ.get("PATH", "/usr/bin:/bin")}
                cap = 120
            print(f"pair {index+1} {arm} starting at {parent.now()}", flush=True)
            usage = parent.one_run(folder, command, environment, cap, group)
            print(f"pair {index+1} {arm}: {usage['status']} after "
                  f"{usage['wall_seconds']:.3f}s", flush=True)
        parent.write_json(pair / "system_after.json", system_counters())
        parent.write_json(pair / "calibration_after.json", calibrate())
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
