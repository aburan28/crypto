#!/usr/bin/env python3
"""Screen exact root11/root22 index layouts at the N83 cold frontier."""

import argparse
import importlib.util
import json
import os
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    value = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(value)
    return value


holdout = module("n83_crosslayout_holdout_run",
                 ROOT / "research/n83_e2e_holdout_20261016/run.py")
screen = holdout.screen
dynamic = holdout.dynamic
cold = holdout.cold
parent = holdout.parent
SOURCE_SHA256 = screen.SOURCE_SHA256
TARGET_SHA256 = "2c068d2f659588edcdbaa685fc73bb146e2bfee3be123d377d537658f8021c6d"
RHO_SEED = cold.RHO_SEED
IC_CAP_SECONDS = 480
RHO_CAP_SECONDS = 120
ARMS = (
    ("k1104root11_a", 1104, "root11_dynamic"),
    ("k1104root22_a", 1104, "root22_dynamic"),
    ("k1120root22", 1120, "root22_dynamic"),
    ("k1120root11", 1120, "root11_dynamic"),
    ("k1104root22_b", 1104, "root22_dynamic"),
    ("k1104root11_b", 1104, "root11_dynamic"),
)


def frozen_inputs():
    hashes = screen.frozen_inputs()
    if parent.digest(cold.TARGETS) != TARGET_SHA256:
        raise ValueError("exposed screen target differs")
    return hashes


def run(binary, rho_binary, output):
    hashes = frozen_inputs()
    if output.exists():
        raise ValueError("output directory already exists")
    output.mkdir(parents=True)
    group = parent.cgroup_root()
    host = parent.host_receipt(group)
    host["schema"] = "n83.e2e-crosslayout-host/v1"
    host["topology"] = holdout.topology_receipt(host["visible_affinity"] or [])
    host["cpu_telemetry_before"] = holdout.cpu_telemetry()
    parent.write_json(output / "host.json", host)
    memory_max = (host["cgroup_before"]["memory.max"] or "").strip()
    memory_too_small = memory_max.isdigit() and int(memory_max) < cold.RSS_CAP_BYTES
    if (host["mem_total_kib"] is None
            or host["mem_total_kib"] * 1024 < cold.RSS_CAP_BYTES + 256 * 1024**2
            or host["disk_free_bytes"] < 2 * 1024**3
            or len(host["visible_affinity"] or []) < 4 or memory_too_small):
        parent.write_json(output / "preflight.json", {
            "status": "REJECTED_HOST_CAPACITY", "host": host})
        return 2
    parent.write_json(output / "preflight.json", {
        "status": "PASS_CAPACITY_PREFLIGHT", "rss_cap_bytes": cold.RSS_CAP_BYTES})
    parent.write_json(output / "source.json", {
        "schema": "n83.e2e-crosslayout-source/v1", "hashes": hashes,
        "source_sha256": SOURCE_SHA256,
        "binary_sha256": parent.digest(binary),
        "rho_binary_sha256": parent.digest(rho_binary),
        "runner_sha256": parent.digest(Path(__file__)),
        "target_sha256": TARGET_SHA256,
        "arm_order": [name for name, _, _ in ARMS],
        "arm_layouts": {name: layout for name, _, layout in ARMS},
        "rho_seed": RHO_SEED,
    })
    for name, k, layout in ARMS:
        folder = output / "runs" / name
        command = [str(binary), f"construct:83:1:{k}", str(cold.TARGETS),
                   "29", str(folder / "targets_result.jsonl")]
        print(f"starting {name} at {parent.now()}", flush=True)
        usage = holdout.run_with_telemetry(
            folder, command, dynamic.environment(folder, layout),
            IC_CAP_SECONDS, group)
        print(f"{name}: {usage['status']} after {usage['wall_seconds']:.3f}s",
              flush=True)
    folder = output / "runs" / "rho"
    command = [str(rho_binary), str(cold.TARGETS), RHO_SEED,
               str(cold.RHO_CONFIG), str(folder / "receipt.json")]
    print(f"starting same-point rho at {parent.now()}", flush=True)
    usage = holdout.run_with_telemetry(
        folder, command, {"PATH": os.environ.get("PATH", "/usr/bin:/bin")},
        RHO_CAP_SECONDS, group)
    print(f"rho: {usage['status']} after {usage['wall_seconds']:.3f}s", flush=True)
    parent.write_json(output / "host_telemetry_after.json",
                      holdout.cpu_telemetry())
    return 0


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--preflight", action="store_true")
    parser.add_argument("--binary", type=Path)
    parser.add_argument("--rho-binary", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.preflight:
        print(json.dumps({"hashes": frozen_inputs(), "arms": ARMS,
                          "target_sha256": TARGET_SHA256}, sort_keys=True))
        return 0
    if not args.binary or not args.rho_binary or not args.output:
        parser.error("--binary, --rho-binary and --output are required")
    return run(args.binary.resolve(), args.rho_binary.resolve(),
               args.output.resolve())


if __name__ == "__main__":
    raise SystemExit(main())
