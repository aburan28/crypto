#!/usr/bin/env python3
"""Compare root11 and root22 at the N83 K1107 and K1112 screen leaders."""

import argparse
import importlib.util
import json
import os
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
ARMS = (
    ("k1107root11_a", 1107, "root11_dynamic"),
    ("k1107root22_a", 1107, "root22_dynamic"),
    ("k1112root11_a", 1112, "root11_dynamic"),
    ("k1112root22_a", 1112, "root22_dynamic"),
    ("k1112root22_b", 1112, "root22_dynamic"),
    ("k1112root11_b", 1112, "root11_dynamic"),
    ("k1107root22_b", 1107, "root22_dynamic"),
    ("k1107root11_b", 1107, "root11_dynamic"),
)
IC_CAP_SECONDS = 480
RHO_CAP_SECONDS = 120


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    value = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(value)
    return value


finegrid = module("n83_root22_neighbors_prior",
                  ROOT / "research/n83_e2e_finegrid_20261016/run.py")
parent = finegrid.parent
holdout = finegrid.holdout
dynamic = finegrid.dynamic
cold = finegrid.cold
SOURCE_SHA256 = finegrid.SOURCE_SHA256
TARGET_SHA256 = finegrid.TARGET_SHA256
RHO_SEED = finegrid.RHO_SEED


def frozen_inputs():
    return finegrid.frozen_inputs()


def run(binary, rho_binary, output):
    hashes = frozen_inputs()
    if output.exists():
        raise ValueError("output directory already exists")
    output.mkdir(parents=True)
    group = parent.cgroup_root()
    host = parent.host_receipt(group)
    host["schema"] = "n83.e2e-root22-neighbors-host/v1"
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
        "schema": "n83.e2e-root22-neighbors-source/v1", "hashes": hashes,
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
