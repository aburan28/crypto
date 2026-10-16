#!/usr/bin/env python3
"""Bounded one-point online K-range screen for exact N83 root indexes."""

import argparse
import importlib.util
import json
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
DYNAMIC_RUNNER = ROOT / "research/n83_root11_dynamic_20261016/run.py"
SOURCE_SHA256 = "3ee6bc11399b382e72789b152043e5a67abc5b783c4e998ad8c534d60baebc03"
ARMS = (
    ("k2543root22_a", 2543, "root22_dynamic", 900),
    ("k1024root22", 1024, "root22_dynamic", 900),
    ("k1152root22", 1152, "root22_dynamic", 900),
    ("k1152root11", 1152, "root11_dynamic", 900),
    ("k1536root22", 1536, "root22_dynamic", 900),
    ("k1536root11", 1536, "root11_dynamic", 900),
    ("k2048root22", 2048, "root22_dynamic", 900),
    ("k2543root11", 2543, "root11_dynamic", 900),
    ("k3072root11", 3072, "root11_dynamic", 900),
    ("k3450root11fixed", 3450, "root11_fixed", 900),
    ("k3584root11fixed", 3584, "root11_fixed", 900),
    ("k2543root22_b", 2543, "root22_dynamic", 900),
)

spec = importlib.util.spec_from_file_location("n83_online_dynamic_runner", DYNAMIC_RUNNER)
dynamic = importlib.util.module_from_spec(spec)
spec.loader.exec_module(dynamic)
cold = dynamic.cold
parent = dynamic.parent
dynamic.SOURCE_SHA256 = SOURCE_SHA256
parent.EXPECTED[parent.SOURCE] = SOURCE_SHA256


def frozen_inputs():
    hashes = dynamic.frozen_inputs()
    if parent.digest(dynamic.SOURCE) != SOURCE_SHA256:
        raise ValueError("executed merged-main N83 source differs")
    return hashes


def run(binary, rho_binary, output):
    hashes = frozen_inputs()
    if output.exists():
        raise ValueError("output directory already exists")
    output.mkdir(parents=True)
    group = parent.cgroup_root()
    host = parent.host_receipt(group)
    host["schema"] = "n83.online-range-host/v1"
    parent.write_json(output / "host.json", host)
    parent.write_json(output / "source.json", {
        "schema": "n83.online-range-source/v1", "hashes": hashes,
        "source_sha256": SOURCE_SHA256,
        "binary_sha256": parent.digest(binary),
        "rho_binary_sha256": parent.digest(rho_binary),
        "runner_sha256": parent.digest(Path(__file__)),
        "arm_order": [arm for arm, _, _, _ in ARMS],
        "rho_seed": cold.RHO_SEED,
    })
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
    for arm, k, layout, cap in ARMS:
        folder = output / "runs" / arm
        command = [str(binary), f"construct:83:1:{k}", str(cold.TARGETS),
                   "29", str(folder / "targets_result.jsonl")]
        print(f"starting {arm} at {parent.now()}", flush=True)
        usage = parent.one_run(folder, command, dynamic.environment(folder, layout),
                               cap, group)
        print(f"{arm}: {usage['status']} after {usage['wall_seconds']:.3f}s",
              flush=True)
    folder = output / "runs" / "rho"
    command = [str(rho_binary), str(cold.TARGETS), cold.RHO_SEED,
               str(cold.RHO_CONFIG), str(folder / "receipt.json")]
    print(f"starting same-point rho at {parent.now()}", flush=True)
    usage = parent.one_run(folder, command,
                           {"PATH": os.environ.get("PATH", "/usr/bin:/bin")},
                           cold.RHO_CAP_SECONDS, group)
    print(f"rho: {usage['status']} after {usage['wall_seconds']:.3f}s", flush=True)
    return 0


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--preflight", action="store_true")
    parser.add_argument("--binary", type=Path)
    parser.add_argument("--rho-binary", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.preflight:
        print(json.dumps({"hashes": frozen_inputs(), "arms": ARMS}, sort_keys=True))
        return 0
    if not args.binary or not args.rho_binary or not args.output:
        parser.error("--binary, --rho-binary and --output are required")
    return run(args.binary.resolve(), args.rho_binary.resolve(),
               args.output.resolve())


if __name__ == "__main__":
    raise SystemExit(main())
