#!/usr/bin/env python3
"""Bounded same-point N83 one-root dynamic-index comparison."""

import argparse
import importlib.util
import json
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
COLD_RUNNER = ROOT / "research/n83_cold_factor_base_20261016/run.py"
SOURCE = ROOT / "examples/koblitz_orbit_dlp_fast.rs"
SOURCE_SHA256 = "d094345474605075566d646ddecef8f694cdd644b0f9190a6b1ef1551fd06888"
ARMS = (
    ("k1271fixed_a", 1271, "root11_fixed", 900),
    ("k1271dynamic_a", 1271, "root11_dynamic", 900),
    ("k1271root22", 1271, "root22_dynamic", 900),
    ("k1024dynamic", 1024, "root11_dynamic", 900),
    ("k1536dynamic", 1536, "root11_dynamic", 900),
    ("k1271dynamic_b", 1271, "root11_dynamic", 900),
    ("k1271fixed_b", 1271, "root11_fixed", 900),
)

spec = importlib.util.spec_from_file_location("n83_cold_runner", COLD_RUNNER)
cold = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cold)
parent = cold.parent
parent.EXPECTED[parent.SOURCE] = SOURCE_SHA256


def frozen_inputs():
    hashes = cold.frozen_inputs()
    if parent.digest(SOURCE) != SOURCE_SHA256:
        raise ValueError("dynamic one-root source differs")
    return hashes


def environment(folder, layout):
    flags = cold.environment(folder, "root22", "off")
    if layout == "root22_dynamic":
        return flags
    if layout not in ("root11_fixed", "root11_dynamic"):
        raise ValueError("unsupported one-root layout")
    flags["KIC_ROOT_ONLY22"] = "0"
    flags["KIC_ROOT_ONLY11"] = "1"
    if layout == "root11_fixed":
        flags["KIC_INDEXED5_LOAD75"] = "0"
        flags["KIC_INDEXED5_FIXED30"] = "1"
    else:
        flags["KIC_INDEXED5_LOAD75"] = "1"
        flags["KIC_INDEXED5_FIXED30"] = "0"
    return flags


def run(binary, rho_binary, output):
    hashes = frozen_inputs()
    if output.exists():
        raise ValueError("output directory already exists")
    output.mkdir(parents=True)
    group = parent.cgroup_root()
    host = parent.host_receipt(group)
    host["schema"] = "n83.root11-dynamic-host/v1"
    parent.write_json(output / "host.json", host)
    parent.write_json(output / "source.json", {
        "schema": "n83.root11-dynamic-source/v1", "hashes": hashes,
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
        command = [str(binary), f"construct:83:1:{k}", str(cold.TARGETS), "29",
                   str(folder / "targets_result.jsonl")]
        print(f"starting {arm} at {parent.now()}", flush=True)
        usage = parent.one_run(folder, command, environment(folder, layout), cap, group)
        print(f"{arm}: {usage['status']} after {usage['wall_seconds']:.3f}s", flush=True)
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
    return run(args.binary.resolve(), args.rho_binary.resolve(), args.output.resolve())


if __name__ == "__main__":
    raise SystemExit(main())
