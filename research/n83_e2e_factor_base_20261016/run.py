#!/usr/bin/env python3
"""Frozen current-source one-target cold K-grid screen for N83."""

import argparse
import importlib.util
import json
import os
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SOURCE = ROOT / "examples/koblitz_orbit_dlp_fast.rs"
SOURCE_SHA256 = "8a2fe55dd90dc33892c197d27260182bf783fc8212fb485d8297daf2bbba8e12"
PRIOR = ROOT / "research/n83_root11_dynamic_20261016/run.py"
spec = importlib.util.spec_from_file_location("n83_dynamic_run_e2e", PRIOR)
dynamic = importlib.util.module_from_spec(spec)
spec.loader.exec_module(dynamic)
dynamic.SOURCE_SHA256 = SOURCE_SHA256
dynamic.parent.EXPECTED[SOURCE] = SOURCE_SHA256
parent = dynamic.parent
cold = dynamic.cold

K_GRID = tuple(sorted(set(range(896, 1313, 32)) |
                      {1096, 1100, 1101, 1102, 1104, 1108}))
CONTROL_K = 1152
ARMS = (("k1152_a", CONTROL_K),
        *((f"k{k}", k) for k in K_GRID if k != CONTROL_K),
        ("k1152_b", CONTROL_K))
RHO_SEED = "2026101601"
IC_CAP_SECONDS = 360
RHO_CAP_SECONDS = 120


def frozen_inputs():
    hashes = dynamic.frozen_inputs()
    if parent.digest(SOURCE) != SOURCE_SHA256:
        raise ValueError("current IC source differs")
    if len(cold.TARGETS.read_text().splitlines()) != 1:
        raise ValueError("one frozen public target required")
    return hashes


def run(binary, rho_binary, output):
    hashes = frozen_inputs()
    if output.exists():
        raise ValueError("output directory already exists")
    output.mkdir(parents=True)
    group = parent.cgroup_root()
    host = parent.host_receipt(group)
    host["schema"] = "n83.e2e-factor-base-host/v1"
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
        "schema": "n83.e2e-factor-base-source/v1", "hashes": hashes,
        "source_sha256": SOURCE_SHA256,
        "binary_sha256": parent.digest(binary),
        "rho_binary_sha256": parent.digest(rho_binary),
        "runner_sha256": parent.digest(Path(__file__)),
        "target_sha256": parent.digest(cold.TARGETS),
        "arm_order": [arm for arm, _ in ARMS], "k_grid": list(K_GRID),
        "rho_seed": RHO_SEED,
    })
    for arm, k in ARMS:
        folder = output / "runs" / arm
        command = [str(binary), f"construct:83:1:{k}", str(cold.TARGETS),
                   "29", str(folder / "targets_result.jsonl")]
        print(f"starting {arm} at {parent.now()}", flush=True)
        usage = parent.one_run(folder, command,
                               dynamic.environment(folder, "root11_dynamic"),
                               IC_CAP_SECONDS, group)
        print(f"{arm}: {usage['status']} after {usage['wall_seconds']:.3f}s",
              flush=True)
    folder = output / "runs" / "rho"
    command = [str(rho_binary), str(cold.TARGETS), RHO_SEED,
               str(cold.RHO_CONFIG), str(folder / "receipt.json")]
    usage = parent.one_run(folder, command,
                           {"PATH": os.environ.get("PATH", "/usr/bin:/bin")},
                           RHO_CAP_SECONDS, group)
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
