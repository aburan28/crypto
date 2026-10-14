#!/usr/bin/env python3
"""Frozen secondary N83 factor-base panel on 100 archived public points."""

import argparse
import importlib.util
import json
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PARENT = ROOT / "research/n83_root4_20261014/run.py"
HOLDOUT = ROOT / "research/n83_factor_base_100_20261014/targets.jsonl"
HOLDOUT_SHA256 = "fa0044f7cad85185812ac4c79f4feabb300ddcdc1732b6a0890778d5b8e72af6"
ARMS = ((3240, "root21", 1300), (3400, "root4", 1500),
        (3300, "root11", 1200), (3500, "root4", 1500),
        (3220, "root21", 1300), (3550, "root4", 1500),
        (3114, "root22", 1100), (3596, "root4", 1500),
        (3200, "root11", 1100), (3400, "root5", 1300),
        (3190, "root21", 1300), (2940, "root22", 1000))

spec = importlib.util.spec_from_file_location("n83_root11_runner", PARENT)
parent = importlib.util.module_from_spec(spec)
spec.loader.exec_module(parent)


def frozen_inputs():
    hashes = parent.check_inputs()
    if parent.digest(HOLDOUT) != HOLDOUT_SHA256:
        raise ValueError("archived 100-point selection fixture changed")
    hashes[str(HOLDOUT.relative_to(ROOT))] = HOLDOUT_SHA256
    return hashes


def run(binary, output):
    hashes = frozen_inputs()
    if output.exists():
        raise ValueError("output directory already exists")
    output.mkdir(parents=True)
    group = parent.cgroup_root()
    host = parent.host_receipt(group)
    parent.write_json(output / "host.json", host)
    parent.write_json(output / "source.json", {
        "schema": "n83.factor-base-100-source/v1",
        "hashes": hashes,
        "binary_sha256": parent.digest(binary),
        "panel_runner_sha256": parent.digest(Path(__file__)),
    })
    if (host["mem_total_kib"] is None or host["mem_total_kib"] < 15_000_000
            or host["disk_free_bytes"] < 2 * 1024**3
            or len(host["visible_affinity"] or []) < 4):
        parent.write_json(output / "preflight.json", {
            "status": "REJECTED_HOST_CAPACITY", "host": host})
        return 2
    parent.write_json(output / "preflight.json", {
        "status": "PASS_CAPACITY_PREFLIGHT", "rss_cap_bytes": parent.RSS_CAP_BYTES})
    for k, variant, cap in ARMS:
        arm = f"k{k}_{variant}"
        folder = output / "runs" / arm
        env = {
            "PATH": os.environ.get("PATH", "/usr/bin:/bin"),
            "KIC_PACKED_ROOT_TABLE": "1",
            "KIC_PACKED_ROOT_TABLE_DENSE": "0",
            "KIC_PACKED_ROOT_TABLE_LOAD75": "0",
            "KIC_STATE28": "0", "KIC_STATE26": "0",
            "KIC_ROOT_ONLY22": "1" if variant == "root22" else "0",
            "KIC_ROOT_ONLY21": "1" if variant == "root21" else "0",
            "KIC_ROOT_ONLY11": "1" if variant in ("root11", "root4", "root5") else "0",
            "KIC_INDEXED_ROOT5": "0" if variant == "root4" else "1",
            "KIC_INDEXED_ROOT4": "1" if variant == "root4" else "0",
            "KIC_ROOT_BLOOM": "0",
            "KIC_SYMMETRIC_S3": "1",
            "KIC_INDEXED5_LOAD75": "1" if variant == "root22" else "0",
            "KIC_INDEXED5_FIXED30": "1" if variant in ("root11", "root21", "root5") else "0",
            "KIC_INDEXED4_FIXED31": "1" if variant == "root4" else "0",
            "KIC_SPARSE_CANONICAL": "0", "KIC_RANK_THREADS": "4",
            "KIC_TARGET_THREADS": "1", "KIC_DUMP_BASE": str(folder / "base_dump.jsonl"),
        }
        command = [str(binary), f"construct:83:1:{k}", str(HOLDOUT), "29",
                   str(folder / "targets_result.jsonl")]
        print(f"starting {arm} at {parent.now()}", flush=True)
        usage = parent.one_run(folder, command, env, cap, group)
        print(f"{arm}: {usage['status']} after {usage['wall_seconds']:.3f}s",
              flush=True)
    return 0


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--preflight", action="store_true")
    parser.add_argument("--binary", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.preflight:
        print(json.dumps({"hashes": frozen_inputs(),
                          "arms": ARMS, "host": parent.host_receipt(parent.cgroup_root())},
                         sort_keys=True))
        return 0
    if not args.binary or not args.output:
        parser.error("--binary and --output are required")
    return run(args.binary.resolve(), args.output.resolve())


if __name__ == "__main__":
    raise SystemExit(main())
