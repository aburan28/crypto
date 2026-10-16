#!/usr/bin/env python3
"""Bounded same-base N83 root11 canonical-selector screen."""

import argparse
import importlib.util
import json
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PARENT = ROOT / "research/n83_root4_20261014/run.py"
POINTS = ROOT / "research/n83_factor_base_refine_20261014/targets.jsonl"
POINTS_SHA256 = "35e34151dcc1ad6048939560adf28fc7c1b1017833fe005a656c3f4e2bc07de8"
SOURCE_SHA256 = "c9824acebda0360416c3c507594bb07c17ca798c32cbf8f84bba95ca2c44d5b1"
PARENT_SOURCE_SHA256 = "5e9fd826d22816b429a82e658c517d9387916aabe7660f708e4c29901cce647e"
ARMS = (("canonical_scan", False, 1200), ("canonical_sparse", True, 1200))

spec = importlib.util.spec_from_file_location("n83_sparse_canonical_parent", PARENT)
parent = importlib.util.module_from_spec(spec)
spec.loader.exec_module(parent)
if parent.EXPECTED[parent.SOURCE] != PARENT_SOURCE_SHA256:
    raise ValueError("parent source pin changed")
parent.EXPECTED[parent.SOURCE] = SOURCE_SHA256


def frozen_inputs():
    hashes = parent.check_inputs()
    if parent.digest(POINTS) != POINTS_SHA256:
        raise ValueError("frozen 1000-point panel changed")
    hashes[str(POINTS.relative_to(ROOT))] = POINTS_SHA256
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
        "schema": "n83.sparse-canonical-source/v1", "hashes": hashes,
        "binary_sha256": parent.digest(binary),
        "screen_runner_sha256": parent.digest(Path(__file__)),
    })
    if (host["mem_total_kib"] is None or host["mem_total_kib"] < 15_000_000
            or host["disk_free_bytes"] < 2 * 1024**3
            or len(host["visible_affinity"] or []) < 4):
        parent.write_json(output / "preflight.json", {
            "status": "REJECTED_HOST_CAPACITY", "host": host})
        return 2
    parent.write_json(output / "preflight.json", {
        "status": "PASS_CAPACITY_PREFLIGHT", "rss_cap_bytes": parent.RSS_CAP_BYTES})
    for arm, sparse, cap in ARMS:
        folder = output / "runs" / arm
        env = {
            "PATH": os.environ.get("PATH", "/usr/bin:/bin"),
            "KIC_PACKED_ROOT_TABLE": "1", "KIC_PACKED_ROOT_TABLE_DENSE": "0",
            "KIC_PACKED_ROOT_TABLE_LOAD75": "0", "KIC_STATE28": "0",
            "KIC_STATE26": "0", "KIC_ROOT_ONLY22": "0",
            "KIC_ROOT_ONLY21": "0", "KIC_ROOT_ONLY11": "1",
            "KIC_INDEXED_ROOT5": "1", "KIC_INDEXED_ROOT4": "0",
            "KIC_ROOT_BLOOM": "0",
            "KIC_SYMMETRIC_S3": "1", "KIC_INDEXED5_LOAD75": "0",
            "KIC_INDEXED5_FIXED30": "1", "KIC_INDEXED4_FIXED31": "0",
            "KIC_SPARSE_CANONICAL": "1" if sparse else "0",
            "KIC_RANK_THREADS": "4",
            "KIC_TARGET_THREADS": "1", "KIC_DUMP_BASE": str(folder / "base_dump.jsonl"),
        }
        command = [str(binary), "construct:83:1:3300", str(POINTS), "29",
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
        print(json.dumps({"hashes": frozen_inputs(), "arms": ARMS,
                          "host": parent.host_receipt(parent.cgroup_root())},
                         sort_keys=True))
        return 0
    if not args.binary or not args.output:
        parser.error("--binary and --output are required")
    return run(args.binary.resolve(), args.output.resolve())


if __name__ == "__main__":
    raise SystemExit(main())
