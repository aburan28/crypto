#!/usr/bin/env python3
"""Frozen N83 four-byte fixed30 table with 11-byte reconstructing states."""

import argparse
import importlib.util
import json
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "research/n83_root4_fixed30_20261014/run.py"
SOURCE_SHA256 = "5281bfc1c370255d8c03851fe85b4b35499ee55c066ba460b5de414b59d202a3"
ARMS = ((3300, "root4f30r11", 1200), (3300, "root11", 1200),
        (3300, "root4f30", 1200), (3350, "root4f30r11", 1200),
        (3400, "root4f30r11", 1200), (3450, "root4f30r11", 1200),
        (3500, "root4f30r11", 1200))

spec = importlib.util.spec_from_file_location("n83_fixed30_parent", BASE)
base = importlib.util.module_from_spec(spec)
spec.loader.exec_module(base)
base.parent.EXPECTED[base.parent.SOURCE] = SOURCE_SHA256


def frozen_inputs():
    return base.frozen_inputs()


def run(binary, output):
    hashes = frozen_inputs()
    if output.exists():
        raise ValueError("output directory already exists")
    output.mkdir(parents=True)
    group = base.parent.cgroup_root()
    host = base.parent.host_receipt(group)
    base.parent.write_json(output / "host.json", host)
    base.parent.write_json(output / "source.json", {
        "schema": "n83.root4-fixed30-root11-source/v1",
        "hashes": hashes,
        "binary_sha256": base.parent.digest(binary),
        "runner_sha256": base.parent.digest(Path(__file__)),
        "parent_runner_sha256": base.parent.digest(BASE),
    })
    if (host["mem_total_kib"] is None or host["mem_total_kib"] < 15_000_000
            or host["disk_free_bytes"] < 2 * 1024**3
            or len(host["visible_affinity"] or []) < 4):
        base.parent.write_json(output / "preflight.json", {
            "status": "REJECTED_HOST_CAPACITY", "host": host})
        return 2
    base.parent.write_json(output / "preflight.json", {
        "status": "PASS_CAPACITY_PREFLIGHT",
        "rss_cap_bytes": base.parent.RSS_CAP_BYTES})
    for k, variant, cap in ARMS:
        arm = f"k{k}_{variant}"
        folder = output / "runs" / arm
        four = variant.startswith("root4f30")
        eleven = variant in ("root4f30r11", "root11")
        env = {
            "PATH": os.environ.get("PATH", "/usr/bin:/bin"),
            "KIC_PACKED_ROOT_TABLE": "1",
            "KIC_PACKED_ROOT_TABLE_DENSE": "0",
            "KIC_PACKED_ROOT_TABLE_LOAD75": "0",
            "KIC_STATE28": "0", "KIC_STATE26": "0",
            "KIC_ROOT_ONLY22": "0",
            "KIC_ROOT_ONLY21": "0" if eleven else "1",
            "KIC_ROOT_ONLY11": "1" if eleven else "0",
            "KIC_INDEXED_ROOT5": "0" if four else "1",
            "KIC_INDEXED_ROOT4": "1" if four else "0",
            "KIC_ROOT_BLOOM": "0", "KIC_SYMMETRIC_S3": "1",
            "KIC_INDEXED5_LOAD75": "0",
            "KIC_INDEXED5_FIXED30": "0" if four else "1",
            "KIC_INDEXED4_FIXED30": "1" if four else "0",
            "KIC_INDEXED4_FIXED31": "0",
            "KIC_SPARSE_CANONICAL": "0", "KIC_RANK_THREADS": "4",
            "KIC_TARGET_THREADS": "1",
            "KIC_DUMP_BASE": str(folder / "base_dump.jsonl"),
        }
        command = [str(binary), f"construct:83:1:{k}", str(base.HOLDOUT),
                   "29", str(folder / "targets_result.jsonl")]
        print(f"starting {arm} at {base.parent.now()}", flush=True)
        usage = base.parent.one_run(folder, command, env, cap, group)
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
                          "host": base.parent.host_receipt(base.parent.cgroup_root())},
                         sort_keys=True))
        return 0
    if not args.binary or not args.output:
        parser.error("--binary and --output are required")
    return run(args.binary.resolve(), args.output.resolve())


if __name__ == "__main__":
    raise SystemExit(main())
