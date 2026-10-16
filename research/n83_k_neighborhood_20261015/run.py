#!/usr/bin/env python3
"""Frozen four-arm K-neighborhood screen on 4000 N83 public points."""

import argparse
import importlib.util
import json
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PARENT = ROOT / "research/n83_root4_20261014/run.py"
HOLDOUT = ROOT / "research/n83_factor_base_confirm_20261015/targets.jsonl"
HOLDOUT_SHA256 = "451df8524b2c356dd4a523260620c4dff406975a122762b98ddddb7987b75c2c"
ARMS = ((3200, "root11", 1200), (3250, "root11", 1200),
        (3300, "root11", 1200), (3350, "root11", 1200))
SOURCE_SHA256 = "a1eb3b937106bb7cfdee6c54d6422ceb126b628462361c9ff7a127127f230c89"
PARENT_SOURCE_SHA256 = "5e9fd826d22816b429a82e658c517d9387916aabe7660f708e4c29901cce647e"

spec = importlib.util.spec_from_file_location("n83_root4_runner", PARENT)
parent = importlib.util.module_from_spec(spec)
spec.loader.exec_module(parent)
if parent.EXPECTED[parent.SOURCE] != PARENT_SOURCE_SHA256:
    raise ValueError("parent runner source pin changed")
parent.EXPECTED[parent.SOURCE] = SOURCE_SHA256


def frozen_inputs():
    hashes = parent.check_inputs()
    if parent.digest(HOLDOUT) != HOLDOUT_SHA256:
        raise ValueError("frozen 4000-point neighborhood fixture changed")
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
        "schema": "n83.k-neighborhood-source/v1",
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
