#!/usr/bin/env python3
"""Balanced N83 root4 Bloom factor-base screen."""

import argparse
import importlib.util
import json
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
PARENT = ROOT / "research/n83_root4_20261014/run.py"
HOLDOUT = HERE / "targets.jsonl"
ANSWERS = HERE / "answers.jsonl"
HOLDOUT_SHA256 = "71b93fd3a397ee2c1204dc14fc451894eee7537b99f50d75c98d65dee6b1e3d1"
ANSWERS_SHA256 = "5417c1a10a4c2240eaec304627eddfbd287e1d005f674b1a2b17b3b957aa53e5"
ARMS = (("root4k3350std_a", 3350, "root4", "standard", 1200),
        ("root4k3400std_a", 3400, "root4", "standard", 1200),
        ("root4k3425std_a", 3425, "root4", "standard", 1200),
        ("root11k3300wide3_a", 3300, "root11", "wide3", 1200),
        ("root11k3300wide3_b", 3300, "root11", "wide3", 1200),
        ("root4k3425std_b", 3425, "root4", "standard", 1200),
        ("root4k3400std_b", 3400, "root4", "standard", 1200),
        ("root4k3350std_b", 3350, "root4", "standard", 1200))
SOURCE_SHA256 = "f2745227d95e36ef45c73c5104a23ada6cebb0794405368581d2be17044845dc"
PARENT_SOURCE_SHA256 = "5e9fd826d22816b429a82e658c517d9387916aabe7660f708e4c29901cce647e"

spec = importlib.util.spec_from_file_location("n83_root4_runner", PARENT)
parent = importlib.util.module_from_spec(spec)
spec.loader.exec_module(parent)
if parent.EXPECTED[parent.SOURCE] != PARENT_SOURCE_SHA256:
    raise ValueError("parent runner source pin changed")
parent.EXPECTED[parent.SOURCE] = SOURCE_SHA256
parent.RSS_CAP_BYTES = 29 * 1024**3 // 2


def frozen_inputs():
    hashes = parent.check_inputs()
    if (parent.digest(HOLDOUT) != HOLDOUT_SHA256
            or parent.digest(ANSWERS) != ANSWERS_SHA256):
        raise ValueError("frozen balanced-order fixture changed")
    hashes[str(HOLDOUT.relative_to(ROOT))] = HOLDOUT_SHA256
    hashes[str(ANSWERS.relative_to(ROOT))] = ANSWERS_SHA256
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
        "schema": "n83.root4-bloom-source/v1",
        "hashes": hashes,
        "binary_sha256": parent.digest(binary),
        "panel_runner_sha256": parent.digest(Path(__file__)),
        "arm_order": [arm for arm, _, _, _, _ in ARMS],
    })
    if (host["mem_total_kib"] is None
            or host["mem_total_kib"] * 1024 < parent.RSS_CAP_BYTES + 256 * 1024**2
            or host["disk_free_bytes"] < 2 * 1024**3
            or len(host["visible_affinity"] or []) < 4):
        parent.write_json(output / "preflight.json", {
            "status": "REJECTED_HOST_CAPACITY", "host": host})
        return 2
    parent.write_json(output / "preflight.json", {
        "status": "PASS_CAPACITY_PREFLIGHT", "rss_cap_bytes": parent.RSS_CAP_BYTES})
    for arm, k, variant, profile, cap in ARMS:
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
            "KIC_ROOT_BLOOM": "1" if profile != "off" else "0",
            "KIC_SYMMETRIC_S3": "1",
            "KIC_INDEXED5_LOAD75": "1" if variant == "root22" else "0",
            "KIC_INDEXED5_FIXED30": "1" if variant in ("root11", "root21", "root5") else "0",
            "KIC_INDEXED4_FIXED31": "1" if variant == "root4" else "0",
            "KIC_SPARSE_CANONICAL": "1", "KIC_RANK_THREADS": "4",
            "KIC_TARGET_THREADS": "1", "KIC_DUMP_BASE": str(folder / "base_dump.jsonl"),
        }
        if profile != "off":
            env["KIC_ROOT_BLOOM_PROFILE"] = profile
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
