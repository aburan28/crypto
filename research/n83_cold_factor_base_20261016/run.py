#!/usr/bin/env python3
"""Bounded one-target cold N83 factor-base geometry screen."""

import argparse
import importlib.util
import json
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
BASE_RUNNER = ROOT / "research/n83_compact4_crosslayout_20261016/run.py"
TARGETS = HERE / "targets.jsonl"
ANSWERS = HERE / "answers.jsonl"
RESERVE_TARGETS = HERE / "reserve_targets.jsonl"
RESERVE_ANSWERS = HERE / "reserve_answers.jsonl"
ADMISSION = ROOT / "research/n83_compact4_crosslayout_results_20261016/selection.json"
RHO_CONFIG = ROOT / "research/n83_dense_root_20261013/rho_config.json"
FROZEN_HASHES = {
    "research/n83_cold_factor_base_20261016/generate.py":
        "1807dadb4dc423850f32742faa03beb5774dacb4d3c17618502b6debb623ca58",
    "research/n83_cold_factor_base_20261016/prior_manifest.json":
        "f0b29cb6e7d0081711ccb4b5ff54455f06198a2a3562d4654fbe7e81d07b1c8b",
    "research/n83_cold_factor_base_20261016/targets.jsonl":
        "2c068d2f659588edcdbaa685fc73bb146e2bfee3be123d377d537658f8021c6d",
    "research/n83_cold_factor_base_20261016/answers.jsonl":
        "b71c5a8bfdae603931fa6c03fc6f6faa3b261288cf787ac682bf8a3e3afe262d",
    "research/n83_cold_factor_base_20261016/reserve_targets.jsonl":
        "6ab08db41c5e86c96162f538a0c45d53331c2f9f823a66eea9bf0b458dd3554d",
    "research/n83_cold_factor_base_20261016/reserve_answers.jsonl":
        "b21d3a1d54c701255a079a3f1d88940f46563a2f857303b4e3bb98ac40483f39",
    "research/n83_compact4_crosslayout_results_20261016/selection.json":
        "75a8f9f7638839d1d3d37e91576a490597f97fb0dfa5cab8b7526d38176415be",
}
ARMS = (
    ("k2543root22_a", 2543, "root22", "off", 900),
    ("k256root22", 256, "root22", "off", 900),
    ("k512root22", 512, "root22", "off", 900),
    ("k768root22", 768, "root22", "off", 900),
    ("k1024root22", 1024, "root22", "off", 900),
    ("k1536root22", 1536, "root22", "off", 900),
    ("k2048root22", 2048, "root22", "off", 900),
    ("k3450root11wide3", 3450, "root11", "wide3", 900),
    ("k2543root22_b", 2543, "root22", "off", 900),
)
RHO_SEED = "2026101601"
RHO_CAP_SECONDS = 120
RSS_CAP_BYTES = 61 * 1024**3 // 4

spec = importlib.util.spec_from_file_location("n83_compact4_runner", BASE_RUNNER)
crossrun = importlib.util.module_from_spec(spec)
spec.loader.exec_module(crossrun)
parent = crossrun.parent
parent.RSS_CAP_BYTES = RSS_CAP_BYTES


def frozen_inputs():
    hashes = parent.check_inputs()
    for relative, expected in FROZEN_HASHES.items():
        path = ROOT / relative
        if parent.digest(path) != expected:
            raise ValueError(f"frozen cold-sweep input differs: {relative}")
        hashes[relative] = expected
    if (len(TARGETS.read_text().splitlines()) != 1
            or len(ANSWERS.read_text().splitlines()) != 1
            or len(RESERVE_TARGETS.read_text().splitlines()) != 4
            or len(RESERVE_ANSWERS.read_text().splitlines()) != 4):
        raise ValueError("one primary and four reserve target/answer rows required")
    admission = json.loads(ADMISSION.read_text())
    if (admission["status"] != "BOUNDED_UNKNOWN"
            or admission["descriptive_lowest_pipeline"] != "root11k3450wide3"):
        raise ValueError("prior cross-layout decision differs")
    return hashes


def environment(folder, variant, profile):
    root22 = variant == "root22"
    if not root22 and (variant, profile) != ("root11", "wide3"):
        raise ValueError("unsupported cold-sweep pipeline")
    flags = {
        "KIC_PACKED_ROOT_TABLE": "1", "KIC_PACKED_ROOT_TABLE_DENSE": "0",
        "KIC_PACKED_ROOT_TABLE_LOAD75": "0", "KIC_STATE28": "0",
        "KIC_STATE26": "0", "KIC_ROOT_ONLY22": "1" if root22 else "0",
        "KIC_ROOT_ONLY21": "0", "KIC_ROOT_ONLY11": "0" if root22 else "1",
        "KIC_INDEXED_ROOT5": "1", "KIC_INDEXED_ROOT4": "0",
        "KIC_ROOT_BLOOM": "0" if root22 else "1",
        "KIC_SYMMETRIC_S3": "1", "KIC_INDEXED5_LOAD75": "1" if root22 else "0",
        "KIC_INDEXED5_FIXED30": "0" if root22 else "1",
        "KIC_INDEXED4_FIXED31": "0", "KIC_SPARSE_CANONICAL": "1",
        "KIC_RANK_THREADS": "4", "KIC_TARGET_THREADS": "1",
        "KIC_DUMP_BASE": str(folder / "base_dump.jsonl"),
    }
    if profile != "off":
        flags["KIC_ROOT_BLOOM_PROFILE"] = profile
    return {"PATH": os.environ.get("PATH", "/usr/bin:/bin"), **flags}


def run(binary, rho_binary, output):
    hashes = frozen_inputs()
    if output.exists():
        raise ValueError("output directory already exists")
    output.mkdir(parents=True)
    group = parent.cgroup_root()
    host = parent.host_receipt(group)
    host["schema"] = "n83.cold-factor-base-host/v1"
    parent.write_json(output / "host.json", host)
    parent.write_json(output / "source.json", {
        "schema": "n83.cold-factor-base-source/v1", "hashes": hashes,
        "binary_sha256": parent.digest(binary),
        "rho_binary_sha256": parent.digest(rho_binary),
        "runner_sha256": parent.digest(Path(__file__)),
        "arm_order": [arm for arm, _, _, _, _ in ARMS],
        "rho_seed": RHO_SEED,
    })
    memory_max = (host["cgroup_before"]["memory.max"] or "").strip()
    memory_too_small = memory_max.isdigit() and int(memory_max) < RSS_CAP_BYTES
    if (host["mem_total_kib"] is None
            or host["mem_total_kib"] * 1024 < RSS_CAP_BYTES + 256 * 1024**2
            or host["disk_free_bytes"] < 2 * 1024**3
            or len(host["visible_affinity"] or []) < 4 or memory_too_small):
        parent.write_json(output / "preflight.json", {
            "status": "REJECTED_HOST_CAPACITY", "host": host})
        return 2
    parent.write_json(output / "preflight.json", {
        "status": "PASS_CAPACITY_PREFLIGHT", "rss_cap_bytes": RSS_CAP_BYTES})
    for arm, k, variant, profile, cap in ARMS:
        folder = output / "runs" / arm
        command = [str(binary), f"construct:83:1:{k}", str(TARGETS), "29",
                   str(folder / "targets_result.jsonl")]
        print(f"starting {arm} at {parent.now()}", flush=True)
        usage = parent.one_run(folder, command, environment(folder, variant, profile), cap,
                               group)
        print(f"{arm}: {usage['status']} after {usage['wall_seconds']:.3f}s", flush=True)
    folder = output / "runs" / "rho"
    command = [str(rho_binary), str(TARGETS), RHO_SEED, str(RHO_CONFIG),
               str(folder / "receipt.json")]
    print(f"starting same-point rho at {parent.now()}", flush=True)
    usage = parent.one_run(folder, command, {"PATH": os.environ.get("PATH", "/usr/bin:/bin")},
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
    return run(args.binary.resolve(), args.rho_binary.resolve(), args.output.resolve())


if __name__ == "__main__":
    raise SystemExit(main())
