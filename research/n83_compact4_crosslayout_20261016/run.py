#!/usr/bin/env python3
"""Fresh 8,000-point six-arm compact Bloom and cross-layout comparison."""

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
HOLDOUT_SHA256 = "6c5558613d4f1115dc7bdb09b74d569645c0e85328ddad6555a2484317803a08"
ANSWERS_SHA256 = "82decefca05b868e2e9753de3f4892e9bcccfab97ed2c8f6cfbdfd871e13d847"
ADMISSION = ROOT / "research/n83_root4_highk_confirm_results_20261016/selection.json"
ADMISSION_SHA256 = "aaccd4b4cc8140bdae5ba1193c459c01780aef19a58486905664f3005f6af141"
ARMS = (("root4k3500std_a", 3500, "root4", "standard", 1500),
        ("root4k3596compact4_a", 3596, "root4", "compact4", 1500),
        ("root11k3450wide3_a", 3450, "root11", "wide3", 1500),
        ("root11k3450wide3_b", 3450, "root11", "wide3", 1500),
        ("root4k3596compact4_b", 3596, "root4", "compact4", 1500),
        ("root4k3500std_b", 3500, "root4", "standard", 1500))
SOURCE_SHA256 = "f78382ada0ee4eabfa6820d805d0e51598dbed82154b43ea309bcd489305e2aa"
PARENT_SOURCE_SHA256 = "5e9fd826d22816b429a82e658c517d9387916aabe7660f708e4c29901cce647e"

spec = importlib.util.spec_from_file_location("n83_root4_runner", PARENT)
parent = importlib.util.module_from_spec(spec)
spec.loader.exec_module(parent)
if parent.EXPECTED[parent.SOURCE] != PARENT_SOURCE_SHA256:
    raise ValueError("parent runner source pin changed")
parent.EXPECTED[parent.SOURCE] = SOURCE_SHA256
parent.RSS_CAP_BYTES = 61 * 1024**3 // 4


def frozen_inputs():
    hashes = parent.check_inputs()
    if (parent.digest(HOLDOUT) != HOLDOUT_SHA256
            or parent.digest(ANSWERS) != ANSWERS_SHA256
            or parent.digest(ADMISSION) != ADMISSION_SHA256):
        raise ValueError("frozen cross-layout fixture changed")
    admission = json.loads(ADMISSION.read_text())
    if (admission["selected_pipeline"] != "root4k3500std"
            or set(admission["balanced_comparison"]["pipeline_mean_ms"])
               != {"root4k3500std", "root4k3596std"}):
        raise ValueError("source panel selection differs")
    hashes[str(HOLDOUT.relative_to(ROOT))] = HOLDOUT_SHA256
    hashes[str(ANSWERS.relative_to(ROOT))] = ANSWERS_SHA256
    hashes[str(ADMISSION.relative_to(ROOT))] = ADMISSION_SHA256
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
        "schema": "n83.compact4-crosslayout-source/v1",
        "hashes": hashes,
        "binary_sha256": parent.digest(binary),
        "panel_runner_sha256": parent.digest(Path(__file__)),
        "arm_order": [arm for arm, _, _, _, _ in ARMS],
    })
    memory_max = (host["cgroup_before"]["memory.max"] or "").strip()
    memory_too_small = memory_max.isdigit() and int(memory_max) < parent.RSS_CAP_BYTES
    if (host["mem_total_kib"] is None
            or host["mem_total_kib"] * 1024 < parent.RSS_CAP_BYTES + 256 * 1024**2
            or host["disk_free_bytes"] < 2 * 1024**3
            or len(host["visible_affinity"] or []) < 4 or memory_too_small):
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
