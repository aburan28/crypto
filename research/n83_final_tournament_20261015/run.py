#!/usr/bin/env python3
"""Run the frozen N83 factor-base finalists on one fresh public-point panel."""

import argparse
import hashlib
import importlib.util
import json
import os
import resource
import statistics
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
PARENT = ROOT / "research/n83_root4_20261014/run.py"
ADMISSION = HERE / "admission.json"
TARGETS = HERE / "targets.jsonl"
ANSWERS = HERE / "answers.jsonl"
TARGET_SHA256 = "cfa943ae538fcc41d2af8bb8fd1c8fc0f29adc340b55644e90382da25b213424"
ANSWER_SHA256 = "1dbd4250796de4434badc2458089d7ad9fe63cc59d7b9d3f6825ad431476a13b"
SOURCE_SHA256 = "f2745227d95e36ef45c73c5104a23ada6cebb0794405368581d2be17044845dc"
PARENT_SHA256 = "5e9fd826d22816b429a82e658c517d9387916aabe7660f708e4c29901cce647e"
RSS_CAP_BYTES = 29 * 1024**3 // 2
CAP_SECONDS = 1500

spec = importlib.util.spec_from_file_location("n83_tournament_parent", PARENT)
parent = importlib.util.module_from_spec(spec)
spec.loader.exec_module(parent)
if parent.EXPECTED[parent.SOURCE] != PARENT_SHA256:
    raise ValueError("parent runner source pin changed")
parent.EXPECTED[parent.SOURCE] = SOURCE_SHA256
parent.RSS_CAP_BYTES = RSS_CAP_BYTES


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def frozen():
    hashes = parent.check_inputs()
    if digest(TARGETS) != TARGET_SHA256 or digest(ANSWERS) != ANSWER_SHA256:
        raise ValueError("final public target or answer sidecar changed")
    hashes[str(TARGETS.relative_to(ROOT))] = TARGET_SHA256
    hashes[str(ANSWERS.relative_to(ROOT))] = ANSWER_SHA256
    if not ADMISSION.is_file():
        raise ValueError("source-screen admission must be sealed before execution")
    admission = json.loads(ADMISSION.read_text())
    candidates = admission["candidate_order"]
    if (admission["schema"] != "n83.final-tournament-admission/v1"
            or admission["status"] != "FROZEN_ADMISSION"
            or not 2 <= len(candidates) <= 6
            or admission["target_sha256"] != TARGET_SHA256
            or admission["answer_sidecar_sha256"] != ANSWER_SHA256
            or [row["position"] for row in candidates] != list(range(len(candidates)))
            or len({row["pipeline"] for row in candidates}) != len(candidates)
            or admission["mirrored_arm_order"] !=
            [row["pipeline"] for row in candidates] +
            [row["pipeline"] for row in reversed(candidates)]):
        raise ValueError("frozen candidate admission malformed")
    for row in candidates:
        origin = "direct" if row["origin"] == "control" else row["origin"]
        source_screen = next(x for x in admission["source_screens"] if x["screen"] == origin)
        result_dir = {
            "direct": "n83_best_e2e_results_20261015",
            "legacy": "n83_legacy_challenger_results_20261015",
            "root4_bloom": "n83_root4_bloom_results_20261015",
            "root11_highk": "n83_root11_highk_results_20261015",
        }[origin]
        manifest = ROOT / "research" / result_dir / (
            f"candidate_K{row['k']}_{row['variant']}_{row['bloom_profile']}.json")
        record = json.loads(manifest.read_text())
        if (row["source_sha256"] != SOURCE_SHA256
                or row["candidate_id"] != record["candidate_id"]
                or row["candidate_manifest_sha256"] != digest(manifest)
                or row["base_hash"] != record["factor_base"]["source_base_hash_blake3"]
                or row["source_sha256"] != record["implementation"]["source_sha256"]["ic"]
                or source_screen["source_commit"] != record["implementation"]["source_commit"]):
            raise ValueError(f"{row['pipeline']}: source manifest changed")
    hashes[str(ADMISSION.relative_to(ROOT))] = digest(ADMISSION)
    return admission, hashes


def arms(candidates):
    first = [(f"{row['pipeline']}_a", row) for row in candidates]
    second = [(f"{row['pipeline']}_b", row) for row in reversed(candidates)]
    return first + second


def environment(folder, row):
    variant, profile = row["variant"], row["bloom_profile"]
    env = {
        "PATH": os.environ.get("PATH", "/usr/bin:/bin"),
        "KIC_PACKED_ROOT_TABLE": "1", "KIC_PACKED_ROOT_TABLE_DENSE": "0",
        "KIC_PACKED_ROOT_TABLE_LOAD75": "0", "KIC_STATE28": "0",
        "KIC_STATE26": "0",
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
    return env


def calibration():
    durations = []
    digests = []
    for repetition in range(3):
        state = (b"n83-final-tournament-calibration-v1" +
                 repetition.to_bytes(4, "little"))
        before_cpu = resource.getrusage(resource.RUSAGE_SELF)
        started = time.monotonic_ns()
        for _ in range(1_000_000):
            state = hashlib.sha256(state).digest()
        durations.append((time.monotonic_ns() - started) / 1_000_000)
        after_cpu = resource.getrusage(resource.RUSAGE_SELF)
        digests.append({"sha256": state.hex(),
                        "user_cpu_seconds": after_cpu.ru_utime - before_cpu.ru_utime,
                        "system_cpu_seconds": after_cpu.ru_stime - before_cpu.ru_stime})
    return {"sample_ms": durations, "median_ms": statistics.median(durations),
            "checksums_and_cpu": digests, "recorded_utc": parent.now()}


def system_counters():
    return {"recorded_utc": parent.now(),
            "proc_stat": parent.read_text(Path("/proc/stat")),
            "proc_interrupts": parent.read_text(Path("/proc/interrupts"))}


def run(binary, output):
    admission, hashes = frozen()
    if output.exists():
        raise ValueError("output directory already exists")
    output.mkdir(parents=True)
    group = parent.cgroup_root()
    host = parent.host_receipt(group)
    parent.write_json(output / "host.json", host)
    order = arms(admission["candidate_order"])
    parent.write_json(output / "source.json", {
        "schema": "n83.final-tournament-source/v1",
        "hashes": hashes, "binary_sha256": digest(binary),
        "panel_runner_sha256": digest(Path(__file__)),
        "admission_sha256": hashes[str(ADMISSION.relative_to(ROOT))],
        "arm_order": [arm for arm, _ in order],
        "candidate_ids": [row["candidate_id"] for row in admission["candidate_order"]],
    })
    if (host["mem_total_kib"] is None
            or host["mem_total_kib"] * 1024 < RSS_CAP_BYTES + 256 * 1024**2
            or host["disk_free_bytes"] < 2 * 1024**3
            or len(host["visible_affinity"] or []) < 4):
        parent.write_json(output / "preflight.json", {
            "status": "REJECTED_HOST_CAPACITY", "host": host})
        return 2
    parent.write_json(output / "preflight.json", {
        "status": "PASS_CAPACITY_PREFLIGHT", "rss_cap_bytes": RSS_CAP_BYTES})
    parent.write_json(output / "calibration_before.json", calibration())
    for arm, row in order:
        folder = output / "runs" / arm
        command = [str(binary), f"construct:83:1:{row['k']}", str(TARGETS),
                   "29", str(folder / "targets_result.jsonl")]
        print(f"starting {arm} at {parent.now()}", flush=True)
        parent.write_json(output / f"{arm}_system_before.json", system_counters())
        usage = parent.one_run(folder, command, environment(folder, row),
                               CAP_SECONDS, group)
        parent.write_json(output / f"{arm}_system_after.json", system_counters())
        print(f"{arm}: {usage['status']} after {usage['wall_seconds']:.3f}s",
              flush=True)
    parent.write_json(output / "calibration_after.json", calibration())
    return 0


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--preflight", action="store_true")
    parser.add_argument("--binary", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.preflight:
        admission, hashes = frozen()
        print(json.dumps({"hashes": hashes,
                          "arms": [(arm, row["candidate_id"]) for arm, row in
                                   arms(admission["candidate_order"])],
                          "host": parent.host_receipt(parent.cgroup_root())},
                         sort_keys=True))
        return 0
    if not args.binary or not args.output:
        parser.error("--binary and --output are required")
    return run(args.binary.resolve(), args.output.resolve())


if __name__ == "__main__":
    raise SystemExit(main())
