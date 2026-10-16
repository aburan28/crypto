#!/usr/bin/env python3
"""Run the frozen fresh-point confirmation of unresolved N83 finalists."""

import argparse
import hashlib
import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
SOURCE = REPO / "research/n83_final_tournament_20261015/run.py"
ADMIT = HERE / "admit.py"
ADMISSION = HERE / "admission.json"
TARGETS = HERE / "targets.jsonl"
ANSWERS = HERE / "answers.jsonl"
TARGET_SHA256 = "66c4683e252f1cd321b4b6c7c9d06744c82140222ddce160faa4e7a0e0984b7e"
ANSWER_SHA256 = "e1207d1be9a7689c42209ed50f048325812436d5f7e6f313b5e5e8a04a5730b1"
RSS_CAP_BYTES = 29 * 1024**3 // 2
CAP_SECONDS = 1500

rspec = importlib.util.spec_from_file_location("n83_confirmation_source", SOURCE)
source = importlib.util.module_from_spec(rspec)
rspec.loader.exec_module(source)
aspec = importlib.util.spec_from_file_location("n83_confirmation_admit", ADMIT)
admit = importlib.util.module_from_spec(aspec)
aspec.loader.exec_module(admit)
parent = source.parent


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def frozen():
    original, hashes = source.frozen()
    if (digest(TARGETS) != TARGET_SHA256
            or digest(ANSWERS) != ANSWER_SHA256
            or not ADMISSION.is_file()
            or json.loads(ADMISSION.read_text()) != admit.admission()):
        raise ValueError("confirmation fixture or sealed admission differs")
    admission = json.loads(ADMISSION.read_text())
    original_by_id = {row["candidate_id"]: row for row in original["candidate_order"]}
    for row in admission["candidate_order"]:
        expected = original_by_id[row["candidate_id"]].copy()
        expected["position"] = row["position"]
        expected["arm_names"] = row["arm_names"]
        if row != expected:
            raise ValueError("confirmed candidate differs from source manifest")
    hashes[str(TARGETS.relative_to(REPO))] = TARGET_SHA256
    hashes[str(ANSWERS.relative_to(REPO))] = ANSWER_SHA256
    hashes[str(ADMISSION.relative_to(REPO))] = digest(ADMISSION)
    hashes[str(ADMIT.relative_to(REPO))] = digest(ADMIT)
    return admission, hashes


def arms(candidates):
    return source.arms(candidates)


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
        "schema": "n83.edge-confirmation-source/v1",
        "hashes": hashes, "binary_sha256": digest(binary),
        "panel_runner_sha256": digest(Path(__file__)),
        "admission_sha256": digest(ADMISSION),
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
    parent.write_json(output / "calibration_before.json", source.calibration())
    for arm, row in order:
        folder = output / "runs" / arm
        command = [str(binary), f"construct:83:1:{row['k']}", str(TARGETS),
                   "29", str(folder / "targets_result.jsonl")]
        print(f"starting {arm} at {parent.now()}", flush=True)
        parent.write_json(output / f"{arm}_system_before.json", source.system_counters())
        usage = parent.one_run(folder, command, source.environment(folder, row),
                               CAP_SECONDS, group)
        parent.write_json(output / f"{arm}_system_after.json", source.system_counters())
        print(f"{arm}: {usage['status']} after {usage['wall_seconds']:.3f}s", flush=True)
    parent.write_json(output / "calibration_after.json", source.calibration())
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
