#!/usr/bin/env python3
"""Five-point cold N83 holdout; finalists derive only from replayed grid."""

import argparse
import importlib.util
import json
import os
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SCREEN = ROOT / "research/n83_e2e_factor_base_20261016/run.py"
TARGETS = HERE / "targets.jsonl"
ANSWERS = HERE / "answers.jsonl"
MANIFEST = HERE / "fixture_manifest.json"
SCREEN_ANALYSIS = HERE / "screen_analysis.json"
SELECTION = HERE / "selection.json"
TARGET_SHA256 = "722bc60cd010e533cd715a9bfb321fba52d950352c99684cd422b23af7f93d26"
ANSWER_SHA256 = "b4c2ae858770b7c8842703558984326bc0002b0f1deb3c898d5ac30562de8b96"
SCREEN_RUN_ID = 38940108434
SCREEN_HEAD_SHA = "551b9dfc632dce7a5ab8218b1b1dff6ae4b4d785"
RHO_SEED_BASE = 2026101700
IC_CAP_SECONDS = 480
RHO_CAP_SECONDS = 120

spec = importlib.util.spec_from_file_location("n83_e2e_screen", SCREEN)
screen = importlib.util.module_from_spec(spec)
spec.loader.exec_module(screen)
dynamic = screen.dynamic
parent = screen.parent
cold = screen.cold


def fixture():
    hashes = screen.frozen_inputs()
    if (parent.digest(TARGETS) != TARGET_SHA256
            or parent.digest(ANSWERS) != ANSWER_SHA256):
        raise ValueError("five-point holdout fixture differs")
    manifest = json.loads(MANIFEST.read_text())
    if (manifest["targets_sha256"] != TARGET_SHA256
            or manifest["answers_sha256"] != ANSWER_SHA256
            or manifest["target_count"] != 5
            or manifest["domain"] != "n83-e2e-cold-holdout-five-20261016-v1"
            or len(TARGETS.read_text().splitlines()) != 5
            or len(ANSWERS.read_text().splitlines()) != 5):
        raise ValueError("five-point holdout manifest differs")
    return hashes


def finalists():
    analysis = json.loads(SCREEN_ANALYSIS.read_text())
    selection = json.loads(SELECTION.read_text())
    if (analysis["schema"] != "n83.e2e-factor-base-analysis/v1"
            or analysis["status"] != "SCREEN_COMPLETE"
            or analysis["evidence_status"] != "PASS_INDEPENDENT_REPLAY"
            or analysis["aa_control"]["status"] != "PASS_10_PERCENT"
            or selection["schema"] != "n83.e2e-holdout-selection/v1"
            or selection["screen_analysis_sha256"] != parent.digest(SCREEN_ANALYSIS)
            or selection["screen_source_sha256"] != screen.SOURCE_SHA256
            or selection["screen_run_id"] != SCREEN_RUN_ID
            or selection["screen_source_commit"] != SCREEN_HEAD_SHA
            or not selection["screen_archive_sha256"]):
        raise ValueError("exposed grid did not admit a frozen holdout selection")
    expected = [int(k) for k in analysis["descriptive_order"]
                if int(k) != screen.CONTROL_K][:3]
    if (len(expected) != 3 or selection["challenger_ks"] != expected
            or any(k not in screen.K_GRID for k in expected)):
        raise ValueError("holdout challengers differ from prespecified rule")
    return expected, selection


def schedule(challengers):
    rows = []
    for i in range(1, 6):
        rotated = challengers[(i - 1) % 3:] + challengers[:(i - 1) % 3]
        rows.append({"point_index": i,
                     "ic_order": ["k1152_a", *(f"k{k}" for k in rotated),
                                  "k1152_b"],
                     "rho_first": i % 2 == 0})
    return rows


def point_inputs(output):
    folder = output / "inputs"
    folder.mkdir()
    paths = {}
    for i, line in enumerate(TARGETS.read_text().splitlines(), start=1):
        path = folder / f"point_{i}.jsonl"
        path.write_text(line + "\n")
        paths[i] = path
    return paths


def run_ic(output, binary, point_index, point_file, arm, group):
    k = screen.CONTROL_K if arm.startswith("k1152_") else int(arm[1:])
    name = f"p{point_index}_{arm}"
    folder = output / "runs" / name
    command = [str(binary), f"construct:83:1:{k}", str(point_file),
               "29", str(folder / "targets_result.jsonl")]
    print(f"starting {name} at {parent.now()}", flush=True)
    usage = parent.one_run(folder, command,
                           dynamic.environment(folder, "root11_dynamic"),
                           IC_CAP_SECONDS, group)
    print(f"{name}: {usage['status']} after {usage['wall_seconds']:.3f}s",
          flush=True)


def run_rho(output, binary, point_index, point_file, group):
    name = f"p{point_index}_rho"
    folder = output / "runs" / name
    command = [str(binary), str(point_file), str(RHO_SEED_BASE + point_index),
               str(cold.RHO_CONFIG), str(folder / "receipt.json")]
    usage = parent.one_run(folder, command,
                           {"PATH": os.environ.get("PATH", "/usr/bin:/bin")},
                           RHO_CAP_SECONDS, group)
    print(f"{name}: {usage['status']} after {usage['wall_seconds']:.3f}s",
          flush=True)


def run(binary, rho_binary, output):
    hashes = fixture()
    challengers, selection = finalists()
    if output.exists():
        raise ValueError("output directory already exists")
    output.mkdir(parents=True)
    group = parent.cgroup_root()
    host = parent.host_receipt(group)
    host["schema"] = "n83.e2e-holdout-host/v1"
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
    paths = point_inputs(output)
    parent.write_json(output / "source.json", {
        "schema": "n83.e2e-holdout-source/v1", "hashes": hashes,
        "source_sha256": screen.SOURCE_SHA256,
        "binary_sha256": parent.digest(binary),
        "rho_binary_sha256": parent.digest(rho_binary),
        "runner_sha256": parent.digest(Path(__file__)),
        "selection_sha256": parent.digest(SELECTION),
        "screen_analysis_sha256": parent.digest(SCREEN_ANALYSIS),
        "challenger_ks": challengers, "screen_archive_sha256":
            selection["screen_archive_sha256"],
        "schedule": schedule(challengers),
        "point_input_sha256": {str(i): parent.digest(path)
                               for i, path in paths.items()},
        "rho_seed_base": RHO_SEED_BASE,
    })
    for row in schedule(challengers):
        i = row["point_index"]
        if row["rho_first"]:
            run_rho(output, rho_binary, i, paths[i], group)
        for arm in row["ic_order"]:
            run_ic(output, binary, i, paths[i], arm, group)
        if not row["rho_first"]:
            run_rho(output, rho_binary, i, paths[i], group)
    return 0


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--preflight", action="store_true")
    parser.add_argument("--binary", type=Path)
    parser.add_argument("--rho-binary", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.preflight:
        print(json.dumps({"hashes": fixture(), "selection_pending":
                          not (SELECTION.exists() and SCREEN_ANALYSIS.exists())},
                         sort_keys=True))
        return 0
    if not args.binary or not args.rho_binary or not args.output:
        parser.error("--binary, --rho-binary and --output are required")
    return run(args.binary.resolve(), args.rho_binary.resolve(),
               args.output.resolve())


if __name__ == "__main__":
    raise SystemExit(main())
