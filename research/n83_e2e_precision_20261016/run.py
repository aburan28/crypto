#!/usr/bin/env python3
"""Run the sealed sixteen-point N83 K1104/K1120 precision panel."""

import argparse
import contextlib
import importlib.util
import io
import json
import os
import subprocess
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
HOLDOUT_RUN = ROOT / "research/n83_e2e_holdout_20261016/run.py"
GENERATOR = HERE / "generate.py"
TARGETS = HERE / "targets.jsonl"
ANSWERS = HERE / "answers.jsonl"
MANIFEST = HERE / "fixture_manifest.json"
RHO_SEED_BASE = 2026101800
IC_CAP_SECONDS = 480
RHO_CAP_SECONDS = 120
K_VALUES = (1104, 1120)


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    value = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(value)
    return value


holdout = module("n83_precision_parent_run", HOLDOUT_RUN)
generator = module("n83_precision_fixture", GENERATOR)
screen = holdout.screen
parent = holdout.parent
dynamic = holdout.dynamic
cold = holdout.cold


def fixture():
    with contextlib.redirect_stdout(io.StringIO()):
        generator.main(True)
    manifest = json.loads(MANIFEST.read_text())
    if (manifest["target_count"] != 16
            or parent.digest(TARGETS) != manifest["targets_sha256"]
            or parent.digest(ANSWERS) != manifest["answers_sha256"]):
        raise ValueError("precision fixture hashes differ")
    return screen.frozen_inputs(), manifest


def schedule():
    rows = []
    for point_index in range(1, 17):
        order = (["k1104_a", "k1120_a", "k1120_b", "k1104_b"]
                 if point_index % 2 else
                 ["k1120_a", "k1104_a", "k1104_b", "k1120_b"])
        rows.append({"point_index": point_index, "ic_order": order,
                     "rho_first": point_index % 2 == 0})
    return rows


def point_inputs(output):
    folder = output / "inputs"
    folder.mkdir()
    paths = {}
    for i, line in enumerate(TARGETS.read_text().splitlines(), start=1):
        path = folder / f"point_{i}.jsonl"
        path.write_text(line + "\n")
        paths[i] = path
    if len(paths) != 16:
        raise ValueError("precision fixture must contain sixteen points")
    return paths


def run_ic(output, binary, point_index, point_file, arm, group):
    k = int(arm.split("_", 1)[0][1:])
    if k not in K_VALUES:
        raise ValueError(f"unexpected precision K: {k}")
    name = f"p{point_index}_{arm}"
    folder = output / "runs" / name
    command = [str(binary), f"construct:83:1:{k}", str(point_file),
               "29", str(folder / "targets_result.jsonl")]
    print(f"starting {name} at {parent.now()}", flush=True)
    usage = holdout.run_with_telemetry(
        folder, command, dynamic.environment(folder, "root11_dynamic"),
        IC_CAP_SECONDS, group)
    print(f"{name}: {usage['status']} after {usage['wall_seconds']:.3f}s",
          flush=True)


def run_rho(output, binary, point_index, point_file, group):
    name = f"p{point_index}_rho"
    folder = output / "runs" / name
    command = [str(binary), str(point_file), str(RHO_SEED_BASE + point_index),
               str(cold.RHO_CONFIG), str(folder / "receipt.json")]
    print(f"starting {name} at {parent.now()}", flush=True)
    usage = holdout.run_with_telemetry(
        folder, command, {"PATH": os.environ.get("PATH", "/usr/bin:/bin")},
        RHO_CAP_SECONDS, group)
    print(f"{name}: {usage['status']} after {usage['wall_seconds']:.3f}s",
          flush=True)


def run(binary, rho_binary, output):
    hashes, manifest = fixture()
    if output.exists():
        raise ValueError("output directory already exists")
    output.mkdir(parents=True)
    group = parent.cgroup_root()
    host = parent.host_receipt(group)
    host["schema"] = "n83.e2e-precision-host/v1"
    host["topology"] = holdout.topology_receipt(host["visible_affinity"] or [])
    host["cpu_telemetry_before"] = holdout.cpu_telemetry()
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
        "schema": "n83.e2e-precision-source/v1", "hashes": hashes,
        "source_sha256": screen.SOURCE_SHA256,
        "binary_sha256": parent.digest(binary),
        "rho_binary_sha256": parent.digest(rho_binary),
        "runner_sha256": parent.digest(Path(__file__)),
        "generator_sha256": parent.digest(GENERATOR),
        "fixture_manifest_sha256": parent.digest(MANIFEST),
        "targets_sha256": manifest["targets_sha256"],
        "answers_sha256": manifest["answers_sha256"],
        "k_values": list(K_VALUES), "schedule": schedule(),
        "point_input_sha256": {str(i): parent.digest(path)
                               for i, path in paths.items()},
        "rho_seed_base": RHO_SEED_BASE,
    })
    for row in schedule():
        i = row["point_index"]
        if row["rho_first"]:
            run_rho(output, rho_binary, i, paths[i], group)
        for arm in row["ic_order"]:
            run_ic(output, binary, i, paths[i], arm, group)
        if not row["rho_first"]:
            run_rho(output, rho_binary, i, paths[i], group)
    parent.write_json(output / "host_telemetry_after.json",
                      holdout.cpu_telemetry())
    return 0


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--preflight", action="store_true")
    parser.add_argument("--binary", type=Path)
    parser.add_argument("--rho-binary", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.preflight:
        hashes, manifest = fixture()
        print(json.dumps({"hashes": hashes, "fixture": manifest,
                          "schedule": schedule()}, sort_keys=True))
        return 0
    if not args.binary or not args.rho_binary or not args.output:
        parser.error("--binary, --rho-binary and --output are required")
    return run(args.binary.resolve(), args.rho_binary.resolve(),
               args.output.resolve())


if __name__ == "__main__":
    raise SystemExit(main())
