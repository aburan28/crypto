#!/usr/bin/env python3
"""Fresh cold processes for four frozen N83 reserve points and four finalists."""

import argparse
import importlib.util
import json
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
DYNAMIC_RUNNER = ROOT / "research/n83_root11_dynamic_20261016/run.py"

spec = importlib.util.spec_from_file_location("n83_dynamic_runner", DYNAMIC_RUNNER)
dynamic = importlib.util.module_from_spec(spec)
spec.loader.exec_module(dynamic)
parent = dynamic.parent
cold = dynamic.cold
SOURCE_SHA256 = dynamic.SOURCE_SHA256
TARGETS = cold.RESERVE_TARGETS
RHO_SEED_BASE = 2026101601
IC_CAP_SECONDS = 900
RHO_CAP_SECONDS = 120

PIPELINES = {
    "k1024root22": (1024, "root22_dynamic"),
    "k1152root22": (1152, "root22_dynamic"),
    "k1152root11": (1152, "root11_dynamic"),
    "k1536root22": (1536, "root22_dynamic"),
}
MIDDLE = {
    1: ("k1152root22", "k1152root11", "k1536root22"),
    2: ("k1152root11", "k1536root22", "k1152root22"),
    3: ("k1536root22", "k1152root22", "k1152root11"),
    4: ("k1152root22", "k1536root22", "k1152root11"),
}


def schedule():
    return [{"point_index": i,
             "ic_order": ["k1024root22_a", *MIDDLE[i], "k1024root22_b"],
             "rho_first": i % 2 == 0}
            for i in range(1, 5)]


def frozen_inputs():
    hashes = dynamic.frozen_inputs()
    if parent.digest(dynamic.SOURCE) != SOURCE_SHA256:
        raise ValueError("executed IC source differs")
    if len(TARGETS.read_text().splitlines()) != 4:
        raise ValueError("four frozen reserve points required")
    return hashes


def make_point_inputs(output):
    directory = output / "inputs"
    directory.mkdir()
    inputs = {}
    for index, line in enumerate(TARGETS.read_text().splitlines(), start=1):
        path = directory / f"point_{index}.jsonl"
        path.write_text(line + "\n")
        inputs[index] = path
    return inputs


def run_ic(output, binary, point_index, point_file, arm, group):
    pipeline = "k1024root22" if arm.startswith("k1024root22_") else arm
    k, layout = PIPELINES[pipeline]
    name = f"p{point_index}_{arm}"
    folder = output / "runs" / name
    command = [str(binary), f"construct:83:1:{k}", str(point_file), "29",
               str(folder / "targets_result.jsonl")]
    print(f"starting {name} at {parent.now()}", flush=True)
    usage = parent.one_run(folder, command, dynamic.environment(folder, layout),
                           IC_CAP_SECONDS, group)
    print(f"{name}: {usage['status']} after {usage['wall_seconds']:.3f}s", flush=True)


def run_rho(output, binary, point_index, point_file, group):
    name = f"p{point_index}_rho"
    folder = output / "runs" / name
    command = [str(binary), str(point_file), str(RHO_SEED_BASE + point_index),
               str(cold.RHO_CONFIG), str(folder / "receipt.json")]
    print(f"starting {name} at {parent.now()}", flush=True)
    usage = parent.one_run(folder, command,
                           {"PATH": os.environ.get("PATH", "/usr/bin:/bin")},
                           RHO_CAP_SECONDS, group)
    print(f"{name}: {usage['status']} after {usage['wall_seconds']:.3f}s", flush=True)


def run(binary, rho_binary, output):
    hashes = frozen_inputs()
    if output.exists():
        raise ValueError("output directory already exists")
    output.mkdir(parents=True)
    group = parent.cgroup_root()
    host = parent.host_receipt(group)
    host["schema"] = "n83.cold-reserve-host/v1"
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
    inputs = make_point_inputs(output)
    parent.write_json(output / "source.json", {
        "schema": "n83.cold-reserve-source/v1", "hashes": hashes,
        "source_sha256": SOURCE_SHA256,
        "binary_sha256": parent.digest(binary),
        "rho_binary_sha256": parent.digest(rho_binary),
        "runner_sha256": parent.digest(Path(__file__)),
        "schedule": schedule(),
        "point_input_sha256": {str(i): parent.digest(path)
                               for i, path in inputs.items()},
        "rho_seed_base": RHO_SEED_BASE,
    })
    for row in schedule():
        i = row["point_index"]
        if row["rho_first"]:
            run_rho(output, rho_binary, i, inputs[i], group)
        for arm in row["ic_order"]:
            run_ic(output, binary, i, inputs[i], arm, group)
        if not row["rho_first"]:
            run_rho(output, rho_binary, i, inputs[i], group)
    return 0


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--preflight", action="store_true")
    parser.add_argument("--binary", type=Path)
    parser.add_argument("--rho-binary", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.preflight:
        print(json.dumps({"hashes": frozen_inputs(), "schedule": schedule()},
                         sort_keys=True))
        return 0
    if not args.binary or not args.rho_binary or not args.output:
        parser.error("--binary, --rho-binary and --output are required")
    return run(args.binary.resolve(), args.rho_binary.resolve(),
               args.output.resolve())


if __name__ == "__main__":
    raise SystemExit(main())
