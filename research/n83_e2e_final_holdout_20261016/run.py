#!/usr/bin/env python3
"""Run one eight-point stratum of the frozen N83 finalist comparison."""

import argparse
import importlib.util
import json
import os
import shutil
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
TARGETS_SHA256 = "407cd248797100214ff039193bb04bf89b5c2708a8cd339ca98f7aca98174062"
ANSWERS_SHA256 = "fe771284a3b0e5c0301309f678751f0630f0fa23844dd0ce1dbc9fa47f243076"
SOURCE_SHA256 = "8a2fe55dd90dc33892c197d27260182bf783fc8212fb485d8297daf2bbba8e12"
ARMS = (
    ("k1104root11", 1104, "root11_dynamic", "IC1N83Ckb1a1fb183264PDP4rootRCguidedLAgaussTDdirectISO0h6f30a77abf9e"),
    ("k1104root22", 1104, "root22_dynamic", "IC1N83Ckb1a1fb183264PDP4rootRCguidedLAgaussTDdirectISO0hce7cffedf9c0"),
    ("k1107root11", 1107, "root11_dynamic", "IC1N83Ckb1a1fb183762PDP4rootRCguidedLAgaussTDdirectISO0h8c736ae9c7b4"),
    ("k1107root22", 1107, "root22_dynamic", "IC1N83Ckb1a1fb183762PDP4rootRCguidedLAgaussTDdirectISO0hcb73613ebc5a"),
    ("k1112root11", 1112, "root11_dynamic", "IC1N83Ckb1a1fb184592PDP4rootRCguidedLAgaussTDdirectISO0h49683406b22f"),
    ("k1112root22", 1112, "root22_dynamic", "IC1N83Ckb1a1fb184592PDP4rootRCguidedLAgaussTDdirectISO0hbc5ec3072dc4"),
    ("k1120root11", 1120, "root11_dynamic", "IC1N83Ckb1a1fb185920PDP4rootRCguidedLAgaussTDdirectISO0hc922b499a234"),
    ("k1120root22", 1120, "root22_dynamic", "IC1N83Ckb1a1fb185920PDP4rootRCguidedLAgaussTDdirectISO0h82b74c23b7bc"),
)
IC_CAP_SECONDS = 480
RHO_CAP_SECONDS = 120
GLOBAL_PRODUCER_CAP_SECONDS = 190 * 60
CONTROL_KEYS = ("k1107root22", "k1104root11")
PAIR_RHO_SEED_BASE = 2026101700


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    value = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(value)
    return value


finegrid = module("n83_finalist_finegrid_run",
                  ROOT / "research/n83_e2e_finegrid_20261016/run.py")
holdout = finegrid.holdout
dynamic = finegrid.dynamic
cold = finegrid.cold
parent = finegrid.parent
CALIBRATION_RHO_SEED = int(cold.RHO_SEED)


def fixture(fixture_dir):
    hashes = finegrid.frozen_inputs()
    targets = fixture_dir / "targets.jsonl"
    answers = fixture_dir / "answers.jsonl"
    manifest_path = fixture_dir / "fixture_manifest.json"
    if (parent.digest(targets) != TARGETS_SHA256
            or parent.digest(answers) != ANSWERS_SHA256):
        raise ValueError("sealed finalist target or answer file differs")
    manifest = json.loads(manifest_path.read_text())
    if (manifest.get("schema") != "n83.e2e-final-holdout-fixture/v1"
            or manifest.get("targets_sha256") != TARGETS_SHA256
            or manifest.get("answers_sha256") != ANSWERS_SHA256
            or manifest.get("target_count") != 16
            or manifest.get("domain") != "n83-e2e-final-holdout-sixteen-20261016-v1"):
        raise ValueError("finalist fixture manifest differs")
    lines = targets.read_text().splitlines(keepends=True)
    if len(lines) != 16 or any(not line.endswith("\n") for line in lines):
        raise ValueError("finalist public-point count or encoding differs")
    return hashes, lines


def schedule(half):
    if half not in (0, 1):
        raise ValueError("half must be zero or one")
    events = []
    for local in range(8):
        index = half * 8 + local
        rotated = ARMS[local:] + ARMS[:local]
        for position, arm in enumerate(rotated):
            name = arm[0]
            events.append({"point_index": index, "name": f"cal_before_{name}",
                           "candidate": name, "kind": "calibration"})
            fixed_index = next(i for i, item in enumerate(ARMS) if item[0] == name)
            pair = ("ic", "rho") if (local + fixed_index) % 2 == 0 else ("rho", "ic")
            for kind in pair:
                events.append({"point_index": index,
                               "name": name if kind == "ic" else f"rho_{name}",
                               "candidate": name, "kind": kind})
            events.append({"point_index": index, "name": f"cal_after_{name}",
                           "candidate": name, "kind": "calibration"})
    events.extend({"point_index": half * 8, "name": f"control_{name}_b",
                   "candidate": name, "kind": "control"}
                  for name in CONTROL_KEYS)
    return events


def run(binary, rho_binary, output, half, fixture_dir):
    hashes, lines = fixture(fixture_dir)
    if output.exists():
        raise ValueError("output directory already exists")
    output.mkdir(parents=True)
    events = schedule(half)
    group = parent.cgroup_root()
    host = parent.host_receipt(group)
    host["schema"] = "n83.e2e-finalist-host/v1"
    host["parent_instance_id"] = host["instance_id"]
    host["instance_id"] = f"{host['instance_id']}-half{half}"
    host["stratum_half"] = half
    host["topology"] = holdout.topology_receipt(host["visible_affinity"] or [])
    host["cpu_telemetry_before"] = holdout.cpu_telemetry()
    parent.write_json(output / "host.json", host)
    memory_max = (host["cgroup_before"]["memory.max"] or "").strip()
    memory_too_small = memory_max.isdigit() and int(memory_max) < cold.RSS_CAP_BYTES
    if (host["mem_total_kib"] is None
            or host["mem_total_kib"] * 1024 < cold.RSS_CAP_BYTES + 256 * 1024**2
            or host["disk_free_bytes"] < 3 * 1024**3
            or len(host["visible_affinity"] or []) < 4 or memory_too_small):
        parent.write_json(output / "preflight.json", {
            "status": "REJECTED_HOST_CAPACITY", "host": host})
        return 2
    parent.write_json(output / "preflight.json", {
        "status": "PASS_CAPACITY_PREFLIGHT", "rss_cap_bytes": cold.RSS_CAP_BYTES})
    executed_bins = output / "executed_bins"
    executed_bins.mkdir()
    binary_copy = executed_bins / "koblitz_orbit_dlp_fast"
    rho_copy = executed_bins / "n83_rho_optimization_eval"
    shutil.copy2(binary, binary_copy)
    shutil.copy2(rho_binary, rho_copy)
    binary_copy.chmod(0o755)
    rho_copy.chmod(0o755)
    parent.write_json(output / "source.json", {
        "schema": "n83.e2e-finalist-source/v1", "hashes": hashes,
        "source_sha256": SOURCE_SHA256,
        "binary_sha256": parent.digest(binary_copy),
        "rho_binary_sha256": parent.digest(rho_copy),
        "runner_sha256": parent.digest(Path(__file__)),
        "targets_sha256": TARGETS_SHA256,
        "answers_sha256": ANSWERS_SHA256,
        "fixture_manifest_sha256": parent.digest(fixture_dir / "fixture_manifest.json"),
        "half": half, "arms": ARMS, "events": events,
        "rho_seed_rule": "2026101700+8*point_index+fixed_candidate_index",
        "calibration_seed": CALIBRATION_RHO_SEED,
        "calibration_target_sha256": parent.digest(cold.TARGETS),
        "calibration_answer_sha256": parent.digest(cold.ANSWERS),
        "global_producer_cap_seconds": GLOBAL_PRODUCER_CAP_SECONDS,
    })
    by_name = {arm[0]: arm for arm in ARMS}
    started = time.monotonic()
    for event in events:
        index, name, candidate, kind = (event[key] for key in
                                         ("point_index", "name", "candidate", "kind"))
        point_output = output / "points" / f"p{index:02d}"
        point_output.mkdir(parents=True, exist_ok=True)
        target = point_output / "target.jsonl"
        if target.exists():
            if target.read_text() != lines[index]:
                raise ValueError("per-point target file differs")
        else:
            target.write_text(lines[index])
        folder = point_output / "runs" / name
        if kind in ("rho", "calibration"):
            rho_target = cold.TARGETS if kind == "calibration" else target
            seed = (CALIBRATION_RHO_SEED if kind == "calibration" else
                    PAIR_RHO_SEED_BASE + 8 * index +
                    next(i for i, arm in enumerate(ARMS) if arm[0] == candidate))
            command = [str(rho_copy), str(rho_target),
                       str(seed), str(cold.RHO_CONFIG),
                       str(folder / "receipt.json")]
            environment = {"PATH": os.environ.get("PATH", "/usr/bin:/bin")}
            cap = RHO_CAP_SECONDS
        else:
            _, k, layout, _ = by_name[candidate]
            command = [str(binary_copy), f"construct:83:1:{k}", str(target),
                       "29", str(folder / "targets_result.jsonl")]
            environment = dynamic.environment(folder, layout)
            cap = IC_CAP_SECONDS
        elapsed = time.monotonic() - started
        if elapsed + cap + 30 > GLOBAL_PRODUCER_CAP_SECONDS:
            parent.write_json(output / "producer_stop.json", {
                "status": "GLOBAL_CAP_BEFORE_NEXT_ARM",
                "elapsed_seconds": elapsed,
                "next_event": event,
                "global_producer_cap_seconds": GLOBAL_PRODUCER_CAP_SECONDS,
            })
            break
        print(f"starting p{index:02d}/{name} at {parent.now()}", flush=True)
        usage = holdout.run_with_telemetry(folder, command, environment,
                                          cap, group)
        print(f"p{index:02d}/{name}: {usage['status']} after "
              f"{usage['wall_seconds']:.3f}s", flush=True)
    parent.write_json(output / "host_telemetry_after.json",
                      holdout.cpu_telemetry())
    return 0


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--preflight", action="store_true")
    parser.add_argument("--half", type=int, choices=(0, 1), required=True)
    parser.add_argument("--fixture-dir", type=Path, default=HERE)
    parser.add_argument("--binary", type=Path)
    parser.add_argument("--rho-binary", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    fixture_dir = args.fixture_dir.resolve()
    if args.preflight:
        hashes, _ = fixture(fixture_dir)
        print(json.dumps({"hashes": hashes, "half": args.half,
                          "arms": ARMS, "events": schedule(args.half),
                          "targets_sha256": TARGETS_SHA256}, sort_keys=True))
        return 0
    if not args.binary or not args.rho_binary or not args.output:
        parser.error("--binary, --rho-binary and --output are required")
    return run(args.binary.resolve(), args.rho_binary.resolve(),
               args.output.resolve(), args.half, fixture_dir)


if __name__ == "__main__":
    raise SystemExit(main())
