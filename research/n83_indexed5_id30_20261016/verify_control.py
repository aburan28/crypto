#!/usr/bin/env python3
"""Independent N83 K1024 exact-base, relation and scalar replay for ID30."""

import argparse
import importlib.util
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RUNNER = HERE / "run_control.py"
DYNAMIC_VERIFIER = ROOT / "research/n83_root11_dynamic_20261016/verify.py"


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    value = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(value)
    return value


runner = module("n83_id30_control_runner", RUNNER)
dynamic = module("n83_id30_control_dynamic_verifier", DYNAMIC_VERIFIER)
group = dynamic.group


def read(path):
    return json.loads(path.read_text())


def verify_arm(output, name, point, scalar):
    folder = output / "runs" / name
    usage = read(folder / "usage.json")
    if usage["status"] != "PASS_PRODUCER":
        return {"arm": name, "status": usage["status"], "usage": usage}
    summary = read(folder / "stdout.jsonl")
    dump = read(folder / "base_dump.jsonl")
    targets = [json.loads(line) for line in (
        folder / "targets_result.jsonl").read_text().splitlines() if line]
    id30 = name == "opt_in_id30"
    states = 83 * runner.K * (runner.K - 1) // 2 + 42 * runner.K
    request = (4 * states + 2) // 3
    slots = max(16, 1 << (2 * request - 1).bit_length())
    if (summary["n"] != 83 or summary["a"] != 1
            or summary["targets"] != 1 or summary["targets_solved"] != 1
            or summary["targets_failed"] != 0
            or summary["orbit_columns"] != runner.K
            or summary["rank"] != runner.K
            or summary["factor_base_points"] != 166 * runner.K
            or summary["regular_states"] != states
            or summary["root_state_slot_bytes"] != 11
            or summary["root_table_layout"] != (
                "indexed83_5byte_exact_id30" if id30 else "indexed83_5byte_exact")
            or summary["root_table_slot_bytes"] != 5
            or summary["root_table_request_entries"] != request
            or summary["root_table_allocated_bytes"] != 5 * slots
            or summary["root_table_entries"] >= slots
            or summary["root_table_indexed5_enabled"] is not True
            or summary["root_table_indexed5_id30_enabled"] is not id30
            or summary["indexed5_load75_enabled"] is not True
            or summary["indexed5_fixed30_enabled"] is not False
            or summary["root_only11_enabled"] is not True
            or summary["root_bloom_enabled"] is not False
            or summary["sparse_canonical_enabled"] is not True
            or len(targets) != 1):
        raise ValueError(f"{name}: exact ID, layout, base or target receipt differs")
    base = group.replay_base(dump, summary)
    target = targets[0]
    indices = target["point_indices"]
    recovered = int(target["recovered_scalar"])
    if (len(base) != 166 * runner.K or target["fixture_index"] != 0
            or group.point_from_json(target["target"]) != point
            or recovered != scalar or target["group_verified"] is not True
            or group.multiply(group.G, recovered) != point
            or len(indices) != 4
            or any(not isinstance(index, int) or index < 0 or index >= len(base)
                   for index in indices)):
        raise ValueError(f"{name}: independent factor-base or scalar replay failed")
    relation = None
    for index in indices:
        relation = group.add(relation, base[index])
    if relation != point:
        raise ValueError(f"{name}: four-point target relation differs")
    phase_keys = ("target_query_ms", "target_pdp_ms",
                  "target_relation_check_ms", "target_descent_ms",
                  "target_recovery_check_ms")
    online = float(target["target_ms"])
    phases = {key: float(target[key]) for key in phase_keys}
    if (not math.isfinite(online) or online <= 0
            or any(not math.isfinite(value) or value < 0 for value in phases.values())
            or abs(online - sum(phases.values())) > 1.0):
        raise ValueError(f"{name}: online phase accounting differs")
    return {"arm": name, "status": "PASS_INDEPENDENT_REPLAY",
            "base_hash": summary["base_hash"],
            "root_table_entries": summary["root_table_entries"],
            "recovered_scalar": recovered, "point_indices": indices,
            "actual_factor_base_points": len(base), "online_ms_exploratory": online,
            "online_phases_ms": phases, "cold_wall_seconds_exploratory":
                usage["wall_seconds"]}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    output = args.output.resolve()
    source = read(output / "source.json")
    dynamic_runner = runner.dynamic_runner()
    if (source["schema"] != "n83.id30-control-source/v1"
            or source["source_sha256"] != runner.SOURCE_SHA256
            or runner.digest(runner.SOURCE) != runner.SOURCE_SHA256
            or source["runner_sha256"] != runner.digest(RUNNER)
            or source["binary_sha256"] != runner.digest(
                output / "koblitz_orbit_dlp_fast")
            or source["public_input_sha256"] != runner.digest(
                dynamic_runner.cold.TARGETS)
            or source["hashes"] != dynamic_runner.frozen_inputs()
            or tuple(source["arm_order"]) != runner.ARMS
            or source["k"] != runner.K):
        raise ValueError("source, binary or public input differs")
    public = read(dynamic_runner.cold.TARGETS)
    answer = read(dynamic_runner.cold.ANSWERS)
    point = group.point_from_json(public)
    scalar = int(answer["scalar"])
    if (answer["point"] != public or not group.on_curve(point)
            or group.multiply(point, group.R) is not None
            or group.multiply(group.G, scalar) != point):
        raise ValueError("frozen point or answer fails subgroup replay")
    arms = []
    for name in runner.ARMS:
        try:
            arms.append(verify_arm(output, name, point, scalar))
        except Exception as error:
            arms.append({"arm": name, "status": "VERIFICATION_FAILED",
                         "error": str(error)})
    accepted = all(row["status"] == "PASS_INDEPENDENT_REPLAY" for row in arms)
    same = (accepted
            and arms[0]["base_hash"] == arms[1]["base_hash"]
            and arms[0]["root_table_entries"] == arms[1]["root_table_entries"]
            and arms[0]["recovered_scalar"] == arms[1]["recovered_scalar"])
    analysis = {"schema": "n83.id30-control-analysis/v1",
                "hardware_tier": "local_correctness_control",
                "status": "PASS_EXACT_PARITY" if same else "BOUNDED_UNKNOWN",
                "selected_factor_base": None, "arms": arms,
                "same_base_index_and_scalar": same,
                "performance_claim": None}
    (output / "analysis.json").write_text(json.dumps(
        analysis, sort_keys=True, indent=2) + "\n")
    print(json.dumps(analysis, sort_keys=True))
    return 0 if same else 1


if __name__ == "__main__":
    raise SystemExit(main())
