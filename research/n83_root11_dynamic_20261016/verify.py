#!/usr/bin/env python3
"""Independent base, relation, scalar, and resource replay for one-root layouts."""

import argparse
import importlib.util
import json
import math
import statistics
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RUNNER = HERE / "run.py"
COLD_VERIFIER = ROOT / "research/n83_cold_factor_base_20261016/verify.py"

rspec = importlib.util.spec_from_file_location("n83_dynamic_runner", RUNNER)
runner = importlib.util.module_from_spec(rspec)
rspec.loader.exec_module(runner)
vspec = importlib.util.spec_from_file_location("n83_cold_verifier", COLD_VERIFIER)
coldver = importlib.util.module_from_spec(vspec)
vspec.loader.exec_module(coldver)
group = coldver.group
coldver.runner.parent.EXPECTED[coldver.runner.parent.SOURCE] = runner.SOURCE_SHA256


def one_json(path):
    return json.loads(path.read_text())


def expected_layout(k, layout):
    states = 83 * k * (k - 1) // 2 + 42 * k
    fixed = layout == "root11_fixed"
    one_root = layout != "root22_dynamic"
    request = 1 << 29 if fixed else (4 * states + 2) // 3
    slots = max(16, 1 << (2 * request - 1).bit_length())
    return states, request, slots, one_root, fixed


def verify_ic(output, arm, k, layout, point, scalar):
    folder = output / "runs" / arm
    if not (folder / "usage.json").exists():
        return {"arm": arm, "k": k, "layout": layout, "status": "NOT_RUN"}
    usage = one_json(folder / "usage.json")
    if usage["status"] != "PASS_PRODUCER":
        return {"arm": arm, "k": k, "layout": layout,
                "status": usage["status"], "usage": usage}
    summary = one_json(folder / "stdout.jsonl")
    targets = [json.loads(line) for line in
               (folder / "targets_result.jsonl").read_text().splitlines() if line]
    dump = one_json(folder / "base_dump.jsonl")
    states, request, slots, one_root, fixed = expected_layout(k, layout)
    if (summary["n"] != 83 or summary["a"] != 1
            or summary["targets"] != 1 or summary["targets_solved"] != 1
            or summary["targets_failed"] != 0
            or summary["orbit_columns"] != k or summary["rank"] != k
            or summary["factor_base_points"] != 166 * k
            or summary["regular_states"] != states
            or summary["root_state_layout"] != (
                "root_only83_11byte_states_reconstructed_second_root" if one_root
                else "root_only83_22byte_states")
            or summary["root_state_slot_bytes"] != (11 if one_root else 22)
            or summary["root_table_layout"] != "indexed83_5byte_exact"
            or summary["root_table_slot_bytes"] != 5
            or summary["root_table_request_entries"] != request
            or summary["root_table_allocated_bytes"] != 5 * slots
            or summary["root_table_entries"] >= slots
            or summary["root_table_indexed5_enabled"] is not True
            or summary["root_table_indexed4_enabled"] is not False
            or summary["indexed5_load75_enabled"] is fixed
            or summary["indexed5_fixed30_enabled"] is not fixed
            or summary["root_only11_enabled"] is not one_root
            or summary["symmetric_s3_states_enabled"] is not True
            or summary["sparse_canonical_enabled"] is not True
            or summary["root_bloom_enabled"] is not False
            or (summary["root_bloom_profile"], summary["root_bloom_hashes"],
                summary["root_bloom_bits_per_request"], summary["root_bloom_bytes"])
               != ("none", 0, 0, 0)
            or len(targets) != 1):
        raise ValueError(f"{arm}: source layout, base, rank or target count differs")
    base = group.replay_base(dump, summary)
    if len(base) != 166 * k:
        raise ValueError(f"{arm}: independent factor-base replay failed")
    target = targets[0]
    indices = target["point_indices"]
    recovered = int(target["recovered_scalar"])
    if (target["fixture_index"] != 0
            or group.point_from_json(target["target"]) != point
            or recovered != scalar or target["group_verified"] is not True
            or group.multiply(group.G, recovered) != point
            or len(indices) != 4
            or any(not isinstance(index, int) or index < 0 or index >= len(base)
                   for index in indices)):
        raise ValueError(f"{arm}: scalar or factor-base witness differs")
    relation = None
    for index in indices:
        relation = group.add(relation, base[index])
    if relation != point:
        raise ValueError(f"{arm}: four-point relation replay failed")
    phase_keys = ("target_query_ms", "target_pdp_ms", "target_relation_check_ms",
                  "target_descent_ms", "target_recovery_check_ms")
    online_ms = float(target["target_ms"])
    phases = {key: float(target[key]) for key in phase_keys}
    if (not math.isfinite(online_ms) or online_ms <= 0
            or any(not math.isfinite(value) or value < 0 for value in phases.values())
            or abs(online_ms - sum(phases.values())) > 1.0):
        raise ValueError(f"{arm}: exclusive online phases differ")
    timing = summary["timing_ms"]
    cold_ms = float(timing["process_total"])
    setup_ms = float(timing["setup_total"])
    rank_ms = float(timing["rank_stage"])
    targets_ms = float(timing["targets_total"])
    residual = cold_ms - setup_ms - rank_ms - targets_ms
    if (not all(math.isfinite(x) and x >= 0 for x in
                (cold_ms, setup_ms, rank_ms, targets_ms))
            or not 0 <= residual <= 1_000
            or not 0 <= targets_ms - online_ms <= 1_000):
        raise ValueError(f"{arm}: cold process interval or phase accounting differs")
    return {"arm": arm, "k": k, "layout": layout,
            "status": "PASS_INDEPENDENT_REPLAY",
            "actual_factor_base_points": len(base), "effective_columns": k,
            "base_hash": summary["base_hash"],
            "root_table_entries": summary["root_table_entries"],
            "root_table_allocated_bytes": summary["root_table_allocated_bytes"],
            "root_state_allocated_bytes": summary["root_state_allocated_bytes"],
            "rank_attempts": summary["rank_attempts"],
            "rank_relations": summary["rank_relations"],
            "cold_process_ms_exploratory": cold_ms,
            "setup_ms": setup_ms, "rank_stage_ms": rank_ms,
            "targets_total_ms": targets_ms, "cold_unattributed_ms": residual,
            "one_target_online_ms": online_ms,
            "online_phases_ms": phases, "target_probes": target["probes"],
            "full_child_wall_seconds": usage["wall_seconds"],
            "resources": coldver.resource(usage)}


def repeatability(rows, a_name, b_name):
    a = rows.get(a_name)
    b = rows.get(b_name)
    if a is None or b is None:
        return None
    if (a["base_hash"], a["root_table_entries"]) != (
            b["base_hash"], b["root_table_entries"]):
        return {"status": "BASE_OR_INDEX_MISMATCH"}
    times = (a["cold_process_ms_exploratory"],
             b["cold_process_ms_exploratory"])
    delta = abs(times[0] - times[1]) / statistics.mean(times)
    return {"cold_process_ms": list(times), "relative_difference": delta,
            "status": "PASS_10_PERCENT" if delta <= 0.10 else "NOISY"}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    output = args.output.resolve()
    hashes, point, scalar = coldver.frozen_point()
    if hashes != runner.frozen_inputs():
        raise ValueError("frozen source or target changed")
    host = one_json(output / "host.json")
    source = one_json(output / "source.json")
    if (host["schema"] != "n83.root11-dynamic-host/v1"
            or source["schema"] != "n83.root11-dynamic-source/v1"
            or source["source_sha256"] != runner.SOURCE_SHA256
            or source["hashes"] != hashes
            or source["runner_sha256"] != runner.parent.digest(RUNNER)
            or source["arm_order"] != [arm for arm, _, _, _ in runner.ARMS]
            or source["rho_seed"] != runner.cold.RHO_SEED
            or host["git_sha"] is None):
        raise ValueError("executed source, host, or arm order differs")
    arms = []
    for arm, k, layout, _ in runner.ARMS:
        try:
            arms.append(verify_ic(output, arm, k, layout, point, scalar))
        except Exception as error:
            arms.append({"arm": arm, "k": k, "layout": layout,
                         "status": "VERIFICATION_FAILED", "error": str(error)})
    try:
        rho = coldver.verify_rho(output, point, scalar)
    except Exception as error:
        rho = {"status": "VERIFICATION_FAILED", "error": str(error)}
    accepted = [row for row in arms if row["status"] == "PASS_INDEPENDENT_REPLAY"
                and row["resources"]["status"] == "PASS"]
    by_arm = {row["arm"]: row for row in accepted}
    fixed_aa = repeatability(by_arm, "k1271fixed_a", "k1271fixed_b")
    dynamic_aa = repeatability(by_arm, "k1271dynamic_a", "k1271dynamic_b")
    complete = (len(accepted) == len(runner.ARMS)
                and rho["status"] == "PASS_INDEPENDENT_REPLAY"
                and rho["resources"]["status"] == "PASS"
                and fixed_aa is not None and fixed_aa["status"] == "PASS_10_PERCENT"
                and dynamic_aa is not None and dynamic_aa["status"] == "PASS_10_PERCENT")
    cold_by_pipeline = {}
    if len(accepted) == len(runner.ARMS):
        for row in accepted:
            key = f"k{row['k']}{row['layout']}"
            cold_by_pipeline.setdefault(key, []).append(
                row["cold_process_ms_exploratory"])
        cold_by_pipeline = {key: statistics.mean(times)
                            for key, times in cold_by_pipeline.items()}
    analysis = {
        "schema": "n83.root11-dynamic-analysis/v1",
        "hardware_tier": "cloud_instance_exploratory",
        "workload_role": "one_fresh_public_point_root_layout_screen",
        "target_count": 1, "arm_order": source["arm_order"],
        "status": "SCREEN_COMPLETE" if complete else "BOUNDED_UNKNOWN",
        "selected_factor_base": None,
        "descriptive_cold_process_ms_by_pipeline": cold_by_pipeline,
        "fixed_aa_control": fixed_aa, "dynamic_aa_control": dynamic_aa,
        "arms": arms, "rho": rho, "controlled_online_speedup": None,
    }
    (output / "analysis.json").write_text(json.dumps(analysis, indent=2,
                                                    sort_keys=True) + "\n")
    print(json.dumps(analysis, sort_keys=True))
    return 0 if complete else 1


if __name__ == "__main__":
    raise SystemExit(main())
