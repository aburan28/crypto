#!/usr/bin/env python3
"""Independent base, target, resource and same-point rho replay for the cold screen."""

import argparse
import importlib.util
import json
import math
import statistics
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RUNNER = HERE / "run.py"
GROUP = ROOT / "research/n83_indexed5_20261013/analyze.py"
EXPECTED_BASE = {
    (2543, "root22"): ("eb0a315471ae0d01598fc24d10b8dae41873cc3fab556fc91f85a584880ccf58",
                       536737856),
    (3450, "root11"): ("ae75c86565f5fa12fb6f30a5539f6f4638bcf1ca8b347a5cefa20f5073a9f572",
                       987887478),
}

rspec = importlib.util.spec_from_file_location("n83_cold_runner", RUNNER)
runner = importlib.util.module_from_spec(rspec)
rspec.loader.exec_module(runner)
gspec = importlib.util.spec_from_file_location("n83_cold_group", GROUP)
group = importlib.util.module_from_spec(gspec)
gspec.loader.exec_module(group)


def one_json(path):
    return json.loads(path.read_text())


def counter(snapshot, section, key):
    for line in (snapshot.get(section) or "").splitlines():
        fields = line.split()
        if len(fields) == 2 and fields[0] == key:
            return int(fields[1])
    return None


def counter_delta(usage, section, key):
    before = counter(usage["cgroup_before"], section, key)
    after = counter(usage["cgroup_after"], section, key)
    return after - before if before is not None and after is not None else None


def resource(usage):
    throttled = counter_delta(usage, "cpu.stat", "nr_throttled")
    oom = counter_delta(usage, "memory.events", "oom")
    oom_kills = counter_delta(usage, "memory.events", "oom_kill")
    return {"peak_rss_bytes": usage["peak_rss_bytes"],
            "major_faults": usage["major_faults"],
            "cpu_throttled_periods_delta": throttled,
            "oom_events_delta": oom, "oom_kills_delta": oom_kills,
            "status": ("PASS" if usage["peak_rss_bytes"] <= runner.RSS_CAP_BYTES
                       and (throttled, oom, oom_kills) == (0, 0, 0)
                       else "NOISY_OR_OVER_CAP")}


def frozen_point():
    hashes = runner.frozen_inputs()
    point_rows = [json.loads(line) for line in runner.TARGETS.read_text().splitlines() if line]
    answer_rows = [json.loads(line) for line in runner.ANSWERS.read_text().splitlines() if line]
    if len(point_rows) != 1 or len(answer_rows) != 1:
        raise ValueError("one primary public point and answer required")
    point = group.point_from_json(point_rows[0])
    answer = answer_rows[0]
    scalar = int(answer["scalar"])
    if (answer["fixture_index"] != 0
            or group.point_from_json(answer["point"]) != point
            or not group.on_curve(point)
            or group.multiply(point, group.R) is not None
            or group.multiply(group.G, scalar) != point):
        raise ValueError("frozen one-target scalar/subgroup replay failed")
    return hashes, point, scalar


def verify_ic(output, arm, k, variant, profile, point, scalar):
    folder = output / "runs" / arm
    if not (folder / "usage.json").exists():
        return {"arm": arm, "status": "NOT_RUN"}
    usage = one_json(folder / "usage.json")
    if usage["status"] != "PASS_PRODUCER":
        return {"arm": arm, "k": k, "variant": variant,
                "status": usage["status"], "usage": usage}
    summary = one_json(folder / "stdout.jsonl")
    targets = [json.loads(line) for line in
               (folder / "targets_result.jsonl").read_text().splitlines() if line]
    dump = one_json(folder / "base_dump.jsonl")
    root22 = variant == "root22"
    state_bound = 83 * k * (k - 1) // 2 + 42 * k
    request = (4 * state_bound + 2) // 3 if root22 else 1 << 29
    slots = max(16, 1 << (2 * request - 1).bit_length())
    bloom = ("none", 0, 0, 0) if root22 else ("wide3", 3, 16, 1 << 30)
    if (summary["n"] != 83 or summary["a"] != 1
            or summary["targets"] != 1 or summary["targets_solved"] != 1
            or summary["targets_failed"] != 0
            or summary["orbit_columns"] != k or summary["rank"] != k
            or summary["factor_base_points"] != 166 * k
            or summary["regular_states"] != state_bound
            or summary["root_table_layout"] != "indexed83_5byte_exact"
            or summary["root_table_slot_bytes"] != 5
            or summary["root_table_indexed5_enabled"] is not True
            or summary["root_table_indexed4_enabled"] is not False
            or summary["root_state_layout"] != (
                "root_only83_22byte_states" if root22 else
                "root_only83_11byte_states_reconstructed_second_root")
            or summary["root_state_slot_bytes"] != (22 if root22 else 11)
            or summary["root_table_request_entries"] != request
            or summary["root_table_allocated_bytes"] != 5 * slots
            or summary["root_table_entries"] >= slots
            or summary["indexed5_load75_enabled"] is not root22
            or summary["indexed5_fixed30_enabled"] is root22
            or summary["root_only11_enabled"] is root22
            or summary["symmetric_s3_states_enabled"] is not True
            or summary["sparse_canonical_enabled"] is not True
            or summary["root_bloom_enabled"] is root22
            or (summary["root_bloom_profile"], summary["root_bloom_hashes"],
                summary["root_bloom_bits_per_request"], summary["root_bloom_bytes"])
               != bloom or len(targets) != 1):
        raise ValueError(f"{arm}: source layout, base, rank or target count differs")
    pinned = EXPECTED_BASE.get((k, variant))
    if pinned and (summary["base_hash"], summary["root_table_entries"]) != pinned:
        raise ValueError(f"{arm}: admitted base or exact index changed")
    base = group.replay_base(dump, summary)
    if len(base) != 166 * k:
        raise ValueError(f"{arm}: independent factor-base replay failed")
    target = targets[0]
    recovered = int(target["recovered_scalar"])
    indices = target["point_indices"]
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
        raise ValueError(f"{arm}: exclusive one-target online phases differ")
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
    resources = resource(usage)
    return {"arm": arm, "k": k, "variant": variant, "bloom_profile": profile,
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
            "resources": resources}


def verify_rho(output, point, scalar):
    folder = output / "runs" / "rho"
    if not (folder / "usage.json").exists():
        return {"status": "NOT_RUN"}
    usage = one_json(folder / "usage.json")
    if usage["status"] != "PASS_PRODUCER":
        return {"status": usage["status"], "usage": usage}
    receipt = one_json(folder / "receipt.json")
    recovered = int(receipt["recovered_scalar"])
    online_ms = float(receipt["online_ms"])
    if (receipt["status"] != "PASS_VERIFIED_PUBLIC_TARGET"
            or group.point_from_json(receipt["public_target_q"]) != point
            or recovered != scalar
            or receipt["generic_scalar_replay_verified"] is not True
            or group.multiply(group.G, recovered) != point
            or not math.isfinite(online_ms) or online_ms <= 0):
        raise ValueError("same-point rho scalar replay failed")
    return {"status": "PASS_INDEPENDENT_REPLAY", "online_ms": online_ms,
            "resources": resource(usage),
            "full_child_wall_seconds": usage["wall_seconds"]}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    output = args.output.resolve()
    hashes, point, scalar = frozen_point()
    host = one_json(output / "host.json")
    source = one_json(output / "source.json")
    if (host["schema"] != "n83.cold-factor-base-host/v1"
            or source["schema"] != "n83.cold-factor-base-source/v1"
            or source["hashes"] != hashes
            or source["runner_sha256"] != runner.parent.digest(RUNNER)
            or source["arm_order"] != [arm for arm, _, _, _, _ in runner.ARMS]
            or source["rho_seed"] != runner.RHO_SEED
            or host["git_sha"] is None):
        raise ValueError("executed source, host, or arm order differs")
    arms = []
    for arm, k, variant, profile, _ in runner.ARMS:
        try:
            arms.append(verify_ic(output, arm, k, variant, profile, point, scalar))
        except Exception as error:
            arms.append({"arm": arm, "k": k, "variant": variant,
                         "status": "VERIFICATION_FAILED", "error": str(error)})
    try:
        rho = verify_rho(output, point, scalar)
    except Exception as error:
        rho = {"status": "VERIFICATION_FAILED", "error": str(error)}
    accepted = [row for row in arms if row["status"] == "PASS_INDEPENDENT_REPLAY"
                and row["resources"]["status"] == "PASS"]
    by_arm = {row["arm"]: row for row in accepted}
    aa = None
    if "k2543root22_a" in by_arm and "k2543root22_b" in by_arm:
        a = by_arm["k2543root22_a"]
        b = by_arm["k2543root22_b"]
        if (a["base_hash"], a["root_table_entries"]) != (
                b["base_hash"], b["root_table_entries"]):
            a["status"] = b["status"] = "REPEATED_BASE_OR_INDEX_MISMATCH"
        else:
            times = (a["cold_process_ms_exploratory"], b["cold_process_ms_exploratory"])
            aa = {"cold_process_ms": list(times),
                  "relative_difference": abs(times[0] - times[1]) / statistics.mean(times),
                  "status": ("PASS_10_PERCENT" if
                             abs(times[0] - times[1]) / statistics.mean(times) <= 0.10
                             else "NOISY")}
    complete = (len(accepted) == len(runner.ARMS)
                and all(row["status"] == "PASS_INDEPENDENT_REPLAY" for row in accepted)
                and rho["status"] == "PASS_INDEPENDENT_REPLAY"
                and rho["resources"]["status"] == "PASS"
                and aa is not None and aa["status"] == "PASS_10_PERCENT")
    evidence_valid = (len(accepted) > 0
                      and all(row["status"] in ("PASS_INDEPENDENT_REPLAY", "TIMEOUT")
                              for row in arms)
                      and rho["status"] == "PASS_INDEPENDENT_REPLAY")
    cold_by_pipeline = {}
    if len(accepted) == len(runner.ARMS):
        for row in accepted:
            key = f"k{row['k']}{row['variant']}" + (
                "wide3" if row["bloom_profile"] == "wide3" else "")
            cold_by_pipeline.setdefault(key, []).append(row["cold_process_ms_exploratory"])
        cold_by_pipeline = {key: statistics.mean(times)
                            for key, times in cold_by_pipeline.items()}
    ranked = sorted(cold_by_pipeline, key=cold_by_pipeline.get)
    analysis = {
        "schema": "n83.cold-factor-base-analysis/v1",
        "hardware_tier": "cloud_instance_exploratory",
        "workload_role": "one_fresh_public_point_cold_factor_base_screen",
        "target_count": 1, "arm_order": source["arm_order"],
        "status": "SCREEN_COMPLETE" if complete else "BOUNDED_UNKNOWN",
        "evidence_status": ("PASS_INDEPENDENT_REPLAY" if evidence_valid
                            else "VERIFICATION_FAILED"),
        "selected_factor_base": None,
        "descriptive_cold_process_ms_by_pipeline": cold_by_pipeline,
        "confirmation_candidates": ranked[:2] if complete else [],
        "aa_control": aa, "arms": arms, "rho": rho,
        "controlled_online_speedup": None,
    }
    (output / "analysis.json").write_text(json.dumps(analysis, indent=2,
                                                    sort_keys=True) + "\n")
    print(json.dumps(analysis, sort_keys=True))
    return 0 if evidence_valid else 1


if __name__ == "__main__":
    raise SystemExit(main())
