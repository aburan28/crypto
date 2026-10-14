#!/usr/bin/env python3
"""Independent five-pair N83 Tier 1 correctness and noise gate."""

import argparse
import importlib.util
import json
import math
import random
import statistics
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
ANALYZER = ROOT / "research/n83_indexed5_20261013/analyze.py"
RUNNER = HERE / "run.py"
CALIBRATION_SHA256 = (
    "27f6d3422dbe05f9c859931528c860b2005763e60fc43d61f9f9e27cd9ce73d7",
    "7df4ffe1d660e67ab1ddbfdf0fdb2d29d6158565da1a25b3fd6d8fff118b00e8",
    "4aac5ed01887596307333d7ff92134cf542cdd7307d3ac4ceab6c79968b28824",
)

spec = importlib.util.spec_from_file_location("n83_final_group", ANALYZER)
group = importlib.util.module_from_spec(spec)
spec.loader.exec_module(group)
runner_spec = importlib.util.spec_from_file_location("n83_final_runner", RUNNER)
runner = importlib.util.module_from_spec(runner_spec)
runner_spec.loader.exec_module(runner)


def one_json(path):
    return json.loads(path.read_text())


def json_lines(path):
    return [json.loads(line) for line in path.read_text().splitlines() if line]


def counters(text):
    result = {}
    for line in (text or "").splitlines():
        fields = line.split()
        if len(fields) == 2 and fields[1].isdigit():
            result[fields[0]] = int(fields[1])
    return result


def counter_delta(before, after, key):
    a, b = counters(before).get(key), counters(after).get(key)
    return None if a is None or b is None or b < a else b - a


def proc_steal_ticks(text, cpus):
    rows = {}
    for line in (text or "").splitlines():
        fields = line.split()
        if fields and fields[0].startswith("cpu") and fields[0][3:].isdigit():
            if len(fields) >= 9:
                rows[int(fields[0][3:])] = int(fields[8])
    return sum(rows[cpu] for cpu in cpus) if all(cpu in rows for cpu in cpus) else None


def fresh_points():
    runner.frozen_inputs()
    points = [group.point_from_json(row) for row in json_lines(runner.TARGETS)]
    answers = json_lines(runner.ANSWERS)
    if len(points) != 5 or len(answers) != 5:
        raise ValueError("fresh target panel is incomplete")
    for i, (point, answer) in enumerate(zip(points, answers)):
        scalar = int(answer["scalar"])
        if (answer["fixture_index"] != i
                or group.point_from_json(answer["point"]) != point
                or not group.on_curve(point)
                or group.multiply(point, group.R) is not None
                or group.multiply(group.G, scalar) != point):
            raise ValueError(f"fresh target {i}: independent scalar replay failed")
    return points, answers


def verify_ic(folder, point, scalar, selected):
    usage = one_json(folder / "usage.json")
    if usage["status"] != "PASS_PRODUCER":
        return {"status": usage["status"], "usage": usage}
    summary = one_json(folder / "stdout.jsonl")
    target = one_json(folder / "targets_result.jsonl")
    dump = one_json(folder / "base_dump.jsonl")
    k, variant = selected["k"], selected["variant"]
    layout = {"root11": "root_only83_11byte_states_reconstructed_second_root",
              "root21": "root_only83_21byte_states",
              "root22": "root_only83_22byte_states"}[variant]
    state_bytes = {"root11": 11, "root21": 21, "root22": 22}[variant]
    state_bound = 83 * k * (k - 1) // 2 + 42 * k
    if (summary["n"] != 83 or summary["a"] != 1
            or summary["targets"] != 1 or summary["targets_solved"] != 1
            or summary["targets_failed"] != 0 or summary["orbit_columns"] != k
            or summary["rank"] != k or summary["regular_states"] != state_bound
            or summary["root_table_layout"] != "indexed83_5byte_exact"
            or summary["root_table_slot_bytes"] != 5
            or summary["root_table_indexed5_enabled"] is not True
            or summary["root_state_layout"] != layout
            or summary["root_state_slot_bytes"] != state_bytes
            or summary["symmetric_s3_states_enabled"] is not True
            or summary["indexed5_load75_enabled"] is not (variant == "root22")
            or summary["indexed5_fixed30_enabled"] is not (variant != "root22")
            or summary["root_only11_enabled"] is not (variant == "root11")):
        raise ValueError("selected IC source layout, base or rank differs")
    if variant != "root22" and (summary["root_table_request_entries"] != 1 << 29
                                or summary["root_table_allocated_bytes"] != 5 * (1 << 30)):
        raise ValueError("selected IC fixed table allocation differs")
    base = group.replay_base(dump, summary)
    if len(base) != 166 * k:
        raise ValueError("selected IC factor base has wrong usable count")
    recovered = int(target["recovered_scalar"])
    if (target["fixture_index"] != 0
            or group.point_from_json(target["target"]) != point
            or recovered != scalar
            or target["group_verified"] is not True
            or group.multiply(group.G, recovered) != point):
        raise ValueError("selected IC target scalar replay failed")
    indices = target["point_indices"]
    if len(indices) != 4 or any(not isinstance(i, int) or i < 0 or i >= len(base)
                                for i in indices):
        raise ValueError("selected IC four-point witness indices invalid")
    relation = None
    for index in indices:
        relation = group.add(relation, base[index])
    if relation != point:
        raise ValueError("selected IC four-point relation failed")
    phase_ms = sum(float(target[key]) for key in
                   ("target_query_ms", "target_pdp_ms", "target_relation_check_ms",
                    "target_descent_ms", "target_recovery_check_ms"))
    online_ms = float(target["target_ms"])
    if abs(phase_ms - online_ms) > 1.0 or not math.isfinite(online_ms) or online_ms <= 0:
        raise ValueError("selected IC exclusive online interval invalid")
    return {"status": "PASS_INDEPENDENT_REPLAY", "online_ms": online_ms,
            "phase_sum_ms": phase_ms, "base_hash": summary["base_hash"],
            "root_table_entries": summary["root_table_entries"],
            "actual_factor_base_points": len(base), "setup_ms": summary["timing_ms"]["setup_total"],
            "peak_rss_bytes": usage["peak_rss_bytes"], "usage": usage}


def verify_rho(folder, point, scalar):
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
        raise ValueError("strong rho same-point scalar replay failed")
    return {"status": "PASS_INDEPENDENT_REPLAY", "online_ms": online_ms,
            "peak_rss_bytes": usage["peak_rss_bytes"], "usage": usage}


def pair_noise(pair, host, ic, rho):
    before = one_json(pair / "calibration_before.json")
    after = one_json(pair / "calibration_after.json")
    valid = (all(0 < float(x) < 10_000 for x in before["sample_ms"] + after["sample_ms"])
             and len(before["sample_ms"]) == len(after["sample_ms"]) == 3
             and tuple(x["sha256"] for x in before["checksums_and_cpu"]) ==
             CALIBRATION_SHA256
             and tuple(x["sha256"] for x in after["checksums_and_cpu"]) ==
             CALIBRATION_SHA256)
    checks = {}
    for arm, row in (("ic", ic), ("rho", rho)):
        if row["status"] != "PASS_INDEPENDENT_REPLAY":
            checks[arm] = {"status": "ARM_INCOMPLETE"}
            valid = False
            continue
        usage = row["usage"]
        cgroup_before = usage["cgroup_before"]
        cgroup_after = usage["cgroup_after"]
        throttle = counter_delta(cgroup_before["cpu.stat"],
                                 cgroup_after["cpu.stat"], "nr_throttled")
        oom = counter_delta(cgroup_before["memory.events"],
                            cgroup_after["memory.events"], "oom")
        oom_kill = counter_delta(cgroup_before["memory.events"],
                                 cgroup_after["memory.events"], "oom_kill")
        arm_ok = (throttle == 0 and oom == 0 and oom_kill == 0
                  and row["peak_rss_bytes"] <= runner.parent.RSS_CAP_BYTES)
        checks[arm] = {"status": "PASS" if arm_ok else "NOISE_OR_RESOURCE_FAILURE",
                       "cgroup_throttle_events": throttle, "oom_events": oom,
                       "oom_kills": oom_kill, "peak_rss_bytes": row["peak_rss_bytes"]}
        valid &= arm_ok
    system_before = one_json(pair / "system_before.json")
    system_after = one_json(pair / "system_after.json")
    cpus = host["execution_cpus"] or []
    steal_before = proc_steal_ticks(system_before["proc_stat"], cpus)
    steal_after = proc_steal_ticks(system_after["proc_stat"], cpus)
    steal_ticks = (steal_after - steal_before if steal_before is not None
                   and steal_after is not None and steal_after >= steal_before else None)
    # The whole pair includes IC setup; this is a conservative allocation-level guard.
    pair_seconds = sum(row.get("usage", {}).get("wall_seconds", 0.0)
                       for row in (ic, rho))
    ticks_per_second = host.get("clock_ticks_per_second")
    steal_seconds = (steal_ticks / ticks_per_second
                     if steal_ticks is not None and ticks_per_second else None)
    mean_steal_seconds = (steal_seconds / len(cpus)
                          if steal_seconds is not None and cpus else None)
    steal_ok = (mean_steal_seconds is not None and pair_seconds > 0
                and mean_steal_seconds <= 0.01 * pair_seconds)
    valid &= steal_ok
    return {"status": "PASS" if valid else "NOISE_GATE_FAILED",
            "calibration_before_ms": before["median_ms"],
            "calibration_after_ms": after["median_ms"],
            "steal_ticks_visible_cpus": steal_ticks,
            "steal_seconds_visible_cpus": steal_seconds,
            "mean_steal_seconds_per_execution_cpu": mean_steal_seconds,
            "steal_limit_fraction_of_pair_wall": 0.01,
            "arm_checks": checks}


def verify_pair(output, index, point, answer, host, selected):
    pair = output / "pairs" / f"p{index+1}"
    if not (pair / "order.json").exists():
        return {"pair": index + 1, "status": "NOT_RUN"}
    order = one_json(pair / "order.json")
    expected_order = ["ic", "rho"] if index % 2 == 0 else ["rho", "ic"]
    point_file = output / "inputs" / f"point_{index+1}.jsonl"
    if (order["target_index"] != index or order["order"] != expected_order
            or order["point_sha256"] != runner.parent.digest(point_file)
            or json_lines(point_file) != [json_lines(runner.TARGETS)[index]]):
        raise ValueError("paired point or alternating order changed")
    scalar = int(answer["scalar"])
    ic = verify_ic(pair / "ic", point, scalar, selected)
    rho = verify_rho(pair / "rho", point, scalar)
    noise = pair_noise(pair, host, ic, rho)
    status = ("PASS_TIER1_PAIR" if ic["status"] == rho["status"] ==
              "PASS_INDEPENDENT_REPLAY" and noise["status"] == "PASS"
              else "PAIR_INCOMPLETE_OR_NOISY")
    return {"pair": index + 1, "status": status, "order": expected_order,
            "ic": ic, "rho": rho, "noise": noise,
            "ratio_rho_over_ic": (rho["online_ms"] / ic["online_ms"]
                                  if status == "PASS_TIER1_PAIR" else None)}


def paired_bootstrap_interval(ratios):
    rng = random.Random(2026101405)
    resamples = []
    for _ in range(10_000):
        resamples.append(statistics.median(ratios[rng.randrange(len(ratios))]
                                           for _ in ratios))
    resamples.sort()
    return [resamples[249], resamples[9749]]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    output = args.output.resolve()
    points, answers = fresh_points()
    selected = runner.selected_config()
    host = one_json(output / "host.json")
    source = one_json(output / "source.json")
    if (source["selection_sha256"] != runner.parent.digest(runner.SELECTION)
            or source["selection"] != selected
            or source["hashes"][str(runner.TARGETS.relative_to(ROOT))] !=
            runner.TARGETS_SHA256):
        raise ValueError("executed source or selection receipt changed")
    rows = []
    for index, (point, answer) in enumerate(zip(points, answers)):
        try:
            rows.append(verify_pair(output, index, point, answer, host, selected))
        except Exception as error:
            rows.append({"pair": index + 1, "status": "VERIFICATION_FAILED",
                         "error": str(error)})
    valid = all(row["status"] == "PASS_TIER1_PAIR" for row in rows)
    calibrations = ([row["noise"][key] for row in rows for key in
                     ("calibration_before_ms", "calibration_after_ms")]
                    if valid else [])
    calibration_ratio = max(calibrations) / min(calibrations) if calibrations else None
    calibration_ok = calibration_ratio is not None and calibration_ratio <= 1.15
    ratios = [row["ratio_rho_over_ic"] for row in rows] if valid else []
    all_faster = valid and all(ratio > 1 for ratio in ratios)
    base_hashes = {row["ic"]["base_hash"] for row in rows} if valid else set()
    rank_base_consistent = len(base_hashes) == 1
    promote = valid and calibration_ok and all_faster and rank_base_consistent
    result = {"schema": "n83.tier1-final-analysis/v1",
              "hardware_tier": "tier1_cloud_instance",
              "status": "PASS_TIER1_PAIRED_SPEEDUP" if promote else "TIER1_GATE_INCOMPLETE",
              "target_count": len(points), "candidate_id": selected["candidate_id"],
              "pairs": rows, "all_five_correct_faster": all_faster,
              "base_hash_consistent": rank_base_consistent,
              "calibration_max_over_min": calibration_ratio,
              "calibration_bound": 1.15,
              "median_ratio_rho_over_ic": statistics.median(ratios) if promote else None,
              "median_ratio_target_block_bootstrap_95pct": (
                  paired_bootstrap_interval(ratios) if promote else None),
              "bootstrap_seed": 2026101405,
              "mean_ic_online_ms": (statistics.mean(row["ic"]["online_ms"] for row in rows)
                                    if promote else None),
              "mean_rho_online_ms": (statistics.mean(row["rho"]["online_ms"] for row in rows)
                                     if promote else None)}
    (output / "analysis.json").write_text(json.dumps(result, indent=2,
                                                    sort_keys=True) + "\n")
    print(json.dumps(result, sort_keys=True))
    return 0 if promote else 1


if __name__ == "__main__":
    raise SystemExit(main())
