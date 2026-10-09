"""Validate a Modal C1200 holdout ledger and score its guest timing fidelity."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
from pathlib import Path
import random
import re
import statistics


HERE = Path(__file__).resolve().parent
PHASES = ("target_query", "target_pdp_charged", "target_relation_check",
          "target_descent", "target_recovery_check")
WEIGHTS = {"topology": 10, "affinity": 15, "numa_binding": 10,
           "cpu_quota": 10, "no_throttling": 15, "no_steal": 10,
           "cpu_psi": 10, "memory_psi": 5, "interrupts": 5,
           "perf_ipc_available": 5, "frequency_stability": 5}
assert sum(WEIGHTS.values()) == 100


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def decide(values: list[bool | None]) -> str:
    return "fail" if False in values else "unknown" if None in values else "pass"


def cpu_set(spec: str) -> set[int]:
    result: set[int] = set()
    for token in spec.split(","):
        if token:
            ends = token.split("-")
            result.update(range(int(ends[0]), int(ends[-1]) + 1))
    return result


def counter(raw: str | None, name: str) -> int | None:
    if not raw:
        return None
    match = re.search(rf"^{re.escape(name)}\s+(\d+)$", raw, re.M)
    return int(match.group(1)) if match else None


def psi_total(raw: str | None, line: str = "some") -> int | None:
    if not raw:
        return None
    match = re.search(rf"^{line}\s+[^\n]*\btotal=(\d+)", raw, re.M)
    return int(match.group(1)) if match else None


def delta(before: int | None, after: int | None) -> int | None:
    return after - before if before is not None and after is not None and after >= before else None


def whole_process_perf(doc: dict) -> dict:
    result = {}
    for entry in doc.get("perf_whole_process", []):
        probe = entry["result"]
        counts = {}
        for line in probe.get("stderr", "").splitlines():
            parts = line.split(",")
            if len(parts) >= 3 and parts[2].strip() in ("cycles", "instructions"):
                value = parts[0].strip()
                if re.fullmatch(r"\d+(?:\.\d+)?", value):
                    counts[parts[2].strip()] = float(value)
        cycles = counts.get("cycles")
        instructions = counts.get("instructions")
        result[entry["arm"]] = {
            "scope": probe.get("scope"), "exit_code": probe.get("exit_code"),
            "cycles": cycles, "instructions": instructions,
            "instructions_per_cycle": (instructions / cycles if cycles and instructions else None),
            "error": probe.get("error") or probe.get("stderr", "")[-2000:] if
                     probe.get("exit_code") != 0 else None}
    return result


def fidelity(doc: dict, rows: list[dict]) -> dict:
    selection = doc["host"]["selection"]
    selected = set(selection["selected"])
    siblings = [set(group) for group in selection["selected_siblings"]]
    checks: dict[str, str] = {}
    checks["topology"] = decide([
        len(selected) == 15 and selection["distinct_physical_cpus"] == 15,
        selection["sampler_cpu"] not in selected and selection["node"] is not None,
        all(not left.intersection(right) for i, left in enumerate(siblings)
            for right in siblings[i + 1:])])
    checks["affinity"] = decide([
        bool(row["thread_masks_seen"]) and all(
            cpu_set(mask) and cpu_set(mask) <= selected
            for mask in row["thread_masks_seen"])
        for row in rows])
    checks["numa_binding"] = decide([doc.get("numa_bind_used") is True,
                                      doc.get("numa_bind_probe", {}).get("exit_code") == 0])
    quota = doc["host"]["cgroup"].get("cpu.max")
    if quota:
        fields = quota.split()
        permitted = (True if fields[0] == "max" else
                     int(fields[0]) / int(fields[1]) >= 16 if len(fields) == 2 else False)
    else:
        permitted = None
    checks["cpu_quota"] = decide([permitted])
    throttling, stealing, cpu_psi, memory_psi, interrupts, frequency = [], [], [], [], [], []
    for row in rows:
        before, after = row["before"], row["after"]
        elapsed_us = row["outer_wall_ns"] / 1000
        stat_before, stat_after = before["cgroup"].get("cpu.stat"), after["cgroup"].get("cpu.stat")
        throttled = delta(counter(stat_before, "nr_throttled"),
                          counter(stat_after, "nr_throttled"))
        throttling.append(throttled == 0 if throttled is not None else None)
        steal_before, steal_after = before.get("steal_ticks", {}), after.get("steal_ticks", {})
        stealing.append(all(steal_after.get(f"cpu{cpu}") == steal_before.get(f"cpu{cpu}")
                            for cpu in selected) if all(f"cpu{cpu}" in steal_before and
                            f"cpu{cpu}" in steal_after for cpu in selected) else None)
        for name, checks_list, limit in (("cpu", cpu_psi, .005), ("memory", memory_psi, 0.0)):
            pressure = delta(psi_total(before["cgroup"].get(f"{name}.pressure")),
                             psi_total(after["cgroup"].get(f"{name}.pressure")))
            checks_list.append(pressure / elapsed_us <= limit if pressure is not None and
                               elapsed_us > 0 else None)
        before_irq, after_irq = before["interrupts"], after["interrupts"]
        if before_irq.get("available") and after_irq.get("available"):
            changes = [delta(before_irq["totals"].get(str(cpu)),
                             after_irq["totals"].get(str(cpu))) for cpu in selected]
            interrupts.append(sum(changes) / (elapsed_us / 1e6) <= 50 * len(selected)
                              if all(value is not None for value in changes) and elapsed_us > 0
                              else None)
        else:
            interrupts.append(None)
        changes = []
        for cpu in selected:
            start = before["cpufreq_cur_khz"].get(str(cpu), {}).get("scaling_cur_freq")
            end = after["cpufreq_cur_khz"].get(str(cpu), {}).get("scaling_cur_freq")
            if not start or not end or not start.strip().isdigit() or not end.strip().isdigit():
                changes.append(None)
            else:
                changes.append(abs(int(end) - int(start)) / int(start) <= .05)
        frequency.append({"pass": True, "fail": False, "unknown": None}[decide(changes)])
    checks["no_throttling"] = decide(throttling)
    checks["no_steal"] = decide(stealing)
    checks["cpu_psi"] = decide(cpu_psi)
    checks["memory_psi"] = decide(memory_psi)
    checks["interrupts"] = decide(interrupts)
    checks["frequency_stability"] = decide(frequency)
    perf = whole_process_perf(doc)
    checks["perf_ipc_available"] = decide([
        True if perf.get(arm, {}).get("instructions_per_cycle") is not None else None
        for arm in ("baseline", "expanded")])
    lower = sum(WEIGHTS[key] for key, status in checks.items() if status == "pass")
    upper = lower + sum(WEIGHTS[key] for key, status in checks.items() if status == "unknown")
    return {"guest_score_lower": lower, "guest_score_upper": upper,
            "component_weights": WEIGHTS, "checks": checks,
            "perf_whole_process": perf,
            "scope": "guest observations around whole solver processes; online work is timed by the solver",
            "host_exclusive_partition_verified": False}


def bootstrap(blocks: list[dict]) -> list[float]:
    rng = random.Random(20261009)
    values = []
    for _ in range(20_000):
        sampled = []
        for block in rng.choices(blocks, k=len(blocks)):
            logs = {arm: list(map(math.log, block[arm])) for arm in ("baseline", "expanded")}
            sampled.append(statistics.mean(rng.choices(logs["baseline"], k=2)) -
                           statistics.mean(rng.choices(logs["expanded"], k=2)))
        values.append(math.exp(statistics.mean(sampled)))
    values.sort()
    return [values[499], values[19_499]]


def repeat_plan(blocks: list[dict]) -> dict:
    effects, within = [], {arm: [] for arm in ("baseline", "expanded")}
    for block in blocks:
        logs = {arm: list(map(math.log, block[arm])) for arm in within}
        effects.append(statistics.mean(logs["baseline"]) - statistics.mean(logs["expanded"]))
        for arm in within:
            within[arm].append(statistics.variance(logs[arm]))
    observed = statistics.variance(effects)
    within_by_arm = {arm: statistics.mean(values) for arm, values in within.items()}
    within_sum = sum(within_by_arm.values())
    between = max(0.0, observed - within_sum / 2)
    margin = statistics.mean(effects) - math.log(2)
    choices = {}
    for repeats in (2, 4, 8, 16):
        required = max(24, math.ceil((1.96 + .84) ** 2 *
                                     (between + within_sum / repeats) / margin ** 2)) if margin > 0 else None
        choices[str(repeats)] = {"independent_target_pairs_required": required,
                                 "single_target_solves_per_curve": (2 * repeats * required
                                                                    if required is not None else None)}
    return {"margin_above_log_2": margin,
            "assumption": "independent targets and repetitions under a stable host and target law",
            "two_sided_95_percent_interval_80_percent_power": choices}


def analyze(raw_path: Path, output: Path) -> dict:
    if output.exists():
        raise FileExistsError(output)
    payload = raw_path.read_bytes()
    doc = json.loads(gzip.decompress(payload) if payload[:2] == b"\x1f\x8b" else payload)
    plan_bytes = (HERE / "MODAL_HOLDOUT_FREEZE.json").read_bytes()
    plan = json.loads(plan_bytes)
    manifest_bytes = (HERE / "MODAL_UPLOAD_MANIFEST.json").read_bytes()
    assert doc["kind"] == "modal_s3_fullrank_panel_v1"
    assert doc["holdout_freeze_sha256"] == sha(plan_bytes)
    assert doc["upload_manifest_sha256"] == sha(manifest_bytes)
    assert doc["strict_host_isolation_receipt"] is None and doc["controlled_speedup"] is None
    assert len(doc["runs"]) == 98 and len(doc["blocks"]) == 24
    failures, timing_blocks = [], []
    warmup_failures = [{"arm": row["arm"], "error": row.get("error")}
                       for row in doc["runs"][:2] if not row["verified_native_and_fixture"]]
    for index, target in enumerate(plan["targets"]):
        rows = doc["runs"][2 + 4 * index:2 + 4 * (index + 1)]
        expected_order = (["baseline", "expanded", "expanded", "baseline"] if index % 2
                          else ["expanded", "baseline", "baseline", "expanded"])
        assert [row["arm"] for row in rows] == expected_order
        assert all(row["label"] == target["label"] for row in rows)
        by_arm = {arm: [] for arm in ("baseline", "expanded")}
        for row in rows:
            arm = row["arm"]
            frozen = plan["arms"][arm]
            assert row["candidate_id"] == frozen["candidate_id"]
            assert row["workload_id"] == target["workload_id"]
            assert row["point_sha256"] == target["point_sha256"]
            assert row["source_sha256"] == frozen["source_sha256"]
            assert row["run_id"] == f"{row['candidate_id']}W{row['workload_id']}R{row['repetition']}"
            assert row["repetition"] == 1 + len(by_arm[arm])
            report = row.get("report")
            if not row["verified_native_and_fixture"] or report is None or row.get("error"):
                failures.append({"label": target["label"], "arm": arm,
                                 "repetition": row["repetition"], "error": row.get("error")})
                by_arm[arm].append(None)
                continue
            timing = report["timing_ms"]
            assert row["phase_ok"] and all(timing[key] >= 0 for key in PHASES)
            assert abs(sum(timing[key] for key in PHASES) -
                       timing["target_online_after_reusable_setup"]) < .02
            assert report["target"] == target["public_point"]
            assert report["recovered_scalar"] == target["fixture_scalar"]
            assert report["factor_base_digest"] == frozen["factor_base_digest"]
            if arm == "expanded":
                assert report["prepared_rank"] == 1200
                assert report["target_full_rank_direct_recovery"] is True
                assert report["target_span_checks"] == 0
            by_arm[arm].append(timing["target_online_after_reusable_setup"])
        timing_blocks.append(by_arm)
    assert all(len(block[arm]) == 2 for block in timing_blocks for arm in block)
    observed_rows = doc["runs"][2:]
    result = {"kind": "modal_s3_fullrank_analysis_v1", "raw_sha256": sha(payload),
              "holdout_freeze_sha256": sha(plan_bytes), "run_count": len(observed_rows),
              "failure_rows": failures, "warmup_failure_rows": warmup_failures,
              "all_native_verified": not failures and not warmup_failures,
              "controlled_speedup": None, "fidelity": fidelity(doc, observed_rows)}
    if not failures:
        ratios = [math.exp(statistics.mean(map(math.log, block["baseline"])) -
                           statistics.mean(map(math.log, block["expanded"])))
                  for block in timing_blocks]
        result["exploratory_online_ratio"] = math.exp(statistics.mean(map(math.log, ratios)))
        result["hierarchical_bootstrap_95"] = bootstrap(timing_blocks)
        result["points_at_least_2x"] = sum(value >= 2 for value in ratios)
        result["target_ratios"] = dict(zip((target["label"] for target in plan["targets"]), ratios))
        result["repeat_count_plan"] = repeat_plan(timing_blocks)
        result["aa_p95_abs_log_ratio_by_arm"] = {
            arm: sorted(abs(math.log(block[arm][0] / block[arm][1]))
                        for block in timing_blocks)[math.ceil(.95 * len(timing_blocks)) - 1]
            for arm in ("baseline", "expanded")}
        result["guest_noise_gate_pass"] = all(value <= .10 for value in
                                               result["aa_p95_abs_log_ratio_by_arm"].values())
        result["guest_high_fidelity"] = (not warmup_failures and result["guest_noise_gate_pass"] and
                                         result["fidelity"]["guest_score_lower"] >= 90 and
                                         all(value != "fail" for value in
                                             result["fidelity"]["checks"].values()))
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    summary = analyze(args.raw, args.output)
    print(json.dumps({key: summary.get(key) for key in
                      ("run_count", "failure_rows", "exploratory_online_ratio",
                       "hierarchical_bootstrap_95", "guest_high_fidelity",
                       "controlled_speedup")}, sort_keys=True))
