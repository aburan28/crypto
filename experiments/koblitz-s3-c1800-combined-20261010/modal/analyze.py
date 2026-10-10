"""Verify the Modal receipt against local scalar fixtures and summarize pairs."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import re


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
PARENT = ROOT / "experiments/koblitz-s3-c1800-combined-20261010/CONFIRM_FREEZE.json"
PLAN = HERE / "PILOT_FREEZE.json"
BUNDLE = HERE / "public_bundle.json"
PHASES = ("target_query", "target_pdp_charged", "target_relation_check",
          "target_descent", "target_recovery_check")
QUALITY_CHECKS = ("verified_scalar_and_fixture", "exclusive_phase_sum",
                  "eight_physical_cores_selected", "whole_smt_siblings_selected",
                  "observed_thread_affinity", "observed_worker_threads",
                  "one_numa_node", "numa_pages_on_selected_node",
                  "unlimited_cpu_quota", "zero_cpu_throttling",
                  "zero_guest_steal", "zero_selected_cpu_interrupts",
                  "cpu_psi_at_most_1_percent",
                  "memory_io_psi_at_most_0_1_percent",
                  "zero_oom_events", "fixed_cpu_frequency",
                  "perf_counters_unmultiplexed")


def sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def cpu_list(raw: str | None) -> list[int]:
    result = set()
    for item in (raw or "").strip().split(","):
        if not item:
            continue
        pieces = item.split("-", 1)
        if all(piece.isdigit() for piece in pieces):
            result.update(range(int(pieces[0]), int(pieces[-1]) + 1))
    return sorted(result)


def named_ints(raw: str | None) -> dict:
    result = {}
    for line in (raw or "").splitlines():
        pieces = line.split()
        if len(pieces) >= 2 and pieces[1].isdigit():
            result[pieces[0]] = int(pieces[1])
    return result


def counter_delta(row: dict, name: str, field: str) -> int | None:
    before = named_ints(row["before"]["cgroup"].get(name))
    after = named_ints(row["after"]["cgroup"].get(name))
    if field not in before or field not in after:
        return None
    return after[field] - before[field]


def psi_fraction(row: dict, name: str) -> float | None:
    def total(raw: str | None) -> int | None:
        for line in (raw or "").splitlines():
            if line.startswith("some "):
                match = re.search(r"\btotal=(\d+)\b", line)
                if match:
                    return int(match.group(1))
        return None
    old = total(row["before"]["cgroup"].get(name))
    new = total(row["after"]["cgroup"].get(name))
    elapsed_us = row["outer_wall_ns"] / 1000
    return ((new - old) / elapsed_us if old is not None and new is not None
            and elapsed_us > 0 else None)


def difference(row: dict, category: str, keys: list[str], field: str | None = None) -> int | None:
    before, after = row["before"].get(category), row["after"].get(category)
    if not isinstance(before, dict) or not isinstance(after, dict):
        return None
    try:
        return sum((after[key][field] - before[key][field]) if field else
                   (after[key] - before[key]) for key in keys)
    except (KeyError, TypeError):
        return None


def guest_quality(row: dict, host: dict, verified: bool) -> dict:
    selected = host["selected_cpus"]
    entries = [host["topology"]["cpus"].get(str(cpu), {}) for cpu in selected]
    identities = [(entry.get("physical_package_id"), entry.get("core_id"))
                  for entry in entries]
    physical = (len(selected) == 16 and all(all(part is not None for part in pair)
                                           for pair in identities)
                and len(set(identities)) == 8)
    siblings = [cpu_list(entry.get("thread_siblings_list")) for entry in entries]
    sibling_complete = (all(siblings) and all(set(group).issubset(selected)
                                               for group in siblings)
                        if physical else None)
    nodes = [name for name, cpus in host["topology"]["numa_nodes"].items()
             if set(selected).issubset(cpu_list(cpus))]
    one_node = len(nodes) == 1 if host["topology"]["numa_nodes"] else None
    observed_nodes = set(row["sampled_numa_pages"])
    pages_one_node = (observed_nodes == {nodes[0].replace("node", "N")}
                      if one_node and observed_nodes else None)
    masks = row["thread_affinity_masks_seen"]
    affinity = (all(cpu_list(mask) == selected for mask in masks) if masks else None)
    quota = row["before"]["cgroup"].get("cpu.max")
    throttled = counter_delta(row, "cpu.stat", "nr_throttled")
    throttle_us = counter_delta(row, "cpu.stat", "throttled_usec")
    steal = difference(row, "cpu_ticks", [str(cpu) for cpu in selected], "steal")
    irq = difference(row, "interrupts", [str(cpu) for cpu in selected])
    cpu_psi = psi_fraction(row, "cpu.pressure")
    mem_psi = psi_fraction(row, "memory.pressure")
    io_psi = psi_fraction(row, "io.pressure")
    oom = counter_delta(row, "memory.events", "oom")
    oom_kill = counter_delta(row, "memory.events", "oom_kill")
    frequency = []
    for cpu in selected:
        governor = (host["cpu_governor"].get(str(cpu)) or "").strip()
        limits = host["cpu_frequency_limits_khz"].get(str(cpu), {})
        low = (limits.get("scaling_min_freq") or "").strip()
        high = (limits.get("scaling_max_freq") or "").strip()
        frequency.append((governor, low, high))
    fixed = (all(governor == "performance" and low and low == high
                 for governor, low, high in frequency)
             if all(governor and low and high for governor, low, high in frequency)
             else None)
    perf = row["perf"]
    perf_coverage = (None if perf["cycles_error"] or perf["instructions_error"] else
                     bool(perf["cycles"] and perf["instructions"] and
                          perf["both_events_at_least_99_percent_running"]))
    checks = {
        "verified_scalar_and_fixture": verified,
        "exclusive_phase_sum": row["phase_ok"],
        "eight_physical_cores_selected": physical,
        "whole_smt_siblings_selected": sibling_complete,
        "observed_thread_affinity": affinity,
        "observed_worker_threads": row["max_threads_seen"] >= 14,
        "one_numa_node": one_node,
        "numa_pages_on_selected_node": pages_one_node,
        "unlimited_cpu_quota": quota.split()[0] == "max" if quota else None,
        "zero_cpu_throttling": (throttled == 0 and throttle_us == 0
                                if throttled is not None and throttle_us is not None else None),
        "zero_guest_steal": steal == 0 if steal is not None else None,
        "zero_selected_cpu_interrupts": irq == 0 if irq is not None else None,
        "cpu_psi_at_most_1_percent": (0 <= cpu_psi <= .01 if cpu_psi is not None else None),
        "memory_io_psi_at_most_0_1_percent": (
            0 <= mem_psi <= .001 and 0 <= io_psi <= .001
            if mem_psi is not None and io_psi is not None else None),
        "zero_oom_events": (oom == 0 and oom_kill == 0
                            if oom is not None and oom_kill is not None else None),
        "fixed_cpu_frequency": fixed,
        "perf_counters_unmultiplexed": perf_coverage,
    }
    assert tuple(checks) == QUALITY_CHECKS
    return {
        "guest_score_0_to_100": round(100 * sum(value is True for value in checks.values()) /
                                      len(checks), 2),
        "available_check_count": sum(value is not None for value in checks.values()),
        "checks": checks,
        "failed_checks": [key for key, value in checks.items() if value is False],
        "unavailable_checks": [key for key, value in checks.items() if value is None],
        "counters": {"cpu_throttled_events": throttled,
                     "cpu_throttled_us": throttle_us, "selected_cpu_steal_ticks": steal,
                     "selected_cpu_interrupts": irq, "cgroup_cpu_psi_fraction": cpu_psi,
                     "cgroup_memory_psi_fraction": mem_psi,
                     "cgroup_io_psi_fraction": io_psi,
                     "oom_events": oom, "oom_kills": oom_kill,
                     "cycles_running_fraction": (perf["cycles_event"] or {}).get("running_fraction"),
                     "instructions_running_fraction":
                         (perf["instructions_event"] or {}).get("running_fraction")},
    }


def idle_noise(host: dict) -> dict:
    before, after = host["idle_before"], host["idle_after"]
    row = {"before": before, "after": after,
           "outer_wall_ns": after["unix_ns"] - before["unix_ns"]}
    selected = [str(cpu) for cpu in host["selected_cpus"]]
    return {"elapsed_ns": row["outer_wall_ns"],
            "cpu_psi_fraction": psi_fraction(row, "cpu.pressure"),
            "memory_psi_fraction": psi_fraction(row, "memory.pressure"),
            "io_psi_fraction": psi_fraction(row, "io.pressure"),
            "cpu_throttled_events": counter_delta(row, "cpu.stat", "nr_throttled"),
            "selected_cpu_steal_ticks": difference(row, "cpu_ticks", selected, "steal"),
            "selected_cpu_interrupts": difference(row, "interrupts", selected)}


def analyze(receipt_path: Path) -> dict:
    plan = json.loads(PLAN.read_bytes())
    assert plan["analyzer_sha256"] == sha(Path(__file__).read_bytes())
    bundle = json.loads(BUNDLE.read_bytes())
    parent = json.loads(PARENT.read_bytes())
    receipt = json.loads(receipt_path.read_bytes())
    assert receipt["kind"] == "s3_modal_c1800_launch_receipt_v1"
    assert receipt["freeze_sha256"] == sha(PLAN.read_bytes())
    inventory = receipt["transfer_inventory"]
    assert len(inventory) == plan["transferred_file_count"] == 3
    assert [row["remote"] for row in inventory] == [
        "/opt/s3/baseline", "/opt/s3/expanded", "/opt/s3/public_bundle.json"]
    assert [row["sha256"] for row in inventory] == [
        plan["binary_sha256"]["baseline"], plan["binary_sha256"]["expanded"],
        sha(BUNDLE.read_bytes())]
    assert sum(row["bytes"] for row in inventory) == plan["transferred_bytes"]
    assert inventory[-1]["bytes"] == BUNDLE.stat().st_size
    guest = receipt["guest"]
    if guest is None:
        return {"kind": "s3_modal_c1800_pilot_analysis_v1",
                "receipt_sha256": sha(receipt_path.read_bytes()),
                "status": "producer_failure", "error": receipt["error"],
                "rows": [], "controlled_speedup": None}
    assert guest["bundle_sha256"] == sha(BUNDLE.read_bytes())
    assert guest["binary_sha256"] == plan["binary_sha256"]
    assert guest["run_order"] == bundle["run_order"]
    cases = {case["n"]: case for case in bundle["curves"]}
    fixtures = {(item["n"], item["label"]): item for item in parent["targets"]
                if item["label"] == "T205" and item["n"] == 53}
    summaries = []
    for index, row in enumerate(guest["rows"]):
        n, arm = bundle["run_order"][index]
        case = cases[n]
        frozen = case["arms"][arm]
        target = fixtures[n, "T205"]
        report = row.get("report") or {}
        timing = report.get("timing_ms") or {}
        point = ROOT / target["point_file_relative"]
        point_ok = (sha(point.read_bytes()) == target["point_sha256"]
                    and json.loads(point.read_bytes()) == case["public_point"])
        phases_ok = (row["phase_ok"] and all(timing.get(k, -1) >= 0 for k in PHASES)
                     and abs(sum(timing[k] for k in PHASES) -
                             timing["target_online_after_reusable_setup"]) < .02)
        verified = bool(row["native_verified"] and point_ok and phases_ok
                        and report.get("recovered_scalar") == target["fixture_scalar"]
                        and report.get("group_verified") is True
                        and report.get("target") == case["public_point"]
                        and report.get("factor_base_digest") == frozen["factor_base_digest"]
                        and row["candidate_id"] == frozen["candidate_id"]
                        and row["workload_id"] == case["workload_id"])
        quality = guest_quality(row, guest["host"], verified)
        summaries.append({
            "n": n, "label": "T205", "arm": arm,
            "run_id": row["run_id"], "verified_native_and_fixture": verified,
            "online_ms": timing.get("target_online_after_reusable_setup"),
            "target_s3_calls": report.get("target_s3_calls"),
            "setup_ms": timing.get("reusable_setup_total"),
            "sampled_peak_rss_bytes": row["sampled_peak_rss_bytes"],
            "ipc": row["perf"]["ipc"],
            "ipc_scaled_estimate": row["perf"]["ipc_scaled_estimate"],
            "guest_quality": quality,
            "error": row["error"],
        })
    complete = len(summaries) == len(bundle["run_order"])
    verified = complete and all(row["verified_native_and_fixture"] for row in summaries)
    curves = {}
    if verified:
        for n in (53,):
            arms = {arm: [row["online_ms"] for row in summaries
                          if row["n"] == n and row["arm"] == arm]
                    for arm in ("baseline", "expanded")}
            assert all(len(times) == 2 and all(value > 0 for value in times)
                       for times in arms.values())
            aa = {arm: max(times) / min(times) for arm, times in arms.items()}
            curves[str(n)] = {
                "baseline_online_ms": arms["baseline"],
                "expanded_online_ms": arms["expanded"],
                "paired_exploratory_online_ratio": math.sqrt(math.prod(arms["baseline"])) /
                                                    math.sqrt(math.prod(arms["expanded"])),
                "aa_max_min_ratio": aa,
                "aa_within_10_percent": all(value <= 1.10 for value in aa.values()),
            }
    return {"kind": "s3_modal_c1800_pilot_analysis_v1",
            "receipt_sha256": sha(receipt_path.read_bytes()),
            "status": ("verified_complete" if verified else
                       "incomplete_or_invalid" if guest["rows"] else "no_rows"),
            "planned_rows": len(bundle["run_order"]),
            "completed_rows": len(summaries),
            "guest_stop_reason": guest["stop_reason"],
            "idle_noise": idle_noise(guest["host"]),
            "rows": summaries, "curves": curves,
            "guest_score_min": min((item["guest_quality"]["guest_score_0_to_100"]
                                    for item in summaries), default=None),
            "guest_score_mean": (sum(item["guest_quality"]["guest_score_0_to_100"]
                                     for item in summaries) / len(summaries)
                                 if summaries else None),
            "host_exclusive_partition_verified": False,
            "controlled_speedup": None}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("receipt", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = analyze(args.receipt)
    if args.output.exists():
        raise FileExistsError(args.output)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": result["status"], "completed_rows": result.get("completed_rows"),
                      "curves": result.get("curves"), "controlled_speedup": None},
                     sort_keys=True))
