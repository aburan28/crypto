"""Verify the Modal cross-curve pilot against local fixtures and score telemetry."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import re


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
PARENT = ROOT / "experiments/koblitz-s3-crosscurve-20261009/FREEZE.json"
PLAN = HERE / "PILOT_FREEZE.json"
BUNDLE = HERE / "public_bundle.json"
PHASES = ("target_query", "target_pdp_charged", "target_relation_check",
          "target_descent", "target_recovery_check")
QUALITY_CHECKS = ("verified_scalar_and_fixture", "exclusive_phase_sum",
                  "sixteen_physical_cores_selected", "smt_siblings_available",
                  "observed_thread_affinity", "observed_worker_threads",
                  "one_numa_node", "numa_pages_on_selected_node",
                  "unlimited_cpu_quota", "zero_cpu_throttling",
                  "zero_guest_steal", "zero_selected_cpu_interrupts",
                  "cpu_psi_at_most_1_percent",
                  "memory_io_psi_at_most_0_1_percent",
                  "zero_oom_events", "fixed_cpu_frequency",
                  "perf_cycles_and_instructions")


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
                and len(set(identities)) == 16)
    siblings = [cpu_list(entry.get("thread_siblings_list")) for entry in entries]
    sibling_available = (all(siblings) and all(set(group).issubset(host["allowed_cpus"])
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
    checks = {
        "verified_scalar_and_fixture": verified,
        "exclusive_phase_sum": row["phase_ok"],
        "sixteen_physical_cores_selected": physical,
        "smt_siblings_available": sibling_available,
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
        "perf_cycles_and_instructions": bool(row["perf"]["cycles"] and
                                             row["perf"]["instructions"]),
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
                     "oom_events": oom, "oom_kills": oom_kill},
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


def analyze(receipt_path: Path) -> tuple[dict, list[dict]]:
    plan = json.loads(PLAN.read_bytes())
    assert plan["analyzer_sha256"] == sha(Path(__file__).read_bytes())
    bundle = json.loads(BUNDLE.read_bytes())
    parent = json.loads(PARENT.read_bytes())
    receipt = json.loads(receipt_path.read_bytes())
    assert receipt["kind"] == "s3_crosscurve_modal_launch_receipt_v1"
    assert receipt["freeze_sha256"] == sha(PLAN.read_bytes())
    inventory = receipt["transfer_inventory"]
    assert len(inventory) == plan["transferred_file_count"] == 3
    assert [row["remote"] for row in inventory] == [
        "/opt/s3/reference", "/opt/s3/prepared", "/opt/s3/public_bundle.json"]
    assert [row["sha256"] for row in inventory] == [
        plan["binary_sha256"]["reference"], plan["binary_sha256"]["prepared"],
        sha(BUNDLE.read_bytes())]
    assert sum(row["bytes"] for row in inventory) == plan["transferred_bytes"]
    guest = receipt["guest"]
    if guest is None:
        return ({"kind": "s3_crosscurve_modal_pilot_analysis_v1",
                 "receipt_sha256": sha(receipt_path.read_bytes()),
                 "status": "producer_failure", "error": receipt["error"],
                 "rows": [], "controlled_speedup": None}, [])
    assert guest["kind"] == "s3_crosscurve_modal_pilot_guest_receipt_v1"
    assert guest["bundle_sha256"] == sha(BUNDLE.read_bytes())
    assert guest["binary_sha256"] == plan["binary_sha256"]
    assert guest["run_order"] == bundle["run_order"]
    cases = {(case["n"], case["label"]): case for case in bundle["cases"]}
    fixtures = {(item["n"], item["label"]): item for item in parent["targets"]}
    summaries, replay_rows = [], []
    for index, row in enumerate(guest["rows"]):
        frozen = bundle["run_order"][index]
        n, label, arm = frozen["n"], frozen["label"], frozen["arm"]
        case, target = cases[n, label], fixtures[n, label]
        report = row.get("report") or {}
        timing = report.get("timing_ms") or {}
        point = PARENT.parent / target["point_file_relative"]
        point_ok = (sha(point.read_bytes()) == target["point_sha256"]
                    and json.loads(point.read_bytes()) == case["public_point"])
        phases_ok = bool(row["phase_ok"] and all(timing.get(k, -1) >= 0 for k in PHASES)
                         and abs(sum(timing[k] for k in PHASES) -
                                 timing["target_online_after_reusable_setup"]) < .02)
        id_ok = all(row.get(key) == frozen[key] for key in
                    ("n", "label", "arm", "repetition", "candidate_id",
                     "workload_id", "run_id"))
        verified = bool(row["native_verified"] and point_ok and phases_ok and id_ok
                        and report.get("recovered_scalar") == target["fixture_scalar"]
                        and report.get("group_verified") is True
                        and report.get("target") == case["public_point"]
                        and report.get("factor_base_points") == case["factor_base_points"])
        quality = guest_quality(row, guest["host"], verified)
        summaries.append({
            "n": n, "label": label, "arm": arm, "run_id": row["run_id"],
            "verified_native_and_fixture": verified, "phase_ok": phases_ok,
            "online_ms": timing.get("target_online_after_reusable_setup"),
            "target_s3_calls": report.get("target_s3_calls"),
            "setup_ms": timing.get("reusable_setup_total"),
            "fallback_attempts": report.get("fallback_attempts", 0),
            "sampled_peak_rss_bytes": row["sampled_peak_rss_bytes"],
            "ipc": row["perf"]["ipc"], "guest_quality": quality,
            "error": row["error"],
        })
        replay_rows.append({
            "kind": "s3_crosscurve_single_target_run_v1",
            "n": n, "label": label, "arm": arm, "run_id": row["run_id"],
            "candidate_id": row["candidate_id"], "workload_id": row["workload_id"],
            "target_sha256": target["point_sha256"],
            "binary_sha256": plan["binary_sha256"][arm],
            "phase_ok": phases_ok, "verified_native_and_fixture": verified,
            "report": row.get("report"),
        })
    complete = len(summaries) == len(bundle["run_order"])
    verified = complete and all(row["verified_native_and_fixture"] for row in summaries)
    pairs = []
    if verified:
        for case in bundle["cases"]:
            selected = [row for row in summaries if
                        (row["n"], row["label"]) == (case["n"], case["label"])]
            arms = {arm: [row for row in selected if row["arm"] == arm]
                    for arm in ("reference", "prepared")}
            assert all(len(rows) == 2 for rows in arms.values())
            times = {arm: [row["online_ms"] for row in rows]
                     for arm, rows in arms.items()}
            calls = {arm: [row["target_s3_calls"] for row in rows]
                     for arm, rows in arms.items()}
            assert all(value > 0 for rows in (*times.values(), *calls.values())
                       for value in rows)
            aa = {arm: max(values) / min(values) for arm, values in times.items()}
            pairs.append({"n": case["n"], "label": case["label"],
                          "online_wall_ratio_exploratory":
                              math.sqrt(math.prod(times["reference"]) /
                                        math.prod(times["prepared"])),
                          "online_root_call_ratio":
                              math.sqrt(math.prod(calls["reference"]) /
                                        math.prod(calls["prepared"])),
                          "aa_max_min_ratio": aa,
                          "aa_within_10_percent": all(value <= 1.10 for value in aa.values())})
    return ({"kind": "s3_crosscurve_modal_pilot_analysis_v1",
             "receipt_sha256": sha(receipt_path.read_bytes()),
             "status": ("verified_complete" if verified else
                        "incomplete_or_invalid" if guest["rows"] else "no_rows"),
             "planned_rows": len(bundle["run_order"]),
             "completed_rows": len(summaries),
             "guest_stop_reason": guest["stop_reason"],
             "idle_noise": idle_noise(guest["host"]),
             "rows": summaries, "pairs": pairs,
             "guest_score_min": min((item["guest_quality"]["guest_score_0_to_100"]
                                     for item in summaries), default=None),
             "guest_score_mean": (sum(item["guest_quality"]["guest_score_0_to_100"]
                                      for item in summaries) / len(summaries)
                                  if summaries else None),
             "host_exclusive_partition_verified": False,
             "controlled_speedup": None}, replay_rows)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("receipt", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--replay-ledger", type=Path)
    args = parser.parse_args()
    result, replay_rows = analyze(args.receipt)
    if args.output.exists():
        raise FileExistsError(args.output)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    if args.replay_ledger:
        if args.replay_ledger.exists():
            raise FileExistsError(args.replay_ledger)
        args.replay_ledger.write_text("".join(
            json.dumps(row, sort_keys=True, separators=(",", ":")) + "\n"
            for row in replay_rows))
    print(json.dumps({"status": result["status"], "completed_rows": result.get("completed_rows"),
                      "pairs": result.get("pairs"), "controlled_speedup": None},
                     sort_keys=True))
