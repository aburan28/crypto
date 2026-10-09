"""Score a Modal rank-243 panel's guest noise and one-target pairing.

The score is a transparent fraction of passed guest checks, with unavailable
checks scored zero and listed separately. It never certifies host exclusivity.
"""

from __future__ import annotations

import argparse
from collections import Counter
import gzip
import hashlib
import json
import math
from pathlib import Path
import random
import re
import statistics


WORKERS = 14
CPU_PSI_FRACTION_CAP = 0.01
MEMORY_IO_PSI_FRACTION_CAP = 0.001
AA_ABS_LOG_CAP = math.log(1.10)
CHECKS = (
    "verified_scalar_and_fixture", "exclusive_phase_sum", "physical_core_selection",
    "observed_thread_affinity", "observed_worker_threads", "numa_memory_bind",
    "unlimited_cpu_quota", "zero_cpu_throttling", "zero_guest_steal",
    "zero_guest_interrupts", "cpu_psi_at_most_1_percent",
    "memory_io_psi_at_most_0_1_percent", "zero_oom_events", "fixed_cpu_frequency",
)


def parse_cpu_list(spec: str) -> list[int]:
    result = set()
    for token in spec.strip().split(","):
        if not token:
            continue
        bounds = token.split("-")
        result.update(range(int(bounds[0]), int(bounds[-1]) + 1))
    return sorted(result)


def parse_named_ints(value: str | None) -> dict[str, int] | None:
    if value is None:
        return None
    result = {}
    for line in value.splitlines():
        fields = line.split()
        if len(fields) >= 2 and fields[1].isdigit():
            result[fields[0]] = int(fields[1])
    return result


def counter_delta(before: dict, after: dict, category: str, filename: str,
                  field: str) -> int | None:
    old = parse_named_ints(before.get(category, {}).get(filename))
    new = parse_named_ints(after.get(category, {}).get(filename))
    if old is None or new is None or field not in old or field not in new:
        return None
    return new[field] - old[field]


def psi_total(value: str | None, kind: str = "some") -> int | None:
    if value is None:
        return None
    for line in value.splitlines():
        if line.startswith(kind + " "):
            match = re.search(r"\btotal=(\d+)\b", line)
            if match:
                return int(match.group(1))
    return None


def psi_fraction(before: dict, after: dict, category: str, name: str,
                 elapsed_us: float) -> float | None:
    old = psi_total(before.get(category, {}).get(name))
    new = psi_total(after.get(category, {}).get(name))
    if old is None or new is None or elapsed_us <= 0:
        return None
    return (new - old) / elapsed_us


def dict_counter_delta(before: dict, after: dict, category: str,
                       keys: list[str]) -> int | None:
    old = before.get(category, {})
    new = after.get(category, {})
    if not isinstance(old, dict) or not isinstance(new, dict):
        return None
    if not all(key in old and key in new for key in keys):
        return None
    return sum(new[key] - old[key] for key in keys)


def fixed_frequency(host: dict, selected: list[int]) -> bool | None:
    by_cpu = host.get("frequency_by_cpu")
    if not isinstance(by_cpu, dict):
        return None
    for cpu in selected:
        item = by_cpu.get(str(cpu), {})
        governor = (item.get("scaling_governor") or "").strip()
        minimum = (item.get("scaling_min_freq") or "").strip()
        maximum = (item.get("scaling_max_freq") or "").strip()
        if not governor or not minimum or not maximum:
            return None
        if governor != "performance" or minimum != maximum:
            return False
    return True


def physical_selection(selection: dict) -> bool:
    selected = selection.get("selected", [])
    siblings = [tuple(group) for group in selection.get("selected_siblings", [])]
    sampler = tuple(selection.get("sampler_siblings", []))
    return (len(selected) == WORKERS + 1 and len(siblings) == len(selected)
            and len(set(siblings)) == len(siblings) and len(sampler) > 0
            and sampler not in siblings and all(0 not in group for group in siblings)
            and all(cpu in group for cpu, group in zip(selected, siblings)))


def run_quality(row: dict, document: dict) -> dict:
    host = document["host"]
    selection = host["selection"]
    selected = sorted(selection["selected"])
    before, after = row["before"], row["after"]
    elapsed_us = row["outer_wall_ns"] / 1000
    cgroup_cpu_psi = psi_fraction(before, after, "cgroup", "cpu.pressure", elapsed_us)
    cgroup_mem_psi = psi_fraction(before, after, "cgroup", "memory.pressure", elapsed_us)
    cgroup_io_psi = psi_fraction(before, after, "cgroup", "io.pressure", elapsed_us)
    throttled = counter_delta(before, after, "cgroup", "cpu.stat", "nr_throttled")
    throttle_us = counter_delta(before, after, "cgroup", "cpu.stat", "throttled_usec")
    oom = counter_delta(before, after, "cgroup", "memory.events", "oom")
    oom_kill = counter_delta(before, after, "cgroup", "memory.events", "oom_kill")
    steal = dict_counter_delta(before, after, "steal_ticks", [f"cpu{cpu}" for cpu in selected])
    old_irq = before.get("interrupts", {})
    new_irq = after.get("interrupts", {})
    irq = (dict_counter_delta(old_irq, new_irq, "totals", [str(cpu) for cpu in selected])
           if old_irq.get("available") and new_irq.get("available") else None)
    masks = row.get("thread_masks_seen", [])
    affinity = bool(masks) and all(parse_cpu_list(mask) == selected for mask in masks)
    quota = host.get("cgroup", {}).get("cpu.max")
    checks = {
        "verified_scalar_and_fixture": row.get("verified_native_and_fixture") is True,
        "exclusive_phase_sum": row.get("phase_ok") is True,
        "physical_core_selection": physical_selection(selection),
        "observed_thread_affinity": affinity if masks else None,
        "observed_worker_threads": row.get("max_threads_seen", 0) >= WORKERS,
        "numa_memory_bind": (document.get("numa_bind_used") is True
                             and selection.get("node") is not None),
        "unlimited_cpu_quota": quota.split()[0] == "max" if quota else None,
        "zero_cpu_throttling": (throttled == 0 and throttle_us == 0
                                if throttled is not None and throttle_us is not None else None),
        "zero_guest_steal": steal == 0 if steal is not None else None,
        "zero_guest_interrupts": irq == 0 if irq is not None else None,
        "cpu_psi_at_most_1_percent": (0 <= cgroup_cpu_psi <= CPU_PSI_FRACTION_CAP
                                      if cgroup_cpu_psi is not None else None),
        "memory_io_psi_at_most_0_1_percent": (
            0 <= cgroup_mem_psi <= MEMORY_IO_PSI_FRACTION_CAP
            and 0 <= cgroup_io_psi <= MEMORY_IO_PSI_FRACTION_CAP
            if cgroup_mem_psi is not None and cgroup_io_psi is not None else None),
        "zero_oom_events": (oom == 0 and oom_kill == 0
                            if oom is not None and oom_kill is not None else None),
        "fixed_cpu_frequency": fixed_frequency(host, selected),
    }
    assert tuple(checks) == CHECKS
    return {
        "tag": row["tag"], "family": row["family"], "label": row["label"],
        "arm": row["arm"], "checks": checks,
        "guest_score_0_to_100": round(100 * sum(value is True for value in checks.values()) / len(CHECKS), 2),
        "available_check_count": sum(value is not None for value in checks.values()),
        "failed_checks": [name for name, value in checks.items() if value is False],
        "unavailable_checks": [name for name, value in checks.items() if value is None],
        "counters": {"cpu_throttled_events": throttled, "cpu_throttled_us": throttle_us,
                     "selected_cpu_steal_ticks": steal, "selected_cpu_interrupts": irq,
                     "cgroup_cpu_psi_fraction": cgroup_cpu_psi,
                     "cgroup_memory_psi_fraction": cgroup_mem_psi,
                     "cgroup_io_psi_fraction": cgroup_io_psi,
                     "oom_events": oom, "oom_kills": oom_kill,
                     "memory_peak_after": after.get("cgroup", {}).get("memory.peak")},
    }


def percentile(sorted_values: list[float], p: float) -> float:
    return sorted_values[max(0, math.ceil(p * len(sorted_values)) - 1)]


def wilson(successes: int, total: int) -> list[float] | None:
    if total == 0:
        return None
    z = 1.959963984540054
    p = successes / total
    denominator = 1 + z * z / total
    midpoint = (p + z * z / (2 * total)) / denominator
    radius = z * math.sqrt(p * (1 - p) / total + z * z / (4 * total * total)) / denominator
    return [max(0, midpoint - radius), min(1, midpoint + radius)]


def summarize_blocks(document: dict) -> tuple[list[dict], dict]:
    by_tag = {row["tag"]: row for row in document["runs"]}
    assert len(by_tag) == len(document["runs"])
    blocks = []
    for block in document["blocks"]:
        rows = [by_tag[tag] for tag in block["run_tags"]]
        assert len(rows) == 5 and Counter(row["arm"] for row in rows) == {
            "reference": 2, "deep": 2, "shallow": 1}
        assert len({tuple(row["point"]) for row in rows}) == 1
        valid = all(row["verified_native_and_fixture"] and row["phase_ok"] for row in rows)
        ratio = aa = None
        if valid:
            times = [(row["arm"], json.loads(row["raw_text"])["timing_ms"][
                "target_online_after_reusable_setup"]) for row in rows]
            by_arm = {arm: [value for name, value in times if name == arm]
                      for arm in ("reference", "deep")}
            assert all(value > 0 for values in by_arm.values() for value in values)
            ratio = math.sqrt(math.prod(by_arm["reference"]) / math.prod(by_arm["deep"]))
            aa = abs(math.log(by_arm["deep"][0] / by_arm["deep"][1]))
            assert math.isclose(ratio, block["exploratory_online_ratio"], rel_tol=1e-10)
            assert math.isclose(aa, block["deep_aa_abs_log_ratio"], rel_tol=1e-10)
        blocks.append({"label": block["label"], "family": block["family"],
                       "all_verified": valid, "online_ratio": ratio,
                       "deep_aa_abs_log_ratio": aa})
    summaries = {}
    for family in ("fresh", "full_rank_stress_control"):
        selected = [item for item in blocks if item["family"] == family]
        ratios = [item["online_ratio"] for item in selected if item["online_ratio"] is not None]
        aas = sorted(item["deep_aa_abs_log_ratio"] for item in selected
                     if item["deep_aa_abs_log_ratio"] is not None)
        interval = None
        if len(ratios) >= 4:
            rng = random.Random(20261009)
            logs = list(map(math.log, ratios))
            samples = sorted(math.exp(statistics.mean(rng.choices(logs, k=len(logs))))
                             for _ in range(20_000))
            interval = [samples[499], samples[19_499]]
        required = None
        if len(ratios) >= 4 and statistics.mean(map(math.log, ratios)) > math.log(2):
            margin = statistics.mean(map(math.log, ratios)) - math.log(2)
            required = max(24, math.ceil(((1.6448536269514722 + 0.8416212335729143)
                                           * statistics.stdev(map(math.log, ratios)) / margin) ** 2))
        summaries[family] = {
            "target_count": len(selected), "verified_count": len(ratios),
            "minimum_online_ratio": min(ratios) if ratios else None,
            "at_least_2_count": sum(value >= 2 for value in ratios),
            "at_least_2_wilson_95": wilson(sum(value >= 2 for value in ratios), len(selected)),
            "geometric_mean_online_ratio": math.exp(statistics.mean(map(math.log, ratios))) if ratios else None,
            "target_bootstrap_95_descriptive": interval,
            "deep_aa_p95_abs_log_ratio": percentile(aas, .95) if aas else None,
            "aa_10_percent_gate": percentile(aas, .95) <= AA_ABS_LOG_CAP if aas else None,
            "pilot_pairs_per_curve_for_geomean_above_2_at_one_sided_95_and_80_percent_power": required,
        }
    return blocks, summaries


def perf_ipc(document: dict) -> list[dict]:
    result = []
    for item in document.get("perf_whole_process", []):
        counters = {}
        run = item.get("result", {})
        for line in (run.get("stderr") or "").splitlines():
            parts = [part.strip() for part in line.split(",")]
            for event in ("cycles", "instructions"):
                if event in parts:
                    try:
                        counters[event] = float(parts[0])
                    except ValueError:
                        pass
        cycles = counters.get("cycles")
        instructions = counters.get("instructions")
        result.append({"arm": item["arm"], "scope": run.get("scope"),
                       "perf_exit_code": run.get("exit_code"),
                       "cycles": cycles, "instructions": instructions,
                       "whole_process_ipc": instructions / cycles if cycles and instructions else None})
    return result


def main(input_path: Path, output_path: Path) -> None:
    raw = input_path.read_bytes()
    document = json.loads(gzip.decompress(raw) if raw[:2] == b"\x1f\x8b" else raw)
    assert document["kind"] == "modal_vm_s3_rank243_single_target_panel_v1"
    assert document["strict_host_isolation_receipt"] is None
    assert document["controlled_speedup"] is None
    runs = [run_quality(row, document) for row in document["runs"]]
    blocks, summaries = summarize_blocks(document)
    aas = [value["aa_10_percent_gate"] for value in summaries.values()
           if value["aa_10_percent_gate"] is not None]
    score = min(row["guest_score_0_to_100"] for row in runs) if runs else 0
    result = {
        "kind": "modal_s3_rank243_guest_fidelity_analysis_v1",
        "source_panel": str(input_path),
        "source_panel_sha256": hashlib.sha256(raw).hexdigest(),
        "score_definition": "100 times passed checks divided by 14; missing checks score zero; panel score is the minimum run score",
        "thresholds": {"cpu_psi_fraction": CPU_PSI_FRACTION_CAP,
                       "memory_io_psi_fraction": MEMORY_IO_PSI_FRACTION_CAP,
                       "aa_p95_abs_log_ratio": AA_ABS_LOG_CAP},
        "guest_fidelity_score_0_to_100": score,
        "all_guest_checks_pass": bool(runs) and bool(aas) and all(
            all(value is True for value in row["checks"].values()) for row in runs) and all(aas),
        "strict_host_isolation_receipt": None,
        "controlled_speedup": None,
        "run_count": len(runs), "block_count": len(blocks),
        "run_quality": runs, "blocks": blocks, "summary": summaries,
        "whole_process_perf": perf_ipc(document),
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"output": str(output_path), "guest_fidelity_score": score,
                      "run_count": len(runs), "block_count": len(blocks),
                      "fresh": summaries["fresh"], "stress": summaries["full_rank_stress_control"]},
                     sort_keys=True))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    main(args.input, args.output)
