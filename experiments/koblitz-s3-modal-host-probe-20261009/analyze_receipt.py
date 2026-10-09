"""Score guest telemetry availability without promoting a CPU speedup claim."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def measure(receipt: dict) -> dict:
    remote = receipt.get("remote") or {}
    selected = str(remote.get("selected_cpu"))
    before = remote.get("before") or {}
    idle = remote.get("after_idle") or {}
    pin = remote.get("after_pin") or {}
    busy = remote.get("after_busy") or {}

    def psi_total(snapshot: dict, path: str) -> int | None:
        raw = snapshot.get(path)
        if not raw:
            return None
        for line in raw.splitlines():
            if line.startswith("some "):
                for field in line.split():
                    if field.startswith("total="):
                        return int(field.split("=", 1)[1])
        return None

    def cpu_stat(snapshot: dict, key: str) -> int | None:
        raw = (snapshot.get("cgroup") or {}).get("cpu.stat") or ""
        for line in raw.splitlines():
            parts = line.split()
            if len(parts) == 2 and parts[0] == key:
                return int(parts[1])
        return None

    def delta(a: int | None, b: int | None) -> int | None:
        return b - a if a is not None and b is not None else None

    def per_cpu_delta(field: str, a: dict, b: dict) -> int | None:
        return delta((a.get(field) or {}).get(selected),
                     (b.get(field) or {}).get(selected))

    perf = remote.get("perf") or {}
    cg = busy.get("cgroup") or {}
    checks = {
        "guest_pin_held": remote.get("pin_error") is None
                          and pin.get("affinity") == [remote.get("selected_cpu")]
                          and busy.get("affinity") == [remote.get("selected_cpu")],
        "guest_cpu_topology_visible": bool(remote.get("cpu_model"))
                                      and bool(remote.get("topology")),
        "guest_numa_policy_visible": bool((busy.get("process_numa_maps") or {})
                                          .get("available"))
                                     and bool(remote.get("numa_nodes")),
        "guest_cgroup_cpu_visible": cg.get("cpu.max") is not None
                                    and cg.get("cpu.stat") is not None,
        "guest_psi_visible": all(
            psi_total(s, "system_cpu_psi") is not None
            and psi_total(s, "system_memory_psi") is not None
            and psi_total({"system_cpu_psi": (s.get("cgroup") or {}).get("cpu.pressure")},
                          "system_cpu_psi") is not None
            and psi_total({"system_cpu_psi": (s.get("cgroup") or {}).get("memory.pressure")},
                          "system_cpu_psi") is not None
            for s in (before, idle, pin, busy)),
        "guest_steal_irq_visible": all(
            (s.get("steal_ticks") or {}).get(selected) is not None
            and (s.get("interrupt_totals") or {}).get(selected) is not None
            for s in (before, idle, pin, busy)),
        "process_perf_ipc_visible": perf.get("cycles") is not None
                                    and perf.get("instructions") is not None
                                    and perf.get("ipc") is not None,
    }
    result = {
        "kind": "source_free_modal_vm_cpu_probe_analysis_v1",
        "sandbox_id": receipt.get("sandbox_id"),
        "remote_exit_code": receipt.get("remote_exit_code"),
        "guest_observability_checks": checks,
        "guest_observability_score": sum(checks.values()),
        "guest_observability_max": len(checks),
        "guest_observability_fraction": sum(checks.values()) / len(checks),
        "idle_system_cpu_psi_some_us": delta(
            psi_total(before, "system_cpu_psi"),
            psi_total(idle, "system_cpu_psi")),
        "busy_system_cpu_psi_some_us": delta(
            psi_total(pin, "system_cpu_psi"),
            psi_total(busy, "system_cpu_psi")),
        "idle_cgroup_throttled_us": delta(
            cpu_stat(before, "throttled_usec"), cpu_stat(idle, "throttled_usec")),
        "busy_cgroup_throttled_us": delta(
            cpu_stat(pin, "throttled_usec"), cpu_stat(busy, "throttled_usec")),
        "idle_selected_cpu_steal_ticks": per_cpu_delta("steal_ticks", before, idle),
        "busy_selected_cpu_steal_ticks": per_cpu_delta("steal_ticks", pin, busy),
        "idle_selected_cpu_interrupts": per_cpu_delta("interrupt_totals", before, idle),
        "busy_selected_cpu_interrupts": per_cpu_delta("interrupt_totals", pin, busy),
        "busy_involuntary_context_switches": delta(
            pin.get("involuntary_context_switches"),
            busy.get("involuntary_context_switches")),
        "perf_ipc": perf.get("ipc"),
        "host_exclusive_partition_verified": False,
        "controlled_speedup": None,
    }
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("receipt", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = measure(json.loads(args.receipt.read_text()))
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(rendered)
    print(rendered, end="")
