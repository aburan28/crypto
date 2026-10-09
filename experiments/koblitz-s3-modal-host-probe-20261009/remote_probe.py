"""Generic, source-free Linux VM CPU telemetry probe for Modal."""

from __future__ import annotations

import ctypes
import errno
import json
import os
from pathlib import Path
import platform
import resource
import struct
import time


def read(path: str, limit: int = 2048) -> str | None:
    try:
        return Path(path).read_text(errors="replace")[:limit]
    except OSError:
        return None


def cpu_list(spec: str | None) -> list[int]:
    result: set[int] = set()
    for item in (spec or "").strip().split(","):
        if not item:
            continue
        ends = item.split("-", 1)
        if all(value.isdigit() for value in ends):
            result.update(range(int(ends[0]), int(ends[-1]) + 1))
    return sorted(result)


def cpu_model() -> str | None:
    raw = read("/proc/cpuinfo", 50_000) or ""
    for line in raw.splitlines():
        if line.startswith(("model name", "Hardware", "Processor")) and ":" in line:
            return line.split(":", 1)[1].strip()
    return None


def cgroup_files() -> dict:
    root = Path("/sys/fs/cgroup")
    names = ("cpu.max", "cpu.stat", "cpu.pressure", "memory.pressure",
             "memory.current", "memory.max", "memory.events",
             "cpuset.cpus", "cpuset.cpus.effective", "cpuset.cpus.exclusive.effective",
             "cpuset.cpus.partition", "cpuset.mems", "cpuset.mems.effective")
    return {name: read(str(root / name)) for name in names}


def numa_maps_summary() -> dict:
    """Summarize process memory policy and observed NUMA page placement."""
    raw = read("/proc/self/numa_maps", 1_000_000)
    if raw is None:
        return {"available": False}
    policies: dict[str, int] = {}
    node_pages: dict[str, int] = {}
    for line in raw.splitlines():
        fields = line.split()
        if len(fields) < 2:
            continue
        policies[fields[1]] = policies.get(fields[1], 0) + 1
        for field in fields[2:]:
            if field.startswith("N") and "=" in field:
                node, count = field.split("=", 1)
                if node[1:].isdigit() and count.isdigit():
                    node_pages[node] = node_pages.get(node, 0) + int(count)
    return {"available": True, "policies": policies, "node_pages": node_pages}


def steal_ticks(cpus: list[int]) -> dict[str, int]:
    raw = read("/proc/stat", 200_000) or ""
    selected = set(cpus)
    result = {}
    for line in raw.splitlines():
        fields = line.split()
        if fields and fields[0].startswith("cpu") and fields[0][3:].isdigit():
            cpu = int(fields[0][3:])
            if cpu in selected and len(fields) > 8:
                result[str(cpu)] = int(fields[8])
    return result


def interrupt_totals(cpus: list[int]) -> dict[str, int] | None:
    raw = read("/proc/interrupts", 1_000_000)
    if not raw:
        return None
    lines = raw.splitlines()
    if not lines:
        return None
    headers = [int(item[3:]) for item in lines[0].split()
               if item.startswith("CPU") and item[3:].isdigit()]
    totals = {str(cpu): 0 for cpu in cpus if cpu in headers}
    for line in lines[1:]:
        parts = line.split()
        if len(parts) < len(headers) + 1:
            continue
        counts = parts[1:1 + len(headers)]
        if not all(item.isdigit() for item in counts):
            continue
        for cpu, count in zip(headers, counts):
            if str(cpu) in totals:
                totals[str(cpu)] += int(count)
    return totals


def snapshot(cpus: list[int]) -> dict:
    usage = resource.getrusage(resource.RUSAGE_SELF)
    return {
        "unix_ns": time.time_ns(),
        "affinity": sorted(os.sched_getaffinity(0)),
        "cgroup": cgroup_files(),
        "system_cpu_psi": read("/proc/pressure/cpu"),
        "system_memory_psi": read("/proc/pressure/memory"),
        "steal_ticks": steal_ticks(cpus),
        "interrupt_totals": interrupt_totals(cpus),
        "frequency_khz": {str(cpu): read(
            f"/sys/devices/system/cpu/cpu{cpu}/cpufreq/scaling_cur_freq")
            for cpu in cpus[:64]},
        "involuntary_context_switches": usage.ru_nivcsw,
        "voluntary_context_switches": usage.ru_nvcsw,
        "max_rss_raw": usage.ru_maxrss,
        "process_schedstat": read("/proc/self/schedstat"),
        "process_numa_maps": numa_maps_summary(),
    }


def perf_open(config: int) -> tuple[int | None, str | None]:
    syscall = {"x86_64": 298, "aarch64": 241}.get(platform.machine())
    if syscall is None:
        return None, "unknown_perf_event_open_syscall_for_architecture"
    attr = ctypes.create_string_buffer(64)
    struct.pack_into("IIQ", attr, 0, 0, 64, config)
    libc = ctypes.CDLL(None, use_errno=True)
    libc.syscall.restype = ctypes.c_long
    fd = libc.syscall(syscall, ctypes.byref(attr), 0, -1, -1, 0)
    if fd < 0:
        code = ctypes.get_errno()
        return None, f"{errno.errorcode.get(code, code)}:{os.strerror(code)}"
    return int(fd), None


def perf_read(fd: int | None) -> int | None:
    if fd is None:
        return None
    try:
        raw = os.read(fd, 8)
        return struct.unpack("Q", raw)[0] if len(raw) == 8 else None
    finally:
        os.close(fd)


def main() -> None:
    allowed = sorted(os.sched_getaffinity(0))
    assert allowed
    topology = {}
    for cpu in allowed[:128]:
        root = f"/sys/devices/system/cpu/cpu{cpu}/topology"
        topology[str(cpu)] = {
            "thread_siblings_list": read(root + "/thread_siblings_list"),
            "core_id": read(root + "/core_id"),
            "physical_package_id": read(root + "/physical_package_id"),
        }
    nodes = {}
    for path in sorted(Path("/sys/devices/system/node").glob("node[0-9]*")):
        nodes[path.name] = read(str(path / "cpulist"))
    before = snapshot(allowed)
    time.sleep(3)
    after_idle = snapshot(allowed)
    pin_error = None
    selected_cpu = allowed[0]
    try:
        os.sched_setaffinity(0, {selected_cpu})
    except OSError as exc:
        pin_error = f"{type(exc).__name__}:{exc.errno}"
    after_pin = snapshot(allowed)
    cycles_fd, cycles_error = perf_open(0)
    instructions_fd, instructions_error = perf_open(1)
    started = time.perf_counter_ns()
    value = 1
    while time.perf_counter_ns() - started < 5_000_000_000:
        value = (value * 1664525 + 1013904223) & ((1 << 64) - 1)
    busy_elapsed_ns = time.perf_counter_ns() - started
    cycles = perf_read(cycles_fd)
    instructions = perf_read(instructions_fd)
    after_busy = snapshot(allowed)
    try:
        os.sched_setaffinity(0, set(allowed))
    except OSError:
        pass
    result = {
        "kind": "source_free_modal_vm_cpu_fidelity_probe_v1",
        "platform": platform.platform(), "machine": platform.machine(),
        "cpu_model": cpu_model(), "cpu_count_visible": os.cpu_count(),
        "allowed_cpus": allowed, "selected_cpu": selected_cpu,
        "pin_error": pin_error, "topology": topology,
        "numa_nodes": nodes,
        "proc_self_status_cpu_mems": [line for line in
            (read("/proc/self/status", 20_000) or "").splitlines()
            if line.startswith(("Cpus_allowed_list", "Mems_allowed_list"))],
        "nohz_full": read("/sys/devices/system/cpu/nohz_full"),
        "isolated_cpus": read("/sys/devices/system/cpu/isolated"),
        "perf_event_paranoid": read("/proc/sys/kernel/perf_event_paranoid"),
        "perf": {"cycles": cycles, "instructions": instructions,
                 "ipc": instructions / cycles if cycles and instructions else None,
                 "cycles_error": cycles_error,
                 "instructions_error": instructions_error},
        "before": before, "after_idle": after_idle,
        "after_pin": after_pin, "after_busy": after_busy,
        "busy_elapsed_ns": busy_elapsed_ns, "busy_checksum": value,
        "host_exclusive_partition_verified": False,
        "controlled_speedup": None,
    }
    print(json.dumps(result, sort_keys=True, separators=(",", ":")), flush=True)


if __name__ == "__main__":
    main()
