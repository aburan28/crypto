"""Run the frozen S3 point pair inside one bounded Modal VM Sandbox."""

from __future__ import annotations

import base64
import ctypes
import errno
import gzip
import hashlib
import json
import os
from pathlib import Path
import platform
import signal
import struct
import subprocess
import sys
import time


BASE = Path("/opt/s3")
PHASES = ("target_query", "target_pdp_charged", "target_relation_check",
          "target_descent", "target_recovery_check")
CGROUP = ("cpu.max", "cpu.stat", "cpu.pressure", "memory.pressure",
          "io.pressure", "memory.current", "memory.max", "memory.peak",
          "memory.events", "memory.numa_stat",
          "cpuset.cpus", "cpuset.cpus.effective",
          "cpuset.cpus.exclusive.effective", "cpuset.cpus.partition",
          "cpuset.mems", "cpuset.mems.effective")


def sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def read(path: str, limit: int = 200_000) -> str | None:
    try:
        return Path(path).read_text(errors="replace")[:limit]
    except OSError:
        return None


def model() -> str | None:
    for line in (read("/proc/cpuinfo") or "").splitlines():
        if line.startswith(("model name", "Hardware", "Processor")) and ":" in line:
            return line.split(":", 1)[1].strip()
    return None


def topology(allowed: list[int]) -> dict:
    cores = {}
    for cpu in allowed:
        prefix = f"/sys/devices/system/cpu/cpu{cpu}/topology/"
        cores[str(cpu)] = {key: read(prefix + key, 80) for key in
                           ("core_id", "physical_package_id", "thread_siblings_list")}
    nodes = {path.name: read(str(path / "cpulist"), 1000) for path in
             sorted(Path("/sys/devices/system/node").glob("node[0-9]*"))}
    return {"cpus": cores, "numa_nodes": nodes}


def choose_cpus(allowed: list[int], topo: dict, count: int) -> list[int]:
    selected = []
    seen = set()
    for cpu in allowed:
        entry = topo["cpus"][str(cpu)]
        identity = (entry["physical_package_id"], entry["core_id"])
        if identity not in seen:
            seen.add(identity)
            selected.append(cpu)
        if len(selected) == count:
            return selected
    for cpu in allowed:
        if cpu not in selected:
            selected.append(cpu)
        if len(selected) == count:
            return selected
    return selected


def per_cpu_ticks(cpus: list[int]) -> dict:
    selected = set(cpus)
    result = {}
    for line in (read("/proc/stat") or "").splitlines():
        parts = line.split()
        if parts and parts[0].startswith("cpu") and parts[0][3:].isdigit():
            cpu = int(parts[0][3:])
            if cpu in selected and len(parts) >= 9:
                result[str(cpu)] = {"user": int(parts[1]), "system": int(parts[3]),
                                    "idle": int(parts[4]), "steal": int(parts[8])}
    return result


def interrupts(cpus: list[int]) -> dict | None:
    raw = read("/proc/interrupts", 1_000_000)
    if not raw:
        return None
    lines = raw.splitlines()
    columns = [int(field[3:]) for field in lines[0].split()
               if field.startswith("CPU") and field[3:].isdigit()]
    counts = {str(cpu): 0 for cpu in cpus if cpu in columns}
    for line in lines[1:]:
        parts = line.split()
        if len(parts) < len(columns) + 1:
            continue
        values = parts[1:1 + len(columns)]
        if not all(value.isdigit() for value in values):
            continue
        for cpu, value in zip(columns, values):
            if str(cpu) in counts:
                counts[str(cpu)] += int(value)
    return counts


def irq_routing(cpus: list[int]) -> dict:
    raw = read("/proc/interrupts", 1_000_000) or ""
    lines = raw.splitlines()
    if not lines:
        return {}
    columns = [int(field[3:]) for field in lines[0].split()
               if field.startswith("CPU") and field[3:].isdigit()]
    selected = set(cpus)
    ranked = []
    for line in lines[1:]:
        parts = line.split()
        if len(parts) < len(columns) + 1 or not parts[0].rstrip(":").isdigit():
            continue
        values = parts[1:1 + len(columns)]
        if not all(value.isdigit() for value in values):
            continue
        count = sum(int(value) for cpu, value in zip(columns, values)
                    if cpu in selected)
        if count:
            ranked.append((count, parts[0].rstrip(":")))
    return {irq: {"selected_cpu_interrupts": count,
                  "smp_affinity_list": read(f"/proc/irq/{irq}/smp_affinity_list", 100)}
            for count, irq in sorted(ranked, reverse=True)[:64]}


def snapshot(cpus: list[int]) -> dict:
    return {
        "unix_ns": time.time_ns(),
        "cpu_ticks": per_cpu_ticks(cpus),
        "interrupts": interrupts(cpus),
        "cpu_psi": read("/proc/pressure/cpu", 4096),
        "memory_psi": read("/proc/pressure/memory", 4096),
        "cgroup": {key: read("/sys/fs/cgroup/" + key, 4096) for key in CGROUP},
        "frequency_khz": {str(cpu): read(
            f"/sys/devices/system/cpu/cpu{cpu}/cpufreq/scaling_cur_freq", 80)
            for cpu in cpus},
    }


def perf_open(pid: int, config: int) -> tuple[int | None, str | None]:
    syscall = {"x86_64": 298, "aarch64": 241}.get(platform.machine())
    if syscall is None:
        return None, "unsupported_architecture"
    attr = ctypes.create_string_buffer(64)
    struct.pack_into("IIQ", attr, 0, 0, 64, config)
    struct.pack_into("Q", attr, 40, 2)  # inherit future threads and children
    libc = ctypes.CDLL(None, use_errno=True)
    libc.syscall.restype = ctypes.c_long
    fd = libc.syscall(syscall, ctypes.byref(attr), pid, -1, -1, 0)
    if fd < 0:
        code = ctypes.get_errno()
        return None, f"{errno.errorcode.get(code, code)}:{os.strerror(code)}"
    return int(fd), None


def perf_close(fd: int | None) -> int | None:
    if fd is None:
        return None
    try:
        raw = os.read(fd, 8)
        return struct.unpack("Q", raw)[0] if len(raw) == 8 else None
    finally:
        os.close(fd)


def process_rss(pid: int) -> int | None:
    for line in (read(f"/proc/{pid}/status", 20_000) or "").splitlines():
        if line.startswith("VmRSS:"):
            fields = line.split()
            if len(fields) >= 2 and fields[1].isdigit():
                return int(fields[1]) * 1024
    return None


def thread_affinity(pid: int) -> tuple[int, list[str]]:
    task_dir = Path(f"/proc/{pid}/task")
    masks = set()
    tids = list(task_dir.glob("[0-9]*"))
    for task in tids:
        for line in (read(str(task / "status"), 20_000) or "").splitlines():
            if line.startswith("Cpus_allowed_list:"):
                masks.add(line.split(":", 1)[1].strip())
                break
    return len(tids), sorted(masks)


def numa_pages(pid: int) -> dict:
    pages = {}
    for line in (read(f"/proc/{pid}/numa_maps", 1_000_000) or "").splitlines():
        for field in line.split():
            if field.startswith("N") and "=" in field:
                node, count = field.split("=", 1)
                if node[1:].isdigit() and count.isdigit():
                    pages[node] = pages.get(node, 0) + int(count)
    return pages


def run_one(bundle: dict, case: dict, arm: str, repetition: int,
            cpus: list[int], seconds: float) -> dict:
    frozen = case["arms"][arm]
    point = Path(f"/tmp/s3_n{case['n']}_T181.json")
    point.write_text(json.dumps(case["public_point"], separators=(",", ":")) + "\n")
    output = Path(f"/tmp/s3_n{case['n']}_{arm}_R{repetition}.json")
    binary = BASE / ("baseline" if arm == "baseline" else "expanded")
    argv = [str(binary), str(case["n"]), "0", str(frozen["columns"]),
            str(bundle["relation_seed"]), str(point), str(output),
            str(bundle["relation_workers"]), str(frozen["prepare_rank"]),
            str(bundle["target_workers"])]
    before = snapshot(cpus)
    started = time.time_ns()
    error = None
    stdout = stderr = b""
    peak_rss = 0
    observed_numa = {}
    observed_schedstat = None
    max_threads = 0
    thread_masks = set()
    polls = 0
    child = subprocess.Popen(argv, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                             start_new_session=True,
                             preexec_fn=lambda: os.sched_setaffinity(0, set(cpus)))
    actual_affinity = sorted(os.sched_getaffinity(child.pid))
    cycles_fd, cycles_error = perf_open(child.pid, 0)
    instructions_fd, instructions_error = perf_open(child.pid, 1)
    deadline = time.monotonic() + seconds
    while True:
        remaining = deadline - time.monotonic()
        if remaining <= 0:
            error = "run_timeout"
            os.killpg(child.pid, signal.SIGKILL)
            stdout, stderr = child.communicate()
            break
        try:
            stdout, stderr = child.communicate(timeout=min(0.25, remaining))
            break
        except subprocess.TimeoutExpired:
            polls += 1
            peak_rss = max(peak_rss, process_rss(child.pid) or 0)
            observed_schedstat = read(f"/proc/{child.pid}/schedstat", 200)
            if polls % 4 == 0:
                count, masks = thread_affinity(child.pid)
                max_threads = max(max_threads, count)
                thread_masks.update(masks)
            if polls % 8 == 0:
                observed_numa = numa_pages(child.pid)
    finished = time.time_ns()
    after = snapshot(cpus)
    cycles = perf_close(cycles_fd)
    instructions = perf_close(instructions_fd)
    raw = output.read_bytes() if output.exists() else b""
    try:
        report = json.loads(raw) if raw else None
    except ValueError:
        report, error = None, "invalid_json"
    timing = report.get("timing_ms", {}) if report else {}
    phase_ok = bool(all(isinstance(timing.get(k), (int, float)) and timing[k] >= 0
                        for k in PHASES)
                    and isinstance(timing.get("target_online_after_reusable_setup"),
                                   (int, float))
                    and abs(sum(timing[k] for k in PHASES) -
                            timing["target_online_after_reusable_setup"]) < .02)
    native_verified = bool(child.returncode == 0 and stdout == raw and report
                           and phase_ok and report.get("group_verified") is True
                           and report.get("target_count") == 1
                           and report.get("target") == case["public_point"]
                           and report.get("factor_base_digest") == frozen["factor_base_digest"]
                           and report.get("factor_base_points") == frozen["factor_base_points"]
                           and report.get("orbit_columns") == frozen["columns"]
                           and report.get("prepared_rank") == frozen["prepare_rank"]
                           and report.get("root_table_slots") == frozen["root_table_slots"]
                           and report.get("target_query_workers") == 1)
    return {
        "n": case["n"], "label": "T181", "arm": arm, "repetition": repetition,
        "candidate_id": frozen["candidate_id"], "workload_id": case["workload_id"],
        "run_id": f"{frozen['candidate_id']}W{case['workload_id']}R{repetition}",
        "argv": argv, "affinity": actual_affinity,
        "started_unix_ns": started, "finished_unix_ns": finished,
        "outer_wall_ns": finished - started, "exit_code": child.returncode,
        "error": error, "stderr": stderr.decode(errors="replace")[-4000:],
        "raw_sha256": sha(raw) if raw else None,
        "native_verified": native_verified, "phase_ok": phase_ok,
        "sampled_peak_rss_bytes": peak_rss,
        "sampled_numa_pages": observed_numa,
        "sampled_process_schedstat": observed_schedstat,
        "max_threads_seen": max_threads,
        "thread_affinity_masks_seen": sorted(thread_masks),
        "perf": {"cycles": cycles, "instructions": instructions,
                 "ipc": instructions / cycles if cycles and instructions else None,
                 "cycles_error": cycles_error,
                 "instructions_error": instructions_error,
                 "boundary": "child counters opened immediately after process start; inheriting later threads"},
        "before": before, "after": after, "report": report,
    }


def main() -> None:
    expected = json.loads(sys.argv[1])
    bundle_raw = (BASE / "public_bundle.json").read_bytes()
    assert sha(bundle_raw) == expected["bundle_sha256"]
    bundle = json.loads(bundle_raw)
    assert bundle["kind"] == "s3_modal_combined_public_t181_bundle_v1"
    for name in ("baseline", "expanded"):
        binary = BASE / name
        assert sha(binary.read_bytes()) == expected["binary_sha256"][name]
        binary.chmod(0o755)
    allowed = sorted(os.sched_getaffinity(0))
    topo = topology(allowed)
    cpus = choose_cpus(allowed, topo, expected["logical_execution_cpus"])
    assert len(cpus) == expected["logical_execution_cpus"]
    host = {
        "platform": platform.platform(), "cpu_model": model(),
        "cpu_count_visible": os.cpu_count(), "allowed_cpus": allowed,
        "selected_cpus": cpus, "topology": topo,
        "irq_routing_top64": irq_routing(cpus),
        "cpu_governor": {str(cpu): read(
            f"/sys/devices/system/cpu/cpu{cpu}/cpufreq/scaling_governor", 100)
            for cpu in cpus},
        "cpu_frequency_limits_khz": {str(cpu): {
            key: read(f"/sys/devices/system/cpu/cpu{cpu}/cpufreq/{key}", 100)
            for key in ("scaling_min_freq", "scaling_max_freq")}
            for cpu in cpus},
        "process_memory_allowed": [line for line in
                                   (read("/proc/self/status", 20_000) or "").splitlines()
                                   if line.startswith("Mems_allowed_list")],
        "perf_event_paranoid": read("/proc/sys/kernel/perf_event_paranoid", 80),
        "host_exclusive_cpu_partition_receipt": None,
        "controlled_speedup": None,
    }
    host["idle_before"] = snapshot(cpus)
    time.sleep(2)
    host["idle_after"] = snapshot(cpus)
    cases = {case["n"]: case for case in bundle["curves"]}
    repetitions = {(n, arm): 0 for n in cases for arm in ("baseline", "expanded")}
    rows = []
    stop_reason = None
    deadline = time.monotonic() + expected["workload_deadline_seconds"]
    for n, arm in bundle["run_order"]:
        remaining = deadline - time.monotonic()
        if remaining < expected["minimum_run_seconds"]:
            stop_reason = "workload_deadline"
            break
        repetitions[n, arm] += 1
        row = run_one(bundle, cases[n], arm, repetitions[n, arm], cpus,
                      min(expected["per_run_timeout_seconds"], remaining - 5))
        rows.append(row)
        if not row["native_verified"]:
            stop_reason = "invalid_or_timeout"
            break
    result = {
        "kind": "s3_modal_combined_pilot_guest_receipt_v1",
        "bundle_sha256": expected["bundle_sha256"],
        "binary_sha256": expected["binary_sha256"],
        "run_order": bundle["run_order"], "rows": rows,
        "stop_reason": stop_reason, "host": host,
    }
    encoded = base64.b64encode(gzip.compress(json.dumps(
        result, sort_keys=True, separators=(",", ":")).encode(), mtime=0)).decode()
    print(json.dumps({"kind": "s3_modal_combined_compressed_receipt_v1",
                      "encoding": "gzip_base64", "sha256": sha(encoded.encode()),
                      "payload": encoded}, separators=(",", ":")), flush=True)


if __name__ == "__main__":
    main()
