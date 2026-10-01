"""The host capsule: everything about the machine a measurement ran on.

Two parts, with different jobs.

``stable``
    Facts that must be identical for two wall-clock measurements to be paired:
    CPU model, microcode, feature flags, cache and core topology, the kernel
    release and boot command line, frequency policy (governor, turbo/boost),
    SMT state, transparent huge pages, the sysctls in ``SYSCTLS``,
    virtualisation, cgroup limits and the toolchain.  The subset in
    ``env_class()`` is hashed into ``env_class_id``.

``volatile``
    Facts about the moment: load, pressure stall information (PSI), per-CPU
    jiffies including hypervisor steal, frequencies and temperatures.  Never
    hashed; sampled before and after every run.

Every probe returns ``None`` for something the kernel does not expose.  An
absent fact is unknown, never a pass: the isolation gate treats ``None`` as a
failed check for the level that needs it.

AGENTS.md section 10 asks every baseline to record the commit, ``rustc
--version``, the CPU model and its relevant features (popcnt, avx2, avx512*,
pclmulqdq, NEON/aes), the core count, memory, OS, architecture and the exact
command.  ``stable_capsule()`` records all of these and more; the command is
in the run record.
"""
from __future__ import annotations

import glob
import os
import platform
import re
import shutil
import subprocess
import sys
import time
from typing import Any

from .canonical import env_class_id as _env_class_id, sha256_bytes, sha256_file

CAPSULE_SCHEMA = "icms.capsule/v1"

SYSCTLS = (
    "kernel/randomize_va_space",
    "kernel/nmi_watchdog",
    "kernel/watchdog",
    "kernel/perf_event_paranoid",
    "kernel/sched_autogroup_enabled",
    "kernel/numa_balancing",
    "kernel/sched_rt_runtime_us",
    "kernel/sched_rt_period_us",
    "kernel/timer_migration",
    "kernel/sched_child_runs_first",
    "kernel/sched_util_clamp_min",
    "kernel/sched_schedstats",
    "vm/swappiness",
    "vm/overcommit_memory",
    "vm/zone_reclaim_mode",
    "vm/stat_interval",
    "vm/dirty_ratio",
    "vm/dirty_background_ratio",
    "vm/nr_hugepages",
    "vm/compaction_proactiveness",
)

CMDLINE_KEYS = (
    "isolcpus", "nohz_full", "rcu_nocbs", "irqaffinity", "nosmt", "mitigations",
    "intel_pstate", "amd_pstate", "processor.max_cstate", "intel_idle.max_cstate",
    "idle", "psi", "transparent_hugepage", "skew_tick", "tsc", "clocksource",
    "nmi_watchdog", "nowatchdog", "cpufreq.default_governor",
)

# The features AGENTS.md section 10 names, plus the ones the binary-field and
# big-integer code dispatches on.  Recorded explicitly; the whole flag set is
# recorded as a hash.
FEATURES = (
    "popcnt", "bmi1", "bmi2", "adx", "avx", "avx2", "fma", "avx512f", "avx512bw",
    "avx512dq", "avx512vl", "avx512ifma", "avx512vbmi", "vpclmulqdq", "pclmulqdq",
    "gfni", "aes", "vaes", "sha_ni", "rdrand", "hypervisor",
    # Arm64 names as /proc/cpuinfo spells them
    "asimd", "pmull", "sha2", "sve", "sve2",
)


def _read(path: str) -> str | None:
    try:
        with open(path, "r", encoding="utf-8", errors="replace") as fh:
            return fh.read().strip()
    except OSError:
        return None


def _run(argv: list[str], timeout: float = 15.0, cwd: str | None = None) -> str | None:
    exe = argv[0] if os.path.isabs(argv[0]) else shutil.which(argv[0])
    if exe is None or not os.path.exists(exe):
        return None
    try:
        out = subprocess.run([exe, *argv[1:]], capture_output=True, text=True,
                             timeout=timeout, check=False, cwd=cwd)
    except (OSError, subprocess.SubprocessError):
        return None
    text = (out.stdout or "").strip() or (out.stderr or "").strip()
    return text or None


def parse_cpu_list(spec: str | None) -> list[int]:
    out: list[int] = []
    if not spec:
        return out
    for part in spec.replace("\n", ",").split(","):
        part = part.strip()
        if not part:
            continue
        if "-" in part:
            a, b = part.split("-", 1)
            out.extend(range(int(a), int(b) + 1))
        else:
            out.append(int(part))
    return sorted(set(out))


# --------------------------------------------------------------- stable ------

def cpu_info() -> dict[str, Any]:
    text = _read("/proc/cpuinfo") or ""
    first: dict[str, str] = {}
    for line in text.splitlines():
        if not line.strip():
            if first:
                break
            continue
        if ":" in line:
            k, v = line.split(":", 1)
            first.setdefault(k.strip(), v.strip())
    flags = sorted(set((first.get("flags") or first.get("Features") or "").split()))
    return {
        "vendor_id": first.get("vendor_id") or first.get("CPU implementer"),
        "model_name": first.get("model name") or first.get("Model") or platform.processor() or None,
        "family": first.get("cpu family"),
        "model": first.get("model") or first.get("CPU part"),
        "stepping": first.get("stepping") or first.get("CPU revision"),
        "microcode": first.get("microcode"),
        "cache_size": first.get("cache size"),
        "flags_count": len(flags),
        "flags_sha256": sha256_bytes(" ".join(flags).encode()) if flags else None,
        "features": {f: (f in flags) for f in FEATURES} if flags else None,
        "logical_cpus": os.cpu_count(),
    }


def topology() -> dict[str, Any]:
    cpus: dict[str, Any] = {}
    for d in sorted(glob.glob("/sys/devices/system/cpu/cpu[0-9]*"), key=lambda p: int(re.search(r"(\d+)$", p).group(1))):
        n = re.search(r"cpu(\d+)$", d).group(1)
        online = _read(f"{d}/online")
        cpus[n] = {
            "online": online in (None, "1"),
            "core_id": _read(f"{d}/topology/core_id"),
            "package_id": _read(f"{d}/topology/physical_package_id"),
            "thread_siblings": _read(f"{d}/topology/thread_siblings_list"),
            "governor": _read(f"{d}/cpufreq/scaling_governor"),
            "driver": _read(f"{d}/cpufreq/scaling_driver"),
            "min_khz": _read(f"{d}/cpufreq/scaling_min_freq"),
            "max_khz": _read(f"{d}/cpufreq/scaling_max_freq"),
            "epp": _read(f"{d}/cpufreq/energy_performance_preference"),
        }
    caches = []
    for d in sorted(glob.glob("/sys/devices/system/cpu/cpu0/cache/index*")):
        caches.append({k: _read(f"{d}/{k}") for k in ("level", "type", "size", "ways_of_associativity", "shared_cpu_list")})
    return {
        "online": _read("/sys/devices/system/cpu/online"),
        "isolated": _read("/sys/devices/system/cpu/isolated"),
        "nohz_full": _read("/sys/devices/system/cpu/nohz_full"),
        "smt_control": _read("/sys/devices/system/cpu/smt/control"),
        "smt_active": _read("/sys/devices/system/cpu/smt/active"),
        "intel_pstate_status": _read("/sys/devices/system/cpu/intel_pstate/status"),
        "intel_pstate_no_turbo": _read("/sys/devices/system/cpu/intel_pstate/no_turbo"),
        "cpufreq_boost": _read("/sys/devices/system/cpu/cpufreq/boost"),
        "cpuidle_driver": _read("/sys/devices/system/cpu/cpuidle/current_driver"),
        "numa_nodes": sorted(os.path.basename(p) for p in glob.glob("/sys/devices/system/node/node[0-9]*")),
        "caches_cpu0": caches,
        "cpus": cpus,
    }


def kernel() -> dict[str, Any]:
    u = platform.uname()
    cmdline = _read("/proc/cmdline") or ""
    flags: dict[str, Any] = {}
    for tok in cmdline.split():
        if tok == "--":
            break
        key = tok.split("=", 1)[0]
        if key in CMDLINE_KEYS:
            flags[key] = tok.split("=", 1)[1] if "=" in tok else True
    vulns = {os.path.basename(p): _read(p) for p in sorted(glob.glob("/sys/devices/system/cpu/vulnerabilities/*"))}
    return {
        "sysname": u.system,
        "release": u.release,
        "version": u.version,
        "machine": u.machine,
        "cmdline": cmdline,
        "cmdline_sha256": sha256_bytes(cmdline.encode()),
        "cmdline_flags": flags,
        "clocksource": _read("/sys/devices/system/clocksource/clocksource0/current_clocksource"),
        "clocksources_available": _read("/sys/devices/system/clocksource/clocksource0/available_clocksource"),
        "sysctl": {k: _read(f"/proc/sys/{k}") for k in SYSCTLS},
        "thp_enabled": _read("/sys/kernel/mm/transparent_hugepage/enabled"),
        "thp_defrag": _read("/sys/kernel/mm/transparent_hugepage/defrag"),
        "psi_available": os.path.exists("/proc/pressure/cpu"),
        "vulnerabilities": vulns,
    }


def virtualization() -> dict[str, Any]:
    return {
        "detect_virt": _run(["systemd-detect-virt"]),
        "detect_vm": _run(["systemd-detect-virt", "--vm"]),
        "detect_container": _run(["systemd-detect-virt", "--container"]),
        "dmi_sys_vendor": _read("/sys/class/dmi/id/sys_vendor"),
        "dmi_product_name": _read("/sys/class/dmi/id/product_name"),
        "dockerenv": os.path.exists("/.dockerenv"),
        "containerenv": os.path.exists("/run/.containerenv"),
    }


def cgroup() -> dict[str, Any]:
    limits = {}
    for path in ("/sys/fs/cgroup/cpu.max", "/sys/fs/cgroup/cpuset.cpus.effective",
                 "/sys/fs/cgroup/memory.max", "/sys/fs/cgroup/cpu/cpu.cfs_quota_us",
                 "/sys/fs/cgroup/cpu/cpu.cfs_period_us", "/sys/fs/cgroup/cpuset/cpuset.cpus",
                 "/sys/fs/cgroup/memory/memory.limit_in_bytes"):
        v = _read(path)
        if v is not None:
            limits[path] = v
    try:
        affinity = sorted(os.sched_getaffinity(0))
    except (AttributeError, OSError):
        affinity = None
    return {"self_cgroup": (_read("/proc/self/cgroup") or "").splitlines(), "limits": limits,
            "process_affinity": affinity}


def memory() -> dict[str, Any]:
    info = {}
    for line in (_read("/proc/meminfo") or "").splitlines():
        k, _, v = line.partition(":")
        if k in ("MemTotal", "SwapTotal", "HugePages_Total", "Hugepagesize"):
            info[k] = v.strip()
    return info


def os_release() -> dict[str, Any]:
    out = {}
    for line in (_read("/etc/os-release") or "").splitlines():
        if "=" in line:
            k, v = line.split("=", 1)
            if k in ("ID", "VERSION_ID", "PRETTY_NAME"):
                out[k] = v.strip('"')
    out["libc"] = "-".join(x for x in platform.libc_ver() if x) or None
    return out


def toolchain(extra: dict[str, list[str]] | None = None) -> dict[str, Any]:
    probes = {
        "rustc": ["rustc", "-Vv"],
        "cargo": ["cargo", "-V"],
        "gcc": ["gcc", "--version"],
        "clang": ["clang", "--version"],
        "cc": ["cc", "--version"],
        "valgrind": ["valgrind", "--version"],
        "msolve": ["msolve", "-V"],
    }
    probes.update(extra or {})
    out: dict[str, Any] = {
        "python": sys.version.split()[0],
        "python_implementation": platform.python_implementation(),
        "python_executable": sys.executable,
        "python_executable_sha256": sha256_file(sys.executable),
    }
    for mod in ("yaml", "numpy"):
        try:
            out["pyyaml" if mod == "yaml" else mod] = __import__(mod).__version__
        except ImportError:
            out["pyyaml" if mod == "yaml" else mod] = None
    # Build flags the measured binaries inherit (a -march=native kernel or a
    # RUSTFLAGS target-cpu changes the code that runs, not just its speed).
    out["build_env"] = {k: os.environ.get(k) for k in ("RUSTFLAGS", "CARGO_BUILD_RUSTFLAGS", "CARGO_PROFILE_RELEASE_LTO",
                                                       "CC", "CFLAGS", "CXX", "CXXFLAGS", "LDFLAGS")}
    for name, argv in probes.items():
        text = _run(argv)
        out[name] = text.splitlines()[0] if text else None
        if name == "rustc" and text:
            out["rustc_verbose"] = text
    return out


def git_state(path: str) -> dict[str, Any] | None:
    top = _run(["git", "-C", path, "rev-parse", "--show-toplevel"])
    if not top:
        return None
    head = _run(["git", "-C", top, "rev-parse", "HEAD"])
    status = _run(["git", "-C", top, "status", "--porcelain", "--untracked-files=no"]) or ""
    diff = subprocess.run(["git", "-C", top, "diff", "HEAD"], capture_output=True).stdout if status else b""
    remote = _run(["git", "-C", top, "config", "--get", "remote.origin.url"])
    return {"toplevel": top, "commit": head, "dirty": bool(status.strip()),
            "dirty_paths": status.splitlines()[:50], "diff_sha256": sha256_bytes(diff) if diff else None,
            "remote": remote}


def stable_capsule(extra_tools: dict[str, list[str]] | None = None) -> dict[str, Any]:
    return {
        "schema": CAPSULE_SCHEMA,
        "captured_unix": time.time(),
        "hostname_sha256": sha256_bytes(platform.node().encode()),
        "cpu": cpu_info(),
        "topology": topology(),
        "kernel": kernel(),
        "virtualization": virtualization(),
        "cgroup": cgroup(),
        "memory": memory(),
        "os": os_release(),
        "toolchain": toolchain(extra_tools),
    }


def env_class(stable: dict[str, Any]) -> dict[str, Any]:
    """The facts two wall-clock measurements must share to be paired."""
    t, k, c = stable["topology"], stable["kernel"], stable["cpu"]
    governors = sorted({str(v.get("governor")) for v in t["cpus"].values()}) if t["cpus"] else ["None"]
    return {
        "cpu": {x: c.get(x) for x in ("vendor_id", "model_name", "family", "model", "stepping", "microcode",
                                       "flags_sha256", "logical_cpus")},
        "frequency": {"governors": governors, "no_turbo": t["intel_pstate_no_turbo"],
                      "boost": t["cpufreq_boost"], "pstate": t["intel_pstate_status"]},
        "smt": [t["smt_control"], t["smt_active"]],
        "cpus": {"online": t["online"], "isolated": t["isolated"], "nohz_full": t["nohz_full"]},
        "kernel": {"release": k["release"], "cmdline_flags": k["cmdline_flags"], "clocksource": k["clocksource"],
                   "sysctl": k["sysctl"], "thp": [k["thp_enabled"], k["thp_defrag"]]},
        "virtualization": stable["virtualization"]["detect_virt"],
        "memory_total": stable["memory"].get("MemTotal"),
        "os": {x: stable["os"].get(x) for x in ("ID", "VERSION_ID", "libc")},
        "toolchain": {x: stable["toolchain"].get(x) for x in ("python", "rustc", "gcc", "clang", "cc", "numpy", "build_env")},
    }


def capsule(extra_tools: dict[str, list[str]] | None = None) -> dict[str, Any]:
    st = stable_capsule(extra_tools)
    cls = env_class(st)
    return {"schema": CAPSULE_SCHEMA, "env_class_id": _env_class_id(cls), "env_class": cls, "stable": st}


# --------------------------------------------------------------- volatile ----

def psi() -> dict[str, Any]:
    out: dict[str, Any] = {}
    for res in ("cpu", "memory", "io"):
        text = _read(f"/proc/pressure/{res}")
        if text is None:
            out[res] = None
            continue
        res_out = {}
        for line in text.splitlines():
            kind, *pairs = line.split()
            res_out[kind] = {p.split("=")[0]: float(p.split("=")[1]) for p in pairs}
        out[res] = res_out
    return out


JIFFY_FIELDS = ("user", "nice", "system", "idle", "iowait", "irq", "softirq", "steal", "guest", "guest_nice")


def proc_stat() -> dict[str, Any]:
    out: dict[str, Any] = {}
    for line in (_read("/proc/stat") or "").splitlines():
        if line.startswith("cpu"):
            name, *vals = line.split()
            out[name] = dict(zip(JIFFY_FIELDS, (int(v) for v in vals)))
        elif line.startswith(("ctxt ", "procs_running ", "procs_blocked ")):
            k, v = line.split()
            out[k] = int(v)
    return out


def volatile_snapshot() -> dict[str, Any]:
    return {
        "t_unix": time.time(),
        "t_monotonic_ns": time.monotonic_ns(),
        "loadavg": _read("/proc/loadavg"),
        "psi": psi(),
        "proc_stat": proc_stat(),
        "cur_khz": {re.search(r"cpu(\d+)", p).group(1): _read(p)
                    for p in sorted(glob.glob("/sys/devices/system/cpu/cpu[0-9]*/cpufreq/scaling_cur_freq"))},
        "thermal_mC": {os.path.basename(os.path.dirname(p)): _read(p)
                       for p in sorted(glob.glob("/sys/class/thermal/thermal_zone*/temp"))},
    }
