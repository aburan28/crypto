"""Host inventory: what this machine is, read from /sys and /proc.

Every reader takes the roots it reads from (``Paths``) and the process
runner it shells out with, so the parsers run unchanged on a synthetic tree
in the tests and on a live host in the worker. Nothing here changes the
host; see :mod:`isolab.isolation` for that.
"""
from __future__ import annotations

import json
import os
import platform
import re
import shutil
import socket
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable

Runner = Callable[[list[str], float], tuple[int, str]]


@dataclass(frozen=True)
class Paths:
    sys: str = "/sys"
    proc: str = "/proc"
    etc: str = "/etc"


def default_run(argv: list[str], timeout: float = 20.0) -> tuple[int, str]:
    try:
        p = subprocess.run(argv, capture_output=True, text=True, timeout=timeout)
    except (OSError, subprocess.TimeoutExpired) as err:
        return 127, str(err)
    return p.returncode, p.stdout + (p.stderr if p.returncode else "")


def read_text(path: str | Path) -> str | None:
    try:
        return Path(path).read_text().strip()
    except OSError:
        return None


def read_int(path: str | Path) -> int | None:
    text = read_text(path)
    if text is None:
        return None
    try:
        return int(text.split()[0])
    except (ValueError, IndexError):
        return None


def parse_cpulist(text: str | None) -> list[int]:
    cpus: set[int] = set()
    for part in (text or "").replace("\n", "").split(","):
        part = part.strip()
        if not part:
            continue
        if "-" in part:
            lo, hi = part.split("-", 1)
            cpus.update(range(int(lo), int(hi) + 1))
        else:
            cpus.add(int(part))
    return sorted(cpus)


def format_cpulist(cpus) -> str:
    cpus = sorted(set(cpus))
    if not cpus:
        return ""
    out, start, prev = [], cpus[0], cpus[0]
    for c in cpus[1:]:
        if c == prev + 1:
            prev = c
            continue
        out.append(f"{start}-{prev}" if start != prev else str(start))
        start = prev = c
    out.append(f"{start}-{prev}" if start != prev else str(start))
    return ",".join(out)


def kernel_version(release: str | None = None) -> tuple[int, ...]:
    release = release or platform.release()
    m = re.match(r"(\d+)\.(\d+)(?:\.(\d+))?", release)
    return tuple(int(x) for x in m.groups() if x is not None) if m else (0,)


# -- topology ----------------------------------------------------------------

def cpu_topology(paths: Paths = Paths()) -> dict[str, Any]:
    base = Path(paths.sys) / "devices/system/cpu"
    present = parse_cpulist(read_text(base / "present")) or list(range(os.cpu_count() or 1))
    online = set(parse_cpulist(read_text(base / "online")) or present)
    nodes: dict[int, dict[str, Any]] = {}
    node_base = Path(paths.sys) / "devices/system/node"
    cpu_node: dict[int, int] = {}
    if node_base.is_dir():
        for nd in sorted(node_base.glob("node[0-9]*"), key=lambda p: int(p.name[4:])):
            nid = int(nd.name[4:])
            ncpus = parse_cpulist(read_text(nd / "cpulist"))
            mem = {}
            for line in (read_text(nd / "meminfo") or "").splitlines():
                m = re.match(r"Node \d+ (\w+):\s+(\d+) kB", line)
                if m:
                    mem[m.group(1)] = int(m.group(2))
            dist = read_text(nd / "distance")
            nodes[nid] = {"cpus": ncpus, "mem_total_kb": mem.get("MemTotal"),
                          "mem_free_kb": mem.get("MemFree"),
                          "distances": [int(x) for x in dist.split()] if dist else None}
            for c in ncpus:
                cpu_node[c] = nid
    if not nodes:
        nodes[0] = {"cpus": present, "mem_total_kb": None, "mem_free_kb": None, "distances": None}
        cpu_node = {c: 0 for c in present}
    cpus: dict[int, dict[str, Any]] = {}
    for c in present:
        d = base / f"cpu{c}" / "topology"
        sib = parse_cpulist(read_text(d / "thread_siblings_list")) or [c]
        cpus[c] = {"core_id": read_int(d / "core_id") if (d / "core_id").exists() else c,
                   "package_id": read_int(d / "physical_package_id") or 0,
                   "node": cpu_node.get(c, 0), "siblings": sib, "online": c in online}
    smt_control = read_text(base / "smt/control")
    caches = []
    for idx in sorted((base / "cpu0/cache").glob("index[0-9]*")) if (base / "cpu0/cache").is_dir() else []:
        caches.append({"level": read_int(idx / "level"), "type": read_text(idx / "type"),
                       "size": read_text(idx / "size"),
                       "shared_cpu_list": read_text(idx / "shared_cpu_list")})
    cores = {(v["package_id"], v["core_id"]) for v in cpus.values()}
    packages = {v["package_id"] for v in cpus.values()}
    threads = max((len(v["siblings"]) for v in cpus.values()), default=1)
    return {"cpus": cpus, "nodes": nodes, "present": present, "online": sorted(online),
            "smt_control": smt_control, "caches": caches,
            "summary": {"sockets": len(packages), "cores": len(cores),
                        "threads_per_core": threads, "numa_nodes": len(nodes)}}


def node_of(topology: dict[str, Any], cpu: int) -> int:
    return topology["cpus"].get(cpu, {}).get("node", 0)


def cores_of(topology: dict[str, Any], cpus) -> dict[tuple[int, int], list[int]]:
    out: dict[tuple[int, int], list[int]] = {}
    for c in sorted(cpus):
        info = topology["cpus"].get(c)
        if not info:
            continue
        out.setdefault((info["package_id"], info["core_id"]), []).append(c)
    return out


# -- cpu / kernel facts ------------------------------------------------------

def cpuinfo(paths: Paths = Paths()) -> dict[str, Any]:
    text = read_text(Path(paths.proc) / "cpuinfo") or ""
    first = text.split("\n\n", 1)[0]
    fields: dict[str, str] = {}
    for line in first.splitlines():
        if ":" in line:
            k, v = line.split(":", 1)
            fields.setdefault(k.strip(), v.strip())
    flags = (fields.get("flags") or fields.get("Features") or "").split()
    model = fields.get("model name") or fields.get("Model name")
    if not model:
        dt = read_text(Path(paths.sys) / "firmware/devicetree/base/model")
        if dt:
            model = dt.replace("\x00", "")
        elif fields.get("CPU implementer"):
            model = f"implementer {fields['CPU implementer']} part {fields.get('CPU part')}"
    return {"model": model, "vendor": fields.get("vendor_id"), "flags": flags,
            "microcode": fields.get("microcode"),
            "hypervisor_flag": "hypervisor" in flags,
            "logical_cpus": text.count("processor\t:") or text.count("processor :") or os.cpu_count()}


_ISOL_FLAGS = {"domain", "nohz", "managed_irq"}


def kernel_cmdline(paths: Paths = Paths()) -> dict[str, Any]:
    raw = read_text(Path(paths.proc) / "cmdline") or ""
    out: dict[str, Any] = {"raw": raw, "isolcpus": [], "isolcpus_flags": [], "nohz_full": [],
                           "rcu_nocbs": [], "mitigations": None, "params": {}}
    for tok in raw.split():
        k, _, v = tok.partition("=")
        out["params"][k] = v
        if k == "isolcpus":
            parts = v.split(",")
            flags = [p for p in parts if p in _ISOL_FLAGS]
            rest = ",".join(p for p in parts if p not in _ISOL_FLAGS)
            out["isolcpus_flags"] = flags
            try:
                out["isolcpus"] = parse_cpulist(rest)
            except ValueError:
                out["isolcpus"] = []
        elif k in ("nohz_full", "rcu_nocbs"):
            try:
                out[k] = parse_cpulist(v)
            except ValueError:
                pass
        elif k == "mitigations":
            out["mitigations"] = v
    return out


def cpufreq(paths: Paths, cpus: list[int]) -> dict[str, Any]:
    base = Path(paths.sys) / "devices/system/cpu"
    c0 = base / f"cpu{cpus[0] if cpus else 0}" / "cpufreq"
    if not c0.is_dir():
        return {"available": False, "driver": None, "governors": {}, "cur_khz": {},
                "min_khz": None, "max_khz": None, "turbo": {"control": None, "enabled": None}}
    governors, cur = {}, {}
    for c in cpus:
        d = base / f"cpu{c}" / "cpufreq"
        governors[c] = read_text(d / "scaling_governor")
        cur[c] = read_int(d / "scaling_cur_freq")
    turbo: dict[str, Any] = {"control": None, "enabled": None}
    no_turbo = base / "intel_pstate/no_turbo"
    boost = base / "cpufreq/boost"
    if no_turbo.exists():
        turbo = {"control": str(no_turbo), "enabled": read_int(no_turbo) == 0, "inverted": True}
    elif boost.exists():
        turbo = {"control": str(boost), "enabled": read_int(boost) == 1, "inverted": False}
    return {"available": True, "driver": read_text(c0 / "scaling_driver"),
            "governors": governors, "cur_khz": cur,
            "available_governors": (read_text(c0 / "scaling_available_governors") or "").split(),
            "min_khz": read_int(c0 / "cpuinfo_min_freq"), "max_khz": read_int(c0 / "cpuinfo_max_freq"),
            "turbo": turbo}


def memory(paths: Paths = Paths()) -> dict[str, Any]:
    out: dict[str, Any] = {}
    for line in (read_text(Path(paths.proc) / "meminfo") or "").splitlines():
        k, _, v = line.partition(":")
        try:
            out[k.strip()] = int(v.split()[0])
        except (ValueError, IndexError):
            pass
    return {"total_kb": out.get("MemTotal"), "available_kb": out.get("MemAvailable"),
            "swap_total_kb": out.get("SwapTotal"), "swap_free_kb": out.get("SwapFree"),
            "hugepages_total": out.get("HugePages_Total")}


def _bracketed(text: str | None) -> str | None:
    if not text:
        return None
    m = re.search(r"\[(\w+)\]", text)
    return m.group(1) if m else text.split()[0]


def kernel_knobs(paths: Paths = Paths()) -> dict[str, Any]:
    sysk = Path(paths.proc) / "sys/kernel"
    sysvm = Path(paths.proc) / "sys/vm"
    thp = Path(paths.sys) / "kernel/mm/transparent_hugepage"
    return {"thp_enabled": _bracketed(read_text(thp / "enabled")),
            "thp_defrag": _bracketed(read_text(thp / "defrag")),
            "aslr": read_int(sysk / "randomize_va_space"),
            "perf_event_paranoid": read_int(sysk / "perf_event_paranoid"),
            "numa_balancing": read_int(sysk / "numa_balancing"),
            "nmi_watchdog": read_int(sysk / "nmi_watchdog"),
            "timer_migration": read_int(sysk / "timer_migration"),
            "sched_autogroup": read_int(sysk / "sched_autogroup_enabled"),
            "ksm_run": read_int(Path(paths.sys) / "kernel/mm/ksm/run"),
            "swappiness": read_int(sysvm / "swappiness"),
            "stat_interval": read_int(sysvm / "stat_interval")}


def cgroups(paths: Paths = Paths()) -> dict[str, Any]:
    root = Path(paths.sys) / "fs/cgroup"
    controllers = (read_text(root / "cgroup.controllers") or "").split()
    v2 = bool(controllers) or (root / "cgroup.controllers").exists()
    kv = kernel_version()
    return {"v2": v2, "root": str(root), "controllers": controllers,
            "subtree_control": (read_text(root / "cgroup.subtree_control") or "").split(),
            "cpuset_partition": v2 and "cpuset" in controllers and kv >= (5, 15),
            "isolated_partition": v2 and "cpuset" in controllers and kv >= (6, 1)}


def virtualization(paths: Paths = Paths(), run: Runner = default_run) -> dict[str, Any]:
    out: dict[str, Any] = {"detect": None, "hypervisor": None, "is_vm": False}
    if shutil.which("systemd-detect-virt"):
        rc, text = run(["systemd-detect-virt"], 5)
        out["detect"] = text.strip().splitlines()[0] if text.strip() else None
        out["is_vm"] = rc == 0 and out["detect"] not in (None, "none")
    hv = read_text(Path(paths.sys) / "hypervisor/type")
    if hv:
        out["hypervisor"] = hv
        out["is_vm"] = True
    dmi = Path(paths.sys) / "class/dmi/id"
    out["dmi"] = {"vendor": read_text(dmi / "sys_vendor"), "product": read_text(dmi / "product_name"),
                  "bios": read_text(dmi / "bios_version")}
    if cpuinfo(paths)["hypervisor_flag"]:
        out["is_vm"] = True
        out.setdefault("hypervisor", "cpuid-flag")
    prod = (out["dmi"]["product"] or "").lower()
    if any(w in prod for w in ("kvm", "qemu", "virtual", "vmware", "virtualbox")):
        out["is_vm"] = True
    return out


def psi_available(paths: Paths = Paths()) -> bool:
    return (Path(paths.proc) / "pressure/cpu").exists()


def resctrl_available(paths: Paths = Paths()) -> bool:
    return (Path(paths.sys) / "fs/resctrl/info").is_dir()


def thermal(paths: Paths = Paths()) -> dict[str, Any]:
    zones = list((Path(paths.sys) / "class/thermal").glob("thermal_zone*"))
    throttle = (Path(paths.sys) / "devices/system/cpu/cpu0/thermal_throttle").is_dir()
    return {"zones": len(zones), "throttle_counters": throttle}


def irq_count(paths: Paths = Paths()) -> int:
    return len([p for p in (Path(paths.proc) / "irq").glob("[0-9]*")]) if (Path(paths.proc) / "irq").is_dir() else 0


def os_release(paths: Paths = Paths()) -> str | None:
    for line in (read_text(Path(paths.etc) / "os-release") or "").splitlines():
        if line.startswith("PRETTY_NAME="):
            return line.split("=", 1)[1].strip().strip('"')
    return None


# -- tools, images, gpus -----------------------------------------------------

_VERSION_ARGS = {"podman": ["--version"], "docker": ["--version"], "runsc": ["--version"],
                 "crun": ["--version"], "runc": ["--version"], "perf": ["--version"],
                 "numactl": ["--version"], "nvidia-smi": ["--version"], "git": ["--version"],
                 "cc": ["--version"], "cpupower": ["--version"], "systemd-run": ["--version"]}


def tools(which: Callable[[str], str | None] = shutil.which, run: Runner = default_run) -> dict[str, Any]:
    out: dict[str, Any] = {}
    for name, args in _VERSION_ARGS.items():
        path = which(name)
        if not path:
            out[name] = None
            continue
        rc, text = run([path, *args], 10)
        first = next((ln.strip() for ln in text.splitlines() if ln.strip()), "")
        out[name] = {"path": path, "version": first[:120]}
    return out


def images(run: Runner = default_run, backend: str = "podman") -> list[dict[str, Any]]:
    rc, text = run([backend, "images", "--format", "json"], 30)
    if rc != 0 or not text.strip():
        return []
    try:
        rows = json.loads(text)
        if isinstance(rows, dict):
            rows = [rows]
    except ValueError:
        rows = []
        for line in text.splitlines():  # docker prints one JSON object per line
            try:
                rows.append(json.loads(line))
            except ValueError:
                pass
    out = []
    for r in rows:
        names = r.get("Names") or r.get("RepoTags") or []
        if not names and r.get("Repository"):
            names = [f"{r['Repository']}:{r.get('Tag', 'latest')}"]
        digests = r.get("Digest") or (r.get("RepoDigests") or [None])[0]
        out.append({"names": names, "id": (r.get("Id") or r.get("ID") or "")[:64],
                    "digest": digests if isinstance(digests, str) else None,
                    "size": r.get("Size")})
    return out


def gpus(run: Runner = default_run, which: Callable[[str], str | None] = shutil.which) -> list[dict[str, Any]]:
    if not which("nvidia-smi"):
        return []
    fields = ["index", "name", "uuid", "memory.total", "driver_version", "persistence_mode",
              "compute_mode", "clocks.max.sm", "clocks.max.memory", "pci.bus_id"]
    rc, text = run(["nvidia-smi", f"--query-gpu={','.join(fields)}", "--format=csv,noheader,nounits"], 20)
    if rc != 0:
        return []
    out = []
    for line in text.strip().splitlines():
        vals = [v.strip() for v in line.split(",")]
        if len(vals) != len(fields):
            continue
        row = dict(zip(fields, vals))
        out.append({"index": int(row["index"]), "name": row["name"], "uuid": row["uuid"],
                    "memory_mb": int(row["memory.total"]) if row["memory.total"].isdigit() else None,
                    "driver": row["driver_version"], "persistence_mode": row["persistence_mode"],
                    "compute_mode": row["compute_mode"], "max_sm_mhz": row["clocks.max.sm"],
                    "max_mem_mhz": row["clocks.max.memory"], "pci": row["pci.bus_id"]})
    return out


def cdi_available(paths: Paths = Paths()) -> bool:
    return any((Path(paths.etc) / "cdi").glob("*.yaml")) or any(Path("/var/run/cdi").glob("*.yaml"))


# -- the whole thing ---------------------------------------------------------

def host_inventory(paths: Paths = Paths(), run: Runner = default_run,
                   which: Callable[[str], str | None] = shutil.which,
                   lab_cpus: list[int] | None = None, image_backend: str | None = "podman") -> dict[str, Any]:
    topo = cpu_topology(paths)
    cpus = lab_cpus if lab_cpus else topo["online"]
    tl = tools(which, run)
    inv = {
        "collected_at": time.time(),
        "hostname": socket.gethostname(),
        "arch": {"arm64": "aarch64", "AMD64": "x86_64", "amd64": "x86_64"}.get(platform.machine(), platform.machine()),
        "kernel": platform.release(),
        "os": os_release(paths) or platform.platform(),
        "boot_id": read_text(Path(paths.proc) / "sys/kernel/random/boot_id"),
        "cpu": cpuinfo(paths),
        "topology": topo,
        "cmdline": kernel_cmdline(paths),
        "cpufreq": cpufreq(paths, cpus),
        "memory": memory(paths),
        "knobs": kernel_knobs(paths),
        "cgroups": cgroups(paths),
        "virtualization": virtualization(paths, run),
        "psi": psi_available(paths),
        "resctrl": resctrl_available(paths),
        "thermal": thermal(paths),
        "irqs": irq_count(paths),
        "tools": tl,
        "gpus": gpus(run, which),
        "cdi": cdi_available(paths),
        "images": images(run, image_backend) if image_backend and tl.get(image_backend) else [],
        "python": platform.python_version(),
    }
    return inv


def inventory_summary(inv: dict[str, Any]) -> str:
    t = inv["topology"]["summary"]
    cpu = inv["cpu"].get("model") or "unknown cpu"
    mem = (inv["memory"].get("total_kb") or 0) // 1024
    vm = inv["virtualization"].get("detect") or ("vm" if inv["virtualization"].get("is_vm") else "bare-metal")
    gpu = f", {len(inv['gpus'])} gpu" if inv.get("gpus") else ""
    return (f"{cpu} ({inv['arch']}): {t['sockets']}s/{t['cores']}c/{t['threads_per_core']}t, "
            f"{t['numa_nodes']} numa, {mem} MiB{gpu}, kernel {inv['kernel']}, {vm}")
