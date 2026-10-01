"""Readers for the conditions around a run, and the thread that samples them.

Everything is a delta of a kernel counter: PSI stall totals, per-CPU tick
accounting on the job's CPUs (idle, irq, softirq, steal), interrupts that
landed on them, frequency, thermal throttling, other processes' CPU time,
network and disk bytes. The summary a run carries is computed from the
samples, and the samples themselves are kept as an artifact.
"""
from __future__ import annotations

import os
import re
import statistics
import threading
import time
from pathlib import Path
from typing import Any, Callable, Iterable

from .inventory import Paths, read_int, read_text

TICK = os.sysconf("SC_CLK_TCK") if hasattr(os, "sysconf") else 100
_STAT_FIELDS = ("user", "nice", "system", "idle", "iowait", "irq", "softirq", "steal", "guest", "guest_nice")


def read_psi(kind: str, paths: Paths = Paths()) -> dict[str, dict[str, float]] | None:
    text = read_text(Path(paths.proc) / "pressure" / kind)
    if text is None:
        return None
    out: dict[str, dict[str, float]] = {}
    for line in text.splitlines():
        name, *fields = line.split()
        out[name] = {k: float(v) for k, v in (f.split("=") for f in fields)}
    return out


def read_cgroup_psi(cgroup: Path | None, kind: str) -> dict[str, dict[str, float]] | None:
    if cgroup is None:
        return None
    text = read_text(cgroup / f"{kind}.pressure")
    if text is None:
        return None
    out: dict[str, dict[str, float]] = {}
    for line in text.splitlines():
        name, *fields = line.split()
        out[name] = {k: float(v) for k, v in (f.split("=") for f in fields)}
    return out


def read_cpu_stat(paths: Paths = Paths()) -> dict[int, dict[str, int]]:
    out: dict[int, dict[str, int]] = {}
    for line in (read_text(Path(paths.proc) / "stat") or "").splitlines():
        if line.startswith("cpu") and line[3:4].isdigit():
            parts = line.split()
            cpu = int(parts[0][3:])
            vals = [int(x) for x in parts[1:11]]
            vals += [0] * (10 - len(vals))
            out[cpu] = dict(zip(_STAT_FIELDS, vals))
    return out


def _per_cpu_table(text: str | None) -> dict[int, int]:
    """Sum the per-CPU columns of /proc/interrupts or /proc/softirqs."""
    if not text:
        return {}
    lines = text.splitlines()
    header = lines[0].split()
    cpus = [int(h[3:]) for h in header if h.startswith("CPU") and h[3:].isdigit()]
    totals = {c: 0 for c in cpus}
    for line in lines[1:]:
        parts = line.split()
        if not parts:
            continue
        cols = parts[1:1 + len(cpus)]
        for c, v in zip(cpus, cols):
            if v.isdigit():
                totals[c] += int(v)
    return totals


def read_interrupts(paths: Paths = Paths()) -> dict[int, int]:
    return _per_cpu_table(read_text(Path(paths.proc) / "interrupts"))


def read_softirqs(paths: Paths = Paths()) -> dict[int, int]:
    return _per_cpu_table(read_text(Path(paths.proc) / "softirqs"))


def read_freq(cpus: Iterable[int], paths: Paths = Paths()) -> dict[int, int | None]:
    base = Path(paths.sys) / "devices/system/cpu"
    out = {}
    for c in cpus:
        d = base / f"cpu{c}" / "cpufreq"
        out[c] = read_int(d / "scaling_cur_freq") or read_int(d / "cpuinfo_cur_freq")
    return out


def read_throttle(cpus: Iterable[int], paths: Paths = Paths()) -> dict[int, int]:
    base = Path(paths.sys) / "devices/system/cpu"
    out = {}
    for c in cpus:
        d = base / f"cpu{c}" / "thermal_throttle"
        core = read_int(d / "core_throttle_count")
        pkg = read_int(d / "package_throttle_count")
        if core is not None or pkg is not None:
            out[c] = (core or 0) + (pkg or 0)
    return out


def read_thermal(paths: Paths = Paths()) -> dict[str, float]:
    out = {}
    for zone in sorted((Path(paths.sys) / "class/thermal").glob("thermal_zone*")):
        temp = read_int(zone / "temp")
        if temp is not None:
            out[(read_text(zone / "type") or zone.name)[:32] + f"#{zone.name[12:]}"] = temp / 1000.0
    return out


def read_meminfo(paths: Paths = Paths()) -> dict[str, int]:
    out = {}
    for line in (read_text(Path(paths.proc) / "meminfo") or "").splitlines():
        k, _, v = line.partition(":")
        try:
            out[k.strip()] = int(v.split()[0])
        except (ValueError, IndexError):
            pass
    return out


def read_netdev(paths: Paths = Paths()) -> dict[str, int]:
    rx = tx = 0
    for line in (read_text(Path(paths.proc) / "net/dev") or "").splitlines()[2:]:
        name, _, rest = line.partition(":")
        if name.strip() == "lo":
            continue
        f = rest.split()
        if len(f) >= 9:
            rx += int(f[0])
            tx += int(f[8])
    return {"rx_bytes": rx, "tx_bytes": tx}


def read_diskstats(paths: Paths = Paths()) -> dict[str, int]:
    """Sectors moved by whole devices (partitions would double count)."""
    rows = [ln.split() for ln in (read_text(Path(paths.proc) / "diskstats") or "").splitlines()]
    rows = [f for f in rows if len(f) >= 10]
    names = [f[2] for f in rows]

    def is_partition(name: str) -> bool:
        return any(other != name and name.startswith(other) and re.fullmatch(r"p?\d+", name[len(other):])
                   for other in names)
    rd = wr = 0
    for f in rows:
        if not is_partition(f[2]):
            rd += int(f[5])
            wr += int(f[9])
    return {"read_sectors": rd, "write_sectors": wr}


def read_loadavg(paths: Paths = Paths()) -> list[float]:
    text = read_text(Path(paths.proc) / "loadavg")
    return [float(x) for x in text.split()[:3]] if text else [0.0, 0.0, 0.0]


def cpu_snapshot(paths: Paths = Paths()) -> dict[int, tuple[str, int, int]]:
    """Every process: (command, utime+stime ticks, parent pid)."""
    out = {}
    root = Path(paths.proc)
    if not root.is_dir():
        return out
    for entry in root.iterdir():
        if not entry.name.isdigit():
            continue
        try:
            stat = (entry / "stat").read_text()
        except OSError:
            continue
        try:
            name = stat[stat.index("(") + 1:stat.rindex(")")]
            fields = stat[stat.rindex(")") + 2:].split()
            out[int(entry.name)] = (name, int(fields[11]) + int(fields[12]), int(fields[1]))
        except (ValueError, IndexError):
            continue
    return out


def descendants(root: int, snapshot: dict[int, tuple[str, int, int]]) -> set[int]:
    children: dict[int, list[int]] = {}
    for pid, (_, _, ppid) in snapshot.items():
        children.setdefault(ppid, []).append(pid)
    out, stack = set(), [root]
    while stack:
        pid = stack.pop()
        for child in children.get(pid, []):
            if child not in out:
                out.add(child)
                stack.append(child)
    return out


def other_use(before: dict[int, tuple[str, int, int]], after: dict[int, tuple[str, int, int]],
              exclude: set[int]) -> dict[str, float]:
    used = {}
    for pid, (name, ticks, _) in after.items():
        if pid in exclude:
            continue
        delta = ticks - before.get(pid, (name, 0, 0))[1]
        if delta > 0:
            used[f"{name}[{pid}]"] = delta / TICK
    return dict(sorted(used.items(), key=lambda kv: -kv[1]))


def cgroup_procs(cgroup: Path | None) -> set[int]:
    if cgroup is None:
        return set()
    try:
        return {int(x) for x in (cgroup / "cgroup.procs").read_text().split()}
    except (OSError, ValueError):
        return set()


def _stats(xs: list[float]) -> dict[str, float] | None:
    xs = [x for x in xs if x is not None]
    if not xs:
        return None
    mean = statistics.fmean(xs)
    sd = statistics.pstdev(xs) if len(xs) > 1 else 0.0
    return {"min": min(xs), "max": max(xs), "mean": mean, "cv": (sd / mean) if mean else 0.0, "n": len(xs)}


class Sampler(threading.Thread):
    """Sample the conditions around a run every ``period`` seconds.

    ``exclude`` returns the pids whose CPU time is the job's own (so that it
    is not charged as contention). ``cgroup`` is the job's cgroup directory,
    read for its own CPU usage and pressure when present.
    """

    def __init__(self, job_cpus: list[int], period: float = 1.0, paths: Paths = Paths(),
                 exclude: Callable[[], set[int]] | None = None, cgroup: Path | None = None,
                 other_cpu_threshold: float = 0.05):
        super().__init__(daemon=True, name="isolab-sampler")
        self.cpus = sorted(job_cpus)
        self.period = period
        self.paths = paths
        self.exclude = exclude or (lambda: set())
        self.cgroup = cgroup
        self.threshold = other_cpu_threshold
        self.samples: list[dict[str, Any]] = []
        self._halt = threading.Event()
        self.started_at: float | None = None
        self.stopped_at: float | None = None
        self._first: dict[str, Any] = {}
        self._last: dict[str, Any] = {}
        self._my_pid = os.getpid()

    def _snapshot(self) -> dict[str, Any]:
        stat = read_cpu_stat(self.paths)
        job = {k: sum(stat.get(c, {}).get(k, 0) for c in self.cpus) for k in _STAT_FIELDS}
        irqs = read_interrupts(self.paths)
        sirqs = read_softirqs(self.paths)
        return {
            "t": time.time(), "mono": time.monotonic(),
            "psi": {k: read_psi(k, self.paths) for k in ("cpu", "memory", "io")},
            "job_ticks": job,
            "irqs_job": sum(irqs.get(c, 0) for c in self.cpus) if irqs else None,
            "softirqs_job": sum(sirqs.get(c, 0) for c in self.cpus) if sirqs else None,
            "freq_khz": read_freq(self.cpus, self.paths),
            "throttle": read_throttle(self.cpus, self.paths),
            "thermal_c": read_thermal(self.paths),
            "loadavg": read_loadavg(self.paths),
            "mem_available_kb": read_meminfo(self.paths).get("MemAvailable"),
            "net": read_netdev(self.paths), "disk": read_diskstats(self.paths),
            "procs": cpu_snapshot(self.paths),
            "cg_cpu_usage_us": _cg_int(self.cgroup, "cpu.stat", "usage_usec"),
            "cg_mem_current": read_int(self.cgroup / "memory.current") if self.cgroup else None,
            "cg_psi": {k: read_cgroup_psi(self.cgroup, k) for k in ("cpu", "memory", "io")},
        }

    def run(self) -> None:
        self.started_at = time.time()
        self._first = prev = self._snapshot()
        while not self._halt.wait(self.period):
            cur = self._snapshot()
            self.samples.append(self._delta(prev, cur))
            prev = cur
        cur = self._snapshot()
        if cur["mono"] - prev["mono"] > 0.05:
            self.samples.append(self._delta(prev, cur))
        self._last = cur
        self.stopped_at = time.time()

    def stop(self) -> None:
        self._halt.set()
        self.join(timeout=max(5.0, self.period * 3))

    def _delta(self, a: dict[str, Any], b: dict[str, Any]) -> dict[str, Any]:
        dt = max(b["mono"] - a["mono"], 1e-6)
        ticks = {k: b["job_ticks"][k] - a["job_ticks"][k] for k in _STAT_FIELDS}
        total_ticks = sum(ticks.values())
        exclude = self.exclude() | {self._my_pid}
        others = other_use(a["procs"], b["procs"], exclude)
        other_s = sum(others.values())
        if total_ticks > 0:
            job_cpu = {"busy_pct": 100.0 * (1 - ticks["idle"] / total_ticks),
                       "irq_pct": 100.0 * (ticks["irq"] + ticks["softirq"]) / total_ticks,
                       "steal_pct": 100.0 * ticks["steal"] / total_ticks,
                       "user_pct": 100.0 * (ticks["user"] + ticks["nice"]) / total_ticks,
                       "system_pct": 100.0 * ticks["system"] / total_ticks, "ticks": ticks}
        else:  # no per-cpu accounting here (not Linux, or a window shorter than a tick)
            job_cpu = {"busy_pct": None, "irq_pct": None, "steal_pct": None, "user_pct": None, "system_pct": None, "ticks": None}
        return {
            "t": b["t"], "dt": dt,
            "psi_delta_us": _psi_delta(a["psi"], b["psi"]),
            "psi_avg10": {k: (b["psi"][k] or {}).get("some", {}).get("avg10") for k in ("cpu", "memory", "io")},
            "job_cpu": job_cpu,
            "irqs_job": None if a["irqs_job"] is None else b["irqs_job"] - a["irqs_job"],
            "softirqs_job": None if a["softirqs_job"] is None else b["softirqs_job"] - a["softirqs_job"],
            "freq_khz": b["freq_khz"],
            "throttle_delta": sum(b["throttle"].values()) - sum(a["throttle"].values()),
            "thermal_max_c": max(b["thermal_c"].values()) if b["thermal_c"] else None,
            "loadavg": b["loadavg"], "mem_available_kb": b["mem_available_kb"],
            "net_bytes": (b["net"]["rx_bytes"] - a["net"]["rx_bytes"]) + (b["net"]["tx_bytes"] - a["net"]["tx_bytes"]),
            "disk_sectors": (b["disk"]["read_sectors"] - a["disk"]["read_sectors"]) + (b["disk"]["write_sectors"] - a["disk"]["write_sectors"]),
            "other_cpu_s": other_s, "other_top": dict(list(others.items())[:5]),
            "contended": other_s > self.threshold * dt,
            "cg_cpu_us": None if a["cg_cpu_usage_us"] is None or b["cg_cpu_usage_us"] is None else b["cg_cpu_usage_us"] - a["cg_cpu_usage_us"],
            "cg_mem_current": b["cg_mem_current"],
            "cg_psi_delta_us": _psi_delta(a["cg_psi"], b["cg_psi"]),
        }

    def summary(self) -> dict[str, Any]:
        s = self.samples
        if not s:
            return {"n_samples": 0}
        wall = sum(x["dt"] for x in s)
        freq_by_cpu: dict[int, list[float]] = {c: [] for c in self.cpus}
        for x in s:
            for c, khz in x["freq_khz"].items():
                if khz:
                    freq_by_cpu[c].append(khz / 1000.0)
        freq = {str(c): _stats(v) for c, v in freq_by_cpu.items() if v}
        other = sum(x["other_cpu_s"] for x in s)
        tops: dict[str, float] = {}
        for x in s:
            for k, v in x["other_top"].items():
                tops[k] = tops.get(k, 0.0) + v
        return {
            "n_samples": len(s), "wall_s": wall,
            "psi_delta_us": _sum_psi(x["psi_delta_us"] for x in s),
            "psi_some_avg10_max": {k: max((x["psi_avg10"][k] or 0.0) for x in s) for k in ("cpu", "memory", "io")},
            "job_cpu_busy_pct": _wmean(s, lambda x: x["job_cpu"]["busy_pct"]),
            "job_cpu_irq_pct": _wmean(s, lambda x: x["job_cpu"]["irq_pct"]),
            "job_cpu_steal_pct": _wmean(s, lambda x: x["job_cpu"]["steal_pct"]),
            "irqs_on_job_cpus": _sum_opt(x["irqs_job"] for x in s),
            "softirqs_on_job_cpus": _sum_opt(x["softirqs_job"] for x in s),
            "other_cpu_s": other, "other_cpu_ratio": other / wall if wall else 0.0,
            "other_top": dict(sorted(tops.items(), key=lambda kv: -kv[1])[:8]),
            "contended_samples": sum(1 for x in s if x["contended"]),
            "freq_mhz": freq,
            "freq_cv_max": max((v["cv"] for v in freq.values()), default=None) if freq else None,
            "throttle_events": sum(x["throttle_delta"] for x in s),
            "thermal_max_c": max((x["thermal_max_c"] for x in s if x["thermal_max_c"] is not None), default=None),
            "loadavg_max": max(x["loadavg"][0] for x in s),
            "mem_available_min_kb": min((x["mem_available_kb"] for x in s if x["mem_available_kb"] is not None), default=None),
            "net_bytes": sum(x["net_bytes"] for x in s), "disk_sectors": sum(x["disk_sectors"] for x in s),
            "cg_cpu_s": (_sum_opt(x["cg_cpu_us"] for x in s) or 0) / 1e6 if any(x["cg_cpu_us"] is not None for x in s) else None,
            "cg_mem_peak_sampled": max((x["cg_mem_current"] for x in s if x["cg_mem_current"] is not None), default=None),
            "cg_psi_delta_us": _sum_psi(x["cg_psi_delta_us"] for x in s),
        }


def _cg_int(cgroup: Path | None, file: str, key: str) -> int | None:
    if cgroup is None:
        return None
    for line in (read_text(cgroup / file) or "").splitlines():
        k, _, v = line.partition(" ")
        if k == key:
            try:
                return int(v)
            except ValueError:
                return None
    return None


def _psi_delta(a: dict[str, Any], b: dict[str, Any]) -> dict[str, dict[str, float] | None]:
    out: dict[str, dict[str, float] | None] = {}
    for kind in ("cpu", "memory", "io"):
        pa, pb = a.get(kind), b.get(kind)
        if not pa or not pb:
            out[kind] = None
            continue
        out[kind] = {lvl: pb[lvl]["total"] - pa[lvl]["total"] for lvl in pb if lvl in pa}
    return out


def _sum_psi(deltas) -> dict[str, dict[str, float] | None]:
    acc: dict[str, dict[str, float] | None] = {"cpu": None, "memory": None, "io": None}
    for d in deltas:
        for kind, v in d.items():
            if v is None:
                continue
            cur = acc[kind] or {}
            for lvl, us in v.items():
                cur[lvl] = cur.get(lvl, 0.0) + us
            acc[kind] = cur
    return acc


def _wmean(samples, f) -> float | None:
    """Duration-weighted mean of f over the samples where it is known."""
    known = [x for x in samples if f(x) is not None]
    w = sum(x["dt"] for x in known)
    return sum(f(x) * x["dt"] for x in known) / w if w else None


def _sum_opt(values) -> int | None:
    vals = [v for v in values if v is not None]
    return sum(vals) if vals else None
