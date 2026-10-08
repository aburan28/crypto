"""``perf stat`` over exactly the job's CPUs, and its CSV output parsed.

Counting system-wide on the reserved CPUs (``-a -C``) is deliberate: once
the CPUs are exclusive, every instruction on them is the job's or is noise
that the fidelity checks should see, and it works the same for a container,
a gVisor sandbox and a bare process.
"""
from __future__ import annotations

import shutil
from typing import Any, Callable

Runner = Callable[[list[str], float], tuple[int, str]]

NAME_MAP = {
    "task-clock": "task_clock_ms", "context-switches": "context_switches", "cs": "context_switches",
    "cpu-migrations": "cpu_migrations", "migrations": "cpu_migrations",
    "page-faults": "page_faults", "faults": "page_faults", "cycles": "cycles",
    "instructions": "instructions", "branches": "branches", "branch-misses": "branch_misses",
    "cache-misses": "cache_misses", "cache-references": "cache_references",
    "LLC-load-misses": "llc_load_misses", "stalled-cycles-frontend": "stalled_cycles_frontend",
    "stalled-cycles-backend": "stalled_cycles_backend", "cpu-clock": "cpu_clock_ms",
}
HARDWARE = {"cycles", "instructions", "branches", "branch-misses", "cache-misses", "cache-references"}


def find_perf(which: Callable[[str], str | None] = shutil.which) -> str | None:
    return which("perf")


def probe(perf: str | None, run: Runner, cpus: list[int] | None = None) -> dict[str, Any]:
    """Can perf count at all here, and are hardware events exposed?"""
    if not perf:
        return {"available": False, "hw_events": False, "detail": "perf not installed"}
    argv = [perf, "stat", "-x,", "-e", "task-clock,cycles,instructions"]
    if cpus:
        argv += ["-a", "-C", ",".join(map(str, cpus))]
    argv += ["--", "sleep", "0.05"]
    rc, text = run(argv, 20)
    if rc != 0:
        return {"available": False, "hw_events": False, "detail": text.strip()[-300:]}
    parsed = parse_csv(text)
    hw = any(parsed["events"].get(e, {}).get("status") == "ok" for e in ("cycles", "instructions"))
    sw = parsed["events"].get("task-clock", {}).get("status") == "ok"
    return {"available": sw, "hw_events": hw, "detail": "ok" if sw else text.strip()[-300:]}


def argv_prefix(perf: str, events: list[str], cpus: list[int], outfile: str,
                system_wide: bool = True) -> list[str]:
    argv = [perf, "stat", "-x,", "-o", outfile, "-e", ",".join(events)]
    if system_wide:
        argv += ["-a", "-C", ",".join(map(str, cpus))]
    return argv + ["--"]


def parse_csv(text: str) -> dict[str, Any]:
    events: dict[str, dict[str, Any]] = {}
    elapsed = None
    for line in text.splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        parts = line.split(",")
        if len(parts) < 3:
            continue
        raw, unit, name = parts[0].strip(), parts[1].strip(), parts[2].strip()
        if not name:
            continue
        base = name.split(":")[0]
        status, value = "ok", None
        if raw.startswith("<not supported>"):
            status = "not supported"
        elif raw.startswith("<not counted>"):
            status = "not counted"
        else:
            try:
                value = float(raw)
            except ValueError:
                status = "unparsed"
        run_pct = None
        if len(parts) >= 5:
            try:
                run_pct = float(parts[4])
            except ValueError:
                pass
        events[base] = {"value": value, "unit": unit, "status": status, "run_pct": run_pct}
    return {"events": events, "elapsed_s": elapsed}


def counters(parsed: dict[str, Any]) -> dict[str, Any]:
    """Flatten parsed events into the result's counter block."""
    out: dict[str, Any] = {}
    status: dict[str, str] = {}
    for name, ev in parsed["events"].items():
        key = NAME_MAP.get(name, name.replace("-", "_"))
        v = ev["value"]
        if v is not None and key not in ("task_clock_ms", "cpu_clock_ms"):
            v = int(v)
        out[key] = v
        status[key] = ev["status"]
    ins, cyc = out.get("instructions"), out.get("cycles")
    out["ipc"] = (ins / cyc) if ins and cyc else None
    out["_status"] = status
    out["hardware_counters"] = any(parsed["events"].get(e, {}).get("status") == "ok" for e in HARDWARE)
    return out
