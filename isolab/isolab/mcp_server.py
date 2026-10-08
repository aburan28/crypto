"""MCP server: lets an agent run experiments on the lab and read them back.

The server runs on the agent's machine. It uploads local input files from
there, validates a job against the live roster before enqueueing it, and
writes fetched artifacts to local paths. It never runs anything itself.
Supports mcp 1.x (FastMCP) and 2.x (MCPServer).
"""
from __future__ import annotations

import asyncio
import json
import os
from pathlib import Path
from typing import Any

try:  # mcp >= 2
    from mcp.server.mcpserver import MCPServer as FastMCP
except ImportError:  # mcp 1.x
    from mcp.server.fastmcp import FastMCP

from . import protocol
from .blobs import resolve_local_inputs, sha256_file
from .fabric import Fabric
from .planner import match_worker
from .views import result_view, worker_line as _worker_line

mcp = FastMCP("isolab", instructions=(
    "isolab runs experiments on dedicated Linux lab hosts under enforced measurement fidelity "
    "(exclusive pinned CPUs, one NUMA node, interrupts moved away, PSI-gated start, perf counters) "
    "and returns results with the evidence of isolation attached. Start with isolab_overview to see "
    "which workers are online and what they can do. Submit with isolab_run (or isolab_submit for a full "
    "isolab.job/v1 document), then isolab_wait or isolab_status, then isolab_result. "
    "A result whose status is timeout, cancelled or infra_error (outcome_class not_completed) is never "
    "evidence about what was measured. fidelity.grade is the isolation tier the run actually got (A best), "
    "fidelity.violations lists every check that failed, and contended repeats are excluded from summary. "
    "Instructions and cycles are the primary metric; wall time is secondary. Call isolab_schema for the schemas."))

_fabric: Fabric | None = None
_lock = asyncio.Lock()


async def fabric() -> Fabric:
    global _fabric
    async with _lock:
        if _fabric is None or _fabric.nc is None or _fabric.nc.is_closed:
            _fabric = await Fabric(name="isolab-mcp").connect()
    return _fabric


async def _eligible(spec: dict[str, Any]) -> tuple[list[str], dict[str, list[str]]]:
    f = await fabric()
    workers = await f.workers()
    ok, why = [], {}
    for w in workers:
        reasons = match_worker(spec, w)
        if reasons:
            why[w["id"]] = reasons
        else:
            ok.append(w["id"])
    return ok, why


async def _submit(spec: dict[str, Any], base_dir: str | None = None) -> dict[str, Any]:
    f = await fabric()
    norm = protocol.normalize_spec(spec, allow_local_files=True)
    pending: list[tuple[str, Path]] = []
    resolve_local_inputs(norm, lambda d, p: pending.append((d, p)), lambda d: False, base_dir)
    uploads = []
    for digest, path in pending:
        if not await f.has_blob(digest):
            await f.put_blob(digest, path)
        uploads.append({"sha256": digest, "bytes": path.stat().st_size, "file": str(path)})
    norm = protocol.normalize_spec(norm)
    ok, why = await _eligible(norm)
    if not ok:
        detail = "; ".join(f"{w}: {', '.join(r)}" for w, r in why.items()) or "no workers online"
        raise ValueError(f"no online worker can take this job: {detail}")
    out = await f.submit(norm)
    out.update(eligible_workers=ok, ineligible=why, uploaded=uploads)
    return out


@mcp.tool()
async def isolab_overview() -> dict[str, Any]:
    """Hub reachability, online workers (one line each), queue depths and running jobs. Call this first."""
    f = await fabric()
    workers = await f.workers()
    running = await f.list_jobs(limit=20, state="running")
    queued = await f.list_jobs(limit=20, state="queued")
    return {"hub": f.urls, "namespace": f.ns, "workers": [_worker_line(w) for w in workers],
            "queue": await f.queue_depths(), "running": running, "queued": queued}


@mcp.tool()
async def isolab_workers(include_stale: bool = False) -> list[dict[str, Any]]:
    """The roster: each worker's shape, pools, labels, images, backends and the isolation tier it can reach."""
    f = await fabric()
    return [_worker_line(w) for w in await f.workers(include_stale)]


@mcp.tool()
async def isolab_worker(worker_id: str) -> dict[str, Any] | None:
    """One worker's full inventory: topology, NUMA nodes, cpufreq, kernel knobs, capabilities, images, GPUs."""
    f = await fabric()
    for w in await f.workers(include_stale=True):
        if w["id"] == worker_id or w.get("hostname") == worker_id:
            return w
    return None


@mcp.tool()
async def isolab_doctor(worker_id: str) -> dict[str, Any] | None:
    """A worker's readiness for strict runs: each check with its value and the command that fixes it."""
    from .doctor import as_json, diagnose
    w = await isolab_worker(worker_id)
    if w is None:
        return None
    inv = {**w, "topology": w["topology"], "tools": w.get("tools") or {}}
    linfo = w.get("launcher") or {}
    items = diagnose(inv, w["capabilities"], w["lab_cpus"], launcher=linfo.get("spin"),
                     launcher_present=bool(linfo.get("spin")) and linfo.get("error") is None)
    return {"worker": w["id"], "max_tier": w.get("max_tier"), "items": as_json(items)}


@mcp.tool()
async def isolab_run(command: list[str], image: str | None = None, cpus: int = 1, memory_mb: int | None = None,
                     repeats: int = 1, warmups: int = 0, policy: str = "standard", worker: str | None = None,
                     pool: str = "default", inputs: list[dict[str, Any]] | None = None,
                     build: list[list[str]] | None = None, env: dict[str, str] | None = None,
                     cwd: str = ".", timeout_s: float = 3600, backend: str = "auto", oci_runtime: str = "auto",
                     network: str = "none", numa_node: str | int = "single", smt: str = "isolate", gpus: int = 0,
                     require: list[str] | None = None, cooldown_s: float = 0, drop_caches: bool = False,
                     labels: dict[str, str] | None = None, name: str | None = None,
                     verify_certificate: bool = False, idempotency_key: str | None = None,
                     input_base_dir: str | None = None) -> dict[str, Any]:
    """Run a command on a lab worker under fidelity enforcement and return its job_id.

    `command` is argv (no shell), run in the working directory where `inputs` were placed. Each input is
    {"path": "where/in/workdir", ...} with exactly one of: "content" (inline text; add "encoding": "base64"
    for binary), "local_file" (a path on this machine, uploaded for you; executable bit preserved), or
    "git": {"url", "commit" (40 hex), "sparse_paths"}. `build` steps run untimed before the measured repeats.
    `image` names a container image the worker has (see isolab_workers); omit it for the worker's default,
    or set backend="direct" to run on the host without a container. `policy` is strict (tier A, bare metal,
    counters, refuses otherwise), standard, or best_effort. `repeats` timed runs after `warmups`; strict
    re-runs contended repeats. `require` adds fidelity requirements such as perf_counters, bare_metal,
    cpu_partition, numa_bind, irq_moved, evicted, frequency_pinned, thp_never, smt_off, no_steal.
    The command may write metrics.json and other outputs into $ISOLAB_OUTPUT_DIR; they come back as artifacts.
    Set verify_certificate when it writes certificate.json (a discrete_log claim) to have it independently checked.
    The job is validated against the live roster first and refused, with reasons per worker, if nobody can take it."""
    spec: dict[str, Any] = {
        "schema": protocol.SPEC_SCHEMA_ID, "name": name, "pool": pool,
        "runtime": {"backend": backend, "oci_runtime": oci_runtime, "image": image, "network": network},
        "inputs": inputs or [],
        "build": [{"argv": b} for b in (build or [])],
        "command": {"argv": command, "cwd": cwd, "env": env or {}},
        "measure": {"warmups": warmups, "repeats": repeats, "cooldown_s": cooldown_s, "drop_caches": drop_caches},
        "resources": {"cpus": cpus, "memory_mb": memory_mb, "numa_node": numa_node, "smt": smt, "gpus": gpus},
        "fidelity": {"policy": policy, **({"require": require} if require else {})},
        "placement": {"worker": worker} if worker else {},
        "limits": {"timeout_s": timeout_s},
        "labels": labels or {},
    }
    if name is None:
        spec.pop("name")
    if verify_certificate:
        spec["verify"] = {"builtin": "taskq-certificate"}
    if idempotency_key:
        spec["idempotency_key"] = idempotency_key
    return await _submit(spec, input_base_dir)


@mcp.tool()
async def isolab_submit(spec: dict[str, Any], input_base_dir: str | None = None) -> dict[str, Any]:
    """Enqueue a full isolab.job/v1 document (see isolab_schema). local_file inputs are uploaded first."""
    return await _submit(spec, input_base_dir)


@mcp.tool()
async def isolab_status(job_id: str) -> dict[str, Any] | None:
    """State, worker, attempt, phase and repeat, elapsed, log tail, and why a job is waiting if it is."""
    f = await fabric()
    rec = await f.get_job(job_id)
    if rec is None:
        return None
    out = {k: v for k, v in rec.items() if k not in ("spec", "_revision")}
    out["progress"] = await f.get_progress(job_id)
    if rec["state"] in protocol.TERMINAL_STATES:
        res = await f.get_result(job_id)
        if res:
            out["result_status"] = res["status"]
            out["grade"] = res["fidelity"]["grade"]
    return out


@mcp.tool()
async def isolab_wait(job_id: str, timeout_s: float = 60) -> dict[str, Any] | None:
    """Block up to timeout_s (max 300) for a terminal state; returns the status."""
    f = await fabric()
    await f.wait(job_id, min(max(timeout_s, 0), 300))
    return await isolab_status(job_id)


@mcp.tool()
async def isolab_result(job_id: str, section: str = "summary") -> dict[str, Any] | None:
    """The write-once result. section: summary (default, compact), fidelity, runs, host, placement, full."""
    f = await fabric()
    res = await f.get_result(job_id)
    return None if res is None else result_view(res, section)


@mcp.tool()
async def isolab_logs(job_id: str, tail_bytes: int = 8192) -> dict[str, Any]:
    """stdout and stderr tails: live from progress while running, from the stored logs when finished."""
    f = await fabric()
    res = await f.get_result(job_id)
    out: dict[str, Any] = {"job_id": job_id}
    if res:
        for name in ("stdout", "stderr"):
            try:
                out[name] = (await f.get_artifact_bytes(job_id, f"_isolab_logs/{name}.log", tail_bytes)).decode(errors="replace")
            except Exception:  # noqa: BLE001
                out[name] = (res.get("logs") or {}).get(f"{name}_tail")
    else:
        prog = await f.get_progress(job_id)
        out["stdout"] = (prog or {}).get("log_tail")
        out["progress"] = prog
    return out


@mcp.tool()
async def isolab_artifacts(job_id: str) -> list[dict[str, Any]] | None:
    """Every file the job wrote, with sha256 and size, and whether it was stored."""
    f = await fabric()
    res = await f.get_result(job_id)
    return None if res is None else res.get("artifacts")


@mcp.tool()
async def isolab_fetch(job_id: str, path: str | None = None, dest_dir: str = ".") -> dict[str, Any]:
    """Download one artifact (path) or all stored artifacts into dest_dir on this machine; hashes verified."""
    f = await fabric()
    res = await f.get_result(job_id)
    if res is None:
        raise ValueError(f"{job_id} has no result yet")
    wanted = [a for a in res.get("artifacts") or [] if a.get("stored") and (path is None or a["path"] == path)]
    if path and not wanted:
        raise ValueError(f"{path!r} is not a stored artifact of {job_id}")
    dest = Path(dest_dir).expanduser() / job_id
    got = []
    for a in wanted:
        target = await f.get_artifact(job_id, a["path"], dest / a["path"])
        ok = sha256_file(target) == a["sha256"]
        got.append({"path": str(target), "bytes": a["bytes"], "sha256_ok": ok})
    return {"dest": str(dest), "files": got}


@mcp.tool()
async def isolab_jobs(state: str | None = None, pool: str | None = None, worker: str | None = None,
                      labels: dict[str, str] | None = None, limit: int = 50) -> list[dict[str, Any]]:
    """Newest first, filtered by state (queued, claimed, running, succeeded, failed, timeout, cancelled, infra_error, dead), pool, worker, labels."""
    f = await fabric()
    return await f.list_jobs(min(limit, 500), state, pool, worker, labels)


@mcp.tool()
async def isolab_cancel(job_id: str) -> str:
    """Cancel a queued job outright or ask the running worker to stop it. Returns the resulting state."""
    f = await fabric()
    return await f.cancel(job_id)


@mcp.tool()
async def isolab_calibrate(worker: str | None = None, repeats: int = 10, cpus: int = 1, policy: str = "standard",
                           iterations: int = 2_000_000_000, backend: str = "auto", image: str | None = None) -> dict[str, Any]:
    """Measure a worker's noise floor: run the fixed-work isolab-spin kernel `repeats` times under isolation.
    The result's summary.wall_s.cv and summary.instructions are the A/A spread every A/B difference on that host is read against."""
    f = await fabric()
    workers = await f.workers()
    if worker:
        workers = [w for w in workers if w["id"] == worker or w.get("hostname") == worker]
    cands = [w for w in workers if w.get("calibration")]
    if not cands:
        raise ValueError("no online worker advertises a calibration kernel (the worker needs a C compiler at start)")
    w = cands[0]
    cal = w["calibration"]
    spec = {"schema": protocol.SPEC_SCHEMA_ID, "name": f"calibrate-{w['id']}", "pool": w["pools"][0],
            "runtime": {"backend": backend, "image": image},
            "inputs": [{"path": "isolab-spin", "sha256": cal["sha256"], "bytes": cal["bytes"], "mode": "0755"}],
            "command": {"argv": ["./isolab-spin", str(iterations)]},
            "measure": {"warmups": 1, "repeats": repeats},
            "resources": {"cpus": cpus},
            "fidelity": {"policy": policy},
            "placement": {"worker": w["id"]},
            "labels": {"isolab": "calibration"}}
    if image is None:
        spec["runtime"].pop("image")
    return await _submit(spec)


@mcp.tool()
async def isolab_schema() -> dict[str, Any]:
    """The isolab.job/v1 and isolab.result/v1 JSON schemas."""
    return {"job_spec": protocol.SPEC_SCHEMA, "result": protocol.RESULT_SCHEMA}


def main() -> None:
    mcp.run(os.environ.get("ISOLAB_MCP_TRANSPORT", "stdio"))


if __name__ == "__main__":
    main()
