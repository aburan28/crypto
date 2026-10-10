"""The worker daemon: inventory, heartbeat, claim, execute, complete.

One job at a time per worker. A big host can run several *slots*
(``isolab worker --slots N`` or ``--slot CPULIST`` repeated): each slot is a
worker of its own with a disjoint set of whole cores, its own cgroup
partition and its own lock, and they share one host lock. Jobs in different
slots run at the same time, each one's CPU is excluded from the other's
contention checks and recorded as a co-tenant instead, and a co-tenanted
result is graded B at best. A ``strict`` job takes the host lock exclusively,
so it waits for the other slots to drain and none start until it is done.
SIGTERM stops the running child and hands the job back without counting an
attempt, so a restart never records a result.
"""
from __future__ import annotations

import asyncio
import logging
import os
import platform
import shutil
import signal
import socket
import subprocess
import sys
import time
from importlib import resources
from pathlib import Path
from typing import Any

from . import __version__, perfstat, protocol
from .blobs import LocalBlobCache, sha256_file
from .execute import Fenced, InfraError, WorkerContext, execute
from .fabric import Fabric, FencedError
from .inventory import Paths, default_run, host_inventory, inventory_summary
from .isolation import DEFAULT_LOCK, max_tier, probe_capabilities
from .planner import capacity, match_worker
from .runners import make_backends, oci_runtimes

log = logging.getLogger("isolab.worker")


class SlotRegistry:
    """The jobs running in this process's slots, so each can tell the others' CPU from contention."""

    def __init__(self) -> None:
        import threading
        self._lock = threading.Lock()
        self._jobs: dict[str, Any] = {}

    def enter(self, job_id: str, pids_fn) -> None:
        with self._lock:
            self._jobs[job_id] = pids_fn

    def leave(self, job_id: str) -> None:
        with self._lock:
            self._jobs.pop(job_id, None)

    def others(self, job_id: str) -> dict[str, set[int]]:
        with self._lock:
            items = [(j, f) for j, f in self._jobs.items() if j != job_id]
        out = {}
        for j, f in items:
            try:
                out[j] = set(f())
            except Exception:  # noqa: BLE001
                out[j] = set()
        return out


def build_launcher(dest_dir: Path, static: bool = True) -> tuple[Path | None, dict[str, Any]]:
    """Compile the in-container launcher and the calibration kernel into dest_dir."""
    dest_dir.mkdir(parents=True, exist_ok=True)
    cc = shutil.which("cc") or shutil.which("gcc") or shutil.which("clang")
    info: dict[str, Any] = {"cc": cc, "static": None, "spin": None, "error": None}
    if not cc:
        info["error"] = "no C compiler (install gcc or build-essential)"
        return None, info
    launcher = dest_dir / "isolab-launch"
    spin = dest_dir / "isolab-spin"
    src_l = resources.files("isolab").joinpath("launcher.c").read_text()
    src_s = resources.files("isolab").joinpath("spin.c").read_text()
    for target, src, name in ((launcher, src_l, "launcher"), (spin, src_s, "spin")):
        tmp = dest_dir / f"{name}.c"
        tmp.write_text(src)
        built = False
        for flags in (["-static"] if static else [], []):
            p = subprocess.run([cc, "-O2", *flags, "-o", str(target), str(tmp)], capture_output=True, text=True)
            if p.returncode == 0:
                built = True
                info["error"] = None
                if name == "launcher":
                    info["static"] = bool(flags)
                break
            info["error"] = p.stderr.strip()[-300:]
        tmp.unlink(missing_ok=True)
        if not built:
            return (None, info) if name == "launcher" else (launcher, info)
    info["spin"] = str(spin)
    return launcher, info


class Worker:
    def __init__(self, fabric: Fabric, *, worker_id: str, pools: list[str], labels: dict[str, str],
                 lab_cpus: list[int], state_dir: Path, backend: str = "auto", default_image: str | None = None,
                 oci_runtime: str | None = None, pull_policy: str = "missing", userns: str | None = None,
                 keep_job_dirs: bool = False, lock_path: str = DEFAULT_LOCK, cgroup_root: str = "/sys/fs/cgroup",
                 paths: Paths = Paths(), block_s: float = 5.0, heartbeat_s: float = 10.0,
                 slot: str | None = None, host_lab_cpus: list[int] | None = None,
                 registry: SlotRegistry | None = None):
        self.fabric = fabric
        self.id = worker_id
        self.pools = pools or ["default"]
        self.labels = labels
        self.lab_cpus = sorted(lab_cpus)
        self.state_dir = state_dir
        self.backend_pref = backend
        self.default_image = default_image
        self.oci_pref = oci_runtime
        self.pull_policy, self.userns, self.keep_job_dirs = pull_policy, userns, keep_job_dirs
        self.lock_path, self.cgroup_root, self.paths = lock_path, cgroup_root, paths
        self.block_s, self.heartbeat_s = block_s, heartbeat_s
        # slots: every slot's housekeeping is the host minus all slots' cpus
        self.slot, self.host_lab_cpus, self.registry = slot, sorted(host_lab_cpus or []), registry
        self._shutdown = asyncio.Event()
        self.current: str | None = None
        self.ctx: WorkerContext | None = None
        self.inventory: dict[str, Any] = {}
        self.started_at = time.time()
        self.loop: asyncio.AbstractEventLoop | None = None
        self.calibration: dict[str, Any] | None = None

    # ---- setup
    async def start(self) -> None:
        self.loop = asyncio.get_running_loop()
        state = self.state_dir
        for d in ("bin", "blobs", "git", "jobs"):
            (state / d).mkdir(parents=True, exist_ok=True)
        inv = host_inventory(self.paths, lab_cpus=self.lab_cpus or None, image_backend=self._image_backend())
        online = inv["topology"]["online"]
        if not self.lab_cpus:
            self.lab_cpus = online[1:] if len(online) > 1 else online
            log.warning("no --cpus given; lab cpus default to %s", self.lab_cpus)
        bad = sorted(set(self.lab_cpus) - set(online))
        if bad:
            raise SystemExit(f"--cpus names offline or absent cpus {bad}; online: {online}")
        housekeeping = sorted(set(online) - set(self.lab_cpus) - set(self.host_lab_cpus)) or online
        if set(housekeeping) == set(online):
            log.warning("lab cpus cover every cpu; the worker and the OS will share them (tier C at best)")
        if hasattr(os, "sched_setaffinity") and set(housekeeping) != set(online):
            try:
                os.sched_setaffinity(0, housekeeping)
            except OSError as err:
                log.warning("cannot pin the worker to housekeeping cpus %s: %s", housekeeping, err)
        launcher, linfo = build_launcher(state / "bin")  # slots start one after another, so this never races
        if launcher is None:
            log.warning("launcher unavailable: %s; container backends are off, direct uses host wait4", linfo.get("error"))
        perf = perfstat.find_perf()
        perf_probe = perfstat.probe(perf, default_run, self.lab_cpus[:1])
        privileged = hasattr(os, "geteuid") and os.geteuid() == 0
        backends = make_backends(privileged=privileged)
        if launcher is None:
            backends = {k: v for k, v in backends.items() if k == "direct"}
        if self.backend_pref != "auto":
            if self.backend_pref not in backends:
                raise SystemExit(f"backend {self.backend_pref!r} is not available; have {sorted(backends)}")
            default_backend = self.backend_pref
        else:
            default_backend = next((b for b in ("podman", "docker", "direct") if b in backends), "direct")
        runtimes = oci_runtimes()
        oci_default = self.oci_pref or next((r for r in ("crun", "runc") if r in runtimes), None)
        caps = probe_capabilities(self.paths, self.cgroup_root, perf_probe, numactl=bool(shutil.which("numactl")),
                                  container_backend=default_backend != "direct")
        self.caps, self.perf_probe = caps, perf_probe
        self.backends, self.default_backend, self.oci_default, self.runtimes = backends, default_backend, oci_default, runtimes
        self.housekeeping, self.launcher, self.launcher_info = housekeeping, launcher, linfo
        if linfo.get("spin"):
            spin = Path(linfo["spin"])
            digest = sha256_file(spin)
            try:
                if not await self.fabric.has_blob(digest):
                    await self.fabric.put_blob(digest, spin)
                self.calibration = {"sha256": digest, "bytes": spin.stat().st_size, "arch": inv["arch"], "path": "isolab-spin"}
            except Exception as err:  # noqa: BLE001
                log.warning("could not publish the calibration kernel: %s", err)
        self.inventory = inv
        self.ctx = WorkerContext(
            worker_id=self.id, hostname=socket.gethostname(), labels=self.labels, version=__version__,
            inventory=inv, caps=caps, lab_cpus=self.lab_cpus, housekeeping=housekeeping, backends=backends,
            default_backend=default_backend, default_image=self.default_image, oci_default=oci_default,
            launcher=str(launcher) if launcher else None, perf=perf, perf_probe=perf_probe,
            blob_cache=LocalBlobCache(state / "blobs"), fetch_blob=self._fetch_blob_sync,
            put_artifact=self._put_artifact_sync, git_cache_dir=state / "git",
            jobs_dir=state / "jobs" / self.slot if self.slot else state / "jobs",
            paths=self.paths, cgroup_root=self.cgroup_root, lock_path=self.lock_path, pull_policy=self.pull_policy,
            userns=self.userns, keep_job_dirs=self.keep_job_dirs, progress=self._progress_sync,
            slot=self.slot, slot_lock_path=f"{self.lock_path}.{self.slot}" if self.slot else None,
            lab_cgroup=f"isolab.lab.{self.slot}" if self.slot else "isolab.lab", registry=self.registry)
        log.info("worker %s: %s", self.id, inventory_summary(inv))
        log.info("lab cpus %s, housekeeping %s, backend %s (%s), tier up to %s, perf hw=%s",
                 self.lab_cpus, housekeeping, default_backend, oci_default, max_tier(caps), perf_probe.get("hw_events"))
        await self.fabric.register_worker(self.record())

    def _image_backend(self) -> str | None:
        if self.backend_pref in ("podman", "docker"):
            return self.backend_pref
        return "podman" if shutil.which("podman") else ("docker" if shutil.which("docker") else None)

    def record(self) -> dict[str, Any]:
        inv = self.inventory
        topo = inv["topology"]
        return {
            "id": self.id, "hostname": inv["hostname"], "version": __version__, "pools": self.pools,
            "labels": self.labels, "arch": inv["arch"], "kernel": inv["kernel"], "os": inv["os"],
            "summary": inventory_summary(inv), "cpu": inv["cpu"], "topology": topo,
            "lab_cpus": self.lab_cpus, "housekeeping_cpus": self.housekeeping,
            "capacity": capacity(topo, self.lab_cpus), "memory": inv["memory"], "gpus": inv["gpus"],
            "images": inv["images"], "backends": sorted(self.backends), "default_backend": self.default_backend,
            "default_image": self.default_image, "oci_runtimes": self.runtimes, "oci_default": self.oci_default,
            "capabilities": self.caps, "max_tier": max_tier(self.caps), "perf": self.perf_probe,
            "cpufreq": inv["cpufreq"], "knobs": inv["knobs"], "virtualization": inv["virtualization"],
            "cmdline": inv["cmdline"], "cgroups": inv["cgroups"], "psi": inv["psi"], "tools": inv["tools"],
            "launcher": self.launcher_info, "calibration": self.calibration,
            "slot": self.slot, "busy": self.current, "started_at": self.started_at, "pid": os.getpid(),
        }

    # ---- sync bridges for the execution thread
    def _fetch_blob_sync(self, digest: str, dest: Path) -> None:
        fut = asyncio.run_coroutine_threadsafe(self.fabric.get_blob(digest, dest), self.loop)
        fut.result(timeout=3600)

    def _put_artifact_sync(self, job_id: str, rel: str, path: Path) -> str:
        fut = asyncio.run_coroutine_threadsafe(self.fabric.put_artifact(job_id, rel, path), self.loop)
        return fut.result(timeout=3600)

    def _progress_sync(self, progress: dict[str, Any]) -> None:
        if self.current and self.loop:
            self._live_progress = {**getattr(self, "_live_progress", {}), **progress}
            asyncio.run_coroutine_threadsafe(self.fabric.put_progress(self.current, self._live_progress), self.loop)

    # ---- main loop
    def request_shutdown(self, *_: Any) -> None:
        log.info("shutdown requested")
        self._shutdown.set()

    async def run_forever(self, max_jobs: int | None = None) -> int:
        done = 0
        hb = asyncio.create_task(self._heartbeat_roster())
        try:
            while not self._shutdown.is_set():
                if await self.run_once():
                    done += 1
                    if max_jobs is not None and done >= max_jobs:
                        break
        finally:
            hb.cancel()
            await self.fabric.unregister_worker(self.id)
        return done

    async def _heartbeat_roster(self) -> None:
        while True:
            try:
                await self.fabric.register_worker(self.record())
            except Exception as err:  # noqa: BLE001
                log.warning("roster heartbeat failed: %s", err)
            await asyncio.sleep(self.heartbeat_s)

    async def run_once(self) -> bool:
        per_pool = max(0.5, self.block_s / len(self.pools))
        for pool in self.pools:
            if self._shutdown.is_set():
                return False
            try:
                msg = await self.fabric.fetch(pool, per_pool)
            except Exception as err:  # noqa: BLE001
                log.warning("fetch from %s failed: %s", pool, err)
                await asyncio.sleep(1.0)
                continue
            if msg is None:
                continue
            return await self._handle(msg)
        return False

    async def _handle(self, msg) -> bool:
        job_id = msg.data.decode(errors="replace")
        rec = await self.fabric.get_job(job_id)
        if rec is None:
            log.warning("%s: message without a record; retrying later", job_id)
            await msg.nak(delay=5)
            return False
        if rec["state"] in protocol.TERMINAL_STATES:
            await self.fabric.drop(msg)
            return False
        try:
            spec = protocol.normalize_spec(rec["spec"])
        except protocol.SpecError as err:
            log.warning("%s: declining, spec unsupported here: %s", job_id, err)
            await self.fabric.decline(job_id, self.id, [f"spec unsupported by worker {__version__}: {err}"], msg)
            return False
        reasons = match_worker(spec, self.record())
        if reasons:
            log.info("%s: declining: %s", job_id, "; ".join(reasons))
            await self.fabric.decline(job_id, self.id, reasons, msg)
            return False
        delivered = msg.metadata.num_delivered if msg.metadata else 1
        claimed = await self.fabric.claim(job_id, self.id, delivered)
        if claimed is None:
            await self.fabric.drop(msg)
            return False
        rec, fence, attempt = claimed
        self.current = job_id
        self._live_progress = {}
        try:
            await self._run(rec, spec, msg, attempt)
        finally:
            self.current = None
        return True

    async def _run(self, rec: dict[str, Any], spec: dict[str, Any], msg, attempt: int) -> None:
        job_id = rec["job_id"]
        state = {"stop": None}
        stop_event = asyncio.Event()

        async def beat() -> None:
            period = max(1.0, self.fabric.lease_s / 4)
            n = 0
            while not stop_event.is_set():
                try:
                    await asyncio.wait_for(stop_event.wait(), timeout=period if n else 2.0)
                    break
                except asyncio.TimeoutError:
                    pass
                n += 1
                try:
                    r = await self.fabric.heartbeat(rec, msg)
                except Exception as err:  # noqa: BLE001
                    log.warning("%s: heartbeat failed: %s", job_id, err)
                    continue
                if r != "ok":
                    state["stop"] = r
                live = self.ctx.live.get("stdout_log") if self.ctx else None
                tail = ""
                if live and Path(live).exists():
                    try:
                        with open(live, "rb") as fh:
                            fh.seek(max(0, Path(live).stat().st_size - 2048))
                            tail = fh.read().decode(errors="replace")
                    except OSError:
                        pass
                self._live_progress = {**getattr(self, "_live_progress", {}), "log_tail": tail, "attempt": attempt,
                                       "elapsed_s": round(time.time() - started, 1), "worker": self.id}
                await self.fabric.put_progress(job_id, self._live_progress)

        def stop() -> str | None:
            if self._shutdown.is_set():
                return "shutdown"
            return state["stop"]

        started = time.time()
        try:
            await self.fabric.mark_running(rec)
        except FencedError:
            log.warning("%s: lost before it started", job_id)
            return
        hb = asyncio.create_task(beat())
        log.info("%s: running (attempt %d, fence %d): %s", job_id, attempt, rec["fence"], spec["command"]["argv"])
        try:
            body = await asyncio.to_thread(execute, spec, job_id, attempt, self.ctx, stop)
        except Fenced:
            log.warning("%s: fenced mid-run; abandoning without writing", job_id)
            return
        except InfraError as err:
            await self._infra(rec, msg, attempt, started, str(err))
            return
        except Exception as err:  # noqa: BLE001
            log.exception("%s: worker error", job_id)
            await self._infra(rec, msg, attempt, started, f"worker error: {type(err).__name__}: {err}")
            return
        finally:
            stop_event.set()
            await hb
        if state["stop"] == "fenced":
            return
        await self._complete(rec, msg, attempt, started, body)

    async def _infra(self, rec, msg, attempt, started, reason) -> None:
        job_id = rec["job_id"]
        log.warning("%s: infrastructure failure: %s", job_id, reason)
        try:
            if self._shutdown.is_set():
                await self.fabric.requeue(rec, msg, f"worker shutdown: {reason}", count_attempt=False, delay=2.0)
            elif attempt < rec["max_attempts"]:
                await self.fabric.requeue(rec, msg, reason, count_attempt=True, delay=5.0)
            else:
                body = {"host": self.inventory, "placement": {}, "runtime": {}, "build": [], "runs": [], "summary": None,
                        "artifacts": [], "verification": None, "inputs": [], "notes": [], "phases": {}, "logs": {},
                        "fidelity": {"policy": rec["spec"].get("fidelity", {}).get("policy", "standard"), "grade": "none",
                                     "contended": False, "checks": [], "violations": [], "tier": None},
                        "status": "infra_error", "error": reason}
                await self._complete(rec, msg, attempt, started, body)
        except FencedError as err:
            log.warning("%s: %s", job_id, err)

    async def _complete(self, rec, msg, attempt, started, body) -> None:
        job_id = rec["job_id"]
        result = {
            "schema": protocol.RESULT_SCHEMA_ID, "job_id": job_id, "spec_sha256": rec["spec_sha256"],
            "attempt": attempt, "fence": rec["fence"], "status": body["status"],
            "outcome_class": protocol.outcome_class(body["status"]), "error": body["error"],
            "worker": {"id": self.id, "hostname": socket.gethostname(), "labels": self.labels, "version": __version__},
            "host": body["host"], "placement": body["placement"], "runtime": body["runtime"], "fidelity": body["fidelity"],
            "timing": {"queued_at": rec["created_at"], "claimed_at": rec.get("claimed_at"), "started_at": started,
                       "finished_at": time.time(), "phases": body.get("phases", {})},
            "build": body["build"], "runs": body["runs"], "summary": body["summary"], "artifacts": body["artifacts"],
            "verification": body.get("verification"), "inputs": body.get("inputs", []), "notes": body.get("notes", []),
            "logs": body.get("logs", {}), "labels": rec["labels"],
        }
        mirror = self.state_dir / "results" / job_id
        mirror.mkdir(parents=True, exist_ok=True)
        import json
        (mirror / f"attempt-{attempt}-fence-{rec['fence']}.json").write_text(json.dumps(result, indent=1, sort_keys=True, default=str))
        try:
            await self.fabric.complete(rec, result, msg)
            log.info("%s: %s (grade %s)", job_id, body["status"], body["fidelity"].get("grade"))
        except FencedError as err:
            log.warning("%s: result refused: %s", job_id, err)


def install_signal_handlers(*workers: Worker) -> None:
    """SIGTERM and SIGINT stop every given worker (a handler per signal, so slots share one)."""
    loop = asyncio.get_running_loop()

    def stop_all(*_: Any) -> None:
        for w in workers:
            w.request_shutdown()
    for sig in (signal.SIGTERM, signal.SIGINT):
        try:
            loop.add_signal_handler(sig, stop_all)
        except NotImplementedError:
            signal.signal(sig, stop_all)


def default_worker_id() -> str:
    return f"{socket.gethostname()}"
