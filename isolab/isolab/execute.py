"""Run one job on this host: fetch, plan, reserve, build, measure, collect.

Nothing here talks to the hub. ``execute()`` takes a normalised spec, a
worker context (what this host has) and a ``stop()`` callback polled while
processes run, and returns the body of an ``isolab.result/v1`` document; the
worker adds identity and fence. Infrastructure failures raise
:class:`InfraError` (retried); losing the job raises :class:`Fenced`.
"""
from __future__ import annotations

import hashlib
import json
import logging
import os
import shutil
import statistics
import subprocess
import sys
import time
from contextlib import contextmanager
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable

from . import perfstat, protocol
from .blobs import LocalBlobCache, decode_content, elf_info, sha256_file
from .inventory import Paths
from .isolation import (Check, FidelityUnmet, QuiesceError, Reservation, grade, post_checks,
                        tier_from_mechanisms)
from .planner import PlanError, plan_cpus
from .runners import Backend, BackendError, ExecOutcome, JobContext, JobDirs
from .samplers import Sampler, cgroup_procs, cpu_snapshot, descendants

log = logging.getLogger("isolab.execute")
StopFn = Callable[[], str | None]


class InfraError(RuntimeError):
    """The environment failed, not the command. Retried; never evidence."""


class Fenced(RuntimeError):
    """This worker lost the job mid-run; nothing may be written."""


@dataclass
class WorkerContext:
    worker_id: str
    hostname: str
    labels: dict[str, str]
    version: str
    inventory: dict[str, Any]
    caps: dict[str, Any]
    lab_cpus: list[int]
    housekeeping: list[int]
    backends: dict[str, Backend]
    default_backend: str
    default_image: str | None
    oci_default: str | None
    launcher: str | None
    perf: str | None
    perf_probe: dict[str, Any]
    blob_cache: LocalBlobCache
    fetch_blob: Callable[[str, Path], None]
    put_artifact: Callable[[str, str, Path], str]
    git_cache_dir: Path
    jobs_dir: Path
    paths: Paths = field(default_factory=Paths)
    cgroup_root: str = "/sys/fs/cgroup"
    lock_path: str = "/tmp/isolab.lock"
    pull_policy: str = "missing"
    userns: str | None = None
    keep_job_dirs: bool = False
    progress: Callable[[dict[str, Any]], None] = lambda p: None
    live: dict[str, Any] = field(default_factory=dict)


class _Phases:
    def __init__(self) -> None:
        self.phases: dict[str, float] = {}

    @contextmanager
    def phase(self, name: str):
        t0 = time.perf_counter()
        try:
            yield
        finally:
            self.phases[name] = round(self.phases.get(name, 0.0) + time.perf_counter() - t0, 4)


# -- inputs ------------------------------------------------------------------

def _git(args: list[str], cwd: Path | None = None, timeout: float = 1800) -> str:
    try:
        p = subprocess.run(["git", *args], cwd=cwd, capture_output=True, text=True, timeout=timeout)
    except (subprocess.TimeoutExpired, OSError) as err:
        raise InfraError(f"git {' '.join(args[:2])}: {err}") from err
    if p.returncode != 0:
        raise InfraError(f"git {' '.join(args[:3])}: {p.stderr.strip()[-400:]}")
    return p.stdout.strip()


def export_git(url: str, commit: str, dest: Path, cache_dir: Path, sparse: list[str] | None) -> str:
    cache_dir.mkdir(parents=True, exist_ok=True)
    mirror = cache_dir / (hashlib.sha256(url.encode()).hexdigest()[:16] + ".git")
    if not mirror.exists():
        _git(["clone", "--mirror", "--filter=blob:none", url, str(mirror)])
    has = subprocess.run(["git", "cat-file", "-e", f"{commit}^{{commit}}"], cwd=mirror, capture_output=True).returncode == 0
    if not has:
        _git(["remote", "update", "--prune"], cwd=mirror)
        has = subprocess.run(["git", "cat-file", "-e", f"{commit}^{{commit}}"], cwd=mirror, capture_output=True).returncode == 0
    if not has:
        try:
            _git(["fetch", "--filter=blob:none", "origin", commit], cwd=mirror)
        except InfraError:
            pass
    resolved = _git(["rev-parse", f"{commit}^{{commit}}"], cwd=mirror)
    if resolved != commit:
        raise InfraError(f"commit {commit} resolved to {resolved}")
    dest.mkdir(parents=True, exist_ok=True)
    args = ["archive", "--format=tar", commit, *(sparse or [])]
    try:
        with subprocess.Popen(["git", *args], cwd=mirror, stdout=subprocess.PIPE) as ga:
            tar = subprocess.run(["tar", "-x", "-C", str(dest)], stdin=ga.stdout, capture_output=True, timeout=1800)
        if tar.returncode != 0:
            raise InfraError(f"tar extract: {tar.stderr.decode(errors='replace')[-300:]}")
    except (OSError, subprocess.TimeoutExpired) as err:
        raise InfraError(f"git archive: {err}") from err
    return resolved


def stage_inputs(spec: dict[str, Any], dirs: JobDirs, ctx: WorkerContext, notes: list[str]) -> list[dict[str, Any]]:
    staged = []
    host_arch = ctx.inventory.get("arch")
    for inp in spec["inputs"]:
        target = dirs.work / inp["path"]
        target.parent.mkdir(parents=True, exist_ok=True)
        entry: dict[str, Any] = {"path": inp["path"]}
        if "content" in inp:
            data = decode_content(inp)
            target.write_bytes(data)
            entry.update(kind="content", bytes=len(data))
        elif "sha256" in inp:
            digest = inp["sha256"]
            cached = ctx.blob_cache.path_for(digest)
            if not ctx.blob_cache.has(digest):
                try:
                    ctx.fetch_blob(digest, cached)
                except Exception as err:
                    raise InfraError(f"fetch blob {digest[:12]} for {inp['path']}: {err}") from err
                if sha256_file(cached) != digest:
                    cached.unlink(missing_ok=True)
                    raise InfraError(f"blob {digest[:12]} failed verification after download")
            shutil.copyfile(cached, target)
            entry.update(kind="blob", sha256=digest, bytes=target.stat().st_size)
        elif "git" in inp:
            g = inp["git"]
            resolved = export_git(g["url"], g["commit"], target, ctx.git_cache_dir, g.get("sparse_paths"))
            entry.update(kind="git", url=g["url"], commit=resolved)
        if inp.get("mode") and target.is_file():
            os.chmod(target, int(inp["mode"], 8))
            if int(inp["mode"], 8) & 0o111:
                info = elf_info(target)
                if info:
                    entry["elf"] = info
                    if host_arch and info["arch"] != host_arch:
                        notes.append(f"input {inp['path']} is an {info['arch']} ELF binary; this worker is {host_arch}")
                    if not info["static"]:
                        notes.append(f"input {inp['path']} is dynamically linked ({info['interpreter']}); it needs a matching libc in the image")
        staged.append(entry)
    return staged


# -- logs and artifacts ------------------------------------------------------

def _cap_log(path: Path, cap: int) -> bool:
    if not path.exists():
        return False
    size = path.stat().st_size
    if size <= cap:
        return False
    half = cap // 2
    with open(path, "rb") as fh:
        head = fh.read(half)
        fh.seek(size - half)
        tail = fh.read()
    path.write_bytes(head + f"\n[isolab: {size - 2 * half} bytes elided]\n".encode() + tail)
    return True


def _tail(path: Path, n: int = 2048) -> str:
    try:
        size = path.stat().st_size
        with open(path, "rb") as fh:
            fh.seek(max(0, size - n))
            return fh.read().decode(errors="replace")
    except OSError:
        return ""


def collect_artifacts(job_id: str, src: Path, ctx: WorkerContext, max_bytes: int, max_files: int,
                      notes: list[str]) -> list[dict[str, Any]]:
    out = []
    total = 0
    files = sorted(p for p in src.rglob("*") if p.is_file())
    for p in files:
        rel = p.relative_to(src).as_posix()
        size = p.stat().st_size
        entry: dict[str, Any] = {"path": rel, "sha256": sha256_file(p), "bytes": size}
        if len(out) >= max_files or total + size > max_bytes:
            entry["stored"] = False
            entry["reason"] = "artifact budget exceeded"
        else:
            try:
                entry["key"] = ctx.put_artifact(job_id, rel, p)
                entry["stored"] = True
                total += size
            except Exception as err:
                entry["stored"] = False
                entry["reason"] = f"upload failed: {err}"
        out.append(entry)
    if any(not e.get("stored") for e in out):
        notes.append("some artifacts were not stored; see artifacts[].reason")
    return out


# -- stats -------------------------------------------------------------------

def _stats(xs: list[float]) -> dict[str, float] | None:
    xs = [x for x in xs if x is not None]
    if not xs:
        return None
    mean = statistics.fmean(xs)
    sd = statistics.stdev(xs) if len(xs) > 1 else 0.0
    return {"n": len(xs), "min": min(xs), "median": statistics.median(xs), "mean": mean, "max": max(xs),
            "stdev": sd, "cv": sd / mean if mean else 0.0}


def summarize_runs(runs: list[dict[str, Any]]) -> dict[str, Any] | None:
    timed = [r for r in runs if not r["warmup"] and r["exit_code"] == 0]
    if not timed:
        return None
    clean = [r for r in timed if not r.get("contended")]
    used = clean or timed
    out: dict[str, Any] = {"n": len(used), "n_timed": len(timed), "n_contended": len(timed) - len(clean),
                           "used_indices": [r["index"] for r in used], "uses_contended": not clean,
                           "wall_s": _stats([r["wall_s"] for r in used]),
                           "cpu_s": _stats([(r.get("rusage") or {}).get("user_s", 0) + (r.get("rusage") or {}).get("sys_s", 0)
                                            for r in used if r.get("rusage")])}
    for key in ("instructions", "cycles", "task_clock_ms", "context_switches", "page_faults"):
        vals = [(r.get("counters") or {}).get(key) for r in used]
        out[key] = _stats([v for v in vals if v is not None]) if any(v is not None for v in vals) else None
    ipc = [(r.get("counters") or {}).get("ipc") for r in used]
    out["ipc"] = _stats([v for v in ipc if v is not None]) if any(v is not None for v in ipc) else None
    peaks = [(r.get("cgroup") or {}).get("memory_peak_bytes") for r in used]
    out["memory_peak_bytes_max"] = max([p for p in peaks if p is not None], default=None)
    return out


# -- verification ------------------------------------------------------------

_VERIFY_EXIT = {0: "verified", 1: "refuted", 2: "no_claim"}


def verify_run(vspec: dict[str, Any], backend: Backend, jctx: JobContext, run_out: Path, cwd: str,
               env: dict[str, str], logs: Path, index: int, stop: StopFn) -> dict[str, Any]:
    cert = run_out / vspec["certificate_file"]
    if vspec.get("builtin") == "taskq-certificate":
        if not cert.exists():
            return {"status": "no_claim", "verifier": "taskq.verify", "detail": f"no {vspec['certificate_file']} written"}
        import importlib.util
        if importlib.util.find_spec("taskq") is None:
            return {"status": "error", "verifier": "taskq.verify", "detail": "taskq is not installed on this worker"}
        t0 = time.perf_counter()
        try:
            p = subprocess.run([sys.executable, "-m", "taskq.verify", str(cert)], capture_output=True, text=True,
                               timeout=vspec["timeout_s"], stdin=subprocess.DEVNULL)
        except subprocess.TimeoutExpired:
            return {"status": "error", "verifier": "taskq.verify", "detail": f"exceeded {vspec['timeout_s']}s"}
        last = (p.stdout.strip().splitlines() or [""])[-1]
        try:
            detail: Any = json.loads(last)
        except ValueError:
            detail = (p.stdout[-1000:] + p.stderr[-1000:]) or None
        return {"status": _VERIFY_EXIT.get(p.returncode, "error"), "verifier": "taskq.verify", "exit_code": p.returncode,
                "wall_s": time.perf_counter() - t0, "certificate_sha256": sha256_file(cert), "detail": detail}
    venv = {**env, "ISOLAB_CERTIFICATE": backend.paths_for(jctx, cert)}
    out = backend.exec(jctx, vspec["argv"], cwd, venv, vspec["timeout_s"], logs / f"verify-{index}.stdout",
                       logs / f"verify-{index}.stderr", stop)
    status = "error" if out.timed_out or out.exit_code is None else _VERIFY_EXIT.get(out.exit_code, "error")
    return {"status": status, "verifier": " ".join(vspec["argv"]), "exit_code": out.exit_code, "wall_s": out.wall_host_s,
            "certificate_sha256": sha256_file(cert) if cert.exists() else None,
            "detail": _tail(logs / f"verify-{index}.stdout", 1000) or None}


def rollup_verification(runs: list[dict[str, Any]]) -> dict[str, Any]:
    counts = {"verified": 0, "refuted": 0, "no_claim": 0, "error": 0}
    for r in runs:
        v = r.get("verification")
        if v:
            counts[v["status"]] = counts.get(v["status"], 0) + 1
    status = "refuted" if counts["refuted"] else "error" if counts["error"] else "verified" if counts["verified"] else "no_claim"
    return {"status": status, "counts": counts}


# -- the job -----------------------------------------------------------------

def execute(spec: dict[str, Any], job_id: str, attempt: int, ctx: WorkerContext, stop: StopFn) -> dict[str, Any]:
    ph = _Phases()
    notes: list[str] = []
    res_block: dict[str, Any] = {"backend": None, "oci_runtime": None, "image": None, "image_digest": None}
    body: dict[str, Any] = {"host": ctx.inventory, "placement": {}, "runtime": res_block,
                            "fidelity": {"policy": spec["fidelity"]["policy"], "grade": "none", "contended": False,
                                         "checks": [], "violations": [], "tier": None},
                            "build": [], "runs": [], "summary": None, "artifacts": [], "verification": None,
                            "inputs": [], "notes": notes, "phases": ph.phases, "logs": {}}
    dirs = JobDirs.create(ctx.jobs_dir / job_id / f"attempt-{attempt}")
    status, error = "succeeded", None
    fid = spec["fidelity"]
    limits, measure, resources = spec["limits"], spec["measure"], spec["resources"]
    try:
        ctx.progress({"phase": "fetch"})
        with ph.phase("fetch"):
            body["inputs"] = stage_inputs(spec, dirs, ctx, notes)
        _check_stop(stop)

        backend_name = spec["runtime"]["backend"] if spec["runtime"]["backend"] != "auto" else ctx.default_backend
        backend = ctx.backends.get(backend_name)
        if backend is None:
            raise InfraError(f"backend {backend_name!r} is not available on {ctx.worker_id}")
        oci = spec["runtime"]["oci_runtime"]
        if oci == "auto":
            oci = ctx.oci_default if backend.container else None
        image_info: dict[str, Any] = {}
        image = None
        if backend.container:
            image = spec["runtime"]["image"] or ctx.default_image
            if not image:
                raise InfraError("no image in the spec and the worker has no default image")
            ctx.progress({"phase": "image", "image": image})
            with ph.phase("image"):
                try:
                    image_info = backend.resolve_image(image, ctx.pull_policy)
                except BackendError as err:
                    raise InfraError(str(err)) from err
        res_block.update(backend=backend.name, oci_runtime=oci, image=image, image_digest=image_info.get("digest"),
                         image_id=image_info.get("id"), backend_version=getattr(backend, "version", lambda: None)(),
                         launcher=bool(ctx.launcher), perf=ctx.perf_probe)
        _check_stop(stop)

        topo = ctx.inventory["topology"]
        try:
            placement = plan_cpus(topo, ctx.lab_cpus, resources["cpus"], resources["numa_node"], resources["smt"])
        except PlanError as err:
            raise InfraError(f"placement: {err}") from err
        gpu_ids = [g["index"] for g in (ctx.inventory.get("gpus") or [])][: resources["gpus"]]
        if len(gpu_ids) < resources["gpus"]:
            raise InfraError(f"wants {resources['gpus']} gpu(s); worker has {len(gpu_ids)}")
        placement.gpus = gpu_ids

        ctx.progress({"phase": "reserve", "cpus": placement.cpus})
        with Reservation(job_id, placement, fid, ctx.caps, ctx.housekeeping, resources["memory_mb"], resources["pids"],
                         paths=ctx.paths, cgroup_root=ctx.cgroup_root, lock_path=ctx.lock_path) as res:
            with ph.phase("reserve"):
                pre = res.pre_checks(ctx.inventory)
                try:
                    res.enforce_requirements(pre)
                except FidelityUnmet as err:
                    raise InfraError(f"fidelity requirements unmet on {ctx.worker_id}: {err}") from err
            body["placement"] = {**placement.to_dict(), "housekeeping_cpus": ctx.housekeeping, "isolation": res.mechanisms}
            tier = res.mechanisms.get("tier") or tier_from_mechanisms(res.mechanisms, ctx.caps)
            body["fidelity"]["tier"] = tier
            ctx.progress({"phase": "quiesce", "settle_s": fid["settle_s"]})
            with ph.phase("quiesce"):
                try:
                    settle = res.quiesce(fid["settle_s"], retries=5 if fid["policy"] != "best_effort" else 1)
                except QuiesceError as err:
                    raise InfraError(str(err)) from err
            body["fidelity"]["checks"] = [c.to_dict() for c in pre + settle]

            jctx = JobContext(job_id=job_id, dirs=dirs, cpus=placement.cpus, mems=placement.mems,
                              cgroup_dir=res.job_cgroup, cgroup_parent=res.mechanisms.get("cgroup_parent"),
                              image=image, oci_runtime=oci, network=spec["runtime"]["network"],
                              memory_mb=resources["memory_mb"], pids=resources["pids"], scratch_mb=resources["scratch_mb"],
                              user=spec["runtime"]["user"], extra_args=spec["runtime"]["extra_args"], gpus=gpu_ids,
                              launcher=ctx.launcher, aslr_off=fid.get("aslr") == "off", userns=ctx.userns,
                              env_base=protocol.host_env_passthrough() if not backend.container else {})
            ctx.progress({"phase": "start"})
            with ph.phase("start"):
                try:
                    start_info = backend.start(jctx)
                except BackendError as err:
                    raise InfraError(str(err)) from err
            res_block["container"] = {k: v for k, v in start_info.items() if k != "run_argv"}
            res_block["run_argv"] = start_info.get("run_argv")
            controls = start_info.get("cgroup_controls")
            if controls is not None and not controls.get("cpuset"):
                notes.append("container engine has no cpuset controller here (rootless?): the container was not pinned "
                             "to the reserved cpus by cgroup; only the host-side reservation applies")
                body["fidelity"]["tier"] = "D"
            try:
                # ---- build (untimed, measured, reported separately)
                with ph.phase("build"):
                    for i, step in enumerate(spec["build"]):
                        ctx.progress({"phase": "build", "step": i})
                        env = {**step["env"], **_contract(backend, jctx, job_id, attempt, None, None, dirs.out / "build")}
                        (dirs.out / "build").mkdir(exist_ok=True)
                        ctx.live["stdout_log"] = dirs.logs / f"build-{i}.stdout"
                        out = backend.exec(jctx, step["argv"], step["cwd"], env, step["timeout_s"] or limits["build_timeout_s"],
                                           dirs.logs / f"build-{i}.stdout", dirs.logs / f"build-{i}.stderr", stop)
                        body["build"].append({"index": i, "argv": step["argv"], "exit_code": out.exit_code, "signal": out.signal,
                                              "timed_out": out.timed_out, "wall_s": (out.launcher or {}).get("wall_s", out.wall_host_s),
                                              "wall_host_s": out.wall_host_s, "rusage": out.launcher})
                        status, error = _outcome_status(out, f"build step {i}", step["timeout_s"] or limits["build_timeout_s"])
                        if status != "succeeded":
                            break
                # ---- measure
                if status == "succeeded":
                    status, error = _measure(spec, job_id, attempt, ctx, backend, jctx, res, dirs, body, ph, stop)
            finally:
                with ph.phase("stop"):
                    backend.stop(jctx)
        # reservation released: the host is free while we hash and upload
        ctx.progress({"phase": "collect"})
        with ph.phase("collect"):
            body["summary"] = summarize_runs(body["runs"])
            if spec.get("verify"):
                body["verification"] = rollup_verification(body["runs"])
            for p in list(dirs.logs.iterdir()):
                if _cap_log(p, limits["max_log_bytes"]):
                    p.rename(p.with_name(p.name + ".truncated"))
            body["logs"] = {"stdout_tail": _tail(dirs.logs / "stdout.log"), "stderr_tail": _tail(dirs.logs / "stderr.log")}
            shutil.move(str(dirs.logs), str(dirs.out / "_isolab_logs"))
            body["artifacts"] = collect_artifacts(job_id, dirs.out, ctx, limits["max_artifact_bytes"],
                                                  limits["max_artifact_files"], notes)
        fidelity = body["fidelity"]
        violations: set[str] = {c["name"] for c in fidelity["checks"] if c["status"] == "fail"}
        used = set((body["summary"] or {}).get("used_indices") or [])
        any_contended_used = False
        for r in body["runs"]:
            for c in r.get("checks") or []:
                if c["status"] == "fail":
                    violations.add(f"run{r['index']}:{c['name']}")
            if r["index"] in used and r.get("contended"):
                any_contended_used = True
        fidelity["violations"] = sorted(violations)
        fidelity["contended"] = any(r.get("contended") for r in body["runs"] if not r["warmup"])
        fidelity["grade"] = grade(fidelity.get("tier") or "D", any_contended_used)
    finally:
        if not ctx.keep_job_dirs:
            shutil.rmtree(dirs.root, ignore_errors=True)
    return {**body, "status": status, "error": error}


def _check_stop(stop: StopFn) -> None:
    s = stop()
    if s == "fenced":
        raise Fenced()
    if s == "shutdown":
        raise InfraError("worker shutdown")


def _contract(backend: Backend, jctx: JobContext, job_id: str, attempt: int, repeat: int | None,
              warmup: bool | None, out_dir: Path) -> dict[str, str]:
    out_dir.mkdir(parents=True, exist_ok=True)
    os.chmod(out_dir, 0o1777)
    return protocol.env_contract(job_id, attempt, repeat, warmup, backend.paths_for(jctx, out_dir),
                                 backend.paths_for(jctx, jctx.dirs.work), backend.paths_for(jctx, jctx.dirs.scratch))


def _outcome_status(out: ExecOutcome, what: str, timeout: float) -> tuple[str, str | None]:
    if out.stopped == "fenced":
        raise Fenced()
    if out.stopped == "shutdown":
        raise InfraError(f"worker shutdown interrupted {what}")
    if out.stopped == "cancel":
        return "cancelled", f"cancelled during {what}"
    if out.timed_out:
        return "timeout", f"{what} exceeded {timeout}s"
    if out.exit_code != 0:
        return "failed", f"{what} exited {out.exit_code}" + (f" (signal {out.signal})" if out.signal else "")
    return "succeeded", None


def _measure(spec, job_id, attempt, ctx, backend, jctx, res, dirs, body, ph, stop) -> tuple[str, str | None]:
    fid, measure, limits, cmd = spec["fidelity"], spec["measure"], spec["limits"], spec["command"]
    plan = [True] * measure["warmups"] + [False] * measure["repeats"]
    strict = fid["policy"] == "strict"
    extras_allowed = measure["repeats"] if strict else 0
    i, clean, extras = 0, 0, 0
    perf_ok = bool(ctx.perf and ctx.perf_probe.get("available"))
    status, error = "succeeded", None
    with ph.phase("measure"):
        while i < len(plan) or (strict and clean < measure["repeats"] and extras < extras_allowed):
            warm = plan[i] if i < len(plan) else False
            is_extra = i >= len(plan)
            run_out = dirs.out / f"run-{i}"
            meta = run_out / "_isolab"
            meta.mkdir(parents=True, exist_ok=True)
            os.chmod(run_out, 0o1777)
            os.chmod(meta, 0o1777)
            env = {**cmd["env"], **_contract(backend, jctx, job_id, attempt, i, warm, run_out)}
            ctx.progress({"phase": "measure", "repeat": i, "of": len(plan), "warmup": warm})
            ctx.live["stdout_log"] = dirs.logs / "stdout.log"
            if i > 0 and measure["cooldown_s"]:
                time.sleep(measure["cooldown_s"])
            if measure["drop_caches"] and ctx.caps.get("drop_caches"):
                try:
                    os.sync()
                    Path("/proc/sys/vm/drop_caches").write_text("3")
                except OSError as err:
                    body["notes"].append(f"drop_caches failed: {err}")
            cg_before = backend.cgroup_stats(jctx) or {}
            perf_out = meta / "perf.csv"
            perf_prefix = perfstat.argv_prefix(ctx.perf, measure["perf_events"], jctx.cpus, str(perf_out)) if perf_ok else None
            exclude = _exclude_fn(jctx)
            sampler = Sampler(res.placement.reserved, period=measure["sample_period_s"], paths=ctx.paths,
                              exclude=exclude, cgroup=jctx.cgroup_dir, other_cpu_threshold=fid["max_other_cpu"])
            sampler.start()
            out = backend.exec(jctx, cmd["argv"], cmd["cwd"], env, limits["timeout_s"], dirs.logs / "stdout.log",
                               dirs.logs / "stderr.log", stop, cmd.get("stdin"), meta / "launch.json", perf_prefix)
            sampler.stop()
            cg_after = backend.cgroup_stats(jctx) or {}
            counters = None
            if perf_out.exists():
                parsed = perfstat.parse_csv(perf_out.read_text())
                counters = perfstat.counters(parsed)
            summ = sampler.summary()
            with open(meta / "samples.jsonl", "w") as fh:
                for s in sampler.samples:
                    fh.write(json.dumps(s, default=str) + "\n")
            cg = _delta(cg_before, cg_after)
            checks = post_checks(summ, counters, cg, out.launcher, fid, res.mechanisms, out.wall_host_s)
            contended = any(c.status == "fail" for c in checks)
            metrics = None
            mfile = run_out / "metrics.json"
            if mfile.exists():
                try:
                    metrics = json.loads(mfile.read_text())
                except ValueError as e:
                    metrics = {"_isolab_error": f"metrics.json unparseable: {e}"}
            launcher = out.launcher or {}
            run = {"index": i, "warmup": warm, "extra": is_extra, "exit_code": out.exit_code, "signal": out.signal,
                   "timed_out": out.timed_out, "wall_s": launcher.get("wall_s", out.wall_host_s),
                   "wall_host_s": out.wall_host_s,
                   "timing_source": "launcher" if out.launcher and launcher.get("source") != "host-wait4" else "host-wait4",
                   "rusage": out.launcher, "cgroup": cg, "counters": counters, "conditions": summ,
                   "checks": [c.to_dict() for c in checks], "contended": contended, "metrics": metrics}
            body["runs"].append(run)
            status, error = _outcome_status(out, f"repeat {i}", limits["timeout_s"])
            if status != "succeeded":
                break
            if spec.get("verify") and not warm:
                run["verification"] = verify_run(spec["verify"], backend, jctx, run_out, cmd["cwd"], env, dirs.logs, i, stop)
            if not warm and not contended:
                clean += 1
            if is_extra:
                extras += 1
            i += 1
    return status, error


def _exclude_fn(jctx: JobContext):
    cg = jctx.cgroup_dir

    def exclude() -> set[int]:
        pids = cgroup_procs(cg) if cg is not None else set()
        if not pids:
            snap = cpu_snapshot()
            pids = descendants(os.getpid(), snap)
        return pids
    return exclude


def _delta(before: dict[str, Any], after: dict[str, Any]) -> dict[str, Any] | None:
    if not after:
        return None
    out: dict[str, Any] = {}
    for k, v in after.items():
        if isinstance(v, int) and isinstance(before.get(k), int) and k not in ("memory_peak_bytes", "pids_peak"):
            out[k] = v - before[k]
        else:
            out[k] = v
    if "usage_usec" in out:
        out["usage_s"] = out["usage_usec"] / 1e6
    if "user_usec" in out:
        out["user_s"] = out["user_usec"] / 1e6
    if "system_usec" in out:
        out["system_s"] = out["system_usec"] / 1e6
    return out
