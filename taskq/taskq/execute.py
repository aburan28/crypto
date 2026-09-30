"""Check out a pinned commit, run a spec, and measure it.

Nothing here talks to Redis. `execute()` takes a normalised spec and a
`stop()` callback (polled about every 100 ms) and returns the body of a
taskq.task-result/v1 document; the worker adds identity and fence.

Measurement is per child process tree via wait4(2): wall clock from the
parent, user/sys CPU and peak RSS from the kernel's rusage for the child and
every descendant it waited for. Setup steps are measured the same way and
reported separately, never folded into the timed runs.
"""
from __future__ import annotations

import hashlib
import json
import os
import platform
import resource
import shutil
import signal
import statistics
import subprocess
import sys
import tempfile
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable

StopFn = Callable[[], str | None]  # None, "cancel", "fenced" or "shutdown"


class InfraError(RuntimeError):
    """The environment failed, not the command. Retried; never evidence."""


class Fenced(RuntimeError):
    """This worker lost the task mid-run; nothing may be written."""


# -- source checkout ---------------------------------------------------------

def _git(*args: str, cwd: str | Path | None = None, timeout: float = 1800) -> str:
    try:
        p = subprocess.run(["git", *args], cwd=cwd, capture_output=True,
                           text=True, timeout=timeout)
    except subprocess.TimeoutExpired as err:
        raise InfraError(f"git {' '.join(args[:2])} timed out") from err
    if p.returncode != 0:
        raise InfraError(f"git {' '.join(args)}: {p.stderr.strip()[-500:]}")
    return p.stdout.strip()


class RepoCache:
    """A blobless mirror per repo alias, plus a git worktree per commit.

    The mirror is a partial clone (`--filter=blob:none`): history is cheap,
    and only the blobs a checked-out commit needs are ever downloaded. Trees
    are worktrees of the mirror, reused across tasks at the same commit (and
    sparse set) so build caches such as cargo `target/` survive, and every use
    starts with `git reset --hard`: tracked files are always exactly the
    commit. The least recently used trees are removed beyond `max_trees`.
    """

    def __init__(self, root: str | Path, repos: dict[str, str], max_trees: int = 8):
        self.root = Path(root)
        self.repos = repos
        self.max_trees = max_trees

    def prepare(self, alias: str, commit: str, sparse: list[str] | None = None) -> Path:
        if alias not in self.repos:
            raise InfraError(f"repo {alias!r} is not in this worker's allowlist "
                             f"({', '.join(sorted(self.repos)) or 'empty'})")
        mirror = self.root / "mirrors" / f"{alias}.git"
        if not mirror.exists():
            mirror.parent.mkdir(parents=True, exist_ok=True)
            _git("clone", "--mirror", "--filter=blob:none", self.repos[alias], str(mirror))
        if not self._has(mirror, commit):
            _git("remote", "update", "--prune", cwd=mirror)
        if not self._has(mirror, commit):
            try:  # an unadvertised commit (e.g. an unmerged PR head)
                _git("fetch", "--filter=blob:none", "origin", commit, cwd=mirror)
            except InfraError:
                pass
        if not self._has(mirror, commit):
            raise InfraError(f"commit {commit} not found in {alias}")
        name = commit
        if sparse:
            name += "-" + hashlib.sha256("\0".join(sorted(sparse)).encode()).hexdigest()[:12]
        tree = self.root / "trees" / alias / name
        if not (tree / ".git").exists():
            tree.parent.mkdir(parents=True, exist_ok=True)
            _git("worktree", "prune", cwd=mirror)
            _git("worktree", "add", "--no-checkout", "--detach", str(tree), commit, cwd=mirror)
            if sparse:
                _git("sparse-checkout", "set", "--cone", *sparse, cwd=tree)
        _git("reset", "--hard", commit, cwd=tree)
        os.utime(tree)
        self._evict(mirror, alias, keep=tree)
        return tree

    def _evict(self, mirror: Path, alias: str, keep: Path) -> None:
        trees = sorted((p for p in (self.root / "trees" / alias).iterdir() if p != keep),
                       key=lambda p: p.stat().st_mtime)
        for old in trees[: max(0, len(trees) + 1 - self.max_trees)]:
            subprocess.run(["git", "worktree", "remove", "--force", str(old)],
                           cwd=mirror, capture_output=True)
            shutil.rmtree(old, ignore_errors=True)

    @staticmethod
    def _has(mirror: Path, commit: str) -> bool:
        return subprocess.run(["git", "cat-file", "-e", f"{commit}^{{commit}}"],
                              cwd=mirror, capture_output=True).returncode == 0


# -- environment -------------------------------------------------------------

def _read(path: str) -> str | None:
    try:
        return Path(path).read_text().strip()
    except OSError:
        return None


def capture_environment() -> dict[str, Any]:
    cpu_model = None
    info = _read("/proc/cpuinfo") or ""
    for line in info.splitlines():
        if line.startswith("model name"):
            cpu_model = line.split(":", 1)[1].strip()
            break
    try:
        affinity = sorted(os.sched_getaffinity(0))
    except AttributeError:
        affinity = None
    return {
        "hostname": platform.node(),
        "platform": platform.platform(),
        "kernel": platform.release(),
        "machine": platform.machine(),
        "python": sys.version.split()[0],
        "cpu_model": cpu_model,
        "cpu_count": os.cpu_count(),
        "cpu_affinity": affinity,
        "loadavg": list(os.getloadavg()),
        "cpu_governor": _read("/sys/devices/system/cpu/cpu0/cpufreq/scaling_governor"),
        "cgroup_cpu_max": _read("/sys/fs/cgroup/cpu.max"),
        "cgroup_memory_max": _read("/sys/fs/cgroup/memory.max"),
        "k8s_node": os.environ.get("NODE_NAME"),
        "k8s_pod": os.environ.get("POD_NAME"),
        "image": os.environ.get("TASKQ_IMAGE"),
    }


# -- one measured process ----------------------------------------------------

@dataclass
class Measured:
    exit_code: int | None
    signal: int | None
    wall_seconds: float
    user_cpu_seconds: float
    sys_cpu_seconds: float
    max_rss_kb: int
    timed_out: bool
    stopped: str | None  # "cancel" / "fenced" when interrupted


def run_measured(argv: list[str], cwd: Path, env: dict[str, str],
                 timeout: float, stdout: Path, stderr: Path, stop: StopFn,
                 cpus: list[int] | None = None,
                 memory_mb: int | None = None) -> Measured:
    def preexec() -> None:
        if cpus:
            os.sched_setaffinity(0, cpus)
        if memory_mb:
            lim = memory_mb * 1024 * 1024
            resource.setrlimit(resource.RLIMIT_AS, (lim, lim))

    with open(stdout, "ab") as out, open(stderr, "ab") as err:
        t0 = time.perf_counter()
        try:
            proc = subprocess.Popen(argv, cwd=cwd, env=env, stdout=out,
                                    stderr=err, stdin=subprocess.DEVNULL,
                                    start_new_session=True, preexec_fn=preexec)
        except OSError as e:
            # e.g. argv[0] not found: the spec is wrong for this commit.
            err.write(f"taskq: cannot exec {argv[0]!r}: {e}\n".encode())
            return Measured(127, None, 0.0, 0.0, 0.0, 0, False, None)
        timed_out, stopped = False, None
        while True:
            pid, status, ru = os.wait4(proc.pid, os.WNOHANG)
            if pid:
                break
            elapsed = time.perf_counter() - t0
            if elapsed > timeout:
                timed_out = True
            else:
                stopped = stop()
            if timed_out or stopped:
                status, ru = _kill_group(proc.pid)
                break
            time.sleep(0.05 if elapsed < 2 else 0.1)
        wall = time.perf_counter() - t0
    proc.returncode = 0  # reaped by wait4; stop Popen from waiting again
    code = os.waitstatus_to_exitcode(status)
    return Measured(
        exit_code=code if code >= 0 else None,
        signal=-code if code < 0 else None,
        wall_seconds=wall, user_cpu_seconds=ru.ru_utime,
        sys_cpu_seconds=ru.ru_stime, max_rss_kb=ru.ru_maxrss,
        timed_out=timed_out, stopped=stopped)


def _kill_group(pid: int):
    """SIGTERM the process group, SIGKILL after a grace period; reap with rusage."""
    for sig, grace in ((signal.SIGTERM, 5.0), (signal.SIGKILL, None)):
        try:
            os.killpg(pid, sig)
        except ProcessLookupError:
            pass
        if grace is None:
            _, status, ru = os.wait4(pid, 0)
            return status, ru
        deadline = time.monotonic() + grace
        while time.monotonic() < deadline:
            got, status, ru = os.wait4(pid, os.WNOHANG)
            if got:
                # the leader is gone; take down any stragglers in its group
                try:
                    os.killpg(pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
                return status, ru
            time.sleep(0.05)


# -- artifacts ---------------------------------------------------------------

def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _cap_log(path: Path, cap: int) -> bool:
    """Keep the head and tail of an oversized log. Returns True if cut."""
    size = path.stat().st_size
    if size <= cap:
        return False
    half = cap // 2
    with open(path, "rb") as fh:
        head = fh.read(half)
        fh.seek(size - half)
        tail = fh.read()
    marker = f"\n[taskq: {size - 2 * half} bytes elided]\n".encode()
    path.write_bytes(head + marker + tail)
    return True


def collect_artifacts(src: Path, dest: Path | None) -> list[dict[str, Any]]:
    out = []
    for p in sorted(x for x in src.rglob("*") if x.is_file()):
        rel = p.relative_to(src).as_posix()
        entry = {"path": rel, "sha256": _sha256(p), "bytes": p.stat().st_size}
        if dest is not None:
            target = dest / rel
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(p, target)
            entry["uri"] = target.as_posix()
        out.append(entry)
    return out


# -- certificate verification ------------------------------------------------

_VERIFY_EXIT = {0: "verified", 1: "refuted", 2: "no_claim"}


def verify_run(vspec: dict[str, Any], run_out: Path, cwd: Path,
               env: dict[str, str]) -> dict[str, Any]:
    """Check one run's certificate in a separate process, after timing stopped.

    Exit 0/1/2 of the verifier mean verified/refuted/no_claim; anything else,
    a timeout or unparseable output is `error`, which is not a refutation.
    """
    cert = run_out / vspec["certificate_file"]
    if "builtin" in vspec:
        if not cert.exists():
            return {"status": "no_claim", "verifier": "taskq.verify",
                    "detail": f"no {vspec['certificate_file']} written"}
        argv = [sys.executable, "-m", "taskq.verify", str(cert)]
        name = "taskq.verify"
    else:
        argv, name = vspec["argv"], " ".join(vspec["argv"])
    venv = {**env, "TASKQ_CERTIFICATE": str(cert)}
    if "builtin" in vspec:  # the built-in must resolve to this worker's taskq
        venv["PYTHONPATH"] = os.pathsep.join(
            [str(Path(__file__).resolve().parent.parent)]
            + ([env["PYTHONPATH"]] if env.get("PYTHONPATH") else []))
    t0 = time.perf_counter()
    try:
        p = subprocess.run(argv, cwd=cwd, env=venv, capture_output=True, text=True,
                           timeout=vspec["timeout_seconds"], stdin=subprocess.DEVNULL)
    except subprocess.TimeoutExpired:
        return {"status": "error", "verifier": name,
                "detail": f"verifier exceeded {vspec['timeout_seconds']}s"}
    except OSError as err:
        return {"status": "error", "verifier": name, "detail": f"cannot exec verifier: {err}"}
    out = {"status": _VERIFY_EXIT.get(p.returncode, "error"), "verifier": name,
           "exit_code": p.returncode, "wall_seconds": time.perf_counter() - t0,
           "certificate_sha256": _sha256(cert) if cert.exists() else None}
    last = (p.stdout.strip().splitlines() or [""])[-1]
    try:
        out["detail"] = json.loads(last)
    except ValueError:
        out["detail"] = (p.stdout[-2000:] + p.stderr[-2000:]) or None
    return out


def rollup_verification(runs: list[dict[str, Any]]) -> dict[str, Any]:
    counts = {"verified": 0, "refuted": 0, "no_claim": 0, "error": 0}
    for r in runs:
        v = r.get("verification")
        if v:
            counts[v["status"]] += 1
    if counts["refuted"]:
        status = "refuted"
    elif counts["error"]:
        status = "error"
    elif counts["verified"]:
        status = "verified"
    else:
        status = "no_claim"
    return {"status": status, "counts": counts}


# -- the whole task ----------------------------------------------------------

def _summary(runs: list[dict[str, Any]]) -> dict[str, Any] | None:
    timed = [r for r in runs if not r["warmup"] and r["exit_code"] == 0]
    if not timed:
        return None
    def stats(field: str) -> dict[str, float]:
        xs = [r[field] for r in timed]
        return {"min": min(xs), "median": statistics.median(xs),
                "mean": statistics.fmean(xs), "max": max(xs),
                "stdev": statistics.stdev(xs) if len(xs) > 1 else 0.0}
    return {"n": len(timed), "wall_seconds": stats("wall_seconds"),
            "user_cpu_seconds": stats("user_cpu_seconds"),
            "sys_cpu_seconds": stats("sys_cpu_seconds"),
            "max_rss_kb": max(r["max_rss_kb"] for r in timed)}


def execute(spec: dict[str, Any], task_id: str, attempt: int,
            repos: RepoCache, artifact_root: Path | None, stop: StopFn,
            worker_cpus: list[int] | None = None) -> dict[str, Any]:
    """Run a normalised spec. Raises InfraError (retry) or Fenced (abandon)."""
    src = spec["source"]
    body: dict[str, Any] = {
        "environment": capture_environment(),
        "source": {"repo": src["repo"], "commit": src["commit"],
                   "sparse_paths": src.get("sparse_paths"), "resolved_commit": None},
        "setup": [], "runs": [], "summary": None, "artifacts": [],
        "error": None,
    }
    tree = repos.prepare(src["repo"], src["commit"], src.get("sparse_paths"))
    resolved = _git("rev-parse", "HEAD", cwd=tree)
    if resolved != src["commit"]:
        raise InfraError(f"checkout resolved to {resolved}, not {src['commit']}")
    body["source"]["resolved_commit"] = resolved

    cmd, limits = spec["command"], spec["limits"]
    cwd = (tree / cmd["cwd"]).resolve()
    if not cwd.is_dir() or tree.resolve() not in (cwd, *cwd.parents):
        return {**body, "status": "failed",
                "error": f"cwd {cmd['cwd']!r} does not exist at this commit"}
    cpus = None
    want = spec["placement"].get("cpus")
    if want:
        pool = worker_cpus or sorted(os.sched_getaffinity(0))
        if want > len(pool):
            raise InfraError(f"task wants {want} cpus; worker has {len(pool)}")
        cpus = pool[:want]

    work = Path(tempfile.mkdtemp(prefix=f"taskq-{task_id}-"))
    status, error = "succeeded", None
    try:
        logs, outputs = work / "logs", work / "outputs"
        logs.mkdir()
        outputs.mkdir()
        base_env = {**_base_env(), **cmd["env"],
                    "TASKQ_TASK_ID": task_id, "TASKQ_ATTEMPT": str(attempt),
                    "TASKQ_SOURCE_DIR": str(tree)}

        t_setup = 0.0
        for i, step in enumerate(cmd["setup"]):
            m = run_measured(step, cwd, base_env, limits["setup_timeout_seconds"],
                             logs / f"setup-{i}.stdout", logs / f"setup-{i}.stderr",
                             stop, cpus, limits["memory_mb"])
            t_setup += m.wall_seconds
            body["setup"].append({"index": i, "argv": step, **m.__dict__})
            if m.stopped == "fenced":
                raise Fenced(task_id)
            if m.stopped == "shutdown":
                raise InfraError("worker shutdown interrupted setup")
            if m.stopped == "cancel":
                status, error = "cancelled", f"cancelled during setup step {i}"
            elif m.timed_out:
                status, error = "timeout", f"setup step {i} exceeded {limits['setup_timeout_seconds']}s"
            elif m.exit_code != 0:
                status, error = "failed", f"setup step {i} exited {m.exit_code} (signal {m.signal})"
            if status != "succeeded":
                break
        body["setup_wall_seconds"] = t_setup

        if status == "succeeded":
            bench = spec.get("benchmark") or {"warmups": 0, "repetitions": 1}
            plan = [True] * bench["warmups"] + [False] * bench["repetitions"]
            for i, warm in enumerate(plan):
                run_out = outputs / (f"run-{i}" if len(plan) > 1 else "")
                run_out.mkdir(parents=True, exist_ok=True)
                env = {**base_env, "TASKQ_OUTPUT_DIR": str(run_out),
                       "TASKQ_REPETITION": str(i), "TASKQ_WARMUP": "1" if warm else "0"}
                m = run_measured(cmd["argv"], cwd, env, limits["timeout_seconds"],
                                 logs / "stdout.log", logs / "stderr.log",
                                 stop, cpus, limits["memory_mb"])
                metrics = None
                mfile = run_out / "metrics.json"
                if mfile.exists():
                    try:
                        metrics = json.loads(mfile.read_text())
                    except ValueError as e:
                        metrics = {"_taskq_error": f"metrics.json unparseable: {e}"}
                rec = {"index": i, "warmup": warm, "metrics": metrics, **m.__dict__}
                body["runs"].append(rec)
                if m.stopped == "fenced":
                    raise Fenced(task_id)
                if m.stopped == "shutdown":
                    raise InfraError("worker shutdown interrupted the run")
                if m.stopped == "cancel":
                    status, error = "cancelled", f"cancelled during run {i}"
                elif m.timed_out:
                    status, error = "timeout", f"run {i} exceeded {limits['timeout_seconds']}s"
                elif m.exit_code != 0:
                    status, error = "failed", f"run {i} exited {m.exit_code} (signal {m.signal})"
                if status != "succeeded":
                    break
                if spec.get("verify") and not warm:
                    rec["verification"] = verify_run(spec["verify"], run_out, cwd, env)
        body["summary"] = _summary(body["runs"])
        if spec.get("verify"):
            body["verification"] = rollup_verification(body["runs"])

        for p in logs.iterdir():
            if _cap_log(p, limits["max_log_bytes"]):
                p.rename(p.with_name(p.name + ".truncated"))
        dest = artifact_root / task_id / f"attempt-{attempt}" if artifact_root else None
        shutil.move(str(logs), str(outputs / "_taskq_logs"))
        body["artifacts"] = collect_artifacts(outputs, dest)
    finally:
        shutil.rmtree(work, ignore_errors=True)
    body.setdefault("setup_wall_seconds", 0.0)
    return {**body, "status": status, "error": error}


def _base_env() -> dict[str, str]:
    keep = ("PATH", "HOME", "LANG", "LC_ALL", "TZ", "TMPDIR", "USER",
            "CARGO_HOME", "RUSTUP_HOME", "CARGO_TARGET_DIR", "SAGE_ROOT",
            "PYTHONPATH", "VIRTUAL_ENV", "CUDA_HOME", "LD_LIBRARY_PATH")
    return {k: os.environ[k] for k in keep if k in os.environ}
