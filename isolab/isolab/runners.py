"""Execution backends: where the user's process runs and how it is timed.

``direct`` runs on the host inside the job cgroup. ``podman`` and ``docker``
hold a container open for the whole job (``isolab-launch --hold`` as PID 1)
and ``exec`` each build step and repeat into it, so container start-up is
never inside a measurement. Any container backend takes ``runsc`` (gVisor)
as its OCI runtime. The static launcher takes the clock and ``wait4`` inside
the box; the host wraps the same exec in ``perf stat`` over the job CPUs.
"""
from __future__ import annotations

import json
import os
import resource
import shlex
import shutil
import signal
import subprocess
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable

from .inventory import default_run, read_text
from .isolation import kill_cgroup

StopFn = Callable[[], str | None]
CONTAINER_WORK = "/isolab/work"
CONTAINER_OUT = "/isolab/out"
CONTAINER_SCRATCH = "/isolab/scratch"
CONTAINER_LAUNCH = "/isolab/bin/launch"


class BackendError(RuntimeError):
    """The backend itself failed (image missing, runtime broken): an infrastructure error."""


@dataclass
class JobDirs:
    root: Path
    work: Path
    out: Path
    logs: Path
    scratch: Path

    @classmethod
    def create(cls, root: Path) -> "JobDirs":
        d = cls(root, root / "work", root / "out", root / "logs", root / "scratch")
        for p in (d.root, d.work, d.out, d.logs, d.scratch):
            p.mkdir(parents=True, exist_ok=True)
        for p in (d.work, d.out, d.scratch):
            os.chmod(p, 0o1777)
        return d


@dataclass
class ExecOutcome:
    exit_code: int | None
    signal: int | None
    timed_out: bool
    stopped: str | None
    wall_host_s: float
    launcher: dict[str, Any] | None
    perf_text: str | None = None
    host_rusage: dict[str, Any] | None = None


@dataclass
class JobContext:
    job_id: str
    dirs: JobDirs
    cpus: list[int]
    mems: list[int]
    cgroup_dir: Path | None
    cgroup_parent: str | None
    image: str | None
    oci_runtime: str | None
    network: str
    memory_mb: int | None
    pids: int
    scratch_mb: int
    user: str | None
    extra_args: list[str]
    gpus: list[int]
    launcher: str | None
    aslr_off: bool = False
    userns: str | None = None
    env_base: dict[str, str] = field(default_factory=dict)


def _translate(dirs: JobDirs, container: bool, path: str) -> str:
    """A host path under the job dirs, as the process will see it."""
    if not container:
        return path
    for host, inside in ((dirs.work, CONTAINER_WORK), (dirs.out, CONTAINER_OUT), (dirs.scratch, CONTAINER_SCRATCH)):
        try:
            rel = Path(path).relative_to(host)
            return str(Path(inside) / rel) if str(rel) != "." else inside
        except ValueError:
            continue
    return path


class Backend:
    name = "abstract"
    container = False

    def __init__(self, run=default_run, which=shutil.which, privileged: bool = False,
                 cgroup_manager: str = "cgroupfs"):
        self.run, self.which, self.privileged = run, which, privileged
        self.cgroup_manager = cgroup_manager

    # ---- lifecycle
    def resolve_image(self, image: str, pull: str = "missing") -> dict[str, Any]:
        return {"image": None, "digest": None, "id": None}

    def start(self, ctx: JobContext) -> dict[str, Any]:
        return {}

    def stop(self, ctx: JobContext) -> None:
        pass

    def paths_for(self, ctx: JobContext, host_path: Path) -> str:
        return _translate(ctx.dirs, self.container, str(host_path))

    def exec(self, ctx: JobContext, argv: list[str], cwd: str, env: dict[str, str], timeout: float,
             stdout: Path, stderr: Path, stop: StopFn, stdin_text: str | None = None,
             record: Path | None = None, perf_prefix: list[str] | None = None) -> ExecOutcome:
        raise NotImplementedError

    def cgroup_stats(self, ctx: JobContext) -> dict[str, Any] | None:
        cg = ctx.cgroup_dir
        if cg is None or not cg.exists():
            return None
        out: dict[str, Any] = {}
        for line in (read_text(cg / "cpu.stat") or "").splitlines():
            k, _, v = line.partition(" ")
            if v.strip().lstrip("-").isdigit():
                out[k] = int(v)
        peak = read_text(cg / "memory.peak")
        out["memory_peak_bytes"] = int(peak) if peak and peak.isdigit() else None
        for line in (read_text(cg / "memory.events") or "").splitlines():
            k, _, v = line.partition(" ")
            if v.strip().isdigit():
                out[k] = int(v)
        pp = read_text(cg / "pids.peak")
        out["pids_peak"] = int(pp) if pp and pp.isdigit() else None
        return out

    # ---- shared process driving
    def _drive(self, popen_argv: list[str], cwd: Path | None, env: dict[str, str] | None, timeout: float,
               stdout: Path, stderr: Path, stop: StopFn, kill: Callable[[], None],
               preexec=None, stdin_text: str | None = None) -> tuple[int | None, int | None, bool, str | None, float, dict]:
        with open(stdout, "ab") as out, open(stderr, "ab") as err:
            stdin = subprocess.PIPE if stdin_text is not None else subprocess.DEVNULL
            t0 = time.perf_counter()
            try:
                proc = subprocess.Popen(popen_argv, cwd=cwd, env=env, stdout=out, stderr=err, stdin=stdin,
                                        start_new_session=True, preexec_fn=preexec)
            except OSError as e:
                err.write(f"isolab: cannot exec {popen_argv[0]!r}: {e}\n".encode())
                return 127, None, False, None, 0.0, {}
            if stdin_text is not None and proc.stdin:
                try:
                    proc.stdin.write(stdin_text.encode())
                    proc.stdin.close()
                except OSError:
                    pass
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
                    kill()
                    status, ru = _kill_group(proc.pid)
                    break
                time.sleep(0.02 if elapsed < 2 else 0.1)
            wall = time.perf_counter() - t0
        proc.returncode = 0
        code = os.waitstatus_to_exitcode(status)
        host_ru = {"user_s": ru.ru_utime, "sys_s": ru.ru_stime, "max_rss_kb": ru.ru_maxrss,
                   "involuntary_switches": ru.ru_nivcsw, "voluntary_switches": ru.ru_nvcsw}
        return (code if code >= 0 else None), (-code if code < 0 else None), timed_out, stopped, wall, host_ru


def _kill_group(pid: int):
    for sig, grace in ((signal.SIGTERM, 3.0), (signal.SIGKILL, None)):
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
                with contextlib_suppress(ProcessLookupError):
                    os.killpg(pid, signal.SIGKILL)
                return status, ru
            time.sleep(0.05)


class contextlib_suppress:
    def __init__(self, *exc):
        self.exc = exc

    def __enter__(self):
        return self

    def __exit__(self, et, ev, tb):
        return et is not None and issubclass(et, self.exc)


def _read_record(path: Path | None) -> dict[str, Any] | None:
    if path is None or not path.exists():
        return None
    try:
        return json.loads(path.read_text())
    except ValueError:
        return None


# -- direct ------------------------------------------------------------------

class DirectBackend(Backend):
    name = "direct"
    container = False

    def exec(self, ctx, argv, cwd, env, timeout, stdout, stderr, stop, stdin_text=None, record=None, perf_prefix=None):
        workdir = (ctx.dirs.work / cwd).resolve()
        full_env = {**ctx.env_base, **env}
        cmd: list[str] = []
        if ctx.mems and self.which("numactl") and ctx.cgroup_dir is None:
            cmd += ["numactl", f"--membind={','.join(map(str, ctx.mems))}"]
        if ctx.aslr_off and self.which("setarch"):
            cmd += ["setarch", os.uname().machine, "-R"]
        if ctx.launcher:
            cmd += [ctx.launcher]
            if record:
                cmd += ["--record", str(record)]
            cmd += ["--"]
        cmd += list(argv)
        if perf_prefix:
            cmd = perf_prefix + cmd
        cpus = ctx.cpus
        cg = ctx.cgroup_dir
        mem = ctx.memory_mb

        def preexec() -> None:
            if cg is not None:
                try:
                    (cg / "cgroup.procs").write_text(str(os.getpid()))
                except OSError:
                    pass
            if cpus and hasattr(os, "sched_setaffinity"):
                try:
                    os.sched_setaffinity(0, cpus)
                except OSError:
                    pass
            if mem and cg is None:
                lim = mem * 1024 * 1024
                try:
                    resource.setrlimit(resource.RLIMIT_AS, (lim, lim))
                except (ValueError, OSError):
                    pass

        def kill() -> None:
            if cg is not None:
                kill_cgroup(cg)

        code, sig, timed_out, stopped, wall, host_ru = self._drive(
            cmd, workdir, full_env, timeout, stdout, stderr, stop, kill, preexec, stdin_text)
        rec = _read_record(record)
        if rec is None and ctx.launcher is None:
            rec = {"wall_s": wall, **host_ru, "exit_code": code, "signal": sig or 0, "source": "host-wait4"}
        return ExecOutcome(code, sig, timed_out, stopped, wall, rec, host_rusage=host_ru)


# -- podman / docker ---------------------------------------------------------

class PodmanBackend(Backend):
    name = "podman"
    container = True
    exe = "podman"

    def __init__(self, *a, **kw):
        super().__init__(*a, **kw)
        self.path = self.which(self.exe)

    def available(self) -> bool:
        if not self.path:
            return False
        rc, _ = self.run([self.path, "info", "--format", "{{.Host.OS}}" if self.exe == "podman" else "{{.OSType}}"], 30)
        return rc == 0

    def _base(self) -> list[str]:
        base = [self.path]
        if self.exe == "podman" and self.privileged and self.cgroup_manager:
            base += ["--cgroup-manager", self.cgroup_manager]
        return base

    def resolve_image(self, image: str, pull: str = "missing") -> dict[str, Any]:
        fmt = "{{.Id}}|{{.Digest}}|{{json .RepoDigests}}"
        rc, text = self.run([*self._base(), "image", "inspect", "--format", fmt, image], 60)
        if rc != 0:
            if pull == "never":
                raise BackendError(f"image {image!r} not present and pull policy is never")
            rc2, text2 = self.run([*self._base(), "pull", "-q", image], 1800)
            if rc2 != 0:
                raise BackendError(f"cannot pull {image!r}: {text2.strip()[-400:]}")
            rc, text = self.run([*self._base(), "image", "inspect", "--format", fmt, image], 60)
            if rc != 0:
                raise BackendError(f"image {image!r} not inspectable after pull: {text.strip()[-300:]}")
        parts = text.strip().split("|", 2)
        digest = parts[1] if len(parts) > 1 and parts[1] not in ("", "<no value>") else None
        if (not digest or not digest.startswith("sha256:")) and len(parts) > 2:
            try:
                rd = json.loads(parts[2]) or []
                digest = rd[0].split("@", 1)[1] if rd and "@" in rd[0] else digest
            except (ValueError, IndexError):
                pass
        return {"image": image, "id": parts[0], "digest": digest}

    def container_name(self, ctx: JobContext) -> str:
        return f"isolab-{ctx.job_id.lower()}"

    def start(self, ctx: JobContext) -> dict[str, Any]:
        if not ctx.launcher:
            raise BackendError("the in-container launcher is missing; build it with `isolab launcher-build`")
        name = self.container_name(ctx)
        self.run([*self._base(), "rm", "-f", "-t", "1", name], 60)
        argv = [*self._base(), "run", "-d", "--name", name, "--init=false",
                "--cpuset-cpus", ",".join(map(str, ctx.cpus)),
                "--pids-limit", str(ctx.pids), "--network", ctx.network,
                "--cap-drop", "ALL", "--security-opt", "no-new-privileges", "--read-only",
                "--tmpfs", "/tmp:rw,exec,nosuid,size=2g",
                "-v", f"{ctx.launcher}:{CONTAINER_LAUNCH}:ro",
                "-v", f"{ctx.dirs.work}:{CONTAINER_WORK}:rw",
                "-v", f"{ctx.dirs.out}:{CONTAINER_OUT}:rw",
                "-w", CONTAINER_WORK,
                "-e", f"ISOLAB_WORK={CONTAINER_WORK}", "-e", f"ISOLAB_SCRATCH={CONTAINER_SCRATCH}"]
        if ctx.mems:
            argv += ["--cpuset-mems", ",".join(map(str, ctx.mems))]
        if ctx.memory_mb:
            argv += ["--memory", f"{ctx.memory_mb}m", "--memory-swap", f"{ctx.memory_mb}m"]
        if ctx.scratch_mb:
            argv += ["--tmpfs", f"{CONTAINER_SCRATCH}:rw,exec,nosuid,size={ctx.scratch_mb}m"]
        else:
            argv += ["-v", f"{ctx.dirs.scratch}:{CONTAINER_SCRATCH}:rw"]
        if ctx.cgroup_parent and self.privileged and self._cgroup_parent_is_path():
            argv += ["--cgroup-parent", ctx.cgroup_parent]
        if ctx.oci_runtime and ctx.oci_runtime != "auto":
            argv += ["--runtime", ctx.oci_runtime]
        if ctx.user:
            argv += ["--user", ctx.user]
        if ctx.userns:
            argv += ["--userns", ctx.userns]
        for g in ctx.gpus:
            argv += self._gpu_args(g)
        argv += list(ctx.extra_args)
        argv += [ctx.image, CONTAINER_LAUNCH, "--hold"]
        rc, text = self.run(argv, 300)
        if rc != 0:
            raise BackendError(f"{self.exe} run failed: {text.strip()[-600:]}\ncommand: {shlex.join(argv)}")
        cid = text.strip().splitlines()[-1] if text.strip() else name
        rc, info = self.run([*self._base(), "inspect", "--format",
                             "{{.Id}}|{{.ImageName}}|{{.Image}}|{{.HostConfig.Runtime}}|{{.HostConfig.CgroupParent}}", name], 60)
        fields = (info.strip().split("|") + [None] * 5)[:5] if rc == 0 else [cid, ctx.image, None, None, None]
        return {"container": name, "id": fields[0], "image_name": fields[1], "image_id": fields[2],
                "runtime": fields[3], "cgroup_parent": fields[4], "run_argv": argv}

    def _cgroup_parent_is_path(self) -> bool:
        return True

    def _gpu_args(self, index: int) -> list[str]:
        return ["--device", f"nvidia.com/gpu={index}"]

    def stop(self, ctx: JobContext) -> None:
        self.run([*self._base(), "rm", "-f", "-t", "2", self.container_name(ctx)], 120)

    def exec(self, ctx, argv, cwd, env, timeout, stdout, stderr, stop, stdin_text=None, record=None, perf_prefix=None):
        name = self.container_name(ctx)
        inside_cwd = str(Path(CONTAINER_WORK) / cwd) if cwd not in (".", "") else CONTAINER_WORK
        cmd = [*self._base(), "exec", "-w", inside_cwd]
        if stdin_text is not None:
            cmd += ["-i"]
        for k, v in env.items():
            cmd += ["-e", f"{k}={v}"]
        cmd += [name, CONTAINER_LAUNCH]
        if record is not None:
            cmd += ["--record", self.paths_for(ctx, record)]
        cmd += ["--", *argv]
        if perf_prefix:
            cmd = perf_prefix + cmd

        def kill() -> None:
            if ctx.cgroup_dir is not None:
                kill_cgroup(ctx.cgroup_dir)
            self.run([*self._base(), "kill", "--signal", "KILL", name], 30)

        code, sig, timed_out, stopped, wall, host_ru = self._drive(
            cmd, None, None, timeout, stdout, stderr, stop, kill, None, stdin_text)
        rec = _read_record(record)
        return ExecOutcome(code, sig, timed_out, stopped, wall, rec, host_rusage=host_ru)

    def version(self) -> str | None:
        rc, text = self.run([self.path, "--version"], 10)
        return text.strip().splitlines()[0] if rc == 0 and text.strip() else None


class DockerBackend(PodmanBackend):
    name = "docker"
    exe = "docker"

    def __init__(self, *a, **kw):
        super().__init__(*a, **kw)
        self._driver = None

    def _base(self) -> list[str]:
        return [self.path]

    def _cgroup_parent_is_path(self) -> bool:
        if self._driver is None:
            rc, text = self.run([self.path, "info", "--format", "{{.CgroupDriver}}"], 30)
            self._driver = text.strip() if rc == 0 else "systemd"
        return self._driver == "cgroupfs"

    def _gpu_args(self, index: int) -> list[str]:
        return ["--gpus", f"device={index}"]


def oci_runtimes(which=shutil.which) -> list[str]:
    return [r for r in ("crun", "runc", "runsc") if which(r)]


def make_backends(which=shutil.which, run=default_run, privileged: bool = False,
                  cgroup_manager: str = "cgroupfs") -> dict[str, Backend]:
    out: dict[str, Backend] = {}
    for cls in (PodmanBackend, DockerBackend):
        b = cls(run=run, which=which, privileged=privileged, cgroup_manager=cgroup_manager)
        if b.available():
            out[b.name] = b
    out["direct"] = DirectBackend(run=run, which=which, privileged=privileged)
    return out
