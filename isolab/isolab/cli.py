"""`isolab` command line: hub, worker, doctor, images, and the client commands."""
from __future__ import annotations

import argparse
import asyncio
import json
import logging
import os
import platform
import shutil
import subprocess
import sys
from importlib import resources
from pathlib import Path
from typing import Any

from . import __version__, protocol
from .inventory import Paths, format_cpulist, parse_cpulist


def _print(doc: Any) -> None:
    print(json.dumps(doc, indent=2, sort_keys=True, default=str))


def _labels(pairs: list[str] | None) -> dict[str, str]:
    out = {}
    for p in pairs or []:
        k, _, v = p.partition("=")
        out[k] = v
    return out


def _default_state_dir() -> Path:
    if hasattr(os, "geteuid") and os.geteuid() == 0 and platform.system() == "Linux":
        return Path("/var/lib/isolab")
    return Path.home() / ".isolab"


async def _with_fabric(fn):
    from .fabric import Fabric
    f = await Fabric(name="isolab-cli").connect()
    try:
        return await fn(f)
    finally:
        await f.close()


def _run(fn) -> Any:
    return asyncio.run(_with_fabric(fn))


# -- lab host commands -------------------------------------------------------

def cmd_hub(args) -> None:
    from .hub import run_hub
    sys.exit(run_hub(args))


def cmd_nats_download(args) -> None:
    from .hub import nats_download
    print(nats_download(Path(args.dest)))


def cmd_launcher_build(args) -> None:
    from .worker import build_launcher
    launcher, info = build_launcher(Path(args.dest), static=not args.dynamic)
    _print({"launcher": str(launcher) if launcher else None, **info})
    if launcher is None:
        sys.exit(1)


def cmd_worker(args) -> None:
    from .fabric import Fabric
    from .worker import Worker, default_worker_id, install_signal_handlers
    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO,
                        format="%(asctime)s %(levelname)s %(name)s %(message)s")
    cpus = parse_cpulist(args.cpus) if args.cpus else []

    async def go() -> None:
        fabric = await Fabric(name=f"isolab-worker-{args.id or default_worker_id()}", lease_s=args.lease_seconds).connect()
        w = Worker(fabric, worker_id=args.id or default_worker_id(), pools=args.pool or ["default"],
                   labels=_labels(args.label), lab_cpus=cpus, state_dir=Path(args.state_dir).expanduser(),
                   backend=args.backend, default_image=args.default_image, oci_runtime=args.oci_runtime,
                   pull_policy=args.pull, userns=args.userns, keep_job_dirs=args.keep_job_dirs,
                   lock_path=args.lock, cgroup_root=args.cgroup_root)
        install_signal_handlers(w)
        try:
            await w.start()
            await w.run_forever(max_jobs=args.max_jobs)
        finally:
            await fabric.close()
    asyncio.run(go())


def cmd_inventory(args) -> None:
    from .inventory import default_run, host_inventory, inventory_summary
    from .isolation import max_tier, probe_capabilities
    from . import perfstat
    inv = host_inventory(lab_cpus=parse_cpulist(args.cpus) if args.cpus else None)
    perf = perfstat.probe(perfstat.find_perf(), default_run)
    caps = probe_capabilities(perf=perf, numactl=bool(shutil.which("numactl")), container_backend=bool(shutil.which("podman") or shutil.which("docker")))
    if args.json:
        _print({"inventory": inv, "capabilities": caps, "max_tier": max_tier(caps)})
    else:
        print(inventory_summary(inv))
        print(f"max isolation tier here: {max_tier(caps)}; perf hw counters: {perf.get('hw_events')}")
        t = inv["topology"]
        for nid, nd in t["nodes"].items():
            print(f"  node {nid}: cpus {format_cpulist(nd['cpus'])}, {(nd.get('mem_total_kb') or 0) // 1024} MiB")


def cmd_doctor(args) -> None:
    from .doctor import apply, as_json, diagnose, render
    from .inventory import default_run, host_inventory
    from .isolation import probe_capabilities
    from . import perfstat
    cpus = parse_cpulist(args.cpus) if args.cpus else []
    inv = host_inventory(lab_cpus=cpus or None)
    if not cpus:
        cpus = inv["topology"]["online"][1:]
    perf = perfstat.probe(perfstat.find_perf(), default_run, cpus[:1])
    caps = probe_capabilities(perf=perf, numactl=bool(shutil.which("numactl")), container_backend=bool(shutil.which("podman") or shutil.which("docker")))
    launcher = Path(args.state_dir).expanduser() / "bin" / "isolab-launch"
    items = diagnose(inv, caps, cpus, launcher=str(launcher))
    if args.apply:
        done = apply(items)
        print("applied: " + (", ".join(done) if done else "nothing to apply"))
        inv = host_inventory(lab_cpus=cpus)
        caps = probe_capabilities(perf=perf, numactl=bool(shutil.which("numactl")), container_backend=True)
        items = diagnose(inv, caps, cpus, launcher=str(launcher))
    if args.json:
        _print(as_json(items))
    else:
        print(render(items))


def _images_dir() -> Path:
    return Path(str(resources.files("isolab").joinpath("images")))


def cmd_images(args) -> None:
    tool = args.tool or ("podman" if shutil.which("podman") else "docker")
    if args.action == "list":
        subprocess.run([tool, "images", "--filter", "reference=*isolab*"], check=False)
        return
    for name in args.names:
        d = _images_dir() / name
        if not (d / "Containerfile").exists():
            sys.exit(f"no image recipe {name!r}; have {[p.name for p in _images_dir().iterdir() if p.is_dir()]}")
        tag = f"localhost/isolab-{name}:latest"
        argv = [tool, "build", "-t", tag, "-f", str(d / "Containerfile")]
        for a in args.build_arg or []:
            argv += ["--build-arg", a]
        argv.append(str(d))
        print("+", " ".join(argv), file=sys.stderr)
        subprocess.run(argv, check=True)
        print(tag)


# -- client commands ---------------------------------------------------------

def _spec_from_run_args(args) -> dict[str, Any]:
    inputs = []
    for item in args.input or []:
        src, _, dest = item.partition(":")
        inputs.append({"path": dest or Path(src).name, "local_file": src})
    spec: dict[str, Any] = {
        "schema": protocol.SPEC_SCHEMA_ID, "pool": args.pool,
        "runtime": {"backend": args.backend, "oci_runtime": args.oci_runtime, "image": args.image, "network": args.network},
        "inputs": inputs, "build": [{"argv": b.split()} for b in args.build or []],
        "command": {"argv": args.argv, "cwd": args.cwd, "env": _labels(args.env)},
        "measure": {"warmups": args.warmups, "repeats": args.repeats, "cooldown_s": args.cooldown},
        "resources": {"cpus": args.cpus, "memory_mb": args.memory_mb, "numa_node": args.numa, "smt": args.smt, "gpus": args.gpus},
        "fidelity": {"policy": args.policy, **({"require": args.require} if args.require else {})},
        "placement": {"worker": args.worker} if args.worker else {},
        "limits": {"timeout_s": args.timeout}, "labels": _labels(args.label),
    }
    if args.name:
        spec["name"] = args.name
    if args.verify:
        spec["verify"] = {"builtin": "taskq-certificate"}
    if args.idempotency_key:
        spec["idempotency_key"] = args.idempotency_key
    if args.numa.isdigit():
        spec["resources"]["numa_node"] = int(args.numa)
    return spec


def cmd_submit(args) -> None:
    from .blobs import resolve_local_inputs
    from .planner import match_worker
    if getattr(args, "spec", None):
        spec = json.loads(Path(args.spec).read_text() if args.spec != "-" else sys.stdin.read())
    else:
        if not args.argv:
            sys.exit("give a command after --")
        spec = _spec_from_run_args(args)

    async def go(f):
        norm = protocol.normalize_spec(spec, allow_local_files=True)
        pending: list[tuple[str, Path]] = []
        resolve_local_inputs(norm, lambda d, p: pending.append((d, p)), lambda d: False)
        for digest, path in pending:
            if not await f.has_blob(digest):
                await f.put_blob(digest, path)
        norm = protocol.normalize_spec(norm)
        workers = await f.workers()
        eligible = [w["id"] for w in workers if not match_worker(norm, w)]
        if not eligible and not args.force:
            why = {w["id"]: match_worker(norm, w) for w in workers}
            sys.exit("no online worker can take this job: " + (json.dumps(why) if why else "none online") +
                     "\n(--force enqueues it anyway)")
        out = await f.submit(norm)
        out["eligible_workers"] = eligible
        if args.wait:
            await f.wait(out["job_id"], args.timeout + 600)
            out["result"] = await f.get_result(out["job_id"])
        return out
    _print(_run(go))


def cmd_fetch(args) -> None:
    from .blobs import sha256_file

    async def go(f):
        res = await f.get_result(args.job_id)
        if res is None:
            sys.exit(f"{args.job_id} has no result yet")
        dest = Path(args.dest) / args.job_id
        got = []
        for a in res.get("artifacts") or []:
            if not a.get("stored") or (args.path and a["path"] != args.path):
                continue
            target = await f.get_artifact(args.job_id, a["path"], dest / a["path"])
            got.append({"path": str(target), "sha256_ok": sha256_file(target) == a["sha256"]})
        return got
    _print(_run(go))


def cmd_logs(args) -> None:
    async def go(f):
        try:
            for name in ("stdout", "stderr"):
                data = await f.get_artifact_bytes(args.job_id, f"_isolab_logs/{name}.log", args.tail)
                print(f"===== {name} =====")
                print(data.decode(errors="replace"))
        except Exception:  # noqa: BLE001
            prog = await f.get_progress(args.job_id)
            print((prog or {}).get("log_tail") or "(no logs yet)")
    _run(go)


def cmd_calibrate(args) -> None:
    from .mcp_server import isolab_calibrate
    async def go(f):
        import isolab.mcp_server as m
        m._fabric = f
        return await isolab_calibrate(args.worker, args.repeats, args.cpus, args.policy, args.iterations, args.backend, args.image)
    _print(_run(go))


def _add_run_args(p: argparse.ArgumentParser) -> None:
    p.add_argument("--name")
    p.add_argument("--pool", default="default")
    p.add_argument("--image")
    p.add_argument("--backend", default="auto", choices=["auto", "podman", "docker", "direct"])
    p.add_argument("--oci-runtime", default="auto", choices=["auto", "crun", "runc", "runsc"])
    p.add_argument("--network", default="none", choices=["none", "host", "bridge"])
    p.add_argument("--input", action="append", metavar="SRC[:DEST]", help="local file to upload; DEST is its path in the workdir")
    p.add_argument("--build", action="append", metavar="CMD", help="untimed build step, whitespace-split")
    p.add_argument("--cwd", default=".")
    p.add_argument("--env", action="append", metavar="K=V")
    p.add_argument("--cpus", type=int, default=1)
    p.add_argument("--memory-mb", type=int)
    p.add_argument("--numa", default="single", help="single | any | <node id>")
    p.add_argument("--smt", default="isolate", choices=["isolate", "allow", "off"])
    p.add_argument("--gpus", type=int, default=0)
    p.add_argument("--repeats", type=int, default=1)
    p.add_argument("--warmups", type=int, default=0)
    p.add_argument("--cooldown", type=float, default=0.0)
    p.add_argument("--policy", default="standard", choices=["strict", "standard", "best_effort"])
    p.add_argument("--require", action="append")
    p.add_argument("--worker")
    p.add_argument("--timeout", type=float, default=3600)
    p.add_argument("--label", action="append", metavar="K=V")
    p.add_argument("--verify", action="store_true", help="independently verify each repeat's certificate.json")
    p.add_argument("--idempotency-key")
    p.add_argument("--force", action="store_true", help="enqueue even if no online worker matches")
    p.add_argument("--wait", action="store_true", help="wait for the result and print it")
    p.add_argument("argv", nargs=argparse.REMAINDER)


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(prog="isolab", description="experiment execution environment with enforced fidelity")
    ap.add_argument("--version", action="version", version=__version__)
    sub = ap.add_subparsers(dest="cmd", required=True)

    h = sub.add_parser("hub", help="run the NATS JetStream hub")
    h.add_argument("--listen", default=os.environ.get("ISOLAB_HUB_LISTEN", "0.0.0.0:4222"))
    h.add_argument("--store", default=os.environ.get("ISOLAB_HUB_STORE", str(_default_state_dir() / "hub")))
    h.add_argument("--token", default=None)
    h.add_argument("--name")
    h.add_argument("--http", help="monitoring endpoint, e.g. 127.0.0.1:8222")
    h.add_argument("--max-memory", default="1GB")
    h.add_argument("--max-file", default="50GB")
    h.add_argument("--cluster-name")
    h.add_argument("--cluster-listen")
    h.add_argument("--routes", help="comma-separated nats://host:6222 of the other hubs")
    h.add_argument("--leaf-remote", help="nats://TOKEN@hub:7422 to join another lab as a leaf")
    h.add_argument("--leaf-listen", help="accept leaf nodes, e.g. 0.0.0.0:7422")
    h.add_argument("--nats-server", help="path to the nats-server binary")
    h.add_argument("--download", action="store_true", help="download nats-server if not found")
    h.add_argument("--print-config", action="store_true")
    h.set_defaults(fn=cmd_hub)

    d = sub.add_parser("nats-download", help="download the nats-server binary for this platform")
    d.add_argument("--dest", default=str(_default_state_dir() / "bin"))
    d.set_defaults(fn=cmd_nats_download)

    lb = sub.add_parser("launcher-build", help="compile the in-container launcher and calibration kernel")
    lb.add_argument("--dest", default=str(_default_state_dir() / "bin"))
    lb.add_argument("--dynamic", action="store_true")
    lb.set_defaults(fn=cmd_launcher_build)

    w = sub.add_parser("worker", help="run a worker daemon")
    w.add_argument("--id")
    w.add_argument("--pool", action="append", help="in priority order; default: default")
    w.add_argument("--label", action="append", metavar="K=V")
    w.add_argument("--cpus", help="lab cpus this worker may hand to jobs, e.g. 4-15")
    w.add_argument("--backend", default="auto", choices=["auto", "podman", "docker", "direct"])
    w.add_argument("--default-image", default=os.environ.get("ISOLAB_DEFAULT_IMAGE"))
    w.add_argument("--oci-runtime", help="default OCI runtime for containers (crun, runc, runsc)")
    w.add_argument("--state-dir", default=os.environ.get("ISOLAB_STATE_DIR", str(_default_state_dir())))
    w.add_argument("--pull", default="missing", choices=["missing", "never", "always"])
    w.add_argument("--userns", help="podman --userns value, e.g. auto")
    w.add_argument("--keep-job-dirs", action="store_true")
    w.add_argument("--lock", default=os.environ.get("ISOLAB_LOCK", "/run/isolab.lock" if os.access("/run", os.W_OK) else "/tmp/isolab.lock"))
    w.add_argument("--cgroup-root", default="/sys/fs/cgroup")
    w.add_argument("--lease-seconds", type=float, default=90)
    w.add_argument("--max-jobs", type=int)
    w.add_argument("-v", "--verbose", action="store_true")
    w.set_defaults(fn=cmd_worker)

    inv = sub.add_parser("inventory", help="what this host is")
    inv.add_argument("--cpus")
    inv.add_argument("--json", action="store_true")
    inv.set_defaults(fn=cmd_inventory)

    doc = sub.add_parser("doctor", help="readiness of this host for strict runs")
    doc.add_argument("--cpus", help="the lab cpus (as given to the worker)")
    doc.add_argument("--apply", action="store_true", help="apply the runtime-tunable fixes")
    doc.add_argument("--state-dir", default=os.environ.get("ISOLAB_STATE_DIR", str(_default_state_dir())))
    doc.add_argument("--json", action="store_true")
    doc.set_defaults(fn=cmd_doctor)

    im = sub.add_parser("images", help="build or list the lab images on this host")
    im.add_argument("action", choices=["build", "list"])
    im.add_argument("names", nargs="*", help="base, sage, cuda")
    im.add_argument("--tool", choices=["podman", "docker"])
    im.add_argument("--build-arg", action="append")
    im.set_defaults(fn=cmd_images)

    r = sub.add_parser("run", help="submit one command (the common case)")
    _add_run_args(r)
    r.set_defaults(fn=cmd_submit)

    s = sub.add_parser("submit", help="submit a full isolab.job/v1 JSON document")
    s.add_argument("spec", help="file, or - for stdin")
    s.add_argument("--force", action="store_true")
    s.add_argument("--wait", action="store_true")
    s.add_argument("--timeout", type=float, default=3600)
    s.set_defaults(fn=cmd_submit, argv=None)

    for name, fn in (("status", lambda a: _print(_run(lambda f: f.get_job(a.job_id)))),
                     ("result", lambda a: _print(_run(lambda f: f.get_result(a.job_id)))),
                     ("cancel", lambda a: print(_run(lambda f: f.cancel(a.job_id))))):
        p = sub.add_parser(name)
        p.add_argument("job_id")
        p.set_defaults(fn=fn)

    p = sub.add_parser("wait")
    p.add_argument("job_id")
    p.add_argument("--timeout", type=float, default=3600)
    p.set_defaults(fn=lambda a: _print(_run(lambda f: f.wait(a.job_id, a.timeout))))

    p = sub.add_parser("logs")
    p.add_argument("job_id")
    p.add_argument("--tail", type=int, default=16384)
    p.set_defaults(fn=cmd_logs)

    p = sub.add_parser("fetch", help="download artifacts")
    p.add_argument("job_id")
    p.add_argument("path", nargs="?")
    p.add_argument("--dest", default=".")
    p.set_defaults(fn=cmd_fetch)

    p = sub.add_parser("jobs")
    p.add_argument("--state")
    p.add_argument("--pool")
    p.add_argument("--worker")
    p.add_argument("--label", action="append", metavar="K=V")
    p.add_argument("--limit", type=int, default=50)
    p.set_defaults(fn=lambda a: _print(_run(lambda f: f.list_jobs(a.limit, a.state, a.pool, a.worker, _labels(a.label)))))

    p = sub.add_parser("workers")
    p.add_argument("--full", action="store_true")
    p.add_argument("--stale", action="store_true")
    p.set_defaults(fn=lambda a: _print(_run(lambda f: f.workers(a.stale))) if a.full else _print(_run(
        lambda f: _summaries(f, a.stale))))

    p = sub.add_parser("calibrate", help="measure a worker's noise floor")
    p.add_argument("--worker")
    p.add_argument("--repeats", type=int, default=10)
    p.add_argument("--cpus", type=int, default=1)
    p.add_argument("--policy", default="standard")
    p.add_argument("--iterations", type=int, default=2_000_000_000)
    p.add_argument("--backend", default="auto")
    p.add_argument("--image")
    p.set_defaults(fn=cmd_calibrate)

    m = sub.add_parser("mcp", help="serve the MCP tools over stdio")
    m.set_defaults(fn=lambda a: __import__("isolab.mcp_server", fromlist=["main"]).main())

    sub.add_parser("schema").set_defaults(fn=lambda a: _print({"job_spec": protocol.SPEC_SCHEMA, "result": protocol.RESULT_SCHEMA}))

    args = ap.parse_args(argv)
    if getattr(args, "argv", None) and args.argv[:1] == ["--"]:
        args.argv = args.argv[1:]
    args.fn(args)


async def _summaries(f, stale: bool):
    from .mcp_server import _worker_line
    return [_worker_line(w) for w in await f.workers(stale)]


if __name__ == "__main__":
    main()
