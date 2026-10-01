"""`taskq` command line: submit and inspect tasks, run a worker or the MCP server."""
from __future__ import annotations

import argparse
import json
import logging
import os
import subprocess
import sys
from pathlib import Path

from . import protocol
from .config import load_repos, store_from_env


def _print(doc) -> None:
    print(json.dumps(doc, indent=2, sort_keys=True, default=str))


def _labels(pairs: list[str]) -> dict[str, str]:
    out = {}
    for p in pairs or []:
        k, _, v = p.partition("=")
        out[k] = v
    return out


def _commit(args) -> str:
    if args.commit:
        return args.commit
    here = args.checkout or "."
    sha = subprocess.run(["git", "rev-parse", "HEAD"], cwd=here, capture_output=True,
                         text=True, check=True).stdout.strip()
    dirty = subprocess.run(["git", "status", "--porcelain", "--untracked-files=no"],
                           cwd=here, capture_output=True, text=True).stdout.strip()
    if dirty:
        sys.exit(f"{here}: tracked changes are uncommitted; the worker would not run them. "
                 "Commit and push, or pass --commit.")
    return sha


def cmd_submit(args) -> None:
    store = store_from_env()
    if args.spec:
        spec = json.loads(Path(args.spec).read_text() if args.spec != "-" else sys.stdin.read())
    else:
        if not args.argv:
            sys.exit("give --spec FILE or a command after --")
        spec = {"schema": protocol.SPEC_SCHEMA_ID, "queue": args.queue,
                "kind": "benchmark" if args.repetitions else "command",
                "source": {"repo": args.repo, "commit": _commit(args)},
                "command": {"argv": args.argv, "cwd": args.cwd,
                            "setup": [s.split() for s in args.setup or []]},
                "limits": {"timeout_seconds": args.timeout},
                "labels": _labels(args.label)}
        if args.repetitions:
            spec["benchmark"] = {"repetitions": args.repetitions, "warmups": args.warmups}
        if args.idempotency_key:
            spec["idempotency_key"] = args.idempotency_key
        if args.verify:
            spec["verify"] = {"builtin": "certificate"}
    _print(store.submit(spec))


def cmd_worker(args) -> None:
    from .execute import RepoCache
    from .store import default_worker_id
    from .worker import Worker, install_signal_handlers
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(name)s %(message)s")
    cpus = None
    if args.cpus:
        cpus = [int(c) for c in args.cpus.split(",")]
    w = Worker(store_from_env(), args.queue or ["cpu"],
               RepoCache(args.cache_dir, load_repos(args.repos), args.max_trees),
               worker_id=args.id or os.environ.get("POD_NAME") or default_worker_id(),
               labels=_labels(args.label),
               artifact_dir=Path(args.artifact_dir) if args.artifact_dir else None,
               result_dir=Path(args.result_dir) if args.result_dir else None,
               lease_seconds=args.lease_seconds, cpus=cpus)
    install_signal_handlers(w)
    w.run_forever(max_tasks=args.max_tasks)


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(prog="taskq")
    sub = ap.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("submit", help="enqueue a task")
    s.add_argument("--spec", help="a taskq.task-spec/v1 JSON file, or - for stdin")
    s.add_argument("--queue", default="cpu")
    s.add_argument("--repo", default="crypto")
    s.add_argument("--commit", help="full sha; default: HEAD of --checkout (must be clean)")
    s.add_argument("--checkout", help="local clone used to resolve HEAD")
    s.add_argument("--cwd", default=".")
    s.add_argument("--setup", action="append", help="build step, whitespace-split")
    s.add_argument("--timeout", type=float, default=3600)
    s.add_argument("--repetitions", type=int, help="make it a benchmark")
    s.add_argument("--warmups", type=int, default=1)
    s.add_argument("--label", action="append", metavar="K=V")
    s.add_argument("--idempotency-key")
    s.add_argument("--verify", action="store_true",
                   help="independently verify each run's certificate.json")
    s.add_argument("argv", nargs=argparse.REMAINDER)
    s.set_defaults(fn=cmd_submit)

    for name, fn in (("status", lambda a: _print(store_from_env().get_task(a.task_id))),
                     ("result", lambda a: _print(store_from_env().get_result(a.task_id))),
                     ("cancel", lambda a: print(store_from_env().cancel(a.task_id)))):
        p = sub.add_parser(name)
        p.add_argument("task_id")
        p.set_defaults(fn=fn)

    p = sub.add_parser("wait")
    p.add_argument("task_id")
    p.add_argument("--timeout", type=float, default=3600)
    p.set_defaults(fn=lambda a: _print(store_from_env().wait(a.task_id, a.timeout)))

    p = sub.add_parser("list")
    p.add_argument("--state")
    p.add_argument("--queue")
    p.add_argument("--label", action="append", metavar="K=V")
    p.add_argument("--limit", type=int, default=50)
    p.set_defaults(fn=lambda a: _print(store_from_env().list_tasks(
        a.limit, a.state, a.queue, _labels(a.label))))

    sub.add_parser("stats").set_defaults(fn=lambda a: _print(store_from_env().queue_stats()))
    sub.add_parser("workers").set_defaults(fn=lambda a: _print(store_from_env().workers()))

    w = sub.add_parser("worker", help="run a worker daemon")
    w.add_argument("--queue", action="append", help="in priority order; default cpu")
    w.add_argument("--id")
    w.add_argument("--label", action="append", metavar="K=V")
    w.add_argument("--repos", help="repo allowlist file (default: TASKQ_REPOS or built-in)")
    w.add_argument("--cache-dir", default=os.environ.get("TASKQ_CACHE_DIR", "/var/cache/taskq"))
    w.add_argument("--artifact-dir", default=os.environ.get("TASKQ_ARTIFACT_DIR"))
    w.add_argument("--result-dir", default=os.environ.get("TASKQ_RESULT_DIR"))
    w.add_argument("--max-trees", type=int, default=8, help="cached checkouts per repo")
    w.add_argument("--lease-seconds", type=float, default=60)
    w.add_argument("--cpus", help="comma-separated CPU ids this worker may pin tasks to")
    w.add_argument("--max-tasks", type=int)
    w.set_defaults(fn=cmd_worker)

    m = sub.add_parser("mcp", help="serve the MCP tools over stdio")
    m.set_defaults(fn=lambda a: __import__("taskq.mcp_server", fromlist=["main"]).main())

    args = ap.parse_args(argv)
    if getattr(args, "argv", None) and args.argv[:1] == ["--"]:
        args.argv = args.argv[1:]
    args.fn(args)


if __name__ == "__main__":
    main()
