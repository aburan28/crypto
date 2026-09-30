"""MCP server: lets an agent submit measurements and read results.

Tools are thin wrappers over `Store`; the server never runs anything itself.
Configure it in `.mcp.json` with TASKQ_REDIS_URL in its environment.
Supports mcp 1.x (FastMCP) and 2.x (MCPServer).
"""
from __future__ import annotations

from typing import Any

try:  # mcp >= 2
    from mcp.server.mcpserver import MCPServer as FastMCP
except ImportError:  # mcp 1.x
    from mcp.server.fastmcp import FastMCP

from . import protocol
from .config import store_from_env

mcp = FastMCP("taskq", instructions=(
    "A persistent, Redis-backed queue of CPU/GPU jobs run by worker pods. "
    "Submit a task pinned to a pushed commit, then poll get_task or "
    "wait_for_task and read get_result. Status timeout/cancelled/infra_error "
    "(outcome_class not_completed) is never evidence about what was measured. "
    "Call get_protocol for the full spec and result schemas."))

_store = None


def store():
    global _store
    if _store is None:
        _store = store_from_env()
    return _store


@mcp.tool()
def get_protocol() -> dict[str, Any]:
    """The taskq.task-spec/v1 and taskq.task-result/v1 JSON schemas."""
    return {"task_spec": protocol.SPEC_SCHEMA, "task_result": protocol.RESULT_SCHEMA}


@mcp.tool()
def submit_task(spec: dict[str, Any]) -> dict[str, Any]:
    """Enqueue a full taskq.task-spec/v1 document. Returns task_id and spec_sha256."""
    return store().submit(spec)


@mcp.tool()
def submit_command(repo: str, commit: str, argv: list[str], queue: str = "cpu",
                   cwd: str = ".", setup: list[list[str]] | None = None,
                   env: dict[str, str] | None = None,
                   repetitions: int | None = None, warmups: int = 1,
                   timeout_seconds: float = 3600, cpus: int | None = None,
                   labels: dict[str, str] | None = None,
                   idempotency_key: str | None = None) -> dict[str, Any]:
    """Enqueue one command at a pushed commit (full 40-hex sha) of an allowlisted repo
    (crypto, crypto-autoresearcher, cryptanalysis). Give `repetitions` to make it a
    benchmark (warmups + timed repetitions). The command may write metrics.json and
    any other outputs into $TASKQ_OUTPUT_DIR; they are hashed and returned."""
    spec: dict[str, Any] = {
        "schema": protocol.SPEC_SCHEMA_ID, "queue": queue,
        "kind": "benchmark" if repetitions else "command",
        "source": {"repo": repo, "commit": commit},
        "command": {"argv": argv, "cwd": cwd, "setup": setup or [], "env": env or {}},
        "limits": {"timeout_seconds": timeout_seconds},
        "placement": {"cpus": cpus},
        "labels": labels or {},
    }
    if repetitions:
        spec["benchmark"] = {"repetitions": repetitions, "warmups": warmups}
    if idempotency_key:
        spec["idempotency_key"] = idempotency_key
    return store().submit(spec)


@mcp.tool()
def get_task(task_id: str) -> dict[str, Any] | None:
    """State, attempts, fence, infrastructure failures and the normalised spec."""
    return store().get_task(task_id)


@mcp.tool()
def get_result(task_id: str, include_runs: bool = True) -> dict[str, Any] | None:
    """The write-once result, or null if the task has not finished."""
    res = store().get_result(task_id)
    if res and not include_runs:
        res = {k: v for k, v in res.items() if k not in ("runs", "setup")}
    return res


@mcp.tool()
def wait_for_task(task_id: str, timeout_seconds: float = 60) -> dict[str, Any] | None:
    """Block up to timeout_seconds (max 300) for a terminal state; returns the task."""
    return store().wait(task_id, min(max(timeout_seconds, 0), 300))


@mcp.tool()
def cancel_task(task_id: str) -> str:
    """Cancel a queued task outright or ask the running worker to stop it."""
    return store().cancel(task_id)


@mcp.tool()
def list_tasks(state: str | None = None, queue: str | None = None,
               labels: dict[str, str] | None = None, limit: int = 50) -> list[dict[str, Any]]:
    """Newest first, filtered by state, queue and exact label matches."""
    return store().list_tasks(min(limit, 500), state, queue, labels)


@mcp.tool()
def queue_stats() -> dict[str, Any]:
    """Per-queue stream length, running (pending) count, unread lag and consumers."""
    return store().queue_stats()


@mcp.tool()
def list_workers() -> list[dict[str, Any]]:
    """Live workers: queues, labels and the task each is running."""
    return store().workers()


def main() -> None:
    import os
    # stdio for a local agent; streamable-http to serve a fleet from the cluster
    mcp.run(os.environ.get("TASKQ_MCP_TRANSPORT", "stdio"))


if __name__ == "__main__":
    main()
