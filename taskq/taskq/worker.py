"""The worker daemon: pop, claim, execute, complete. One task at a time.

One task per process is deliberate. The point of this queue is measurement,
and two CPU-bound tasks sharing a box measure each other. Scale out with more
pods (or `--cpus` slices on a large node), not threads.

SIGTERM (a pod eviction, a rollout) stops the running child and hands the
task back as an infrastructure failure, so another worker picks it up; it is
never recorded as a result.
"""
from __future__ import annotations

import json
import logging
import os
import signal
import threading
import time
from pathlib import Path
from typing import Any

from . import protocol
from .execute import Fenced, InfraError, RepoCache, execute
from .store import FencedError, Store

log = logging.getLogger("taskq.worker")


class Worker:
    def __init__(self, store: Store, queues: list[str], repos: RepoCache,
                 worker_id: str, labels: dict[str, str] | None = None,
                 artifact_dir: Path | None = None, result_dir: Path | None = None,
                 lease_seconds: float = 60.0, block_seconds: float = 5.0,
                 cpus: list[int] | None = None):
        self.store, self.queues, self.repos = store, queues, repos
        self.id = worker_id
        self.labels = labels or {}
        self.artifact_dir, self.result_dir = artifact_dir, result_dir
        self.lease_ms = int(lease_seconds * 1000)
        self.block_ms = int(block_seconds * 1000)
        self.cpus = cpus
        self._shutdown = threading.Event()
        self.current: str | None = None

    def request_shutdown(self, *_: Any) -> None:
        log.info("shutdown requested")
        self._shutdown.set()

    def _register(self) -> None:
        self.store.register_worker(self.id, {
            "hostname": os.uname().nodename, "queues": ",".join(self.queues),
            "labels": self.labels, "current_task": self.current or "",
            "pid": os.getpid()}, ttl=max(30, self.lease_ms // 1000 * 2))

    def run_forever(self, max_tasks: int | None = None) -> int:
        done = 0
        while not self._shutdown.is_set():
            if self.run_once():
                done += 1
                if max_tasks is not None and done >= max_tasks:
                    break
        return done

    def run_once(self) -> bool:
        """Process at most one task. Returns True if one was handled."""
        self._register()
        msg = self.store.next_message(self.queues, self.id, self.lease_ms, self.block_ms)
        if msg is None:
            return False
        queue, msg_id, task_id = msg
        claimed = self.store.claim(queue, msg_id, task_id, self.id)
        if claimed is None:
            return False
        fence, attempt = claimed
        task = self.store.get_task(task_id)
        spec = task["spec"]
        try:
            # A spec written by a newer taskq may carry fields this worker would
            # silently ignore (e.g. `verify`). Hand it to a worker that knows them.
            protocol.normalize_spec(spec)
        except protocol.SpecError as err:
            log.warning("%s: declining, spec unsupported here: %s", task_id, err)
            self.store.requeue(queue, msg_id, task_id, fence,
                               f"declined by {self.id}: spec unsupported by this worker ({err})",
                               count_attempt=False, worker_id=self.id)
            time.sleep(min(1.0, self.block_ms / 1000))
            return False
        missing = {k: v for k, v in spec["placement"]["require_labels"].items()
                   if self.labels.get(k) != v}
        if missing:
            log.info("%s: declining, lacks labels %s", task_id, missing)
            self.store.requeue(queue, msg_id, task_id, fence,
                               f"declined by {self.id}: lacks {missing}",
                               count_attempt=False, worker_id=self.id)
            time.sleep(min(1.0, self.block_ms / 1000))
            return False
        self.current = task_id
        try:
            self._run(queue, msg_id, task, fence, attempt)
        finally:
            self.current = None
        return True

    def _run(self, queue: str, msg_id: str, task: dict[str, Any],
             fence: int, attempt: int) -> None:
        task_id, spec = task["task_id"], task["spec"]
        state = {"stop": None}
        beat_stop = threading.Event()

        def beat() -> None:
            period = max(0.2, self.lease_ms / 4000)
            while not beat_stop.wait(period):
                try:
                    r = self.store.heartbeat(queue, msg_id, task_id, self.id, fence)
                except Exception as err:  # Redis blip: keep running, lease may lapse
                    log.warning("%s: heartbeat failed: %s", task_id, err)
                    continue
                if r != "ok":
                    state["stop"] = r
                self._register()

        def stop() -> str | None:
            if self._shutdown.is_set():
                return "shutdown"
            return state["stop"]

        hb = threading.Thread(target=beat, daemon=True)
        hb.start()
        started = time.time()
        log.info("%s: running (attempt %d, fence %d)", task_id, attempt, fence)
        try:
            body = execute(spec, task_id, attempt, self.repos, self.artifact_dir,
                           stop=stop, worker_cpus=self.cpus)
        except Fenced:
            log.warning("%s: fenced mid-run; abandoning without writing", task_id)
            return
        except InfraError as err:
            self._infra(queue, msg_id, task, fence, attempt, started, str(err))
            return
        except Exception as err:  # a worker bug is an infrastructure failure
            log.exception("%s: worker error", task_id)
            self._infra(queue, msg_id, task, fence, attempt, started,
                        f"worker error: {type(err).__name__}: {err}")
            return
        finally:
            beat_stop.set()
            hb.join()
        if state["stop"] == "fenced":
            return
        self._complete(queue, msg_id, task, fence, attempt, started, body)

    def _requeue(self, queue, msg_id, task_id, fence, reason) -> None:
        try:
            self.store.requeue(queue, msg_id, task_id, fence, reason,
                               count_attempt=True, worker_id=self.id)
        except FencedError:
            log.warning("%s: requeue refused, fence moved", task_id)

    def _infra(self, queue, msg_id, task, fence, attempt, started, reason) -> None:
        log.warning("%s: infrastructure failure: %s", task["task_id"], reason)
        if attempt < task["max_attempts"] and not self._shutdown.is_set():
            self._requeue(queue, msg_id, task["task_id"], fence, reason)
            return
        body = {"status": "infra_error", "error": reason, "environment": {},
                "source": {**task["spec"]["source"], "resolved_commit": None},
                "setup": [], "runs": [], "summary": None, "artifacts": [],
                "setup_wall_seconds": 0.0}
        if self._shutdown.is_set():
            self._requeue(queue, msg_id, task["task_id"], fence, f"worker shutdown: {reason}")
            return
        self._complete(queue, msg_id, task, fence, attempt, started, body)

    def _complete(self, queue, msg_id, task, fence, attempt, started, body) -> None:
        task_id = task["task_id"]
        result = {
            "schema": protocol.RESULT_SCHEMA_ID,
            "task_id": task_id, "attempt": attempt, "fence": fence,
            "spec_sha256": task["spec_sha256"],
            "status": body["status"],
            "outcome_class": protocol.outcome_class(body["status"]),
            "error": body["error"],
            "worker": {"id": self.id, "hostname": os.uname().nodename,
                       "labels": self.labels},
            "environment": body["environment"],
            "source": body["source"],
            "timing": {"queued_at": float(task["created_at"]),
                       "started_at": started, "finished_at": time.time(),
                       "setup_wall_seconds": body.get("setup_wall_seconds", 0.0)},
            "setup": body["setup"], "runs": body["runs"],
            "summary": body["summary"], "artifacts": body["artifacts"],
            "verification": body.get("verification"),
            "labels": task["labels"],
        }
        if self.result_dir:
            # Mirror before commit: a mirrored file with no Redis result is a
            # fenced-out attempt, visible as such by its fence number.
            d = self.result_dir / task_id
            d.mkdir(parents=True, exist_ok=True)
            (d / f"attempt-{attempt}-fence-{fence}.json").write_text(
                json.dumps(result, indent=2, sort_keys=True))
        try:
            self.store.complete(queue, msg_id, task_id, fence, result)
            log.info("%s: %s", task_id, body["status"])
        except FencedError:
            log.warning("%s: result refused, another worker holds the task", task_id)


def install_signal_handlers(worker: Worker) -> None:
    signal.signal(signal.SIGTERM, worker.request_shutdown)
    signal.signal(signal.SIGINT, worker.request_shutdown)
