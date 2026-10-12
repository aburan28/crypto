r"""Redis layout and the state machine, enforced in Lua.

Keys (all under a namespace prefix, default `taskq`):

    {ns}:q:{queue}          stream; one message per enqueue; consumer group `workers`
    {ns}:task:{id}          hash; the task's state (see FIELDS below)
    {ns}:result:{id}        string; the result JSON, written once (SET NX)
    {ns}:attempts:{id}      list; one JSON line per infrastructure failure
    {ns}:tasks              zset; task ids scored by creation time (ms)
    {ns}:idem:{key}         string; idempotency key -> task id
    {ns}:events             stream; audit feed of every transition (capped)
    {ns}:worker:{id}        hash with TTL; live worker roster

State machine:

    queued --claim--> running --complete--> succeeded | failed | timeout | cancelled
       ^                 |  \--infra error, attempts left--> queued
       |                 |  \--infra error, none left-----> infra_error
       |                 \--worker lost (lease idle) --> reclaimed by another worker
       \--placement declined (not an attempt)
    queued --cancel--> cancelled          (claim beyond max_attempts) --> dead

A lease is the stream message's pending-entry idle time: the running worker
resets it with XCLAIM ... JUSTID every heartbeat, and another worker may
XAUTOCLAIM it once it has been idle longer than the lease. Every claim bumps a
**fence**; completion and requeue are Lua scripts that refuse a stale fence,
so a worker that was frozen past its lease cannot overwrite the task after it
wakes up.

Redis is the system of record here, so run it with AOF persistence
(`appendonly yes`). Workers can additionally mirror every result to a
directory (`--result-dir`) so that losing Redis loses the queue, not results.
"""
from __future__ import annotations

import json
import socket
import time
from typing import Any

import redis

from . import protocol

GROUP = "workers"
EVENTS_MAXLEN = 200_000

_CLAIM = """
-- KEYS: task hash, attempts list  ARGV: worker, msg_id, now
local state = redis.call('HGET', KEYS[1], 'state')
if not state then return {-3, 'missing'} end
if redis.call('HGET', KEYS[1], 'cancel_requested') == '1' and state == 'queued' then
  redis.call('HSET', KEYS[1], 'state', 'cancelled', 'finished_at', ARGV[3])
  return {-1, 'cancelled'}
end
if state ~= 'queued' and state ~= 'running' then return {-1, state} end
if state == 'running' then
  -- reclaimed after its lease lapsed: the previous holder was lost
  redis.call('RPUSH', KEYS[2], cjson.encode({at=ARGV[3], worker=redis.call('HGET', KEYS[1], 'worker'),
             fence=tonumber(redis.call('HGET', KEYS[1], 'fence')), reason='worker lost (lease lapsed)'}))
end
local attempt = redis.call('HINCRBY', KEYS[1], 'attempt', 1)
local maxa = tonumber(redis.call('HGET', KEYS[1], 'max_attempts'))
if attempt > maxa then
  redis.call('HSET', KEYS[1], 'state', 'dead', 'finished_at', ARGV[3],
             'error', 'exceeded max_attempts (' .. maxa .. ') without completing')
  return {-2, 'dead'}
end
local fence = redis.call('HINCRBY', KEYS[1], 'fence', 1)
redis.call('HSET', KEYS[1], 'state', 'running', 'worker', ARGV[1],
           'msg_id', ARGV[2], 'started_at', ARGV[3])
return {fence, tostring(attempt)}
"""

_COMPLETE = """
-- KEYS: task hash, result key, queue stream, events stream
-- ARGV: fence, msg_id, status, result_json, now, task_id
if redis.call('HGET', KEYS[1], 'fence') ~= ARGV[1] then return 0 end
if redis.call('HGET', KEYS[1], 'state') ~= 'running' then return 0 end
if redis.call('SET', KEYS[2], ARGV[4], 'NX') == false then return 0 end
redis.call('HSET', KEYS[1], 'state', ARGV[3], 'finished_at', ARGV[5])
redis.call('XACK', KEYS[3], 'workers', ARGV[2])
redis.call('XADD', KEYS[4], 'MAXLEN', '~', '200000', '*',
           'task_id', ARGV[6], 'event', ARGV[3], 'fence', ARGV[1])
return 1
"""

_REQUEUE = """
-- KEYS: task hash, queue stream, events stream, attempts list
-- ARGV: fence, msg_id, now, task_id, reason, count_attempt('1'|'0'), log_line
if redis.call('HGET', KEYS[1], 'fence') ~= ARGV[1] then return 0 end
if redis.call('HGET', KEYS[1], 'state') ~= 'running' then return 0 end
if ARGV[6] == '0' then
  redis.call('HINCRBY', KEYS[1], 'attempt', -1)
  local d = redis.call('HINCRBY', KEYS[1], 'declines', 1)
  if d > 1000 then
    redis.call('HSET', KEYS[1], 'state', 'dead', 'finished_at', ARGV[3],
               'error', 'unplaceable: declined by 1000 workers')
    redis.call('XACK', KEYS[2], 'workers', ARGV[2])
    return 2
  end
else
  redis.call('RPUSH', KEYS[4], ARGV[7])
end
redis.call('HSET', KEYS[1], 'state', 'queued', 'worker', '', 'error', ARGV[5])
redis.call('XACK', KEYS[2], 'workers', ARGV[2])
redis.call('XADD', KEYS[2], '*', 'task_id', ARGV[4])
redis.call('XADD', KEYS[3], 'MAXLEN', '~', '200000', '*',
           'task_id', ARGV[4], 'event', 'requeued', 'reason', ARGV[5])
return 1
"""

_CANCEL = """
-- KEYS: task hash, events stream  ARGV: now, task_id
local state = redis.call('HGET', KEYS[1], 'state')
if not state then return 'missing' end
if state ~= 'queued' and state ~= 'running' then return state end
redis.call('HSET', KEYS[1], 'cancel_requested', '1')
if state == 'queued' then
  redis.call('HSET', KEYS[1], 'state', 'cancelled', 'finished_at', ARGV[1])
end
redis.call('XADD', KEYS[2], 'MAXLEN', '~', '200000', '*',
           'task_id', ARGV[2], 'event', 'cancel_requested')
return redis.call('HGET', KEYS[1], 'state')
"""


class FencedError(RuntimeError):
    """This worker no longer holds the task; its writes were refused."""


def _now() -> str:
    return f"{time.time():.6f}"


class Store:
    def __init__(self, client: redis.Redis, namespace: str = "taskq"):
        self.r = client
        self.ns = namespace
        self._claim = client.register_script(_CLAIM)
        self._complete = client.register_script(_COMPLETE)
        self._requeue = client.register_script(_REQUEUE)
        self._cancel = client.register_script(_CANCEL)
        self._groups: set[str] = set()

    @classmethod
    def from_url(cls, url: str, namespace: str = "taskq") -> "Store":
        return cls(redis.Redis.from_url(url, decode_responses=True), namespace)

    # -- keys --------------------------------------------------------------
    def k_queue(self, q: str) -> str: return f"{self.ns}:q:{q}"
    def k_task(self, t: str) -> str: return f"{self.ns}:task:{t}"
    def k_result(self, t: str) -> str: return f"{self.ns}:result:{t}"
    def k_attempts(self, t: str) -> str: return f"{self.ns}:attempts:{t}"
    def k_worker(self, w: str) -> str: return f"{self.ns}:worker:{w}"
    @property
    def k_tasks(self) -> str: return f"{self.ns}:tasks"
    @property
    def k_events(self) -> str: return f"{self.ns}:events"

    def ensure_group(self, queue: str) -> None:
        if queue in self._groups:
            return
        try:
            self.r.xgroup_create(self.k_queue(queue), GROUP, id="0", mkstream=True)
        except redis.ResponseError as err:
            if "BUSYGROUP" not in str(err):
                raise
        self._groups.add(queue)

    # -- submitter side ----------------------------------------------------
    def submit(self, spec: dict[str, Any]) -> dict[str, Any]:
        """Enqueue a spec. Returns {task_id, spec_sha256, deduplicated}."""
        norm = protocol.normalize_spec(spec)
        digest = protocol.spec_sha256(norm)
        task_id = protocol.new_task_id()
        key = norm.get("idempotency_key")
        if key:
            idem = f"{self.ns}:idem:{key}"
            if not self.r.set(idem, task_id, nx=True):
                existing = self.r.get(idem)
                prior = self.r.hget(self.k_task(existing), "spec_sha256")
                if prior and prior != digest:
                    raise protocol.SpecError(
                        f"idempotency_key {key!r} already names {existing} "
                        "with a different spec")
                return {"task_id": existing, "spec_sha256": prior,
                        "deduplicated": True}
        queue = norm["queue"]
        self.ensure_group(queue)
        now = _now()
        pipe = self.r.pipeline(transaction=True)
        pipe.hset(self.k_task(task_id), mapping={
            "task_id": task_id, "spec": protocol.canonical_json(norm).decode(),
            "spec_sha256": digest, "queue": queue, "state": "queued",
            "attempt": 0, "fence": 0, "declines": 0,
            "max_attempts": norm["retry"]["max_attempts"],
            "cancel_requested": 0, "created_at": now, "worker": "",
            "labels": json.dumps(norm["labels"], sort_keys=True),
        })
        pipe.zadd(self.k_tasks, {task_id: int(float(now) * 1000)})
        pipe.xadd(self.k_queue(queue), {"task_id": task_id})
        pipe.xadd(self.k_events, {"task_id": task_id, "event": "queued"},
                  maxlen=EVENTS_MAXLEN, approximate=True)
        pipe.execute()
        return {"task_id": task_id, "spec_sha256": digest, "deduplicated": False}

    def get_task(self, task_id: str) -> dict[str, Any] | None:
        h = self.r.hgetall(self.k_task(task_id))
        if not h:
            return None
        h["spec"] = json.loads(h["spec"])
        h["labels"] = json.loads(h.get("labels") or "{}")
        for f in ("attempt", "fence", "declines", "max_attempts"):
            h[f] = int(h.get(f, 0))
        h["cancel_requested"] = h.get("cancel_requested") == "1"
        h["infra_failures"] = [json.loads(x) for x in
                               self.r.lrange(self.k_attempts(task_id), 0, -1)]
        h["has_result"] = bool(self.r.exists(self.k_result(task_id)))
        return h

    def get_result(self, task_id: str) -> dict[str, Any] | None:
        raw = self.r.get(self.k_result(task_id))
        return json.loads(raw) if raw else None

    def cancel(self, task_id: str) -> str:
        return self._cancel(keys=[self.k_task(task_id), self.k_events],
                            args=[_now(), task_id])

    def list_tasks(self, limit: int = 50, state: str | None = None,
                   queue: str | None = None,
                   labels: dict[str, str] | None = None,
                   scan: int = 5000) -> list[dict[str, Any]]:
        """Newest first. Filters are applied over the newest `scan` ids."""
        out = []
        for tid in self.r.zrevrange(self.k_tasks, 0, scan - 1):
            h = self.r.hmget(self.k_task(tid), "state", "queue", "labels",
                             "created_at", "finished_at", "attempt", "worker")
            if h[0] is None:
                continue
            lab = json.loads(h[2] or "{}")
            if state and h[0] != state:
                continue
            if queue and h[1] != queue:
                continue
            if labels and any(lab.get(k) != v for k, v in labels.items()):
                continue
            out.append({"task_id": tid, "state": h[0], "queue": h[1],
                        "labels": lab, "created_at": h[3], "finished_at": h[4],
                        "attempt": int(h[5] or 0), "worker": h[6]})
            if len(out) >= limit:
                break
        return out

    def queue_stats(self) -> dict[str, Any]:
        stats = {}
        for key in self.r.scan_iter(match=f"{self.ns}:q:*", _type="stream"):
            q = key.split(":", 2)[2]
            try:
                groups = self.r.xinfo_groups(key)
            except redis.ResponseError:
                groups = []
            g = next((x for x in groups if x["name"] == GROUP), {})
            stats[q] = {"stream_length": self.r.xlen(key),
                        "pending_running": g.get("pending", 0),
                        "lag_unread": g.get("lag"),
                        "consumers": g.get("consumers", 0)}
        return stats

    def workers(self) -> list[dict[str, Any]]:
        out = []
        for key in self.r.scan_iter(match=f"{self.ns}:worker:*"):
            h = self.r.hgetall(key)
            if h:
                h["labels"] = json.loads(h.get("labels") or "{}")
                out.append(h)
        return sorted(out, key=lambda h: h.get("id", ""))

    def wait(self, task_id: str, timeout: float, poll: float = 0.25) -> dict[str, Any] | None:
        deadline = time.monotonic() + timeout
        while True:
            t = self.get_task(task_id)
            if t is None or t["state"] in protocol.TERMINAL_STATES:
                return t
            if time.monotonic() >= deadline:
                return t
            time.sleep(poll)

    # -- worker side -------------------------------------------------------
    def register_worker(self, worker_id: str, info: dict[str, Any], ttl: int) -> None:
        mapping = {**info, "id": worker_id, "seen_at": _now(),
                   "labels": json.dumps(info.get("labels", {}), sort_keys=True)}
        pipe = self.r.pipeline()
        pipe.hset(self.k_worker(worker_id), mapping=mapping)
        pipe.expire(self.k_worker(worker_id), ttl)
        pipe.execute()

    def next_message(self, queues: list[str], worker_id: str,
                     lease_ms: int, block_ms: int) -> tuple[str, str, str] | None:
        """Return (queue, msg_id, task_id): reclaimed work first, then new."""
        for q in queues:
            self.ensure_group(q)
        for q in queues:
            res = self.r.xautoclaim(self.k_queue(q), GROUP, worker_id,
                                    min_idle_time=lease_ms, start_id="0-0", count=1)
            for msg_id, fields in (res[1] if res else []):
                if fields:  # a deleted entry comes back as None
                    return q, msg_id, fields["task_id"]
                self.r.xack(self.k_queue(q), GROUP, msg_id)
        # Strict priority by queue order: a non-blocking pass first.
        for q in queues:
            got = self.r.xreadgroup(GROUP, worker_id, {self.k_queue(q): ">"}, count=1)
            if got:
                msg_id, fields = got[0][1][0]
                return q, msg_id, fields["task_id"]
        got = self.r.xreadgroup(GROUP, worker_id,
                                {self.k_queue(q): ">" for q in queues},
                                count=1, block=block_ms)
        if got:
            stream, msgs = got[0]
            msg_id, fields = msgs[0]
            return stream.split(":", 2)[2], msg_id, fields["task_id"]
        return None

    def claim(self, queue: str, msg_id: str, task_id: str, worker_id: str) -> tuple[int, int] | None:
        """Take the task under a new fence. None means: ack and move on."""
        fence, attempt = self._claim(keys=[self.k_task(task_id), self.k_attempts(task_id)],
                                     args=[worker_id, msg_id, _now()])
        if int(fence) < 0:
            self.r.xack(self.k_queue(queue), GROUP, msg_id)
            if int(fence) == -2:
                self.r.xadd(self.k_events, {"task_id": task_id, "event": "dead"},
                            maxlen=EVENTS_MAXLEN, approximate=True)
            return None
        self.r.xadd(self.k_events, {"task_id": task_id, "event": "running",
                                    "worker": worker_id, "fence": fence},
                    maxlen=EVENTS_MAXLEN, approximate=True)
        return int(fence), int(attempt)

    def heartbeat(self, queue: str, msg_id: str, task_id: str,
                  worker_id: str, fence: int) -> str:
        """Extend the lease. Returns 'ok', 'cancel' or 'fenced'."""
        cur, cancel = self.r.hmget(self.k_task(task_id), "fence", "cancel_requested")
        if cur != str(fence):
            return "fenced"
        self.r.xclaim(self.k_queue(queue), GROUP, worker_id, 0, [msg_id], justid=True)
        return "cancel" if cancel == "1" else "ok"

    def complete(self, queue: str, msg_id: str, task_id: str, fence: int,
                 result: dict[str, Any]) -> None:
        protocol.validate_result(result)
        ok = self._complete(
            keys=[self.k_task(task_id), self.k_result(task_id),
                  self.k_queue(queue), self.k_events],
            args=[fence, msg_id, result["status"],
                  json.dumps(result, sort_keys=True), _now(), task_id])
        if not ok:
            raise FencedError(f"{task_id}: fence {fence} is stale; result refused")

    def requeue(self, queue: str, msg_id: str, task_id: str, fence: int,
                reason: str, count_attempt: bool, worker_id: str = "") -> None:
        line = json.dumps({"at": _now(), "worker": worker_id,
                           "fence": fence, "reason": reason})
        ok = self._requeue(
            keys=[self.k_task(task_id), self.k_queue(queue), self.k_events,
                  self.k_attempts(task_id)],
            args=[fence, msg_id, _now(), task_id, reason,
                  "1" if count_attempt else "0", line])
        if not ok:
            raise FencedError(f"{task_id}: fence {fence} is stale; requeue refused")


def default_worker_id() -> str:
    return f"{socket.gethostname()}-{int(time.time())}"
