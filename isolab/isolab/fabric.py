"""The hub as seen by clients and workers: a NATS JetStream work queue,
fenced job records, write-once results, blobs and artifacts.

Everything a worker writes to a job record is a compare-and-swap against the
revision it last saw, so a worker that lost its lease cannot overwrite the
job when it wakes up. Results are created, never put: the hub refuses a
second one. Execution is at-least-once; the result commit is exactly-once.
"""
from __future__ import annotations

import asyncio
import io
import json
import logging
import os
import time
from pathlib import Path
from typing import Any

import nats
from nats.js import api
from nats.js.errors import (BadRequestError, BucketNotFoundError, KeyNotFoundError,
                            KeyWrongLastSequenceError, NoKeysError, NotFoundError,
                            ObjectNotFoundError)

from . import protocol

log = logging.getLogger("isolab.fabric")
DEFAULT_URL = "nats://127.0.0.1:4222"
WORKER_STALE_S = 60.0


class FencedError(RuntimeError):
    """This worker no longer holds the job; its write was refused."""


class JobNotFound(KeyError):
    pass


def _now() -> float:
    return time.time()


def _loads(b: bytes) -> Any:
    return json.loads(b.decode())


def _dumps(doc: Any) -> bytes:
    return json.dumps(doc, sort_keys=True, separators=(",", ":"), default=str).encode()


class Fabric:
    def __init__(self, url: str | None = None, namespace: str | None = None, token: str | None = None,
                 creds: str | None = None, name: str = "isolab", lease_s: float = 90.0):
        self.urls = [u.strip() for u in (url or os.environ.get("ISOLAB_URL") or DEFAULT_URL).split(",")]
        self.ns = (namespace or os.environ.get("ISOLAB_NAMESPACE") or "ISOLAB").upper()
        self.token = token or os.environ.get("ISOLAB_TOKEN")
        self.creds = creds or os.environ.get("ISOLAB_CREDS")
        self.name = name
        self.lease_s = lease_s
        self.decline_delay_s = 5.0
        self.nc = None
        self.js = None
        self._kv: dict[str, Any] = {}
        self._os: dict[str, Any] = {}
        self._subs: dict[str, Any] = {}

    # ---- names
    @property
    def subj(self) -> str:
        return self.ns.lower()

    def work_subject(self, pool: str) -> str:
        return f"{self.subj}.work.{pool}"

    # ---- connection
    async def connect(self, ensure: bool = True, connect_timeout: float = 10.0) -> "Fabric":
        opts: dict[str, Any] = dict(servers=self.urls, name=self.name, connect_timeout=connect_timeout,
                                    max_reconnect_attempts=-1, reconnect_time_wait=2,
                                    allow_reconnect=True, drain_timeout=3)
        if self.token:
            opts["token"] = self.token
        if self.creds:
            opts["user_credentials"] = self.creds
        self.nc = await nats.connect(**opts)
        self.js = self.nc.jetstream(timeout=15)
        if ensure:
            await self.ensure_resources()
        return self

    async def close(self) -> None:
        if self.nc:
            # close, not drain: a pull subscription mid-fetch makes drain wait out its timeout
            with _suppress():
                await self.nc.flush(timeout=2)
            with _suppress():
                await self.nc.close()
            self.nc = None

    @property
    def max_payload(self) -> int:
        return getattr(self.nc, "max_payload", 1024 * 1024) if self.nc else 1024 * 1024

    async def ensure_resources(self) -> None:
        js = self.js
        for name, subjects, retention, extra in (
                (f"{self.ns}_WORK", [f"{self.subj}.work.*"], api.RetentionPolicy.WORK_QUEUE,
                 {"max_msg_size": 4096, "discard": api.DiscardPolicy.NEW}),
                (f"{self.ns}_EVENTS", [f"{self.subj}.events.>"], api.RetentionPolicy.LIMITS,
                 {"max_msgs": 200_000, "max_age": 30 * 86400})):
            try:
                await js.stream_info(name)
            except NotFoundError:
                await js.add_stream(api.StreamConfig(name=name, subjects=subjects, retention=retention,
                                                     storage=api.StorageType.FILE, **extra))
        for bucket, ttl in (("JOBS", 0), ("RESULTS", 0), ("PROGRESS", 6 * 3600), ("WORKERS", WORKER_STALE_S), ("IDEM", 0)):
            self._kv[bucket] = await self._bucket(f"{self.ns}_{bucket}", ttl)
        for bucket in ("BLOBS", "ARTIFACTS"):
            self._os[bucket] = await self._store(f"{self.ns}_{bucket}")

    async def _bucket(self, name: str, ttl: float):
        try:
            return await self.js.key_value(name)
        except (BucketNotFoundError, NotFoundError):
            cfg = api.KeyValueConfig(bucket=name, history=1, storage=api.StorageType.FILE,
                                     ttl=ttl or None, max_value_size=self.max_payload)
            return await self.js.create_key_value(config=cfg)

    async def _store(self, name: str):
        try:
            return await self.js.object_store(name)
        except (BucketNotFoundError, NotFoundError):
            return await self.js.create_object_store(bucket=name, config=api.ObjectStoreConfig(bucket=name, storage=api.StorageType.FILE))

    def kv(self, bucket: str):
        return self._kv[bucket]

    # ---- job records
    async def get_job(self, job_id: str) -> dict[str, Any] | None:
        try:
            e = await self.kv("JOBS").get(job_id)
        except KeyNotFoundError:
            return None
        rec = _loads(e.value)
        rec["_revision"] = e.revision
        return rec

    async def _cas(self, job_id: str, rec: dict[str, Any], last: int) -> int:
        body = {k: v for k, v in rec.items() if not k.startswith("_")}
        try:
            return await self.kv("JOBS").update(job_id, _dumps(body), last=last)
        except (KeyWrongLastSequenceError, BadRequestError) as err:
            raise FencedError(f"{job_id}: record moved past revision {last}") from err

    async def _event(self, job_id: str, event: str, **fields: Any) -> None:
        with _suppress():
            await self.js.publish(f"{self.subj}.events.{job_id}",
                                  _dumps({"job_id": job_id, "event": event, "at": _now(), **fields}))

    # ---- submitter side
    async def submit(self, spec: dict[str, Any]) -> dict[str, Any]:
        norm = protocol.normalize_spec(spec)
        digest = protocol.spec_sha256(norm)
        job_id = protocol.new_job_id()
        key = norm.get("idempotency_key")
        if key:
            try:
                await self.kv("IDEM").create(key, job_id.encode())
            except (KeyWrongLastSequenceError, BadRequestError):
                existing = (await self.kv("IDEM").get(key)).value.decode()
                prior = await self.get_job(existing)
                if prior and prior["spec_sha256"] != digest:
                    raise protocol.SpecError(f"idempotency_key {key!r} already names {existing} with a different spec")
                return {"job_id": existing, "spec_sha256": digest, "deduplicated": True}
        rec = {"job_id": job_id, "schema": "isolab.job-record/v1", "spec": norm, "spec_sha256": digest,
               "name": norm.get("name"), "pool": norm["pool"], "state": "queued", "attempt": 0, "fence": 0,
               "max_attempts": norm["retry"]["max_attempts"], "worker": None, "created_at": _now(),
               "claimed_at": None, "started_at": None, "finished_at": None, "error": None,
               "cancel_requested": False, "declines": [], "infra_failures": [], "labels": norm["labels"],
               "summary": protocol.summarize_spec(norm)}
        await self.kv("JOBS").create(job_id, _dumps(rec))
        await self.js.publish(self.work_subject(norm["pool"]), job_id.encode(),
                              headers={"Nats-Msg-Id": job_id})
        await self._event(job_id, "queued", pool=norm["pool"])
        return {"job_id": job_id, "spec_sha256": digest, "deduplicated": False}

    async def get_result(self, job_id: str) -> dict[str, Any] | None:
        try:
            e = await self.kv("RESULTS").get(job_id)
        except KeyNotFoundError:
            return None
        return _loads(e.value)

    async def get_progress(self, job_id: str) -> dict[str, Any] | None:
        try:
            e = await self.kv("PROGRESS").get(job_id)
        except KeyNotFoundError:
            return None
        return _loads(e.value)

    async def cancel(self, job_id: str) -> str:
        for _ in range(20):
            rec = await self.get_job(job_id)
            if rec is None:
                return "missing"
            if rec["state"] in protocol.TERMINAL_STATES:
                return rec["state"]
            rev = rec["_revision"]
            if rec["state"] == "queued":
                rec.update(state="cancelled", finished_at=_now(), error="cancelled before it was claimed")
            else:
                rec["cancel_requested"] = True
            try:
                await self._cas(job_id, rec, rev)
            except FencedError:
                await asyncio.sleep(0.05)
                continue
            await self._event(job_id, "cancel_requested")
            return rec["state"]
        return "retry"

    async def list_jobs(self, limit: int = 50, state: str | None = None, pool: str | None = None,
                        worker: str | None = None, labels: dict[str, str] | None = None) -> list[dict[str, Any]]:
        try:
            keys = await self.kv("JOBS").keys()
        except NoKeysError:
            return []
        keys = sorted(keys, reverse=True)
        out = []
        for k in keys:
            rec = await self.get_job(k)
            if rec is None:
                continue
            if state and rec["state"] != state:
                continue
            if pool and rec["pool"] != pool:
                continue
            if worker and rec.get("worker") != worker:
                continue
            if labels and any(rec["labels"].get(a) != b for a, b in labels.items()):
                continue
            out.append({"job_id": rec["job_id"], "name": rec.get("name"), "state": rec["state"], "pool": rec["pool"],
                        "worker": rec.get("worker"), "attempt": rec["attempt"], "created_at": rec["created_at"],
                        "finished_at": rec.get("finished_at"), "error": rec.get("error"),
                        "labels": rec["labels"], "summary": rec.get("summary")})
            if len(out) >= limit:
                break
        return out

    async def wait(self, job_id: str, timeout: float, poll: float = 0.5) -> dict[str, Any] | None:
        deadline = time.monotonic() + timeout
        while True:
            rec = await self.get_job(job_id)
            if rec is None or rec["state"] in protocol.TERMINAL_STATES or time.monotonic() >= deadline:
                return rec
            await asyncio.sleep(min(poll, max(0.0, deadline - time.monotonic())))

    async def workers(self, include_stale: bool = False) -> list[dict[str, Any]]:
        try:
            keys = await self.kv("WORKERS").keys()
        except NoKeysError:
            return []
        out = []
        for k in keys:
            try:
                e = await self.kv("WORKERS").get(k)
            except KeyNotFoundError:
                continue
            inv = _loads(e.value)
            age = _now() - inv.get("heartbeat_at", 0)
            inv["stale"] = age > WORKER_STALE_S
            inv["heartbeat_age_s"] = round(age, 1)
            if include_stale or not inv["stale"]:
                out.append(inv)
        return sorted(out, key=lambda w: w.get("id", ""))

    async def queue_depths(self) -> dict[str, Any]:
        out: dict[str, Any] = {}
        try:
            info = await self.js.stream_info(f"{self.ns}_WORK")
            out["queued_messages"] = info.state.messages
            try:
                consumers = await self.js.consumers_info(f"{self.ns}_WORK")
            except Exception:
                consumers = []
            out["pools"] = {c.config.filter_subject.rsplit(".", 1)[-1]: {"pending": c.num_pending, "running": c.num_ack_pending,
                                                                         "redelivered": c.num_redelivered}
                            for c in consumers if c.config.filter_subject}
        except Exception as err:
            out["error"] = str(err)
        return out

    # ---- blobs and artifacts
    async def has_blob(self, digest: str) -> bool:
        try:
            await self._os["BLOBS"].get_info(f"sha256-{digest}")
            return True
        except (ObjectNotFoundError, NotFoundError):
            return False

    async def put_blob(self, digest: str, source: Path | bytes) -> None:
        data = source if isinstance(source, bytes) else open(source, "rb")
        try:
            await self._os["BLOBS"].put(f"sha256-{digest}", data)
        finally:
            if not isinstance(data, bytes):
                data.close()

    async def get_blob(self, digest: str, dest: Path) -> None:
        dest.parent.mkdir(parents=True, exist_ok=True)
        with open(dest, "wb") as fh:
            await self._os["BLOBS"].get(f"sha256-{digest}", writeinto=fh)

    async def put_artifact(self, job_id: str, relpath: str, source: Path) -> str:
        name = f"{job_id}/{relpath}"
        with open(source, "rb") as fh:
            await self._os["ARTIFACTS"].put(name, fh)
        return name

    async def get_artifact(self, job_id: str, relpath: str, dest: Path) -> Path:
        dest.parent.mkdir(parents=True, exist_ok=True)
        with open(dest, "wb") as fh:
            await self._os["ARTIFACTS"].get(f"{job_id}/{relpath}", writeinto=fh)
        return dest

    async def get_artifact_bytes(self, job_id: str, relpath: str, limit: int | None = None) -> bytes:
        buf = io.BytesIO()
        await self._os["ARTIFACTS"].get(f"{job_id}/{relpath}", writeinto=buf)
        data = buf.getvalue()
        return data[-limit:] if limit else data

    # ---- worker side
    async def register_worker(self, inv: dict[str, Any]) -> None:
        inv = {**inv, "heartbeat_at": _now()}
        await self.kv("WORKERS").put(inv["id"], _dumps(inv))

    async def unregister_worker(self, worker_id: str) -> None:
        with _suppress():
            await self.kv("WORKERS").delete(worker_id)

    async def pull_subscription(self, pool: str):
        if pool not in self._subs:
            subject = self.work_subject(pool)
            cfg = api.ConsumerConfig(durable_name=f"pool-{pool}", ack_policy=api.AckPolicy.EXPLICIT,
                                     ack_wait=self.lease_s, max_deliver=-1, filter_subject=subject,
                                     max_ack_pending=1000, deliver_policy=api.DeliverPolicy.ALL)
            self._subs[pool] = await self.js.pull_subscribe(subject, durable=f"pool-{pool}",
                                                            stream=f"{self.ns}_WORK", config=cfg)
        return self._subs[pool]

    async def fetch(self, pool: str, timeout: float):
        sub = await self.pull_subscription(pool)
        try:
            msgs = await sub.fetch(1, timeout=timeout)
        except (asyncio.TimeoutError, nats.errors.TimeoutError):
            return None
        return msgs[0] if msgs else None

    async def claim(self, job_id: str, worker_id: str, delivered: int) -> tuple[dict[str, Any], int, int] | None:
        """Take the job under a new fence. Returns (record, fence, attempt) or None to drop the message."""
        for _ in range(10):
            rec = await self.get_job(job_id)
            if rec is None:
                return None
            if rec["state"] in protocol.TERMINAL_STATES:
                return None
            rev = rec["_revision"]
            if rec["cancel_requested"] and rec["state"] == "queued":
                rec.update(state="cancelled", finished_at=_now(), error="cancelled before it was claimed")
                try:
                    await self._cas(job_id, rec, rev)
                except FencedError:
                    continue
                await self._event(job_id, "cancelled")
                return None
            if rec["state"] in ("claimed", "running"):
                if delivered <= 1 or rec.get("worker") == worker_id:
                    # a first delivery of a job someone holds, or our own redelivery: not a takeover
                    return None
                rec["infra_failures"].append({"at": _now(), "worker": rec.get("worker"), "fence": rec["fence"],
                                              "reason": "worker lost (lease lapsed)"})
            rec["attempt"] += 1
            if rec["attempt"] > rec["max_attempts"]:
                rec.update(state="dead", finished_at=_now(),
                           error=f"exceeded max_attempts ({rec['max_attempts']}) without completing")
                try:
                    await self._cas(job_id, rec, rev)
                except FencedError:
                    continue
                await self._event(job_id, "dead")
                return None
            rec.update(state="claimed", worker=worker_id, claimed_at=_now(), error=None)
            try:
                new_rev = await self._cas(job_id, rec, rev)
            except FencedError:
                await asyncio.sleep(0.05)
                continue
            rec["fence"] = new_rev
            rec["_revision"] = new_rev
            # record the fence itself (one more CAS so the fence is visible to readers)
            try:
                rec["_revision"] = await self._cas(job_id, rec, new_rev)
            except FencedError:
                continue
            await self._event(job_id, "claimed", worker=worker_id, fence=new_rev, attempt=rec["attempt"])
            return rec, new_rev, rec["attempt"]
        return None

    async def mark_running(self, rec: dict[str, Any]) -> None:
        rec.update(state="running", started_at=_now())
        rec["_revision"] = await self._cas(rec["job_id"], rec, rec["_revision"])
        await self._event(rec["job_id"], "running", worker=rec["worker"])

    async def heartbeat(self, rec: dict[str, Any], msg) -> str:
        """Extend the lease and read the record. Returns ok, cancel or fenced."""
        with _suppress():
            await msg.in_progress()
        cur = await self.get_job(rec["job_id"])
        if cur is None:
            return "fenced"
        if cur.get("worker") != rec["worker"] or cur["fence"] != rec["fence"]:
            return "fenced"
        rec["_revision"] = cur["_revision"]
        rec["cancel_requested"] = cur["cancel_requested"]
        return "cancel" if cur["cancel_requested"] else "ok"

    async def put_progress(self, job_id: str, progress: dict[str, Any]) -> None:
        with _suppress():
            await self.kv("PROGRESS").put(job_id, _dumps({**progress, "at": _now()}))

    async def complete(self, rec: dict[str, Any], result: dict[str, Any], msg) -> None:
        job_id = rec["job_id"]
        protocol.validate_result(result)
        body = _dumps(result)
        if len(body) > self.max_payload - 4096:
            result = compact_result(result, self.max_payload - 8192)
            body = _dumps(result)
        cur = await self.get_job(job_id)
        if cur is None or cur.get("worker") != rec["worker"] or cur["fence"] != rec["fence"]:
            raise FencedError(f"{job_id}: fence {rec['fence']} is stale; result refused")
        try:
            await self.kv("RESULTS").create(job_id, body)
        except (KeyWrongLastSequenceError, BadRequestError) as err:
            raise FencedError(f"{job_id}: a result already exists") from err
        for _ in range(10):
            cur = await self.get_job(job_id)
            if cur is None:
                break
            cur.update(state=result["status"], finished_at=_now(), error=result.get("error"))
            try:
                await self._cas(job_id, cur, cur["_revision"])
                break
            except FencedError:
                await asyncio.sleep(0.05)
        with _suppress():
            await msg.ack_sync()
        await self._event(job_id, result["status"], worker=rec["worker"], fence=rec["fence"])

    async def requeue(self, rec: dict[str, Any], msg, reason: str, count_attempt: bool, delay: float) -> None:
        job_id = rec["job_id"]
        for _ in range(10):
            cur = await self.get_job(job_id)
            if cur is None:
                break
            if cur.get("worker") != rec["worker"] or cur["fence"] != rec["fence"]:
                raise FencedError(f"{job_id}: requeue refused, fence moved")
            cur.update(state="queued", worker=None, error=reason)
            if count_attempt:
                cur["infra_failures"].append({"at": _now(), "worker": rec["worker"], "fence": rec["fence"], "reason": reason})
            else:
                cur["attempt"] = max(0, cur["attempt"] - 1)
            try:
                await self._cas(job_id, cur, cur["_revision"])
                break
            except FencedError:
                await asyncio.sleep(0.05)
        with _suppress():
            await msg.nak(delay=delay)
        await self._event(job_id, "requeued", reason=reason, counted=count_attempt)

    async def decline(self, job_id: str, worker_id: str, reasons: list[str], msg) -> None:
        """Hand back a job this worker cannot place; not an attempt."""
        n = 0
        for _ in range(10):
            cur = await self.get_job(job_id)
            if cur is None:
                break
            cur["declines"] = ([*cur.get("declines", []), {"at": _now(), "worker": worker_id, "reasons": reasons}])[-20:]
            n = len(cur["declines"])
            cur["error"] = f"waiting for a worker: {worker_id} declined ({'; '.join(reasons)})"
            try:
                await self._cas(job_id, cur, cur["_revision"])
                break
            except FencedError:
                await asyncio.sleep(0.05)
        with _suppress():
            await msg.nak(delay=min(self.decline_delay_s * max(1, n), 60.0))

    async def drop(self, msg) -> None:
        with _suppress():
            await msg.ack_sync()


def compact_result(result: dict[str, Any], limit: int) -> dict[str, Any]:
    """Trim the parts of a result that are reproducible from artifacts until it fits."""
    r = json.loads(json.dumps(result, default=str))
    r.setdefault("truncated", [])
    steps = [("host.images", lambda: r["host"].pop("images", None)),
             ("runs[].conditions.other_top", lambda: [run.get("conditions", {}).pop("other_top", None) for run in r["runs"]]),
             ("runs[].conditions.freq_mhz", lambda: [run.get("conditions", {}).pop("freq_mhz", None) for run in r["runs"]]),
             ("fidelity.checks.detail", lambda: [c.pop("detail", None) for c in r["fidelity"]["checks"]]),
             ("host.topology.cpus", lambda: r["host"].get("topology", {}).pop("cpus", None)),
             ("runs[].counters._status", lambda: [(run.get("counters") or {}).pop("_status", None) for run in r["runs"]]),
             ("logs", lambda: r.pop("logs", None))]
    for name, fn in steps:
        if len(_dumps(r)) <= limit:
            break
        fn()
        r["truncated"].append(name)
    return r


class _suppress:
    def __enter__(self):
        return self

    def __exit__(self, et, ev, tb):
        if et is not None:
            log.debug("suppressed %s: %s", et.__name__, ev)
        return et is not None and not issubclass(et, (KeyboardInterrupt, SystemExit))

    async def __aenter__(self):
        return self

    async def __aexit__(self, et, ev, tb):
        return self.__exit__(et, ev, tb)
