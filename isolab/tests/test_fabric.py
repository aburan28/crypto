"""The queue against a real JetStream server: claims, leases, fences, cancels, declines, blobs."""
import asyncio
import time

import pytest

from isolab import protocol
from isolab.fabric import FencedError


def spec(**over):
    s = {"schema": protocol.SPEC_SCHEMA_ID, "command": {"argv": ["true"]}, "runtime": {"backend": "direct"}}
    s.update(over)
    return s


def result_for(rec, status="succeeded"):
    return {"schema": protocol.RESULT_SCHEMA_ID, "job_id": rec["job_id"], "spec_sha256": rec["spec_sha256"],
            "attempt": rec["attempt"], "fence": rec["fence"], "status": status, "outcome_class": protocol.outcome_class(status),
            "error": None, "worker": {"id": rec["worker"], "hostname": "h"}, "host": {}, "placement": {}, "runtime": {},
            "fidelity": {"policy": "standard", "grade": "C", "contended": False, "checks": [], "violations": []},
            "timing": {"queued_at": rec["created_at"], "started_at": 1.0, "finished_at": 2.0}, "build": [], "runs": [],
            "summary": None, "artifacts": [], "labels": rec["labels"]}


async def test_submit_and_record(fabric):
    out = await fabric.submit(spec(labels={"exp": "1"}, name="n1"))
    rec = await fabric.get_job(out["job_id"])
    assert rec["state"] == "queued" and rec["attempt"] == 0 and rec["spec_sha256"] == out["spec_sha256"]
    assert rec["labels"] == {"exp": "1"} and rec["summary"]["argv"] == ["true"]
    jobs = await fabric.list_jobs(labels={"exp": "1"})
    assert [j["job_id"] for j in jobs] == [out["job_id"]]
    assert await fabric.list_jobs(state="running") == []
    assert await fabric.get_result(out["job_id"]) is None


async def test_idempotency_key(fabric):
    a = await fabric.submit(spec(idempotency_key="EXP/RUN-1"))
    b = await fabric.submit(spec(idempotency_key="EXP/RUN-1"))
    assert b["deduplicated"] and b["job_id"] == a["job_id"]
    with pytest.raises(protocol.SpecError, match="different spec"):
        await fabric.submit(spec(idempotency_key="EXP/RUN-1", resources={"cpus": 2}))


async def test_claim_heartbeat_complete(fabric):
    out = await fabric.submit(spec())
    msg = await fabric.fetch("default", 5)
    assert msg is not None and msg.data.decode() == out["job_id"]
    claimed = await fabric.claim(out["job_id"], "w1", msg.metadata.num_delivered)
    assert claimed is not None
    rec, fence, attempt = claimed
    assert rec["state"] == "claimed" and attempt == 1 and fence == rec["fence"] > 0
    await fabric.mark_running(rec)
    assert (await fabric.get_job(out["job_id"]))["state"] == "running"
    assert await fabric.heartbeat(rec, msg) == "ok"
    await fabric.put_progress(out["job_id"], {"phase": "measure"})
    assert (await fabric.get_progress(out["job_id"]))["phase"] == "measure"
    await fabric.complete(rec, result_for(rec), msg)
    final = await fabric.get_job(out["job_id"])
    assert final["state"] == "succeeded" and final["finished_at"]
    assert (await fabric.get_result(out["job_id"]))["status"] == "succeeded"
    # the message is gone from the work queue
    assert await fabric.fetch("default", 0.5) is None
    # a second result is refused
    with pytest.raises(FencedError):
        await fabric.complete(rec, result_for(rec, "failed"), msg)


async def test_cancel_queued_and_running(fabric):
    q = await fabric.submit(spec())
    assert await fabric.cancel(q["job_id"]) == "cancelled"
    msg = await fabric.fetch("default", 5)
    assert await fabric.claim(q["job_id"], "w1", 1) is None  # cancelled: drop it
    await fabric.drop(msg)
    r = await fabric.submit(spec())
    msg = await fabric.fetch("default", 5)
    rec, fence, _ = await fabric.claim(r["job_id"], "w1", 1)
    await fabric.mark_running(rec)
    assert await fabric.cancel(r["job_id"]) == "running"
    assert await fabric.heartbeat(rec, msg) == "cancel"
    await fabric.complete(rec, result_for(rec, "cancelled"), msg)
    assert (await fabric.get_job(r["job_id"]))["state"] == "cancelled"
    assert await fabric.cancel(r["job_id"]) == "cancelled"
    assert await fabric.cancel("J-nope") == "missing"


async def test_decline_redelivers_and_records_reason(fabric):
    out = await fabric.submit(spec())
    msg = await fabric.fetch("default", 5)
    await fabric.decline(out["job_id"], "w1", ["wants 99 cpus"], msg)
    rec = await fabric.get_job(out["job_id"])
    assert rec["state"] == "queued" and rec["attempt"] == 0
    assert rec["declines"][0]["worker"] == "w1" and "wants 99 cpus" in rec["error"]
    again = await fabric.fetch("default", 5)
    assert again is not None and again.data.decode() == out["job_id"] and again.metadata.num_delivered == 2
    await fabric.drop(again)


async def test_requeue_counts_attempt_and_dead_after_max(fabric):
    out = await fabric.submit(spec(retry={"max_attempts": 2}))
    for expected_attempt in (1, 2):
        msg = await fabric.fetch("default", 5)
        rec, _, attempt = await fabric.claim(out["job_id"], "w1", 1)
        assert attempt == expected_attempt
        await fabric.requeue(rec, msg, "disk on fire", count_attempt=True, delay=0.2)
    rec = await fabric.get_job(out["job_id"])
    assert rec["state"] == "queued" and len(rec["infra_failures"]) == 2
    msg = await fabric.fetch("default", 5)
    assert await fabric.claim(out["job_id"], "w1", 1) is None
    assert (await fabric.get_job(out["job_id"]))["state"] == "dead"
    await fabric.drop(msg)


async def test_lost_worker_takeover_fences_the_old_one(fabric, fabric2):
    out = await fabric.submit(spec(retry={"max_attempts": 3}))
    msg = await fabric.fetch("default", 5)
    rec, fence, _ = await fabric.claim(out["job_id"], "w1", 1)
    await fabric.mark_running(rec)
    # w1 stops heartbeating; the lease (5 s) lapses and the hub redelivers to w2
    msg2 = None
    deadline = time.monotonic() + 20
    while msg2 is None and time.monotonic() < deadline:
        msg2 = await fabric2.fetch("default", 2)
    assert msg2 is not None and msg2.metadata.num_delivered >= 2
    claimed = await fabric2.claim(out["job_id"], "w2", msg2.metadata.num_delivered)
    assert claimed is not None
    rec2, fence2, attempt2 = claimed
    assert attempt2 == 2 and fence2 > fence and rec2["infra_failures"][0]["reason"].startswith("worker lost")
    # the old worker wakes up: fenced everywhere
    assert await fabric.heartbeat(rec, msg) == "fenced"
    with pytest.raises(FencedError):
        await fabric.complete(rec, result_for(rec), msg)
    with pytest.raises(FencedError):
        await fabric.requeue(rec, msg, "x", True, 0.1)
    await fabric2.mark_running(rec2)
    await fabric2.complete(rec2, result_for(rec2), msg2)
    final = await fabric.get_result(out["job_id"])
    assert final["fence"] == fence2 and final["attempt"] == 2


async def test_first_delivery_of_a_held_job_is_not_a_takeover(fabric, fabric2):
    out = await fabric.submit(spec())
    msg = await fabric.fetch("default", 5)
    rec, _, _ = await fabric.claim(out["job_id"], "w1", 1)
    assert await fabric2.claim(out["job_id"], "w2", 1) is None
    await fabric.complete(rec, result_for(rec), msg)


async def test_workers_roster_and_blobs_and_artifacts(fabric, tmp_path):
    await fabric.register_worker({"id": "w1", "pools": ["default"]})
    ws = await fabric.workers()
    assert ws[0]["id"] == "w1" and not ws[0]["stale"]
    await fabric.unregister_worker("w1")
    assert await fabric.workers() == []
    digest = "ab" * 32
    assert not await fabric.has_blob(digest)
    await fabric.put_blob(digest, b"payload")
    assert await fabric.has_blob(digest)
    dest = tmp_path / "blob"
    await fabric.get_blob(digest, dest)
    assert dest.read_bytes() == b"payload"
    big = tmp_path / "big.bin"
    big.write_bytes(b"x" * (300 * 1024))
    key = await fabric.put_artifact("J-1", "out/big.bin", big)
    assert key == "J-1/out/big.bin"
    got = await fabric.get_artifact("J-1", "out/big.bin", tmp_path / "got.bin")
    assert got.read_bytes() == big.read_bytes()
    assert await fabric.get_artifact_bytes("J-1", "out/big.bin", 10) == b"x" * 10


async def test_wait_times_out_then_returns(fabric):
    out = await fabric.submit(spec())
    t0 = time.monotonic()
    rec = await fabric.wait(out["job_id"], 0.6)
    assert rec["state"] == "queued" and 0.5 <= time.monotonic() - t0 < 3
    assert "queued_messages" in await fabric.queue_depths()
