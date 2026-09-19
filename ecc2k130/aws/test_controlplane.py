#!/usr/bin/env python3
"""Tests for the RDS/ElastiCache/S3 control plane.

Two decisions shape this suite:

  * **The SQL is executed, not asserted about.**  Every test runs the real
    statements from `slots.py` and `coordinator.py` against SQLite, which
    supports the `ON CONFLICT`, `RETURNING` and `INSERT ... SELECT ... WHERE
    EXISTS` forms they use.  A test that checked the SQL text would pass on a
    claim that cannot win a race; these run the race, with threads.
  * **No network, no AWS SDK, no Redis.**  The cache tests drive a fake that
    implements the handful of commands used, including the failure modes that
    matter -- timeouts, garbage in a key, an evicted lock -- and the last test
    asserts in a clean subprocess that importing the package pulls in neither
    `redis` nor `boto3`.

Run: python3 -m unittest discover -s aws -p test_controlplane.py -v
"""

import json
import os
import subprocess
import sys
import tempfile
import threading
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from controlplane import cache as cachemod            # noqa: E402
from controlplane.cache import Cache, NullCache       # noqa: E402
from controlplane.config import Config                # noqa: E402
from controlplane.coordinator import Coordinator, FenceError, workerInfo  # noqa: E402
from controlplane.db import Database, isTransient, sqliteDatabase  # noqa: E402
from controlplane.schema import MIGRATIONS, migrate, pending  # noqa: E402
from controlplane.slots import PostgresSlots          # noqa: E402
from controlplane.store import LocalObjectStore, ObjectExists  # noqa: E402


class Clock:
    """A clock the tests move by hand; leases are time, and time is the test."""

    def __init__(self, now=1_800_000_000):
        self.now = float(now)

    def __call__(self):
        return self.now

    def advance(self, seconds):
        self.now += seconds
        return self.now


def database(path=None):
    db = sqliteDatabase(path or ":memory:")
    migrate(db)
    return db


def gpuInfo(instance="i-1", gpu=0, name="NVIDIA RTX PRO 6000 Blackwell Server Edition",
            family="g7e", instanceType="g7e.2xlarge"):
    return workerInfo(gpu=gpu, gpuName=name, gpuFamily=family, instance=instance,
                      instanceType=instanceType)


# ---------------------------------------------------------------------------
class SchemaTests(unittest.TestCase):
    def test_migrate_is_idempotent(self):
        db = sqliteDatabase()
        first = migrate(db)
        self.assertEqual(first, [ident for ident, _ in MIGRATIONS])
        self.assertEqual(migrate(db), [])
        self.assertEqual(pending(db), [])

    def test_migrate_resumes_after_a_partial_apply(self):
        """A migration that died mid-list re-runs without complaining.

        Every statement is `IF NOT EXISTS` and the marker row is written last,
        so the failure mode a deploy actually has -- the host going away
        between two `CREATE INDEX` statements -- costs a repeat, not a manual
        repair.
        """
        db = sqliteDatabase()
        ident, statements = MIGRATIONS[0]
        db.run(statements[0])                      # as if it had died here
        self.assertIn(ident, migrate(db))
        self.assertEqual(pending(db), [])


class DatabaseTests(unittest.TestCase):
    def test_transient_errors_are_recognised_by_shape(self):
        class OperationalError(Exception):
            pass

        self.assertTrue(isTransient(OperationalError("server closed the connection")))
        self.assertTrue(isTransient(RuntimeError("terminating connection due to failover")))
        self.assertFalse(isTransient(RuntimeError("syntax error at or near")))

    def test_a_failover_is_retried_on_a_fresh_connection(self):
        """One dead connection must not surface as a dead campaign."""
        attempts = {"n": 0}
        real = []

        def connect():
            import sqlite3

            attempts["n"] += 1
            conn = sqlite3.connect(":memory:")
            if attempts["n"] == 1:
                class Dying:
                    def cursor(self):
                        raise _Operational("server closed the connection unexpectedly")

                    def rollback(self):
                        raise _Operational("connection already closed")

                    def close(self):
                        pass

                return Dying()
            real.append(conn)
            return conn

        db = Database(connect, paramstyle="qmark", dialect="sqlite", sleep=lambda _s: None)
        self.assertEqual(db.fetchOne("SELECT 1"), (1,))
        self.assertEqual(attempts["n"], 2)

    def test_a_real_error_is_not_retried(self):
        db = sqliteDatabase()
        with self.assertRaises(Exception):
            db.run("SELECT * FROM a_table_that_does_not_exist")


class _Operational(Exception):
    pass


# ---------------------------------------------------------------------------
class SlotTests(unittest.TestCase):
    def setUp(self):
        self.clock = Clock()
        self.db = database()
        self.slots = PostgresSlots(self.db, campaign="test", clock=self.clock)

    def test_first_claim_creates_slot_zero_with_fence_one(self):
        slot, fence = self.slots.claimLease("worker-a", gpuInfo())
        self.assertEqual((slot, fence), (0, 1))
        item = self.slots.get(0)
        self.assertEqual(item["owner"], "worker-a")
        self.assertEqual(item["state"], "active")
        self.assertEqual(item["leaseUntil"], int(self.clock.now) + 180)

    def test_a_live_lease_is_never_handed_to_a_second_worker(self):
        self.slots.claimLease("worker-a", gpuInfo())
        slot, _ = self.slots.claimLease("worker-b", gpuInfo(instance="i-2"))
        self.assertEqual(slot, 1)              # a new run id, not a's
        self.assertEqual(self.slots.get(0)["owner"], "worker-a")

    def test_an_expired_lease_is_resumed_and_the_fence_moves(self):
        slot, fence = self.slots.claimLease("worker-a", gpuInfo())
        self.clock.advance(181)
        resumed, newFence = self.slots.claimLease("worker-b", gpuInfo(instance="i-2"))
        self.assertEqual(resumed, slot)
        self.assertEqual(newFence, fence + 1)
        self.assertEqual(self.slots.get(slot)["owner"], "worker-b")

    def test_a_retired_slot_keeps_its_run_id_forever(self):
        slot, _ = self.slots.claimLease("worker-a", gpuInfo())
        self.slots.retire(slot, "checkpoint refused by the client")
        self.clock.advance(10_000)
        again, _ = self.slots.claimLease("worker-b", gpuInfo(instance="i-2"))
        self.assertNotEqual(again, slot)
        self.assertEqual(self.slots.get(slot)["state"], "retired")

    def test_heartbeat_requires_the_lease_and_the_fence(self):
        slot, fence = self.slots.claimLease("worker-a", gpuInfo())
        self.assertTrue(self.slots.heartbeat(slot, "worker-a", {"iters": 5}, fence=fence))
        self.assertEqual(self.slots.get(slot)["iters"], 5)
        self.assertFalse(self.slots.heartbeat(slot, "worker-b", {}, fence=fence))
        self.clock.advance(181)
        self.assertFalse(self.slots.heartbeat(slot, "worker-a", {}, fence=fence))

    def test_heartbeat_drops_fields_that_are_not_columns(self):
        """A registry field has to exist in the schema to be written.

        The other registries take arbitrary attributes, so a caller may pass
        one; it is dropped rather than interpolated into the statement.
        """
        slot, fence = self.slots.claimLease("worker-a", gpuInfo())
        self.assertTrue(self.slots.heartbeat(
            slot, "worker-a", {"iters": 7, "owner": "worker-b", "fence": 99,
                               "nonsense; DROP TABLE rho_slots": 1}, fence=fence))
        item = self.slots.get(slot)
        self.assertEqual((item["owner"], item["fence"], item["iters"]), ("worker-a", fence, 7))

    def test_release_hands_back_the_run_id_but_not_somebody_elses(self):
        slot, fence = self.slots.claimLease("worker-a", gpuInfo())
        self.assertFalse(self.slots.release(slot, "worker-b"))
        self.assertTrue(self.slots.release(slot, "worker-a", state="idle",
                                           extra={"reason": "spot interruption"}))
        item = self.slots.get(slot)
        self.assertEqual((item["state"], item["leaseUntil"]), ("idle", 0))
        again, newFence = self.slots.claimLease("worker-b", gpuInfo(instance="i-2"))
        self.assertEqual((again, newFence), (slot, fence + 1))

    def test_cpu_workers_never_resume_a_gpu_checkpoint(self):
        """The rule lives in worker.py; this checks it is the one in force.

        Loading a packed GPU checkpoint into the host client is exit 6, and
        the supervisor then retires the run id -- so a registry that let a CPU
        worker resume a GPU slot would burn run ids one interruption at a
        time.
        """
        slot, _ = self.slots.claimLease("worker-gpu", gpuInfo())
        self.clock.advance(181)
        cpu = workerInfo(gpu=0, gpuName="cpu/8", gpuFamily="cpu", instance="i-cpu")
        got, _ = self.slots.claimLease("worker-cpu", cpu, newOnly=True)
        self.assertNotEqual(got, slot)

    def test_an_ada_worker_does_not_take_an_untagged_blackwell_slot(self):
        self.db.run("INSERT INTO rho_slots (campaign_id, slot, state, lease_until, fence) "
                    "VALUES ('test', 0, 'idle', 0, 3)")
        ada = workerInfo(gpu=0, gpuName="NVIDIA L40S", gpuFamily="g6e", instance="i-ada")
        slot, _ = self.slots.claimLease("worker-ada", ada)
        self.assertNotEqual(slot, 0)

    def test_concurrent_claims_never_share_a_run_id(self):
        """The property the whole registry exists for, run as a race.

        Sixteen threads claim at once against one database, twice: once with
        no slots to resume (so they race on creation) and once with every slot
        expired (so they race on the conditional UPDATE).
        """
        path = os.path.join(tempfile.mkdtemp(), "slots.db")
        db = database(path)
        slots = PostgresSlots(db, campaign="test", clock=self.clock)
        got, errors = [], []
        lock = threading.Lock()

        def claim(name):
            try:
                slot, fence = slots.claimLease(name, gpuInfo(instance=name))
                with lock:
                    got.append((slot, fence))
            except Exception as exc:                       # pragma: no cover
                with lock:
                    errors.append(exc)

        threads = [threading.Thread(target=claim, args=("w%d" % i,)) for i in range(16)]
        for t in threads:
            t.start()
        for t in threads:
            t.join()
        self.assertEqual(errors, [])
        self.assertEqual(len({slot for slot, _ in got}), 16)

        self.clock.advance(1000)                            # every lease expires
        got2 = []

        def reclaim(name):
            try:
                slot, fence = slots.claimLease(name, gpuInfo(instance=name))
                with lock:
                    got2.append((slot, fence))
            except Exception as exc:                        # pragma: no cover
                with lock:
                    errors.append(exc)

        threads = [threading.Thread(target=reclaim, args=("r%d" % i,)) for i in range(16)]
        for t in threads:
            t.start()
        for t in threads:
            t.join()
        self.assertEqual(errors, [])
        self.assertEqual(len({slot for slot, _ in got2}), 16)
        # Sixteen workers, sixteen live leases, and every fence strictly
        # greater than the one the slot was first issued under.
        live = [it for it in slots.scan(fresh=True) if it["leaseUntil"] >= self.clock.now]
        self.assertEqual(len(live), 16)
        self.assertTrue(all(it["fence"] >= 2 for it in live))
        db.close()

    def test_events_record_every_transition(self):
        slot, fence = self.slots.claimLease("worker-a", gpuInfo())
        self.slots.release(slot, "worker-a", state="idle")
        kinds = [e["kind"] for e in self.slots.events(slot)]
        self.assertEqual(sorted(kinds), ["create", "release"])

    def test_expired_lists_what_a_spot_interruption_took(self):
        self.slots.claimLease("worker-a", gpuInfo())
        self.assertEqual(self.slots.expired(), [])
        self.clock.advance(181)
        self.assertEqual([it["slot"] for it in self.slots.expired()], [0])


# ---------------------------------------------------------------------------
class FakeRedis:
    """Just enough Redis: strings with TTLs, a hash, a sorted set, EVAL.

    TTLs are evaluated against the test's clock, so an expiring lock or a
    fallen-out heartbeat is a clock move rather than a sleep.
    """

    def __init__(self, clock):
        self.clock = clock
        self.values = {}       # key -> (value, expiresAt or None)
        self.hashes = {}
        self.zsets = {}
        self.fail = None       # set to an exception to make every call raise
        self.calls = 0

    # -- helpers
    def _check(self):
        self.calls += 1
        if self.fail:
            raise self.fail

    def _live(self, key):
        item = self.values.get(key)
        if item is None:
            return None
        value, expires = item
        if expires is not None and self.clock() >= expires:
            del self.values[key]
            return None
        return value

    # -- commands
    def get(self, key):
        self._check()
        return self._live(key)

    def set(self, key, value, ex=None, px=None, nx=False):
        self._check()
        if nx and self._live(key) is not None:
            return None
        expires = None
        if ex:
            expires = self.clock() + ex
        elif px:
            expires = self.clock() + px / 1000.0
        self.values[key] = (value if isinstance(value, (str, bytes)) else str(value), expires)
        return True

    def delete(self, *keys):
        self._check()
        n = 0
        for key in keys:
            n += 1 if self.values.pop(key, None) is not None else 0
        return n

    def hset(self, name, field, value):
        self._check()
        self.hashes.setdefault(name, {})[field] = value
        return 1

    def hdel(self, name, field):
        self._check()
        return 1 if self.hashes.get(name, {}).pop(field, None) is not None else 0

    def hmget(self, name, fields):
        self._check()
        table = self.hashes.get(name, {})
        return [table.get(f) for f in fields]

    def zadd(self, name, mapping):
        self._check()
        self.zsets.setdefault(name, {}).update(mapping)
        return len(mapping)

    def zrem(self, name, member):
        self._check()
        return 1 if self.zsets.get(name, {}).pop(member, None) is not None else 0

    def zrangebyscore(self, name, low, high):
        self._check()
        top = float("inf") if high in ("+inf", b"+inf") else float(high)
        return [m for m, s in sorted(self.zsets.get(name, {}).items(), key=lambda kv: kv[1])
                if float(low) <= s <= top]

    def expire(self, name, ttl):
        self._check()
        return True

    def eval(self, script, numkeys, *args):
        self._check()
        key, token = args[0], args[1]
        if self._live(key) != token:
            return 0
        if "del" in script:
            del self.values[key]
            return 1
        value, _ = self.values[key]
        self.values[key] = (value, self.clock() + float(args[2]) / 1000.0)
        return 1

    def pipeline(self, transaction=False):
        return _Pipeline(self)

    def close(self):
        pass


class _Pipeline:
    def __init__(self, client):
        self.client = client
        self.queued = []

    def __getattr__(self, name):
        def queue(*args, **kw):
            self.queued.append((name, args, kw))
            return self
        return queue

    def execute(self):
        return [getattr(self.client, name)(*args, **kw) for name, args, kw in self.queued]


class CacheTests(unittest.TestCase):
    def setUp(self):
        self.clock = Clock()
        self.redis = FakeRedis(self.clock)
        self.cache = Cache(self.redis, campaign="test", clock=self.clock)

    def test_keys_carry_a_hash_tag_so_one_campaign_lands_on_one_shard(self):
        self.assertTrue(self.cache.key("fleet").startswith("rho:{test}:"))

    def test_a_dead_cluster_stops_being_called_after_a_few_failures(self):
        self.redis.fail = TimeoutError("timed out")
        for _ in range(3):
            self.assertIsNone(self.cache.getJson("anything"))
        calls = self.redis.calls
        for _ in range(10):
            self.assertIsNone(self.cache.getJson("anything"))
        self.assertEqual(self.redis.calls, calls)           # skipped, not retried
        self.assertTrue(self.cache.stats()["cooling"])
        # And it comes back on its own once the cooldown passes.
        self.redis.fail = None
        self.clock.advance(31)
        self.cache.setJson("anything", {"a": 1})
        self.assertEqual(self.cache.getJson("anything"), {"a": 1})

    def test_undecodable_bytes_are_a_miss_not_an_error(self):
        self.redis.values[self.cache.key("junk")] = ("not json", None)
        self.assertIsNone(self.cache.getJson("junk"))

    def test_fleet_distinguishes_empty_from_unknown(self):
        self.assertEqual(self.cache.fleet(), [])            # healthy, nobody beating
        self.cache.beat("w1", {"slot": 3})
        self.assertEqual([w["owner"] for w in self.cache.fleet()], ["w1"])
        self.clock.advance(10_000)                          # w1 stops beating
        self.assertEqual(self.cache.fleet(), [])
        self.redis.fail = TimeoutError("timed out")
        self.assertIsNone(self.cache.fleet())               # unknown, ask the database

    def test_a_lock_is_released_only_by_its_holder(self):
        first = self.cache.lock("merge", ttl=60)
        self.assertTrue(first.acquire())
        second = self.cache.lock("merge", ttl=60)
        self.assertFalse(second.acquire())
        self.assertFalse(second.release())
        self.assertTrue(first.release())
        self.assertTrue(second.acquire())

    def test_a_stalled_holder_cannot_delete_the_next_holders_lock(self):
        stalled = self.cache.lock("merge", ttl=30)
        self.assertTrue(stalled.acquire())
        self.clock.advance(31)                              # the lock expires
        nextHolder = self.cache.lock("merge", ttl=30)
        self.assertTrue(nextHolder.acquire())
        self.assertFalse(stalled.release())                 # token mismatch
        self.assertIsNotNone(self.redis.get(self.cache.key("lock", "merge")))

    def test_a_lock_is_granted_when_the_cache_cannot_answer(self):
        """Work continues when the cache is down; it is only ever advisory."""
        self.redis.fail = TimeoutError("timed out")
        self.assertTrue(self.cache.lock("merge").acquire())

    def test_point_hints_are_hints(self):
        self.assertFalse(self.cache.seenPoint("abc"))       # first sighting
        self.assertTrue(self.cache.seenPoint("abc"))        # offered again
        self.clock.advance(3601)
        self.assertFalse(self.cache.seenPoint("abc"))       # the window is an hour
        self.redis.fail = TimeoutError("timed out")
        self.assertIsNone(self.cache.seenPoint("abc"))      # unknown, never False

    def test_null_cache_answers_unknown_to_everything(self):
        null = NullCache()
        self.assertIsNone(null.fleet())
        self.assertIsNone(null.seenPoint("abc"))
        self.assertIsNone(null.slotSnapshot())
        self.assertTrue(null.lock("merge").acquire())
        self.assertFalse(null.stats()["enabled"])


# ---------------------------------------------------------------------------
class CoordinatorTests(unittest.TestCase):
    def setUp(self):
        self.clock = Clock()
        self.db = database()
        self.root = tempfile.mkdtemp()
        self.store = LocalObjectStore(os.path.join(self.root, "s3"))
        self.redis = FakeRedis(self.clock)
        self.cache = Cache(self.redis, campaign="test", clock=self.clock)
        self.coord = Coordinator(self.db, self.store, self.cache, campaign="test",
                                 owner="worker-a", clock=self.clock)

    def other(self, owner="worker-b"):
        return Coordinator(self.db, self.store, self.cache, campaign="test",
                           owner=owner, clock=self.clock)

    def test_an_object_is_in_s3_before_it_is_in_the_ledger(self):
        lease = self.coord.acquire(gpuInfo())
        row = self.coord.publishObject(lease, "dp/slot-00000/a-0.bin", b"x" * 64,
                                       streamId="a", offset=0)
        self.assertEqual(row["records"], 2)
        self.assertEqual(self.store.getBytes("dp/slot-00000/a-0.bin"), b"x" * 64)
        ledger = self.coord.objects(slot=lease.slot)
        self.assertEqual([o["objectKey"] for o in ledger], ["dp/slot-00000/a-0.bin"])

    def test_re_uploading_the_same_object_does_not_double_count(self):
        lease = self.coord.acquire(gpuInfo())
        for _ in range(3):
            self.coord.publishObject(lease, "dp/slot-00000/a-0.bin", b"x" * 64)
        rows = self.coord.objects(slot=lease.slot)
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["records"], 2)

    def test_different_bytes_at_one_key_are_refused(self):
        lease = self.coord.acquire(gpuInfo())
        self.coord.publishObject(lease, "dp/slot-00000/a-0.bin", b"x" * 64)
        with self.assertRaises(ObjectExists):
            self.coord.publishObject(lease, "dp/slot-00000/a-0.bin", b"y" * 64)

    def test_a_zombie_worker_cannot_write_to_the_ledger(self):
        """The reason the fence exists.

        Worker A is frozen past its lease; B claims the slot and walks it.
        A wakes up still holding a lease object that looks valid to A, and
        every write it attempts is rejected -- not because A noticed, but
        because the database did.
        """
        lease = self.coord.acquire(gpuInfo())
        self.clock.advance(181)
        b = self.other()
        b.acquire(gpuInfo(instance="i-2"))
        with self.assertRaises(FenceError):
            self.coord.publishObject(lease, "dp/slot-00000/zombie-0.bin", b"x" * 32)
        with self.assertRaises(FenceError):
            self.coord.checkFence(lease)
        self.assertFalse(self.coord.renew(lease))

    def test_a_zombie_worker_cannot_pull_back_the_checkpoint_pointer(self):
        lease = self.coord.acquire(gpuInfo())
        self.coord.setCheckpoint(lease, "ck/slot-00000/1.ck", iters=1000)
        self.clock.advance(181)
        b = self.other()
        newLease = b.acquire(gpuInfo(instance="i-2"))
        b.setCheckpoint(newLease, "ck/slot-00000/2.ck", iters=2000)
        with self.assertRaises(FenceError):
            self.coord.setCheckpoint(lease, "ck/slot-00000/1.ck", iters=1000)
        pointer = b.checkpoint(newLease.slot)
        self.assertEqual((pointer["objectKey"], pointer["iters"]), ("ck/slot-00000/2.ck", 2000))

    def test_a_checkpoint_moves_forward_under_the_same_lease(self):
        lease = self.coord.acquire(gpuInfo())
        self.coord.setCheckpoint(lease, "ck/1.ck", iters=10)
        self.coord.setCheckpoint(lease, "ck/2.ck", iters=20)
        self.assertEqual(self.coord.checkpoint(lease.slot)["iters"], 20)

    def test_heartbeat_keeps_the_lease_and_publishes_the_fleet(self):
        lease = self.coord.acquire(gpuInfo())
        self.clock.advance(120)
        self.assertTrue(self.coord.heartbeat(lease, {"iters": 10 ** 12, "points": 40,
                                                     "rate": 14.1e9}))
        fleet = self.coord.fleet()
        self.assertEqual(fleet["source"], "cache")
        self.assertEqual(fleet["workers"][0]["slot"], lease.slot)
        self.assertEqual(self.coord.slots.get(lease.slot)["iters"], 10 ** 12)

    def test_the_fleet_falls_back_to_the_database_when_the_cache_is_down(self):
        lease = self.coord.acquire(gpuInfo())
        self.coord.heartbeat(lease, {"iters": 5})
        self.redis.fail = TimeoutError("timed out")
        fleet = self.coord.fleet()
        self.assertEqual(fleet["source"], "database")
        self.assertEqual([w["slot"] for w in fleet["workers"]], [lease.slot])

    def test_the_walk_survives_a_cache_outage_entirely(self):
        """A cache failure may cost a dashboard panel and nothing else."""
        self.redis.fail = TimeoutError("timed out")
        lease = self.coord.acquire(gpuInfo())
        self.assertTrue(self.coord.heartbeat(lease, {"iters": 1}))
        self.coord.publishObject(lease, "dp/slot-00000/a-0.bin", b"x" * 32)
        self.assertEqual(len(self.coord.objects(slot=lease.slot)), 1)

    def test_status_reports_every_store(self):
        lease = self.coord.acquire(gpuInfo())
        self.coord.publishObject(lease, "dp/slot-00000/a-0.bin", b"x" * 64)
        self.coord.heartbeat(lease, {"iters": 7})
        status = self.coord.status()
        self.assertEqual(status["slots"], {"total": 1, "leased": 1,
                                           "byState": {"active": 1}, "expired": 0})
        self.assertEqual(status["corpus"]["records"], 2)
        self.assertEqual(status["fleet"]["live"], 1)
        self.clock.advance(181)
        self.assertEqual(self.coord.status()["slots"]["expired"], 1)

    def test_prune_keeps_the_record_of_which_seeds_were_walked(self):
        lease = self.coord.acquire(gpuInfo())
        self.coord.publishObject(lease, "dp/slot-00000/a-0.bin", b"x" * 32)
        self.coord.heartbeat(lease, {"iters": 1})
        self.clock.advance(30 * 86400)
        pruned = self.coord.prune()
        self.assertEqual(pruned["workers"], 1)
        self.assertGreaterEqual(pruned["events"], 1)
        self.assertEqual(len(self.coord.slots.scan(fresh=True)), 1)
        self.assertEqual(len(self.coord.objects()), 1)

    def test_a_released_slot_is_resumed_by_the_next_worker(self):
        lease = self.coord.acquire(gpuInfo())
        self.coord.setCheckpoint(lease, "ck/1.ck", iters=500)
        self.assertTrue(self.coord.release(lease, state="idle",
                                           extra={"reason": "spot interruption"}))
        b = self.other()
        resumed = b.acquire(gpuInfo(instance="i-2"))
        self.assertEqual(resumed.slot, lease.slot)
        self.assertEqual(resumed.fence, lease.fence + 1)
        # The checkpoint the new owner must load is still there, and the new
        # fence is high enough to move it.
        self.assertEqual(b.checkpoint(resumed.slot)["iters"], 500)
        b.setCheckpoint(resumed, "ck/2.ck", iters=900)

    def test_run_id_is_the_slot_plus_one(self):
        lease = self.coord.acquire(gpuInfo())
        self.assertEqual(lease.runId, lease.slot + 1)

    def test_point_hints_reach_the_cache(self):
        self.assertFalse(self.coord.offerPoint("deadbeef"))
        self.assertTrue(self.coord.offerPoint("deadbeef"))


# ---------------------------------------------------------------------------
class ConfigTests(unittest.TestCase):
    def test_database_url_is_read_under_the_name_the_ingest_uses(self):
        config = Config.fromEnv({"DATABASE_URL": "postgresql://h/db", "RHO_CAMPAIGN": "c"})
        self.assertEqual(config.databaseUrl, "postgresql://h/db")
        self.assertEqual(config.campaign, "c")

    def test_no_redis_url_means_no_cache_rather_than_an_error(self):
        self.assertFalse(cachemod.Cache.fromEnv(Config.fromEnv({})).enabled)

    def test_cache_can_be_switched_off_to_get_the_control(self):
        config = Config.fromEnv({"RHO_REDIS_URL": "rediss://x", "RHO_CACHE": "0"})
        self.assertFalse(config.cacheEnabled)
        self.assertFalse(cachemod.Cache.fromEnv(config).enabled)

    def test_object_keys_go_under_the_prefix(self):
        config = Config.fromEnv({"ECC_PREFIX": "campaign-2"})
        self.assertEqual(config.objectKey("dp", "slot-00001"), "campaign-2/dp/slot-00001")


class ImportPurityTests(unittest.TestCase):
    def test_importing_the_package_pulls_in_neither_redis_nor_boto3(self):
        """The lazy imports are a property of the package, asserted as one.

        A control-plane host, a CI runner and a developer laptop all import
        this package; only the first has the AWS SDKs.  An eager import here
        would be found in production, by the fleet.
        """
        code = ("import sys; sys.path.insert(0, %r); import controlplane; "
                "print('redis' in sys.modules, 'boto3' in sys.modules)"
                % os.path.dirname(os.path.abspath(__file__)))
        out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)
        self.assertEqual(out.returncode, 0, out.stderr)
        self.assertEqual(out.stdout.strip(), "False False")


class WorkerContractTests(unittest.TestCase):
    """The registry has to stay a drop-in for the three worker.py already has.

    `worker.py` calls `scan()`, `claim(owner, info, newOnly)`,
    `heartbeat(slot, owner, fields)` and `release(slot, owner, state, extra)`
    positionally.  A signature that drifts here is not caught by any test of
    this package on its own -- it is caught on a GPU, at 3 a.m., by a
    supervisor that cannot renew a lease.
    """

    def test_the_four_methods_take_what_the_supervisor_passes(self):
        import inspect

        import worker

        for name in ("scan", "claim", "heartbeat", "release"):
            reference = inspect.signature(getattr(worker.DynamoSlots, name))
            ours = inspect.signature(getattr(PostgresSlots, name))
            required = [p for p in reference.parameters.values()
                        if p.default is inspect.Parameter.empty]
            have = list(ours.parameters.values())
            self.assertGreaterEqual(len(have), len(required), name)
            for expected, actual in zip(required, have):
                self.assertEqual(expected.name, actual.name, name)

    def test_the_supervisor_can_be_pointed_at_the_control_plane(self):
        clock = Clock()
        db = database()
        slots = PostgresSlots(db, campaign="test", clock=clock)
        # Exactly the call worker.py makes, positionally.
        slot = slots.claim("i-1:gpu0:abc", gpuInfo(), False)
        self.assertTrue(slots.heartbeat(slot, "i-1:gpu0:abc", {"points": 3}))
        self.assertTrue(slots.release(slot, "i-1:gpu0:abc", "idle", {"reason": "done"}))
        self.assertEqual([it["slot"] for it in slots.scan(fresh=True)], [slot])


class CliTests(unittest.TestCase):
    def test_doctor_reports_each_store_separately(self):
        root = tempfile.mkdtemp()
        env = dict(os.environ, ECC_LOCAL_STORE=root, RHO_CAMPAIGN="test")
        env.pop("RHO_REDIS_URL", None)
        env.pop("DATABASE_URL", None)
        here = os.path.dirname(os.path.abspath(__file__))
        run = subprocess.run([sys.executable, "-m", "controlplane", "migrate"],
                             capture_output=True, text=True, cwd=here, env=env)
        self.assertEqual(run.returncode, 0, run.stderr)
        self.assertEqual(json.loads(run.stdout)["pending"], [])
        run = subprocess.run([sys.executable, "-m", "controlplane", "doctor"],
                             capture_output=True, text=True, cwd=here, env=env)
        report = json.loads(run.stdout)
        self.assertTrue(report["ok"], report)
        self.assertTrue(report["database"]["ok"])
        self.assertIn("no cluster configured", report["cache"]["note"])

    def test_status_runs_against_a_local_rehearsal_store(self):
        root = tempfile.mkdtemp()
        env = dict(os.environ, ECC_LOCAL_STORE=root, RHO_CAMPAIGN="test")
        here = os.path.dirname(os.path.abspath(__file__))
        subprocess.run([sys.executable, "-m", "controlplane", "migrate"],
                       capture_output=True, text=True, cwd=here, env=env, check=True)
        run = subprocess.run([sys.executable, "-m", "controlplane", "status"],
                             capture_output=True, text=True, cwd=here, env=env)
        self.assertEqual(run.returncode, 0, run.stderr)
        status = json.loads(run.stdout)
        self.assertEqual(status["slots"]["total"], 0)
        self.assertFalse(status["cache"]["enabled"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
