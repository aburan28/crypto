#!/usr/bin/env python3
"""The control plane against real Postgres, real Redis and a real S3 endpoint.

`test_controlplane.py` runs the package's SQL against SQLite, its cache logic
against a fake, and its object store against a directory.  That is the right
trade for a suite that must pass on a laptop with no services, and there are
things it structurally cannot see:

  * SQLite has no `pg_advisory_lock`, so nothing there covers what two hosts
    running `migrate` in the same rollout do to each other.
  * SQLite is dynamically typed: a `bigint` that must hold 2**62 and an
    `integer` that cannot are the same column to it.
  * A fake Redis expires keys against the test's clock.  A real one expires
    them against its own, runs the unlock script itself, and refuses the
    connection when it is not there.
  * A directory has POSIX rename semantics, no HTTP in front of it, and no
    1,000-key page limit on a listing.

So this file does two things.  It re-runs the unit suite's own test bodies
against the real services -- the same assertions with different fixtures,
which is what writing them against interfaces was for -- and it adds the
cases only a real server can answer.

Nothing here is skipped quietly.  With `ECC_INTEGRATION` unset the module
skips as a whole and says why; with it set, a service that cannot be reached
is an error rather than a skip, because a green run that tested nothing is
worse than a red one.

Run:  make -C ecc2k130 test-integration     (brings the services up first)
"""

import json
import os
import subprocess
import sys
import threading
import time
import unittest
import uuid
from concurrent.futures import ThreadPoolExecutor

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from controlplane.cache import Cache                  # noqa: E402
from controlplane.config import Config                # noqa: E402
from controlplane.coordinator import Coordinator      # noqa: E402
from controlplane.db import Database                  # noqa: E402
from controlplane.schema import MIGRATIONS, migrate, pending  # noqa: E402
from controlplane.slots import PostgresSlots          # noqa: E402
from controlplane.store import ObjectStore, ObjectExists      # noqa: E402

import test_controlplane as unit                      # noqa: E402

ENABLED = (os.environ.get("ECC_INTEGRATION") or "").strip().lower() in ("1", "true", "yes", "on")
if not ENABLED:
    raise unittest.SkipTest(
        "integration services not requested: set ECC_INTEGRATION=1 with DATABASE_URL, "
        "RHO_REDIS_URL and ECC_BUCKET, or run `make -C ecc2k130 test-integration`")

CONFIG = None
DB = None
CLIENT = None          # one boto3 S3 client, shared: botocore clients are thread-safe
LOCAL_PREFIX = "it"    # every key this suite writes lives under it


def setUpModule():
    """Connect to all three stores, or fail loudly.

    Every failure here is a real one: the caller asked for the integration
    suite, so a store that is not there is the finding, not a reason to pass.
    """
    global CONFIG, DB, CLIENT
    CONFIG = Config.fromEnv()
    missing = [name for name, value in (("DATABASE_URL", CONFIG.databaseUrl),
                                        ("RHO_REDIS_URL", CONFIG.redisUrl),
                                        ("ECC_BUCKET", CONFIG.bucket)) if not value]
    if missing:
        raise RuntimeError("ECC_INTEGRATION is set but %s unset" % ", ".join(missing))

    DB = Database.fromEnv(CONFIG)
    DB.fetchOne("SELECT 1")
    if DB.dialect != "postgres":
        raise RuntimeError("DATABASE_URL did not resolve to Postgres")
    migrate(DB)

    client = liveCache().client
    client.ping()

    store = ObjectStore.fromEnv(CONFIG)
    if CONFIG.s3Endpoint:
        # A bucket is made by Terraform on a deployed campaign; on a throwaway
        # endpoint the suite makes its own rather than needing a setup step
        # somebody can forget.
        store.ensureBucket()
    store.list(limit=1)
    CLIENT = store.client


def tearDownModule():
    if DB is not None:
        DB.close()


def campaign():
    """A campaign nobody else is using.

    Every table is keyed by campaign, so a unique name isolates a test from
    the rest of the suite and from whatever else is in the database -- which
    is the same property the schema gives a rehearsal running beside a live
    campaign.
    """
    return "it-%s" % uuid.uuid4().hex[:12]


def liveDatabase():
    db = Database.fromEnv(CONFIG)
    db.fetchOne("SELECT 1")
    return db


def liveCache():
    """The cache production builds: same timeouts, same retry policy, no retry."""
    cache = Cache.fromEnv(CONFIG)
    if not cache.enabled:
        raise RuntimeError("RHO_REDIS_URL is set but Cache.fromEnv returned a NullCache")
    return cache


def liveStore(prefix=""):
    return ObjectStore(CONFIG.bucket, prefix, client=CLIENT,
                       endpoint=CONFIG.s3Endpoint, region=CONFIG.s3Region)


def deadCache(name, clock):
    """A Cache over a client built the way production builds one, pointed at a
    port nobody answers.  Not a fake failure: a real socket, really refused."""
    config = Config(campaign=name, redisUrl="redis://127.0.0.1:1/0")
    cache = Cache.fromEnv(config)
    return Cache(cache.client, campaign=name, clock=clock)


class LiveRedis:
    """The real client, with the two seams the unit suite's bodies reach for.

    `fail` injects the outage a cache must survive; `calls` counts what
    actually went to the server, which is how the cooldown test proves a dead
    cluster stops being called; `values[key] = (blob, ttl)` plants a value the
    way the fake's dict did, and here that is a real `SET`.  Everything else
    is forwarded, so a test that is not injecting anything is talking to
    Redis over a socket.
    """

    def __init__(self, client):
        self.client = client
        self.fail = None
        self.calls = 0
        self.values = _Values(client)

    def __getattr__(self, name):
        attr = getattr(self.client, name)
        if not callable(attr):
            return attr

        def call(*args, **kw):
            self.calls += 1
            if self.fail is not None:
                raise self.fail
            return attr(*args, **kw)

        return call


class _Values:
    def __init__(self, client):
        self.client = client

    def __setitem__(self, key, value):
        blob, ttl = value if isinstance(value, tuple) else (value, None)
        self.client.set(key, blob, ex=ttl)


# ---------------------------------------------------------------------------
# The unit suite's bodies, against the real thing.
# ---------------------------------------------------------------------------
class SlotTestsOnPostgres(unit.SlotTests):
    """Every assertion in the unit suite's SlotTests, run against Postgres.

    Including the sixteen-thread claim race, which is the one the registry
    exists for and the one whose result on SQLite is only an analogy.
    """

    def setUp(self):
        self.clock = unit.Clock()
        self.campaign = campaign()
        self.db = DB
        self.slots = PostgresSlots(self.db, campaign=self.campaign, clock=self.clock)

    def freshDatabase(self):
        return liveDatabase()


class CacheTestsOnRedis(unit.CacheTests):
    """The unit suite's CacheTests against Redis.

    Two bodies are replaced rather than inherited: they expire a key by moving
    the test's clock, and a real server does not read it.  The replacements
    test the same property against a TTL the server can actually reach, which
    makes them stronger, not weaker.
    """

    def setUp(self):
        self.clock = unit.Clock()
        self.campaign = campaign()
        self.redis = LiveRedis(liveCache().client)
        self.cache = Cache(self.redis, campaign=self.campaign, clock=self.clock)

    def test_a_stalled_holder_cannot_delete_the_next_holders_lock(self):
        stalled = self.cache.lock("merge", ttl=1)
        self.assertTrue(stalled.acquire())
        time.sleep(1.3)                                     # the server expires it
        nextHolder = self.cache.lock("merge", ttl=30)
        self.assertTrue(nextHolder.acquire())
        self.assertFalse(stalled.release())                 # token mismatch
        self.assertIsNotNone(self.redis.get(self.cache.key("lock", "merge")))
        self.assertTrue(nextHolder.release())

    def test_point_hints_are_hints(self):
        self.assertFalse(self.cache.seenPoint("abc", ttl=1))    # first sighting
        self.assertTrue(self.cache.seenPoint("abc", ttl=1))     # offered again
        time.sleep(1.3)                                         # the window passes
        self.assertFalse(self.cache.seenPoint("abc", ttl=1))
        self.redis.fail = TimeoutError("timed out")
        self.assertIsNone(self.cache.seenPoint("abc"))          # unknown, never False


class CoordinatorTestsOnRealStores(unit.CoordinatorTests):
    """The unit suite's CoordinatorTests with Postgres, Redis and S3 under it."""

    def setUp(self):
        self.clock = unit.Clock()
        self.campaign = campaign()
        self.db = DB
        self.store = liveStore(prefix="%s/%s" % (LOCAL_PREFIX, self.campaign))
        self.redis = LiveRedis(liveCache().client)
        self.cache = Cache(self.redis, campaign=self.campaign, clock=self.clock)
        self.coord = Coordinator(self.db, self.store, self.cache, campaign=self.campaign,
                                 owner="worker-a", clock=self.clock)


# ---------------------------------------------------------------------------
# What only Postgres can answer.
# ---------------------------------------------------------------------------
class PostgresTests(unittest.TestCase):
    def test_two_hosts_migrating_at_once_apply_each_migration_once(self):
        """A rollout runs `migrate` on every host at the same moment.

        Without the advisory lock two of them interleave `CREATE INDEX` and
        both write the marker; with it, one applies the list and the rest find
        nothing to do.  Run on a database of its own, because the point is the
        migration from empty.
        """
        import psycopg
        from psycopg.conninfo import conninfo_to_dict, make_conninfo

        def url(dbname):
            info = conninfo_to_dict(CONFIG.databaseUrl)
            info["dbname"] = dbname
            return make_conninfo(**info)

        name = "it_migrate_%s" % uuid.uuid4().hex[:10]
        admin = psycopg.connect(url("postgres"), autocommit=True)
        try:
            admin.execute('CREATE DATABASE "%s"' % name)
        finally:
            admin.close()

        target = url(name)
        handles = [Database(lambda t=target: psycopg.connect(t, autocommit=False))
                   for _ in range(6)]
        applied, errors = [], []
        lock = threading.Lock()

        def run(db):
            try:
                got = migrate(db)
                with lock:
                    applied.extend(got)
            except Exception as exc:                        # pragma: no cover
                with lock:
                    errors.append(exc)

        threads = [threading.Thread(target=run, args=(h,)) for h in handles]
        for t in threads:
            t.start()
        for t in threads:
            t.join()
        try:
            self.assertEqual(errors, [])
            # Exactly once each: the lock serialises the appliers, and the
            # marker row keeps the rest from repeating the work.
            self.assertEqual(sorted(applied), sorted(ident for ident, _ in MIGRATIONS))
            self.assertEqual(pending(handles[0]), [])
        finally:
            for h in handles:
                h.close()
            admin = psycopg.connect(url("postgres"), autocommit=True)
            try:
                admin.execute('DROP DATABASE IF EXISTS "%s"' % name)
            finally:
                admin.close()

    def test_a_terminated_backend_is_retried_on_a_fresh_connection(self):
        """The failover path, with the failover actually happening.

        `isTransient` is unit-tested against exception shapes; this kills the
        connection under the handle with `pg_terminate_backend` and requires
        the next call to succeed on a new one.  That is what an RDS failover
        looks like from a worker's side.
        """
        import psycopg

        db = liveDatabase()
        try:
            before = db.fetchOne("SELECT pg_backend_pid()")[0]
            killer = psycopg.connect(CONFIG.databaseUrl, autocommit=True)
            try:
                killer.execute("SELECT pg_terminate_backend(%s)", (before,))
            finally:
                killer.close()
            self.assertEqual(db.fetchOne("SELECT 1")[0], 1)
            after = db.fetchOne("SELECT pg_backend_pid()")[0]
            self.assertNotEqual(before, after)
        finally:
            db.close()

    def test_every_connection_carries_the_statement_timeout(self):
        """A stuck claim holds a GPU idle, so the timeout is not advisory.

        SQLite has no such setting, so this is the only place the `SET` in
        `Database.fromEnv`'s connect hook is known to have run.
        """
        db = liveDatabase()
        try:
            self.assertEqual(db.fetchOne("SHOW statement_timeout")[0], "15s")
            self.assertEqual(db.fetchOne("SHOW idle_in_transaction_session_timeout")[0], "30s")
        finally:
            db.close()

    def test_the_counters_hold_a_real_campaigns_iteration_count(self):
        """2**62 iterations is inside this campaign's range and outside int4.

        On SQLite every integer column takes it, so a schema that said
        `integer` would pass the unit suite and overflow in the fleet's first
        week.
        """
        slots = PostgresSlots(DB, campaign=campaign(), clock=unit.Clock())
        slot, fence = slots.claimLease("worker-a", unit.gpuInfo())
        big = 2 ** 62
        self.assertTrue(slots.heartbeat(slot, "worker-a", {"iters": big}, fence=fence))
        self.assertEqual(slots.get(slot)["iters"], big)


# ---------------------------------------------------------------------------
# What only Redis can answer.
# ---------------------------------------------------------------------------
class RedisTests(unittest.TestCase):
    def setUp(self):
        self.clock = unit.Clock()
        self.campaign = campaign()
        self.client = liveCache().client
        self.cache = Cache(self.client, campaign=self.campaign, clock=self.clock)

    def test_the_unlock_script_runs_on_the_server(self):
        """Compare-and-delete is a Lua script, and here Redis is the one running it."""
        key = self.cache.key("lock", "merge")
        holder = self.cache.lock("merge", ttl=60)
        self.assertTrue(holder.acquire())

        impostor = self.cache.lock("merge", ttl=60)
        impostor.held = True                                # as if it believed it held
        self.assertFalse(impostor.release())                # the script says otherwise
        self.assertIsNotNone(self.client.get(key))

        self.assertTrue(holder.extend(ttl=45))
        self.assertGreater(self.client.pttl(key), 30_000)
        self.assertTrue(holder.release())
        self.assertIsNone(self.client.get(key))

    def test_a_heartbeat_writes_both_halves_of_the_fleet_view(self):
        """The pipeline is four commands, and all four have to land.

        A fake that implements `pipeline` as a list of calls cannot show that
        the hash, the sorted set and both expiries survived one round trip to
        a real server.
        """
        self.cache.beat("w1", {"slot": 3})
        blob = self.client.hget(self.cache.key("fleet"), "w1")
        self.assertEqual(json.loads(blob)["slot"], 3)
        self.assertEqual(int(self.client.zscore(self.cache.key("fleet", "seen"), "w1")),
                         int(self.clock.now))
        self.assertGreater(self.client.ttl(self.cache.key("fleet")), 0)
        self.assertGreater(self.client.ttl(self.cache.key("fleet", "seen")), 0)

    def test_a_cluster_that_is_not_there_degrades_instead_of_raising(self):
        """Not a fake failure: a real client pointed at a port nobody answers.

        Every call must come back with the "cannot say" answer its caller
        knows how to handle, the lock must be granted anyway, and after three
        failures the cache must stop being called at all.
        """
        dead = deadCache(self.campaign, self.clock)
        self.assertIsNone(dead.getJson("anything"))
        self.assertIsNone(dead.fleet())                     # None, not []
        self.assertIsNone(dead.seenPoint("abc"))            # None, never False
        self.assertTrue(dead.stats()["cooling"])
        self.assertTrue(dead.lock("merge").acquire())       # the work goes on
        self.assertFalse(dead.setJson("anything", {"a": 1}))

    def test_an_unreachable_cluster_costs_one_timeout_and_not_ten(self):
        """The cache's promise is that it can never be slower than the database.

        redis-py 6 and later retry a connection failure ten times with backoff
        unless told not to, which measured 4.5 seconds per call here -- nine
        times the timeout this package sets, on the claim path, before the
        cooldown has even seen three failures.  The client `Cache.fromEnv`
        builds has that turned off; this is what says so.
        """
        dead = deadCache(self.campaign, self.clock)
        started = time.time()
        for _ in range(3):
            self.assertIsNone(dead.getJson("anything"))
        elapsed = time.time() - started
        # Three calls against a refused port. With the default retry policy
        # this is upwards of thirteen seconds.
        self.assertLess(elapsed, 2.0, "a dead cluster cost %.1fs for three calls" % elapsed)


# ---------------------------------------------------------------------------
# What only an object store behind HTTP can answer.
# ---------------------------------------------------------------------------
class ObjectStoreTests(unittest.TestCase):
    def setUp(self):
        self.prefix = "%s/%s" % (LOCAL_PREFIX, campaign())
        self.store = liveStore(prefix=self.prefix)

    def test_bytes_and_their_digest_survive_the_round_trip(self):
        data = b"a corpus object, thirty-two byt"
        put = self.store.putBytes("dp/slot-00000/a-0.bin", data)
        self.assertFalse(put["existed"])
        head = self.store.head("dp/slot-00000/a-0.bin")
        # The digest is user metadata on the object, which is what `putBytes`
        # compares a re-upload against.  A store that dropped it would make
        # every retry look like a conflict.
        self.assertEqual(head["sha256"], put["sha256"])
        self.assertEqual(head["bytes"], len(data))
        self.assertEqual(self.store.getBytes("dp/slot-00000/a-0.bin"), data)

    def test_a_retry_is_not_a_conflict_and_different_bytes_are(self):
        self.store.putBytes("dp/a.bin", b"x" * 64)
        again = self.store.putBytes("dp/a.bin", b"x" * 64)
        self.assertTrue(again["existed"])
        with self.assertRaises(ObjectExists):
            self.store.putBytes("dp/a.bin", b"y" * 64)
        self.assertEqual(self.store.getBytes("dp/a.bin"), b"x" * 64)

    def test_a_missing_object_is_none_rather_than_an_exception(self):
        self.assertIsNone(self.store.head("dp/never-written.bin"))

    def test_the_listing_pages_past_the_first_thousand_keys(self):
        """`list_objects_v2` returns 1,000 keys and a token, not the corpus.

        The pagination in `store.py` is written by hand; against a directory
        there is nothing to paginate, so this is the only place the loop that
        follows the continuation token is exercised.  A corpus silently
        truncated at 1,000 objects is a merge that quietly misses points.
        """
        count = 1002
        with ThreadPoolExecutor(max_workers=16) as pool:
            list(pool.map(lambda i: self.store.putBytes("many/%04d.bin" % i, b"\x01"),
                          range(count)))
        listed = self.store.list("many/")
        self.assertEqual(len(listed), count)
        self.assertEqual(len({item["key"] for item in listed}), count)

    def test_a_listing_is_scoped_to_its_prefix(self):
        self.store.putBytes("dp/one.bin", b"1")
        self.store.putBytes("ck/two.bin", b"2")
        keys = [item["key"] for item in self.store.list("dp/")]
        self.assertEqual(keys, ["%s/dp/one.bin" % self.prefix])


# ---------------------------------------------------------------------------
# All three at once.
# ---------------------------------------------------------------------------
class CrossStoreTests(unittest.TestCase):
    def setUp(self):
        self.clock = unit.Clock()
        self.campaign = campaign()
        self.store = liveStore(prefix="%s/%s" % (LOCAL_PREFIX, self.campaign))
        self.cache = Cache(liveCache().client, campaign=self.campaign, clock=self.clock)
        self.coord = Coordinator(DB, self.store, self.cache, campaign=self.campaign,
                                 owner="worker-a", clock=self.clock)

    def test_an_upload_that_fails_leaves_no_ledger_row(self):
        """The ordering invariant, from the failing side.

        The happy path is covered by the unit suite and re-run above; what
        matters here is that a store that refuses the bytes cannot leave a row
        claiming points that are not in the bucket.
        """
        lease = self.coord.acquire(unit.gpuInfo())

        class Refusing:
            def __init__(self, inner):
                self.inner = inner

            def __getattr__(self, name):
                return getattr(self.inner, name)

            def putBytes(self, *args, **kw):
                raise RuntimeError("the bucket said no")

        self.coord.store = Refusing(self.store)
        with self.assertRaises(RuntimeError):
            self.coord.publishObject(lease, "dp/slot-00000/a-0.bin", b"x" * 64)
        self.assertEqual(self.coord.objects(slot=lease.slot), [])
        self.assertIsNone(self.store.head("dp/slot-00000/a-0.bin"))

    def test_a_walk_runs_end_to_end_across_the_three_stores(self):
        """One worker's life, with nothing faked underneath it.

        Claim, upload, checkpoint, heartbeat, lose the lease to a second
        worker, and have every stale write refused -- against Postgres, Redis
        and S3 at once.  This is the shape of what a supervisor does between
        two spot interruptions.
        """
        lease = self.coord.acquire(unit.gpuInfo())
        for i in range(3):
            self.coord.publishObject(lease, "dp/slot-%05d/a-%d.bin" % (lease.slot, i),
                                     b"p" * 64, streamId="a", offset=i * 64)
        self.coord.setCheckpoint(lease, "ck/slot-%05d-1.ck" % lease.slot, iters=10 ** 9)
        self.assertTrue(self.coord.heartbeat(lease, {"iters": 10 ** 9, "points": 6}))

        fleet = self.coord.fleet()
        self.assertEqual(fleet["source"], "cache")
        self.assertEqual([w["owner"] for w in fleet["workers"]], ["worker-a"])

        status = self.coord.status()
        self.assertEqual(status["corpus"]["objects"], 3)
        self.assertEqual(status["corpus"]["records"], 6)

        # The spot interruption: A freezes past its lease, B takes the slot.
        self.clock.advance(181)
        b = Coordinator(DB, self.store, self.cache, campaign=self.campaign,
                        owner="worker-b", clock=self.clock)
        resumed = b.acquire(unit.gpuInfo(instance="i-2"))
        self.assertEqual(resumed.slot, lease.slot)
        self.assertGreater(resumed.fence, lease.fence)
        self.assertEqual(b.checkpoint(resumed.slot)["iters"], 10 ** 9)

        # A wakes up. Postgres refuses it, whatever A believes it still holds.
        with self.assertRaises(Exception):
            self.coord.setCheckpoint(lease, "ck/stale.ck", iters=10 ** 10)
        self.assertEqual(b.checkpoint(resumed.slot)["iters"], 10 ** 9)
        self.assertFalse(self.coord.heartbeat(lease, {"iters": 10 ** 10}))


class CliTests(unittest.TestCase):
    def test_doctor_reports_all_three_stores_from_a_clean_process(self):
        """`doctor` in a subprocess is the whole stack with nothing injected.

        It resolves its own config, imports boto3 and redis for real, and has
        to reach all three services; the unit suite runs it with no cluster
        and no bucket, which is the one path that cannot fail this way.
        """
        here = os.path.dirname(os.path.abspath(__file__))
        env = dict(os.environ, RHO_CAMPAIGN=campaign())
        env.pop("ECC_LOCAL_STORE", None)                    # the real stores, not a directory
        run = subprocess.run([sys.executable, "-m", "controlplane", "doctor"],
                             capture_output=True, text=True, cwd=here, env=env)
        self.assertEqual(run.returncode, 0, run.stderr)
        report = json.loads(run.stdout)
        self.assertTrue(report["ok"], report)
        self.assertTrue(report["database"]["ok"], report["database"])
        self.assertEqual(report["database"]["dialect"], "postgres")
        self.assertEqual(report["database"]["pending"], [])
        self.assertTrue(report["cache"]["ok"], report["cache"])
        self.assertTrue(report["corpus"]["ok"], report["corpus"])
        self.assertEqual(report["corpus"]["bucket"], CONFIG.bucket)


if __name__ == "__main__":
    unittest.main(verbosity=2)
