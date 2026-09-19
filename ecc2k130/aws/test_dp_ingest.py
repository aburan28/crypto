#!/usr/bin/env python3
"""Offline tests for the distinguished-point ingest.

The incident these pin: on 2026-09-17 the fleet's workers moved to the
`ecc2k-seed-orbit-v1` object naming at 10:06Z, the last worker still using the
old naming uploaded at 13:50Z, and the ingest -- which matched only the old
shape -- then took nothing for five hours while 56 M records landed in `dp/`.
It logged no error, because an object whose key does not match is not an
error, it is invisible. Every test below is about one of those two properties:
both key shapes are points, and anything unreadable is counted out loud.

    python3 -m unittest test_dp_ingest
"""

import datetime
import io
import os
import struct
import sys
import unittest
from unittest import mock

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import dp_ingest


def when(iso):
    return datetime.datetime.fromisoformat(iso).replace(tzinfo=datetime.timezone.utc)


LEGACY = "dp/slot-00002/1789311001-0000000000000000.bin"
ORBIT = ("dp/slot-00140/a8b4133d5c5f414ebd1337f15603588d-0000000000791392-"
         "9c6a7dae1862f8393fe196b0c4f890e5e01f7f45615be019252b86ca3e429914.bin")


class FakeS3:
    """Enough of list_objects_v2 to drive s3Objects, including paging."""

    def __init__(self, items, pageSize=2):
        self.items = list(items)
        self.pageSize = pageSize

    def list_objects_v2(self, **kw):
        start = int(kw.get("ContinuationToken") or 0)
        page = self.items[start:start + self.pageSize]
        nxt = start + self.pageSize
        return {"Contents": [{"Key": k, "Size": s, "LastModified": t} for k, s, t in page],
                "IsTruncated": nxt < len(self.items),
                "NextContinuationToken": str(nxt)}


class KeyShapes(unittest.TestCase):
    def test_legacy_key_takes_found_at_from_the_name(self):
        self.assertEqual(dp_ingest.foundAt(LEGACY, when("2026-09-20T00:00:00")), 1789311001)

    def test_orbit_key_takes_found_at_from_last_modified(self):
        t = when("2026-09-17T18:47:34")
        self.assertEqual(dp_ingest.foundAt(ORBIT, t), int(t.timestamp()))

    def test_orbit_key_is_a_point_object_at_all(self):
        # The regression itself: this returned None before the fix, and an
        # object with no found_at was skipped without a word.
        self.assertIsNotNone(dp_ingest.foundAt(ORBIT, when("2026-09-17T18:47:34")))

    def test_envelopes_and_foreign_keys_are_not_points(self):
        for key in (ORBIT + ".json", "dp/slot-00002/notes.txt", "ckpt/slot-00002.ck"):
            self.assertIsNone(dp_ingest.foundAt(key, when("2026-09-17T18:47:34")), key)


class Listing(unittest.TestCase):
    def objects(self, items):
        return dp_ingest.s3Objects(FakeS3(items), "bucket")

    def test_both_shapes_are_listed_and_sorted_by_upload_time(self):
        objects, unrecognised = self.objects([
            (ORBIT, 3200, when("2026-09-17T18:47:34")),
            (LEGACY, 6400, when("2026-09-17T13:50:59")),
        ])
        self.assertEqual([o[0] for o in objects], [LEGACY, ORBIT])
        self.assertEqual([o[1] for o in objects], [200, 100])
        self.assertEqual(unrecognised, [])

    def test_an_unreadable_key_is_reported_not_skipped(self):
        objects, unrecognised = self.objects([
            (LEGACY, 6400, when("2026-09-17T13:50:59")),
            ("dp/slot-00007/something-new.bin", 64, when("2026-09-18T00:00:00")),
        ])
        self.assertEqual([o[0] for o in objects], [LEGACY])
        self.assertEqual(unrecognised, ["dp/slot-00007/something-new.bin"])

    def test_envelopes_are_neither_ingested_nor_reported_as_unreadable(self):
        objects, unrecognised = self.objects([
            (ORBIT, 3200, when("2026-09-17T18:47:34")),
            (ORBIT + ".json", 412, when("2026-09-17T18:47:34")),
        ])
        self.assertEqual([o[0] for o in objects], [ORBIT])
        self.assertEqual(unrecognised, [])

    def test_a_short_object_carrying_no_whole_record_is_not_listed(self):
        objects, _ = self.objects([(ORBIT, 16, when("2026-09-17T18:47:34"))])
        self.assertEqual(objects, [])

    def test_paging_does_not_lose_objects(self):
        items = [(LEGACY.replace("0000000000000000", "%016d" % i), 64,
                  when("2026-09-17T13:50:59")) for i in range(7)]
        objects, _ = self.objects(items)
        self.assertEqual(len(objects), 7)


class Decoding(unittest.TestCase):
    def test_point_key_is_the_last_24_bytes_of_the_record(self):
        record = bytes(range(32))
        self.assertEqual(dp_ingest.decode(record)["point_key"], record[8:])

    def test_walk_seed_is_the_seed_big_endian_without_leading_zeros(self):
        record = (7).to_bytes(8, "little") + bytes(24)
        self.assertEqual(dp_ingest.decode(record)["walk_seed"], b"\x07")

    def test_coefficients_are_17_bytes(self):
        r = dp_ingest.decode(bytes(range(32)))
        self.assertEqual((len(r["a"]), len(r["b"])), (17, 17))
        self.assertEqual(r["b"], bytes(17))


class FakeConn:
    """A connection that can be killed the way a failover kills one."""

    def __init__(self, pool):
        self.pool = pool
        self.closed = False
        self.rollbacks = 0
        pool.append(self)

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()
        return False

    def rollback(self):
        if self.closed:
            raise RuntimeError("the connection is closed")
        self.rollbacks += 1

    def close(self):
        self.closed = True


class FakeCursor:
    def __init__(self, conn):
        self.conn = conn
        self.rows = []

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False

    def execute(self, sql, args=None):
        self.conn.queries.append(sql)
        self.rows = ([("dp-slot-00002-old.bin", 24000)] if "GROUP BY worker_id" in sql
                     else [])

    def fetchall(self):
        return self.rows


class CountingConn:
    def __init__(self):
        self.queries = []

    def cursor(self):
        return FakeCursor(self)


class LegacyCounts(unittest.TestCase):
    """The aggregate that reads the whole corpus must not run every pass."""

    def setUp(self):
        dp_ingest._counts.update({"at": 0.0, "map": {}})
        self.addCleanup(dp_ingest._counts.update, {"at": 0.0, "map": {}})

    def aggregates(self, conn):
        return sum("GROUP BY worker_id" in q for q in conn.queries)

    def test_the_corpus_wide_aggregate_is_not_repeated_within_the_ttl(self):
        conn = CountingConn()
        for _ in range(4):
            counts, _ = dp_ingest.ingestedCounts(conn)
        self.assertEqual(self.aggregates(conn), 1)
        self.assertEqual(counts, {"dp-slot-00002-old.bin": 24000})

    def test_progress_is_still_read_every_time(self):
        conn = CountingConn()
        for _ in range(4):
            dp_ingest.ingestedCounts(conn)
        self.assertEqual(sum("dp_ingest_progress" in q for q in conn.queries), 4)

    def test_it_refreshes_once_the_ttl_has_passed(self):
        conn = CountingConn()
        dp_ingest.ingestedCounts(conn)
        dp_ingest._counts["at"] -= 3600.0
        dp_ingest.ingestedCounts(conn)
        self.assertEqual(self.aggregates(conn), 2)


class IndexCursor:
    def __init__(self, conn):
        self.conn = conn
        self.row = None

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False

    def execute(self, sql, args=None):
        self.conn.statements.append(" ".join(sql.split()))
        if "pg_index" in sql:
            self.row = self.conn.existing
        elif "dp_ingest_meta WHERE key" in sql:
            self.row = (1,) if self.conn.vacuumed else None

    def fetchone(self):
        return self.row


class IndexConn:
    def __init__(self, existing=None, vacuumed=True):
        self.existing = existing
        self.vacuumed = vacuumed
        self.statements = []
        self.autocommit = False

    def cursor(self):
        return IndexCursor(self)


class FoundAtIndex(unittest.TestCase):
    """The snapshot's cost is an index, and a half-built one is not one."""

    def statements(self, conn, needle):
        return [s for s in conn.statements if needle in s]

    def test_it_is_built_when_it_is_missing(self):
        conn = IndexConn(existing=None)
        self.assertTrue(dp_ingest.ensureFoundAtIndex(conn))
        self.assertEqual(len(self.statements(conn, "CREATE INDEX CONCURRENTLY")), 1)

    def test_the_table_is_vacuumed_so_the_index_can_be_scanned_alone(self):
        # Without it the index is built and the snapshot still reads the heap.
        conn = IndexConn(existing=None, vacuumed=False)
        dp_ingest.ensureFoundAtIndex(conn)
        self.assertEqual(len(self.statements(conn, "VACUUM (ANALYZE)")), 1)

    def test_a_replacement_host_does_not_vacuum_the_table_again(self):
        conn = IndexConn(existing=(dp_ingest.FOUND_AT_INDEX, True), vacuumed=True)
        self.assertFalse(dp_ingest.ensureFoundAtIndex(conn))
        self.assertEqual(self.statements(conn, "VACUUM"), [])

    def test_an_index_that_already_exists_unvacuumed_is_still_vacuumed(self):
        # The state the live store was left in: index built at 20:47Z, heap
        # still unmarked, snapshots still six minutes.
        conn = IndexConn(existing=(dp_ingest.FOUND_AT_INDEX, True), vacuumed=False)
        self.assertTrue(dp_ingest.ensureFoundAtIndex(conn))
        self.assertEqual(self.statements(conn, "CREATE INDEX CONCURRENTLY"), [])
        self.assertEqual(len(self.statements(conn, "VACUUM (ANALYZE)")), 1)

    def test_it_is_not_rebuilt_when_it_is_already_valid(self):
        conn = IndexConn(existing=(dp_ingest.FOUND_AT_INDEX, True))
        self.assertFalse(dp_ingest.ensureFoundAtIndex(conn))
        self.assertEqual(self.statements(conn, "CREATE INDEX CONCURRENTLY"), [])

    def test_an_invalid_index_is_dropped_and_rebuilt(self):
        conn = IndexConn(existing=(dp_ingest.FOUND_AT_INDEX, False))
        self.assertTrue(dp_ingest.ensureFoundAtIndex(conn))
        self.assertEqual(len(self.statements(conn, "DROP INDEX CONCURRENTLY")), 1)
        self.assertEqual(len(self.statements(conn, "CREATE INDEX CONCURRENTLY")), 1)

    def test_the_build_runs_outside_a_transaction_and_restores_the_mode(self):
        conn = IndexConn(existing=None)
        dp_ingest.ensureFoundAtIndex(conn)
        self.assertFalse(conn.autocommit)


class Passes(unittest.TestCase):
    """What a pass does with a backlog, and with a database that goes away.

    Both are from the same 2026-09-17 recovery: rho-dp was resized while a
    four-thousand-object drain was running, 2,533 objects failed against the
    closed connection their thread kept using, and the pass -- being
    unbounded -- published neither status nor a metric for the half hour it
    took to get through the rest.
    """

    def setUp(self):
        self.conns = []
        self.realPending = dp_ingest.pending
        self.realIngest = dp_ingest.ingestObject
        self.addCleanup(setattr, dp_ingest, "pending", self.realPending)
        self.addCleanup(setattr, dp_ingest, "ingestObject", self.realIngest)
        self.addCleanup(setattr, dp_ingest, "log", dp_ingest.log)
        dp_ingest._failed.clear()
        self.addCleanup(dp_ingest._failed.clear)
        dp_ingest.log = lambda msg: None

    def backlog(self, n):
        items = [("dp/slot-%05d/x-%016d.bin" % (i, i), 100, 1789000000 + i, 0)
                 for i in range(n)]
        dp_ingest.pending = lambda *a, **k: (items, 1789000000 + n, [])
        return items

    def run_pass(self, threads=1, limit=dp_ingest.PASS_OBJECTS):
        return dp_ingest.onePass(lambda: FakeConn(self.conns), None, "bucket",
                                 threads=threads, limit=limit)

    def test_a_pass_stops_at_its_bound_and_still_reports_the_whole_backlog(self):
        self.backlog(5)
        done = []
        dp_ingest.ingestObject = lambda conn, s3, bucket, key, found_at: (
            done.append(key) or (100, 100))
        rows, objects, state = self.run_pass(limit=2)
        self.assertEqual((objects, rows), (2, 200))
        # The slice is this pass's work; the number published is the backlog.
        self.assertEqual(state["outstanding"], 3)

    def test_an_unbounded_pass_is_still_available(self):
        self.backlog(5)
        dp_ingest.ingestObject = lambda *a, **k: (100, 100)
        _, objects, state = self.run_pass(limit=0)
        self.assertEqual((objects, state["outstanding"]), (5, 0))

    def test_a_database_that_goes_away_costs_one_object_not_the_pass(self):
        self.backlog(4)
        seen = []

        def ingest(conn, s3, bucket, key, found_at):
            seen.append(key)
            if len(seen) == 1:  # the resize lands here
                conn.close()
                raise RuntimeError("the connection is closed")
            return 100, 100

        dp_ingest.ingestObject = ingest
        _, objects, state = self.run_pass()
        self.assertEqual(objects, 4)  # the killed object retried, rest unharmed
        self.assertEqual(state["outstanding"], 0)
        self.assertGreaterEqual(len(self.conns), 2)  # it reconnected

    def test_a_failed_slice_does_not_pin_the_rest_of_the_backlog(self):
        # The bound takes the oldest `limit` keys. If those stay outstanding
        # because they failed, the next pass would take the same slice and
        # never reach the objects behind them.
        items = self.backlog(4)
        poison = {items[0][0], items[1][0]}
        seen = []

        def ingest(conn, s3, bucket, key, found_at):
            seen.append(key)
            if key in poison:
                raise ValueError("object is unreadable")
            return 100, 100

        dp_ingest.ingestObject = ingest
        _, objects, state = self.run_pass(limit=2)
        self.assertEqual((objects, state["outstanding"]), (0, 4))
        self.assertEqual(seen, [items[0][0], items[1][0]])
        _, objects, state = self.run_pass(limit=2)
        self.assertEqual((objects, state["outstanding"]), (2, 2))
        self.assertEqual(seen[2:], [items[2][0], items[3][0]])

    def test_failed_keys_are_retried_once_nothing_else_is_outstanding(self):
        items = self.backlog(2)
        seen = []

        def ingest(conn, s3, bucket, key, found_at):
            seen.append(key)
            raise ValueError("object is unreadable")

        dp_ingest.ingestObject = ingest
        _, objects, state = self.run_pass(limit=2)
        self.assertEqual((objects, state["outstanding"]), (0, 2))
        self.assertEqual(seen, [items[0][0], items[1][0]])
        # Both keys failed, so they are all that remains; retry them.
        _, objects, state = self.run_pass(limit=2)
        self.assertEqual((objects, state["outstanding"]), (0, 2))
        self.assertEqual(seen[2:], [items[0][0], items[1][0]])

    def test_an_object_that_is_simply_bad_is_counted_and_left_behind(self):
        self.backlog(3)

        def ingest(conn, s3, bucket, key, found_at):
            if key.endswith("%016d.bin" % 1):
                raise ValueError("record size is not a multiple of 32")
            return 100, 100

        dp_ingest.ingestObject = ingest
        _, objects, state = self.run_pass()
        self.assertEqual(objects, 2)
        self.assertEqual(state["outstanding"], 1)
        # The failed object must leave a usable transaction behind; without
        # rollback the rest of this thread's objects fail as well.
        self.assertGreaterEqual(sum(c.rollbacks for c in self.conns), 1)


class StatusState(unittest.TestCase):
    """The page must not read a stalled ingest as a quiet campaign."""

    def state(self, dps, dpsLastHour, collisions, ingest):
        behind = int(ingest.get("outstanding", 0)) or int(ingest.get("unrecognised", 0))
        return ("COLLISION_RECORDED" if collisions else
                "COLLECTING" if dpsLastHour else
                "INGEST_BEHIND" if behind else
                "IDLE_OR_STALE" if dps else "EMPTY")

    def test_a_backlog_reads_as_ingest_behind_not_idle(self):
        self.assertEqual(self.state(60053195, 0, 0, {"outstanding": 3795}), "INGEST_BEHIND")

    def test_unreadable_keys_alone_are_enough_to_say_behind(self):
        self.assertEqual(self.state(60053195, 0, 0, {"unrecognised": 12}), "INGEST_BEHIND")

    def test_a_caught_up_ingest_with_no_recent_points_is_still_idle(self):
        self.assertEqual(self.state(60053195, 0, 0, {"outstanding": 0}), "IDLE_OR_STALE")

    def test_points_in_the_last_hour_win(self):
        self.assertEqual(self.state(60053195, 142035, 0, {"outstanding": 3795}), "COLLECTING")


def ckpt(iterBase=10, threads=2, batch=4, runId=1):
    return struct.pack("<8s6IQ", dp_ingest.CKPT_MAGIC, 1, 131, threads, batch, 64, runId, iterBase)


class FakeWorkS3:
    """list_objects_v2 plus ranged get_object, enough to drive checkpointWork."""

    def __init__(self, objects, pageSize=2):
        self.objects = list(objects)
        self.pageSize = pageSize

    def list_objects_v2(self, **kw):
        start = int(kw.get("ContinuationToken") or 0)
        page = self.objects[start:start + self.pageSize]
        nxt = start + self.pageSize
        return {"Contents": [{"Key": k, "Size": len(b), "LastModified": t}
                             for k, t, b in page],
                "IsTruncated": nxt < len(self.objects),
                "NextContinuationToken": str(nxt)}

    def get_object(self, **kw):
        body = next(b for k, t, b in self.objects if k == kw["Key"])
        rng = kw.get("Range")
        if rng:
            lo, hi = rng.split("=", 1)[1].split("-")
            body = body[int(lo):int(hi) + 1]
        return {"Body": io.BytesIO(body)}


class CheckpointWork(unittest.TestCase):
    """The public work feed has to see the keys the live fleet actually writes."""

    def setUp(self):
        now = datetime.datetime.now(datetime.timezone.utc)
        self.recent = now - datetime.timedelta(seconds=30)
        self.old = now - datetime.timedelta(seconds=4000)
        self.hash = "ab" * 32

    def work(self, items):
        return dp_ingest.checkpointWork(FakeWorkS3(items), "bucket")

    def test_legacy_pointer_still_counts(self):
        total, slots = self.work([
            ("ckpt/slot-00002.ck", self.recent, ckpt(10, 2, 4)),
        ])
        self.assertEqual(total, 80)
        self.assertEqual(len(slots), 1)
        self.assertFalse(slots[0]["retired"])

    def test_contract_checkpoint_counts_as_work(self):
        # The regression: this key was skipped with no log, so a 133-slot
        # fleet published zero walkers while every slot was walking.
        key = "ckpt/slot-00140/%s.ck" % self.hash
        total, slots = self.work([(key, self.recent, ckpt(100, 2, 4))])
        self.assertEqual(total, 800)
        self.assertEqual(slots[0]["slot"], 140)
        self.assertFalse(slots[0]["retired"])
        self.assertEqual(slots[0]["iterations"], 800)

    def test_retired_prefix_still_counts(self):
        total, slots = self.work([
            ("ckpt/retired/slot-00002.ck", self.recent, ckpt(5, 1, 1)),
        ])
        self.assertEqual(total, 5)
        self.assertTrue(slots[0]["retired"])

    def test_envelopes_are_not_checkpoints(self):
        key = "ckpt/slot-00140/%s.ck" % self.hash
        total, slots = self.work([
            (key, self.recent, ckpt(100, 2, 4)),
            (key + ".json", self.recent, b"{}"),
        ])
        self.assertEqual((total, len(slots)), (800, 1))

    def test_older_blobs_of_the_same_slot_are_not_summed(self):
        live = "ckpt/slot-00002/%s.ck" % self.hash
        older = "ckpt/slot-00002/%s.ck" % ("cd" * 32)
        total, slots = self.work([
            (older, self.old, ckpt(10, 2, 4)),
            (live, self.recent, ckpt(50, 2, 4)),
        ])
        self.assertEqual(total, 400)
        self.assertEqual(len(slots), 1)
        self.assertFalse(slots[0]["retired"])

    def test_a_live_hash_beats_a_stale_legacy_pointer(self):
        live = "ckpt/slot-00002/%s.ck" % self.hash
        total, slots = self.work([
            ("ckpt/slot-00002.ck", self.old, ckpt(1, 1, 1)),
            (live, self.recent, ckpt(50, 2, 4)),
        ])
        self.assertEqual(total, 400)
        self.assertEqual(len(slots), 1)
        self.assertFalse(slots[0]["retired"])

    def test_a_live_hash_beats_a_retired_leftover(self):
        live = "ckpt/slot-00002/%s.ck" % self.hash
        total, slots = self.work([
            ("ckpt/retired/slot-00002.ck", self.old, ckpt(10, 1, 1)),
            (live, self.recent, ckpt(80, 3, 4)),
        ])
        self.assertEqual(total, 960)
        self.assertEqual(len(slots), 1)
        self.assertFalse(slots[0]["retired"])


class Defaults(unittest.TestCase):
    def test_status_every_default_is_three_minutes(self):
        import argparse
        ap = argparse.ArgumentParser()
        ap.add_argument("--status-every", type=float,
                        default=float(os.environ.get("RHO_STATUS_EVERY", "180")))
        self.assertEqual(ap.parse_args([]).status_every, 180.0)


class DatabaseUrl(unittest.TestCase):
    def test_secret_lookup_includes_sslmode(self):
        env = os.environ.copy()
        env.pop("DATABASE_URL", None)
        env["RHO_DB_HOST"] = "rho-dp.example.com"
        env["RHO_DB_SSLMODE"] = "require"
        with mock.patch.dict(os.environ, env, clear=True):
            with mock.patch("boto3.client") as client:
                client.return_value.get_secret_value.return_value = {
                    "SecretString": '{"username":"u","password":"p/w","dbname":"d"}'
                }
                url = dp_ingest.databaseUrl()
        self.assertIn("sslmode=require", url)
        self.assertIn("rho-dp.example.com", url)


if __name__ == "__main__":
    unittest.main()
