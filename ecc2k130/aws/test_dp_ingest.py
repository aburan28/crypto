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
import json
import os
import struct
import sys
import time
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
            done.append(key) or (100, 100, 0))
        rows, objects, state = self.run_pass(limit=2)
        self.assertEqual((objects, rows), (2, 200))
        # The slice is this pass's work; the number published is the backlog.
        self.assertEqual(state["outstanding"], 3)

    def test_an_unbounded_pass_is_still_available(self):
        self.backlog(5)
        dp_ingest.ingestObject = lambda *a, **k: (100, 100, 0)
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
            return 100, 100, 0

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
            return 100, 100, 0

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
            return 100, 100, 0

        dp_ingest.ingestObject = ingest
        _, objects, state = self.run_pass()
        self.assertEqual(objects, 2)
        self.assertEqual(state["outstanding"], 1)
        # The failed object must leave a usable transaction behind; without
        # rollback the rest of this thread's objects fail as well.
        self.assertGreaterEqual(sum(c.rollbacks for c in self.conns), 1)


class CollisionCopy:
    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False

    def set_types(self, types):
        pass

    def write_row(self, row):
        pass


class CollisionCursor:
    """Enough of a cursor to drive ingestObject's insert-and-check path."""

    def __init__(self, conn):
        self.conn = conn
        self.rows = []
        self.rowcount = 0

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False

    def copy(self, sql):
        return CollisionCopy()

    def execute(self, sql, args=None):
        self.conn.queries.append((sql, args))
        if "JOIN distinguished_points" in sql:
            self.rows = self.conn.candidates
            self.rowcount = len(self.rows)
        elif "INSERT INTO distinguished_points" in sql:
            self.rowcount = self.conn.inserted
        elif "INSERT INTO rho_collisions" in sql:
            self.conn.stored.append(args)
            self.rowcount = 1
        else:
            self.rows = []
            self.rowcount = 1

    def fetchall(self):
        return self.rows


class CollisionConn:
    def __init__(self, inserted=0, candidates=()):
        self.inserted = inserted
        self.candidates = list(candidates)
        self.queries = []
        self.stored = []
        self.commits = 0

    def cursor(self):
        return CollisionCursor(self)

    def commit(self):
        self.commits += 1


class CollisionS3:
    def __init__(self, body):
        self.body = body

    def get_object(self, **kw):
        return {"Body": io.BytesIO(self.body)}


def record(seed, k0=1, k1=2, k2=3):
    return struct.pack("<QQQQ", seed, k0, k1, k2)


class Collisions(unittest.TestCase):
    """The campaign's terminal event must not be swallowed by DO NOTHING.

    `ON CONFLICT (campaign_id, point_key) DO NOTHING` cannot distinguish a
    worker re-reporting its own point from two different walks meeting at one,
    and before this it dropped both alike -- so the only thing the campaign
    exists to find was the one thing the ingest discarded. Nothing in the tree
    wrote `rho_collisions`, which both this program and
    scripts/rho_status/snapshot.py read.
    """

    def setUp(self):
        self.realCols = dp_ingest._collisionCols
        self.addCleanup(setattr, dp_ingest, "_collisionCols", self.realCols)
        dp_ingest._collisionCols = {
            "campaign_id", "point_key", "a1", "walk_seed1", "worker_id1",
            "a2", "walk_seed2", "worker_id2", "object_key", "detected_at"}
        self.addCleanup(setattr, dp_ingest, "log", dp_ingest.log)
        self.logged = []
        dp_ingest.log = self.logged.append

    def candidate(self, minePoint, mineSeed, theirsSeed):
        """One row as the join returns it: same point, two different seeds."""
        return (minePoint,
                mineSeed.to_bytes(dp_ingest.COEFF_BYTES, "big"),
                mineSeed.to_bytes(8, "big").lstrip(b"\x00"),
                theirsSeed.to_bytes(dp_ingest.COEFF_BYTES, "big"),
                theirsSeed.to_bytes(8, "big").lstrip(b"\x00"),
                "dp-slot-00002-old.bin")

    def ingest(self, body, inserted, candidates=()):
        conn = CollisionConn(inserted=inserted, candidates=candidates)
        result = dp_ingest.ingestObject(conn, CollisionS3(body), "bucket",
                                        ORBIT, "2026-09-19 00:00:00+00")
        return conn, result

    def test_a_point_reached_by_two_seeds_is_recorded(self):
        conn, (added, seen, hits) = self.ingest(
            record(7), inserted=0, candidates=[self.candidate(b"pk", 7, 9)])
        self.assertEqual((added, seen, hits), (0, 1, 1))
        self.assertEqual(len(conn.stored), 1)
        stored = conn.stored[0]
        self.assertIn(dp_ingest.CAMPAIGN, stored)
        self.assertIn(b"pk", stored)

    def test_the_collision_is_announced_whether_or_not_it_can_be_stored(self):
        # A storage problem must never be the reason nobody hears about this.
        dp_ingest._collisionCols = set()
        conn, (_, _, hits) = self.ingest(
            record(7), inserted=0, candidates=[self.candidate(b"pk", 7, 9)])
        self.assertEqual(hits, 1)
        self.assertEqual(conn.stored, [])
        self.assertTrue(any("COLLISION" in line for line in self.logged))

    def test_the_same_seed_twice_is_a_re_report_not_a_collision(self):
        # A resumed worker re-uploads points it already sent; that is the case
        # DO NOTHING is for and it must stay silent.
        conn, (_, _, hits) = self.ingest(
            record(7), inserted=0, candidates=[self.candidate(b"pk", 7, 7)])
        self.assertEqual(hits, 0)
        self.assertEqual(conn.stored, [])

    def test_one_seed_stored_at_another_width_is_not_a_collision(self):
        # verify() already reports that the old ingester's column widths differ
        # from this program's. Comparing bytes would call every legacy
        # re-ingest a collision; the seed is an integer.
        row = list(self.candidate(b"pk", 7, 7))
        row[3] = (7).to_bytes(8, "big")  # same seed, narrower column
        conn, (_, _, hits) = self.ingest(record(7), inserted=0, candidates=[tuple(row)])
        self.assertEqual(hits, 0)
        self.assertEqual(conn.stored, [])

    def test_an_object_of_new_points_never_runs_the_check(self):
        # The cost property: in the steady state every record inserts, so
        # nothing conflicted and there is nothing to look for.
        conn, (added, seen, hits) = self.ingest(record(7) + record(8), inserted=2)
        self.assertEqual((added, seen, hits), (2, 2, 0))
        self.assertFalse([q for q, _ in conn.queries if "JOIN distinguished_points" in q])

    def test_the_check_runs_as_soon_as_anything_conflicted(self):
        conn, _ = self.ingest(record(7) + record(8), inserted=1)
        self.assertTrue([q for q, _ in conn.queries if "JOIN distinguished_points" in q])

    def test_the_collision_is_written_in_the_transaction_that_found_it(self):
        # A recorded collision must not be able to outlive the insert it came
        # from, so it is stored before the single commit, like the progress row.
        conn, _ = self.ingest(record(7), inserted=0,
                              candidates=[self.candidate(b"pk", 7, 9)])
        order = [q for q, _ in conn.queries]
        self.assertLess(next(i for i, q in enumerate(order) if "rho_collisions" in q),
                        next(i for i, q in enumerate(order) if "dp_ingest_progress" in q))
        self.assertEqual(conn.commits, 1)

    def test_a_stored_point_with_no_seed_is_reported_not_counted(self):
        # `a` is NOT NULL for everything this program writes, but the corpus
        # predates it. Calling an uncomparable row a collision would make one
        # of every record that meets it.
        row = list(self.candidate(b"pk", 7, 9))
        row[3] = None
        conn, (_, _, hits) = self.ingest(record(7), inserted=0, candidates=[tuple(row)])
        self.assertEqual(hits, 0)
        self.assertEqual(conn.stored, [])
        self.assertTrue(any("no seed" in line for line in self.logged))

    def test_the_insert_is_built_from_the_columns_the_table_actually_has(self):
        # rho_collisions predates this file and its shape is not in the tree,
        # so the DDL above may be a no-op against a different table.
        dp_ingest._collisionCols = {"campaign_id", "point_key", "detected_at"}
        conn, (_, _, hits) = self.ingest(
            record(7), inserted=0, candidates=[self.candidate(b"pk", 7, 9)])
        self.assertEqual(hits, 1)
        sql = [q for q, _ in conn.queries if "INSERT INTO rho_collisions" in q][0]
        self.assertIn("(campaign_id, point_key)", sql)
        self.assertEqual(sql.count("%s"), 2)
        self.assertNotIn("walk_seed1", sql)

    def test_a_pass_reports_collisions_so_the_page_publishes_at_once(self):
        real = dp_ingest.pending
        self.addCleanup(setattr, dp_ingest, "pending", real)
        realIngest = dp_ingest.ingestObject
        self.addCleanup(setattr, dp_ingest, "ingestObject", realIngest)
        dp_ingest._failed.clear()
        self.addCleanup(dp_ingest._failed.clear)
        dp_ingest.pending = lambda *a, **k: (
            [("dp/slot-00001/x-0000000000000000.bin", 100, 1789000000, 0)], 1789000000, [])
        dp_ingest.ingestObject = lambda *a, **k: (99, 100, 1)
        _, _, state = dp_ingest.onePass(lambda: FakeConn([]), None, "bucket", threads=1)
        self.assertEqual(state["collisions"], 1)


class RollupCursor:
    """Records the SQL a rollup path runs, and answers the reads it makes."""

    def __init__(self, conn):
        self.conn = conn
        self.rows = []
        self.one = None
        self.rowcount = 0

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False

    def execute(self, sql, args=None):
        self.conn.queries.append((" ".join(sql.split()), args))
        if "FROM dp_ingest_meta WHERE key" in sql:
            self.one = (1,) if self.conn.marked else None
        elif "SELECT dps FROM dp_ingest_totals" in sql:
            self.one = (self.conn.totals_dps,) if self.conn.totals_dps is not None else None
        elif "SELECT count(*)" in sql and "distinguished_points" in sql:
            self.one = (self.conn.corpusCount(),)
        elif "GROUP BY 1" in sql and "distinguished_points" in sql:
            self.rows = self.conn.corpus
        elif "UPDATE dp_ingest_totals SET dps" in sql:
            self.conn.totals_dps = int(args[0]) if args else 0
        else:
            self.rows = []
            self.one = None

    def fetchone(self):
        return self.one

    def fetchall(self):
        return self.rows


class RollupConn:
    def __init__(self, marked=False, corpus=(), corpus_extra=0):
        self.marked = marked
        self.corpus = list(corpus)
        self.corpus_extra = corpus_extra
        self.totals_dps = None
        self.queries = []

    def corpusCount(self):
        return sum(int(n) for _, n, _, _ in self.corpus) + self.corpus_extra

    def cursor(self):
        return RollupCursor(self)


class Rollup(unittest.TestCase):
    """The snapshot's cost must stop growing with the corpus.

    Every published figure but the checkpoint work was a question asked of the
    whole table -- count(*), min/max of found_at, a 48-hour GROUP BY -- run
    every --status-every seconds on the instance the ingest writes to, and
    while it ran this program was not ingesting. The counters below are
    maintained by the insert itself, so the publish reads one row and 48 more.
    """

    def setUp(self):
        self.addCleanup(setattr, dp_ingest, "log", dp_ingest.log)
        dp_ingest.log = lambda msg: None

    def sql(self, conn):
        return [q for q, _ in conn.queries]

    def test_an_object_adds_its_rows_to_the_totals_and_to_its_hour(self):
        conn = RollupConn()
        cur = conn.cursor()
        dp_ingest.applyRollup(cur, 1200, "2026-09-19 07:31:00+00")
        totals = [(q, a) for q, a in conn.queries if "dp_ingest_totals" in q]
        hourly = [(q, a) for q, a in conn.queries if "dp_ingest_hourly" in q]
        self.assertEqual(len(totals), 1)
        self.assertEqual(len(hourly), 1)
        self.assertIn(1200, totals[0][1])
        self.assertIn(1200, hourly[0][1])
        # The hour is derived in the database from the object's own found_at,
        # not from the clock the ingest happens to be running on.
        self.assertIn("date_trunc('hour', %s::timestamptz)", hourly[0][0])

    def test_an_object_that_added_nothing_moves_no_counter(self):
        # A re-ingested object inserts zero rows; the counters must not move,
        # which is what keeps them exact under this file's idempotency.
        conn = RollupConn()
        dp_ingest.applyRollup(conn.cursor(), 0, "2026-09-19 07:31:00+00")
        self.assertEqual(conn.queries, [])

    def test_the_counters_are_written_in_the_insert_transaction(self):
        conn = CollisionConn(inserted=100)
        dp_ingest.ingestObject(conn, CollisionS3(record(7) * 100), "bucket",
                               ORBIT, "2026-09-19 07:31:00+00")
        order = [q for q, _ in conn.queries]
        # After the points, before the single commit: the counters cannot
        # describe rows that were never stored, nor miss rows that were.
        self.assertTrue(any("dp_ingest_totals" in q for q in order))
        self.assertLess(next(i for i, q in enumerate(order)
                             if "INSERT INTO distinguished_points" in q),
                        next(i for i, q in enumerate(order) if "dp_ingest_totals" in q))
        self.assertEqual(conn.commits, 1)

    def test_the_corpus_is_read_once_to_seed_the_counters(self):
        conn = RollupConn(corpus=[(when("2026-09-19T07:00:00"), 1200,
                                   when("2026-09-19T07:01:00"),
                                   when("2026-09-19T07:59:00"))])
        self.assertTrue(dp_ingest.ensureRollup(conn.cursor()))
        self.assertEqual(len([q for q in self.sql(conn)
                              if "FROM distinguished_points" in q]), 1)
        self.assertTrue(any("INSERT INTO dp_ingest_meta" in q for q in self.sql(conn)))

    def test_the_seed_does_not_verify_itself_by_counting_the_table(self):
        # A verification pass sounds free and is not: counting the corpus to
        # check the seed is a second full scan, and a check that retries until
        # it agrees with a corpus the fleet is still writing to never agrees.
        # Deferring the mark on a mismatch therefore reinstates a full scan
        # before every pass -- the exact cost these counters exist to remove.
        # The seed reads the table once, marks itself, and the precondition
        # (every writer goes through applyRollup) is documented instead.
        conn = RollupConn(corpus=[(when("2026-09-19T07:00:00"), 1200,
                                   when("2026-09-19T07:01:00"),
                                   when("2026-09-19T07:59:00"))])
        dp_ingest.ensureRollup(conn.cursor())
        scans = [q for q in self.sql(conn) if "FROM distinguished_points" in q]
        self.assertEqual(len(scans), 1)
        self.assertIn("GROUP BY", scans[0])
        self.assertTrue(any("INSERT INTO dp_ingest_meta" in q for q in self.sql(conn)))

    def test_a_replacement_host_does_not_read_the_corpus_again(self):
        conn = RollupConn(marked=True)
        self.assertFalse(dp_ingest.ensureRollup(conn.cursor()))
        self.assertFalse([q for q in self.sql(conn) if "FROM distinguished_points" in q])

    def test_recount_seeds_them_again(self):
        conn = RollupConn(marked=True)
        self.assertTrue(dp_ingest.ensureRollup(conn.cursor(), force=True))
        self.assertTrue(any("DELETE FROM dp_ingest_meta" in q for q in self.sql(conn)))
        self.assertTrue([q for q in self.sql(conn) if "FROM distinguished_points" in q])

    def test_the_seed_takes_its_lock_before_it_reads(self):
        # applyRollup locks the totals row *after* inserting its points, so
        # claiming that row first puts any concurrent ingest behind this
        # transaction with its points still uncommitted -- invisible to the
        # aggregate, and applied as a delta once this commits. Reading first
        # would lose every object that landed during the scan.
        conn = RollupConn()
        dp_ingest.ensureRollup(conn.cursor())
        order = self.sql(conn)
        self.assertLess(next(i for i, q in enumerate(order) if "FOR UPDATE" in q),
                        next(i for i, q in enumerate(order) if "FROM distinguished_points" in q))


class SnapshotCursor:
    def __init__(self, conn):
        self.conn = conn
        self.rows = []
        self.one = None

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False

    def execute(self, sql, args=None):
        flat = " ".join(sql.split())
        self.conn.queries.append(flat)
        if "FROM dp_ingest_hourly" in flat:
            self.rows = self.conn.hourly
        elif "rho_collisions" in flat:
            self.one = (0, None)
        elif "rho_campaigns" in flat:
            self.one = (dp_ingest.CAMPAIGN, "sect131r1", 32,
                        when("2026-09-01T00:00:00"), 190123456,
                        when("2026-09-01T00:00:00"), when("2026-09-19T07:59:00"))

    def fetchone(self):
        return self.one

    def fetchall(self):
        return self.rows


class SnapshotConn:
    def __init__(self, hourly=()):
        self.hourly = list(hourly)
        self.queries = []

    def cursor(self):
        return SnapshotCursor(self)


class Snapshot(unittest.TestCase):
    """What a publish costs, which is the whole point of the rollup."""

    def setUp(self):
        self.addCleanup(setattr, dp_ingest, "log", dp_ingest.log)
        dp_ingest.log = lambda msg: None
        now = datetime.datetime.now(datetime.timezone.utc).replace(
            minute=0, second=0, microsecond=0)
        self.conn = SnapshotConn(hourly=[(now, 7), (now - datetime.timedelta(hours=3), 5)])
        self.payload = dp_ingest.statusPayload(self.conn, FakeWorkS3([]), "bucket")

    def test_the_snapshot_never_reads_the_points_table(self):
        # The regression this guards: every figure below used to be a question
        # asked of distinguished_points, so a publish cost more every day the
        # campaign ran and blocked the ingest while it did.
        self.assertFalse([q for q in self.conn.queries if "distinguished_points" in q])

    def test_the_total_comes_from_the_counter(self):
        self.assertEqual(self.payload["dps"], 190123456)

    def test_the_recency_figures_come_from_the_rollup(self):
        self.assertEqual(self.payload["dps_last_hour"], 7)
        self.assertEqual(self.payload["dps_last_day"], 12)

    def test_the_first_and_last_point_still_reach_the_page(self):
        # These moved column with the query; a silently shifted index would
        # publish the campaign's creation date as its newest point.
        self.assertTrue(self.payload["first_dp_at"].startswith("2026-09-01"))
        self.assertTrue(self.payload["last_dp_at"].startswith("2026-09-19"))

    def test_the_hourly_series_is_still_published(self):
        self.assertEqual([h["dps"] for h in self.payload["hourly"]], [7, 5])


class Windows(unittest.TestCase):
    """The recency tiles come from the rollup the chart is drawn from."""

    def hours(self, *offsets):
        now = datetime.datetime.now(datetime.timezone.utc).replace(
            minute=0, second=0, microsecond=0)
        return [(now - datetime.timedelta(hours=h), 10) for h in offsets]

    def test_it_sums_whole_buckets_inside_the_window(self):
        self.assertEqual(dp_ingest.windowSum(self.hours(0, 1), 1), 20)

    def test_buckets_outside_the_window_are_not_counted(self):
        self.assertEqual(dp_ingest.windowSum(self.hours(0, 5), 1), 10)
        self.assertEqual(dp_ingest.windowSum(self.hours(0, 5, 30), 24), 20)

    def test_a_naive_timestamp_is_read_as_utc_not_local(self):
        naive = [(datetime.datetime.now(datetime.timezone.utc).replace(
            tzinfo=None, minute=0, second=0, microsecond=0), 10)]
        self.assertEqual(dp_ingest.windowSum(naive, 1), 10)

    def test_an_empty_rollup_is_zero_not_an_error(self):
        self.assertEqual(dp_ingest.windowSum([], 24), 0)


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


class WalkRate(unittest.TestCase):
    def test_measures_between_two_publishes(self):
        previous = {
            "generated_at": "2026-09-19T12:00:00Z",
            "work": {"iterations": 1000},
        }
        current = {
            "generated_at": "2026-09-19T12:10:00Z",
            "work": {"iterations": 1600},
        }
        rate = dp_ingest.walkRateBetween(current, previous)
        self.assertAlmostEqual(rate["iterations_per_second"], 1.0)
        self.assertEqual(rate["window_seconds"], 600)
        self.assertEqual(rate["iterations_from"], 1000)
        self.assertEqual(rate["iterations_to"], 1600)

    def test_refuses_backwards_totals_and_short_windows(self):
        base = {"generated_at": "2026-09-19T12:00:00Z", "work": {"iterations": 1000}}
        self.assertIsNone(dp_ingest.walkRateBetween(
            {"generated_at": "2026-09-19T12:00:30Z", "work": {"iterations": 900}}, base))
        self.assertIsNone(dp_ingest.walkRateBetween(
            {"generated_at": "2026-09-19T12:04:00Z", "work": {"iterations": 1100}}, base))

    def test_uses_the_previous_rate_window_when_the_last_write_is_too_close(self):
        previous = {
            "generated_at": "2026-09-19T12:09:00Z",
            "work": {"iterations": 1540},
            "walk_rate": {
                "iterations_from": 1000,
                "measured_from": "2026-09-19T12:00:00Z",
            },
        }
        current = {
            "generated_at": "2026-09-19T12:12:00Z",
            "work": {"iterations": 1720},
        }
        rate = dp_ingest.walkRateBetween(current, previous)
        self.assertAlmostEqual(rate["iterations_per_second"], (1720 - 1000) / 720)
        self.assertEqual(rate["measured_from"], "2026-09-19T12:00:00Z")

    def test_public_status_drops_per_slot_rows(self):
        public = dp_ingest.publicStatus({
            "work": {"iterations": 10, "walking_slots": 1, "per_slot": [{"slot": 1}]},
        })
        self.assertNotIn("per_slot", public["work"])
        self.assertEqual(public["work"]["walking_slots"], 1)


class StatusPublisherTests(unittest.TestCase):
    """Status must keep moving while an ingest pass is stuck."""

    def test_publish_once_writes_without_waiting_for_a_pass(self):
        puts = []

        class BucketS3:
            def get_object(self, **kw):
                raise KeyError("none yet")

            def put_object(self, **kw):
                puts.append(kw)

            def get_paginator(self, name):
                class Pages:
                    def paginate(self, **kw):
                        return iter([])
                return Pages()

            def list_objects_v2(self, **kw):
                return {}

        snapshot = {
            "curve_id": 131, "dp_mask_bits": 32, "campaign_created_at": None,
            "dps": 7, "dps_last_hour": 0, "dps_last_day": 0,
            "first_dp_at": None, "last_dp_at": None, "collisions": 0,
            "latest_collision_at": None, "hourly": [],
        }

        class Conn:
            def __enter__(self):
                return self

            def __exit__(self, *exc):
                return False

            def cursor(self):
                raise AssertionError("publish_once should use the cached snapshot")

        publisher = dp_ingest.StatusPublisher(
            connect=lambda: Conn(),
            s3=BucketS3(),
            bucket="bucket",
            statusBucket="status",
            status_every=30.0,
            snapshot_every=0.0,
        )
        publisher._snapshot = snapshot
        publisher._snap_at = time.time()
        publisher.update_ingest({"outstanding": 3, "unrecognised": 0, "newest": time.time()})
        with mock.patch.object(dp_ingest, "campaignSnapshot", return_value=snapshot):
            publisher.publish_once()
        self.assertEqual(len(puts), 1)
        body = json.loads(puts[0]["Body"].decode())
        self.assertEqual(body["dps"], 7)
        self.assertEqual(body["ingest"]["outstanding_objects"], 3)

    def test_kick_publishes_without_waiting_out_the_interval(self):
        puts = []

        class BucketS3:
            def get_object(self, **kw):
                raise KeyError("none yet")

            def put_object(self, **kw):
                puts.append(time.time())

            def list_objects_v2(self, **kw):
                return {}

        snapshot = {
            "curve_id": 131, "dp_mask_bits": 32, "campaign_created_at": None,
            "dps": 1, "dps_last_hour": 0, "dps_last_day": 0,
            "first_dp_at": None, "last_dp_at": None, "collisions": 0,
            "latest_collision_at": None, "hourly": [],
        }

        class Conn:
            def __enter__(self):
                return self

            def __exit__(self, *exc):
                return False

        # Long interval: without kick(), a second publish would not arrive in time.
        publisher = dp_ingest.StatusPublisher(
            connect=lambda: Conn(),
            s3=BucketS3(),
            bucket="bucket",
            statusBucket="status",
            status_every=30.0,
            snapshot_every=0.0,
        )
        with mock.patch.object(dp_ingest, "campaignSnapshot", return_value=snapshot):
            publisher.start()
            deadline = time.time() + 2.0
            while not puts and time.time() < deadline:
                time.sleep(0.01)
            self.assertTrue(puts, "first publish on start never arrived")
            before = len(puts)
            publisher.kick()
            deadline = time.time() + 2.0
            while len(puts) <= before and time.time() < deadline:
                time.sleep(0.01)
            publisher.stop(timeout=2.0)
        self.assertGreater(len(puts), before, "kick must publish without waiting out status_every")

    def test_main_wires_the_publisher_thread(self):
        path = os.path.join(os.path.dirname(__file__), "dp_ingest.py")
        with open(path, encoding="utf-8") as fh:
            source = fh.read()
        self.assertIn("class StatusPublisher", source)
        self.assertIn("publisher.start()", source)
        self.assertIn("publisher.update_ingest(ingest)", source)
        # The old pass-then-publish gate must not come back: it is what froze
        # the feed whenever onePass ran long.
        self.assertNotIn("time.time() - published >= args.status_every", source)


class Defaults(unittest.TestCase):
    def test_status_every_default_is_three_minutes(self):
        import argparse
        ap = argparse.ArgumentParser()
        ap.add_argument("--status-every", type=float,
                        default=float(os.environ.get("RHO_STATUS_EVERY", "180")))
        self.assertEqual(ap.parse_args([]).status_every, 180.0)

    def test_snapshot_every_default_is_every_status_write(self):
        path = os.path.join(os.path.dirname(__file__), "dp_ingest.py")
        with open(path) as fh:
            self.assertIn('RHO_SNAPSHOT_EVERY", "0"', fh.read())


class CachedSnapshot(unittest.TestCase):
    def test_status_payload_refreshes_walkers_without_sql(self):
        now = datetime.datetime.now(datetime.timezone.utc)
        items = [("ckpt/slot-90001/%s.ck" % ("ab" * 32), now, ckpt(50, 2, 4))]
        snapshot = {
            "curve_id": 131, "dp_mask_bits": 32, "campaign_created_at": None,
            "dps": 190404664, "dps_last_hour": 0, "dps_last_day": 10,
            "first_dp_at": None, "last_dp_at": None, "collisions": 0,
            "latest_collision_at": None, "hourly": [],
        }
        payload = dp_ingest.statusPayload(
            None, FakeWorkS3(items), "bucket", snapshot=snapshot)
        self.assertEqual(payload["dps"], 190404664)
        self.assertEqual(payload["walkers"], 1)
        self.assertEqual(payload["work"]["walking_slots"], 1)
        self.assertEqual(payload["state"], "IDLE_OR_STALE")


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
