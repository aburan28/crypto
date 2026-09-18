#!/usr/bin/env python3
"""Run snapshot.py's rollup against a real Postgres.

The unit tests pin the *shape* of the backfill; the claims that matter are
about what Postgres does with it. Two of them cost the campaign a day:

  * the first backfill held SHARE on distinguished_points and never finished
    inside its timeout, so it blocked the ingest and published nothing, and
  * a backfill that runs beside the insert trigger must count each row once,
    which is true only because each hour is a single absolute upsert.

Neither is visible without a server. This module skips unless one is reachable,
so it is a no-op in a CI container and a real check anywhere it can connect:

    RHO_TEST_DATABASE_URL=postgresql://... python3 -m unittest test_rollup_postgres

If RHO_TEST_DATABASE_URL is unset it will start a throwaway cluster with
pg_ctl when initdb is on PATH (see PG_BIN), and stop it afterwards.
"""

import os
import shutil
import subprocess
import tempfile
import time
import unittest
from datetime import datetime

import snapshot

PG_BIN = os.environ.get("PG_BIN", "/usr/lib/postgresql/16/bin")
CAMPAIGN = "ecc2k-test"

SCHEMA = """
DROP TABLE IF EXISTS distinguished_points, rho_campaigns, rho_collisions,
                     rho_dp_hour, rho_dp_recent, rho_dp_meta CASCADE;
CREATE TABLE rho_campaigns (
  campaign_id  text PRIMARY KEY,
  curve_id     text,
  dp_mask_bits integer,
  created_at   timestamptz NOT NULL DEFAULT now(),
  meta         jsonb
);
CREATE TABLE rho_collisions (
  campaign_id text NOT NULL,
  detected_at timestamptz NOT NULL DEFAULT now()
);
CREATE TABLE distinguished_points (
  campaign_id text NOT NULL,
  worker_id   text NOT NULL,
  point_key   bytea NOT NULL,
  found_at    timestamptz NOT NULL,
  PRIMARY KEY (campaign_id, point_key)
);
CREATE INDEX distinguished_points_campaign_found_at
  ON distinguished_points (campaign_id, found_at);
INSERT INTO rho_campaigns (campaign_id, curve_id, dp_mask_bits)
VALUES ('%s', 'certicom-ecc2k-130', 32);
""" % CAMPAIGN

# Six hours of four workers, plus a two-hour gap the backfill must step over in
# one index probe rather than a transaction per empty hour.
FIXTURE = """
INSERT INTO distinguished_points (campaign_id, worker_id, point_key, found_at)
SELECT '%s',
       'w' || (g %% 4),
       decode(lpad(to_hex(g), 16, '0'), 'hex'),
       date_trunc('hour', now()) - interval '9 hours'
         + (CASE WHEN g %% 6 < 3 THEN g %% 6 ELSE (g %% 6) + 2 END) * interval '1 hour'
         + (g %% 47) * interval '7 seconds'
FROM generate_series(1, 20000) g;
""" % CAMPAIGN


def instant(text):
    """Compare timestamps as moments; psql and the document spell them apart."""
    text = text.strip().replace(" ", "T").replace("Z", "+00:00")
    if text.endswith("+00"):
        text += ":00"
    return datetime.fromisoformat(text)


def pg_available():
    return bool(os.environ.get("RHO_TEST_DATABASE_URL")) or os.path.exists(
        os.path.join(PG_BIN, "initdb"))


@unittest.skipUnless(pg_available(), "no Postgres: set RHO_TEST_DATABASE_URL or install initdb")
class RollupAgainstPostgres(unittest.TestCase):
    cluster = None
    url = None

    @classmethod
    def setUpClass(cls):
        cls.url = os.environ.get("RHO_TEST_DATABASE_URL")
        if cls.url:
            return
        cls.cluster = tempfile.mkdtemp(prefix="rho-pg-")
        data = os.path.join(cls.cluster, "data")
        env = dict(os.environ, PATH=PG_BIN + os.pathsep + os.environ.get("PATH", ""))
        subprocess.run([os.path.join(PG_BIN, "initdb"), "-D", data, "-U", "postgres",
                        "--auth=trust", "-E", "UTF8"],
                       check=True, capture_output=True, text=True, env=env)
        # -l matters: without it the postmaster inherits the captured stdout and
        # holds the pipe open, so subprocess.run waits for a server that is
        # already up.
        subprocess.run([os.path.join(PG_BIN, "pg_ctl"), "-D", data, "-w",
                        "-l", os.path.join(cls.cluster, "server.log"),
                        "-o", "-p 5439 -k %s -c listen_addresses=" % cls.cluster, "start"],
                       check=True, capture_output=True, text=True, env=env)
        cls.url = "postgresql://postgres@/postgres?host=%s&port=5439" % cls.cluster

    @classmethod
    def tearDownClass(cls):
        if cls.cluster:
            data = os.path.join(cls.cluster, "data")
            subprocess.run([os.path.join(PG_BIN, "pg_ctl"), "-D", data, "-m", "immediate", "stop"],
                           capture_output=True, text=True)
            shutil.rmtree(cls.cluster, ignore_errors=True)

    def sql(self, text, tuples_only=True):
        return snapshot.psql_script(self.url, text, tuples_only=tuples_only)

    def scalar(self, text):
        return self.sql(text).strip().splitlines()[-1].strip()

    def setUp(self):
        self.sql(SCHEMA, tuples_only=False)
        self.sql(FIXTURE, tuples_only=False)
        self.sql(snapshot.ENSURE_SQL, tuples_only=False)

    def total(self):
        return int(self.scalar(
            "SELECT count(*) FROM distinguished_points WHERE campaign_id = '%s'" % CAMPAIGN))

    def rolled(self):
        return int(self.scalar(
            "SELECT COALESCE(sum(dps),0) FROM rho_dp_hour WHERE campaign_id = '%s'" % CAMPAIGN))

    def backfill(self, budget=600):
        return snapshot.run_backfill(self.url, CAMPAIGN, budget)

    def test_the_backfill_reproduces_the_table_and_marks_itself_ready(self):
        self.assertTrue(self.backfill())
        self.assertEqual(self.rolled(), self.total())
        self.assertEqual(self.scalar(
            "SELECT ready FROM rho_dp_meta WHERE campaign_id = '%s'" % CAMPAIGN), "t")
        hours = int(self.scalar(
            "SELECT count(DISTINCT hour) FROM rho_dp_hour WHERE campaign_id = '%s'" % CAMPAIGN))
        self.assertEqual(hours, 6, "six populated hours, and the gap contributes none")

    def test_it_holds_no_lock_on_the_table_the_ingest_writes(self):
        """An INSERT must go through while the backfill is between hours.

        The shape this replaces took SHARE on distinguished_points for the
        whole attempt, so this is the regression that matters most: run the
        backfill, and assert nothing is left holding a conflicting mode.
        """
        self.assertTrue(self.backfill())
        modes = self.sql(
            "SELECT mode FROM pg_locks l JOIN pg_class c ON c.oid = l.relation "
            "WHERE c.relname = 'distinguished_points'")
        self.assertNotIn("ShareLock", modes)
        self.assertNotIn("ExclusiveLock", modes)
        self.sql(
            "INSERT INTO distinguished_points (campaign_id, worker_id, point_key, found_at) "
            "VALUES ('%s', 'w0', '\\xdeadbeef', now())" % CAMPAIGN, tuples_only=False)

    def test_the_backfill_runs_while_a_writer_holds_the_insert_lock(self):
        """The regression, stated as Postgres states it.

        An INSERT holds ROW EXCLUSIVE, which conflicts with SHARE and not with
        ACCESS SHARE. Hold ROW EXCLUSIVE in another session and run the
        backfill: the shape that shipped on 2026-09-18 would wait here until
        its 800 s timeout, which is the whole failure. This one finishes.
        """
        holder = subprocess.Popen(
            ["psql", self.url, "-v", "ON_ERROR_STOP=1", "-At", "-c",
             "BEGIN; LOCK TABLE distinguished_points IN ROW EXCLUSIVE MODE; "
             "SELECT pg_sleep(30); COMMIT;"],
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        try:
            # Let the lock be taken before asking for the table.
            deadline = time.monotonic() + 10
            while time.monotonic() < deadline:
                held = self.sql(
                    "SELECT mode FROM pg_locks l JOIN pg_class c ON c.oid = l.relation "
                    "WHERE c.relname = 'distinguished_points' AND l.mode = 'RowExclusiveLock'")
                if "RowExclusiveLock" in held:
                    break
                time.sleep(0.2)
            self.assertIn("RowExclusiveLock", held, "no writer to contend with")
            started = time.monotonic()
            self.assertTrue(self.backfill(budget=20))
            self.assertLess(time.monotonic() - started, 20,
                            "the backfill waited on a writer's lock")
            self.assertEqual(self.rolled(), self.total())
        finally:
            holder.kill()
            holder.wait(timeout=10)

    def test_the_trigger_and_the_backfill_together_count_a_row_once(self):
        """Insert after the backfill, then rerun it: the total must not move.

        The trigger adds and the backfill sets. Rerunning an hour the trigger
        has since touched is the interleaving that a DELETE-then-INSERT pair
        would get wrong, and it is free here.
        """
        self.assertTrue(self.backfill())
        before = self.rolled()
        self.sql(
            "INSERT INTO distinguished_points (campaign_id, worker_id, point_key, found_at) "
            "SELECT '%s', 'w1', decode(lpad(to_hex(900000+g),16,'0'),'hex'), now() "
            "FROM generate_series(1,500) g" % CAMPAIGN, tuples_only=False)
        self.assertEqual(self.rolled(), before + 500, "the trigger maintains the rollup")
        self.assertEqual(self.rolled(), self.total())
        # Rewind the cursor and run every hour again over the top.
        self.sql(
            "UPDATE rho_dp_meta SET ready = false, backfill_cursor = "
            "(SELECT date_trunc('hour', min(found_at)) FROM distinguished_points "
            " WHERE campaign_id = '%s') WHERE campaign_id = '%s'" % (CAMPAIGN, CAMPAIGN),
            tuples_only=False)
        self.assertTrue(self.backfill())
        self.assertEqual(self.rolled(), self.total(), "a second pass double-counted nothing")

    def test_a_spent_budget_leaves_a_cursor_and_the_next_run_finishes(self):
        self.assertFalse(self.backfill(budget=0), "a zero budget cannot finish six hours")
        cursor, until = snapshot.backfill_progress(self.url, CAMPAIGN)
        self.assertTrue(cursor, "an unfinished backfill records where to resume")
        self.assertTrue(until)
        partial = self.rolled()
        self.assertLess(partial, self.total())
        self.assertEqual(self.scalar(
            "SELECT ready FROM rho_dp_meta WHERE campaign_id = '%s'" % CAMPAIGN), "f")
        self.assertTrue(self.backfill())
        self.assertEqual(self.rolled(), self.total())

    def test_the_published_document_matches_the_table_it_summarises(self):
        snap = snapshot.take_snapshot(self.url, CAMPAIGN, source="test")
        self.assertEqual(snap["dps"], self.total())
        self.assertEqual(snap["workers"], 4)
        self.assertEqual(snap["state"], "COLLECTING" if snap["dps_last_hour"] else "IDLE_OR_STALE")
        per_worker = {w["worker_id"]: w["dps"] for w in snap["per_worker"]}
        self.assertEqual(sum(per_worker.values()), self.total())
        rows = self.sql(
            "SELECT date_trunc('hour', found_at)::text, count(*) FROM distinguished_points "
            "WHERE campaign_id = '%s' GROUP BY 1 ORDER BY 1" % CAMPAIGN).strip().splitlines()
        table_hours = {}
        for line in rows:
            hour, count = line.split("|")
            table_hours[instant(hour)] = int(count)
        self.assertEqual(
            {instant(h["hour"]): h["dps"] for h in snap["hourly"]},
            table_hours,
            "the hourly series is the table's own histogram")

    def test_an_unfinished_backfill_publishes_nothing(self):
        with self.assertRaises(RuntimeError) as caught:
            snapshot.take_snapshot(self.url, CAMPAIGN, source="test", budget=0)
        self.assertIn("rerun to resume", str(caught.exception))


if __name__ == "__main__":
    unittest.main()
