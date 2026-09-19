#!/usr/bin/env python3
"""Run snapshot.py's rollup backfill against a real Postgres.

The rest of the suite pins the *shape* of this SQL. The claims that cost the
campaign a day are not about its shape, they are about what Postgres does with
it concurrently:

  * the first backfill held SHARE on distinguished_points and never finished
    inside its 800 s timeout, so every 15-minute run blocked the ingest, rolled
    back, and left ready=false for the next one to repeat -- the store went 76
    minutes without ingesting an object while the page stayed on an old
    snapshot, and
  * an hour rebuilt beside a live insert trigger must count each row once,
    which is a statement about isolation levels and row locks and nothing else.

Neither is visible without a server, so this module skips unless one is
reachable, which makes it a no-op in a CI container and a real check anywhere
it can connect:

    RHO_TEST_DATABASE_URL=postgresql://... python3 -m unittest test_rollup_postgres

If RHO_TEST_DATABASE_URL is unset it starts a throwaway cluster with pg_ctl
when initdb is on PATH (see PG_BIN), and stops it afterwards.
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

# Six populated hours of four workers inside the last nine, with a two-hour
# gap: hours_to_backfill walks every hour in the span, so the empty ones have
# to cost nothing and contribute nothing.
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
        # holds the pipe open, so subprocess.run waits for a server already up.
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

    def ready(self):
        return self.scalar(
            "SELECT COALESCE((SELECT ready FROM rho_dp_meta WHERE campaign_id = '%s'), false)"
            % CAMPAIGN)

    def backfill(self, budget=600):
        return snapshot.backfill_hours(self.url, CAMPAIGN, time.monotonic(), budget)

    def test_the_backfill_reproduces_the_table_and_marks_itself_ready(self):
        self.assertTrue(self.backfill())
        self.assertEqual(self.rolled(), self.total())
        self.assertEqual(self.ready(), "t")
        hours = int(self.scalar(
            "SELECT count(DISTINCT hour) FROM rho_dp_hour WHERE campaign_id = '%s'" % CAMPAIGN))
        self.assertEqual(hours, 6, "six populated hours, and the gap contributes none")

    def test_the_backfill_runs_while_a_writer_holds_the_insert_lock(self):
        """The regression, stated the way Postgres states it.

        An INSERT holds ROW EXCLUSIVE, which conflicts with SHARE and not with
        ACCESS SHARE. Hold ROW EXCLUSIVE in another session and run the
        backfill: the shape that shipped on 2026-09-18 waits here until its
        800 s timeout, which is the whole failure. This one finishes.
        """
        holder = subprocess.Popen(
            ["psql", self.url, "-v", "ON_ERROR_STOP=1", "-At", "-c",
             "BEGIN; LOCK TABLE distinguished_points IN ROW EXCLUSIVE MODE; "
             "SELECT pg_sleep(30); COMMIT;"],
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        try:
            held = ""
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
            self.assertTrue(self.backfill(budget=600))
            self.assertLess(time.monotonic() - started, 25,
                            "the backfill waited on a writer's lock")
            self.assertEqual(self.rolled(), self.total())
        finally:
            holder.kill()
            holder.wait(timeout=10)

    def test_no_lock_on_the_table_survives_the_backfill(self):
        self.assertTrue(self.backfill())
        modes = self.sql(
            "SELECT mode FROM pg_locks l JOIN pg_class c ON c.oid = l.relation "
            "WHERE c.relname = 'distinguished_points'")
        self.assertNotIn("ShareLock", modes)
        self.assertNotIn("ExclusiveLock", modes)
        self.sql(
            "INSERT INTO distinguished_points (campaign_id, worker_id, point_key, found_at) "
            "VALUES ('%s', 'w0', '\\xdeadbeef', now())" % CAMPAIGN, tuples_only=False)

    def test_the_trigger_keeps_the_rollup_current_after_the_backfill(self):
        self.assertTrue(self.backfill())
        before = self.rolled()
        self.sql(
            "INSERT INTO distinguished_points (campaign_id, worker_id, point_key, found_at) "
            "SELECT '%s', 'w1', decode(lpad(to_hex(900000+g),16,'0'),'hex'), now() "
            "FROM generate_series(1,500) g" % CAMPAIGN, tuples_only=False)
        self.assertEqual(self.rolled(), before + 500, "the trigger maintains the rollup")
        self.assertEqual(self.rolled(), self.total())

    def test_rebuilding_an_hour_the_trigger_has_touched_counts_each_row_once(self):
        """The double-count the DELETE and the isolation level exist to prevent.

        The trigger adds to a bucket and the backfill rebuilds it. Rewind
        backfill_through so every hour, including the one the trigger has just
        written, is done again over the top of its own work.
        """
        self.assertTrue(self.backfill())
        self.sql(
            "INSERT INTO distinguished_points (campaign_id, worker_id, point_key, found_at) "
            "SELECT '%s', 'w2', decode(lpad(to_hex(800000+g),16,'0'),'hex'), now() "
            "FROM generate_series(1,700) g" % CAMPAIGN, tuples_only=False)
        self.sql(
            "UPDATE rho_dp_meta SET ready = false, backfill_through = NULL "
            "WHERE campaign_id = '%s'" % CAMPAIGN, tuples_only=False)
        self.assertTrue(self.backfill())
        self.assertEqual(self.rolled(), self.total(), "a second pass double-counted nothing")

    def test_a_spent_budget_leaves_a_cursor_and_the_next_run_finishes(self):
        self.assertFalse(self.backfill(budget=0), "a zero budget cannot finish six hours")
        self.assertEqual(self.ready(), "f")
        self.assertTrue(self.backfill(budget=600))
        self.assertEqual(self.ready(), "t")
        self.assertEqual(self.rolled(), self.total())

    def test_the_published_document_matches_the_table_it_summarises(self):
        snap = snapshot.take_snapshot(self.url, CAMPAIGN, source="test")
        self.assertEqual(snap["dps"], self.total())
        self.assertEqual(snap["workers"], 4)
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

    def test_the_fallback_count_agrees_with_the_rollup(self):
        """The two paths the page can publish from must not disagree.

        load_fallback counts the corpus an hour at a time without the rollup;
        a reader cannot tell which one produced a figure, so they have to match.
        """
        fallback = snapshot.load_fallback(self.url, CAMPAIGN, time.monotonic(), 600)
        self.assertTrue(self.backfill())
        rollup = snapshot.take_snapshot(self.url, CAMPAIGN, source="test")
        self.assertEqual(fallback["dps"], rollup["dps"])
        self.assertEqual(fallback["dps"], self.total())


if __name__ == "__main__":
    unittest.main()
