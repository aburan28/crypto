#!/usr/bin/env python3
"""Offline tests for the ECC2K-130 status snapshot and history merge."""

from __future__ import annotations

import json
import os
import re
import tempfile
import unittest
from datetime import datetime, timedelta, timezone

from render import (
    HISTORY_LIMIT,
    MAX_RATE_SPAN_S,
    MIN_RATE_SPAN_S,
    RATE_WINDOW_S,
    apply_rate,
    measure_rate,
    merge_history,
)
from snapshot import (
    BACKFILL_SQL,
    CLAIM_BOUNDARY,
    ENSURE_SQL,
    FORBIDDEN_PUBLIC_KEYS,
    READ_SQL,
    assert_public,
    campaign_state,
    normalize,
    render_sql,
    rollup_ready_from_text,
    sql_campaign,
)
from work_feed import MAX_FEED_AGE_S, feed_url, ingest_block, merge_work, work_block


HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))


class SnapshotTests(unittest.TestCase):
    def test_fixture_normalizes_to_public_document(self):
        path = os.path.join(HERE, "testdata", "status.json")
        with open(path, encoding="utf-8") as fh:
            fixture = json.load(fh)
        snapshot = normalize(
            {
                "campaign_id": fixture["campaign_id"],
                "curve_id": fixture["curve_id"],
                "dp_mask_bits": fixture["dp_mask_bits"],
                "campaign_created_at": fixture["campaign_created_at"],
                "dps": fixture["dps"],
                "workers": fixture["workers"],
                "first_dp_at": fixture["first_dp_at"],
                "last_dp_at": fixture["last_dp_at"],
                "dps_last_hour": fixture["dps_last_hour"],
                "dps_last_day": fixture["dps_last_day"],
                "collisions": fixture["collisions"],
                "latest_collision_at": None,
                "per_worker": fixture["per_worker"],
                "hourly": fixture["hourly"],
            },
            "ecc2k-130",
            "fixture",
        )
        self.assertEqual(snapshot["state"], "COLLECTING")
        self.assertEqual(snapshot["claim_boundary"], CLAIM_BOUNDARY)
        self.assertNotIn("point_key", json.dumps(snapshot))
        assert_public(snapshot)

    def test_collision_state(self):
        self.assertEqual(campaign_state({"collisions": 1, "dps": 10, "dps_last_hour": 0}), "COLLISION_RECORDED")
        self.assertEqual(campaign_state({"collisions": 0, "dps": 10, "dps_last_hour": 0}), "IDLE_OR_STALE")
        self.assertEqual(campaign_state({"collisions": 0, "dps": 0, "dps_last_hour": 0}), "EMPTY")
        self.assertEqual(
            campaign_state({
                "collisions": 0, "dps": 10, "dps_last_hour": 0,
                "ingest_outstanding": 12,
            }),
            "INGEST_BEHIND",
        )
        self.assertEqual(
            campaign_state({
                "collisions": 0, "dps": 10, "dps_last_hour": 5,
                "ingest_outstanding": 12,
            }),
            "COLLECTING",
        )

    def test_refuses_walk_secrets(self):
        with self.assertRaises(SystemExit):
            assert_public({"point_key": "00"})
        for key in FORBIDDEN_PUBLIC_KEYS:
            self.assertTrue(key)


class WorkflowTests(unittest.TestCase):
    def test_publish_job_uses_access_keys_and_deploy_key(self):
        path = os.path.join(ROOT, ".github", "workflows", "ecc2k130-status.yml")
        with open(path, encoding="utf-8") as fh:
            text = fh.read()
        self.assertIn("aws-access-key-id: ${{ secrets.AWS_ACCESS_KEY_ID }}", text)
        self.assertIn("aws-secret-access-key: ${{ secrets.AWS_SECRET_ACCESS_KEY }}", text)
        self.assertIn("secrets.RHO_WALKER_SSH_KEY", text)
        self.assertNotIn("role-to-assume", text)

    def test_history_window_matches_the_publish_cadence(self):
        # history.json is trimmed to a snapshot COUNT, so the window it covers
        # is that count divided by the cron rate. The two live in different
        # files, and speeding the cron up without raising the count silently
        # shortens the published series instead of failing, so pin both to the
        # same seven days here.
        path = os.path.join(ROOT, ".github", "workflows", "ecc2k130-status.yml")
        with open(path, encoding="utf-8") as fh:
            text = fh.read()
        crons = re.findall(r'- cron: "([^"]+)"', text)
        self.assertEqual(len(crons), 1, "expected one schedule, got %s" % (crons,))
        minute, rest = crons[0].split(" ", 1)
        self.assertEqual(rest, "* * * *", "not an every-N-minutes cron: %s" % crons[0])
        step = int(minute[2:]) if minute.startswith("*/") else 60
        self.assertEqual(
            HISTORY_LIMIT,
            7 * (24 * 60 // step),
            "cron runs every %d min, so seven days is %d snapshots, but "
            "render.HISTORY_LIMIT is %d" % (step, 7 * (24 * 60 // step), HISTORY_LIMIT),
        )

    def test_dashboard_stale_banner_allows_for_scheduler_jitter(self):
        # The banner must not cry stale on one late run: GitHub delays
        # scheduled workflows under load. Keep the threshold at several
        # missed runs, and never below two.
        workflow = os.path.join(ROOT, ".github", "workflows", "ecc2k130-status.yml")
        with open(workflow, encoding="utf-8") as fh:
            crons = re.findall(r'- cron: "([^"]+)"', fh.read())
        minute = crons[0].split(" ", 1)[0]
        step = int(minute[2:]) if minute.startswith("*/") else 60
        page = os.path.join(ROOT, "docs", "ecc2k130-status", "index.html")
        with open(page, encoding="utf-8") as fh:
            text = fh.read()
        match = re.search(r"var STALE_AFTER_MS = ([^;]+);", text)
        self.assertIsNotNone(match, "STALE_AFTER_MS not found on the dashboard")
        expression = match.group(1).replace("HOUR_MS", str(3600 * 1000)).strip()
        self.assertRegex(expression, r"^[\d\s*]+$", "unexpected STALE_AFTER_MS: %s" % expression)
        stale_minutes = eval(expression) / 1000 / 60  # noqa: S307 - digits and * only
        self.assertGreaterEqual(stale_minutes, 2 * step, "stale banner will flap on a late run")
        self.assertIn("INGEST_BEHIND", text)
        self.assertIn(
            "The fleet is walking, but the store is behind on ingest",
            text,
        )

    def test_snapshot_reads_the_rollup_and_scans_the_heap_only_to_backfill(self):
        # Three separate LATERAL scans of a 30M-row table is what pushed the
        # walker hop from 41s to 510s and broke publication; grouping the heap
        # once by worker and hour then pushed it to 524s at 187M rows, one
        # timeout tick under 540s, and the next scheduled run died. The
        # published document is assembled from rho_dp_hour / rho_dp_recent.
        # distinguished_points is read only by the one-time backfill.
        self.assertEqual(BACKFILL_SQL.count("FROM distinguished_points"), 1)
        self.assertNotIn("LOCK TABLE", BACKFILL_SQL)
        self.assertIn("date_trunc('hour', found_at)", BACKFILL_SQL)
        self.assertIn("ON CONFLICT (campaign_id, hour, worker_id)", BACKFILL_SQL)
        self.assertNotIn("FROM distinguished_points", READ_SQL)
        self.assertIn("FROM rho_dp_hour", READ_SQL)
        self.assertIn("FROM rho_dp_recent", READ_SQL)
        self.assertIn("CREATE TRIGGER", ENSURE_SQL)
        self.assertIn("rho_dp_rollup_insert", ENSURE_SQL)
        self.assertIn("ON distinguished_points", ENSURE_SQL)
        # Ingest INSERTs fire this as the rho/dp-rds role, which is not the
        # walker DATABASE_URL role that creates the rollup tables.
        self.assertIn("SECURITY DEFINER", ENSURE_SQL)
        self.assertGreaterEqual(BACKFILL_SQL.count("{{campaign}}"), 1)
        self.assertIn("{{campaign}}", READ_SQL)
        with open(os.path.join(HERE, "snapshot.py"), encoding="utf-8") as fh:
            source = fh.read()
        self.assertNotIn('replace("%(campaign)s", "%s")', source)
        self.assertIn("if rollup_ready_from_text", source)

    def test_campaign_id_is_a_literal_not_concatenated_sql(self):
        self.assertEqual(sql_campaign("ecc2k-130"), "ecc2k-130")
        self.assertEqual(sql_campaign("ecc2k-130'"), "ecc2k-130''")
        with self.assertRaises(SystemExit):
            sql_campaign("ecc2k-130;drop")
        self.assertIn("'ecc2k-130'", render_sql("x = '{{campaign}}'", "ecc2k-130"))

    def test_rollup_ready_parses_psql_boolean_output(self):
        self.assertTrue(rollup_ready_from_text("t\n"))
        self.assertTrue(rollup_ready_from_text("true"))
        self.assertFalse(rollup_ready_from_text("f"))
        self.assertFalse(rollup_ready_from_text(""))

    def test_fetch_script_punches_runner_ip_not_launch_key(self):
        path = os.path.join(HERE, "fetch_via_walker.sh")
        with open(path, encoding="utf-8") as fh:
            text = fh.read()
        self.assertIn("authorize-security-group-ingress", text)
        self.assertIn("revoke-security-group-ingress", text)
        self.assertIn("checkip.amazonaws.com", text)
        self.assertIn("gha-ecc2k130-status", text)
        self.assertIn("gha_walker.pub", text)
        self.assertIn("InvalidPermission.Duplicate", text)
        self.assertIn("ecc2k-dp-walker", text)
        # The snapshot query runs silently and is getting slower as the DP
        # table grows; without keepalives a slow query is indistinguishable
        # from a dead connection and the hop dies on a broken pipe.
        self.assertIn("ServerAliveInterval", text)
        self.assertIn("ServerAliveCountMax", text)
        self.assertIn("RHO_REMOTE_TIMEOUT:-1200", text)
        self.assertIn("snapshot exceeded", text)
        self.assertIn("exit 124", text)
        self.assertNotIn("meow34", text)
        pub = os.path.join(HERE, "gha_walker.pub")
        with open(pub, encoding="utf-8") as fh:
            line = fh.read().strip()
        self.assertTrue(line.startswith("ssh-ed25519 "))
        self.assertIn("ecc2k130-status-gha", line)

    def test_gha_iam_user_scopes_sg_mutate_to_walker_tags(self):
        path = os.path.join(HERE, "gha_iam_user.sh")
        with open(path, encoding="utf-8") as fh:
            text = fh.read()
        self.assertIn("PunchWalkerSshOnly", text)
        self.assertIn('ec2:ResourceTag/Name": "rho-ecc2k-walker"', text)
        self.assertIn('ec2:ResourceTag/Purpose": "ecc2k-dp-walker"', text)
        self.assertIn("security-group/*", text)
        authorize_star = (
            '"ec2:AuthorizeSecurityGroupIngress"' in text
            and '"Resource": "*"' in text.split("PunchWalkerSshOnly", 1)[1]
        )
        self.assertFalse(authorize_star)

    def test_publish_job_timeout_covers_the_one_time_backfill(self):
        path = os.path.join(ROOT, ".github", "workflows", "ecc2k130-status.yml")
        with open(path, encoding="utf-8") as fh:
            text = fh.read()
        self.assertIn("timeout-minutes: 30", text)
        fetch = os.path.join(HERE, "fetch_via_walker.sh")
        with open(fetch, encoding="utf-8") as fh:
            script = fh.read()
        self.assertIn("RHO_REMOTE_TIMEOUT:-1200", script)
        self.assertIn("SET LOCAL statement_timeout = '1200s'", BACKFILL_SQL)
        self.assertIn("SET LOCAL lock_timeout = '15s'", BACKFILL_SQL)
        with open(os.path.join(HERE, "snapshot.py"), encoding="utf-8") as fh:
            source = fh.read()
        self.assertIn("lost a lock race; retrying", source)


class HistoryTests(unittest.TestCase):
    def test_merges_and_dedups(self):
        previous = {
            "points": [
                {"generated_at": "2026-09-12T09:00:00Z", "dps": 100, "dps_last_hour": 10, "collisions": 0, "workers": 1, "state": "COLLECTING"},
                {"generated_at": "2026-09-12T10:00:00Z", "dps": 200, "dps_last_hour": 20, "collisions": 0, "workers": 1, "state": "COLLECTING"},
            ]
        }
        snapshot = {
            "generated_at": "2026-09-12T10:00:00Z",
            "campaign_id": "ecc2k-130",
            "dps": 210,
            "dps_last_hour": 30,
            "collisions": 0,
            "workers": 2,
            "state": "COLLECTING",
        }
        history = merge_history(previous, snapshot, limit=2)
        self.assertEqual(len(history["points"]), 2)
        self.assertEqual(history["points"][-1]["dps"], 210)
        self.assertEqual(history["campaign_id"], "ecc2k-130")

    def test_render_cli(self):
        import render

        with tempfile.TemporaryDirectory() as tmp:
            status = os.path.join(tmp, "status.json")
            history = os.path.join(tmp, "history.json")
            with open(status, "w", encoding="utf-8") as fh:
                json.dump({
                    "generated_at": "2026-09-12T11:00:00Z",
                    "campaign_id": "ecc2k-130",
                    "dps": 1,
                    "dps_last_hour": 1,
                    "collisions": 0,
                    "workers": 1,
                    "state": "COLLECTING",
                }, fh)
            self.assertEqual(render.main(["--status", status, "--history-out", history]), 0)
            with open(history, encoding="utf-8") as fh:
                blob = json.load(fh)
            self.assertEqual(blob["points"][0]["dps"], 1)

    def test_history_point_carries_the_iteration_total(self):
        snapshot = {
            "generated_at": "2026-09-15T13:00:00Z",
            "campaign_id": "ecc2k-130",
            "dps": 1,
            "dps_last_hour": 1,
            "collisions": 0,
            "workers": 1,
            "state": "COLLECTING",
            "work": {"iterations": 7817055361302528},
        }
        point = merge_history({}, snapshot)["points"][-1]
        self.assertEqual(point["iterations"], 7817055361302528)
        del snapshot["work"]
        self.assertIsNone(merge_history({}, snapshot)["points"][-1]["iterations"])


def rate_history(samples):
    return [
        {
            "generated_at": at,
            "dps": 0,
            "dps_last_hour": 0,
            "collisions": 0,
            "workers": 1,
            "iterations": iterations,
            "state": "COLLECTING",
        }
        for at, iterations in samples
    ]


class WalkRateTests(unittest.TestCase):
    # The measured rate is the point of the feature, so these are the numbers
    # it was verified against: four consecutive samples of the live campaign
    # feed on 2026-09-15, whose checkpoint totals differ by 109,081,968,771,072
    # iterations over 1504 s = 72.5 B it/s. status.py's own sum of what the
    # walkers reported themselves over the same period was 86-101 B it/s; the
    # difference is each slot's walked-but-not-yet-checkpointed tail, which
    # this figure deliberately does not count.
    LIVE = (
        ("2026-09-15T13:15:37Z", 7817055361302528),
        ("2026-09-15T13:25:26Z", 7860096436535296),
        ("2026-09-15T13:33:51Z", 7888628575371264),
        ("2026-09-15T13:40:41Z", 7926137330073600),
    )

    def test_measures_the_live_campaign_samples(self):
        rate = measure_rate(rate_history(self.LIVE))
        self.assertEqual(rate["window_seconds"], 1504)
        self.assertAlmostEqual(rate["iterations_per_second"] / 1e9, 72.528, places=2)
        self.assertEqual(rate["measured_from"], "2026-09-15T13:15:37Z")
        self.assertEqual(rate["measured_to"], "2026-09-15T13:40:41Z")
        self.assertEqual(rate["iterations_to"] - rate["iterations_from"], 109081968771072)

    def test_window_bounds_the_smoothing_not_the_reach(self):
        # Inside the window the oldest sample wins, so the difference spans
        # about an hour of checkpoints rather than one publish interval.
        samples = [("2026-09-15T%02d:00:00Z" % hour, 1000 + hour) for hour in range(10, 15)]
        rate = measure_rate(rate_history(samples), window_s=RATE_WINDOW_S)
        self.assertEqual(rate["window_seconds"], RATE_WINDOW_S)
        self.assertEqual(rate["measured_from"], "2026-09-15T13:00:00Z")

    def test_reaches_past_the_window_after_an_outage(self):
        # One sample inside the window and one three hours before it: the
        # difference is still worth taking, and the span it used says so.
        rate = measure_rate(rate_history([
            ("2026-09-15T10:00:00Z", 1_000_000_000_000),
            ("2026-09-15T13:00:00Z", 1_360_000_000_000),
        ]))
        self.assertEqual(rate["window_seconds"], 3 * 3600)
        self.assertAlmostEqual(rate["iterations_per_second"], 360e9 / (3 * 3600))

    def test_declines_to_guess(self):
        first = self.LIVE[:1]
        # One sample has nothing to difference against.
        self.assertIsNone(measure_rate(rate_history(first)))
        # Neither does a history whose points predate the iteration total.
        pointless = rate_history(self.LIVE)
        for point in pointless:
            point.pop("iterations")
        self.assertIsNone(measure_rate(pointless))
        # Two samples inside the checkpoint granularity are not a measurement.
        self.assertIsNone(measure_rate(rate_history([
            ("2026-09-15T13:00:00Z", 1_000_000_000_000),
            ("2026-09-15T13:05:00Z", 1_000_500_000_000),
        ])))
        self.assertLess(5 * 60, MIN_RATE_SPAN_S)
        # A total that went backwards is a reset checkpoint, not a negative
        # rate: a slot restarted from an earlier base contributes less than it
        # did, and the page says it has no rate rather than showing that.
        self.assertIsNone(measure_rate(rate_history([
            ("2026-09-15T12:00:00Z", 2_000_000_000_000),
            ("2026-09-15T13:00:00Z", 1_000_000_000_000),
        ])))
        # Nothing within reach of the newest sample at all.
        self.assertIsNone(measure_rate(rate_history([
            ("2026-09-14T13:00:00Z", 1_000_000_000_000),
            ("2026-09-15T13:00:00Z", 2_000_000_000_000),
        ]), max_span_s=MAX_RATE_SPAN_S))

    def test_a_halted_fleet_measures_zero_and_says_so(self):
        # Zero is a measurement: the totals stopped moving. It must not become
        # None (which the page renders as "no rate yet") or an ETA.
        rate = measure_rate(rate_history([
            ("2026-09-15T12:00:00Z", 2_000_000_000_000),
            ("2026-09-15T13:00:00Z", 2_000_000_000_000),
        ]))
        self.assertEqual(rate["iterations_per_second"], 0.0)

    def test_apply_rate_sets_and_clears(self):
        snapshot = {"campaign_id": "ecc2k-130"}
        self.assertTrue(apply_rate(snapshot, {"iterations_per_second": 1.0}))
        self.assertEqual(snapshot["walk_rate"]["iterations_per_second"], 1.0)
        self.assertFalse(apply_rate(snapshot, None))
        self.assertNotIn("walk_rate", snapshot)

    def test_render_cli_publishes_the_rate_into_the_snapshot(self):
        import render

        with tempfile.TemporaryDirectory() as tmp:
            status = os.path.join(tmp, "status.json")
            history_in = os.path.join(tmp, "history-prev.json")
            history_out = os.path.join(tmp, "history.json")
            with open(history_in, "w", encoding="utf-8") as fh:
                json.dump({"points": rate_history(self.LIVE[:-1])}, fh)
            with open(status, "w", encoding="utf-8") as fh:
                json.dump({
                    "generated_at": self.LIVE[-1][0],
                    "campaign_id": "ecc2k-130",
                    "dps": 31635015,
                    "dps_last_hour": 1,
                    "collisions": 0,
                    "workers": 1,
                    "state": "COLLECTING",
                    "work": {"iterations": self.LIVE[-1][1]},
                }, fh)
            self.assertEqual(render.main([
                "--status", status,
                "--history-in", history_in,
                "--history-out", history_out,
                "--status-out", status,
            ]), 0)
            with open(status, encoding="utf-8") as fh:
                published = json.load(fh)
            self.assertAlmostEqual(
                published["walk_rate"]["iterations_per_second"] / 1e9, 72.528, places=2)
            assert_public(published)

    def test_render_cli_does_not_publish_a_historical_rate_without_work(self):
        import render

        with tempfile.TemporaryDirectory() as tmp:
            status = os.path.join(tmp, "status.json")
            history_in = os.path.join(tmp, "history-prev.json")
            history_out = os.path.join(tmp, "history.json")
            with open(history_in, "w", encoding="utf-8") as fh:
                json.dump({"points": rate_history(self.LIVE[:-1])}, fh)
            with open(status, "w", encoding="utf-8") as fh:
                json.dump({
                    "generated_at": self.LIVE[-1][0],
                    "campaign_id": "ecc2k-130",
                    "dps": 31635015,
                    "dps_last_hour": 1,
                    "collisions": 0,
                    "workers": 1,
                    "state": "COLLECTING",
                }, fh)
            self.assertEqual(render.main([
                "--status", status,
                "--history-in", history_in,
                "--history-out", history_out,
                "--status-out", status,
            ]), 0)
            with open(status, encoding="utf-8") as fh:
                published = json.load(fh)
            self.assertNotIn("walk_rate", published)


def work_feed(iterations=7817055361302528, generated_at=None, campaign="ecc2k-130"):
    generated_at = generated_at or datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    return {
        "campaign_id": campaign,
        "generated_at": generated_at,
        "source": "dp_ingest.py (live, direct from rho-dp + slot checkpoints)",
        "work": {
            "iterations": iterations,
            "method": "sum over slot checkpoints of iterBase * threads * batch",
            "per_slot": [
                {"slot": 2, "iterations": 100, "checkpoint_age_s": 206, "retired": False},
                {"slot": 3, "iterations": 100, "checkpoint_age_s": 44, "retired": False},
                {"slot": 4, "iterations": 100, "checkpoint_age_s": 90000, "retired": True},
                {"slot": 5, "iterations": 100, "checkpoint_age_s": 9000, "retired": False},
            ],
        },
    }


class WorkFeedTests(unittest.TestCase):
    def test_takes_the_iteration_total_and_counts_who_is_walking(self):
        block = work_block(work_feed(), "ecc2k-130")
        self.assertEqual(block["iterations"], 7817055361302528)
        self.assertAlmostEqual(block["iterations_log2"], 52.796, places=2)
        self.assertEqual(block["slots"], 4)
        # Retired, and stale by three missed checkpoints, are both not walking.
        self.assertEqual(block["walking_slots"], 2)

    def test_merges_into_a_snapshot_without_publishing_anything_private(self):
        snapshot = {"campaign_id": "ecc2k-130", "dps": 1}
        self.assertTrue(merge_work(snapshot, work_feed(), "ecc2k-130"))
        self.assertEqual(snapshot["work"]["iterations"], 7817055361302528)
        assert_public(snapshot)
        # No per-slot rows: the page needs the total and a count, and each row
        # names a walker's slot and run.
        self.assertNotIn("per_slot", snapshot["work"])

    def test_ingest_health_overrides_a_quiet_store_when_the_feed_is_behind(self):
        feed = work_feed()
        feed["ingest"] = {
            "outstanding_objects": 42,
            "unrecognised_objects": 0,
            "newest_object_at": "2026-09-18T11:00:00Z",
            "lag_seconds": 90,
        }
        snapshot = {"campaign_id": "ecc2k-130", "dps": 10, "state": "IDLE_OR_STALE"}
        self.assertTrue(merge_work(snapshot, feed, "ecc2k-130"))
        self.assertEqual(snapshot["state"], "INGEST_BEHIND")
        self.assertEqual(snapshot["ingest"]["outstanding_objects"], 42)
        assert_public(snapshot)
        self.assertIsNone(ingest_block("not a feed"))
        self.assertIsNone(ingest_block({"campaign_id": "ecc2k-130"}))

    def test_collecting_is_not_relabelled_as_ingest_behind(self):
        feed = work_feed()
        feed["ingest"] = {"outstanding_objects": 1, "unrecognised_objects": 0}
        snapshot = {"campaign_id": "ecc2k-130", "state": "COLLECTING"}
        self.assertTrue(merge_work(snapshot, feed, "ecc2k-130"))
        self.assertEqual(snapshot["state"], "COLLECTING")
        self.assertEqual(snapshot["ingest"]["outstanding_objects"], 1)

    def test_refuses_a_stale_feed_rather_than_freezing_the_total(self):
        # A frozen total does not read as a dead ingest host, it reads as a
        # dead campaign, and it drags the measured rate to zero as the
        # snapshots catch up with it.
        old = datetime.now(timezone.utc) - timedelta(seconds=MAX_FEED_AGE_S + 60)
        feed = work_feed(generated_at=old.strftime("%Y-%m-%dT%H:%M:%SZ"))
        self.assertIsNone(work_block(feed, "ecc2k-130"))
        snapshot = {"campaign_id": "ecc2k-130"}
        self.assertFalse(merge_work(snapshot, feed, "ecc2k-130"))
        self.assertNotIn("work", snapshot)

    def test_refuses_another_campaign_and_unusable_totals(self):
        self.assertIsNone(work_block(work_feed(campaign="eccp-131"), "ecc2k-130"))
        self.assertIsNone(work_block(work_feed(iterations=0), "ecc2k-130"))
        self.assertIsNone(work_block({"campaign_id": "ecc2k-130"}, "ecc2k-130"))
        self.assertIsNone(work_block({"campaign_id": "ecc2k-130", "work": {}}, "ecc2k-130"))
        self.assertIsNone(work_block("not a feed", "ecc2k-130"))
        no_time = work_feed()
        del no_time["generated_at"]
        self.assertIsNone(work_block(no_time, "ecc2k-130"))

    def test_feed_url_is_derived_and_never_hardcoded(self):
        # The bucket is named for the AWS account, which this public tree does
        # not carry; the URL is built from the calling identity at run time.
        url = feed_url("ecc2k130", "us-west-2", account="123456789012")
        self.assertEqual(
            url, "https://ecc2k130-status-123456789012.s3.us-west-2.amazonaws.com/status.json")
        self.assertTrue(url.startswith("https://"))
        with open(os.path.join(HERE, "work_feed.py"), encoding="utf-8") as fh:
            source = fh.read()
        self.assertNotRegex(source, r"\b\d{12}\b", "an account number is hardcoded in work_feed.py")

    def test_workflow_merges_the_feed_before_it_renders_history(self):
        # The rate is a difference between snapshots, so the iteration total
        # has to be in the snapshot before render.py records the point. In the
        # other order every point publishes without a total and no rate is
        # ever measured.
        path = os.path.join(ROOT, ".github", "workflows", "ecc2k130-status.yml")
        with open(path, encoding="utf-8") as fh:
            text = fh.read()
        self.assertIn("--status-out docs/ecc2k130-status/status.json", text)
        self.assertLess(
            text.index("python3 scripts/rho_status/work_feed.py"),
            text.index("python3 scripts/rho_status/render.py"),
        )


if __name__ == "__main__":
    unittest.main()
