#!/usr/bin/env python3
"""Offline tests for the ECC2K-130 status snapshot and history merge."""

from __future__ import annotations

import http.server
import json
import os
import re
import shutil
import socket
import stat
import subprocess
import tempfile
import threading
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
    stamp_published,
)
from snapshot import (
    BACKFILL_SQL,
    BOUNDS_SQL,
    CAMPAIGN_SQL,
    CLAIM_BOUNDARY,
    ENSURE_SQL,
    FALLBACK_SQL,
    FORBIDDEN_PUBLIC_KEYS,
    HOUR_COUNT_SQL,
    MARK_READY_SQL,
    READ_SQL,
    SESSION_PREAMBLE,
    assert_public,
    campaign_state,
    hours_to_backfill,
    normalize,
    parse_bounds,
    parse_hour_counts,
    render_sql,
    rollup_ready_from_text,
    sql_campaign,
    sql_hour,
)
from check_feed_age import describe as describe_feed_age
from work_feed import (
    MAX_FEED_AGE_S,
    feed_url,
    ingest_block,
    merge_work,
    snapshot_from_feed,
    work_block,
)


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

    def test_publish_job_does_not_depend_on_aws_or_the_walker_to_publish(self):
        # 2026-09-18T19:15Z to 2026-09-19T04:41Z: every scheduled publish
        # died with no running walker and Pages froze on 11:25Z. The feed is
        # the floor the job always has; the AWS step and the hop are how it
        # does better, and neither may end the job when it fails.
        path = os.path.join(ROOT, ".github", "workflows", "ecc2k130-status.yml")
        with open(path, encoding="utf-8") as fh:
            text = fh.read()
        aws_step = text.split("- name: Configure AWS", 1)[1].split("- name:", 1)[0]
        self.assertIn("continue-on-error: true", aws_step)
        self.assertIn("aws-actions/configure-aws-credentials", aws_step)
        # The feed URL secret lets the feed be read with no AWS identity, in
        # the snapshot step and again in the iteration-total merge.
        self.assertGreaterEqual(
            text.count("RHO_WORK_FEED_URL: ${{ secrets.RHO_WORK_FEED_URL }}"), 2)
        snapshot_step = text.split("- name: Snapshot", 1)[1].split("- name:", 1)[0]
        self.assertIn("fetch_via_walker.sh", snapshot_step)
        self.assertIn("RHO_WORK_FEED_URL: ${{ secrets.RHO_WORK_FEED_URL }}", snapshot_step)
        # Losing the live history.json loses seven days of points and the
        # measured rate; one Pages hiccup must not be allowed to.
        history_step = text.split("PAGES_HISTORY_URL:", 1)[1].split("render.py", 1)[0]
        self.assertIn("--retry", history_step)
        self.assertIn("--retry-all-errors", history_step)
        self.assertIn("--max-time", history_step)

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
        # When the walker hop is down the snapshot is the ingest host's feed,
        # rewritten every --status-every seconds; the banner must not call a
        # document stale that is simply one feed interval old at publish time
        # plus one cron interval by the next publish.
        ingest = os.path.join(ROOT, "ecc2k130", "aws", "dp_ingest.py")
        with open(ingest, encoding="utf-8") as fh:
            source = fh.read()
        every = re.search(r'RHO_STATUS_EVERY", "([\d.]+)"', source)
        if every is None:
            every = re.search(r'"--status-every", type=float, default=([\d.]+)', source)
        self.assertIsNotNone(every, "--status-every default not found in dp_ingest.py")
        feed_minutes = float(every.group(1)) / 60
        self.assertGreaterEqual(
            stale_minutes, feed_minutes + step,
            "a feed-sourced snapshot is %d min old at publish and %d by the next; "
            "STALE_AFTER_MS of %d min will flap" % (feed_minutes, feed_minutes + step, stale_minutes))
        self.assertGreaterEqual(stale_minutes, 2 * feed_minutes, "one missed feed write is not stale")
        # Two clocks on the page: the source's generated_at and the
        # publisher's published_at, so a stopped publisher and a stopped
        # source do not share one sentence.
        self.assertIn("status.published_at", text)
        self.assertIn("Published to Pages", text)
        self.assertIn("The publisher ran ", text)
        self.assertIn("neither the walker hop nor the ingest feed had anything newer", text)
        self.assertIn("INGEST_BEHIND", text)
        self.assertIn(
            "The fleet is walking, but the store is behind on ingest",
            text,
        )
        # The banner's "last hour" is wall-clock. A frozen snapshot that still
        # says COLLECTING is what made the page claim new points while also
        # saying the job had not run in 16 hours.
        self.assertIn("function campaignDisplayState", text)
        self.assertIn("This snapshot is ", text)
        self.assertNotIn("The publishing job last ran ", text)
        self.assertIn("Per-worker counts need the walker hop", text)

    def test_snapshot_reads_the_rollup_and_scans_the_heap_only_to_backfill(self):
        # Three separate LATERAL scans of a 30M-row table is what pushed the
        # walker hop from 41s to 510s and broke publication; grouping the heap
        # once by worker and hour then pushed it to 524s at 187M rows, one
        # timeout tick under 540s, and the next scheduled run died. The
        # ready document is assembled from rho_dp_hour / rho_dp_recent.
        # distinguished_points is heap-fetched only by one hour of backfill
        # at a time; the publish-first fallback groups by hour, not worker.
        self.assertEqual(BACKFILL_SQL.count("FROM distinguished_points"), 1)
        self.assertIn("ISOLATION LEVEL REPEATABLE READ", BACKFILL_SQL)
        self.assertNotIn("LOCK TABLE", BACKFILL_SQL)
        self.assertIn("date_trunc('hour', found_at)", BACKFILL_SQL)
        self.assertIn("ON CONFLICT (campaign_id, hour, worker_id)", BACKFILL_SQL)
        self.assertIn("found_at >= '{{hour}}'::timestamptz", BACKFILL_SQL)
        self.assertIn("found_at < '{{hour}}'::timestamptz + interval '1 hour'", BACKFILL_SQL)
        recent = BACKFILL_SQL.split("INSERT INTO rho_dp_recent", 1)[1].split("UPDATE rho_dp_meta", 1)[0]
        self.assertIn("worker_id, found_at, count(*)", recent)
        self.assertNotIn("worker_id, hour, dps", recent)
        self.assertNotIn("FROM distinguished_points", READ_SQL)
        self.assertIn("FROM rho_dp_hour", READ_SQL)
        self.assertIn("FROM rho_dp_recent", READ_SQL)
        self.assertIn("CREATE TRIGGER", ENSURE_SQL)
        self.assertIn("rho_dp_rollup_insert", ENSURE_SQL)
        self.assertIn("ON distinguished_points", ENSURE_SQL)
        self.assertIn("backfill_through", ENSURE_SQL)
        # Ingest INSERTs fire this as the rho/dp-rds role, which is not the
        # walker DATABASE_URL role that creates the rollup tables.
        self.assertIn("SECURITY DEFINER", ENSURE_SQL)
        self.assertGreaterEqual(BACKFILL_SQL.count("{{campaign}}"), 1)
        self.assertIn("{{campaign}}", READ_SQL)
        with open(os.path.join(HERE, "snapshot.py"), encoding="utf-8") as fh:
            source = fh.read()
        self.assertNotIn('replace("%(campaign)s", "%s")', source)
        self.assertIn("if rollup_ready_from_text", source)
        self.assertIn("publishing index-only fallback first", source)
        self.assertIn("SET statement_timeout = 0;", SESSION_PREAMBLE)

    def test_fallback_counts_one_hour_without_worker_id(self):
        # The 18:56Z publish died ~200s into a full-table GROUP BY hour:
        # walker SSH reset and Pages stayed on 11:25Z. Fallback is a count
        # of one hour through the found_at index, no worker_id, no GROUP BY.
        self.assertEqual(FALLBACK_SQL, HOUR_COUNT_SQL)
        self.assertIn("FROM distinguished_points", HOUR_COUNT_SQL)
        self.assertNotIn("worker_id", HOUR_COUNT_SQL)
        self.assertNotIn("GROUP BY", HOUR_COUNT_SQL)
        self.assertIn("found_at >= '{{hour}}'::timestamptz", HOUR_COUNT_SQL)
        self.assertIn("found_at < '{{hour}}'::timestamptz + interval '1 hour'", HOUR_COUNT_SQL)
        self.assertIn("SET statement_timeout = '60s'", HOUR_COUNT_SQL)
        self.assertIn("rho_campaigns", CAMPAIGN_SQL)
        self.assertNotIn("distinguished_points", CAMPAIGN_SQL)
        self.assertIn("backfill_through", MARK_READY_SQL)
        self.assertIn("min(found_at)", BOUNDS_SQL)
        self.assertEqual(
            parse_hour_counts("123|2026-09-18 11:00:00+00|2026-09-18 11:59:00+00|4|50"),
            (123, "2026-09-18 11:00:00+00", "2026-09-18 11:59:00+00", 4, 50),
        )
        self.assertEqual(parse_hour_counts("0||||"), (0, None, None, 0, 0))

    def test_backfill_resumes_after_the_hour_already_written(self):
        start = datetime(2026, 9, 11, 17, tzinfo=timezone.utc)
        now = datetime(2026, 9, 11, 20, tzinfo=timezone.utc)
        self.assertEqual(
            [h.hour for h in hours_to_backfill(start, now, None)],
            [17, 18, 19, 20],
        )
        through = datetime(2026, 9, 11, 18, tzinfo=timezone.utc)
        self.assertEqual(
            [h.hour for h in hours_to_backfill(start, now, through)],
            [19, 20],
        )
        self.assertEqual(hours_to_backfill(start, now, now), [])
        start_h, now_h, through_h = parse_bounds(
            "2026-09-11T17:00:00Z|2026-09-18T18:00:00Z|2026-09-12T00:00:00Z"
        )
        self.assertEqual(start_h, datetime(2026, 9, 11, 17, tzinfo=timezone.utc))
        self.assertEqual(now_h, datetime(2026, 9, 18, 18, tzinfo=timezone.utc))
        self.assertEqual(through_h, datetime(2026, 9, 12, 0, tzinfo=timezone.utc))
        self.assertEqual(sql_hour(start_h), "2026-09-11T17:00:00Z")
        with self.assertRaises(SystemExit):
            sql_hour("2026-09-11 17:00:00")
        self.assertIn("'2026-09-11T17:00:00Z'", render_sql(
            "x = '{{hour}}'::timestamptz", "ecc2k-130", hour=start_h
        ))

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
        self.assertIn("no running instance tagged Name=rho-ecc2k-walker", text)
        self.assertIn("tagged walkers (any state)", text)
        # The snapshot query runs silently and is getting slower as the DP
        # table grows; without keepalives a slow query is indistinguishable
        # from a dead connection and the hop dies on a broken pipe.
        self.assertIn("ServerAliveInterval", text)
        self.assertIn("ServerAliveCountMax", text)
        self.assertIn("PYTHONUNBUFFERED=1", text)
        self.assertIn("ServerAliveInterval=5", text)
        self.assertIn("RHO_REMOTE_TIMEOUT:-1500", text)
        self.assertIn("snapshot exceeded", text)
        self.assertIn("GITHUB_RUN_ID", text)
        self.assertIn("REMOTE_JSON=", text)
        # A 124 after the fallback write must still copy this run's file, and
        # must not copy /tmp/rho_status.json leftover from the 11:25Z hop.
        self.assertNotIn('"$USER@$HOST:/tmp/rho_status.json"', text)
        self.assertIn('"$USER@$HOST:$REMOTE_JSON" "$out"', text)
        self.assertNotIn("meow34", text)
        pub = os.path.join(HERE, "gha_walker.pub")
        with open(pub, encoding="utf-8") as fh:
            line = fh.read().strip()
        self.assertTrue(line.startswith("ssh-ed25519 "))
        self.assertIn("ecc2k130-status-gha", line)

    def test_fetch_script_writes_the_feed_first_and_hops_in_a_subshell(self):
        # The order and the shape are the fix. Every earlier repair hardened
        # one step of the hop and the next step froze Pages anyway, because
        # the hop ran first under `set -e` and any failure ended the script
        # before the fallback. Now the feed copy exists before the hop is
        # tried, the hop is a subshell that returns instead of exiting, and
        # the script only exits 1 when both sources came up empty.
        path = os.path.join(HERE, "fetch_via_walker.sh")
        with open(path, encoding="utf-8") as fh:
            text = fh.read()
        self.assertIn("--as-snapshot", text)
        feed_at = text.index("write_from_ingest_feed\nFEED_RC=$?")
        hop_def_at = text.index("walker_hop() (")
        hop_call_at = text.index('walker_hop "$HOP_OUT"')
        self.assertLess(feed_at, hop_def_at)
        self.assertLess(hop_def_at, hop_call_at)
        body = text[hop_def_at:hop_call_at]
        # A subshell body with its own strict mode, an EXIT trap that revokes
        # the punch-hole on every path out, and no `exit` anywhere inside.
        self.assertIn("set -euo pipefail", body)
        self.assertIn("trap cleanup EXIT", body)
        self.assertNotRegex(body, r"^\s*exit \d", "the hop must return, never exit")
        self.assertGreaterEqual(body.count("return 1"), 5)
        # Skips that keep the feed: no key, or asked not to.
        self.assertIn("RHO_SKIP_WALKER", text)
        self.assertIn('[ ! -s "$KEY" ]', text)
        # The hop wins when it produced parseable JSON; the feed otherwise;
        # nothing at all is the only exit 1.
        tail = text[hop_call_at:]
        self.assertIn('mv -f "$HOP_OUT" "$OUT"', tail)
        self.assertIn('mv -f "$FEED_OUT" "$OUT"', tail)
        self.assertLess(tail.index('mv -f "$HOP_OUT" "$OUT"'), tail.index('mv -f "$FEED_OUT" "$OUT"'))
        self.assertEqual(tail.count("exit 1"), 1)
        self.assertIn("neither the ingest status feed nor the walker hop", tail)

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
        self.assertIn("timeout-minutes: 35", text)
        fetch = os.path.join(HERE, "fetch_via_walker.sh")
        with open(fetch, encoding="utf-8") as fh:
            script = fh.read()
        self.assertIn("RHO_REMOTE_TIMEOUT:-1500", script)
        self.assertIn("SET LOCAL statement_timeout = '300s'", BACKFILL_SQL)
        self.assertIn("SET LOCAL lock_timeout = '15s'", BACKFILL_SQL)
        with open(os.path.join(HERE, "snapshot.py"), encoding="utf-8") as fh:
            source = fh.read()
        self.assertIn("lost a lock race; retrying", source)
        self.assertIn("RHO_SNAPSHOT_BUDGET_S", source)


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


class PublishedStampTests(unittest.TestCase):
    def test_stamp_is_the_publish_time_not_the_source_time(self):
        snapshot = {"campaign_id": "ecc2k-130", "generated_at": "2026-09-19T04:26:48Z"}
        now = datetime(2026, 9, 19, 4, 42, 3, 512000, tzinfo=timezone.utc)
        self.assertEqual(stamp_published(snapshot, now), "2026-09-19T04:42:03Z")
        self.assertEqual(snapshot["published_at"], "2026-09-19T04:42:03Z")
        self.assertEqual(snapshot["generated_at"], "2026-09-19T04:26:48Z")

    def test_render_cli_stamps_the_document_it_publishes(self):
        import render

        with tempfile.TemporaryDirectory() as tmp:
            status = os.path.join(tmp, "status.json")
            history = os.path.join(tmp, "history.json")
            with open(status, "w", encoding="utf-8") as fh:
                json.dump({"campaign_id": "ecc2k-130", "generated_at": "2026-09-19T04:26:48Z",
                           "dps": 5}, fh)
            before = datetime.now(timezone.utc) - timedelta(seconds=2)
            self.assertEqual(render.main([
                "--status", status, "--history-out", history, "--status-out", status,
            ]), 0)
            with open(status, encoding="utf-8") as fh:
                published = json.load(fh)
            self.assertEqual(published["generated_at"], "2026-09-19T04:26:48Z")
            stamp = datetime.fromisoformat(published["published_at"].replace("Z", "+00:00"))
            self.assertGreaterEqual(stamp, before.replace(microsecond=0))
            # History points key on generated_at; the publish stamp is not one
            # of them, so republishing the same source document does not
            # multiply history rows.
            with open(history, encoding="utf-8") as fh:
                points = json.load(fh)["points"]
            self.assertEqual(len(points), 1)
            self.assertNotIn("published_at", points[0])


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


def ingest_campaign_feed(generated_at=None, campaign="ecc2k-130",
                         dps=187000000, dps_last_hour=5120):
    generated_at = generated_at or datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    feed = work_feed(generated_at=generated_at, campaign=campaign)
    feed.update({
        "curve_id": "certicom-ecc2k-130",
        "dp_mask_bits": None,
        "campaign_created_at": "2026-09-11T17:51:35Z",
        "dps": dps,
        "dps_last_hour": dps_last_hour,
        "dps_last_day": 120000,
        "collisions": 0,
        "walkers": 13,
        "first_dp_at": "2026-09-11T17:51:35Z",
        "last_dp_at": generated_at,
        "hourly": [{"hour": "2026-09-19T03:00:00Z", "dps": dps_last_hour}],
        "ingest": {
            "outstanding_objects": 0,
            "unrecognised_objects": 0,
            "newest_object_at": generated_at,
            "lag_seconds": 12,
        },
    })
    return feed


class WorkFeedTests(unittest.TestCase):
    def test_takes_the_iteration_total_and_counts_who_is_walking(self):
        block = work_block(work_feed(), "ecc2k-130")
        self.assertEqual(block["iterations"], 7817055361302528)
        self.assertAlmostEqual(block["iterations_log2"], 52.796, places=2)
        self.assertEqual(block["slots"], 4)
        # Retired, and stale by three missed checkpoints, are both not walking.
        self.assertEqual(block["walking_slots"], 2)

    def test_takes_walking_slots_when_the_feed_has_already_stripped_per_slot(self):
        feed = work_feed()
        del feed["work"]["per_slot"]
        feed["work"]["slots"] = 4
        feed["work"]["walking_slots"] = 2
        block = work_block(feed, "ecc2k-130")
        self.assertEqual(block["walking_slots"], 2)
        self.assertEqual(block["slots"], 4)

    def test_counts_slots_off_the_campaign_cutoff_from_the_rows(self):
        # dp_ingest marks each slot against the campaign's weight; a slot at
        # another cutoff mostly cannot collide with the rest, so the page gets
        # it as a separate count and never folded into the walkers.
        feed = work_feed()
        rows = feed["work"]["per_slot"]
        rows[0]["dp_weight_ok"] = True
        rows[1]["dp_weight_ok"] = False      # walking, off weight
        rows[2]["dp_weight_ok"] = False      # retired, off weight
        rows[3]["dp_weight_ok"] = None       # too few records to judge
        feed["work"]["campaign_dp_weight"] = 32
        feed["work"]["iterations_per_dp_log2_expected"] = 28.41
        block = work_block(feed, "ecc2k-130")
        self.assertEqual(block["off_weight_slots"], 2)
        self.assertEqual(block["off_weight_walking_slots"], 1)
        self.assertEqual(block["campaign_dp_weight"], 32)
        self.assertEqual(block["iterations_per_dp_log2_expected"], 28.41)

    def test_takes_the_off_weight_counts_when_the_feed_has_stripped_per_slot(self):
        feed = work_feed()
        del feed["work"]["per_slot"]
        feed["work"].update(slots=10, walking_slots=10, off_weight_slots=5,
                            off_weight_walking_slots=4)
        block = work_block(feed, "ecc2k-130")
        self.assertEqual(block["off_weight_slots"], 5)
        self.assertEqual(block["off_weight_walking_slots"], 4)
        # An older ingest host that does not publish them reads as zero, not
        # as a missing key the page has to special-case.
        del feed["work"]["off_weight_slots"]
        del feed["work"]["off_weight_walking_slots"]
        block = work_block(feed, "ecc2k-130")
        self.assertEqual(block["off_weight_slots"], 0)
        self.assertNotIn("campaign_dp_weight", block)

    def test_re_reports_are_copied_from_the_ingest_block(self):
        feed = work_feed()
        feed["ingest"] = {
            "outstanding_objects": 0, "unrecognised_objects": 0,
            "newest_object_at": "2026-09-21T00:00:00Z", "lag_seconds": 30,
            "duplicate_records": 813112, "duplicate_records_last_day": 4000,
            "duplicate_slots_last_day": 2,
        }
        block = ingest_block(feed)
        self.assertEqual(block["duplicate_records"], 813112)
        self.assertEqual(block["duplicate_records_last_day"], 4000)
        self.assertEqual(block["duplicate_slots_last_day"], 2)
        snapshot = {"campaign_id": "ecc2k-130", "dps": 1, "state": "COLLECTING"}
        self.assertTrue(merge_work(snapshot, feed, "ecc2k-130"))
        self.assertEqual(snapshot["ingest"]["duplicate_records_last_day"], 4000)
        assert_public(snapshot)
        # Absent on an older feed: absent here too, so the page can say so.
        self.assertNotIn("duplicate_records", ingest_block(ingest_campaign_feed()))

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

    def test_snapshot_from_feed_is_public_and_keeps_the_ingest_timestamp(self):
        stamped = "2026-09-19T03:44:00Z"
        feed = ingest_campaign_feed(generated_at=stamped)
        snapshot = snapshot_from_feed(feed, "ecc2k-130")
        self.assertEqual(snapshot["generated_at"], stamped)
        self.assertEqual(snapshot["source"], "ingest-status-feed")
        self.assertEqual(snapshot["query_driver"], "ingest-feed")
        self.assertEqual(snapshot["state"], "COLLECTING")
        self.assertEqual(snapshot["dps"], 187000000)
        self.assertEqual(snapshot["dps_last_hour"], 5120)
        self.assertEqual(snapshot["per_worker"], [])
        self.assertEqual(snapshot["workers"], 0)
        self.assertEqual(snapshot["work"]["walking_slots"], 2)
        self.assertNotIn("per_slot", snapshot["work"])
        assert_public(snapshot)
        dumped = json.dumps(snapshot)
        self.assertNotIn("run_id", dumped)
        self.assertNotIn('"slot"', dumped)

    def test_snapshot_from_feed_keeps_work_on_a_stale_ingest_document(self):
        # Overlaying a stale total onto a fresh walker snapshot is refused;
        # here the ingest document *is* the snapshot, so dropping the total
        # would publish counts with no walk rate for no reason.
        old = datetime.now(timezone.utc) - timedelta(seconds=MAX_FEED_AGE_S + 3600)
        stamped = old.strftime("%Y-%m-%dT%H:%M:%SZ")
        feed = ingest_campaign_feed(generated_at=stamped)
        snapshot = snapshot_from_feed(feed, "ecc2k-130")
        self.assertEqual(snapshot["generated_at"], stamped)
        self.assertEqual(snapshot["work"]["iterations"], 7817055361302528)
        self.assertGreater(snapshot["work"]["feed_age_seconds"], MAX_FEED_AGE_S)
        self.assertIsNone(snapshot_from_feed(work_feed(), "ecc2k-130"))
        self.assertIsNone(snapshot_from_feed(ingest_campaign_feed(campaign="other"), "ecc2k-130"))

    def test_as_snapshot_cli_writes_the_file_and_fails_closed(self):
        import work_feed as wf

        feed = ingest_campaign_feed(generated_at="2026-09-19T03:44:00Z")
        original = wf.fetch
        wf.fetch = lambda url, timeout=20: feed
        try:
            with tempfile.TemporaryDirectory() as tmp:
                out = os.path.join(tmp, "status.json")
                self.assertEqual(wf.main([
                    "--as-snapshot",
                    "--status", out,
                    "--feed-url", "https://example.invalid/status.json",
                ]), 0)
                with open(out, encoding="utf-8") as fh:
                    published = json.load(fh)
                self.assertEqual(published["generated_at"], "2026-09-19T03:44:00Z")
                self.assertEqual(published["dps"], 187000000)
                assert_public(published)
                wf.fetch = lambda url, timeout=20: work_feed()
                self.assertEqual(wf.main([
                    "--as-snapshot",
                    "--status", out,
                    "--feed-url", "https://example.invalid/status.json",
                ]), 1)
        finally:
            wf.fetch = original

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


FAKE_AWS = r"""#!/usr/bin/env bash
# Fake `aws` for the fetch-script tests. FAKE_AWS_MODE: fail | nowalker | walker.
echo "aws $*" >> "$FAKE_LOG"
case "$FAKE_AWS_MODE" in
  fail)
    echo "Unable to locate credentials. You can configure credentials by running \"aws configure\"." >&2
    exit 255 ;;
esac
if [[ "$*" == *"sts get-caller-identity"* ]]; then echo 123456789012; exit 0; fi
if [[ "$*" == *"describe-instances"* ]]; then
  if [[ "$FAKE_AWS_MODE" == walker && "$*" == *"instance-state-name,Values=running"* ]]; then
    printf 'i-0fake\t203.0.113.7\tsg-0fake\n'
  elif [[ "$*" == *"instance-state-name,Values=running"* ]]; then
    printf 'None\tNone\tNone\n'
  else
    printf 'i-0fake\tstopped\tNone\t2026-09-18T19:15:00Z\n'
  fi
  exit 0
fi
if [[ "$*" == *"describe-security-groups"* ]]; then
  if [[ "$*" == *"Purpose"* ]]; then echo ecc2k-dp-walker; else echo rho-ecc2k-walker; fi
  exit 0
fi
if [[ "$*" == *"authorize-security-group-ingress"* ]]; then exit 0; fi
if [[ "$*" == *"revoke-security-group-ingress"* ]]; then exit 0; fi
echo "fake aws: unexpected $*" >&2
exit 2
"""

FAKE_SSH = r"""#!/usr/bin/env bash
echo "ssh $*" >> "$FAKE_LOG"
exit "${FAKE_SSH_RC:-0}"
"""

FAKE_SCP = r"""#!/usr/bin/env bash
# Fake `scp`. A destination without ':' is the remote->local copy of the
# snapshot; FAKE_SCP_MODE=ok drops the fixture there, anything else fails.
echo "scp $*" >> "$FAKE_LOG"
dest="${@: -1}"
if [[ "$dest" == *:* ]]; then exit 0; fi
if [[ "${FAKE_SCP_MODE:-fail}" == ok ]]; then cp "$FAKE_HOP_JSON" "$dest"; exit 0; fi
echo "Connection reset by peer" >&2
exit 1
"""

FAKE_CURL = r"""#!/usr/bin/env bash
echo "curl $*" >> "$FAKE_LOG"
echo 198.51.100.9
"""

FAKE_SLEEP = "#!/usr/bin/env bash\nexit 0\n"


class _FeedHandler(http.server.SimpleHTTPRequestHandler):
    def log_message(self, *args):  # quiet
        pass


def _free_port():
    with socket.socket() as sock:
        sock.bind(("127.0.0.1", 0))
        return sock.getsockname()[1]


class FetchScriptTests(unittest.TestCase):
    """Run fetch_via_walker.sh for real against fake aws/ssh/scp/curl and a local feed.

    These are the failure modes that froze Pages between 2026-09-18T11:25Z
    and 2026-09-19T04:41Z, each of which must now end with a snapshot on
    disk and exit 0, plus the one case that legitimately exits 1.
    """

    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp(prefix="rho-fetch-")
        cls.bin = os.path.join(cls.tmp, "bin")
        os.makedirs(cls.bin)
        for name, body in (("aws", FAKE_AWS), ("ssh", FAKE_SSH), ("scp", FAKE_SCP),
                           ("curl", FAKE_CURL), ("sleep", FAKE_SLEEP)):
            path = os.path.join(cls.bin, name)
            with open(path, "w", encoding="utf-8") as fh:
                fh.write(body)
            os.chmod(path, os.stat(path).st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
        cls.hop_json = os.path.join(HERE, "testdata", "status.json")
        cls.key = os.path.join(cls.tmp, "walker.pem")
        with open(cls.key, "w", encoding="utf-8") as fh:
            fh.write("not a real key\n")
        cls.empty_key = os.path.join(cls.tmp, "empty.pem")
        open(cls.empty_key, "w", encoding="utf-8").close()
        # A local HTTP feed standing in for the status bucket.
        cls.feed_dir = os.path.join(cls.tmp, "feed")
        os.makedirs(cls.feed_dir)
        with open(os.path.join(cls.feed_dir, "status.json"), "w", encoding="utf-8") as fh:
            json.dump(ingest_campaign_feed(), fh)
        port = _free_port()
        handler = lambda *a, **k: _FeedHandler(*a, directory=cls.feed_dir, **k)  # noqa: E731
        cls.server = http.server.ThreadingHTTPServer(("127.0.0.1", port), handler)
        threading.Thread(target=cls.server.serve_forever, daemon=True).start()
        cls.feed_url = "http://127.0.0.1:%d/status.json" % port
        cls.dead_url = "http://127.0.0.1:%d/status.json" % _free_port()

    @classmethod
    def tearDownClass(cls):
        cls.server.shutdown()
        cls.server.server_close()
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def run_script(self, aws_mode, feed_url=None, ssh_rc=0, scp_mode="fail", key=None,
                   skip_walker=False):
        work = tempfile.mkdtemp(dir=self.tmp)
        out = os.path.join(work, "status.json")
        log = os.path.join(work, "calls.log")
        open(log, "w", encoding="utf-8").close()
        env = dict(os.environ)
        env.update({
            "PATH": self.bin + os.pathsep + env.get("PATH", ""),
            "FAKE_LOG": log,
            "FAKE_AWS_MODE": aws_mode,
            "FAKE_SSH_RC": str(ssh_rc),
            "FAKE_SCP_MODE": scp_mode,
            "FAKE_HOP_JSON": self.hop_json,
            "RHO_STATUS_OUT": out,
            "RHO_WORK_FEED_URL": feed_url or self.feed_url,
            "RHO_WORK_FEED_ATTEMPTS": "1",
            "RHO_WORK_FEED_BACKOFF_S": "0",
            "RHO_WALKER_SSH_KEY": self.key if key is None else key,
            "RHO_SKIP_WALKER": "1" if skip_walker else "0",
        })
        env.pop("RHO_WALKER_HOST", None)
        proc = subprocess.run(
            ["bash", os.path.join(HERE, "fetch_via_walker.sh")],
            cwd=work, env=env, capture_output=True, text=True, timeout=120, check=False,
        )
        with open(log, encoding="utf-8") as fh:
            calls = fh.read()
        snapshot = None
        if os.path.exists(out):
            with open(out, encoding="utf-8") as fh:
                snapshot = json.load(fh)
        leftovers = [n for n in os.listdir(work) if n.startswith("status.json.")]
        return proc, snapshot, calls, leftovers

    def test_no_aws_credentials_still_publishes_the_feed(self):
        proc, snapshot, calls, leftovers = self.run_script("fail")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual(snapshot["source"], "ingest-status-feed")
        self.assertEqual(snapshot["dps"], 187000000)
        self.assertIn("walker hop did not produce a snapshot", proc.stderr)
        self.assertIn("from the ingest status feed", proc.stdout)
        self.assertEqual(leftovers, [])

    def test_no_running_walker_still_publishes_the_feed(self):
        proc, snapshot, calls, leftovers = self.run_script("nowalker")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual(snapshot["source"], "ingest-status-feed")
        self.assertIn("no running instance tagged Name=rho-ecc2k-walker", proc.stderr)
        self.assertIn("tagged walkers (any state)", proc.stderr)
        self.assertNotIn("authorize-security-group-ingress", calls)
        self.assertEqual(leftovers, [])

    def test_a_reset_ssh_session_publishes_the_feed_and_revokes_the_rule(self):
        proc, snapshot, calls, leftovers = self.run_script("walker", ssh_rc=255, scp_mode="fail")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual(snapshot["source"], "ingest-status-feed")
        self.assertIn("authorize-security-group-ingress", calls)
        self.assertIn("revoke-security-group-ingress", calls)
        self.assertLess(calls.index("authorize-security-group-ingress"),
                        calls.index("revoke-security-group-ingress"))
        self.assertIn("walker hop failed (ssh=255", proc.stderr)
        self.assertEqual(leftovers, [])

    def test_a_timed_out_query_that_left_a_file_publishes_that_file(self):
        proc, snapshot, calls, leftovers = self.run_script("walker", ssh_rc=124, scp_mode="ok")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual(snapshot["source"], "fixture")
        self.assertIn("publishing the index-only fallback it wrote first", proc.stderr)
        self.assertIn("revoke-security-group-ingress", calls)

    def test_a_working_hop_beats_the_feed(self):
        proc, snapshot, calls, leftovers = self.run_script("walker", ssh_rc=0, scp_mode="ok")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual(snapshot["source"], "fixture")
        self.assertEqual(snapshot["per_worker"][0]["worker_id"], "watch-ip-172-31-58-235")
        self.assertIn("from the walker hop", proc.stdout)
        self.assertIn("revoke-security-group-ingress", calls)
        self.assertEqual(leftovers, [])

    def test_a_working_hop_publishes_even_when_the_feed_is_down(self):
        proc, snapshot, calls, leftovers = self.run_script(
            "walker", feed_url=self.dead_url, ssh_rc=0, scp_mode="ok")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual(snapshot["source"], "fixture")
        self.assertIn("ingest status feed did not answer", proc.stderr)

    def test_both_sources_down_is_the_only_failure(self):
        proc, snapshot, calls, leftovers = self.run_script("nowalker", feed_url=self.dead_url)
        self.assertEqual(proc.returncode, 1)
        self.assertIsNone(snapshot)
        self.assertIn("neither the ingest status feed nor the walker hop", proc.stderr)
        self.assertEqual(leftovers, [])

    def test_skip_walker_publishes_the_feed_without_touching_aws(self):
        proc, snapshot, calls, leftovers = self.run_script("walker", skip_walker=True)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual(snapshot["source"], "ingest-status-feed")
        self.assertEqual(calls, "")

    def test_an_empty_deploy_key_skips_the_hop_rather_than_failing_it(self):
        proc, snapshot, calls, leftovers = self.run_script("walker", key=self.empty_key)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual(snapshot["source"], "ingest-status-feed")
        self.assertIn("RHO_WALKER_SSH_KEY is unset or empty", proc.stderr)
        self.assertNotIn("ssh ", calls)

    def test_feed_fetch_retries_before_giving_up(self):
        import work_feed as wf

        calls = []

        class Boom(OSError):
            pass

        original = wf.urllib.request.urlopen

        def flaky(url, timeout=0):
            calls.append(url)
            if len(calls) < 3:
                raise Boom("blip")
            import io

            class Resp(io.BytesIO):
                def __enter__(self):
                    return self

                def __exit__(self, *exc):
                    return False
            return Resp(json.dumps(work_feed()).encode("utf-8"))

        wf.urllib.request.urlopen = flaky
        try:
            feed = wf.fetch("https://example.invalid/status.json", attempts=3, backoff_s=0)
            self.assertEqual(feed["campaign_id"], "ecc2k-130")
            self.assertEqual(len(calls), 3)
            calls.clear()
            with self.assertRaises(OSError):
                wf.fetch("https://example.invalid/status.json", attempts=2, backoff_s=0)
            self.assertEqual(len(calls), 2)
        finally:
            wf.urllib.request.urlopen = original


class FeedAgeCheck(unittest.TestCase):
    """The end-to-end assertion that a green publish means a moving number."""

    def at(self, minutes_ago, now=None):
        now = now or datetime(2026, 9, 20, 8, 0, tzinfo=timezone.utc)
        return (now - timedelta(minutes=minutes_ago)).strftime("%Y-%m-%dT%H:%M:%SZ")

    @property
    def now(self):
        return datetime(2026, 9, 20, 8, 0, tzinfo=timezone.utc)

    def test_a_current_feed_passes(self):
        ok, message = describe_feed_age(
            {"dps": 10, "generated_at": self.at(4), "published_at": self.at(1)}, self.now)
        self.assertTrue(ok, message)
        self.assertIn("current", message)

    def test_the_20260920_outage_fails_and_blames_the_ingest_side(self):
        # The exact shape of that morning: the ingest host stopped writing at
        # 04:53Z, the publisher kept running every three minutes, and every
        # run was green. Both halves are needed -- a stale generated_at with a
        # fresh published_at is what "the workflow is fine, the feed is dead"
        # looks like, and the message has to say so or the next person reads
        # 802 green runs and goes looking in the wrong repository.
        ok, message = describe_feed_age(
            {"dps": 190497264, "generated_at": self.at(187), "published_at": self.at(13)},
            self.now)
        self.assertFalse(ok)
        self.assertIn("stale", message)
        self.assertIn("ingest feed is what stopped", message)

    def test_a_stopped_publisher_is_named_instead(self):
        ok, message = describe_feed_age(
            {"dps": 190497264, "generated_at": self.at(500), "published_at": self.at(480)},
            self.now)
        self.assertFalse(ok)
        self.assertIn("this workflow is not running either", message)

    def test_a_campaign_with_no_points_is_not_an_outage(self):
        # The committed placeholder. A fresh checkout, or a campaign before its
        # first point, must not read as a dead feed.
        with open(os.path.join(ROOT, "docs", "ecc2k130-status", "status.json")) as fh:
            placeholder = json.load(fh)
        ok, message = describe_feed_age(placeholder, self.now)
        self.assertTrue(ok, message)
        self.assertIn("nothing has been counted", message)

    def test_points_without_a_timestamp_is_still_an_outage(self):
        # The placeholder exemption is guarded on dps: a campaign that
        # collected 190 M points and then lost its timestamp is broken, and
        # must not be excused by the same branch that forgives an empty one.
        ok, message = describe_feed_age({"dps": 190497264, "generated_at": None}, self.now)
        self.assertFalse(ok)
        self.assertIn("no generated_at", message)

    def test_a_future_timestamp_is_not_reported_as_fresh(self):
        ok, message = describe_feed_age(
            {"dps": 10, "generated_at": self.at(-120)}, self.now)
        self.assertFalse(ok)
        self.assertIn("future", message)

    def test_the_threshold_boundary(self):
        limit = MAX_FEED_AGE_S / 60
        ok, _ = describe_feed_age({"dps": 10, "generated_at": self.at(limit - 1)}, self.now)
        self.assertTrue(ok)
        ok, _ = describe_feed_age({"dps": 10, "generated_at": self.at(limit + 1)}, self.now)
        self.assertFalse(ok)

    def test_the_run_reds_before_the_page_apologises(self):
        # Two thresholds describe one event. The run must fail no later than
        # the dashboard starts telling readers the counts are not current --
        # otherwise the page is apologising while CI still reads green, which
        # is the state this whole check exists to remove.
        page = os.path.join(ROOT, "docs", "ecc2k130-status", "index.html")
        with open(page, encoding="utf-8") as fh:
            match = re.search(r"var STALE_AFTER_MS = ([^;]+);", fh.read())
        self.assertIsNotNone(match, "STALE_AFTER_MS not found on the dashboard")
        expression = match.group(1).replace("HOUR_MS", str(3600 * 1000)).strip()
        self.assertRegex(expression, r"^[\d\s*]+$", "unexpected STALE_AFTER_MS: %s" % expression)
        stale_s = eval(expression) / 1000  # noqa: S307 - digits and * only
        self.assertLessEqual(
            MAX_FEED_AGE_S, stale_s,
            "MAX_FEED_AGE_S of %d s is above the dashboard's %d s: the page would "
            "call the counts stale while the run still passed" % (MAX_FEED_AGE_S, stale_s))

    def test_the_threshold_does_not_flap_on_one_late_run(self):
        # Same jitter argument as the banner: the feed is rewritten every
        # --status-every seconds and published once per cron interval, so the
        # limit has to clear one of each with room to spare.
        workflow = os.path.join(ROOT, ".github", "workflows", "ecc2k130-status.yml")
        with open(workflow, encoding="utf-8") as fh:
            crons = re.findall(r'- cron: "([^"]+)"', fh.read())
        minute = crons[0].split(" ", 1)[0]
        step_s = (int(minute[2:]) if minute.startswith("*/") else 60) * 60
        ingest = os.path.join(ROOT, "ecc2k130", "aws", "dp_ingest.py")
        with open(ingest, encoding="utf-8") as fh:
            every = re.search(r'RHO_STATUS_EVERY", "([\d.]+)"', fh.read())
        self.assertIsNotNone(every, "--status-every default not found in dp_ingest.py")
        self.assertGreaterEqual(MAX_FEED_AGE_S, 2 * (float(every.group(1)) + step_s))

    def test_cli_exit_codes_and_warn_only(self):
        import check_feed_age

        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "status.json")
            stale = {"dps": 10, "generated_at": self.at(500, datetime.now(timezone.utc))}
            with open(path, "w") as fh:
                json.dump(stale, fh)
            self.assertEqual(check_feed_age.main(["--status", path]), 1)
            # The rollout lever: report the failure, keep the run green.
            self.assertEqual(check_feed_age.main(["--status", path, "--warn-only"]), 0)
            # A raised limit forgives it, so the threshold is really a knob.
            self.assertEqual(
                check_feed_age.main(["--status", path, "--max-age-s", "99999"]), 0)
            # An unreadable snapshot is a failure, not a pass by default.
            self.assertEqual(
                check_feed_age.main(["--status", os.path.join(tmp, "absent.json")]), 1)


if __name__ == "__main__":
    unittest.main()
