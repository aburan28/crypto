#!/usr/bin/env python3
"""Read-only ECC2K-130 distinguished-point campaign snapshot.

Queries the shared Pollard-rho DP store and writes a public JSON document.
The document is aggregates only: no point keys, walk coefficients, seeds,
or connection strings.

Usage:
    DATABASE_URL=postgresql://... python3 scripts/rho_status/snapshot.py
    python3 scripts/rho_status/snapshot.py --database-url URL --out status.json
"""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import time
from datetime import datetime, timedelta, timezone

DEFAULT_CAMPAIGN = "ecc2k-130"
SCHEMA_VERSION = 1
CLAIM_BOUNDARY = (
    "Public research campaign aggregates for Certicom ECC2K-130 distinguished-point "
    "collection. Not a discrete-log recovery, not a verified collision, and not a "
    "statement about deployed-curve security unless collisions>0 and an independent "
    "solver has checked [k]P = Q."
)

# Lifetime aggregates used to come from one GROUP BY over distinguished_points.
# That was affordable at a few million rows and is not at ~187 M: the walker
# hop took 524 s on 2026-09-18T11:25Z and the next run died at the 540 s remote
# timeout, so Pages froze on that snapshot and the dashboard's stale banner is
# what a visitor sees. An index on (campaign_id, found_at) helps the ingest's
# *count* query (#432) and does not help this one, because grouping by
# worker_id still heap-fetches every row.
#
# The durable answer is the one #432 named and deferred: a rollup maintained
# from INSERTs. Another btree on distinguished_points does not get Pages
# under 540s: the ingest already builds (campaign_id, found_at), and this
# query still has to read worker_id for every row. A covering index would
# still be O(corpus) every 15 minutes. snapshot.py creates the rollup,
# backfills it once, and installs a statement-level trigger so the ingest's
# INSERT ... SELECT keeps the buckets current without a code deploy on that
# host. Later snapshots read the rollup.
#
# #448 landed that design and the live page still did not move. The merged
# backfill set statement_timeout to 800s, took SHARE on distinguished_points,
# and grouped by (worker_id, found_at) — one row per object. PostgreSQL
# cancelled it at 800s, rolled back, and left rho_dp_meta.ready false, so
# every later publish retried the same scan. The hop copies this file each
# run, so the fix is here: (1) write an index-only fallback *before* the
# backfill so a later SIGTERM still leaves status.json to scp, (2) backfill
# one hour at a time through (campaign_id, found_at), resume from
# rho_dp_meta.backfill_through, never lock the heap.

ENSURE_SQL = r"""
CREATE TABLE IF NOT EXISTS rho_dp_hour (
  campaign_id text NOT NULL,
  hour        timestamptz NOT NULL,
  worker_id   text NOT NULL,
  dps         bigint NOT NULL,
  first_at    timestamptz NOT NULL,
  last_at     timestamptz NOT NULL,
  PRIMARY KEY (campaign_id, hour, worker_id)
);
CREATE TABLE IF NOT EXISTS rho_dp_recent (
  campaign_id text NOT NULL,
  worker_id   text NOT NULL,
  found_at    timestamptz NOT NULL,
  dps         bigint NOT NULL,
  PRIMARY KEY (campaign_id, worker_id, found_at)
);
CREATE TABLE IF NOT EXISTS rho_dp_meta (
  campaign_id   text PRIMARY KEY,
  ready         boolean NOT NULL DEFAULT false,
  backfilled_at timestamptz
);
ALTER TABLE rho_dp_meta ADD COLUMN IF NOT EXISTS backfill_through timestamptz;
CREATE INDEX IF NOT EXISTS rho_dp_hour_hour
  ON rho_dp_hour (campaign_id, hour);
CREATE INDEX IF NOT EXISTS rho_dp_recent_found_at
  ON rho_dp_recent (campaign_id, found_at);

CREATE OR REPLACE FUNCTION rho_dp_rollup_insert() RETURNS trigger
LANGUAGE plpgsql
SECURITY DEFINER
SET search_path = public
AS $fn$
BEGIN
  INSERT INTO rho_dp_hour (campaign_id, hour, worker_id, dps, first_at, last_at)
  SELECT campaign_id,
         date_trunc('hour', found_at),
         worker_id,
         count(*)::bigint,
         min(found_at),
         max(found_at)
  FROM new_rows
  GROUP BY 1, 2, 3
  ON CONFLICT (campaign_id, hour, worker_id) DO UPDATE SET
    dps = rho_dp_hour.dps + EXCLUDED.dps,
    first_at = LEAST(rho_dp_hour.first_at, EXCLUDED.first_at),
    last_at = GREATEST(rho_dp_hour.last_at, EXCLUDED.last_at);

  INSERT INTO rho_dp_recent (campaign_id, worker_id, found_at, dps)
  SELECT campaign_id, worker_id, found_at, count(*)::bigint
  FROM new_rows
  GROUP BY 1, 2, 3
  ON CONFLICT (campaign_id, worker_id, found_at) DO UPDATE SET
    dps = rho_dp_recent.dps + EXCLUDED.dps;
  RETURN NULL;
END;
$fn$;

DO $tg$
BEGIN
  IF NOT EXISTS (
    SELECT 1
    FROM pg_trigger t
    JOIN pg_class c ON c.oid = t.tgrelid
    WHERE t.tgname = 'rho_dp_rollup_insert'
      AND c.relname = 'distinguished_points'
      AND NOT t.tgisinternal
  ) THEN
    EXECUTE $create$
      CREATE TRIGGER rho_dp_rollup_insert
      AFTER INSERT ON distinguished_points
      REFERENCING NEW TABLE AS new_rows
      FOR EACH STATEMENT
      EXECUTE PROCEDURE rho_dp_rollup_insert()
    $create$;
  END IF;
END
$tg$;
"""

READY_SQL = r"""
SELECT COALESCE(
  (SELECT ready FROM rho_dp_meta WHERE campaign_id = '{{campaign}}'),
  false
);
"""

INIT_META_SQL = r"""
INSERT INTO rho_dp_meta (campaign_id, ready)
VALUES ('{{campaign}}', false)
ON CONFLICT (campaign_id) DO NOTHING;
"""

BOUNDS_SQL = r"""
SET statement_timeout = '30s';
SELECT
  (SELECT to_char(
     date_trunc('hour', min(found_at)) AT TIME ZONE 'UTC',
     'YYYY-MM-DD"T"HH24:MI:SS"Z"'
   )
   FROM distinguished_points
   WHERE campaign_id = '{{campaign}}'),
  to_char(
    date_trunc('hour', now()) AT TIME ZONE 'UTC',
    'YYYY-MM-DD"T"HH24:MI:SS"Z"'
  ),
  (SELECT to_char(
     backfill_through AT TIME ZONE 'UTC',
     'YYYY-MM-DD"T"HH24:MI:SS"Z"'
   )
   FROM rho_dp_meta
   WHERE campaign_id = '{{campaign}}');
"""

# One hour of distinguished_points, using the #432 (campaign_id, found_at)
# index. Peak hours on this campaign are ~10 M rows, not 187 M. No SHARE
# lock: the trigger is already live, so ON CONFLICT adds buckets ingest
# wrote after this transaction's snapshot. DELETE + scan of one closed
# hour is how a mid-hour trigger start gets the pre-trigger rows without
# wiping hours already finished.
#
# REPEATABLE READ: the scan must not see rows the trigger already counted
# after DELETE. READ COMMITTED + ON CONFLICT ADD double-counts ingest that
# commits in that window, and backfill_through would freeze the inflation.
# rho_dp_recent keeps the object's found_at, not the hour bucket: READ_SQL
# slides a 1h/24h window over that column, and an hour timestamp would make
# dps_last_hour drop to the current clock hour (IDLE_OR_STALE after :00).
BACKFILL_SQL = r"""
BEGIN ISOLATION LEVEL REPEATABLE READ;
SET LOCAL statement_timeout = '300s';
SET LOCAL lock_timeout = '15s';
DELETE FROM rho_dp_hour
 WHERE campaign_id = '{{campaign}}'
   AND hour = '{{hour}}'::timestamptz;
DELETE FROM rho_dp_recent
 WHERE campaign_id = '{{campaign}}'
   AND found_at >= '{{hour}}'::timestamptz
   AND found_at < '{{hour}}'::timestamptz + interval '1 hour';
WITH src AS MATERIALIZED (
  SELECT campaign_id, worker_id, found_at
  FROM distinguished_points
  WHERE campaign_id = '{{campaign}}'
    AND found_at >= '{{hour}}'::timestamptz
    AND found_at < '{{hour}}'::timestamptz + interval '1 hour'
), ins_hour AS (
  INSERT INTO rho_dp_hour (campaign_id, hour, worker_id, dps, first_at, last_at)
  SELECT
    campaign_id,
    date_trunc('hour', found_at),
    worker_id,
    count(*)::bigint,
    min(found_at),
    max(found_at)
  FROM src
  GROUP BY 1, 2, 3
  ON CONFLICT (campaign_id, hour, worker_id) DO UPDATE SET
    dps = rho_dp_hour.dps + EXCLUDED.dps,
    first_at = LEAST(rho_dp_hour.first_at, EXCLUDED.first_at),
    last_at = GREATEST(rho_dp_hour.last_at, EXCLUDED.last_at)
), ins_recent AS (
  INSERT INTO rho_dp_recent (campaign_id, worker_id, found_at, dps)
  SELECT campaign_id, worker_id, found_at, count(*)::bigint
  FROM src
  WHERE found_at > now() - interval '7 days'
  GROUP BY 1, 2, 3
  ON CONFLICT (campaign_id, worker_id, found_at) DO UPDATE SET
    dps = rho_dp_recent.dps + EXCLUDED.dps
)
UPDATE rho_dp_meta
SET backfill_through = '{{hour}}'::timestamptz
WHERE campaign_id = '{{campaign}}';
COMMIT;
"""

MARK_READY_SQL = r"""
UPDATE rho_dp_meta
SET ready = true, backfilled_at = now()
WHERE campaign_id = '{{campaign}}'
  AND backfill_through >= date_trunc('hour', now());
"""

# Index-only: group by hour, never by worker_id, so this does not heap-fetch
# 187 M rows. Used to publish a fresh generated_at while the worker rollup
# is still catching up. per_worker is empty until the rollup is ready.
FALLBACK_SQL = r"""
SET statement_timeout = '600s';
SELECT json_build_object(
  'campaign_id', c.campaign_id,
  'curve_id', c.curve_id,
  'dp_mask_bits', c.dp_mask_bits,
  'campaign_created_at', c.created_at,
  'meta', c.meta,
  'dps', COALESCE(s.dps, 0),
  'workers', 0,
  'first_dp_at', s.first_dp,
  'last_dp_at', s.last_dp,
  'dps_last_hour', COALESCE(s.dps_last_hour, 0),
  'dps_last_day', COALESCE(s.dps_last_day, 0),
  'collisions', COALESCE(k.collisions, 0),
  'latest_collision_at', k.latest_collision_at,
  'per_worker', '[]'::json,
  'hourly', COALESCE(s.hourly, '[]'::json),
  'rollup_ready', false
)
FROM rho_campaigns c
LEFT JOIN LATERAL (
  SELECT
    sum(bucket.dps)::bigint AS dps,
    min(bucket.first_at) AS first_dp,
    max(bucket.last_at) AS last_dp,
    COALESCE(sum(bucket.dps_last_hour), 0)::bigint AS dps_last_hour,
    COALESCE(sum(bucket.dps_last_day), 0)::bigint AS dps_last_day,
    json_agg(json_build_object('hour', bucket.hour, 'dps', bucket.dps) ORDER BY bucket.hour)
      FILTER (WHERE bucket.hour > now() - interval '7 days' AND bucket.dps > 0) AS hourly
  FROM (
    SELECT
      date_trunc('hour', found_at) AS hour,
      count(*)::bigint AS dps,
      min(found_at) AS first_at,
      max(found_at) AS last_at,
      count(*) FILTER (WHERE found_at > now() - interval '1 hour')::bigint AS dps_last_hour,
      count(*) FILTER (WHERE found_at > now() - interval '24 hours')::bigint AS dps_last_day
    FROM distinguished_points
    WHERE campaign_id = c.campaign_id
    GROUP BY 1
  ) bucket
) s ON true
LEFT JOIN LATERAL (
  SELECT
    count(*)::bigint AS collisions,
    max(detected_at) AS latest_collision_at
  FROM rho_collisions x
  WHERE x.campaign_id = c.campaign_id
) k ON true
WHERE c.campaign_id = '{{campaign}}';
"""

PRUNE_SQL = r"""
DELETE FROM rho_dp_recent
WHERE found_at < now() - interval '8 days';
"""

READ_SQL = r"""
SELECT json_build_object(
  'campaign_id', c.campaign_id,
  'curve_id', c.curve_id,
  'dp_mask_bits', c.dp_mask_bits,
  'campaign_created_at', c.created_at,
  'meta', c.meta,
  'dps', COALESCE(s.dps, 0),
  'workers', COALESCE(s.workers, 0),
  'first_dp_at', s.first_dp,
  'last_dp_at', s.last_dp,
  'dps_last_hour', COALESCE(r.dps_last_hour, 0),
  'dps_last_day', COALESCE(r.dps_last_day, 0),
  'collisions', COALESCE(k.collisions, 0),
  'latest_collision_at', k.latest_collision_at,
  'per_worker', COALESCE(w.per_worker, '[]'::json),
  'hourly', COALESCE(h.hourly, '[]'::json),
  'rollup_ready', true
)
FROM rho_campaigns c
JOIN rho_dp_meta m ON m.campaign_id = c.campaign_id AND m.ready
LEFT JOIN LATERAL (
  SELECT
    sum(dps)::bigint AS dps,
    count(DISTINCT worker_id)::bigint AS workers,
    min(first_at) AS first_dp,
    max(last_at) AS last_dp
  FROM rho_dp_hour
  WHERE campaign_id = c.campaign_id
) s ON true
LEFT JOIN LATERAL (
  SELECT
    COALESCE(sum(dps) FILTER (WHERE found_at > now() - interval '1 hour'), 0)::bigint AS dps_last_hour,
    COALESCE(sum(dps) FILTER (WHERE found_at > now() - interval '24 hours'), 0)::bigint AS dps_last_day
  FROM rho_dp_recent
  WHERE campaign_id = c.campaign_id
) r ON true
LEFT JOIN LATERAL (
  SELECT
    count(*)::bigint AS collisions,
    max(detected_at) AS latest_collision_at
  FROM rho_collisions x
  WHERE x.campaign_id = c.campaign_id
) k ON true
LEFT JOIN LATERAL (
  SELECT json_agg(row_to_json(p) ORDER BY p.dps DESC) AS per_worker
  FROM (
    SELECT
      worker_id,
      sum(dps)::bigint AS dps,
      min(first_at) AS first_dp_at,
      max(last_at) AS last_dp_at
    FROM rho_dp_hour
    WHERE campaign_id = c.campaign_id
    GROUP BY worker_id
    ORDER BY sum(dps) DESC
    LIMIT 32
  ) p
) w ON true
LEFT JOIN LATERAL (
  SELECT json_agg(row_to_json(y) ORDER BY y.hour) AS hourly
  FROM (
    SELECT hour, sum(dps)::bigint AS dps
    FROM rho_dp_hour
    WHERE campaign_id = c.campaign_id
      AND hour > now() - interval '7 days'
    GROUP BY hour
    HAVING sum(dps) > 0
    ORDER BY hour
  ) y
) h ON true
WHERE c.campaign_id = '{{campaign}}';
"""

# Kept so a grep for the old one-pass scan still finds the backfill, which
# is the remaining distinguished_points read that groups by worker_id.
SNAPSHOT_SQL = BACKFILL_SQL

HOUR_RE = re.compile(r"^\d{4}-\d{2}-\d{2}T\d{2}:00:00Z$")
SESSION_PREAMBLE = "SET statement_timeout = 0;\n"


FORBIDDEN_PUBLIC_KEYS = (
    "point_key",
    "a",
    "b",
    "a1",
    "b1",
    "a2",
    "b2",
    "walk_seed",
    "walk_seed1",
    "walk_seed2",
    "database_url",
    "password",
    "dsn",
)


def utcnow():
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def iso(value):
    if value is None:
        return None
    if hasattr(value, "isoformat"):
        text = value.isoformat()
        if text.endswith("+00:00"):
            return text[:-6] + "Z"
        return text
    return str(value)


def assert_public(blob):
    dumped = json.dumps(blob)
    for key in FORBIDDEN_PUBLIC_KEYS:
        needle = '"%s"' % key
        if needle in dumped:
            raise SystemExit("refusing to publish field %s" % key)


def campaign_state(row):
    collisions = int(row.get("collisions") or 0)
    if collisions > 0:
        return "COLLISION_RECORDED"
    if int(row.get("dps_last_hour") or 0) > 0:
        return "COLLECTING"
    outstanding = int(row.get("ingest_outstanding") or 0)
    unrecognised = int(row.get("ingest_unrecognised") or 0)
    if outstanding > 0 or unrecognised > 0:
        return "INGEST_BEHIND"
    if int(row.get("dps") or 0) > 0:
        return "IDLE_OR_STALE"
    return "EMPTY"


def normalize(row, campaign, source):
    hourly = row.get("hourly") or []
    if isinstance(hourly, str):
        hourly = json.loads(hourly)
    workers = row.get("per_worker") or []
    if isinstance(workers, str):
        workers = json.loads(workers)
    snapshot = {
        "schema_version": SCHEMA_VERSION,
        "generated_at": utcnow(),
        "source": source,
        "campaign_id": row.get("campaign_id") or campaign,
        "curve_id": row.get("curve_id"),
        "dp_mask_bits": row.get("dp_mask_bits"),
        "campaign_created_at": iso(row.get("campaign_created_at")),
        "state": campaign_state(row),
        "dps": int(row.get("dps") or 0),
        "collisions": int(row.get("collisions") or 0),
        "workers": int(row.get("workers") or 0),
        "first_dp_at": iso(row.get("first_dp_at")),
        "last_dp_at": iso(row.get("last_dp_at")),
        "latest_collision_at": iso(row.get("latest_collision_at")),
        "dps_last_hour": int(row.get("dps_last_hour") or 0),
        "dps_last_day": int(row.get("dps_last_day") or 0),
        "per_worker": [
            {
                "worker_id": item.get("worker_id") or "unknown",
                "dps": int(item.get("dps") or 0),
                "first_dp_at": iso(item.get("first_dp_at")),
                "last_dp_at": iso(item.get("last_dp_at")),
            }
            for item in workers
        ],
        "hourly": [
            {
                "hour": iso(item.get("hour")),
                "dps": int(item.get("dps") or 0),
            }
            for item in hourly
        ],
        "claim_boundary": CLAIM_BOUNDARY,
    }
    assert_public(snapshot)
    return snapshot


def sql_campaign(campaign):
    if not campaign or any(ch in campaign for ch in "\n\r;"):
        raise SystemExit("refusing campaign id %r" % campaign)
    return campaign.replace("'", "''")


def sql_hour(hour):
    if isinstance(hour, datetime):
        text = hour.astimezone(timezone.utc).strftime("%Y-%m-%dT%H:00:00Z")
    else:
        text = str(hour)
    if not HOUR_RE.match(text):
        raise SystemExit("refusing hour %r" % text)
    return text


def parse_hour(text):
    text = (text or "").strip()
    if not text:
        return None
    return datetime.strptime(sql_hour(text), "%Y-%m-%dT%H:00:00Z").replace(tzinfo=timezone.utc)


def parse_bounds(text):
    lines = [line for line in (text or "").strip().splitlines() if line]
    if not lines:
        return None, None, None
    parts = lines[-1].split("|")
    while len(parts) < 3:
        parts.append("")
    return parse_hour(parts[0]), parse_hour(parts[1]), parse_hour(parts[2])


def hours_to_backfill(start_hour, now_hour, through):
    if start_hour is None or now_hour is None:
        return []
    hour = start_hour
    if through is not None:
        hour = through + timedelta(hours=1)
    out = []
    while hour <= now_hour:
        out.append(hour)
        hour += timedelta(hours=1)
    return out


def render_sql(template, campaign, hour=None):
    sql = template.replace("{{campaign}}", sql_campaign(campaign))
    if hour is not None:
        sql = sql.replace("{{hour}}", sql_hour(hour))
    return sql


def psql_script(database_url, sql, tuples_only=True):
    cmd = ["psql", database_url, "-v", "ON_ERROR_STOP=1", "-q"]
    if tuples_only:
        cmd += ["-At"]
    proc = subprocess.run(
        cmd + ["-f", "-"],
        input=SESSION_PREAMBLE + sql,
        capture_output=True,
        text=True,
        check=False,
    )
    if proc.stderr:
        sys.stderr.write(proc.stderr)
        if not proc.stderr.endswith("\n"):
            sys.stderr.write("\n")
    if proc.returncode != 0:
        raise RuntimeError(proc.stderr.strip() or proc.stdout.strip() or "psql failed")
    return proc.stdout


def rollup_ready_from_text(text):
    token = (text or "").strip().splitlines()
    if not token:
        return False
    return token[-1].lower() in ("t", "true", "1")


def is_lock_race(err):
    text = str(err).lower()
    return (
        "serialize" in text
        or "deadlock" in text
        or "lock timeout" in text
        or "lock_not_available" in text
    )


def parse_json_row(text, what):
    lines = [line for line in (text or "").strip().splitlines() if line.startswith("{")]
    if not lines:
        raise RuntimeError(what)
    return json.loads(lines[-1])


def snapshot_budget_s():
    raw = os.environ.get("RHO_SNAPSHOT_BUDGET_S", "1400")
    try:
        return max(60.0, float(raw))
    except ValueError:
        return 1400.0


def write_snapshot_file(path, snapshot):
    text = json.dumps(snapshot, indent=2, sort_keys=True) + "\n"
    if path == "-":
        sys.stdout.write(text)
        return
    directory = os.path.dirname(os.path.abspath(path))
    os.makedirs(directory or ".", exist_ok=True)
    tmp = path + ".tmp"
    with open(tmp, "w", encoding="utf-8") as fh:
        fh.write(text)
    os.replace(tmp, path)


def backfill_one_hour(database_url, campaign, hour):
    last_err = None
    for attempt in range(1, 4):
        try:
            psql_script(
                database_url,
                render_sql(BACKFILL_SQL, campaign, hour=hour),
                tuples_only=False,
            )
            return
        except RuntimeError as err:
            last_err = err
            if attempt < 3 and is_lock_race(err):
                sys.stderr.write(
                    "rho_status: backfill %s attempt %d lost a lock race; retrying\n"
                    % (sql_hour(hour), attempt)
                )
                continue
            raise
    if last_err is not None:
        raise last_err


def backfill_hours(database_url, campaign, started, budget):
    psql_script(database_url, render_sql(INIT_META_SQL, campaign), tuples_only=False)
    bounds_text = psql_script(
        database_url, render_sql(BOUNDS_SQL, campaign), tuples_only=True
    )
    start_hour, now_hour, through = parse_bounds(bounds_text)
    hours = hours_to_backfill(start_hour, now_hour, through)
    if not hours:
        sys.stderr.write("rho_status: no distinguished_points hours to backfill\n")
        psql_script(database_url, render_sql(MARK_READY_SQL, campaign), tuples_only=False)
        return True
    sys.stderr.write(
        "rho_status: backfilling %d hour(s) from %s through %s\n"
        % (len(hours), sql_hour(hours[0]), sql_hour(hours[-1]))
    )
    for index, hour in enumerate(hours, 1):
        leftover = budget - (time.monotonic() - started)
        if leftover < 90:
            sys.stderr.write(
                "rho_status: backfill budget reached after %d/%d hours; will resume\n"
                % (index - 1, len(hours))
            )
            return False
        sys.stderr.write(
            "rho_status: backfill hour %s (%d/%d, %.0fs left)\n"
            % (sql_hour(hour), index, len(hours), leftover)
        )
        backfill_one_hour(database_url, campaign, hour)
    psql_script(database_url, render_sql(MARK_READY_SQL, campaign), tuples_only=False)
    return True


def psql_json(database_url, campaign, publish=None):
    started = time.monotonic()
    budget = snapshot_budget_s()
    sys.stderr.write("rho_status: ensuring rho_dp rollup tables and insert trigger\n")
    psql_script(database_url, ENSURE_SQL, tuples_only=False)
    ready_text = psql_script(database_url, render_sql(READY_SQL, campaign), tuples_only=True)
    if rollup_ready_from_text(ready_text):
        sys.stderr.write("rho_status: rollup ready, skipping distinguished_points scan\n")
    else:
        # Publish the index-only document first so a later timeout still
        # leaves a file the hop can scp. Then spend the rest of the budget
        # on hour chunks; the next run resumes.
        sys.stderr.write(
            "rho_status: rollup not ready; publishing index-only fallback first\n"
        )
        fallback = parse_json_row(
            psql_script(
                database_url, render_sql(FALLBACK_SQL, campaign), tuples_only=True
            ),
            "campaign %s fallback snapshot returned no json" % campaign,
        )
        if publish is not None:
            publish(fallback)
        try:
            backfill_hours(database_url, campaign, started, budget)
        except RuntimeError as err:
            sys.stderr.write("rho_status: backfill interrupted: %s\n" % err)
            return fallback
        ready_text = psql_script(
            database_url, render_sql(READY_SQL, campaign), tuples_only=True
        )
        if not rollup_ready_from_text(ready_text):
            sys.stderr.write(
                "rho_status: rollup still not ready; keeping the index-only snapshot\n"
            )
            return fallback
    psql_script(database_url, PRUNE_SQL, tuples_only=False)
    text = psql_script(database_url, render_sql(READ_SQL, campaign), tuples_only=True)
    return parse_json_row(
        text, "campaign %s rollup not ready (backfill did not complete)" % campaign
    )


def take_snapshot(database_url, campaign, source="direct", out=None):
    if not database_url:
        raise SystemExit("DATABASE_URL is required")
    # The walker hop has psql and often no psycopg. Prefer psql so the
    # ensure/backfill statements run as separate sessions with COMMIT
    # between the trigger becoming visible and each hour chunk.
    if not shutil.which("psql"):
        raise SystemExit("psql is required to take a snapshot")

    snapshot = {"query_driver": "psql"}

    def publish(row):
        if isinstance(row, str):
            row = json.loads(row)
        built = normalize(row, campaign, source)
        built["query_driver"] = "psql"
        snapshot.clear()
        snapshot.update(built)
        if out:
            write_snapshot_file(out, built)
        return built

    row = psql_json(database_url, campaign, publish=publish)
    return publish(row)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--database-url", default=os.environ.get("DATABASE_URL"))
    parser.add_argument("--campaign", default=os.environ.get("RHO_CAMPAIGN", DEFAULT_CAMPAIGN))
    parser.add_argument("--out", default="-")
    parser.add_argument("--source", default="direct")
    args = parser.parse_args(argv)
    out = None if args.out == "-" else args.out
    snapshot = take_snapshot(args.database_url, args.campaign, args.source, out=out)
    if args.out == "-":
        write_snapshot_file("-", snapshot)
    return 0


if __name__ == "__main__":
    sys.exit(main())
