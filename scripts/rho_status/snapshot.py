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
import shutil
import subprocess
import sys
from datetime import datetime, timezone

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

# One scan of distinguished_points, and only on a run that saw ready=false.
#
# This is the same GROUP BY (worker_id, hour) that published in 524s at
# 187 M rows. The first merged run grouped by (worker_id, found_at) — one
# row per object — under a SHARE lock, and PostgreSQL cancelled it at 800s
# (`psql:<stdin>:33`). SHARE waits for ingest writers; with six ingest
# threads that wait is unbounded. The finer grain also spilled more groups
# than the query that actually finished. Later snapshots never reach this
# statement. The Python caller must not run it when ready.
#
# No SHARE lock: the trigger is already live, so ON CONFLICT adds any
# buckets ingest wrote after this transaction's snapshot. A few minutes of
# overlap is tens of thousands of points on 187 M, not a second 800s wait.
BACKFILL_SQL = r"""
BEGIN ISOLATION LEVEL REPEATABLE READ;
SET LOCAL statement_timeout = '1200s';
SET LOCAL lock_timeout = '15s';
INSERT INTO rho_dp_meta (campaign_id, ready)
VALUES ('{{campaign}}', false)
ON CONFLICT (campaign_id) DO NOTHING;
DELETE FROM rho_dp_hour WHERE campaign_id = '{{campaign}}';
DELETE FROM rho_dp_recent WHERE campaign_id = '{{campaign}}';
WITH grouped AS (
  SELECT
    campaign_id,
    worker_id,
    date_trunc('hour', found_at) AS hour,
    count(*)::bigint AS dps,
    min(found_at) AS first_at,
    max(found_at) AS last_at
  FROM distinguished_points
  WHERE campaign_id = '{{campaign}}'
  GROUP BY 1, 2, 3
), ins_hour AS (
  INSERT INTO rho_dp_hour (campaign_id, hour, worker_id, dps, first_at, last_at)
  SELECT campaign_id, hour, worker_id, dps, first_at, last_at
  FROM grouped
  ON CONFLICT (campaign_id, hour, worker_id) DO UPDATE SET
    dps = rho_dp_hour.dps + EXCLUDED.dps,
    first_at = LEAST(rho_dp_hour.first_at, EXCLUDED.first_at),
    last_at = GREATEST(rho_dp_hour.last_at, EXCLUDED.last_at)
), ins_recent AS (
  INSERT INTO rho_dp_recent (campaign_id, worker_id, found_at, dps)
  SELECT campaign_id, worker_id, hour, dps
  FROM grouped
  WHERE hour > now() - interval '7 days'
  ON CONFLICT (campaign_id, worker_id, found_at) DO UPDATE SET
    dps = rho_dp_recent.dps + EXCLUDED.dps
)
UPDATE rho_dp_meta
SET ready = true, backfilled_at = now()
WHERE campaign_id = '{{campaign}}';
COMMIT;
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

# Kept only so a grep for the old one-pass scan still finds the backfill, which
# is the one remaining full-table read, and so tests can pin that the published
# document is assembled from the rollup.
SNAPSHOT_SQL = BACKFILL_SQL


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


def render_sql(template, campaign):
    return template.replace("{{campaign}}", sql_campaign(campaign))


def psql_script(database_url, sql, tuples_only=True):
    cmd = ["psql", database_url, "-v", "ON_ERROR_STOP=1", "-q"]
    if tuples_only:
        cmd += ["-At"]
    proc = subprocess.run(
        cmd + ["-f", "-"],
        input=sql,
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


def psql_json(database_url, campaign):
    sys.stderr.write("rho_status: ensuring rho_dp rollup tables and insert trigger\n")
    psql_script(database_url, ENSURE_SQL, tuples_only=False)
    ready_text = psql_script(database_url, render_sql(READY_SQL, campaign), tuples_only=True)
    if rollup_ready_from_text(ready_text):
        sys.stderr.write("rho_status: rollup ready, skipping distinguished_points scan\n")
    else:
        sys.stderr.write(
            "rho_status: backfilling rho_dp rollup from distinguished_points "
            "(worker, hour — same grain as the 524s snapshot)\n"
        )
        last_err = None
        for attempt in range(1, 4):
            try:
                psql_script(
                    database_url, render_sql(BACKFILL_SQL, campaign), tuples_only=False
                )
                last_err = None
                break
            except RuntimeError as err:
                last_err = err
                text = str(err).lower()
                if attempt < 3 and (
                    "serialize" in text
                    or "deadlock" in text
                    or "lock timeout" in text
                    or "lock_not_available" in text
                ):
                    sys.stderr.write(
                        "rho_status: backfill attempt %d lost a lock race; retrying\n"
                        % attempt
                    )
                    continue
                raise
        if last_err is not None:
            raise last_err
    psql_script(database_url, PRUNE_SQL, tuples_only=False)
    text = psql_script(database_url, render_sql(READ_SQL, campaign), tuples_only=True)
    lines = [line for line in text.strip().splitlines() if line.startswith("{")]
    if not lines:
        raise RuntimeError(
            "campaign %s rollup not ready (backfill did not complete)" % campaign
        )
    return json.loads(lines[-1])


def take_snapshot(database_url, campaign, source="direct"):
    if not database_url:
        raise SystemExit("DATABASE_URL is required")
    # The walker hop has psql and often no psycopg. Prefer psql so the
    # multi-statement ensure/backfill script is one session with COMMIT
    # between the trigger becoming visible and the REPEATABLE READ backfill.
    if not shutil.which("psql"):
        raise SystemExit("psql is required to take a snapshot")
    row = psql_json(database_url, campaign)
    if isinstance(row, str):
        row = json.loads(row)
    snapshot = normalize(row, campaign, source)
    snapshot["query_driver"] = "psql"
    return snapshot


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--database-url", default=os.environ.get("DATABASE_URL"))
    parser.add_argument("--campaign", default=os.environ.get("RHO_CAMPAIGN", DEFAULT_CAMPAIGN))
    parser.add_argument("--out", default="-")
    parser.add_argument("--source", default="direct")
    args = parser.parse_args(argv)
    snapshot = take_snapshot(args.database_url, args.campaign, args.source)
    text = json.dumps(snapshot, indent=2, sort_keys=True) + "\n"
    if args.out == "-":
        sys.stdout.write(text)
    else:
        os.makedirs(os.path.dirname(os.path.abspath(args.out)) or ".", exist_ok=True)
        with open(args.out, "w", encoding="utf-8") as fh:
            fh.write(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
