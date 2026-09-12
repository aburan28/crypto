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

SNAPSHOT_SQL = r"""
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
  'dps_last_hour', COALESCE(s.dps_last_hour, 0),
  'dps_last_day', COALESCE(s.dps_last_day, 0),
  'collisions', COALESCE(k.collisions, 0),
  'latest_collision_at', k.latest_collision_at,
  'per_worker', COALESCE(w.per_worker, '[]'::json),
  'hourly', COALESCE(h.hourly, '[]'::json)
)
FROM rho_campaigns c
LEFT JOIN LATERAL (
  SELECT
    count(*)::bigint AS dps,
    count(DISTINCT worker_id)::bigint AS workers,
    min(found_at) AS first_dp,
    max(found_at) AS last_dp,
    count(*) FILTER (WHERE found_at > now() - interval '1 hour')::bigint AS dps_last_hour,
    count(*) FILTER (WHERE found_at > now() - interval '24 hours')::bigint AS dps_last_day
  FROM distinguished_points d
  WHERE d.campaign_id = c.campaign_id
) s ON true
LEFT JOIN LATERAL (
  SELECT
    count(*)::bigint AS collisions,
    max(detected_at) AS latest_collision_at
  FROM rho_collisions r
  WHERE r.campaign_id = c.campaign_id
) k ON true
LEFT JOIN LATERAL (
  SELECT json_agg(row_to_json(x) ORDER BY x.dps DESC) AS per_worker
  FROM (
    SELECT
      worker_id,
      count(*)::bigint AS dps,
      min(found_at) AS first_dp_at,
      max(found_at) AS last_dp_at
    FROM distinguished_points d
    WHERE d.campaign_id = c.campaign_id
    GROUP BY worker_id
    ORDER BY count(*) DESC
    LIMIT 32
  ) x
) w ON true
LEFT JOIN LATERAL (
  SELECT json_agg(row_to_json(y) ORDER BY y.hour) AS hourly
  FROM (
    SELECT
      date_trunc('hour', found_at) AS hour,
      count(*)::bigint AS dps
    FROM distinguished_points d
    WHERE d.campaign_id = c.campaign_id
      AND found_at > now() - interval '7 days'
    GROUP BY 1
    ORDER BY 1
  ) y
) h ON true
WHERE c.campaign_id = %(campaign)s;
"""


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


def psql_json(database_url, campaign):
    sql = SNAPSHOT_SQL.replace("%(campaign)s", "'%s'" % campaign.replace("'", "''"))
    proc = subprocess.run(
        ["psql", database_url, "-v", "ON_ERROR_STOP=1", "-At", "-c", sql],
        capture_output=True,
        text=True,
        check=False,
    )
    if proc.returncode != 0:
        raise RuntimeError(proc.stderr.strip() or proc.stdout.strip() or "psql failed")
    line = proc.stdout.strip().splitlines()
    if not line:
        raise RuntimeError("campaign %s not found" % campaign)
    return json.loads(line[0])


def psycopg_json(database_url, campaign):
    import psycopg

    sql = SNAPSHOT_SQL.replace("%(campaign)s", "%s")
    with psycopg.connect(database_url) as conn:
        with conn.cursor() as cur:
            cur.execute("SET default_transaction_read_only = on")
            cur.execute(sql, (campaign,))
            row = cur.fetchone()
    if not row or row[0] is None:
        raise RuntimeError("campaign %s not found" % campaign)
    return row[0]


def take_snapshot(database_url, campaign, source="direct"):
    if not database_url:
        raise SystemExit("DATABASE_URL is required")
    try:
        import psycopg  # noqa: F401

        row = psycopg_json(database_url, campaign)
        driver = "psycopg"
    except ImportError:
        row = psql_json(database_url, campaign)
        driver = "psql"
    if isinstance(row, str):
        row = json.loads(row)
    snapshot = normalize(row, campaign, source)
    snapshot["query_driver"] = driver
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
