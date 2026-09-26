#!/usr/bin/env python3
"""Ingest campaign distinguished points from S3 into the rho-dp Postgres store.

The store is a *derived view* for the public dashboard. `s3://bucket/dp/` is the
authoritative corpus and `merge.py` is what searches it for collisions, so this
program is allowed to lag, be restarted, or be re-run from scratch; what it may
not do is invent, drop, or duplicate a point.

It replaces a hand-started directory watcher on the walker host that was not in
version control, had no service unit, competed with that host's own walks for
CPU, and did not survive a database restart. Three properties follow from that
history and are the reason this is a file rather than a shell one-liner:

  * **Idempotent.** `distinguished_points` is keyed `(campaign_id, point_key)`,
    so every insert is `ON CONFLICT DO NOTHING`. Re-ingesting an object is
    harmless, which is what makes crash recovery trivial.
  * **No state outside the database it describes.** Progress lives in
    `dp_ingest_progress` and is written in the same transaction as the points,
    so a record of having ingested an object cannot outlive the insert it
    describes. Drop the table and the next pass rebuilds it; there is no file
    on some host that can disagree with the store.
  * **Restartable by the machine, not by a person.** Deployed as a systemd unit
    with `Restart=always`; a database failover costs a retry, not an outage that
    ends when somebody notices.
  * **Its published figures cost the same on any day of the campaign.** The
    snapshot used to ask the whole table for every number it published, every
    few minutes, on the instance the ingest writes to. `dp_ingest_totals` and
    `dp_ingest_hourly` are maintained by the insert itself, in the transaction
    that adds the points, so they cannot drift from what they summarise and a
    publish reads one row and 48 more.
  * **It refuses what does not match its own commit marker.** The worker
    writes `<object>.bin.json` after the payload, carrying a sha256 and a
    record count; `merge.py` verifies it over the same corpus and this program
    did not. A truncated read self-heals, because the object stays short of
    its S3 size and is retried, but a body corrupted in place does not, and
    once stored it is indistinguishable from real points to the collision
    check above. Absent for legacy keys, which is why a missing marker is not
    an error.
  * **It reports the collision it is uniquely placed to see.** The unique key
    that makes re-ingest harmless is the same constraint that defines a
    collision: two walks reaching one distinguished point from different
    seeds. `DO NOTHING` cannot tell that from a worker re-reporting its own
    point, and dropped both alike, so the campaign's terminal event was the
    one thing this program threw away -- and `rho_collisions`, which it and
    scripts/rho_status/snapshot.py both read, had no writer at all. A
    conflicting record is now compared against the stored one and a genuine
    meeting is recorded and published at once. `merge.py` over the S3 corpus
    stays the independent check rather than the only detector.

Record and column encodings are not invented here. They are read back from
`rho_campaigns.meta`, which the original ingest wrote:

    format         ecc2k130-gpu-packed32
    point_key      3xuint64-le orbit representative (x only)
    walk_seed      GPU client starting seed (big-endian integer)
    coeff_encoding 17-byte big-endian mod 2^131; a=start-seed, b=0

so a 32-byte record `(seed, k0, k1, k2)` little-endian (`merge.py`'s `RECORD`)
maps to `point_key = k0||k1||k2` as stored — literally bytes 8..32 of the record
— with `a` the seed in 17 big-endian bytes, `b` seventeen zero bytes, and
`walk_seed` the seed big-endian. `--verify` checks that against points the old
ingester already wrote before this program is trusted to add any.

`found_at` is taken from the object's own time rather than the wall clock at
insert. The old watcher used insert time, which is the same thing while it
keeps pace and a lie when it does not: catching up on a backlog would stack a
million points into the hour the catch-up ran and leave a spike in the
published hourly chart that no GPU ever produced. The object's time is also
stable across re-ingest, so a repeated object cannot move a bucket -- provided
it is the producer's. `producedAt` in the `.bin.json` commit marker is that;
`LastModified` is storage metadata, and a copy, a replication or a lifecycle
transition rewrites it, moving points between hourly buckets without a walker
doing anything. The marker wins where it carries one, the listing time stands
where it does not.

Two object-key shapes live under `dp/`, and both must be read:

    dp/slot-00002/1789311001-0000000000000000.bin
    dp/slot-00140/a8b4133d5c5f414ebd1337f15603588d-0000000000791392-9c6a7dae...e914.bin

The first is the original worker's `<upload epoch>-<offset>`. The second is
`ecc2k-seed-orbit-v1` (worker.py since #338): `<stream id>-<offset>-<sha256>`,
with no clock in the name, so its upload time is the object's `LastModified`.
The ingester matched only the first shape until 2026-09-17: the fleet moved to
the second at 10:06Z that day, the last old-shape slot uploaded at 13:50Z, and
the store then took nothing for five hours while 56 M records landed in `dp/`
-- with no error logged, because an object that fails the key pattern is not
an error, it is invisible.

Objects are ingested by several threads at once. One object at a time was
enough for the 15-slot fleet this was written for and is not enough for 133:
serialised, a pass sustained about 2,450 records/s against a fleet producing
2,600 points/s, so the store could not have kept pace even with every object
recognised.

Usage:
    dp_ingest.py --verify                 # check encodings, write nothing
    dp_ingest.py --once                   # one pass, then exit
    dp_ingest.py                          # poll forever
"""

from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import math
import os
import queue
import re
import struct
import sys
import threading
import time
import urllib.request

RECORD_BYTES = 32
CAMPAIGN = os.environ.get("RHO_CAMPAIGN", "ecc2k-130")
COEFF_BYTES = 17  # 17-byte big-endian mod 2^131, per rho_campaigns.meta
# dp/slot-00002/1789311001-0000000000000000.bin -- the original worker's
# <upload epoch>-<offset>.  The leading integer is the uploading worker's
# clock, which is where found_at comes from.
LEGACY_KEY_RE = re.compile(r"^dp/(slot-\d+)/(\d+)-(\d+)\.bin$")
# dp/slot-00140/<32-hex stream id>-<offset>-<64-hex sha256>.bin -- worker.py's
# ecc2k-seed-orbit-v1 naming.  It carries a content hash instead of a clock, so
# found_at is the object's own LastModified.
ORBIT_KEY_RE = re.compile(r"^dp/(slot-\d+)/([0-9a-f]{32})-(\d+)-([0-9a-f]{64})\.bin$")
# Written beside a record object by the contract path; metadata, not points.
ENVELOPE_SUFFIX = ".bin.json"
# Objects one pass will take before it stops to publish what it knows. Sized
# so that a pass is a couple of minutes at the measured rate (~3 objects/s),
# which is the cadence the dashboard and the alarms want; the backlog itself
# is unbounded and the next pass simply takes the next slice.
PASS_OBJECTS = 256
# Records decoded and COPYed at a time within one object. 250k decoded rows is
# about 140 MB of Python objects; see ingestObject for the object sizes that
# made a bound necessary.
INGEST_CHUNK_RECORDS = int(os.environ.get("RHO_INGEST_CHUNK_RECORDS", "250000"))


def log(msg):
    sys.stderr.write(time.strftime("%Y-%m-%dT%H:%M:%SZ ", time.gmtime()) + msg + "\n")
    sys.stderr.flush()


def workerId(key):
    """The identifier the existing rows use: the object key, slashes to dashes."""
    return "dp-" + key[len("dp/"):].replace("/", "-")


def foundAt(key, lastModified=None):
    """Upload time of one dp object, or None if the key is not a dp object.

    The two key shapes carry it differently and the difference is the whole
    reason this function exists: the legacy name has the uploader's clock in
    it, the orbit name has a content hash instead, and reading only the first
    shape is how five hours of a 133-worker fleet went into `dp/` without a
    single row reaching the store and without one line in the log.
    """
    m = LEGACY_KEY_RE.match(key)
    if m:
        return int(m.group(2))
    if ORBIT_KEY_RE.match(key) and lastModified is not None:
        return int(lastModified.timestamp())
    return None


def decode(record):
    """One 32-byte record -> the column values for one row."""
    seed, k0, k1, k2 = struct.unpack("<QQQQ", record)
    # walk_seed is the seed as a big-endian *integer*, so it carries no leading
    # zero bytes and its width varies with the seed -- the stored rows are 7
    # bytes wide wherever the seed happens to fit in 7.  Zero-padding to 8 would
    # decode to the same integer, but it would make byte-level comparison
    # against the existing corpus report differences that are not there.
    walk_seed = seed.to_bytes(8, "big").lstrip(b"\x00") or b"\x00"
    return {
        # Stored exactly as the client packs it: three little-endian uint64.
        "point_key": record[8:RECORD_BYTES],
        "a": seed.to_bytes(COEFF_BYTES, "big"),
        "b": (0).to_bytes(COEFF_BYTES, "big"),
        "walk_seed": walk_seed,
    }


def databaseUrl():
    url = os.environ.get("DATABASE_URL")
    if url:
        return url
    secret = os.environ.get("RHO_DB_SECRET", "rho/dp-rds")
    host = os.environ["RHO_DB_HOST"]
    import boto3

    blob = boto3.client("secretsmanager").get_secret_value(SecretId=secret)["SecretString"]
    c = json.loads(blob)
    sslmode = os.environ.get("RHO_DB_SSLMODE", "require")
    return "postgresql://%s:%s@%s:5432/%s?sslmode=%s" % (
        c["username"], urllib.request.quote(c["password"], safe=""),
        host, c["dbname"], sslmode)


def s3Objects(s3, bucket, prefix="dp/"):
    """Every dp object with its record count and upload time, oldest first.

    Returns `(objects, unrecognised)`. An object under `dp/` that matches
    neither key shape and is not an envelope is *not* skipped quietly: it is
    returned so the caller can say so. Silence about an unreadable key is the
    exact failure this ingest had, and a skip that logs nothing is
    indistinguishable from a fleet that has stopped.
    """
    out, unrecognised = [], []
    token = None
    while True:
        kw = {"Bucket": bucket, "Prefix": prefix}
        if token:
            kw["ContinuationToken"] = token
        page = s3.list_objects_v2(**kw)
        for item in page.get("Contents", []):
            key = item["Key"]
            when = foundAt(key, item.get("LastModified"))
            if when is None:
                if not key.endswith(ENVELOPE_SUFFIX):
                    unrecognised.append(key)
                continue
            records = item["Size"] // RECORD_BYTES
            if records:
                out.append((key, records, when))
        if not page.get("IsTruncated"):
            break
        token = page["NextContinuationToken"]
    out.sort(key=lambda kv: kv[2])
    return out, unrecognised


PROGRESS_DDL = """
CREATE TABLE IF NOT EXISTS dp_ingest_progress (
    campaign_id text NOT NULL,
    object_key  text NOT NULL,
    records     bigint NOT NULL,
    ingested_at timestamptz NOT NULL DEFAULT now(),
    PRIMARY KEY (campaign_id, object_key)
)
"""
# Records of an object whose point the store already held *under the same
# seed*. ON CONFLICT DO NOTHING drops them, correctly -- a resumed worker
# re-reporting its own point is what that clause is for -- but the same
# clause also drops a whole run that is walking another run's seeds, and did:
# Modal run 3 replayed AWS slot 2 for two days on 2026-09-19/20, 813k of its
# 912k records landed as silent re-reports, and the page showed a new worker
# adding points. This column is how that reads as what it is.
PROGRESS_DUPLICATES_DDL = """
ALTER TABLE dp_ingest_progress ADD COLUMN IF NOT EXISTS duplicates bigint NOT NULL DEFAULT 0
"""


def ensureProgress(conn):
    with conn.cursor() as cur:
        cur.execute(PROGRESS_DDL)
        cur.execute(PROGRESS_DUPLICATES_DDL)
    conn.commit()


FOUND_AT_INDEX = "distinguished_points_campaign_found_at"
META_DDL = """
CREATE TABLE IF NOT EXISTS dp_ingest_meta (
    key   text PRIMARY KEY,
    at    timestamptz NOT NULL DEFAULT now()
)
"""
VACUUM_MARK = "found_at_index_vacuum"


def ensureFoundAtIndex(conn):
    """Give the snapshot an index to read instead of the whole table.

    Every figure on the dashboard except `dps` is a question about `found_at`
    -- the newest point, the first point, the last 48 hours by hour -- and
    with no index on it each one reads all 42 GB of the table. Measured on
    2026-09-17 at 130 M rows: one `publishStatus` held ~12,000 read IOPS for
    over five minutes, on the same instance the ingest writes to, and
    `--status-every` would have started the next one straight after. `dps`
    itself becomes an index-only count of this index rather than a heap scan.

    CONCURRENTLY so the build does not block the ingest, which means it cannot
    run inside a transaction and can leave an invalid index behind if it
    fails; an invalid one is dropped and rebuilt rather than silently used.
    """
    previous = conn.autocommit
    conn.autocommit = True
    built = False
    try:
        with conn.cursor() as cur:
            cur.execute(META_DDL)
            cur.execute(
                "SELECT c.relname, i.indisvalid FROM pg_class c "
                "JOIN pg_index i ON i.indexrelid = c.oid WHERE c.relname = %s",
                (FOUND_AT_INDEX,))
            row = cur.fetchone()
            if row and not row[1]:
                log("index %s exists but is invalid (a previous build failed); "
                    "dropping it" % FOUND_AT_INDEX)
                cur.execute("DROP INDEX CONCURRENTLY IF EXISTS %s" % FOUND_AT_INDEX)
                row = None
            if not row:
                log("building index %s; the ingest keeps running while it does"
                    % FOUND_AT_INDEX)
                started = time.time()
                cur.execute(
                    "CREATE INDEX CONCURRENTLY IF NOT EXISTS %s "
                    "ON distinguished_points (campaign_id, found_at)" % FOUND_AT_INDEX)
                log("built index %s in %.0fs" % (FOUND_AT_INDEX, time.time() - started))
                built = True
            # An index-only scan reads the heap for every page the visibility
            # map does not mark all-visible, and a bulk-loaded table has
            # almost none marked: measured after the index was built, the
            # snapshot still took over six minutes, because it was still
            # reading 42 GB of heap. VACUUM is what sets that map, so the
            # index is only half the fix and this is the other half. Once
            # ever, recorded in dp_ingest_meta, because a replacement host is
            # the deployment mechanism here and a vacuum per deploy is not a
            # cost this table can carry.
            cur.execute("SELECT 1 FROM dp_ingest_meta WHERE key = %s", (VACUUM_MARK,))
            if cur.fetchone() is None:
                started = time.time()
                cur.execute("VACUUM (ANALYZE) distinguished_points")
                cur.execute("INSERT INTO dp_ingest_meta (key) VALUES (%s) "
                            "ON CONFLICT (key) DO NOTHING", (VACUUM_MARK,))
                log("vacuumed distinguished_points in %.0fs; index-only scans "
                    "are available from here" % (time.time() - started))
                return True
            return built
    finally:
        conn.autocommit = previous


PROGRESS_MARK = "progress_backfilled"


def ingestedObjects(conn):
    """Objects the store records as complete: `{object_key: records}`.

    This used to be two sources, because the corpus has two eras. Objects this
    program ingested have a `dp_ingest_progress` row. Objects the previous
    ingester wrote have none, and were recognised instead by counting rows per
    `worker_id` over the whole table -- 124 M rows grouped into 113 k on
    2026-09-17, three minutes before a pass could read its first object, and
    growing with the corpus. A 30-minute cache hid how often it ran; it did not
    change what it cost, and it was the last question this program asked of
    `distinguished_points` outside an insert.

    `backfillProgress` answers that question once and writes the rows the old
    ingester never wrote, so both eras now live in the same small table and a
    pass reads only this.
    """
    with conn.cursor() as cur:
        cur.execute(
            "SELECT object_key, records FROM dp_ingest_progress WHERE campaign_id = %s",
            (CAMPAIGN,))
        return {row[0]: int(row[1]) for row in cur.fetchall()}


def backfillProgress(conn, s3, bucket, prefix="dp/"):
    """Give the previous ingester's objects the progress rows they never had.

    Runs once, recorded in `dp_ingest_meta` beside the vacuum and rollup marks.
    The per-worker aggregate cannot be inverted -- `worker_id` is the object key
    with its slashes turned into dashes, and slot names contain dashes, so the
    key cannot be recovered from it -- and the bucket is the list of keys, so
    this pairs the two the same way `pending` used to on every pass.

    Only objects the store holds *in full* get a row. A partially ingested one
    is left without, so the next pass re-reads it, which is a no-op for the
    records already there. An object whose points were stored under some other
    era's worker_id (`gpu-i-<instance>-slot-N`) never reaches its record count
    and gets no row here either; it is re-ingested once, and that pass writes
    its row. Under the old scheme it was re-ingested on *every* pass forever,
    which is the non-convergence the progress table was added to fix.
    """
    with conn.cursor() as cur:
        cur.execute(META_DDL)
        cur.execute("SELECT 1 FROM dp_ingest_meta WHERE key = %s", (PROGRESS_MARK,))
        if cur.fetchone() is not None:
            return 0
        log("pairing the corpus with the bucket once to record what the previous "
            "ingester already stored; no pass reads the points table after this")
        started = time.time()
        cur.execute(
            "SELECT worker_id, count(*) FROM distinguished_points "
            "WHERE campaign_id = %s AND worker_id LIKE 'dp-slot-%%' GROUP BY worker_id",
            (CAMPAIGN,))
        counts = {row[0]: int(row[1]) for row in cur.fetchall()}
        objects, _ = s3Objects(s3, bucket, prefix)
        written = 0
        for key, records, _ in objects:
            if counts.get(workerId(key), 0) >= records:
                cur.execute(
                    "INSERT INTO dp_ingest_progress (campaign_id, object_key, records) "
                    "VALUES (%s, %s, %s) ON CONFLICT (campaign_id, object_key) "
                    "DO NOTHING", (CAMPAIGN, key, records))
                written += cur.rowcount
        cur.execute("INSERT INTO dp_ingest_meta (key) VALUES (%s) "
                    "ON CONFLICT (key) DO NOTHING", (PROGRESS_MARK,))
    conn.commit()
    log("recorded %d already-stored object(s) in %.0fs from %d worker groups"
        % (written, time.time() - started, len(counts)))
    return written


# Every figure the snapshot publishes except the checkpoint work was a
# question asked of the whole table -- count(*) over the campaign, min and max
# of found_at, and a 48-hour GROUP BY -- so its cost grew with the corpus
# while the answer it produced stayed the same size. It ran every
# --status-every seconds, on the instance the ingest writes to, and while it
# ran this program was not ingesting: the help text for --status-every already
# conceded "at ~190 M rows a publish takes a few minutes". An index and a
# one-off VACUUM bought time; neither changes the shape of the cost.
#
# These two tables are maintained by the insert itself, in the transaction
# that adds the points, so they cannot drift from what they summarise. The
# snapshot then reads one row and 48 more, and its cost stops depending on how
# long the campaign has been running.
#
# Every record of one object shares that object's found_at -- ingestObject
# inserts a single scalar for the whole COPY -- so an object contributes its
# rows to exactly one hour, and the rollup is an addition rather than a
# regrouping.
TOTALS_DDL = """
CREATE TABLE IF NOT EXISTS dp_ingest_totals (
    campaign_id text PRIMARY KEY,
    dps         bigint NOT NULL DEFAULT 0,
    first_dp_at timestamptz,
    last_dp_at  timestamptz
)
"""
HOURLY_DDL = """
CREATE TABLE IF NOT EXISTS dp_ingest_hourly (
    campaign_id text NOT NULL,
    hour        timestamptz NOT NULL,
    dps         bigint NOT NULL,
    PRIMARY KEY (campaign_id, hour)
)
"""
ROLLUP_MARK = "rollup_backfilled"


def applyRollup(cur, added, found_at):
    """Add one object's rows to the totals and to its hour.

    Called inside ingestObject's transaction, after the insert, with the rows
    that were actually inserted -- so a re-ingested object contributes zero and
    the counters stay exact under the idempotency the rest of this file
    depends on.
    """
    if added <= 0:
        return
    cur.execute(
        "INSERT INTO dp_ingest_totals (campaign_id, dps, first_dp_at, last_dp_at) "
        "VALUES (%s, %s, %s::timestamptz, %s::timestamptz) "
        "ON CONFLICT (campaign_id) DO UPDATE SET "
        "  dps = dp_ingest_totals.dps + EXCLUDED.dps, "
        "  first_dp_at = least(dp_ingest_totals.first_dp_at, EXCLUDED.first_dp_at), "
        "  last_dp_at = greatest(dp_ingest_totals.last_dp_at, EXCLUDED.last_dp_at)",
        (CAMPAIGN, added, found_at, found_at))
    cur.execute(
        "INSERT INTO dp_ingest_hourly (campaign_id, hour, dps) "
        "VALUES (%s, date_trunc('hour', %s::timestamptz), %s) "
        "ON CONFLICT (campaign_id, hour) DO UPDATE SET "
        "  dps = dp_ingest_hourly.dps + EXCLUDED.dps",
        (CAMPAIGN, found_at, added))


def ensureRollup(cur, force=False):
    """Seed the counters from the corpus, once, and say whether it ran.

    The corpus predates these tables, so the first pass has to read it -- the
    same aggregate the snapshot used to run, run once instead of every three
    minutes. One scan yields all of it: the hours are the rollup, their sum is
    the total, and the extremes of found_at come out of the same GROUP BY.

    Ordering is what makes this safe while the fleet is being ingested.
    `applyRollup` takes a row lock on this campaign's `dp_ingest_totals` row
    *after* inserting its points; locking that row here first means a
    concurrent ingest blocks before it adds its delta, with its points still
    uncommitted and therefore invisible to the aggregate below. It applies the
    delta once this commits. Neither a lost object nor a double-counted one is
    reachable.

    That argument has a precondition worth stating plainly, because it is the
    limit of what this can promise: it holds for writers that go through
    `applyRollup`. A superseded ingester still writing to the same table takes
    no such lock, so points it commits during or after the seed never reach
    these counters. The remedy is the order of a cutover -- stop the old
    writer, then seed -- and `--recount` for a corpus that was written to
    behind this program's back. It is deliberately not a self-healing check:
    verifying the seed means counting the table, and a verification that runs
    until it agrees with a moving corpus never agrees, so it would reinstate
    the per-pass full scan this whole mechanism exists to remove.
    """
    cur.execute(TOTALS_DDL)
    cur.execute(HOURLY_DDL)
    cur.execute(META_DDL)
    if force:
        cur.execute("DELETE FROM dp_ingest_meta WHERE key = %s", (ROLLUP_MARK,))
    else:
        cur.execute("SELECT 1 FROM dp_ingest_meta WHERE key = %s", (ROLLUP_MARK,))
        if cur.fetchone() is not None:
            return False
    log("seeding the dp counters from the corpus; this reads the table once "
        "and the ingest waits on it, then no snapshot reads it again")
    started = time.time()
    # Claim the lock before the read, so that anything mid-flight is behind us.
    cur.execute("INSERT INTO dp_ingest_totals (campaign_id) VALUES (%s) "
                "ON CONFLICT (campaign_id) DO NOTHING", (CAMPAIGN,))
    cur.execute("SELECT 1 FROM dp_ingest_totals WHERE campaign_id = %s FOR UPDATE",
                (CAMPAIGN,))
    cur.execute("LOCK TABLE dp_ingest_hourly IN EXCLUSIVE MODE")
    cur.execute(
        "SELECT date_trunc('hour', found_at), count(*), min(found_at), max(found_at) "
        "FROM distinguished_points WHERE campaign_id = %s GROUP BY 1", (CAMPAIGN,))
    rows = cur.fetchall()
    total = sum(int(n) for _, n, _, _ in rows)
    cur.execute("DELETE FROM dp_ingest_hourly WHERE campaign_id = %s", (CAMPAIGN,))
    for hour, n, _, _ in rows:
        cur.execute("INSERT INTO dp_ingest_hourly (campaign_id, hour, dps) "
                    "VALUES (%s, %s, %s)", (CAMPAIGN, hour, int(n)))
    cur.execute(
        "UPDATE dp_ingest_totals SET dps = %s, first_dp_at = %s, last_dp_at = %s "
        "WHERE campaign_id = %s",
        (total,
         min((r[2] for r in rows if r[2] is not None), default=None),
         max((r[3] for r in rows if r[3] is not None), default=None),
         CAMPAIGN))
    cur.execute("INSERT INTO dp_ingest_meta (key) VALUES (%s) "
                "ON CONFLICT (key) DO NOTHING", (ROLLUP_MARK,))
    log("seeded the dp counters in %.0fs: %d points across %d hours"
        % (time.time() - started, total, len(rows)))
    return True


# A collision is the campaign's terminal event, and until now the ingest was
# the one place it could be seen and the one place it was thrown away:
# `ON CONFLICT (campaign_id, point_key) DO NOTHING` cannot tell a worker
# re-reporting its own point from two different walks meeting, and dropped
# both alike. `rho_collisions` is read by this program and by
# scripts/rho_status/snapshot.py and was written by nothing in the tree, so
# the dashboard's COLLISION_RECORDED state had no producer at all.
#
# The store already enforces exactly the constraint that defines a collision.
# This makes it report one. merge.py stays the independent check over the S3
# corpus rather than the only detector.
COLLISIONS_DDL = """
CREATE TABLE IF NOT EXISTS rho_collisions (
    campaign_id text NOT NULL,
    point_key   bytea NOT NULL,
    a1          bytea,
    walk_seed1  bytea,
    worker_id1  text,
    a2          bytea,
    walk_seed2  bytea,
    worker_id2  text,
    object_key  text,
    detected_at timestamptz NOT NULL DEFAULT now(),
    PRIMARY KEY (campaign_id, point_key)
)
"""
# Columns the live table actually has. The DDL above is a no-op against the
# table the original out-of-tree ingest created, whose shape is not in this
# repository, so the INSERT is built from what is there rather than from what
# this file would have chosen. Anything but campaign_id and point_key is
# optional; without those two a row cannot be written at all and the finding
# goes to the log instead, which is never silent either way.
_collisionCols = set()


def ensureCollisions(conn):
    """Make sure a collision has somewhere to be recorded, and learn its shape."""
    global _collisionCols
    with conn.cursor() as cur:
        cur.execute(COLLISIONS_DDL)
        cur.execute("SELECT column_name FROM information_schema.columns "
                    "WHERE table_name = 'rho_collisions' "
                    "  AND table_schema = ANY(current_schemas(false))")
        _collisionCols = {row[0] for row in cur.fetchall()}
    conn.commit()
    missing = {"campaign_id", "point_key"} - _collisionCols
    if missing:
        log("WARNING: rho_collisions is missing %s; a collision will be logged "
            "but cannot be stored" % ", ".join(sorted(missing)))
    return _collisionCols


def seedOf(coeff):
    """The walk seed as an integer, however wide the stored bytes happen to be.

    `a` is the seed in 17 big-endian bytes here, but the corpus predates this
    program and verify() already reports that the old ingester's widths differ.
    Comparing bytes would call every legacy re-ingest a collision; comparing
    the integer is what the seed actually means.
    """
    return int.from_bytes(bytes(coeff or b""), "big")


def findCollisions(cur, key):
    """Sort this object's conflicting records into meetings and re-reports.

    Returns `(collisions, duplicates)`. Runs against the temp table after the
    insert, so every record of the object joins to the stored row for its
    point -- its own row when it was just inserted, someone else's when it was
    not. A point stored under a *different* seed is a collision, two walks
    meeting, which is the whole point of the campaign. A point stored under
    the *same* seed by a different object is a re-report: harmless from a
    resumed worker, and the entire output of a run that is walking another
    run's seeds, so it is counted rather than dropped in silence. The
    object's own rows are neither.
    """
    cur.execute(
        "SELECT i.point_key, i.a, i.walk_seed, d.a, d.walk_seed, d.worker_id "
        "FROM dp_in i JOIN distinguished_points d "
        "  ON d.campaign_id = %s AND d.point_key = i.point_key",
        (CAMPAIGN,))
    mine = workerId(key)
    out = []
    duplicates = 0
    for point_key, mineA, mineSeed, theirsA, theirsSeed, worker in cur.fetchall():
        if worker == mine:
            continue
        if theirsA is None:
            # `a` is NOT NULL for everything this program writes, but the
            # corpus predates it. A stored point with no seed cannot be
            # compared and could not be solved from if it were a meeting, so
            # it is said out loud and not counted -- reporting it as a
            # collision would make every record that meets it one.
            log("WARNING: point %s is stored with no seed (%s); cannot tell a "
                "re-report from a collision" % (bytes(point_key).hex(), worker))
            continue
        # Byte inequality would be the cheap filter; equality of the seeds
        # themselves is the question, and a width difference between eras is
        # not a collision.
        if seedOf(mineA) == seedOf(theirsA):
            duplicates += 1
            continue
        out.append({
            "point_key": bytes(point_key),
            "a1": bytes(theirsA),
            "walk_seed1": bytes(theirsSeed) if theirsSeed is not None else None,
            "worker_id1": worker,
            "a2": bytes(mineA),
            "walk_seed2": bytes(mineSeed),
            "worker_id2": mine,
            "object_key": key,
        })
    return out, duplicates


def recordCollisions(cur, found):
    """Store what can be stored, and say the rest out loud.

    Written in the same transaction as the points that revealed it, so a
    recorded collision cannot outlive the insert it came from. ON CONFLICT DO
    NOTHING because re-ingesting the same object must not multiply the row.
    """
    for row in found:
        # Never let a storage problem be the reason nobody hears about this.
        log("COLLISION: point %s seen from seed %d (%s) and seed %d (%s)"
            % (row["point_key"].hex(), seedOf(row["a1"]), row["worker_id1"],
               seedOf(row["a2"]), row["worker_id2"]))
    if not found or not {"campaign_id", "point_key"} <= _collisionCols:
        return 0
    cols = [c for c in ("point_key", "a1", "walk_seed1", "worker_id1",
                        "a2", "walk_seed2", "worker_id2", "object_key")
            if c in _collisionCols]
    sql = ("INSERT INTO rho_collisions (campaign_id, %s) VALUES (%s) "
           "ON CONFLICT DO NOTHING"
           % (", ".join(cols), ", ".join(["%s"] * (len(cols) + 1))))
    stored = 0
    for row in found:
        cur.execute(sql, tuple([CAMPAIGN] + [row[c] for c in cols]))
        stored += cur.rowcount
    return stored


def readEnvelope(s3, bucket, key):
    """The `.bin.json` commit marker beside an object, or None.

    The worker writes it *after* the payload, so its presence means the object
    is complete. It carries the producer's own sha256, record count and (since
    this change) producedAt. Until now the ingest ignored it entirely: merge.py
    verified it over the same corpus and this program did not, so a body that
    came back corrupted was accepted into the store that collision detection
    reads. Absent for the legacy key shape and for anything uploaded before the
    contract path, which is why nothing here treats a missing one as an error.
    """
    try:
        raw = s3.get_object(Bucket=bucket, Key=key + ".json")["Body"].read()
        manifest = json.loads(raw.decode("utf-8"))
        return manifest if isinstance(manifest, dict) else None
    except Exception:
        return None


def checkEnvelope(key, body, manifest):
    """Refuse a body that does not match its own commit marker.

    Truncation already self-heals -- a short read records fewer records than
    the object's S3 size, so the object stays outstanding and is retried -- but
    a body corrupted in place does not, and is indistinguishable from real
    points once stored. Raising here leaves the object outstanding and logged
    rather than letting it into the corpus.
    """
    if not manifest:
        return
    want = manifest.get("sha256")
    if want:
        got = hashlib.sha256(body).hexdigest()
        if got != want:
            raise ValueError("%s: body does not match its manifest sha256 "
                             "(%s != %s)" % (key, got, want))
    want = manifest.get("records")
    if isinstance(want, int) and want != len(body) // RECORD_BYTES:
        raise ValueError("%s: manifest says %d records, body carries %d"
                         % (key, want, len(body) // RECORD_BYTES))


def ingestObject(conn, s3, bucket, key, found_at):
    """Insert every record of one object.

    Returns (rows added, records seen, collisions found).
    """
    body = s3.get_object(Bucket=bucket, Key=key)["Body"].read()
    manifest = readEnvelope(s3, bucket, key)
    checkEnvelope(key, body, manifest)
    # The producer's own clock, where it recorded one. found_at from the
    # caller is the object's S3 LastModified for the orbit key shape, and that
    # is storage metadata: a copy, a replication or a lifecycle transition
    # rewrites it and silently moves points between hourly buckets.
    produced = manifest.get("producedAt") if manifest else None
    if isinstance(produced, int) and produced > 0:
        found_at = time.strftime("%Y-%m-%d %H:%M:%S+00", time.gmtime(produced))
    whole = len(body) - len(body) % RECORD_BYTES
    wid = workerId(key)
    added = 0
    collisions = []
    duplicates = 0
    with conn.cursor() as cur:
        # A temp table plus one INSERT ... ON CONFLICT is both fast and
        # idempotent; row-at-a-time INSERT was the old ingester's ceiling.
        cur.execute("CREATE TEMP TABLE IF NOT EXISTS dp_in "
                    "(point_key bytea, a bytea, b bytea, walk_seed bytea) ON COMMIT DROP")
        # One object, several COPY batches, one transaction. A worker's delta
        # is a few thousand records, but a sync loop whose state was reset
        # re-sends a whole corpus as one object -- 6.5 M records on
        # 2026-09-21, four times over -- and decoding that into Python dicts in
        # one go is ~3.5 GB; six threads of it is how the ingest died at 10:30
        # that day and stalled again at 14:00. The bytes stay whole (210 MB);
        # only the decoded rows are bounded.
        for start in range(0, whole, INGEST_CHUNK_RECORDS * RECORD_BYTES):
            end = min(whole, start + INGEST_CHUNK_RECORDS * RECORD_BYTES)
            cur.execute("TRUNCATE dp_in")
            # Sorted by the key the target is indexed on. point_key is a hash,
            # so an object's records arrive in random index order and each
            # insert walks to a different leaf page: on the 60 M-row table
            # that measured ~2.2 random reads per row inserted. Sorting costs
            # nothing here and turns a batch into a smaller set of pages
            # touched repeatedly.
            rows = sorted((decode(body[off:off + RECORD_BYTES])
                           for off in range(start, end, RECORD_BYTES)),
                          key=lambda r: r["point_key"])
            with cur.copy("COPY dp_in (point_key, a, b, walk_seed) FROM STDIN (FORMAT BINARY)") as cp:
                cp.set_types(["bytea", "bytea", "bytea", "bytea"])
                for r in rows:
                    cp.write_row((r["point_key"], r["a"], r["b"], r["walk_seed"]))
            cur.execute(
                "INSERT INTO distinguished_points "
                "  (campaign_id, point_key, a, b, walk_seed, worker_id, found_at) "
                "SELECT %s, point_key, a, b, walk_seed, %s, %s FROM dp_in "
                "ORDER BY point_key "
                "ON CONFLICT (campaign_id, point_key) DO NOTHING",
                (CAMPAIGN, wid, found_at))
            addedHere = cur.rowcount
            added += addedHere
            # Every record inserted means nothing conflicted, so there is
            # nothing to look for: in the steady state, where each object is
            # new points, the check below never runs and costs nothing. It
            # runs on the re-reports and the resumed walks -- and on the one
            # object that ends the campaign.
            if addedHere < len(rows):
                found, dupes = findCollisions(cur, key)
                recordCollisions(cur, found)
                collisions += found
                duplicates += dupes
            del rows
        if duplicates:
            log("%s: %d of %d records are re-reports of points the store already "
                "holds under the same seed" % (key, duplicates, whole // RECORD_BYTES))
        # Same transaction as the points, so the record of having ingested an
        # object cannot outlive the insert that it describes.
        cur.execute(
            "INSERT INTO dp_ingest_progress (campaign_id, object_key, records, duplicates) "
            "VALUES (%s, %s, %s, %s) ON CONFLICT (campaign_id, object_key) "
            "DO UPDATE SET records = EXCLUDED.records, "
            "  duplicates = EXCLUDED.duplicates, ingested_at = now()",
            (CAMPAIGN, key, whole // RECORD_BYTES, duplicates))
        # Also this transaction, so the counters the snapshot reads cannot
        # drift from the rows they summarise. Last, because every ingest
        # thread contends on the one totals row and the lock it takes is held
        # until the commit below.
        applyRollup(cur, added, found_at)
    conn.commit()
    return added, whole // RECORD_BYTES, len(collisions)


def verify(conn, s3, bucket, sample=64):
    """Check this program's encodings against points the old ingester wrote.

    Not a replication: it decodes the object independently and asks whether the
    resulting point_key is already stored under the same worker_id. Agreement
    means a resumed ingest lands in the same key space, which is the only
    property collision detection depends on.
    """
    done = ingestedObjects(conn)
    objects, _ = s3Objects(s3, bucket)
    complete = [(k, n) for k, n, _ in objects if done.get(k, 0) == n]
    if not complete:
        log("verify: no fully-ingested object to compare against")
        return False
    key, records = complete[len(complete) // 2]
    body = s3.get_object(Bucket=bucket, Key=key,
                         Range="bytes=0-%d" % (sample * RECORD_BYTES - 1))["Body"].read()
    keys = [decode(body[o:o + RECORD_BYTES])["point_key"]
            for o in range(0, len(body) - len(body) % RECORD_BYTES, RECORD_BYTES)]
    with conn.cursor() as cur:
        cur.execute(
            "SELECT point_key, a, b, walk_seed, worker_id FROM distinguished_points "
            "WHERE campaign_id = %s AND point_key = ANY(%s)", (CAMPAIGN, keys))
        rows = {bytes(r[0]): r for r in cur.fetchall()}
    hits = sum(1 for k in keys if k in rows)
    log("verify: %s -- %d/%d decoded point_keys already present"
        % (key, hits, len(keys)))
    if hits != len(keys):
        log("verify: FAILED. point_key encoding does not match the stored corpus; "
            "ingesting now would build a second, unmatchable key space")
        return False
    # The keys agree, so compare the auxiliary columns on one of them and report
    # rather than refuse: a/b/walk_seed are not in the primary key and the
    # collision path recomputes walks from the seed, so a difference here is
    # worth knowing about but does not break matching.
    k = keys[0]
    mine = decode(body[0:RECORD_BYTES])
    stored = rows[k]
    for i, col in enumerate(("a", "b", "walk_seed"), start=1):
        got, want = bytes(stored[i]), mine[col]
        log("verify: %-9s stored %d bytes, computed %d bytes, %s"
            % (col, len(got), len(want), "match" if got == want else "DIFFER"))
    log("verify: worker_id convention %r vs computed %r"
        % (stored[4], workerId(key)))
    return True


CKPT_MAGIC = b"ECC2K130"
CKPT_HEADER = 40  # magic[8] + 6x u32 + u64 iterBase
# Two live key shapes, the same split as dp/:
#   ckpt/slot-00002.ck and ckpt/retired/slot-00002.ck -- the original worker's
#   mutable pointer, overwritten in place.
#   ckpt/slot-00002/<64-hex>.ck -- ecc2k-seed-orbit-v1; the worker never
#   updates the old pointer, it writes a new immutable blob and keeps the
#   live name on the slot record.  Matching only the first shape is how a
#   walking fleet becomes zero walkers on the public feed.
CKPT_RE = re.compile(r"^ckpt/(retired/)?slot-(\d+)(?:/([0-9a-f]{64}))?\.ck$")

# Slots checkpoint every 600 s (campaign.json `checkpointEvery`). Three missed
# checkpoints is a dead worker rather than a slow one — the same rule as
# scripts/rho_status/work_feed.py, which strips per_slot before publishing.
FRESH_CHECKPOINT_S = 1800
# Below this the checkpoint granularity dominates the difference between two
# consecutive status.json writes.
MIN_WALK_RATE_SPAN_S = 600

# The campaign's distinguished-point rule, and what it costs in iterations per
# point: measured on the live fleet from each client's own counters, every GPU
# family agreeing to 0.04 in the exponent (aws/README.md, benchmarks/dp-interval).
# Walks stop at their own distinguished point, so a slot walking at another
# cutoff can meet a campaign walk and neither records the same point -- about
# 12% of meetings survive at weight 34, 4% at 35. The ratio of a slot's
# checkpointed iterations to the records it uploaded gives its cutoff away
# without reading a point, and that is what `dpWeightVerdict` reads. Weights
# 31 and 33 sit 1.6 either side of 32, so a tolerance of 1.0 separates
# neighbours; below DP_RATIO_MIN_RECORDS the ratio is mostly the lag between
# the checkpoint listing and the ingest.
CAMPAIGN_DP_WEIGHT = 32
CAMPAIGN_ITER_PER_DP_LOG2 = 28.41
DP_RATIO_TOLERANCE_LOG2 = 1.0
DP_RATIO_MIN_RECORDS = 50000
FIELD_BITS = 131


def theoreticalIterPerDpLog2(weight, m=FIELD_BITS):
    """Iterations per point for HW(x) <= weight over uniform m-bit strings."""
    total = 0
    for k in range(0, int(weight) + 1):
        total += math.comb(m, k)
    return m - math.log2(total)


def estimateDpWeight(log2IterPerDp, m=FIELD_BITS):
    """The cutoff that gives about this interval.

    The curve's x-coordinates run 0.6 shorter in the exponent than uniform
    strings at both measured weights (32: 28.41 vs 29.01; 34: 25.27 vs 25.84),
    so the binomial tail is shifted by that before the nearest weight is read.
    """
    shift = CAMPAIGN_ITER_PER_DP_LOG2 - theoreticalIterPerDpLog2(CAMPAIGN_DP_WEIGHT, m)
    best, bestGap = None, None
    for k in range(1, m):
        gap = abs(theoreticalIterPerDpLog2(k, m) + shift - log2IterPerDp)
        if bestGap is None or gap < bestGap:
            best, bestGap = k, gap
    return best


def dpWeightVerdict(iterations, records):
    """(log2 iterations per point, estimated weight, at campaign weight?).

    The last is True, False, or None when there is too little to judge.
    """
    if not records or int(records) < DP_RATIO_MIN_RECORDS or not iterations or iterations <= 0:
        return None, None, None
    ratio = math.log2(float(iterations) / float(records))
    ok = abs(ratio - CAMPAIGN_ITER_PER_DP_LOG2) <= DP_RATIO_TOLERANCE_LOG2
    return round(ratio, 3), estimateDpWeight(ratio), ok


SLOT_OF_KEY_RE = re.compile(r"^dp/slot-(\d+)/")


def slotOfKey(key):
    m = SLOT_OF_KEY_RE.match(key or "")
    return int(m.group(1)) if m else None


def coveredRecords(rows):
    """Records each slot has actually produced, from what it uploaded.

    `rows` are `(slot, stream, records, covered)`: the records summed over a
    slot's objects in one stream, and the furthest `offset / 32 + records`
    any of them reaches. Summing object sizes over-counts a stream that was
    uploaded twice -- modal_sync.py re-sent every Modal corpus from offset 0
    twice on 2026-09-21 when its state directory was reset, and the bucket
    then listed each Modal record about three times, which read as one point
    per 2^26.5 iterations and would have called runs at the campaign's own
    cutoff off-weight. Within one stream offsets are a file position and
    never overlap (worker.py starts a new stream when it rotates the file,
    modal_sync.py keeps one per run), so the furthest end is the count; the
    legacy key shape carries no stream and its era had no re-uploads, so it
    is summed.
    """
    out = {}
    for slot, stream, records, covered in rows:
        if slot is None:
            continue
        n = int(covered or 0) if stream else int(records or 0)
        out[int(slot)] = out.get(int(slot), 0) + n
    return out


def coverageRows(objects):
    """`(slot, stream, records, covered)` per stream from a bucket listing."""
    groups = {}
    for key, records, _ in objects:
        slot = slotOfKey(key)
        if slot is None:
            continue
        m = ORBIT_KEY_RE.match(key)
        stream = m.group(2) if m else None
        end = (int(m.group(3)) // RECORD_BYTES + int(records)) if m else 0
        acc = groups.setdefault((slot, stream), [0, 0])
        acc[0] += int(records)
        acc[1] = max(acc[1], end)
    return [(slot, stream, acc[0], acc[1]) for (slot, stream), acc in groups.items()]


def parseIso(value):
    if not value:
        return None
    text = str(value).replace("Z", "+00:00")
    try:
        stamp = datetime.datetime.fromisoformat(text)
    except ValueError:
        return None
    if stamp.tzinfo is None:
        stamp = stamp.replace(tzinfo=datetime.timezone.utc)
    return stamp


def loadPreviousStatus(s3, statusBucket):
    """The last published status.json, or None when this is the first write."""
    try:
        body = s3.get_object(Bucket=statusBucket, Key="status.json")["Body"].read()
        return json.loads(body.decode("utf-8"))
    except Exception:
        return None


def walkRateBetween(current, previous):
    """Iterations per second between two ingest publishes, or None.

    Consecutive writes are often closer together than MIN_WALK_RATE_SPAN_S
    (the dashboard publishes every ~180s; a slot checkpoints every 600s).
    When the previous document already carries a measured_from, difference
    against that older total instead of dropping the rate.
    """
    if not isinstance(previous, dict):
        return None
    bases = [previous]
    prior = previous.get("walk_rate")
    if isinstance(prior, dict) and prior.get("measured_from"):
        try:
            fromIter = int(prior.get("iterations_from"))
        except (TypeError, ValueError):
            fromIter = 0
        if fromIter > 0:
            bases.append({
                "generated_at": prior["measured_from"],
                "work": {"iterations": fromIter},
            })
    for base in bases:
        rate = _walkRate(current, base)
        if rate is not None:
            return rate
    return None


def _walkRate(current, previous):
    if not isinstance(previous, dict):
        return None
    curWork = current.get("work") or {}
    prevWork = previous.get("work") or {}
    try:
        curIter = int(curWork.get("iterations"))
        prevIter = int(prevWork.get("iterations"))
    except (TypeError, ValueError):
        return None
    if curIter <= 0 or prevIter <= 0 or curIter < prevIter:
        return None
    curAt = parseIso(current.get("generated_at"))
    prevAt = parseIso(previous.get("generated_at"))
    if curAt is None or prevAt is None:
        return None
    span = (curAt - prevAt).total_seconds()
    if span < MIN_WALK_RATE_SPAN_S:
        return None
    return {
        "iterations_per_second": (curIter - prevIter) / span,
        "window_seconds": int(round(span)),
        "measured_from": previous.get("generated_at"),
        "measured_to": current.get("generated_at"),
        "iterations_from": prevIter,
        "iterations_to": curIter,
        "method": "difference between this publish and the previous status.json",
    }


def windowSum(hourly, hours):
    """Points in the last `hours`, from the rollup.

    These two figures were a sliding `found_at > now() - interval '1 hour'`
    over the whole table. They are now a sum of whole hour buckets, which is
    the same series the published chart draws, so the tile and the chart can
    no longer disagree -- but it does mean the window is aligned to the hour
    and can reach back up to an hour further than the name suggests. For what
    they are used for, a COLLECTING/IDLE decision and a headline figure, that
    is the more consistent answer rather than a less accurate one.
    """
    now = datetime.datetime.now(datetime.timezone.utc)
    floor = (now - datetime.timedelta(hours=hours)).replace(
        minute=0, second=0, microsecond=0)
    total = 0
    for hour, count in hourly:
        if hour is None:
            continue
        if hour.tzinfo is None:
            hour = hour.replace(tzinfo=datetime.timezone.utc)
        if hour >= floor:
            total += int(count)
    return total


def publicStatus(payload):
    """The browser-safe document: same aggregates, no per-slot rows."""
    out = json.loads(json.dumps(payload))
    work = out.get("work")
    if isinstance(work, dict):
        work = dict(work)
        work.pop("per_slot", None)
        out["work"] = work
    return out


def checkpointWork(s3, bucket, staleAfter=1800.0):
    """Exact iterations walked, read from the slot checkpoints.

    This is the honest measure of work done and it needs no assumption about
    distinguished-point density.  Deriving work from the stored point count
    instead requires a constant iterations-per-point, and that constant changes
    whenever `dpWeight` changes -- so a campaign that has been retuned (as this
    one was, 34 -> 32 on 2026-09-13) has no single correct multiplier, and using
    the old one silently under-reports.

    The client prints "N iterations of M walks": `iterBase` is per-walk steps,
    not a total, so the work of one slot is `iterBase * walks`, and the client's
    `walksPerLaunch()` is `threads * BATCH * LANES` -- one for the packed
    engine, the word width for the bitsliced ones, and the header carries it.
    Until 2026-09-21 this read `threads * batch` and undercounted every
    bitsliced slot by its lane count (slot 195, 256 lanes, by 256x; about
    10^-4 of the campaign total, so an accounting correction and not a
    finding). Every factor comes out of the checkpoint header, so a slot is
    self-describing and a geometry change needs no bookkeeping here.  Retired
    slots are included -- their work is part of the campaign whether or not
    they still run.

    The contract path leaves every earlier blob in place, so a slot may have
    many `.ck` objects.  Only the newest one is the walk; summing the rest
    would count the same slot once per checkpoint.
    """
    latest, token = {}, None
    while True:
        kw = {"Bucket": bucket, "Prefix": "ckpt/"}
        if token:
            kw["ContinuationToken"] = token
        page = s3.list_objects_v2(**kw)
        for item in page.get("Contents", []):
            m = CKPT_RE.match(item["Key"])
            if not m:
                continue
            slot = int(m.group(2))
            modified = item["LastModified"].timestamp()
            prev = latest.get(slot)
            if prev is None or modified > prev[0]:
                latest[slot] = (modified, m, item)
        if not page.get("IsTruncated"):
            break
        token = page["NextContinuationToken"]
    slots = []
    for modified, m, item in latest.values():
        head = s3.get_object(Bucket=bucket, Key=item["Key"],
                             Range="bytes=0-%d" % (CKPT_HEADER - 1))["Body"].read()
        if len(head) < CKPT_HEADER or head[:8] != CKPT_MAGIC:
            log("checkpoint %s: not a checkpoint header, skipped" % item["Key"])
            continue
        version, m131, threads, batch, lanes, runId = struct.unpack_from("<6I", head, 8)
        iterBase, = struct.unpack_from("<Q", head, 32)
        walks = threads * batch * max(1, lanes)
        age = time.time() - modified
        slots.append({
            "slot": int(m.group(2)),
            "run_id": runId,
            # Whether a slot is still walking is decided by how recently its
            # checkpoint moved, not by the ckpt/retired/ prefix.  That prefix
            # is written by an orderly retirement, and a worker that was
            # stopped abruptly never writes it -- slots 0 and 1 sat under the
            # live prefix for hours after their instances were gone.  Asking
            # the artifact when it last changed cannot be fooled that way.
            "retired": bool(m.group(1)) or age > staleAfter,
            "checkpoint_age_s": int(age),
            "walks": walks,
            "per_walk_steps": iterBase,
            "iterations": iterBase * walks,
        })
    slots.sort(key=lambda s: (s["retired"], s["slot"]))
    total = sum(s["iterations"] for s in slots)
    return total, slots


def campaignSnapshot(conn):
    """SQL half of the dashboard: the insert-maintained rollup, not the table."""
    with conn.cursor() as cur:
        # One row from the counters instead of a count over the campaign, and
        # 48 from the rollup instead of a GROUP BY over two days of points.
        # Neither reads distinguished_points, so the cost of a publish no
        # longer grows with the corpus.
        cur.execute("""
            SELECT c.campaign_id, c.curve_id, c.dp_mask_bits, c.created_at,
                   COALESCE(t.dps, 0) AS dps, t.first_dp_at, t.last_dp_at
            FROM rho_campaigns c
            LEFT JOIN dp_ingest_totals t USING (campaign_id)
            WHERE c.campaign_id = %s""", (CAMPAIGN,))
        row = cur.fetchone() or ()
        cur.execute("SELECT count(*), max(detected_at) FROM rho_collisions WHERE campaign_id = %s",
                    (CAMPAIGN,))
        collisions, latest_collision = cur.fetchone()
        cur.execute("""
            SELECT hour, dps FROM dp_ingest_hourly
            WHERE campaign_id = %s
              AND hour >= date_trunc('hour', now() - interval '48 hours')
            ORDER BY hour""", (CAMPAIGN,))
        hourly = cur.fetchall()
        # One row per (slot, upload stream) from the progress table -- one row
        # per object there, tens of thousands, never the points: what each
        # slot uploaded, how far its stream reaches (see coveredRecords for
        # why that and not the sum), and how much of it the store already
        # held. `statusPayload` divides the checkpointed iterations by the
        # covered records to read the slot's cutoff.
        cur.execute(r"""
            SELECT substring(object_key from '^dp/slot-([0-9]+)/') AS slot,
                   substring(object_key from '^dp/slot-[0-9]+/([0-9a-f]{32})-[0-9]+-[0-9a-f]{64}\.bin$') AS stream,
                   sum(records),
                   max(substring(object_key from '^dp/slot-[0-9]+/[0-9a-f]{32}-([0-9]+)-[0-9a-f]{64}\.bin$')::bigint / 32 + records),
                   sum(duplicates),
                   sum(records) FILTER (WHERE ingested_at > now() - interval '24 hours'),
                   sum(duplicates) FILTER (WHERE ingested_at > now() - interval '24 hours')
            FROM dp_ingest_progress
            WHERE campaign_id = %s
            GROUP BY 1, 2""", (CAMPAIGN,))
        streams = cur.fetchall()
        covered = coveredRecords([(s, st, n, c) for s, st, n, c, _, _, _ in streams])
        perSlot = {}
        for slot, _, records, _, duplicates, recordsDay, duplicatesDay in streams:
            if slot is None:
                continue
            entry = perSlot.setdefault(int(slot), {
                "records": covered.get(int(slot), 0), "uploaded": 0, "duplicates": 0,
                "records_last_day": 0, "duplicates_last_day": 0})
            entry["uploaded"] += int(records or 0)
            entry["duplicates"] += int(duplicates or 0)
            entry["records_last_day"] += int(recordsDay or 0)
            entry["duplicates_last_day"] += int(duplicatesDay or 0)

    def iso(v):
        return v.isoformat() if hasattr(v, "isoformat") else v

    return {
        "per_slot_records": perSlot,
        "curve_id": row[1] if row else None,
        "dp_mask_bits": row[2] if row else None,
        "campaign_created_at": iso(row[3]) if row else None,
        "dps": int(row[4] or 0) if row else 0,
        "dps_last_hour": windowSum(hourly, 1),
        "dps_last_day": windowSum(hourly, 24),
        "first_dp_at": iso(row[5]) if row else None,
        "last_dp_at": iso(row[6]) if row else None,
        "collisions": int(collisions or 0),
        "latest_collision_at": iso(latest_collision),
        "hourly": [{"hour": iso(h), "dps": int(n)} for h, n in hourly],
    }


def statusPayload(conn, s3, bucket, ingest=None, snapshot=None):
    """The dashboard snapshot, computed here so the page can read it directly.

    `ingest` is what the last pass knew about itself. It is published because
    `state` alone cannot tell the difference between a fleet that stopped and
    an ingest that stopped: on 2026-09-17 the page read IDLE_OR_STALE for five
    hours while 133 workers were walking and 56 M records were landing in S3.

    `snapshot` is the SQL half (point counts from the insert-maintained
    rollup). Checkpoint work is always re-read from S3, so walker liveness
    can move every --status-every without waiting on the database.
    """
    if snapshot is None:
        snapshot = campaignSnapshot(conn)
    ingest = dict(ingest or {})
    dps = int(snapshot.get("dps") or 0)
    dps_last_hour = int(snapshot.get("dps_last_hour") or 0)
    collisions = int(snapshot.get("collisions") or 0)
    # An ingest that is behind makes every recency figure below a statement
    # about this program, not about the fleet, so it is named as one rather
    # than left to be read as a quiet campaign.
    behind = int(ingest.get("outstanding", 0)) or int(ingest.get("unrecognised", 0))
    state = ("COLLISION_RECORDED" if collisions else
             "COLLECTING" if dps_last_hour else
             "INGEST_BEHIND" if behind else
             "IDLE_OR_STALE" if dps else "EMPTY")
    iterations, per_slot = checkpointWork(s3, bucket)
    walking = [
        slot for slot in per_slot
        if not slot["retired"] and slot["checkpoint_age_s"] <= FRESH_CHECKPOINT_S
    ]
    # Each slot against the campaign's cutoff, and what it re-reported. The
    # per-slot rows stay private (publicStatus strips them); the counts are
    # what the page gets, and they are the two numbers that would have named
    # both 2026-09-20 Modal problems within an hour.
    perSlotRecords = snapshot.get("per_slot_records") or {}
    for slot in per_slot:
        counts = perSlotRecords.get(slot["slot"]) or {}
        slot["records"] = int(counts.get("records") or 0)
        slot["uploaded"] = int(counts.get("uploaded") or slot["records"])
        slot["duplicates"] = int(counts.get("duplicates") or 0)
        slot["duplicates_last_day"] = int(counts.get("duplicates_last_day") or 0)
        ratio, estimate, ok = dpWeightVerdict(slot["iterations"], slot["records"])
        slot["iterations_per_dp_log2"] = ratio
        slot["dp_weight_estimate"] = estimate
        slot["dp_weight_ok"] = ok
    offWeight = [slot for slot in per_slot if slot["dp_weight_ok"] is False]
    offWeightWalking = [slot for slot in offWeight if slot in walking]
    duplicateRecords = sum(c.get("duplicates") or 0 for c in perSlotRecords.values())
    duplicateRecordsDay = sum(c.get("duplicates_last_day") or 0 for c in perSlotRecords.values())
    duplicateSlotsDay = sum(1 for c in perSlotRecords.values() if c.get("duplicates_last_day"))
    payload = {
        "schema_version": 1,
        "generated_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "source": "dp_ingest.py (live, direct from rho-dp + slot checkpoints)",
        "campaign_id": CAMPAIGN,
        "curve_id": snapshot.get("curve_id"),
        "dp_mask_bits": snapshot.get("dp_mask_bits"),
        "campaign_created_at": snapshot.get("campaign_created_at"),
        "state": state,
        "dps": dps,
        "collisions": collisions,
        "first_dp_at": snapshot.get("first_dp_at"),
        "last_dp_at": snapshot.get("last_dp_at"),
        "latest_collision_at": snapshot.get("latest_collision_at"),
        "dps_last_hour": dps_last_hour,
        "dps_last_day": int(snapshot.get("dps_last_day") or 0),
        "hourly": list(snapshot.get("hourly") or []),
        # Exact, and the reason this block exists rather than a dps multiplier:
        # see checkpointWork.  `walkers` is the count of slots still running,
        # which is what "workers" should have meant all along -- the page's own
        # figure counts distinct worker_id values, and those are per-object
        # names, so it grows with every upload.
        "work": {
            "iterations": iterations,
            "iterations_log2": (math.log2(iterations) if iterations else None),
            "method": "sum over slot checkpoints of iterBase * threads * batch",
            "density_independent": True,
            "slots": len(per_slot),
            "walking_slots": len(walking),
            # Slots whose iterations-per-point says they are not walking at
            # the campaign's cutoff: their points mostly cannot meet anyone
            # else's, however many they upload.
            "campaign_dp_weight": CAMPAIGN_DP_WEIGHT,
            "iterations_per_dp_log2_expected": CAMPAIGN_ITER_PER_DP_LOG2,
            "off_weight_slots": len(offWeight),
            "off_weight_walking_slots": len(offWeightWalking),
            "per_slot": per_slot,
        },
        "walkers": sum(1 for s in per_slot if not s["retired"]),
        # What this program knows about its own health, so that a reader of
        # the page can tell a stopped fleet from a stopped ingest.
        "ingest": {
            "outstanding_objects": int(ingest.get("outstanding", 0)),
            "unrecognised_objects": int(ingest.get("unrecognised", 0)),
            "newest_object_at": (
                time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime(ingest["newest"]))
                if ingest.get("newest") else None),
            "lag_seconds": (int(time.time() - ingest["newest"])
                            if ingest.get("newest") else None),
            # Records dropped as the same walk's point again. A resumed
            # worker produces a few; a run walking another run's seeds
            # produces nothing else.
            "duplicate_records": duplicateRecords,
            "duplicate_records_last_day": duplicateRecordsDay,
            "duplicate_slots_last_day": duplicateSlotsDay,
        },
        "claim_boundary": (
            "Public research campaign aggregates for Certicom ECC2K-130 "
            "distinguished-point collection. Not a discrete-log recovery, not a "
            "verified collision, and not a statement about deployed-curve security "
            "unless collisions>0 and an independent solver has checked [k]P = Q."),
    }
    blob = json.dumps(payload, sort_keys=True)
    for banned in ("point_key", "walk_seed", "password", "secret"):
        if banned in blob:
            raise RuntimeError("refusing to publish field %s" % banned)
    return payload


def publishMetrics(namespace, ingest):
    """Put the ingest's own health where an alarm can see it.

    Best effort by design: a metric this program cannot publish is not a
    reason to stop ingesting, and the same numbers are in the log and in
    status.json.
    """
    try:
        import boto3

        data = [{"MetricName": "OutstandingObjects",
                 "Value": float(ingest.get("outstanding", 0)), "Unit": "Count"},
                {"MetricName": "UnrecognisedObjects",
                 "Value": float(ingest.get("unrecognised", 0)), "Unit": "Count"}]
        if ingest.get("newest"):
            data.append({"MetricName": "NewestObjectAgeSeconds",
                         "Value": float(time.time() - ingest["newest"]), "Unit": "Seconds"})
        boto3.client("cloudwatch").put_metric_data(Namespace=namespace, MetricData=data)
    except Exception as exc:
        log("metric publish failed: %s: %s" % (type(exc).__name__, exc))


def publishStatus(conn, s3, bucket, statusBucket, ingest=None, snapshot=None,
                  source=None):
    """Write the snapshot where a browser can read it, no workflow involved.

    Point counts come from `snapshot` (the SQL query). Walker liveness always
    comes from S3 checkpoints, so a cached snapshot still moves `walkers`.
    """
    started = time.time()
    previous = loadPreviousStatus(s3, statusBucket)
    kind = source or ("cached" if snapshot is not None else "sql")
    payload = statusPayload(conn, s3, bucket, ingest, snapshot=snapshot)
    payload["published_at"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    rate = walkRateBetween(payload, previous)
    if rate is not None:
        payload["walk_rate"] = rate
    public = publicStatus(payload)
    s3.put_object(
        Bucket=statusBucket, Key="status.json",
        Body=(json.dumps(public, indent=2, sort_keys=True) + "\n").encode(),
        ContentType="application/json",
        # Short but non-zero: the underlying data only moves when a worker
        # uploads, so caching for less than that buys nothing and costs requests.
        CacheControl="public, max-age=30")
    # The SQL half now reads the counters the insert maintains, so its cost
    # should stay flat; walker refresh is a list of ckpt/ objects. If the
    # elapsed time starts climbing, the counters are not being read.
    log("published status.json in %.1fs (%s): dps=%d state=%s work=2^%.3f walkers=%d "
        "outstanding=%d unreadable=%d off_weight_walking=%d duplicates_24h=%d"
        % (time.time() - started, kind, payload["dps"], payload["state"],
           payload["work"]["iterations_log2"] or 0, payload["walkers"],
           payload["ingest"]["outstanding_objects"],
           payload["ingest"]["unrecognised_objects"],
           payload["work"]["off_weight_walking_slots"],
           payload["ingest"]["duplicate_records_last_day"]))
    for slot in payload["work"]["per_slot"]:
        if slot.get("dp_weight_ok") is False and not slot["retired"]:
            log("WARNING: slot %05d run %d reports one point per 2^%.2f iterations; the "
                "campaign's weight-%d cutoff gives 2^%.2f, so it looks like weight %s and "
                "its meetings with campaign walks mostly go unrecorded"
                % (slot["slot"], slot["run_id"], slot["iterations_per_dp_log2"],
                   CAMPAIGN_DP_WEIGHT, CAMPAIGN_ITER_PER_DP_LOG2, slot["dp_weight_estimate"]))
    return payload


class StatusPublisher:
    """Publish status.json on its own cadence, never blocked by an ingest pass.

    The main loop used to be pass-then-publish. Status and metrics only ran
    after `onePass` returned, so a slow `pending()` aggregate or a stuck
    object drain froze the public feed for as long as the pass took — the
    page sat on IDLE_OR_STALE with a generated_at that did not move, which
    is indistinguishable from a dead ingest host. This thread owns the
    cadence: it re-reads checkpoints from S3 every `--status-every` seconds
    and refreshes the SQL rollup on `--snapshot-every`, using the latest
    ingest health the main loop has handed it.
    """

    def __init__(self, connect, s3, bucket, statusBucket, status_every=180.0,
                 snapshot_every=0.0):
        self.connect = connect
        self.s3 = s3
        self.bucket = bucket
        self.statusBucket = statusBucket
        self.status_every = max(1.0, float(status_every))
        self.snapshot_every = float(snapshot_every)
        self._lock = threading.Lock()
        self._ingest = {}
        self._snapshot = None
        self._snap_at = 0.0
        self._wake = threading.Event()
        self._stop = threading.Event()
        self._thread = None

    def update_ingest(self, ingest):
        with self._lock:
            self._ingest = dict(ingest or {})

    def kick(self):
        """Publish as soon as the thread is free (caught-up / collision)."""
        self._wake.set()

    def start(self):
        if self._thread is not None:
            return
        self._thread = threading.Thread(
            target=self._loop, name="status-publisher", daemon=True)
        self._thread.start()

    def stop(self, timeout=5.0):
        self._stop.set()
        self._wake.set()
        thread = self._thread
        if thread is not None:
            thread.join(timeout=timeout)

    def publish_once(self):
        """One publish on the calling thread (for --once and tests)."""
        self._publish()

    def _loop(self):
        # First publish immediately so a restart is not silent for status_every.
        self._publish()
        while not self._stop.is_set():
            self._wake.wait(timeout=self.status_every)
            self._wake.clear()
            if self._stop.is_set():
                return
            self._publish()

    def _publish(self):
        with self._lock:
            ingest = dict(self._ingest)
            snapshot = self._snapshot
            snap_at = self._snap_at
        try:
            need_sql = (
                snapshot is None
                or self.snapshot_every <= 0
                or time.time() - snap_at >= self.snapshot_every
            )
            if need_sql:
                with self.connect() as conn:
                    snapshot = campaignSnapshot(conn)
                snap_at = time.time()
                with self._lock:
                    self._snapshot = snapshot
                    self._snap_at = snap_at
            publishStatus(
                None, self.s3, self.bucket, self.statusBucket,
                ingest, snapshot=snapshot,
                source="sql" if need_sql else "cached",
            )
        except Exception as exc:
            # Publishing is a view; never let it stop the ingest or the thread.
            log("status publish failed: %s: %s" % (type(exc).__name__, exc))


def pending(conn, s3, bucket, prefix="dp/"):
    """Objects in the bucket that the store does not have, oldest first.

    Also returns the newest upload time seen and any key it could not read, so
    that being behind and being unable to read are two different, reported
    numbers rather than one silence.
    """
    done = ingestedObjects(conn)
    objects, unrecognised = s3Objects(s3, bucket, prefix)
    todo, newest = [], 0
    for key, records, when in objects:
        newest = max(newest, when)
        have = done.get(key, 0)
        if have >= records:
            continue
        todo.append((key, records, when, have))
    if unrecognised:
        log("WARNING: %d object(s) under dp/ match no known key shape and are "
            "NOT being ingested, for example: %s"
            % (len(unrecognised), ", ".join(unrecognised[:3])))
    return todo, newest, unrecognised


# Keys that failed the last bounded pass. They stay at the front of
# pending(), so taking the oldest `limit` again would pin the window on
# them and never reach the objects behind. The next slice skips them;
# they are retried once they no longer fill the window, or once they are
# all that remains.
_failed = set()


def onePass(connect, s3, bucket, threads=6, prefix="dp/", limit=PASS_OBJECTS):
    """Ingest the oldest outstanding objects, several at a time.

    Each thread owns a connection: psycopg connections are not shared, and the
    temp table the COPY lands in is per-session anyway. Order does not matter
    because every insert is idempotent. A failed object is retried on a later
    pass; the next slice skips it so a poison prefix cannot pin the bound
    window on the same oldest keys.

    A pass is *bounded* rather than "everything outstanding". Status used to
    publish only between passes, so an unbounded drain froze the page for as
    long as it took; `StatusPublisher` now owns that cadence on its own
    thread. `outstanding` stays the whole backlog, not this slice, so the
    number a reader sees is the one that matters.
    """
    with connect() as probe:
        todo, newest, unrecognised = pending(probe, s3, bucket, prefix)
    state = {"newest": newest, "unrecognised": len(unrecognised), "outstanding": len(todo)}
    if not todo:
        _failed.clear()
        return 0, 0, state
    if limit and limit > 0:
        rest = [item for item in todo if item[0] not in _failed]
        slice_ = (rest or todo)[:limit]
    else:
        slice_ = todo
    work = queue.Queue()
    for item in slice_:
        work.put(item)
    tally = {"rows": 0, "objects": 0, "failed": 0, "collisions": 0}
    failed = set()
    lock = threading.Lock()

    def drain():
        conn = None
        try:
            while True:
                try:
                    key, records, when, have = work.get_nowait()
                except queue.Empty:
                    return
                # Two attempts, because the failure worth surviving is the
                # database going away mid-pass -- a failover or an instance
                # resize -- after which this thread's connection is closed and
                # every remaining object would fail against the corpse of it.
                # 2,533 objects failed that way on 2026-09-17 when rho-dp was
                # resized under a running drain.
                for attempt in (0, 1):
                    try:
                        if conn is None or conn.closed:
                            conn = connect()
                        n, seen, hits = ingestObject(
                            conn, s3, bucket, key,
                            time.strftime("%Y-%m-%d %H:%M:%S+00", time.gmtime(when)))
                    except Exception as exc:
                        broken = conn is None or conn.closed
                        if not broken:
                            try:
                                conn.rollback()
                            except Exception:
                                conn, broken = None, True
                        else:
                            conn = None
                        if broken and attempt == 0:
                            log("reconnecting after %s on %s"
                                % (type(exc).__name__, key))
                            time.sleep(2.0)
                            continue
                        # One bad object must not take the pass down with it;
                        # it is skipped on the next slice so it cannot stall
                        # the rest of the backlog, then retried later.
                        with lock:
                            tally["failed"] += 1
                            failed.add(key)
                        log("object %s failed (will retry): %s: %s"
                            % (key, type(exc).__name__, exc))
                    else:
                        with lock:
                            tally["rows"] += n
                            tally["objects"] += 1
                            tally["collisions"] += hits
                        log("ingested %s: %d records, %d new (store had %d)%s"
                            % (key, seen, n, have,
                               ", %d COLLISION(S)" % hits if hits else ""))
                    break
        except Exception as exc:  # a connection that will not open at all
            log("ingest thread stopped: %s: %s" % (type(exc).__name__, exc))
        finally:
            if conn is not None and not conn.closed:
                conn.close()

    pool = [threading.Thread(target=drain, daemon=True) for _ in range(max(1, threads))]
    started = time.time()
    for t in pool:
        t.start()
    for t in pool:
        t.join()
    _failed.clear()
    _failed.update(failed)
    elapsed = max(time.time() - started, 1e-9)
    state["outstanding"] = len(todo) - tally["objects"]
    state["collisions"] = tally["collisions"]
    log("pass complete: %d of %d outstanding objects touched, %d rows added, "
        "%d failed, %d still outstanding, %.0f rows/s%s"
        % (tally["objects"], len(todo), tally["rows"], tally["failed"],
           state["outstanding"], tally["rows"] / elapsed,
           ", %d COLLISION(S) RECORDED" % tally["collisions"]
           if tally["collisions"] else ""))
    return tally["rows"], tally["objects"], state


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--bucket", default=os.environ.get("RHO_BUCKET"))
    ap.add_argument("--prefix", default="dp/",
                    help="restrict to one slot ('dp/slot-00140/') to re-ingest or rehearse")
    ap.add_argument("--status-bucket", default=os.environ.get("RHO_STATUS_BUCKET"),
                    help="publish status.json here for the page to poll directly")
    ap.add_argument("--interval", type=float, default=120.0)
    ap.add_argument("--threads", type=int,
                    default=int(os.environ.get("RHO_INGEST_THREADS", "6")),
                    help="objects ingested concurrently, each on its own connection")
    ap.add_argument("--metric-namespace", default=os.environ.get("RHO_METRIC_NAMESPACE"),
                    help="publish outstanding/lag/unrecognised here as CloudWatch metrics")
    ap.add_argument("--pass-objects", type=int,
                    default=int(os.environ.get("RHO_PASS_OBJECTS", PASS_OBJECTS)),
                    help="objects per pass; status publishes on its own thread "
                         "so a long drain no longer freezes the feed "
                         "(0 = unbounded)")
    ap.add_argument("--no-index", dest="index", action="store_false",
                    help="do not create the found_at index at startup")
    ap.add_argument("--status-every", type=float,
                    default=float(os.environ.get("RHO_STATUS_EVERY", "180")),
                    help="seconds between status.json writes on the publisher "
                         "thread. The SQL half reads the maintained counters "
                         "rather than the points table; walker liveness is an "
                         "S3 list. Both keep moving while a pass runs. The "
                         "page's Actions job runs every 3 minutes and reads "
                         "this copy from the status bucket.")
    ap.add_argument("--snapshot-every", type=float,
                    default=float(os.environ.get("RHO_SNAPSHOT_EVERY", "0")),
                    help="seconds between rollup reads. 0 (the default) reads "
                         "them on every status write. Walker liveness still "
                         "refreshes every --status-every.")
    ap.add_argument("--recount", action="store_true",
                    help="re-seed the dp counters from the table even if they "
                         "were seeded before. They are maintained in the same "
                         "transaction as the points and cannot drift from them, "
                         "so this is for a corpus something else has written to "
                         "-- a superseded ingester that wrote during the first "
                         "seed. Run it once, after that writer is stopped; it "
                         "reads the whole table, and it is consumed on the "
                         "first seed rather than repeated on every reconnect")
    ap.add_argument("--once", action="store_true")
    ap.add_argument("--verify", action="store_true")
    ap.add_argument("--work", action="store_true",
                    help="print exact iterations from the checkpoints and exit")
    ap.add_argument("--pending", action="store_true",
                    help="report what the store is missing and exit, writing nothing")
    args = ap.parse_args(argv)
    if not args.bucket:
        raise SystemExit("--bucket or RHO_BUCKET is required")

    import boto3

    s3 = boto3.client("s3")

    if args.work:  # needs no database, so the driver is not imported yet
        total, per_slot = checkpointWork(s3, args.bucket)
        # Records per slot from the bucket listing rather than the store, so
        # this stays a no-database report; the ratio reads the cutoff either way.
        objects, _ = s3Objects(s3, args.bucket, args.prefix)
        recordsBySlot = coveredRecords(coverageRows(objects))
        for s in per_slot:
            ratio, estimate, ok = dpWeightVerdict(s["iterations"], recordsBySlot.get(s["slot"], 0))
            verdict = ("" if ok is None else
                       "  2^%.2f it/point %s" % (ratio, "ok" if ok else "OFF WEIGHT (~%d)" % estimate))
            log("slot %05d run %d %-8s ckpt %5ds old  %d walks x %d steps = %.6g iterations%s"
                % (s["slot"], s["run_id"], "retired" if s["retired"] else "walking",
                   s["checkpoint_age_s"], s["walks"], s["per_walk_steps"], s["iterations"],
                   verdict))
        log("exact total: %d iterations = 2^%.3f (no density assumption)"
            % (total, math.log2(total) if total else 0))
        return 0

    import psycopg

    url = databaseUrl()

    def connect():
        return psycopg.connect(url, connect_timeout=30)

    recount_once = args.recount
    while True:
        try:
            with connect() as conn:
                ensureProgress(conn)
                # --pending promises to write nothing, and the DDL is a write.
                if not (args.verify or args.pending):
                    ensureCollisions(conn)
                    backfillProgress(conn, s3, args.bucket, args.prefix)
                    with conn.cursor() as cur:
                        ensureRollup(cur, force=recount_once)
                    conn.commit()
                    # Consume it: --recount stays set for the life of the
                    # process, and this block is inside the failover retry, so
                    # leaving it true made a database blip re-read the whole
                    # corpus on every reconnect.
                    recount_once = False
                if args.index and not (args.verify or args.pending):
                    try:
                        ensureFoundAtIndex(conn)
                    except Exception as exc:
                        # A missing index is slow, not wrong: say so and walk on.
                        log("index build failed (snapshots stay slow): %s: %s"
                            % (type(exc).__name__, exc))
                if args.verify:
                    return 0 if verify(conn, s3, args.bucket) else 1
                if args.pending:
                    todo, newest, unrecognised = pending(conn, s3, args.bucket, args.prefix)
                    log("outstanding: %d objects, %d records; newest object %s; "
                        "%d unreadable keys"
                        % (len(todo), sum(n for _, n, _, _ in todo),
                           time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime(newest))
                           if newest else "none",
                           len(unrecognised)))
                    return 0

            publisher = None
            if args.status_bucket and not args.once:
                publisher = StatusPublisher(
                    connect, s3, args.bucket, args.status_bucket,
                    status_every=args.status_every,
                    snapshot_every=args.snapshot_every,
                )
                publisher.start()
            wasBehind = False
            try:
                while True:
                    _, objects, ingest = onePass(
                        connect, s3, args.bucket, args.threads,
                        args.prefix, args.pass_objects)
                    if args.metric_namespace:
                        publishMetrics(args.metric_namespace, ingest)
                    behind = bool(ingest.get("outstanding"))
                    caught_up = wasBehind and not behind
                    wasBehind = behind
                    if publisher is not None:
                        publisher.update_ingest(ingest)
                        # A drain that finishes, or a collision, is the
                        # transition a reader is waiting for — do not wait
                        # out the rest of status_every.
                        if caught_up or bool(ingest.get("collisions")):
                            publisher.kick()
                    elif args.status_bucket:
                        # --once: one synchronous publish after the pass.
                        try:
                            with connect() as conn:
                                snapshot = campaignSnapshot(conn)
                            publishStatus(
                                None, s3, args.bucket, args.status_bucket,
                                ingest, snapshot=snapshot, source="sql")
                        except Exception as exc:
                            log("status publish failed: %s: %s"
                                % (type(exc).__name__, exc))
                    if args.once:
                        return 0
                    # A backlog is drained as fast as the database allows rather
                    # than one pass per interval: the interval exists to keep an
                    # idle ingest cheap, not to rate-limit catching up.  A pass
                    # that added nothing is not a backlog -- it is a retry -- and
                    # sleeping 0 while a poison object stays outstanding is how
                    # one failed row busy-loops the host.
                    time.sleep(0.0 if objects else args.interval)
            finally:
                if publisher is not None:
                    publisher.stop()
        except KeyboardInterrupt:
            return 0
        except Exception as exc:  # a failover is a retry, never an outage
            log("error: %s: %s" % (type(exc).__name__, exc))
            if args.once or args.verify or args.pending:
                return 1
            time.sleep(15.0)


if __name__ == "__main__":
    sys.exit(main())
