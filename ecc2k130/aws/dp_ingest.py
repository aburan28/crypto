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

`found_at` is taken from the object's upload time rather than the wall clock at
insert. The old watcher used insert time, which is the same thing while it
keeps pace and a lie when it does not: catching up on a backlog would stack a
million points into the hour the catch-up ran and leave a spike in the
published hourly chart that no GPU ever produced. The upload time is also
stable across re-ingest, so a repeated object cannot move a bucket.

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
    return "postgresql://%s:%s@%s:5432/%s" % (
        c["username"], urllib.request.quote(c["password"], safe=""), host, c["dbname"])


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


def ensureProgress(conn):
    with conn.cursor() as cur:
        cur.execute(PROGRESS_DDL)
    conn.commit()


def ingestedCounts(conn):
    """Rows already present per worker_id, and objects recorded as complete.

    Two sources because the corpus has two eras.  Objects this program ingested
    have a `dp_ingest_progress` row and that is authoritative.  Objects the
    previous ingester wrote have no such row, so they are recognised by their row
    count under the matching worker_id.

    The count alone is not sufficient, which is the whole reason the table
    exists: a point already stored under some *other* worker_id -- an earlier
    era attributed points to `gpu-i-<instance>-slot-N` rather than the object
    name -- is skipped by `ON CONFLICT DO NOTHING`, so the count under this
    object's own worker_id never reaches its record count and the object would
    be re-ingested on every pass forever.  Harmless for correctness, but it
    means the backlog never converges.
    """
    with conn.cursor() as cur:
        cur.execute(
            "SELECT worker_id, count(*) FROM distinguished_points "
            "WHERE campaign_id = %s AND worker_id LIKE 'dp-slot-%%' GROUP BY worker_id",
            (CAMPAIGN,))
        counts = {row[0]: int(row[1]) for row in cur.fetchall()}
        cur.execute(
            "SELECT object_key, records FROM dp_ingest_progress WHERE campaign_id = %s",
            (CAMPAIGN,))
        done = {row[0]: int(row[1]) for row in cur.fetchall()}
    return counts, done


def ingestObject(conn, s3, bucket, key, found_at):
    """Insert every record of one object. Returns rows actually added."""
    body = s3.get_object(Bucket=bucket, Key=key)["Body"].read()
    whole = len(body) - len(body) % RECORD_BYTES
    wid = workerId(key)
    added = 0
    with conn.cursor() as cur:
        # A temp table plus one INSERT ... ON CONFLICT is both fast and
        # idempotent; row-at-a-time INSERT was the old ingester's ceiling.
        cur.execute("CREATE TEMP TABLE IF NOT EXISTS dp_in "
                    "(point_key bytea, a bytea, b bytea, walk_seed bytea) ON COMMIT DROP")
        cur.execute("TRUNCATE dp_in")
        # Sorted by the key the target is indexed on. point_key is a hash, so
        # an object's records arrive in random index order and each insert
        # walks to a different leaf page: on the 60 M-row table that measured
        # ~2.2 random reads per row inserted. Sorting costs nothing here and
        # turns a batch into a smaller set of pages touched repeatedly.
        rows = sorted((decode(body[off:off + RECORD_BYTES])
                       for off in range(0, whole, RECORD_BYTES)),
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
        added = cur.rowcount
        # Same transaction as the points, so the record of having ingested an
        # object cannot outlive the insert that it describes.
        cur.execute(
            "INSERT INTO dp_ingest_progress (campaign_id, object_key, records) "
            "VALUES (%s, %s, %s) ON CONFLICT (campaign_id, object_key) "
            "DO UPDATE SET records = EXCLUDED.records, ingested_at = now()",
            (CAMPAIGN, key, whole // RECORD_BYTES))
    conn.commit()
    return added, whole // RECORD_BYTES


def verify(conn, s3, bucket, sample=64):
    """Check this program's encodings against points the old ingester wrote.

    Not a replication: it decodes the object independently and asks whether the
    resulting point_key is already stored under the same worker_id. Agreement
    means a resumed ingest lands in the same key space, which is the only
    property collision detection depends on.
    """
    counts, _ = ingestedCounts(conn)
    objects, _ = s3Objects(s3, bucket)
    complete = [(k, n) for k, n, _ in objects if counts.get(workerId(k), 0) == n]
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
CKPT_RE = re.compile(r"^ckpt/(retired/)?slot-(\d+)\.ck$")


def checkpointWork(s3, bucket, staleAfter=1800.0):
    """Exact iterations walked, read from the slot checkpoints.

    This is the honest measure of work done and it needs no assumption about
    distinguished-point density.  Deriving work from the stored point count
    instead requires a constant iterations-per-point, and that constant changes
    whenever `dpWeight` changes -- so a campaign that has been retuned (as this
    one was, 34 -> 32 on 2026-09-13) has no single correct multiplier, and using
    the old one silently under-reports.

    The client prints "N iterations of M walks": `iterBase` is per-walk steps,
    not a total, so the work of one slot is `iterBase * threads * batch`.  Every
    factor comes out of the checkpoint header, so a slot is self-describing and
    a geometry change needs no bookkeeping here.  Retired slots are included --
    their work is part of the campaign whether or not they still run.
    """
    slots, token = [], None
    while True:
        kw = {"Bucket": bucket, "Prefix": "ckpt/"}
        if token:
            kw["ContinuationToken"] = token
        page = s3.list_objects_v2(**kw)
        for item in page.get("Contents", []):
            m = CKPT_RE.match(item["Key"])
            if not m:
                continue
            head = s3.get_object(Bucket=bucket, Key=item["Key"],
                                 Range="bytes=0-%d" % (CKPT_HEADER - 1))["Body"].read()
            if len(head) < CKPT_HEADER or head[:8] != CKPT_MAGIC:
                log("checkpoint %s: not a checkpoint header, skipped" % item["Key"])
                continue
            version, m131, threads, batch, lanes, runId = struct.unpack_from("<6I", head, 8)
            iterBase, = struct.unpack_from("<Q", head, 32)
            walks = threads * batch
            age = time.time() - item["LastModified"].timestamp()
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
        if not page.get("IsTruncated"):
            break
        token = page["NextContinuationToken"]
    slots.sort(key=lambda s: (s["retired"], s["slot"]))
    total = sum(s["iterations"] for s in slots)
    return total, slots


def statusPayload(conn, s3, bucket, ingest=None):
    """The dashboard snapshot, computed here so the page can read it directly.

    `ingest` is what the last pass knew about itself. It is published because
    `state` alone cannot tell the difference between a fleet that stopped and
    an ingest that stopped: on 2026-09-17 the page read IDLE_OR_STALE for five
    hours while 133 workers were walking and 56 M records were landing in S3.
    """
    with conn.cursor() as cur:
        cur.execute("""
            SELECT c.campaign_id, c.curve_id, c.dp_mask_bits, c.created_at,
                   count(d.point_key) AS dps,
                   count(d.point_key) FILTER (WHERE d.found_at > now() - interval '1 hour') AS dps_last_hour,
                   count(d.point_key) FILTER (WHERE d.found_at > now() - interval '1 day') AS dps_last_day,
                   min(d.found_at) AS first_dp_at, max(d.found_at) AS last_dp_at
            FROM rho_campaigns c LEFT JOIN distinguished_points d USING (campaign_id)
            WHERE c.campaign_id = %s GROUP BY 1,2,3,4""", (CAMPAIGN,))
        row = cur.fetchone() or ()
        cur.execute("SELECT count(*), max(detected_at) FROM rho_collisions WHERE campaign_id = %s",
                    (CAMPAIGN,))
        collisions, latest_collision = cur.fetchone()
        cur.execute("""
            SELECT date_trunc('hour', found_at) AS hour, count(*)
            FROM distinguished_points
            WHERE campaign_id = %s AND found_at > now() - interval '48 hours'
            GROUP BY 1 ORDER BY 1""", (CAMPAIGN,))
        hourly = cur.fetchall()

    def iso(v):
        return v.isoformat() if hasattr(v, "isoformat") else v

    dps = int(row[4] or 0) if row else 0
    dps_last_hour = int(row[5] or 0) if row else 0
    collisions = int(collisions or 0)
    ingest = dict(ingest or {})
    # An ingest that is behind makes every recency figure below a statement
    # about this program, not about the fleet, so it is named as one rather
    # than left to be read as a quiet campaign.
    behind = int(ingest.get("outstanding", 0)) or int(ingest.get("unrecognised", 0))
    state = ("COLLISION_RECORDED" if collisions else
             "COLLECTING" if dps_last_hour else
             "INGEST_BEHIND" if behind else
             "IDLE_OR_STALE" if dps else "EMPTY")
    iterations, per_slot = checkpointWork(s3, bucket)
    payload = {
        "schema_version": 1,
        "generated_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "source": "dp_ingest.py (live, direct from rho-dp + slot checkpoints)",
        "campaign_id": CAMPAIGN,
        "curve_id": row[1] if row else None,
        "dp_mask_bits": row[2] if row else None,
        "campaign_created_at": iso(row[3]) if row else None,
        "state": state,
        "dps": dps,
        "collisions": collisions,
        "first_dp_at": iso(row[7]) if row else None,
        "last_dp_at": iso(row[8]) if row else None,
        "latest_collision_at": iso(latest_collision),
        "dps_last_hour": dps_last_hour,
        "dps_last_day": int(row[6] or 0) if row else 0,
        "hourly": [{"hour": iso(h), "dps": int(n)} for h, n in hourly],
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


def publishStatus(conn, s3, bucket, statusBucket, ingest=None):
    """Write the snapshot where a browser can read it, no workflow involved."""
    payload = statusPayload(conn, s3, bucket, ingest)
    s3.put_object(
        Bucket=statusBucket, Key="status.json",
        Body=(json.dumps(payload, indent=2, sort_keys=True) + "\n").encode(),
        ContentType="application/json",
        # Short but non-zero: the underlying data only moves when a worker
        # uploads, so caching for less than that buys nothing and costs requests.
        CacheControl="public, max-age=30")
    log("published status.json: dps=%d state=%s work=2^%.3f walkers=%d "
        "outstanding=%d unreadable=%d"
        % (payload["dps"], payload["state"],
           payload["work"]["iterations_log2"] or 0, payload["walkers"],
           payload["ingest"]["outstanding_objects"],
           payload["ingest"]["unrecognised_objects"]))
    return payload


def pending(conn, s3, bucket, prefix="dp/"):
    """Objects in the bucket that the store does not have, oldest first.

    Also returns the newest upload time seen and any key it could not read, so
    that being behind and being unable to read are two different, reported
    numbers rather than one silence.
    """
    counts, done = ingestedCounts(conn)
    objects, unrecognised = s3Objects(s3, bucket, prefix)
    todo, newest = [], 0
    for key, records, when in objects:
        newest = max(newest, when)
        have = counts.get(workerId(key), 0)
        if done.get(key, -1) >= records or have >= records:
            continue
        todo.append((key, records, when, have))
    if unrecognised:
        log("WARNING: %d object(s) under dp/ match no known key shape and are "
            "NOT being ingested, for example: %s"
            % (len(unrecognised), ", ".join(unrecognised[:3])))
    return todo, newest, unrecognised


def onePass(connect, s3, bucket, threads=6, prefix="dp/"):
    """Ingest everything outstanding, several objects at a time.

    Each thread owns a connection: psycopg connections are not shared, and the
    temp table the COPY lands in is per-session anyway. Order does not matter
    because every insert is idempotent, so a failed object is simply retried
    on the next pass instead of stopping the ones behind it.
    """
    with connect() as probe:
        todo, newest, unrecognised = pending(probe, s3, bucket, prefix)
    state = {"newest": newest, "unrecognised": len(unrecognised), "outstanding": len(todo)}
    if not todo:
        return 0, 0, state
    work = queue.Queue()
    for item in todo:
        work.put(item)
    tally = {"rows": 0, "objects": 0, "failed": 0}
    lock = threading.Lock()

    def drain():
        try:
            with connect() as conn:
                while True:
                    try:
                        key, records, when, have = work.get_nowait()
                    except queue.Empty:
                        return
                    try:
                        n, seen = ingestObject(
                            conn, s3, bucket, key,
                            time.strftime("%Y-%m-%d %H:%M:%S+00", time.gmtime(when)))
                    except Exception as exc:
                        # One bad object must not take the pass down with it;
                        # it keeps its place in the backlog and is retried.
                        with lock:
                            tally["failed"] += 1
                        log("object %s failed (will retry): %s: %s"
                            % (key, type(exc).__name__, exc))
                        continue
                    with lock:
                        tally["rows"] += n
                        tally["objects"] += 1
                    log("ingested %s: %d records, %d new (store had %d)"
                        % (key, seen, n, have))
        except Exception as exc:  # a connection that will not open at all
            log("ingest thread stopped: %s: %s" % (type(exc).__name__, exc))

    pool = [threading.Thread(target=drain, daemon=True) for _ in range(max(1, threads))]
    started = time.time()
    for t in pool:
        t.start()
    for t in pool:
        t.join()
    elapsed = max(time.time() - started, 1e-9)
    state["outstanding"] = len(todo) - tally["objects"]
    log("pass complete: %d objects touched, %d rows added, %d failed, %d still "
        "outstanding, %.0f rows/s"
        % (tally["objects"], tally["rows"], tally["failed"],
           state["outstanding"], tally["rows"] / elapsed))
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
        for s in per_slot:
            log("slot %05d run %d %-8s ckpt %5ds old  %d walks x %d steps = %.6g iterations"
                % (s["slot"], s["run_id"], "retired" if s["retired"] else "walking",
                   s["checkpoint_age_s"], s["walks"], s["per_walk_steps"], s["iterations"]))
        log("exact total: %d iterations = 2^%.3f (no density assumption)"
            % (total, math.log2(total) if total else 0))
        return 0

    import psycopg

    url = databaseUrl()

    def connect():
        return psycopg.connect(url, connect_timeout=30)

    while True:
        try:
            with connect() as conn:
                ensureProgress(conn)
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
            while True:
                _, _, ingest = onePass(connect, s3, args.bucket, args.threads, args.prefix)
                if args.metric_namespace:
                    publishMetrics(args.metric_namespace, ingest)
                if args.status_bucket:
                    try:
                        with connect() as conn:
                            publishStatus(conn, s3, args.bucket, args.status_bucket, ingest)
                    except Exception as exc:
                        # Publishing is a view; never let it stop the ingest.
                        log("status publish failed: %s: %s" % (type(exc).__name__, exc))
                if args.once:
                    return 0
                # A backlog is drained as fast as the database allows rather
                # than one pass per interval: the interval exists to keep an
                # idle ingest cheap, not to rate-limit catching up.
                time.sleep(0.0 if ingest.get("outstanding") else args.interval)
        except KeyboardInterrupt:
            return 0
        except Exception as exc:  # a failover is a retry, never an outage
            log("error: %s: %s" % (type(exc).__name__, exc))
            if args.once or args.verify or args.pending:
                return 1
            time.sleep(15.0)


if __name__ == "__main__":
    sys.exit(main())
