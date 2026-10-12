"""The control-plane schema and its migrations.

Five tables, all keyed by `campaign_id` so one database can carry a running
campaign and a rehearsal without either being able to see the other's run ids:

    rho_slots         one row per run id: owner, lease, state, fence
    rho_slot_events   an append-only log of claims, releases and retirements
    rho_workers       one row per live supervisor, for the dashboard
    rho_objects       the ledger: one row per corpus object accepted into S3
    rho_checkpoints   the fenced pointer to each slot's newest checkpoint

Timestamps are `bigint` epoch seconds rather than `timestamptz`, for one
reason: a lease comparison must be made by the database, against a value the
caller supplies, so that a worker whose clock is wrong loses its lease instead
of holding one forever.  Epoch seconds also happen to make the schema run on
SQLite, which is what lets the suite test the claim logic for real.

The existing dashboard tables (`distinguished_points`, `dp_ingest_*`,
`rho_collisions`) are *not* touched here.  They are a derived view written by
`dp_ingest.py` and they predate this package; migrating them would put two
writers on one table for no gain.
"""

from __future__ import annotations

import time

SCHEMA_VERSION = 1

# Each migration is (id, [statement, ...]).  Migrations are append-only and
# are never edited after they land: a statement that has run on the live
# database cannot be un-run by changing the text here, it can only be
# contradicted by it.  Add a new one instead.
MIGRATIONS = [
    ("0001-slots", [
        """
        CREATE TABLE IF NOT EXISTS rho_slots (
            campaign_id   text    NOT NULL,
            slot          integer NOT NULL,
            owner         text,
            lease_until   bigint  NOT NULL DEFAULT 0,
            fence         bigint  NOT NULL DEFAULT 0,
            state         text    NOT NULL DEFAULT 'idle',
            instance      text,
            gpu           integer,
            gpu_name      text,
            gpu_family    text,
            instance_type text,
            iters         bigint  NOT NULL DEFAULT 0,
            points        bigint  NOT NULL DEFAULT 0,
            reason        text,
            created_at    bigint  NOT NULL DEFAULT 0,
            claimed_at    bigint  NOT NULL DEFAULT 0,
            updated_at    bigint  NOT NULL DEFAULT 0,
            PRIMARY KEY (campaign_id, slot)
        )
        """,
        # The claim path asks for "slots whose lease has expired, oldest
        # first"; without this it reads the whole registry on every claim,
        # and every spot interruption claims at once.
        "CREATE INDEX IF NOT EXISTS rho_slots_lease ON rho_slots (campaign_id, lease_until)",
        """
        CREATE TABLE IF NOT EXISTS rho_slot_events (
            campaign_id text    NOT NULL,
            slot        integer NOT NULL,
            at          bigint  NOT NULL,
            kind        text    NOT NULL,
            owner       text,
            fence       bigint,
            detail      text
        )
        """,
        "CREATE INDEX IF NOT EXISTS rho_slot_events_at ON rho_slot_events (campaign_id, at)",
        "CREATE INDEX IF NOT EXISTS rho_slot_events_slot ON rho_slot_events (campaign_id, slot, at)",
    ]),
    ("0002-workers", [
        """
        CREATE TABLE IF NOT EXISTS rho_workers (
            campaign_id   text NOT NULL,
            owner         text NOT NULL,
            slot          integer,
            instance      text,
            instance_type text,
            gpu           integer,
            gpu_name      text,
            gpu_family    text,
            client_sha    text,
            iters         bigint NOT NULL DEFAULT 0,
            points        bigint NOT NULL DEFAULT 0,
            rate          double precision NOT NULL DEFAULT 0,
            first_seen    bigint NOT NULL DEFAULT 0,
            last_seen     bigint NOT NULL DEFAULT 0,
            PRIMARY KEY (campaign_id, owner)
        )
        """,
        "CREATE INDEX IF NOT EXISTS rho_workers_seen ON rho_workers (campaign_id, last_seen)",
    ]),
    ("0003-objects", [
        # The ledger.  `object_key` is the primary key because S3 keys are
        # content-addressed under the strict protocol (`protocol.py`), so an
        # object re-uploaded after a retry is the same key and the same row:
        # re-ingest is a no-op rather than a double count.
        """
        CREATE TABLE IF NOT EXISTS rho_objects (
            campaign_id text    NOT NULL,
            object_key  text    NOT NULL,
            slot        integer NOT NULL,
            stream_id   text,
            byte_offset bigint  NOT NULL DEFAULT 0,
            bytes       bigint  NOT NULL DEFAULT 0,
            records     bigint  NOT NULL DEFAULT 0,
            sha256      text,
            fence       bigint  NOT NULL DEFAULT 0,
            uploaded_at bigint  NOT NULL DEFAULT 0,
            PRIMARY KEY (campaign_id, object_key)
        )
        """,
        "CREATE INDEX IF NOT EXISTS rho_objects_slot ON rho_objects (campaign_id, slot, byte_offset)",
        "CREATE INDEX IF NOT EXISTS rho_objects_uploaded ON rho_objects (campaign_id, uploaded_at)",
        # One row per slot: where its newest checkpoint is, under which fence.
        # A worker that lost its lease cannot move this pointer, which is the
        # whole reason it is a table and not an S3 object with a timestamp.
        """
        CREATE TABLE IF NOT EXISTS rho_checkpoints (
            campaign_id text    NOT NULL,
            slot        integer NOT NULL,
            object_key  text    NOT NULL,
            iters       bigint  NOT NULL DEFAULT 0,
            fence       bigint  NOT NULL DEFAULT 0,
            sha256      text,
            updated_at  bigint  NOT NULL DEFAULT 0,
            PRIMARY KEY (campaign_id, slot)
        )
        """,
    ]),
]

MIGRATION_TABLE = """
CREATE TABLE IF NOT EXISTS rho_controlplane_migrations (
    id         text   PRIMARY KEY,
    version    integer NOT NULL,
    applied_at bigint  NOT NULL
)
"""

# One arbitrary but fixed key, so that two hosts running `migrate` at the same
# moment -- a fleet rollout does exactly that -- serialise instead of racing
# each other through `CREATE INDEX`.
ADVISORY_LOCK_KEY = 0x5209_3130


def _lostTheCreateRace(exc):
    """True when `CREATE TABLE IF NOT EXISTS` lost a race with another host.

    Postgres checks whether the table exists and inserts into the catalog in
    two steps that are not atomic, so two hosts creating the same table at the
    same moment leave one of them holding a unique violation on a `pg_catalog`
    index instead of the no-op it asked for.  The table is there either way,
    which is all the caller wanted.  Only ever applied to the migration
    table's own statement, where there is no other unique constraint to
    confuse this with.
    """
    if type(exc).__name__ in ("DuplicateTable", "DuplicateObject", "UniqueViolation"):
        return True
    return "already exists" in str(exc).lower()


def _ensureMigrationTable(db):
    try:
        db.run(MIGRATION_TABLE)
    except Exception as exc:
        if not _lostTheCreateRace(exc):
            raise


def migrate(db, now=None, log=None):
    """Apply every migration the database has not seen.  Idempotent.

    Returns the list of ids applied by this call, which is empty on a database
    that is already current -- so a deploy that calls this on every boot is
    free after the first one.
    """
    now = int(now if now is not None else time.time())
    applied = []
    locked = False
    if db.dialect == "postgres":
        # Taken before the migration table is created, not after: a rollout
        # runs this on every host at once, and `CREATE TABLE IF NOT EXISTS`
        # is not itself safe against a concurrent one (see
        # `_lostTheCreateRace`).  Everything that touches the schema has to be
        # inside the lock for the lock to mean anything.
        #
        # Session-level, released by the unlock below or by the connection
        # dying, which is the behaviour wanted if the migrating host is the
        # thing that fails.
        db.run("SELECT pg_advisory_lock(?)", (ADVISORY_LOCK_KEY,), fetch="one")
        locked = True
    try:
        _ensureMigrationTable(db)
        done = {row[0] for row in db.fetchAll("SELECT id FROM rho_controlplane_migrations")}
        for ident, statements in MIGRATIONS:
            if ident in done:
                continue
            for statement in statements:
                # `CREATE INDEX` and `CREATE TABLE` are each their own
                # statement rather than one batch: a partially applied
                # migration re-runs harmlessly because every statement is
                # `IF NOT EXISTS`, and the marker row is only written once
                # the whole list has succeeded.
                db.run(statement)
            db.run("INSERT INTO rho_controlplane_migrations (id, version, applied_at) "
                   "VALUES (?, ?, ?) ON CONFLICT (id) DO NOTHING",
                   (ident, SCHEMA_VERSION, now))
            applied.append(ident)
            if log:
                log("applied migration %s" % ident)
    finally:
        if locked:
            try:
                db.run("SELECT pg_advisory_unlock(?)", (ADVISORY_LOCK_KEY,), fetch="one")
            except Exception:
                pass
    return applied


def pending(db):
    """Migration ids this database has not applied."""
    # Not under the advisory lock: this is a read, and `doctor` calls it on
    # every host.  It still has to survive the create race above, because a
    # rollout is exactly when somebody runs it.
    _ensureMigrationTable(db)
    done = {row[0] for row in db.fetchAll("SELECT id FROM rho_controlplane_migrations")}
    return [ident for ident, _ in MIGRATIONS if ident not in done]
