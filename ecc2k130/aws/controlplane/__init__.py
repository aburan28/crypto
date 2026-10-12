"""Control plane for the distributed rho campaign: RDS, ElastiCache and S3.

The campaign already had three stores and no single place that said how they
relate.  Slots were leased from DynamoDB (or from S3 objects), the dashboard's
derived view lived in RDS Postgres, the corpus lived in S3, and the rules that
keep them consistent -- who owns a run id, which writes a stale owner may still
make, what a Redis miss means -- were spread across `worker.py`, `dp_ingest.py`
and a shell script.  This package is that single place.

    authority       RDS Postgres     slot leases, fencing tokens, the object
                                     ledger and the checkpoint pointers
    cache           ElastiCache      scan snapshots, fleet heartbeats, the
                                     distinguished-point hot set, fleet-wide
                                     locks -- every one of them optional
    corpus          S3               the bytes; immutable, content-addressed

Three invariants hold everything together, and they are the part worth keeping
if the rest is rewritten:

  * **Postgres is the authority and Redis is a hint.**  Nothing here reads a
    Redis miss as evidence, and every cache failure degrades to the database.
    A cache that can change an answer is not a cache.
  * **S3 holds the bytes; the database holds what is true about them.**  A
    ledger row without its object is a lie the merge would believe, so the
    object is uploaded first and the row committed after, never the reverse.
  * **A lease is fenced.**  Every claim bumps a per-slot fence token, and every
    write that a stale owner could still make -- an object row, a checkpoint
    pointer -- carries the fence it was issued under and is rejected if the
    slot has moved on.  A 180-second lease cannot stop a paused process from
    waking up and writing; a fence can.

`psycopg`, `redis` and `boto3` are imported lazily inside the component that
needs them, so this package imports -- and its tests run -- on a machine that
has none of the three.

Under AGENTS.md §3 this is plumbing, not a measured change: it makes no
performance claim and moves no row on any scoreboard.  What it is allowed to
claim is the failure modes it removes, and those are listed in README.md.
"""

from .config import Config
from .db import Database, DatabaseError
from .cache import Cache, NullCache
from .store import ObjectStore, LocalObjectStore, ObjectExists
from .slots import PostgresSlots, LEASE_SECONDS
from .coordinator import Coordinator, SlotLease, FenceError

__all__ = [
    "Config", "Database", "DatabaseError", "Cache", "NullCache",
    "ObjectStore", "LocalObjectStore", "ObjectExists", "PostgresSlots",
    "LEASE_SECONDS", "Coordinator", "SlotLease", "FenceError",
]
