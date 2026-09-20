"""Environment-driven configuration for the control plane.

Every setting has one name, read in one place.  An absent setting means an
absent component: no `RHO_REDIS_URL` is a fleet with no cache, not an error,
and the coordinator then talks to Postgres for everything it would have
cached.  That is also the control to compare a cached run against.
"""

from __future__ import annotations

from dataclasses import dataclass
import os

# A lease is short enough that a lost spot instance frees its run id in
# minutes, and long enough that a worker whose heartbeat is slowed by an S3
# upload does not lose a walk it is still making.  Kept identical to
# worker.py's LEASE_SECONDS: two different lease lengths against one registry
# is a race nobody would find.
LEASE_SECONDS = 180
HEARTBEAT_SECONDS = 60

# Redis operations use this timeout and do not retry.  A cache slower than the
# query it replaces is worse than no cache; the database is right there.
REDIS_TIMEOUT_SECONDS = 0.5
# After this many consecutive failures the cache is skipped entirely for
# REDIS_COOLDOWN_SECONDS, so a dead cluster costs one timeout per cooldown
# rather than one per call.
REDIS_FAILURES_BEFORE_COOLDOWN = 3
REDIS_COOLDOWN_SECONDS = 30

# How long a cached scan of the slot table may be served.  This bounds the
# staleness a *claim* may see, not the claim itself: the claim is a
# conditional UPDATE in Postgres, so a stale snapshot costs a retry and can
# never hand one run id to two workers.
SLOT_SNAPSHOT_TTL_SECONDS = 5
# Heartbeats are the dashboard's view of the fleet.  A worker beats every
# HEARTBEAT_SECONDS, so three missed beats drop it out of the live set.
FLEET_TTL_SECONDS = 3 * HEARTBEAT_SECONDS
# The distinguished-point hot set is a duplicate/collision *hint* for points
# seen in the last hour.  It is never read as evidence that a point is new.
DP_HINT_TTL_SECONDS = 3600


def _int(env, name, default):
    raw = (env.get(name) or "").strip()
    if not raw:
        return default
    return int(raw)


def _flag(env, name, default=True):
    raw = (env.get(name) or "").strip().lower()
    if not raw:
        return default
    return raw in ("1", "true", "yes", "on")


@dataclass(frozen=True)
class Config:
    """Resolved control-plane settings.

    `databaseUrl` is resolved lazily by `db.py` rather than here, because the
    Secrets Manager path needs boto3 and this dataclass must import on a host
    that has none of the AWS SDKs.
    """

    campaign: str = "ecc2k-130"
    bucket: str = ""
    prefix: str = ""
    localStore: str = ""
    databaseUrl: str = ""
    dbSecret: str = "rho/dp-rds"
    dbHost: str = ""
    dbName: str = ""
    dbSslMode: str = "require"
    dbConnectTimeout: int = 30
    dbStatementTimeoutMs: int = 15000
    redisUrl: str = ""
    s3Endpoint: str = ""
    s3Region: str = ""
    cacheEnabled: bool = True
    leaseSeconds: int = LEASE_SECONDS

    @classmethod
    def fromEnv(cls, env=None):
        env = os.environ if env is None else env

        def get(name, default=""):
            return (env.get(name) or default).strip()

        return cls(
            campaign=get("RHO_CAMPAIGN", "ecc2k-130"),
            bucket=get("ECC_BUCKET"),
            prefix=get("ECC_PREFIX"),
            localStore=get("ECC_LOCAL_STORE"),
            # DATABASE_URL is what dp_ingest.py and ingest.sh already use;
            # spelling it differently here would give one campaign two
            # databases the day somebody sets only one of them.
            databaseUrl=get("DATABASE_URL") or get("RHO_DB_URL"),
            dbSecret=get("RHO_DB_SECRET", "rho/dp-rds"),
            dbHost=get("RHO_DB_HOST"),
            dbName=get("RHO_DB_NAME"),
            dbSslMode=get("RHO_DB_SSLMODE", "require"),
            dbConnectTimeout=_int(env, "RHO_DB_CONNECT_TIMEOUT", 30),
            dbStatementTimeoutMs=_int(env, "RHO_DB_STATEMENT_TIMEOUT_MS", 15000),
            # rediss:// means TLS, which is how the ElastiCache cluster is
            # configured; a redis:// URL against it simply will not connect.
            redisUrl=get("RHO_REDIS_URL") or get("INDEXCALC_REDIS_URL"),
            # An S3-compatible endpoint: MinIO in the integration environment,
            # and what an on-premises object store or GCS's S3 interoperability
            # API would be pointed at.  Empty means real S3, which is what
            # every deployed host uses.
            s3Endpoint=get("ECC_S3_ENDPOINT"),
            s3Region=get("ECC_S3_REGION") or get("AWS_REGION") or get("AWS_DEFAULT_REGION"),
            cacheEnabled=_flag(env, "RHO_CACHE", True),
            leaseSeconds=_int(env, "RHO_LEASE_SECONDS", LEASE_SECONDS),
        )

    def objectKey(self, *parts):
        """A bucket key under the configured prefix."""
        path = "/".join(str(p).strip("/") for p in parts if str(p).strip("/"))
        return "%s/%s" % (self.prefix.strip("/"), path) if self.prefix.strip("/") else path
