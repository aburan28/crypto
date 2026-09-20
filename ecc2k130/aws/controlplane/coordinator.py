"""The facade a worker, an ingest or an operator actually calls.

One object holding the three stores and the rules that relate them:

    lease  = coordinator.acquire(info)          # RDS: a fenced run id
    coordinator.heartbeat(lease, metrics)       # RDS + ElastiCache
    coordinator.publishObject(lease, name, data)# S3 first, then the ledger
    coordinator.setCheckpoint(lease, ...)       # fenced pointer
    coordinator.release(lease, state="idle")

Every write that a stale owner could still make takes the lease, not a slot
number, and is rejected if the slot's fence has moved.  That is the property
that makes the campaign safe against a process which was frozen past its
lease and wakes up believing it still owns a run id -- a spot instance under
memory pressure, a host that was paused, a partition that healed.  Without it
a 180-second lease is a hope.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import logging
import time
import uuid

from .cache import NullCache, cacheFromEnv
from .config import Config, FLEET_TTL_SECONDS, HEARTBEAT_SECONDS
from .db import databaseFromEnv
from .schema import migrate
from .slots import PostgresSlots
from .store import ObjectExists, storeFromEnv

logger = logging.getLogger(__name__)


class FenceError(RuntimeError):
    """A write was attempted under a lease the slot has moved past."""


@dataclass
class SlotLease:
    """One worker's right to walk one run id, until `expires`."""

    slot: int
    owner: str
    fence: int
    expires: int
    info: dict = field(default_factory=dict)

    @property
    def runId(self):
        # The client derives its seeds from run id = slot + 1; a slot is
        # never a run id and confusing them walks somebody else's seeds.
        return self.slot + 1

    def alive(self, now=None):
        return int(now if now is not None else time.time()) < self.expires


class Coordinator:
    def __init__(self, db, store, cache=None, campaign="ecc2k-130", owner=None,
                 clock=time.time):
        self.db = db
        self.store = store
        self.cache = cache if cache is not None else NullCache()
        self.campaign = campaign
        self.owner = owner or ("%s:%s" % (_hostTag(), uuid.uuid4().hex))
        self._clock = clock
        self.slots = PostgresSlots(db, campaign=campaign, cache=self.cache, clock=clock)

    @classmethod
    def fromEnv(cls, config=None, owner=None, migrateSchema=False):
        config = config or Config.fromEnv()
        db = databaseFromEnv(config)
        if migrateSchema:
            migrate(db)
        return cls(db, storeFromEnv(config), cacheFromEnv(config),
                   campaign=config.campaign, owner=owner)

    # -- leases -------------------------------------------------------------
    def acquire(self, info, newOnly=False):
        slot, fence = self.slots.claimLease(self.owner, info, newOnly=newOnly)
        lease = SlotLease(slot=slot, owner=self.owner, fence=fence,
                          expires=int(self._clock()) + self.slots.leaseSeconds, info=dict(info))
        logger.info("claimed slot %d fence %d", slot, fence)
        return lease

    def renew(self, lease, fields=None):
        """Extend a lease.  False means it is gone and the walk must stop.

        The fence goes into the WHERE clause, so a worker whose slot was
        re-claimed and released back to it under a *new* fence cannot renew
        the old lease and carry on writing under it.
        """
        ok = self.slots.heartbeat(lease.slot, lease.owner, fields or {}, fence=lease.fence)
        if ok:
            lease.expires = int(self._clock()) + self.slots.leaseSeconds
        return ok

    def release(self, lease, state="idle", extra=None):
        ok = self.slots.release(lease.slot, lease.owner, state=state, extra=extra,
                                fence=lease.fence)
        self.cache.dropWorker(lease.owner)
        return ok

    def checkFence(self, lease):
        """Raise unless the database still agrees this lease owns the slot.

        Read *inside* the caller's transaction where it matters; on its own
        it is a check, not a guarantee, which is why the fenced writes below
        put the fence in their WHERE clause as well.
        """
        item = self.slots.get(lease.slot)
        if item is None or int(item["fence"]) != int(lease.fence) or item["owner"] != lease.owner:
            raise FenceError("slot %d moved on: lease fence %d, registry %s"
                             % (lease.slot, lease.fence,
                                item and item.get("fence")))
        return True

    # -- the corpus ---------------------------------------------------------
    def publishObject(self, lease, name, data, records=None, streamId=None, offset=0,
                      metadata=None):
        """Put bytes in S3, then record them in the ledger.  Never the reverse.

        Returns the ledger row.  An object whose bytes are already there with
        the same hash is accepted as a retry: `putBytes` says so, the row is
        written idempotently, and neither the record count nor the dashboard
        double-counts.
        """
        result = self.store.putBytes(name, data, metadata=metadata)
        if records is None:
            records = len(data) // 32          # 32-byte packed DP records
        row = {
            "objectKey": result["key"], "slot": lease.slot, "streamId": streamId,
            "offset": int(offset), "bytes": result["bytes"], "records": int(records),
            "sha256": result["sha256"], "fence": lease.fence,
            "uploadedAt": int(self._clock()),
        }
        self.recordObject(lease, row)
        return row

    def recordObject(self, lease, row):
        """The ledger insert, fenced.

        The fence check is a subquery in the INSERT rather than a separate
        SELECT: a check that is not in the same statement as the write is a
        race with the next claim, and this one runs often enough to find it.
        """
        inserted = self.db.execute(
            "INSERT INTO rho_objects (campaign_id, object_key, slot, stream_id, byte_offset, "
            "bytes, records, sha256, fence, uploaded_at) "
            "SELECT ?, ?, ?, ?, ?, ?, ?, ?, ?, ? WHERE EXISTS ("
            "  SELECT 1 FROM rho_slots WHERE campaign_id = ? AND slot = ? AND owner = ? "
            "  AND fence = ?) "
            "ON CONFLICT (campaign_id, object_key) DO NOTHING",
            (self.campaign, row["objectKey"], int(row["slot"]), row.get("streamId"),
             int(row.get("offset") or 0), int(row["bytes"]), int(row["records"]),
             row.get("sha256"), int(row["fence"]), int(row["uploadedAt"]),
             self.campaign, int(lease.slot), lease.owner, int(lease.fence)))
        if inserted:
            return True
        # Zero rows is either "already recorded" (a retry, fine) or "the
        # fence moved" (a zombie writer, not fine).  Telling them apart is
        # worth one extra read on a path that is already doing network I/O.
        existing = self.db.fetchOne(
            "SELECT sha256 FROM rho_objects WHERE campaign_id = ? AND object_key = ?",
            (self.campaign, row["objectKey"]))
        if existing is not None:
            if existing[0] and row.get("sha256") and existing[0] != row["sha256"]:
                raise ObjectExists("%s already recorded with a different hash"
                                   % row["objectKey"])
            return False
        raise FenceError("slot %d fence %d is stale; object %s not recorded"
                         % (lease.slot, lease.fence, row["objectKey"]))

    def setCheckpoint(self, lease, objectKey, iters, sha256=None):
        """Move a slot's checkpoint pointer, under the fence that owns it.

        The `WHERE fence <= ?` on the update is what makes a resumed slot
        safe: a worker that lost the lease and comes back cannot pull the
        pointer back to its own older checkpoint, because the new owner's
        fence is higher.
        """
        now = int(self._clock())
        changed = self.db.execute(
            "INSERT INTO rho_checkpoints (campaign_id, slot, object_key, iters, fence, "
            "sha256, updated_at) SELECT ?, ?, ?, ?, ?, ?, ? WHERE EXISTS ("
            "  SELECT 1 FROM rho_slots WHERE campaign_id = ? AND slot = ? AND owner = ? "
            "  AND fence = ?) "
            "ON CONFLICT (campaign_id, slot) DO UPDATE SET object_key = excluded.object_key, "
            "iters = excluded.iters, fence = excluded.fence, sha256 = excluded.sha256, "
            "updated_at = excluded.updated_at WHERE rho_checkpoints.fence <= excluded.fence",
            (self.campaign, int(lease.slot), objectKey, int(iters), int(lease.fence),
             sha256, now, self.campaign, int(lease.slot), lease.owner, int(lease.fence)))
        if not changed:
            raise FenceError("checkpoint pointer for slot %d not moved: fence %d is stale"
                             % (lease.slot, lease.fence))
        return True

    def checkpoint(self, slot):
        row = self.db.fetchOne(
            "SELECT object_key, iters, fence, sha256, updated_at FROM rho_checkpoints "
            "WHERE campaign_id = ? AND slot = ?", (self.campaign, int(slot)))
        if row is None:
            return None
        return {"objectKey": row[0], "iters": row[1], "fence": row[2], "sha256": row[3],
                "updatedAt": row[4]}

    def objects(self, slot=None, since=0, limit=1000):
        if slot is None:
            rows = self.db.fetchAll(
                "SELECT object_key, slot, records, bytes, sha256, uploaded_at FROM rho_objects "
                "WHERE campaign_id = ? AND uploaded_at >= ? ORDER BY uploaded_at, object_key "
                "LIMIT ?", (self.campaign, int(since), int(limit)))
        else:
            rows = self.db.fetchAll(
                "SELECT object_key, slot, records, bytes, sha256, uploaded_at FROM rho_objects "
                "WHERE campaign_id = ? AND slot = ? AND uploaded_at >= ? "
                "ORDER BY byte_offset LIMIT ?",
                (self.campaign, int(slot), int(since), int(limit)))
        return [{"objectKey": r[0], "slot": r[1], "records": r[2], "bytes": r[3],
                 "sha256": r[4], "uploadedAt": r[5]} for r in rows]

    # -- distinguished-point hints -----------------------------------------
    def offerPoint(self, pointKey):
        """Has this point been offered in the last hour?  None if unknown.

        A hint, and nothing more.  A True is worth one indexed lookup in
        `distinguished_points`, which is where a genuine collision is
        confirmed and recorded; a False or a None says nothing at all.  The
        merge over the S3 corpus stays the independent detector either way.
        """
        return self.cache.seenPoint(pointKey)

    # -- the fleet ----------------------------------------------------------
    def heartbeat(self, lease, metrics=None):
        """Renew the lease and publish the worker's telemetry.

        The renewal is the part that matters and it goes first: if the
        database is reachable but the cache is not, the walk continues and
        the dashboard is stale, which is the right way round.
        """
        metrics = dict(metrics or {})
        alive = self.renew(lease, {k: v for k, v in metrics.items()
                                   if k in ("iters", "points", "gpuName", "instanceType")})
        now = int(self._clock())
        payload = dict(metrics, slot=lease.slot, fence=lease.fence,
                       instance=lease.info.get("instance"),
                       instanceType=lease.info.get("instanceType"),
                       gpuName=lease.info.get("gpuName"),
                       gpuFamily=lease.info.get("gpuFamily"), alive=alive)
        self.cache.beat(lease.owner, payload)
        self.db.execute(
            "INSERT INTO rho_workers (campaign_id, owner, slot, instance, instance_type, gpu, "
            "gpu_name, gpu_family, client_sha, iters, points, rate, first_seen, last_seen) "
            "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?) "
            "ON CONFLICT (campaign_id, owner) DO UPDATE SET slot = excluded.slot, "
            "iters = excluded.iters, points = excluded.points, rate = excluded.rate, "
            "gpu_name = excluded.gpu_name, gpu_family = excluded.gpu_family, "
            "instance_type = excluded.instance_type, client_sha = excluded.client_sha, "
            "last_seen = excluded.last_seen",
            (self.campaign, lease.owner, int(lease.slot), lease.info.get("instance"),
             lease.info.get("instanceType"), lease.info.get("gpu"),
             lease.info.get("gpuName"), lease.info.get("gpuFamily"),
             metrics.get("clientSha"), int(metrics.get("iters") or 0),
             int(metrics.get("points") or 0), float(metrics.get("rate") or 0.0), now, now))
        return alive

    def fleet(self, ttl=FLEET_TTL_SECONDS, fresh=False):
        """Live workers.  From the cache when it can answer, else from RDS.

        `source` in the result says which, because a dashboard that cannot
        tell a warm cache from a cold one will eventually report an empty
        fleet as an outage.
        """
        if not fresh:
            rows = self.cache.fleet()
            if rows is not None:
                return {"source": "cache", "workers": rows}
        cutoff = int(self._clock()) - ttl
        rows = self.db.fetchAll(
            "SELECT owner, slot, instance_type, gpu_name, gpu_family, iters, points, rate, "
            "last_seen FROM rho_workers WHERE campaign_id = ? AND last_seen >= ? "
            "ORDER BY slot", (self.campaign, cutoff))
        workers = [{"owner": r[0], "slot": r[1], "instanceType": r[2], "gpuName": r[3],
                    "gpuFamily": r[4], "iters": r[5], "points": r[6], "rate": r[7],
                    "lastSeen": r[8]} for r in rows]
        return {"source": "database", "workers": workers}

    def status(self):
        """One object for the dashboard, the CLI and the alarms.

        Deliberately a handful of aggregate queries rather than a scan: the
        `dp_ingest` note records what a full-table read of the points costs on
        the instance the ingest writes to, and a status page that hurts the
        campaign it reports on does not stay enabled.
        """
        now = int(self._clock())
        slots = self.slots.scan(fresh=True)
        byState = {}
        for item in slots:
            byState[item.get("state") or "idle"] = byState.get(item.get("state") or "idle", 0) + 1
        live = [it for it in slots if int(it.get("leaseUntil") or 0) >= now]
        totals = self.db.fetchOne(
            "SELECT COUNT(*), COALESCE(SUM(records), 0), COALESCE(SUM(bytes), 0), "
            "COALESCE(MAX(uploaded_at), 0) FROM rho_objects WHERE campaign_id = ?",
            (self.campaign,)) or (0, 0, 0, 0)
        fleet = self.fleet()
        return {
            "campaign": self.campaign,
            "at": now,
            "slots": {"total": len(slots), "leased": len(live), "byState": byState,
                      "expired": len(self.slots.expired(now))},
            "corpus": {"objects": totals[0], "records": totals[1], "bytes": totals[2],
                       "lastUploadAt": totals[3]},
            "fleet": {"source": fleet["source"], "live": len(fleet["workers"]),
                      "iters": sum(int(w.get("iters") or 0) for w in fleet["workers"])},
            "cache": self.cache.stats(),
        }

    # -- housekeeping -------------------------------------------------------
    def prune(self, before=None, keepEvents=None):
        """Drop worker rows and slot events older than a cutoff.

        Slots, objects and checkpoints are never pruned: they are the record
        of which seeds have been walked and which bytes exist, and a campaign
        that forgets either has to start over.
        """
        before = int(before if before is not None else self._clock() - 7 * 86400)
        keepEvents = int(keepEvents if keepEvents is not None else before)
        workers = self.db.execute(
            "DELETE FROM rho_workers WHERE campaign_id = ? AND last_seen < ?",
            (self.campaign, before))
        events = self.db.execute(
            "DELETE FROM rho_slot_events WHERE campaign_id = ? AND at < ?",
            (self.campaign, keepEvents))
        return {"workers": workers, "events": events}

    def close(self):
        self.db.close()
        self.cache.close()


def _hostTag():
    import socket

    try:
        return socket.gethostname()
    except Exception:
        return "host"


def workerInfo(gpu=0, gpuName="", gpuFamily="", instance="", instanceType=""):
    """The claimant description the registry's rules read.

    Keys match `worker.py`'s `info` exactly -- `idleSlotClaimable` and
    `slotFamilyCompatible` read `gpuName` and `gpuFamily` from it -- so that
    the two callers cannot drift.
    """
    return {"gpu": int(gpu), "gpuName": gpuName, "gpuFamily": gpuFamily,
            "instance": instance, "instanceType": instanceType}


def heartbeatInterval():
    return HEARTBEAT_SECONDS


__all__ = ["Coordinator", "SlotLease", "FenceError", "workerInfo", "heartbeatInterval"]
