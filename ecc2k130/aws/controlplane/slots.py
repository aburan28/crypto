"""Slot leases in RDS Postgres, with the same contract as `worker.py`'s.

`worker.py` already carries three registries -- DynamoDB, S3 objects and a
locked JSON file -- behind one four-method contract: `scan`, `claim`,
`heartbeat`, `release`.  This is a fourth, on the database the campaign
already runs for the dashboard, and it exists for three things the others
cannot do:

  * **A claim is one statement.**  DynamoDB's claim is a scan followed by a
    conditional update per candidate, and the S3 one is a read-modify-write
    per candidate with an `If-Match`.  Here it is a single conditional
    `UPDATE ... RETURNING`, so the database decides, and a hundred workers
    claiming at once after a spot interruption cost one round trip each
    instead of one per candidate they lose.
  * **A lease has a fence.**  Every successful claim increments the slot's
    fence.  A worker that was paused past its lease -- a stop-the-world GC, a
    hypervisor freeze, a network partition that healed -- comes back holding a
    fence the database has moved past, and every write it attempts is
    rejected.  No lease length can prevent that; only a fence can.
  * **One store fewer.**  The slot registry, the object ledger and the
    dashboard's view are in one database, so "which worker wrote this object,
    under which lease" is a join rather than a correlation by eye across
    DynamoDB and S3.

The claimability rules -- which slots a CPU worker may resume, which GPU
family may take which checkpoint -- are *imported from `worker.py`*, never
restated.  Two copies of the rule that keeps a packed checkpoint away from
the host client is the kind of divergence that retires a run id and loses its
seeds; there is exactly one copy and it lives where the client does.
"""

from __future__ import annotations

import json
import os
import sys
import time

from .config import LEASE_SECONDS

# Run id = slot + 1 must fit in 16 bits; identical to worker.py's MAX_SLOT.
MAX_SLOT = 65534

# States a slot may be claimed from.  'retired', 'solved' and 'error' keep
# their run id forever so that no seed is ever walked twice.
CLAIMABLE_STATES = ("active", "idle")

_worker = None


def workerRules():
    """`worker.py`, imported from the directory above this package.

    Imported lazily and by path because this package is a subdirectory of the
    script it imports from -- `aws/worker.py` is a program, not an installed
    module.  A control-plane host that does not carry the supervisor gets a
    clear ImportError here rather than a subtly different claim rule.
    """
    global _worker
    if _worker is None:
        here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        if here not in sys.path:
            sys.path.insert(0, here)
        import worker  # noqa: E402  (path-dependent by design)

        _worker = worker
    return _worker


SLOT_COLUMNS = ("slot", "owner", "lease_until", "fence", "state", "instance", "gpu",
                "gpu_name", "gpu_family", "instance_type", "iters", "points", "reason",
                "created_at", "claimed_at", "updated_at")

# The contract's records are the camelCase dicts worker.py's other registries
# return; the table is snake_case like the rest of the schema.  One mapping,
# in one place.
_FIELD_TO_COLUMN = {
    "owner": "owner", "leaseUntil": "lease_until", "fence": "fence", "state": "state",
    "instance": "instance", "gpu": "gpu", "gpuName": "gpu_name", "gpuFamily": "gpu_family",
    "instanceType": "instance_type", "iters": "iters", "points": "points",
    "reason": "reason", "createdAt": "created_at", "claimedAt": "claimed_at",
    "updatedAt": "updated_at",
}
_COLUMN_TO_FIELD = {v: k for k, v in _FIELD_TO_COLUMN.items()}
_COLUMN_TO_FIELD["slot"] = "slot"


def rowToItem(row):
    item = {}
    for column, value in zip(SLOT_COLUMNS, row):
        item[_COLUMN_TO_FIELD[column]] = value
    return item


class PostgresSlots:
    """The slot registry.  Same four methods as `worker.DynamoSlots`.

    `cache` is optional and is only ever used to serve `scan()` a few seconds
    stale.  Nothing about a claim, a heartbeat or a release consults it: a
    cache that could decide who owns a run id would be a second authority.
    """

    def __init__(self, db, campaign="ecc2k-130", cache=None,
                 leaseSeconds=LEASE_SECONDS, clock=time.time):
        self.db = db
        self.campaign = campaign
        self.cache = cache
        self.leaseSeconds = leaseSeconds
        self._clock = clock

    # -- reads --------------------------------------------------------------
    def scan(self, fresh=False):
        """Every slot in the registry.

        Served from the cache when one is configured and warm, because the
        claim path calls this first and a spot interruption calls it from a
        hundred workers within a few seconds.  Staleness costs a lost claim
        attempt and nothing else.
        """
        if self.cache is not None and not fresh:
            rows = self.cache.slotSnapshot()
            if rows is not None:
                return rows
        rows = self.db.fetchAll(
            "SELECT %s FROM rho_slots WHERE campaign_id = ? ORDER BY slot"
            % ", ".join(SLOT_COLUMNS), (self.campaign,))
        items = [rowToItem(row) for row in rows]
        if self.cache is not None:
            self.cache.putSlotSnapshot(items)
        return items

    def get(self, slot):
        row = self.db.fetchOne(
            "SELECT %s FROM rho_slots WHERE campaign_id = ? AND slot = ?"
            % ", ".join(SLOT_COLUMNS), (self.campaign, int(slot)))
        return rowToItem(row) if row else None

    def fence(self, slot):
        item = self.get(slot)
        return int(item["fence"]) if item else 0

    # -- writes -------------------------------------------------------------
    def _event(self, tx, slot, kind, owner, fence, detail=None):
        tx.execute(
            "INSERT INTO rho_slot_events (campaign_id, slot, at, kind, owner, fence, detail) "
            "VALUES (?, ?, ?, ?, ?, ?, ?)",
            (self.campaign, int(slot), int(self._clock()), kind, owner, int(fence or 0),
             json.dumps(detail, sort_keys=True) if detail else None))

    def takeSlot(self, slot, owner, info):
        """Try to take one existing slot.  Returns its new fence, or None.

        The condition is the whole design: only an expired lease may be taken,
        and only from a slot still in a claimable state.  Two workers running
        this concurrently both issue the same UPDATE; Postgres serialises
        them, the second one's `lease_until < ?` no longer holds, and it gets
        zero rows back.  There is no window between the read and the write
        because there is no read.
        """
        now = int(self._clock())
        placeholders = ", ".join("?" for _ in CLAIMABLE_STATES)
        with self.db.transaction() as tx:
            row = tx.fetchOne(
                "UPDATE rho_slots SET owner = ?, lease_until = ?, claimed_at = ?, "
                "updated_at = ?, state = 'active', fence = fence + 1, instance = ?, gpu = ?, "
                "gpu_name = ?, gpu_family = ?, instance_type = ?, reason = NULL "
                "WHERE campaign_id = ? AND slot = ? AND lease_until < ? "
                "AND state IN (%s) RETURNING fence" % placeholders,
                (owner, now + self.leaseSeconds, now, now,
                 info.get("instance"), info.get("gpu"), info.get("gpuName"),
                 info.get("gpuFamily") or "", info.get("instanceType") or "",
                 self.campaign, int(slot), now) + CLAIMABLE_STATES)
            if row is None:
                return None
            fence = int(row[0])
            self._event(tx, slot, "claim", owner, fence, {"resumed": True})
        if self.cache is not None:
            self.cache.invalidateSlots()
        return fence

    def createSlot(self, slot, owner, info):
        """Try to create a new slot.  Returns its fence (1), or None if taken.

        `ON CONFLICT DO NOTHING` is the whole race: the run id is the primary
        key, so two workers reaching for the same new slot cannot both get it,
        and the loser simply tries the next one.
        """
        now = int(self._clock())
        with self.db.transaction() as tx:
            row = tx.fetchOne(
                "INSERT INTO rho_slots (campaign_id, slot, owner, lease_until, fence, state, "
                "instance, gpu, gpu_name, gpu_family, instance_type, created_at, claimed_at, "
                "updated_at) VALUES (?, ?, ?, ?, 1, 'active', ?, ?, ?, ?, ?, ?, ?, ?) "
                "ON CONFLICT (campaign_id, slot) DO NOTHING RETURNING fence",
                (self.campaign, int(slot), owner, now + self.leaseSeconds,
                 info.get("instance"), info.get("gpu"), info.get("gpuName"),
                 info.get("gpuFamily") or "", info.get("instanceType") or "", now, now, now))
            if row is None:
                return None
            fence = int(row[0])
            self._event(tx, slot, "create", owner, fence, {"gpuName": info.get("gpuName")})
        if self.cache is not None:
            self.cache.invalidateSlots()
        return fence

    def claim(self, owner, info, newOnly=False):
        """Claim a slot, resuming an expired one if the rules allow.

        Returns the slot number, and stashes the fence issued with it on
        `self.lastFence` -- the four-method contract `worker.py` uses returns
        a bare slot, and changing it would mean changing three other
        registries.  `claimLease` below returns both, and is what the
        coordinator calls.
        """
        slot, _ = self.claimLease(owner, info, newOnly=newOnly)
        return slot

    def claimLease(self, owner, info, newOnly=False):
        """`(slot, fence)`.  Raises when the run ids are exhausted."""
        rules = workerRules()
        now = int(self._clock())
        items = self.scan(fresh=True)
        free = [it for it in items
                if it.get("state") in CLAIMABLE_STATES
                and int(it.get("leaseUntil") or 0) < now
                and rules.idleSlotClaimable(it, info, newOnly)
                and rules.slotFamilyCompatible(it.get("gpuFamily"), info.get("gpuFamily"))]
        for it in sorted(free, key=lambda x: int(x["slot"])):
            fence = self.takeSlot(int(it["slot"]), owner, info)
            if fence is not None:
                self.lastFence = fence
                return int(it["slot"]), fence
        nextSlot = max((int(it["slot"]) for it in items), default=-1) + 1
        for _ in range(64):
            if nextSlot > MAX_SLOT:
                raise RuntimeError("run ids exhausted")
            fence = self.createSlot(nextSlot, owner, info)
            if fence is not None:
                self.lastFence = fence
                return nextSlot, fence
            nextSlot += 1
        raise RuntimeError("could not allocate a slot")

    def heartbeat(self, slot, owner, fields=None, fence=None):
        """Extend a lease.  False when the lease is gone -- and it is gone.

        A false here is not a retryable error: it means another worker holds
        this run id and is walking it right now.  The supervisor's response is
        to stop, which is what `worker.py` already does with `leaseLost`.
        """
        now = int(self._clock())
        fields = dict(fields or {})
        sets = ["lease_until = ?", "updated_at = ?"]
        params = [now + self.leaseSeconds, now]
        for name in sorted(fields):
            column = _FIELD_TO_COLUMN.get(name)
            if column is None or column in ("owner", "fence", "lease_until"):
                # Unknown or protected fields are dropped rather than
                # smuggled into the SQL: a registry field name has to exist
                # in the schema to be written.
                continue
            sets.append("%s = ?" % column)
            params.append(fields[name])
        query = ("UPDATE rho_slots SET %s WHERE campaign_id = ? AND slot = ? AND owner = ? "
                 "AND lease_until >= ?" % ", ".join(sets))
        params += [self.campaign, int(slot), owner, now]
        if fence is not None:
            query += " AND fence = ?"
            params.append(int(fence))
        return bool(self.db.execute(query, tuple(params)))

    def release(self, slot, owner, state="idle", extra=None, fence=None):
        """Give the slot up.  The lease is cleared even if it already expired.

        Deliberately not conditional on the lease still being live: a worker
        shutting down after a long upload should still hand back a run id it
        no longer strictly holds, and the owner check is what stops it from
        handing back somebody else's.
        """
        extra = dict(extra or {})
        sets = ["lease_until = 0", "state = ?", "updated_at = ?"]
        params = [state, int(self._clock())]
        for name in sorted(extra):
            column = _FIELD_TO_COLUMN.get(name)
            if column is None or column in ("owner", "fence", "lease_until"):
                continue
            sets.append("%s = ?" % column)
            params.append(extra[name])
        query = ("UPDATE rho_slots SET %s WHERE campaign_id = ? AND slot = ? AND owner = ?"
                 % ", ".join(sets))
        params += [self.campaign, int(slot), owner]
        if fence is not None:
            query += " AND fence = ?"
            params.append(int(fence))
        with self.db.transaction() as tx:
            changed = tx.execute(query, tuple(params)).rowcount
            if changed:
                self._event(tx, slot, "release", owner, fence, dict(extra, state=state))
        if changed and self.cache is not None:
            self.cache.invalidateSlots()
        return bool(changed)

    # -- operator actions ---------------------------------------------------
    def retire(self, slot, reason, state="retired"):
        """Take a run id out of service, from an operator rather than a worker.

        The run id is kept: its seeds have been walked and must never be
        walked again, so the row stays and only the state changes.
        """
        now = int(self._clock())
        with self.db.transaction() as tx:
            changed = tx.execute(
                "UPDATE rho_slots SET state = ?, reason = ?, lease_until = 0, updated_at = ? "
                "WHERE campaign_id = ? AND slot = ?",
                (state, reason, now, self.campaign, int(slot))).rowcount
            if changed:
                self._event(tx, slot, state, None, None, {"reason": reason})
        if changed and self.cache is not None:
            self.cache.invalidateSlots()
        return bool(changed)

    def events(self, slot=None, limit=100):
        if slot is None:
            rows = self.db.fetchAll(
                "SELECT slot, at, kind, owner, fence, detail FROM rho_slot_events "
                "WHERE campaign_id = ? ORDER BY at DESC, slot DESC LIMIT ?",
                (self.campaign, int(limit)))
        else:
            rows = self.db.fetchAll(
                "SELECT slot, at, kind, owner, fence, detail FROM rho_slot_events "
                "WHERE campaign_id = ? AND slot = ? ORDER BY at DESC LIMIT ?",
                (self.campaign, int(slot), int(limit)))
        return [{"slot": r[0], "at": r[1], "kind": r[2], "owner": r[3], "fence": r[4],
                 "detail": json.loads(r[5]) if r[5] else None} for r in rows]

    def expired(self, now=None):
        """Slots whose lease has run out but whose state still says active.

        This is the fleet's own view of what a spot interruption took, and
        what `status` prints: a slot in this list is not lost -- it is free,
        and the next worker to claim resumes its checkpoint.
        """
        now = int(now if now is not None else self._clock())
        rows = self.db.fetchAll(
            "SELECT %s FROM rho_slots WHERE campaign_id = ? AND state = 'active' "
            "AND lease_until < ? ORDER BY slot" % ", ".join(SLOT_COLUMNS),
            (self.campaign, now))
        return [rowToItem(row) for row in rows]
