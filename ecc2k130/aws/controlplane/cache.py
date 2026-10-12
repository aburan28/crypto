"""ElastiCache (Redis) layer: a hint store, never an authority.

What the cache is for, in the order it earns its keep:

  * **Fleet heartbeats.**  At 130 workers beating every 60 s the dashboard was
    asking Postgres for the fleet on every page load.  The live fleet is a
    hash plus a sorted set here, read in one round trip, and it expires on its
    own if a worker stops beating.  Losing it costs a dashboard panel, not a
    walk.
  * **Slot-scan snapshots.**  Every claim reads the registry first.  A spot
    interruption claims a hundred slots within a few seconds, and those are a
    hundred identical scans.  The snapshot is served for a few seconds; a
    stale snapshot can only cost a claim attempt, because the claim itself is
    a conditional UPDATE the database adjudicates.
  * **Fleet-wide locks.**  One worker, not thirty, rebuilds a shared artifact
    or kicks off a merge.  The lock is advisory: losing it duplicates work,
    never corrupts it, so nothing that must happen exactly once may depend on
    holding it.
  * **The distinguished-point hot set.**  A cheap first look at whether a
    point has been seen in the last hour.  A *hit* is worth checking against
    the database; a *miss* means nothing at all, because the set holds one
    hour of one cluster's memory and the corpus holds billions of records.

Everything here fails open.  A timeout, a dropped connection, a replica
promotion, an `OOM command not allowed` -- all of them are recorded and
treated as a miss, and after a few in a row the cache is skipped entirely for
a cooldown.  A cache that is slower than the database it fronts is worse than
no cache, and 0.5 s per call against a dead cluster is exactly that.

`redis` is imported lazily inside `fromEnv`, so this module imports on a host
that does not have it.
"""

from __future__ import annotations

import json
import logging
import threading
import time
import uuid

from .config import (Config, DP_HINT_TTL_SECONDS, FLEET_TTL_SECONDS,
                     REDIS_COOLDOWN_SECONDS, REDIS_FAILURES_BEFORE_COOLDOWN,
                     REDIS_TIMEOUT_SECONDS, SLOT_SNAPSHOT_TTL_SECONDS)

logger = logging.getLogger(__name__)

# Distinguishes "the cache could not answer" from a command whose own answer
# is None -- `SET NX` returns None when the key was already there, which is
# the single most important not-an-error in this file: it is what a lock
# being held, and a point having been seen, both look like.
_UNKNOWN = object()

# Compare-and-delete: a lock is released by the holder or not at all.  Without
# the token check, a holder whose process stalled past the TTL would delete
# the lock the *next* holder now owns.
_UNLOCK = """
if redis.call('get', KEYS[1]) == ARGV[1] then
    return redis.call('del', KEYS[1])
end
return 0
"""
# Same argument for extending one.
_EXTEND = """
if redis.call('get', KEYS[1]) == ARGV[1] then
    return redis.call('pexpire', KEYS[1], ARGV[2])
end
return 0
"""


class NullCache:
    """The no-cache case, spelled as an object so callers never branch.

    `Coordinator` holds one of these when `RHO_REDIS_URL` is unset or
    `RHO_CACHE=0`.  Reading the database for everything is the control that
    any claim about the cache has to be measured against.
    """

    enabled = False

    def stats(self):
        return {"enabled": False}

    def getJson(self, key):
        return None

    def setJson(self, key, value, ttl=None):
        return False

    def delete(self, *keys):
        return False

    def beat(self, owner, payload, ttl=FLEET_TTL_SECONDS):
        return False

    def fleet(self):
        return None

    def slotSnapshot(self):
        return None

    def putSlotSnapshot(self, rows, ttl=SLOT_SNAPSHOT_TTL_SECONDS):
        return False

    def seenPoint(self, pointKey, ttl=DP_HINT_TTL_SECONDS):
        return None

    def lock(self, name, ttl=60):
        return _NullLock()

    def close(self):
        pass


class _NullLock:
    def acquire(self):
        return True

    def release(self):
        return False

    def extend(self, ttl=None):
        return False

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False


class Cache:
    """A guarded Redis client for one campaign's control-plane keys.

    Keys carry the campaign and the hash tag `{campaign}` so that cluster mode
    keeps one campaign's keys on one shard: the fleet view is read as a unit
    and a cross-slot read would need one round trip per shard.
    """

    enabled = True

    def __init__(self, client, campaign="ecc2k-130", clock=time.time,
                 cooldown=REDIS_COOLDOWN_SECONDS,
                 failuresBeforeCooldown=REDIS_FAILURES_BEFORE_COOLDOWN):
        self.client = client
        self.campaign = campaign
        self._clock = clock
        self._cooldown = cooldown
        self._limit = failuresBeforeCooldown
        self._lock = threading.Lock()
        self._failures = 0
        self._skipUntil = 0.0
        self._stats = {"hits": 0, "misses": 0, "errors": 0, "skipped": 0, "writes": 0}

    # -- construction -------------------------------------------------------
    @classmethod
    def fromEnv(cls, config=None):
        """A Cache, or a NullCache when no cluster is configured.

        Returning a different class rather than an inert Cache is deliberate:
        `cache.enabled` then answers "is there a cluster" honestly, and the
        status output can say "no cache" instead of "no hits".
        """
        config = config or Config.fromEnv()
        if not config.redisUrl or not config.cacheEnabled:
            return NullCache()
        import redis  # lazy

        kw = {}
        try:
            # redis-py 6 and later retry connection failures ten times with
            # backoff *by default*, which turns one unreachable cluster into
            # seconds of stall on a path whose whole promise is one timeout:
            # measured at 4.5s per call against a refused port, against the
            # 0.5s this sets.  `retry_on_timeout=False` does not turn that off
            # and is deprecated; a Retry with no attempts does.
            from redis.backoff import NoBackoff
            from redis.retry import Retry

            kw["retry"] = Retry(NoBackoff(), 0)
        except ImportError:                       # redis-py 4, where it is the default
            kw["retry_on_timeout"] = False

        client = redis.Redis.from_url(
            config.redisUrl,
            socket_timeout=REDIS_TIMEOUT_SECONDS,
            socket_connect_timeout=REDIS_TIMEOUT_SECONDS,
            # No retries and no reconnect-on-timeout storm: the caller has a
            # database that can answer, and it is one hop away.
            health_check_interval=30,
            **kw,
        )
        return cls(client, campaign=config.campaign)

    # -- failure handling ---------------------------------------------------
    def _skipping(self):
        with self._lock:
            if self._skipUntil and self._clock() < self._skipUntil:
                self._stats["skipped"] += 1
                return True
            return False

    def _ok(self):
        with self._lock:
            self._failures = 0
            self._skipUntil = 0.0

    def _failed(self, op, exc):
        with self._lock:
            self._failures += 1
            self._stats["errors"] += 1
            if self._failures >= self._limit:
                self._skipUntil = self._clock() + self._cooldown
                logger.warning("cache disabled for %ds after %d failures (%s: %s)",
                               self._cooldown, self._failures, op, exc)
        logger.debug("cache %s failed: %s", op, exc)

    def _call(self, op, fn, default=None):
        if self._skipping():
            return default
        try:
            out = fn()
        except Exception as exc:          # every Redis failure is a miss
            self._failed(op, exc)
            return default
        self._ok()
        return out

    def stats(self):
        with self._lock:
            out = dict(self._stats)
            out["enabled"] = True
            out["cooling"] = bool(self._skipUntil and self._clock() < self._skipUntil)
            return out

    # -- keys ---------------------------------------------------------------
    def key(self, *parts):
        return "rho:{%s}:%s" % (self.campaign, ":".join(str(p) for p in parts))

    # -- generic ------------------------------------------------------------
    def getJson(self, key):
        raw = self._call("get", lambda: self.client.get(self.key(key)))
        if raw is None:
            self._stats["misses"] += 1
            return None
        try:
            value = json.loads(raw)
        except (ValueError, TypeError) as exc:
            # Bytes the codec cannot read are a miss, never an error the
            # caller has to handle: the database has the same answer.
            self._failed("decode", exc)
            return None
        self._stats["hits"] += 1
        return value

    def setJson(self, key, value, ttl=None):
        blob = json.dumps(value, separators=(",", ":"), sort_keys=True, allow_nan=False)
        self._stats["writes"] += 1
        return bool(self._call(
            "set", lambda: self.client.set(self.key(key), blob, ex=ttl), default=False))

    def delete(self, *keys):
        names = [self.key(k) for k in keys]
        return bool(self._call("del", lambda: self.client.delete(*names), default=False))

    # -- fleet heartbeats ---------------------------------------------------
    def beat(self, owner, payload, ttl=FLEET_TTL_SECONDS):
        """Record one worker's heartbeat.

        A hash field per worker plus a sorted set scored by `last_seen`, so the
        live fleet is a range query and the dead fall out of it without
        anybody sweeping.  Both writes go in one pipeline: a heartbeat is the
        most frequent call in the system and two round trips would double the
        cache's share of it.
        """
        now = int(self._clock())
        blob = json.dumps(dict(payload, owner=owner, lastSeen=now),
                          separators=(",", ":"), sort_keys=True, allow_nan=False)

        def go():
            pipe = self.client.pipeline(transaction=False)
            pipe.hset(self.key("fleet"), owner, blob)
            pipe.zadd(self.key("fleet", "seen"), {owner: now})
            # Expire the whole view rather than each member: a campaign that
            # stops entirely should leave nothing behind, and a live one
            # refreshes this on every beat.
            pipe.expire(self.key("fleet"), ttl * 4)
            pipe.expire(self.key("fleet", "seen"), ttl * 4)
            pipe.execute()
            return True

        self._stats["writes"] += 1
        return bool(self._call("beat", go, default=False))

    def fleet(self, ttl=FLEET_TTL_SECONDS):
        """Live workers, or None when the cache cannot say.

        None and [] are different answers and the caller must be able to tell
        them apart: None means "ask the database", [] means "the cache is
        healthy and the fleet is empty".
        """
        now = int(self._clock())

        def go():
            cutoff = now - ttl
            owners = self.client.zrangebyscore(self.key("fleet", "seen"), cutoff, "+inf")
            if not owners:
                return []
            blobs = self.client.hmget(self.key("fleet"), owners)
            out = []
            for blob in blobs:
                if not blob:
                    continue
                try:
                    out.append(json.loads(blob))
                except ValueError:
                    continue
            return out

        rows = self._call("fleet", go, default=None)
        if rows is None:
            self._stats["misses"] += 1
        else:
            self._stats["hits"] += 1
        return rows

    def dropWorker(self, owner):
        def go():
            pipe = self.client.pipeline(transaction=False)
            pipe.hdel(self.key("fleet"), owner)
            pipe.zrem(self.key("fleet", "seen"), owner)
            pipe.execute()
            return True

        return bool(self._call("dropWorker", go, default=False))

    # -- slot snapshots -----------------------------------------------------
    def slotSnapshot(self):
        return self.getJson("slots:snapshot")

    def putSlotSnapshot(self, rows, ttl=SLOT_SNAPSHOT_TTL_SECONDS):
        return self.setJson("slots:snapshot", rows, ttl=ttl)

    def invalidateSlots(self):
        """Called after a claim or a release.

        Not required for correctness -- the snapshot expires in seconds and the
        database adjudicates every claim -- but it keeps `status` from showing
        a slot as free one second after somebody took it.
        """
        return self.delete("slots:snapshot")

    # -- distinguished-point hints -----------------------------------------
    def seenPoint(self, pointKey, ttl=DP_HINT_TTL_SECONDS):
        """True if this point was offered before, within the hot set's window.

        Returns None when the cache cannot answer.  The caller must treat
        True as "worth checking in Postgres" and must never treat False as
        proof that a point is new: this set remembers one hour, and a
        collision between two walks is exactly the event that could fall
        outside it.
        """
        name = self.key("dp", pointKey)
        # SET NX returns True when it created the key, i.e. when the point was
        # *not* there; invert it.  One round trip for look-and-remember.
        created = self._call("seenPoint", lambda: self.client.set(name, 1, nx=True, ex=ttl),
                             default=_UNKNOWN)
        if created is _UNKNOWN:
            return None
        # `SET NX` returns None when the key already existed; that is the hit,
        # not a failure.  Only the sentinel above means "cannot say".
        return not bool(created)

    # -- locks --------------------------------------------------------------
    def lock(self, name, ttl=60):
        return RedisLock(self, name, ttl)

    def close(self):
        try:
            self.client.close()
        except Exception:
            pass


class RedisLock:
    """An advisory fleet-wide lock.  Losing it duplicates work, never breaks it.

    There is no fencing token here on purpose: a lock that must fence is
    guarding a write, and every write in this package that needs fencing gets
    it from the slot's own fence in Postgres instead.  This exists so that
    thirty workers that miss the same expensive artifact do not all rebuild
    it.
    """

    def __init__(self, cache, name, ttl=60):
        self.cache = cache
        self.name = name
        self.ttl = ttl
        self.token = uuid.uuid4().hex
        self.held = False

    def acquire(self):
        key = self.cache.key("lock", self.name)
        got = self.cache._call(
            "lock", lambda: self.cache.client.set(key, self.token, nx=True, px=int(self.ttl * 1000)),
            default=_UNKNOWN)
        # A cache that cannot answer grants the lock: the alternative is to
        # stop doing the work because the *cache* is down, and the work is
        # what the fleet is for.  Duplicated work is the accepted cost.  A
        # None from Redis is not that case -- it means somebody holds it.
        self.held = True if got is _UNKNOWN else bool(got)
        return self.held

    def extend(self, ttl=None):
        if not self.held:
            return False
        key = self.cache.key("lock", self.name)
        ms = int((ttl or self.ttl) * 1000)
        return bool(self.cache._call(
            "extendLock", lambda: self.cache.client.eval(_EXTEND, 1, key, self.token, ms),
            default=False))

    def release(self):
        if not self.held:
            return False
        self.held = False
        key = self.cache.key("lock", self.name)
        return bool(self.cache._call(
            "unlock", lambda: self.cache.client.eval(_UNLOCK, 1, key, self.token), default=False))

    def __enter__(self):
        self.acquire()
        return self

    def __exit__(self, *exc):
        self.release()
        return False


def cacheFromEnv(config=None):
    """The cache a caller should use; a NullCache when there is no cluster."""
    return Cache.fromEnv(config or Config.fromEnv())
