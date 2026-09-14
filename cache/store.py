"""A tiered read-through store for index-calculus precomputation artifacts.

Three tiers, fastest first:

    L1  in-process LRU of *decoded* objects   process lifetime
    L2  Redis, shared across the fleet        TTL, default 30 days
    L3  S3, the source of truth               no expiry

Redis is a cache and S3 is the authority.  If the cluster is flushed or a
node is replaced mid-run, workers refill from S3 rather than recomputing a
multi-hour Groebner basis; nothing in this module treats a Redis miss as
evidence that an artifact does not exist.

L1 holds objects, not bytes, and that is the point of it.  S_5, S_6 and the
generic Groebner structure are read on every relation attempt, so a network
round trip -- or even a repeated deserialization -- per attempt would
dominate the specialization it is meant to accelerate.  L2 exists to stop
every worker in the fleet hitting S3 on cold start, and to share the
short-lived T3 specialization memos that no one wants in S3 at all.

Every failure path degrades to recomputation and never to a wrong answer.
A tier that raises, times out, or returns bytes the codec cannot decode is
recorded and treated as a miss.  Redis operations use a 0.5 s timeout and
do not retry, and a tier that fails repeatedly is skipped for a cooldown:
a cache that is slower than the work it replaces is worse than no cache,
and 0.5 s per relation attempt against a dead cluster is exactly that.

`redis` and `boto3` are imported lazily, inside the tier that needs them,
so this module imports and its tests run on a machine that has neither.
"""

from __future__ import annotations

from collections import OrderedDict
from dataclasses import dataclass
import json
import logging
import os
import threading
import time
from typing import Any, Callable, Mapping, Optional, Sequence

from .keys import CacheKey

logger = logging.getLogger(__name__)

# Artifacts are immutable and cheap to keep; the TTL exists so that a key we
# stop writing eventually stops costing memory, not to bound staleness.
# Staleness is handled by SCHEMA_VERSION, which changes the key instead.
DEFAULT_TTL_SECONDS = 30 * 24 * 3600

# T3 memoizes individual specialization results during one relation search.
# It is churn: useful across the workers of a single run, worthless a day
# later, and never worth an S3 round trip.
SPECIALIZATION_TTL_SECONDS = 3600

REDIS_TIMEOUT_SECONDS = 0.5

# Above this, a value skips Redis entirely rather than blocking a worker on
# a multi-megabyte SET.  S3 still takes it, so the artifact is not lost --
# the fleet just pays one S3 read per process instead of one per fleet.
MAX_REDIS_VALUE_BYTES = 8 * 1024 * 1024

L1 = "l1"
L2 = "l2"
L3 = "l3"

ARTIFACT_TIERS = (L1, L2, L3)
SPECIALIZATION_TIERS = (L2,)

_BREAKER_THRESHOLD = 3
_BREAKER_COOLDOWN_SECONDS = 30.0


@dataclass(frozen=True)
class Codec:
    """How an artifact becomes bytes, and how it comes back.

    `name` goes into the cache key: two codecs produce different bytes for
    the same object, so they must not share an entry.  The algebra backend
    supplies these -- this package never imports one.
    """

    name: str
    encode: Callable[[Any], bytes]
    decode: Callable[[bytes], Any]


def json_codec(name: str = "json") -> Codec:
    """A codec for plain JSON-able artifacts; handy for tests and metadata."""
    return Codec(
        name=name,
        encode=lambda obj: json.dumps(obj, sort_keys=True, separators=(",", ":")).encode("utf-8"),
        decode=lambda blob: json.loads(blob.decode("utf-8")),
    )


def bytes_codec(name: str = "raw") -> Codec:
    """A codec for callers that already hold the serialized form."""
    def _encode(obj: Any) -> bytes:
        if not isinstance(obj, (bytes, bytearray)):
            raise TypeError(f"bytes_codec expects bytes, got {type(obj).__name__}")
        return bytes(obj)
    return Codec(name=name, encode=_encode, decode=bytes)


class _Breaker:
    """Skip a tier that keeps failing, and try it again after a cooldown."""

    def __init__(self, threshold: int = _BREAKER_THRESHOLD,
                 cooldown: float = _BREAKER_COOLDOWN_SECONDS,
                 clock: Callable[[], float] = time.monotonic):
        self._threshold, self._cooldown, self._clock = threshold, cooldown, clock
        self._failures = 0
        self._open_until = 0.0
        self._lock = threading.Lock()

    def closed(self) -> bool:
        with self._lock:
            return self._clock() >= self._open_until

    def record_success(self) -> None:
        with self._lock:
            self._failures = 0
            self._open_until = 0.0

    def record_failure(self) -> None:
        with self._lock:
            self._failures += 1
            if self._failures >= self._threshold:
                self._open_until = self._clock() + self._cooldown


# -- tiers ----------------------------------------------------------------

class ObjectLRU:
    """L1: a bounded, thread-safe LRU of decoded artifacts."""

    name = L1

    def __init__(self, max_entries: int = 64):
        if max_entries < 1:
            raise ValueError("max_entries must be at least 1")
        self.max_entries = max_entries
        self._entries: "OrderedDict[Any, Any]" = OrderedDict()
        self._lock = threading.Lock()

    def get(self, key: CacheKey) -> Any:
        with self._lock:
            k = key.l1_key()
            if k not in self._entries:
                return None
            self._entries.move_to_end(k)
            return self._entries[k]

    def put(self, key: CacheKey, value: Any) -> None:
        with self._lock:
            k = key.l1_key()
            self._entries[k] = value
            self._entries.move_to_end(k)
            while len(self._entries) > self.max_entries:
                self._entries.popitem(last=False)

    def clear(self) -> None:
        with self._lock:
            self._entries.clear()

    def __len__(self) -> int:
        with self._lock:
            return len(self._entries)


class RedisTier:
    """L2: shared, TTL'd, and never authoritative.

    The client is built lazily so that importing this module does not
    require `redis`, and so that a process which never takes a cache miss
    never opens a connection.
    """

    name = L2

    def __init__(self, url: str, *, client: Any = None,
                 timeout: float = REDIS_TIMEOUT_SECONDS,
                 max_value_bytes: int = MAX_REDIS_VALUE_BYTES):
        self.url = url
        self.timeout = timeout
        self.max_value_bytes = max_value_bytes
        self._client = client
        self._lock = threading.Lock()

    def client(self) -> Any:
        if self._client is None:
            with self._lock:
                if self._client is None:
                    import redis  # lazy: optional dependency

                    # `rediss://` selects TLS through `from_url`; encryption in
                    # transit is a cluster setting and a URL scheme, not a flag
                    # here.  No retry: one slow round trip is the whole budget.
                    self._client = redis.Redis.from_url(
                        self.url,
                        socket_timeout=self.timeout,
                        socket_connect_timeout=self.timeout,
                        retry_on_timeout=False,
                    )
        return self._client

    def get(self, key: CacheKey) -> Optional[bytes]:
        return self.client().get(key.redis_key())

    def put(self, key: CacheKey, blob: bytes, ttl: Optional[int]) -> bool:
        if len(blob) > self.max_value_bytes:
            # Deliberately not chunked.  A caller that wants a 40 MB Groebner
            # structure shared gets it from S3; putting it in Redis buys one
            # avoided S3 read per process and costs every worker the transfer.
            logger.debug("redis: skipping %s, %d bytes over the %d cap",
                         key, len(blob), self.max_value_bytes)
            return False
        self.client().set(key.redis_key(), blob, ex=ttl)
        return True


class S3Tier:
    """L3: the source of truth.  No TTL -- artifacts here are permanent."""

    name = L3

    def __init__(self, bucket: str, prefix: str = "indexcalc", *, client: Any = None):
        self.bucket = bucket
        self.prefix = prefix
        self._client = client
        self._lock = threading.Lock()

    def client(self) -> Any:
        if self._client is None:
            with self._lock:
                if self._client is None:
                    import boto3  # lazy: optional dependency

                    self._client = boto3.client("s3")
        return self._client

    def get(self, key: CacheKey) -> Optional[bytes]:
        client = self.client()
        try:
            response = client.get_object(Bucket=self.bucket,
                                         Key=key.s3_key(self.prefix))
        except Exception as exc:  # NoSuchKey is a miss, not a tier failure
            if _is_missing_object(exc):
                return None
            raise
        return response["Body"].read()

    def put(self, key: CacheKey, blob: bytes, ttl: Optional[int]) -> bool:
        del ttl  # S3 entries do not expire; lifecycle rules are the operator's
        self.client().put_object(Bucket=self.bucket,
                                 Key=key.s3_key(self.prefix),
                                 Body=blob)
        return True


# `NoSuchKey`/`404` mean the artifact is not cached yet, which is an ordinary
# miss.  Anything else -- credentials, throttling, a missing bucket -- is a
# tier failure, and the store must see it as one so the breaker can open.
_S3_MISS_CODES = frozenset({"NoSuchKey", "404"})


def _is_missing_object(exc: BaseException) -> bool:
    """Tell 'not cached' apart from 'S3 is unwell', without importing botocore."""
    response = getattr(exc, "response", None)
    if isinstance(response, Mapping):
        error = response.get("Error")
        if isinstance(error, Mapping):
            return str(error.get("Code")) in _S3_MISS_CODES
    return type(exc).__name__ == "NoSuchKey"


# -- the store ------------------------------------------------------------

@dataclass
class TierStats:
    hit: int = 0
    miss: int = 0
    error: int = 0
    skip: int = 0
    write: int = 0

    def as_dict(self) -> dict:
        return {"hit": self.hit, "miss": self.miss, "error": self.error,
                "skip": self.skip, "write": self.write}


class TieredStore:
    """Read-through across L1/L2/L3, with every tier optional.

    A store with no tiers at all is legal and simply computes; that is what
    `INDEXCALC_CACHE=0` produces, and it is the behaviour every failure
    path degrades to.
    """

    def __init__(self, *, l1: Optional[ObjectLRU] = None,
                 l2: Any = None, l3: Any = None, enabled: bool = True):
        self._tiers = {L1: l1, L2: l2, L3: l3}
        self.enabled = enabled
        self._stats = {name: TierStats() for name in ARTIFACT_TIERS}
        self._stats["compute"] = TierStats()
        self._breakers = {name: _Breaker() for name in (L2, L3)}

    # -- construction

    @classmethod
    def from_env(cls, env: Optional[Mapping[str, str]] = None) -> "TieredStore":
        """Build a store from `INDEXCALC_*`.  Absent config means absent tier."""
        env = os.environ if env is None else env
        enabled = env.get("INDEXCALC_CACHE", "1").strip().lower() not in ("0", "false", "no", "")
        if not enabled:
            return cls(enabled=False)

        # A malformed cache setting degrades like every other cache failure:
        # to the default, never to a dead run.
        raw_entries = env.get("INDEXCALC_L1_ENTRIES", "64")
        try:
            entries = int(raw_entries)
        except (TypeError, ValueError):
            logger.warning("INDEXCALC_L1_ENTRIES=%r is not an integer; using 64", raw_entries)
            entries = 64
        l1 = ObjectLRU(max_entries=entries) if entries > 0 else None

        redis_url = env.get("INDEXCALC_REDIS_URL", "").strip()
        l2 = RedisTier(redis_url) if redis_url else None

        bucket = env.get("INDEXCALC_S3_BUCKET", "").strip()
        prefix = env.get("INDEXCALC_S3_PREFIX", "indexcalc").strip()
        l3 = S3Tier(bucket, prefix) if bucket else None

        return cls(l1=l1, l2=l2, l3=l3, enabled=True)

    # -- introspection

    def tier(self, name: str) -> Any:
        return self._tiers.get(name)

    def stats(self) -> dict:
        return {name: s.as_dict() for name, s in self._stats.items()}

    def reset_stats(self) -> None:
        for s in self._stats.values():
            s.hit = s.miss = s.error = s.skip = s.write = 0

    # -- the read-through

    def get_or_compute(self, key: CacheKey, compute: Callable[[], Any],
                       codec: Codec, *,
                       tiers: Sequence[str] = ARTIFACT_TIERS,
                       ttl: Optional[int] = DEFAULT_TTL_SECONDS) -> Any:
        """Return the artifact at `key`, computing and populating it on a miss.

        `compute` is called at most once per miss and is never called with a
        tier lock held.  Whatever it returns is what the caller gets, even
        if every write below it fails.
        """
        if not self.enabled:
            self._stats["compute"].miss += 1
            return compute()

        active = [name for name in tiers if self._tiers.get(name) is not None]

        # Read fastest-first; remember which tiers missed so a hit further
        # down can backfill them.
        missed: list[str] = []
        for name in active:
            value, blob, found = self._read(name, key, codec)
            if found:
                # Reuse the bytes the hit came from rather than re-encoding:
                # cheaper, and it keeps a non-deterministic codec from writing
                # a different encoding of the same artifact into a faster tier.
                self._write_all(missed, key, value, codec, ttl, encoded=blob)
                return value
            missed.append(name)

        self._stats["compute"].miss += 1
        value = compute()
        self._write_all(active, key, value, codec, ttl)
        return value

    # -- internals

    def _read(self, name: str, key: CacheKey,
              codec: Codec) -> tuple[Any, Optional[bytes], bool]:
        """Return `(value, raw_bytes, found)`; `raw_bytes` is None for L1."""
        tier = self._tiers[name]
        stats = self._stats[name]

        if name == L1:
            value = tier.get(key)
            if value is None:
                stats.miss += 1
                return None, None, False
            stats.hit += 1
            return value, None, True

        breaker = self._breakers[name]
        if not breaker.closed():
            stats.skip += 1
            return None, None, False

        try:
            blob = tier.get(key)
        except Exception:
            stats.error += 1
            breaker.record_failure()
            logger.warning("cache tier %s failed reading %s; recomputing", name, key,
                           exc_info=logger.isEnabledFor(logging.DEBUG))
            return None, None, False
        breaker.record_success()

        if blob is None:
            stats.miss += 1
            return None, None, False

        try:
            value = codec.decode(blob)
        except Exception:
            # An undecodable value is a miss.  We do not delete it: the next
            # write under this key overwrites it, and a SCHEMA_VERSION bump
            # strands it to age out on its own.
            stats.error += 1
            logger.warning("cache tier %s held an undecodable value for %s; recomputing",
                           name, key, exc_info=logger.isEnabledFor(logging.DEBUG))
            return None, None, False

        stats.hit += 1
        return value, blob, True

    def _write_all(self, names: Sequence[str], key: CacheKey, value: Any,
                   codec: Codec, ttl: Optional[int],
                   encoded: Optional[bytes] = None) -> None:
        byte_tiers = [n for n in names if n != L1]

        if L1 in names and self._tiers.get(L1) is not None:
            try:
                self._tiers[L1].put(key, value)
                self._stats[L1].write += 1
            except Exception:
                self._stats[L1].error += 1
                logger.warning("cache tier l1 failed writing %s", key,
                               exc_info=logger.isEnabledFor(logging.DEBUG))

        if not byte_tiers:
            return

        if encoded is None:
            try:
                encoded = codec.encode(value)
            except Exception:
                # Serialization failed: the caller still gets its value, and
                # nothing below L1 is populated.  Charged to the slowest tier
                # we were going to write, since that is where it shows up.
                self._stats[byte_tiers[-1]].error += 1
                logger.warning("codec %s could not encode %s; leaving it uncached",
                               codec.name, key,
                               exc_info=logger.isEnabledFor(logging.DEBUG))
                return

        # Slowest first: S3 is the authority, so if only one write survives
        # it should be the durable one.
        for name in (L3, L2):
            if name not in byte_tiers:
                continue
            tier = self._tiers.get(name)
            if tier is None:
                continue
            breaker = self._breakers[name]
            if not breaker.closed():
                self._stats[name].skip += 1
                continue
            try:
                if tier.put(key, encoded, ttl):
                    self._stats[name].write += 1
                else:
                    self._stats[name].skip += 1
            except Exception:
                self._stats[name].error += 1
                breaker.record_failure()
                logger.warning("cache tier %s failed writing %s", name, key,
                               exc_info=logger.isEnabledFor(logging.DEBUG))
            else:
                breaker.record_success()
