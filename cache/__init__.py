"""Read-through cache for index-calculus precomputation artifacts.

Summation polynomials and the symbolic Groebner structure depend on the
curve and not on the instance, so a fleet recomputes them once per curve
instead of once per run.  Three tiers -- an in-process LRU, Redis, and S3
as the source of truth -- and every failure path in all three degrades to
recomputation rather than to a wrong answer.

    from cache import TieredStore, summation_polynomial, json_codec

    store = TieredStore.from_env()
    s5 = summation_polynomial(store, curve, 5,
                              compute=lambda: backend.symmetrised_s5(curve),
                              codec=backend.codec)

Configuration is read from `INDEXCALC_REDIS_URL`, `INDEXCALC_S3_BUCKET`,
`INDEXCALC_S3_PREFIX`, `INDEXCALC_L1_ENTRIES` and `INDEXCALC_CACHE`; an
absent setting means an absent tier, and a store with no tiers computes.
See `cache/README.md` for the tier rationale and the ElastiCache settings
this schema assumes.
"""

from .artifacts import (
    EXPANDED,
    FORMS,
    SYMMETRIC,
    groebner_structure,
    groebner_structure_key,
    specialization,
    specialization_key,
    summation_polynomial,
    summation_polynomial_key,
)
from .keys import (
    NS_GROEBNER,
    NS_SPECIALIZATION,
    NS_SUMMATION,
    SCHEMA_VERSION,
    CacheKey,
    CacheKeyError,
    canonical_field,
    canonical_int,
    curve_fingerprint,
    inputs_fingerprint,
)
from .store import (
    ARTIFACT_TIERS,
    DEFAULT_TTL_SECONDS,
    MAX_REDIS_VALUE_BYTES,
    REDIS_TIMEOUT_SECONDS,
    SPECIALIZATION_TIERS,
    SPECIALIZATION_TTL_SECONDS,
    Codec,
    ObjectLRU,
    RedisTier,
    S3Tier,
    TieredStore,
    bytes_codec,
    json_codec,
)

__all__ = [
    "ARTIFACT_TIERS", "CacheKey", "CacheKeyError", "Codec", "DEFAULT_TTL_SECONDS",
    "EXPANDED", "FORMS", "MAX_REDIS_VALUE_BYTES", "NS_GROEBNER", "NS_SPECIALIZATION",
    "NS_SUMMATION", "ObjectLRU", "REDIS_TIMEOUT_SECONDS", "RedisTier", "S3Tier",
    "SCHEMA_VERSION", "SPECIALIZATION_TIERS", "SPECIALIZATION_TTL_SECONDS", "SYMMETRIC",
    "TieredStore", "bytes_codec", "canonical_field", "canonical_int",
    "curve_fingerprint", "groebner_structure", "groebner_structure_key",
    "inputs_fingerprint", "json_codec", "specialization", "specialization_key",
    "summation_polynomial", "summation_polynomial_key",
]
