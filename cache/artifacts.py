"""Cached accessors for the artifacts index calculus keeps recomputing.

Three families, and they are cached for different reasons:

`summation_polynomial` is the highest-value entry in the cache.  S_n depends
on the curve equation and the field alone -- not on the instance, the target
point, the factor base or the seed -- so across the whole fleet, forever, it
is computed at most once per curve.  The symmetrised form, expressed in the
elementary symmetric functions, is what gets stored: `RESEARCH_SYMMETRIZED_SEMAEV.md`
records the density reduction it buys (~10x at S_5, ~100x at S_8), and it is
also the form F4 wants, so caching the expanded form would be storing the
larger object in order to throw it away.

`groebner_structure` caches the symbolic shape of the system -- the Macaulay
pattern, the degrees, the monomial order, whatever the backend's "solve this
shape again with different constants" object is.  Read on every relation
attempt, which is why L1 carries it as a decoded object.

`specialization` is the T3 tier: a short-lived memo of one specialization
during relation search.  Redis only, short TTL, no S3.  It is churn -- worth
sharing between the workers of one run and worthless the next day -- and it
is what actually consumes cluster memory, which is why the cluster wants
`volatile-lru` rather than `allkeys-lru`: a relation-search burst must not
evict the artifacts above, which cost hours to rebuild.

No compute function and no serializer is imported here.  The caller injects
both, so this package depends on no algebra backend and the backend is free
to change without this one knowing.  The codec's `name` is part of every key
for exactly that reason: swapping the backend swaps the bytes.
"""

from __future__ import annotations

from typing import Any, Callable, Mapping, Optional

from .keys import (
    NS_GROEBNER,
    NS_SPECIALIZATION,
    NS_SUMMATION,
    CacheKey,
    CacheKeyError,
    curve_fingerprint,
    inputs_fingerprint,
)
from .store import (
    ARTIFACT_TIERS,
    DEFAULT_TTL_SECONDS,
    SPECIALIZATION_TIERS,
    SPECIALIZATION_TTL_SECONDS,
    Codec,
    TieredStore,
)

SYMMETRIC = "symmetric"
EXPANDED = "expanded"
FORMS = (SYMMETRIC, EXPANDED)


def summation_polynomial_key(curve: Mapping[str, Any], n: int, codec_name: str,
                             form: str = SYMMETRIC) -> CacheKey:
    if not isinstance(n, int) or isinstance(n, bool) or n < 2:
        raise CacheKeyError(f"summation polynomial index must be an integer >= 2, got {n!r}")
    if form not in FORMS:
        raise CacheKeyError(f"form must be one of {FORMS}, got {form!r}")
    return CacheKey.build(NS_SUMMATION, curve_fingerprint(curve),
                          n=n, form=form, codec=codec_name)


def summation_polynomial(store: TieredStore, curve: Mapping[str, Any], n: int, *,
                         compute: Callable[[], Any], codec: Codec,
                         form: str = SYMMETRIC,
                         ttl: Optional[int] = DEFAULT_TTL_SECONDS) -> Any:
    """Return S_n for `curve`, computing it once per curve per fleet.

    `curve` is an ic parameter document (`docs/ic/binary-example.json`) or
    any mapping carrying its `field`, `a` and `b`.  `compute` takes no
    arguments and returns the polynomial in whatever form `codec` encodes.
    """
    key = summation_polynomial_key(curve, n, codec.name, form)
    return store.get_or_compute(key, compute, codec, tiers=ARTIFACT_TIERS, ttl=ttl)


def groebner_structure_key(curve: Mapping[str, Any], codec_name: str, *,
                           summands: int, factor_base_dim: int,
                           monomial_order: str, algorithm: str,
                           extra: Optional[Mapping[str, Any]] = None) -> CacheKey:
    if not isinstance(summands, int) or isinstance(summands, bool) or summands < 2:
        raise CacheKeyError(f"summands must be an integer >= 2, got {summands!r}")
    if not isinstance(factor_base_dim, int) or isinstance(factor_base_dim, bool) or factor_base_dim < 1:
        raise CacheKeyError(f"factor_base_dim must be a positive integer, got {factor_base_dim!r}")
    params = dict(m=summands, fbdim=factor_base_dim,
                  order=monomial_order, alg=algorithm, codec=codec_name)
    # `extra` is the escape hatch for backend-specific knobs, and it is
    # fingerprinted into the key rather than dropped: a parameter that
    # changes the bytes and not the key is the one failure mode this whole
    # module is built to avoid.
    if extra:
        params["opts"] = inputs_fingerprint(extra)
    return CacheKey.build(NS_GROEBNER, curve_fingerprint(curve), **params)


def groebner_structure(store: TieredStore, curve: Mapping[str, Any], *,
                       compute: Callable[[], Any], codec: Codec,
                       summands: int, factor_base_dim: int,
                       monomial_order: str, algorithm: str,
                       extra: Optional[Mapping[str, Any]] = None,
                       ttl: Optional[int] = DEFAULT_TTL_SECONDS) -> Any:
    """Return the generic Groebner structure for this system shape.

    **The caller still owes a genericity check.**  The structure returned
    here is the one a *generic* specialization has; a degenerate one -- a
    dropped leading term, a rank collapse, a vanishing resultant -- does not
    have it, and using the cached shape there yields a wrong solve rather
    than a slow one.  Run the cheap genericity test on each specialization
    and fall back to a full solve when it fails.  That check is the caller's
    and is deliberately not in this package, which knows nothing about the
    algebra backend.
    """
    key = groebner_structure_key(curve, codec.name, summands=summands,
                                 factor_base_dim=factor_base_dim,
                                 monomial_order=monomial_order,
                                 algorithm=algorithm, extra=extra)
    return store.get_or_compute(key, compute, codec, tiers=ARTIFACT_TIERS, ttl=ttl)


def specialization_key(curve: Mapping[str, Any], codec_name: str,
                       inputs: Any) -> CacheKey:
    return CacheKey.build(NS_SPECIALIZATION, curve_fingerprint(curve),
                          h=inputs_fingerprint(inputs), codec=codec_name)


def specialization(store: TieredStore, curve: Mapping[str, Any], *,
                   inputs: Any, compute: Callable[[], Any], codec: Codec,
                   ttl: Optional[int] = SPECIALIZATION_TTL_SECONDS) -> Any:
    """Memoize one specialization result across the workers of a run.

    `inputs` is everything the result depends on beyond the curve -- the
    specialised coordinates, the factor-base slice, the target -- and is
    fingerprinted, so it may be any mapping, sequence or scalar with a
    stable `repr`.  Redis only: this tier is churn, and S3 should not carry
    it.
    """
    key = specialization_key(curve, codec.name, inputs)
    return store.get_or_compute(key, compute, codec,
                                tiers=SPECIALIZATION_TIERS, ttl=ttl)
