# Cached index-calculus precomputation artifacts

A read-through cache for the artifacts index calculus recomputes
unnecessarily: summation polynomials and the symbolic Gröbner structure.

**This package makes no performance claim.** Nothing here has been measured
against the frozen regression suite, no end-to-end pipeline comparison has
been run, and no row on `docs/index-calculus-scoreboard.html` moves because
of it. Under AGENTS.md §3 it is not yet even an *engineering* step: it is
plumbing, and it becomes a measurable change only once it is wired into the
relation-search loop, which is not in this package. When that happens the
measurement owes a baseline/candidate comparison under §8 — including the
cache's own cost, which is not zero: key construction, serialization and a
network round trip on every miss all belong in `S`.

## Tiers

| Tier | Contents | Lifetime | Churn |
|---|---|---|---|
| L1 in-process LRU | S₃…S₆, generic GB structure | process | none |
| L2 ElastiCache Redis | same, shared across the fleet | TTL 30d | none |
| L3 S3 | same, source of truth | permanent | none |
| T3 (Redis only) | individual specialization results | TTL 1h | all of it |

### Redis is a cache, not the authority

S3 holds the artifacts. If the cluster is flushed or a node is replaced
mid-run, workers refill from S3 rather than recomputing a multi-hour Gröbner
basis. This is the load-bearing decision, and it is the one thing to
preserve if the rest is rewritten: nothing in `store.py` reads a Redis miss
as evidence that an artifact does not exist.

### L1 matters more than Redis for the immutable artifacts

S₅/S₆ and the generic GB structure are read on every relation attempt. A
network round trip per attempt would dominate the specialization it is meant
to accelerate, so L1 holds *decoded objects* rather than bytes — it skips the
deserialization as well as the hop. Redis exists to stop every worker hitting
S3 on cold start, and to share T3.

### Summation polynomials are the highest-value entry

S_n depends on the curve alone — not the instance, not the target point, not
the factor base — so it is computed at most once per curve across the whole
fleet, forever. `curve_fingerprint` excludes the subgroup order, cofactor,
generator, target and name for exactly this reason; folding any of them in
would quietly turn a fleet-wide artifact into a per-run one.

The symmetrized form, in the elementary symmetric functions, is what gets
cached. It is dramatically smaller than the expanded monomial form and is
what F4 wants anyway — see `RESEARCH_SYMMETRIZED_SEMAEV.md` for the density
reduction that note records. (Those ratios are cited from it, not measured
here.)

## Correctness

**Key completeness is the only real hazard.** A key must encode everything
affecting the bytes, so `keys.py` puts the curve *representation* in the
fingerprint and not merely its isomorphism class. Two spellings of one field
must read as one field, and two genuinely different fields must never read as
one:

- `polynomial_terms` is sorted, because `src/bin/ic/params.rs` emits
  `[low_terms…, degree]` while `docs/ic/binary-example.json` writes it
  ascending. Same field, two orderings, and without sorting, two cache
  entries and a recomputed Gröbner basis every time.
- Coefficients are canonicalised through the repository's own integer
  encoding (`0x`-prefixed hex, else decimal — `number()` in `params.rs`), so
  `"0"` and `"0x0"` name one curve. Prime-field coefficients are reduced mod
  `p`; binary-field coefficients that are out of range are *rejected* rather
  than masked, because an unreduced value there means the caller holds a
  different representation than it thinks it does.
- `binary` and `binary_abstract` never collide. This repository carries both:
  `scripts/ecc2k130_point_decomposition.py` works in a polynomial basis with
  an explicit reduction polynomial, `ecc2k130/codegen/decomp.py` in the
  permuted type-II optimal normal basis. Same field, same curve, different
  bytes. For `binary_abstract` the free-text `representation` string is the
  *only* record of which basis the bytes are in, so it is normalised for
  whitespace and otherwise kept verbatim.
- The codec's name is in every key. Swapping the algebra backend swaps the
  encoding, and a parameter that changes the bytes but not the key is the one
  failure mode this module exists to prevent.

**`SCHEMA_VERSION` is the stale-artifact kill switch.** Bump it when a
serialization or algorithm changes; it feeds both the fingerprint and the key
path, so old keys are then never read and age out. Never mutate an artifact
in place under a live key — a reader already holding the old bytes has no way
to find out.

**Every tier failure degrades to recomputation, never to a wrong answer.** A
tier that raises, times out, or returns bytes the codec cannot decode is
recorded in `store.stats()` and treated as a miss. Redis timeouts are 0.5 s
and non-retrying, and a tier that fails repeatedly is skipped for a 30 s
cooldown: a cache slower than the work it replaces is worse than no cache,
and 0.5 s per relation attempt against a dead cluster is exactly that.

**The generic GB fast path is only valid while specializations stay generic.**
The caller must still run a cheap genericity check and fall back to a full
solve on a degenerate specialization. That check is *not* in this package and
cannot be: it needs the algebra backend, which this package deliberately does
not import.

## Configuration

```
INDEXCALC_REDIS_URL=rediss://<primary-endpoint>:6379   # rediss:// => TLS
INDEXCALC_S3_BUCKET=<bucket>
INDEXCALC_S3_PREFIX=indexcalc
INDEXCALC_L1_ENTRIES=64
INDEXCALC_CACHE=1
```

An absent setting means an absent tier, and a store with no tiers computes.
`INDEXCALC_CACHE=0` disables the cache entirely, which is the behaviour every
failure path above degrades to — so it is also the control to run a candidate
against.

The prefix does **not** carry the schema version: `SCHEMA_VERSION` is emitted
into the key path by `keys.py`, so a code-side bump cannot be forgotten in
configuration. Setting `INDEXCALC_S3_PREFIX=indexcalc/v1` would just produce
`indexcalc/v1/v1/…`.

### ElastiCache settings this schema assumes

- **`maxmemory-policy`: `volatile-lru`, not `allkeys-lru`.** The T1/T2
  artifacts and the T3 churn share a cluster; `allkeys` lets a relation-search
  burst evict the artifacts that cost hours to rebuild. All keys carry TTLs,
  so `volatile-lru` still bounds memory.
- **Cluster mode:** keys carry the hash tag `{curve_fp}`, so one curve's
  artifacts land on one shard and a warming worker does not fan out. S3 keys
  are the same address with the braces dropped.
- **Encryption in transit on**, via the `rediss://` scheme.
- **8 MB value cap.** A larger artifact skips Redis rather than blocking a
  worker on the transfer; S3 still takes it, so the fleet pays one S3 read per
  process instead of one per fleet.
- Sizing is driven by the number of live curves, not request rate — T1/T2 is
  tens of MB per curve at most; T3 is what actually consumes memory.

## Using it

Compute functions and serializers are injected, never imported, so this
package depends on no algebra backend:

```python
from cache import TieredStore, summation_polynomial

store = TieredStore.from_env()
s5 = summation_polynomial(store, curve, 5,
                          compute=lambda: backend.symmetrised_s5(curve),
                          codec=backend.codec)
```

`curve` is an ic parameter document (`docs/ic/binary-example.json`) or any
mapping carrying its `field`, `a` and `b`.

## Tests

```bash
python3 cache/test_cache.py
```

They run without `redis` or `boto3` installed: both are imported lazily
inside the tier that needs them, and the tiers are tested against injected
clients.

## Not here

- Wiring into the relation-search loop (needs the call sites)
- The genericity check and full-solve fallback
- A serializer for the algebra backend (injected; the caller supplies it)
- Terraform for the ElastiCache cluster
- Thundering-herd protection. Two workers that miss the same key at the same
  moment both compute it. That is fine for a specialization and wasteful for
  a multi-hour Gröbner basis, so a fleet-wide warm-up pass — or a lock — is
  worth having before the first cold start at scale.
