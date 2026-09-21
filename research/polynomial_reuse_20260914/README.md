# Polynomial reuse for index calculus

This implementation separates reusable polynomial preprocessing from exact
algebraic reductions. Both are opt-in. Parameterized bases are an explicitly
invoked, bounded experiment, not a default solver path. The earlier unpublished
Redis implementation was removed by workspace maintenance; this change rebuilds
that capability on repository revision `01fa842`.

## What is reusable

For the last binary Semaev link, write

    S3(x,y,r) = (x+y)^2 r^2 + xy r + (xy)^2 + b.

Let A=(x+y)^2, C=xy, D=(xy)^2+b and express the target as
r=sum_k r_k z^k. Since squaring is F2-linear,

    S3(x,y,r) = D + sum_k r_k [ A (z^k)^2 + C z^k ].

The template stores D and each bracketed coefficient as coordinate Boolean
polynomials, plus every preceding target-independent chain link. Instantiation
uses XOR of polynomial lists. It does not repeat symbolic field multiplication.
The same template serves different target coordinates and both the native SAT
and Gröbner frontends. Actual factor-base membership, SAT symmetry/trace
constraints, extraction, point lifting, and relation verification remain in the
existing callers. The separate union-S4 route is not changed.

The preprocessing key includes the complete field multiplication/squaring
structure, ordered subspace basis, b, summand count and source fingerprint.
Curve coefficient a is not part of this key because S3 does not depend on a;
this cache does not certify curve membership. Exact reductions additionally
key the ordered target-specific equations, variable count, engine and degree.

## Cache layers

| Layer | Payload | Reuse scope |
|---|---|---|
| preprocessing | S3-chain template | Fresh targets in the same encoding context |
| exact-reduction | Completed Buchberger output or bounded F4 reduction plus metadata | Identical ordered Boolean system and engine settings |
| parameterized | Boolean Buchberger output with symbolic target bits | Experimental specialization across targets |

Namespaces and source fingerprints separate these payloads. Exact reductions
still consume logical node-budget units when served from cache. A failed
reduction is never cached as UNSAT. SAT uses the preprocessing template;
**SAT answers and learned clauses are not persisted by this implementation**.

Redis stores checksummed JSON with a 24-hour TTL and 4 MiB serialized-entry cap.
Reads are bounded with GETRANGE; writes use atomic SET EX. Checksums detect
accidental corruption, not malicious writers: use a private trusted cache.
Local FIFO storage defaults to 32 MiB per thread, shared among layers; it charges
serialized bytes, key bytes and an estimated per-entry overhead. This does not
bound temporary algebra or serialization allocations. Network failures use local
fallback/recomputation and a five-second reconnect cooldown. Socket timeouts are
100 ms; DNS resolution is not covered by that socket timeout.

```bash
cargo build --release --features redis-cache --bin ic
IC_PREPROCESS_CACHE=redis IC_REDUCTION_CACHE=redis \
  IC_REDIS_URL='rediss://:AUTH_TOKEN@PRIVATE_ENDPOINT:6379/0' \
  IC_REDIS_NAMESPACE=experiment-001 \
  ./target/release/ic run --degree 9 --known-log 5 --seed 101 \
    --solver groebner --batch 1 --max-trials 500 --json
```

Supply the real URL through the worker's secret mechanism; the shown token is a
placeholder. `local` selects a local-only layer; unset/off disables that layer.
`IC_CACHE_LOCAL_BYTES` configures the combined per-thread local allowance, capped
at 1 GiB (zero disables retention). Layer settings are read once per thread.
Without the Redis Cargo feature, redis mode warns and falls back locally.
ElastiCache needs a reachable cluster-mode-disabled endpoint for this client;
AWS provisioning and cloud TLS validation have not occurred.

The CLI reports `counts.algebra_cache_current_thread` in the order
`[preprocessing, exact-reduction, parameterized]`. It is **not a sum across Rayon
workers**; use `--batch 1` for this benchmark. Cache `bytes_written` counts encoded
admissible payload bytes, including local-only writes; it is not network traffic.
`overhead_ns` covers the cache operation itself, excluding caller key preparation
and computation. Full-run elapsed time includes those costs.

## Parameterized experiment

`DecompositionTemplate::parameterized_generators()` appends Boolean target
variables. `parameter_basis_cached(template, &mut AlgebraCache)` computes or
retrieves a basis in its own namespace and refuses more than 16 total variables.
Use `AlgebraCache::local(bytes).with_redis(url)` for Redis-backed offline jobs.
The caller must enforce a timeout and memory cap: the existing Buchberger routine
has no interruption hook. The experiment runner uses 30 seconds and 2 GiB per
process and retains failures.

Specialization is a homomorphism of Boolean rings. If G and F generate the same
ideal before substitution, their images generate the same ideal after assigning
any target bits. This does **not** assert that the specialized list remains a
Gröbner basis for the remaining variables. We run the solver again and check
complete solution sets against the original ANF by exhaustive enumeration.
This finite Boolean construction does not divide by target polynomials and does
not discard exceptional targets. It may, however, be far larger or more
expensive than solving individual targets.

For more general parameter fields, exceptional strata must be handled explicitly;
see the [Singular Gröbner-cover documentation](https://www.singular.uni-kl.de/Manual/4-0-3/sing_956.htm).

## Reproduce

The frozen [contract](contract.json) declares the success/failure criteria and
why a matched native suite replaces the unaffected external WDSat S4 corpus.
Read [RESULTS.md](RESULTS.md) for measured outcomes and limitations.

Build an untouched reference `ic` from `01fa842` in a separate worktree, then:

```bash
cargo build --release --features redis-cache --bin ic --example polynomial_reuse_bench
python3 research/polynomial_reuse_20260914/run.py \
  --output research/polynomial_reuse_20260914/results/run-NEW \
  --reference /absolute/path/to/reference/ic \
  --redis-server /absolute/path/to/redis-server
python3 research/polynomial_reuse_20260914/compare.py \
  research/polynomial_reuse_20260914/results/run-NEW
cargo test --release --features redis-cache --lib polynomial_reuse
cargo test --release --features redis-cache --lib algebra_cache
# Dedicated test Redis only:
IC_TEST_REDIS_URL=redis://127.0.0.1:6379/0 \
  cargo test --release --features redis-cache --lib redis_cross_client -- --ignored
```

The stage benchmark uses baseline/template/parameterized forms on identical
systems, three repetitions, disjoint development/holdout targets, F4 and SAT
complete enumeration, and exhaustive ANF truth. Setup, instantiation, solving and
validation have separate timers. F4 root-elimination XOR counts and SAT conflicts
are partial counters, not total operations. The DLP suite includes an unchanged
reference plus five cache modes, fresh namespaces and new-process replay. Its
incomplete runs remain in the comparison and block unqualified speedup claims.
