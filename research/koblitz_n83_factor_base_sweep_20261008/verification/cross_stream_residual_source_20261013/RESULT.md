# Cross-stream signed-orbit residual constructor: source validation

The paired v3 design covers nested public-x factor-base streams only. The
library now has `build_exact_residual_difference_factor_base` for a second
case: given two already replayed same-curve signed-Frobenius bases, it builds
the residual domain from every **complete** signed orbit in the pool that is
absent from the small base. The small base need not be a prefix of the pool,
and a shared orbit is removed in full. The returned base can be passed to
`MixedOrbitDomains::new` and the existing one/two-large-prime graph path.

The constructor reconstructs both input orbit maps from their point records,
rejects a partial-orbit overlap or changed phase/sign label, reconstructs the
residual map, and checks that the small and residual x-domains are disjoint.
The residual keeps the pool's signed-orbit order, which may differ from the
row order used by the older exact-prefix suffix constructor. A caller must
bind and replay the compressed S3 objects before using this library helper;
the current N83 worker still imports only the two pinned nested pairs.

The focused degree-7 control constructs a non-nested pool with one shared and
two residual signed orbits. It compares the resulting exact point set to an
independently built residual, checks domain admission, compares the nested
case against the older suffix constructor as a set, and rejects a damaged
orbit label, an incomplete point record, and a pool with no residual domain.
The focused release-library test passed: one selected test, zero failures.
Two initial assertions compared row order across different valid orbit
orderings; they were corrected to compare exact point sets and orbit counts.

The earlier [known-witness chain replay](../mixed-orbit-control-chain-preflight-20261013/RESULT.md)
and [paired-capacity screen](../paired-capacity-screen-20261013/RESULT.md)
now hash the worker's **frozen commit blobs** instead of the current working
tree. Their recorded frozen-source hashes and scientific counters are
unchanged, their regenerated receipts pass byte-for-byte checks, and both
focused Python suites pass. This allows later library source work without
mistaking it for a change in the already built K600 worker.

```sh
cargo test --release --lib cross_stream_difference_keeps_only_complete_disjoint_signed_orbits -- --nocapture
```

This is a source-availability step. No cross-stream N83 S3 object pair was
imported, no SAT model or search was launched, and no graph or total-runtime
comparison was measured. It does not change the conservative 17,400-second
local ledger or identify a best factor base. The next empirical gate would
need exact object identities, compressed hashes, replayed overlap counts, a
bounded importer and independent audit before a natural-target search.
