# Original C: full selected-point-order diagnostic

Both native methods verify every one of 3,738 public one-target workloads over
all 1,869 retained C models: **7,476 primary method records**, with an exact
non-clock replay of all records, fixtures and scalar truth. Source-only runs
are separate four-record controls, not extra prefix observations.

The direct control uses EC Pohlig--Hellman / small-prime BSGS. The finite-field
route uses reduced Tate/Miller with the first argument varying, verifies full
prime-power image orders, then uses multiplicative PH and checked projected
factor-base relations. No planted scalar is an argument to either solver.
Negative tests reject wrong orders and a rank-two target outside the actual
cyclic subgroup, even when it passes the annihilation test.

Exact selected point orders range from **2,600,529,297,408** to
**1,637,958,981,148,213,248**. All their 2/3/7 prime-power components are solved.
This is not a maximality or whole-group cyclicity certificate. All orders are
smooth with largest prime 7; their large numeric size is not a large-prime
security level or a new algorithmic advantage.

Native totals: 26,166 setup equations / 26,166 queries; 161,336 target digit
equations / 161,336 queries. Setup includes identity/dependent equations.
Per-component ranks, pair tests, matrix work and all phase ledgers remain in
the records. Dense coverage reflects the tiny digit groups, not a general
large-order index-calculus oracle. This is not NFS/integer-smoothness IC.

[Protocol](PROTOCOL.md), [pre-execution decision](RUN_DECISION.md),
[native implementation](diagnostic.rs), [native aggregation](summarize.rs).
Five validation tests pass. The pinned library source is 2bf7700b7; exact rlib
digests are retained. The initial compile failure is retained separately.

Primary source/prefix and their replays each have an inspected two-page report
under `../../reports/run-20261017T19*Z-C-full-point-*`. CPU 4 was reserved,
PSI was unavailable, and the full-prefix replay was contended. All elapsed
observations remain unqualified; no A/A, speedup, fastest-method or scaling
claim is made. Replays never contribute to primary means. Curve inversions,
sqrt internals and other unpriced work are explicit; primitive categories are
not a combined priced cost.

The C prefix still has only 517 locally exact conductor labels / 24 observed
strata, with 1,352 unknown. Full C remains 1,934,901,248 vertices / 165 strata.
The original full method matrix, transported targets, exact representation
registration, independent arithmetic replay and large-k backends remain open.

Raw runs are retained as one deterministic archive under `evidence/`; its
README gives extraction and hash checks. In an extracted delivery:

```sh
build-final/full-point run inputs/C.sqlite inputs/C-local-conductors.json /fresh/run prefix
build-final/full-point compare prefix /fresh/run /fresh/replay-receipt.json
build-final/full-point-summary /fresh/run /fresh/summary
```

Use `source` instead of `prefix` for the original model alone. Outputs must be
fresh. `BUILD.sh` takes explicit matching Rust dependencies and creates tests,
the diagnostic and the native summarizer. The supplied Linux x86-64 binaries
are exact measured artifacts; byte-reproducible compilation is not asserted.
