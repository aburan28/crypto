# Call-plus-destruction ceiling has room above 2x on frozen discovery

The native Linux x86-64 AVX2 discovery completed all **48 cells and 9,408
timed F5 calls** under one uncontended physical-core reservation. The
registered four n=24, batch-32 primary groups all have a 95% bootstrap
upper bound for the median optimistic ceiling above 2.0. The structural
feasibility gate therefore **PASSES 4/4** and permits the untouched holdout
phase. It does **not** measure a candidate cache or any speedup.

| Seed and affine family | Median `U_call` | 95% upper bound | A/A 97.5% floor | Discovery screen |
|---|---:|---:|---:|---|
| 20261011, independent | 4.692 | 4.714 | 1.008 | PASS |
| 20261011, walk | 4.656 | 4.657 | 1.004 | PASS |
| 3141667, independent | 4.728 | 4.742 | 1.007 | PASS |
| 3141667, walk | 4.774 | 4.777 | 1.003 | PASS |

`U_call = (C+D)/(C+D-build-reduce)` pretends that all matrix building and
reduction work disappears. `C` includes the complete F5 call and its
materialized polynomial return; `D` separately charges returned-object
destruction. Exact output SHA-256 validation was done **outside** those
timers and remains in the whole-process receipt. Across 56 primary batch
observations, build plus reduction account for 78.78% of 142,175,332,640
timed nanoseconds. This gives an optimistic aggregate ceiling near 4.71x,
consistent with the paired groups. Validation itself took 172,800,199,158
untimed nanoseconds and is not an online cost estimate or a hidden speedup.

The inherited full-column guard used direct row packing on 8,512 calls and
the source-defined sorted-row fallback on 896; every route flag and returned
polynomial digest was checked. The CPU receipt was uncontended, with worker
wall 632.012602345 seconds, 6.4 other-process CPU seconds, preflight PSI
some avg10 3.46, SMT siblings 2 and 3 reserved, and peak whole-worker RSS
126,056 KiB. RSS includes all fixtures and references and is not the memory
cost of a proposed cache.

The [sealed discovery bundle](qualified_discovery_01/manifest.json) has
manifest SHA-256 `25ad106c29ffbe06674224ca4739e3d8cba44a84ecbeb7c4821256fc284e6557`,
raw SHA-256 `fc9a89c0d4e1b175fc59d63bd020a04f63e2e90be3edecc4545d31dd6c0af01a`,
and result SHA-256 `52d400f40ead390b2134f70fc51e72a31985e1da7879a47116da1d9aaf0e2563`.
It came from source commit `1dd364de6bdb4df30487af03d084899bb620fd31`
in [workflow 37015858282](https://github.com/aburan28/crypto/actions/runs/37015858282).
All 27 downloaded members match the manifest; native source-bound replay
passed. The new holdout seeds were not yet run when this note was written.

**Decision:** advance to fresh holdouts with unchanged worker/verifier source
and the same call-only metric. A passing upper ceiling will only establish
room for a later actual cache, whose setup, fallback, criterion, output,
memory and full-call timing still need measurement. Natural relation yield,
independent rank, full IC cost and rho ratio remain null. No curve or key
input is present.
