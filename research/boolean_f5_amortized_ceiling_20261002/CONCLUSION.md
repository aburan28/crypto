# Validation-inclusive ceiling rejects the registered gate; call-only ceiling remains open

The second, full-SMT-core Linux x86-64 AVX2 discovery passed correctness,
resource admission, all 48 frozen cells, native result verification and
post-seal replay. It recorded **9,408 timed F5 calls** on generated public
quadratic systems. The registered four-group n=24 optimistic ceiling gate is
**REJECTED: 0/4** groups have a 95% bootstrap upper bound above 2.0.
Fresh holdouts were not run.

| n=24, batch 32 | Median `T/(T-build-reduce)` | 95% upper bound | A/A 97.5% floor | Registered >2 gate |
|---|---:|---:|---:|---|
| Seed 20261008, independent affine | 1.612 | 1.613 | 1.007 | REJECTED |
| Seed 20261008, walk affine | 1.609 | 1.610 | 1.006 | REJECTED |
| Seed 3141627, independent affine | 1.611 | 1.612 | 1.004 | REJECTED |
| Seed 3141627, walk affine | 1.607 | 1.608 | 1.004 | REJECTED |

These are **ceilings**, not measured candidate/reference speedups. They
pretend that the matrix build and reduction phases cost zero. Every row and
phase record was replayed from the source. The conditional conclusion is
specific to the protocol's outer timer `T`, which charges exact SHA-256
returned-row validation and destruction as well as the F5 call.

## Accounting limit on interpretation

The validation choice materially narrows what this result can say. Across
the 56 primary n=24 batch observations, the raw phase sums contain
337,374,139,781 outer nanoseconds. Exact returned-row validation accounts
for 168,718,857,852 ns, about **50.01%** of that total. Build plus reduction
accounts for 127,787,204,917 ns, about **37.88%**. Under the frozen metric,
even removing all of that 37.88% leaves roughly 62.12% and yields an
optimistic aggregate ceiling near 1.61x, consistent with the paired gates.

That does **not** establish a 1.61x ceiling on the F5 API call alone. A
caller need not SHA-256 its entire returned polynomial list on every call.
Moving exact validation outside the timed call, while separately pricing
output destruction, changes the cost boundary and requires a **new frozen
protocol and fresh seeds**. The current holdout seeds cannot be used to
retrofit a new metric after inspecting these discovery phases. No source,
raw sample or original gate is rewritten to disguise this limitation.

The run requested direct row packing and required direct unpacking. The
source-defined full-column guard took its packed-direct path on 8,512 calls
and its sorted-row fallback on 896 calls; all remain in the result. The
whole-worker receipt was uncontended on the two SMT siblings of one physical
core, with 673.2260006 s wall, 5.05 other-process CPU seconds, PSI some
avg10 3.92 and 126,140 KiB peak RSS. Peak RSS includes the fixture corpus,
reference calls, child workers and output validation; it is not a
candidate-specific memory figure.

The [qualified bundle](qualified_discovery_01/manifest.json) has manifest
SHA-256 `f7c177888fbf7c5fca8d173940d098fd037e1269abb826656696e454b8d55468`,
raw SHA-256 `549fa88040223255ab85e8f09259b7d0bc014573221d7d2c374e061679eb156e`,
and result SHA-256 `03c8de9ecdae1a528b2e8728a9e68e78756c265c953816aed193c57863319ff5`.
It came from exact source commit `7404964eb2bb52faa667fce83d3bcc3fb3e9c766` in
[workflow 37006736293](https://github.com/aburan28/crypto/actions/runs/37006736293).
All 27 manifest members match the downloaded files; the workflow's native
replay passed. The first attempt's preflight-only SMT refusal remains sealed
separately in [ATTEMPTS.md](ATTEMPTS.md), with zero admitted timing.

Nine optimized Rust tests pass, including small F4/F5 row-space agreement,
route/fallback checks, phase accounting, deterministic bootstrap, exact
float round trips, CPU-sibling parsing and resource/refusal handling. The
registered route uses the inherited Boolean-safe F5 criterion and preserves
its public Echelon polynomial output. This is a phase-cost diagnostic on
generated systems. No actual graded F5 cache, complete-call speedup,
natural relation-yield gain, calibrated IC cost, Pollard-rho ratio, curve
input or key-related result is established.
