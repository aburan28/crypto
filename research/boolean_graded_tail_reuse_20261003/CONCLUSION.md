# Exact graded reuse is faster at n=24, but the universal discovery gate fails

The third, strongest-control Linux ARM64 discovery is complete, resource
qualified and natively replayed. It verified **200 frozen cells, 14,000 cold
arm batches and 299,600 output matrices**. The candidate reuses the identical
cubic Macaulay block across changing affine tails, and every returned echelon
matched fresh full elimination exactly. The registered all-eight discovery
screen is **REJECTED: 2/8 groups** have a paired 95% lower timing bound above
1.5 and their A/A floor. Both passing groups are at n=24. No fresh holdout
was executed, so the primary holdout 2x gate is unknown.

| Original variables | Independent affine, median [lower] | Walk affine, median [lower] | Groups above 1.5 |
|---:|---:|---:|---:|
| 12 | 1.041 [1.034] | 1.070 [1.065] | 0/2 |
| 16 | 1.336 [1.327] | 1.366 [1.354] | 0/2 |
| 20 | 1.432 [1.423] | 1.471 [1.464] | 0/2 |
| 24 | 1.779 [1.760] | 1.826 [1.819] | 2/2 |

Each dimensionless ratio is the pointwise fastest correct fresh control's
cold batch cost divided by the graded arm's cost for the same fixed
seed/repetition. The table copies `qualified_discovery_01/results.json` and
rounds displayed values only; its gate uses the full-precision paired lower
bounds and A/A 97.5th-percentile floors. The intervals describe repeated
measurements on two fixed discovery seeds, not a population of Boolean
systems. There is **no** universal 2x or full-method claim.

## Matched cold costs at n=24

The following milliseconds are pooled medians of twenty complete batch-64
observations per family. Every arm receives its own cold setup, constructs or
reuses all 64 matrices, eliminates, materializes, verifies and destroys
outputs. The table is descriptive; acceptance uses paired pointwise minima,
not ratios of these pooled medians. Retained bytes count the arm's stored
context capacity and map payload, excluding allocator metadata.

| Arm | Independent affine (ms) | Walk affine (ms) | Packed direct / arm, independent | Retained context bytes | Exact output |
|---|---:|---:|---:|---:|---|
| Packed direct, 16-bit dense lookup | 110.711 | 110.056 | 1.000 | 33,600,796 | PASS |
| Ranked dense intermediate | 180.925 | 177.284 | 0.612 | 46,332 | PASS |
| Ranked sparse intermediate | 217.533 | 213.353 | 0.509 | 46,332 | PASS |
| Completed-matrix cache | 181.150 | 177.783 | 0.611 | 257,436 | PASS |
| **Graded high-block reuse** | **62.250** | **60.218** | **1.779** | **1,831,020** | PASS |

The high-block identity is exact: for quadratic `f_j=q_j+a_j` and any admitted
degree-at-most-one multiplier `t`, the cubic projection of `t*f_j` is
`pi_3(t*q_j)`. The candidate compiles that block once, records row swaps and
XORs, transforms each changing low block and completes low-column elimination.
Its context guard rebuilds when quadratic support changes. All fixed setup,
failed guard work, schedule application, final output and validation are
charged. A finite n=24 stage gain says nothing by itself about the fraction
of complete Boolean solving spent here, natural relation yield, or a
Pollard-rho crossover.

## Provenance and correction chain

The admitted [discovery bundle](qualified_discovery_01/manifest.json) has
manifest SHA-256 `b7a41ed0c59926ce3fa93990c2b4c950865e821fc866a6eb8144b4403965422b`,
raw SHA-256 `aca5d8837a981311d428b607e8a2c35f04d18745b684a31b44dc7c38b789b6b2`,
and result SHA-256 `f94e7e2e865c26020c8cd01b975652d81a669d42f56a98a4bc002d02dadaeee7`.
GitHub [run 36985293072](https://github.com/aburan28/crypto/actions/runs/36985293072)
used exact source commit `0f226534462ee2e0a5deb41b7f8dfe0531ef73d2`.
Every one of its 26 manifest members was checked after download; the native
workflow replay and an untimed cross-host native replay both passed. Its
one-core receipt reports 261.678709939 worker
seconds, 2.78 other-process CPU seconds, PSI some avg10 3.66, no contention,
and whole-worker peak RSS 52,504 KiB. RSS includes every arm and common
references; it is not candidate-specific peak memory.

Two earlier attempts remain unchanged in [ATTEMPTS.md](ATTEMPTS.md) and
[RUNS.json](RUNS.json). The first failed a post-seal one-unit float readback;
the second completed but used a 32-bit dense table that unnecessarily pushed
the n=24 packed control above the unchanged 64 MiB cap. The now included
16-bit table fits that cap. Its earlier apparent n=24 ~2.9x comparison was
not against the strongest applicable constructor and is not reused. Neither
attempt supplied a promoted result or consumed holdouts.

Eleven optimized Rust tests pass: 1,024 small affine assignments, all five
declared sizes and four families, every ambient rank/dense coordinate,
quadratic-support fallback, exact echelon equality, changed-output refusal,
resource contention refusal, float round trips and sealed failure handling.
Strict Clippy and formatting pass. The producer, analysis and replay are Rust;
the existing Python isolation controller is used only to reserve and check
the Linux CPU. No curve or key input is present. Full polynomial solving,
calibrated IC operations, independent relation yield, memory at
cryptanalytic scale and rho ratio remain null. No cryptanalytic breakthrough
has been established.
