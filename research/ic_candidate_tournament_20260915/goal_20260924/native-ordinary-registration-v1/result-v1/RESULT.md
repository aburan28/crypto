# Original audited natural-query panels, 2026-10-05

Both separately sealed, one-use native registrations completed their fixed
512-query panels. Each **original frozen checker** passed source, binary,
publication, worker/role-custody, input-law, geometric, model, row and timing
replay. An audit pass means that the reported outcome is trustworthy; it does
not mean the preparation reached full rank.

| Stage diagnostic on the same ordered n17 points | Bounded MatrixF5 | External CryptoMiniSat |
| --- | ---: | ---: |
| Verified witnesses / 512 | 129 / 512 | 22 / 512 |
| Witness rate; Wilson 95% interval | 25.20%; 21.63–29.13% | 4.30%; 2.85–6.42% |
| Other outcomes | 383 proved geometric negatives | 490 conflict-budget inconclusives |
| Feasible queries left inconclusive | 0 | 107 |
| Accepted relation rows | 127 | 22 |
| Final rank / folded columns | **29 / 29** | **22 / 29** |
| Independently replayed column logs | 29 | none |
| Charged preparation wall | 3,683.65 s | 2,142.39 s |
| Charged preparation wall / accepted row | 29.01 s | 97.38 s |
| Included PDP phase, all attempts | 3,643.84 s | 2,123.23 s |

The independently audited pairing has 22 witnesses in both arms, 107 in F5
only, none in CMS only and 383 in neither. The paired F5-minus-CMS witness
rate is 20.90 percentage points, with a 10,000-resample paired bootstrap 95%
interval of 17.38–24.41 points. CMS-minus-F5 mean observed attempt wall is
-2.969 s (paired bootstrap interval -3.456 to -2.490 s), **an exploratory
stage cost difference, not a speedup**: the registered CMS budget misses 107
feasible queries and cannot solve its 29-column log system. Treating its
shorter incomplete attempts as wins would reverse the scientific conclusion.
The CMS executable reached its frozen 100,000-conflict cap on all 490
inconclusive queries. Its 30-second per-role watchdog did not turn those into
proved negatives. The F5 arm is a degree-3, 8,192-node MatrixF5-criterion
collector, not a claim of complete incremental F5.

The fixed curve is `EC1N17Ce1hdfbf24105ef5`, with 63 geometric base points,
62 subgroup-usable points and 29 folded columns. The two panels use the same
ordered fixture scalars and point law; source-instance nonces differ by arm.
All 512 attempts, failures and costs are in
[`stage-summary.json`](stage-summary.json), produced by the native,
source-bound [paired analysis](../../../../../examples/ic_ordinary_stage_summary.rs).
The report checks each original audit and producer SHA-256, complete attempt
counts, exact trial/scalar pairing, outcome counts and exclusive phase closure.
Wilson intervals describe witness rates under the declared query law; paired
bootstrap resamples whole trial pairs with seed 20261005. The per-attempt wall
comparison comes from an ordinary, unisolated macOS ARM64 host and remains
exploratory. Peak memory is unrecorded, so neither a memory claim nor an
isolation-qualified CPU ratio is available.

The full original execution trees, including all role stdout/stderr, exporter
sources, PID ledgers, progress files, producer and terminal records, are
retained in [F5 execution](f5-execution.tar.gz) and
[CMS execution](cms-execution.tar.gz). Archive SHA-256 values are
`ec670533381cb142da250c94a27a04f9cbd26de29c60f80d00c411be1d5b2c12`
(F5, 1,035 tar entries) and
`d7d4a49332a24be658338d7bd20cef58eedf7a6ff6eddcb867410fd5955c5e7f`
(CMS, 8,716 entries). The original frozen
[F5 audit](f5-original-audit.json) and
[CMS audit](cms-original-audit.json) have SHA-256
`fe8231687eda459a99ae817ee97d683c9fa561457df098bc1e2c9f895faa59a0`
and
`283ada0f6c81d3d5e9df8dc98310e9b0dd10df269e0451dce9865bc4e0fd6f68`,
respectively. The producer SHA-256 values, registration seals, exact phase
costs and all per-trial records are in the machine-readable report.

This is a **target-free preparation-stage result**. Candidate, workload and
run IDs remain null because no complete one-target pipeline is registered.
There is no SAT factor-log table, no new public-target recovery, no paired
rho solve and no online IC speedup. The three historical confirmation sets
and both consumed registrations remain closed. A revised CMS budget or
encoding needs a new source-bound registration and another fixed natural
panel; the 107 missed feasible rows cannot be filled from the geometric
oracle or the F5 arm. The F5 arm may now bind these audited logs to a new,
separate one-target registration.
