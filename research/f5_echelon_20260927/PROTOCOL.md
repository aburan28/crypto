# Matrix-F5 echelon output experiment

## Caller and output contract

The only in-repository calls to `matrix_f5_f2` and `matrix_f5_f2_timed` are the F5 benchmark and module tests. The existing public functions continue returning the same RREF. This experiment adds an explicit `F5OutputForm::Echelon` choice, which skips above-pivot reduction and returns a different polynomial basis of the same row space. The benchmark emits both a raw `rows_fp` and a canonical `row_space_fp` computed by reducing the returned rows again **outside** the timed F5 call. Only the canonical fingerprint is required to match. The changed raw fingerprint is retained in every receipt. The F5 criterion, row construction, rank and pruning counters must match exactly.

## Frozen paired measurement

Use the unchanged seven-case matrix-F5 suite in `examples/f4_f2_bench.rs`, primary `f5_n24_m24_d4`, with frozen seed XOR zero and holdouts `badc0de1` and `5eed2026`. The benchmark selects the explicit output form through `KIC_F5_ECHELON=0/1`. Set `KIC_GF2_WORD_BATCH=0` and `KIC_GF2_DEFER_ABOVE=0` in both arms so this result is independent of PR #883 and the negative deferred-reverse experiment. Build one release binary. Pin one allowed CPU, set `RAYON_NUM_THREADS=1`, warm each arm once per workload, collect five A/A reference pairs and five alternating-order A/B pairs. Preserve every call, process status, full output, binary/source hashes, host/load and exact five-pair bootstrap intervals.

`reduce_ms` covers elimination; `wall_ms` covers the whole matrix-F5 step, including unpacking, but excludes process launch and fixture construction. Canonical row-space fingerprinting is outside both intervals. This is an internal point-decomposition solver-stage measurement, not a one-target IC or DLP speedup.

## Promotion gate

Require exact canonical row-space fingerprints and the same rank, pruning, criterion work and matrix-kernel row-space tests. The frozen primary elimination median reference/candidate ratio must reach 1.20 with a 95% paired lower bound above one; the complete-call ratio and lower bound must exceed one. Both holdout primary ratios must clear their A/A ranges, and no smaller complete-call median may regress below its own A/A range. If one-thread gates pass, run the paired four-thread control before considering a default output change. Otherwise keep echelon as an explicit option and preserve the negative or inconclusive receipt.
