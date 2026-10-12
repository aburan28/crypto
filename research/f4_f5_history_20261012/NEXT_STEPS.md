# F4/F5 follow-up roadmap

2026-10-12 — **Original planning snapshot; current follow-up execution status is below.**

Execution update: the [native hashing study](../f4_f5_hash_20261012/REPORT.md) completed worker-inclusive instruction profiling, tracked live-allocation measurements and a measured F5 optimization, including three fresh-seed holdouts. Complete F5 calls use 3.71% fewer instructions and criterion construction 45.14% fewer under its native controls. Full release-library verification passed with 3,906 tests and 110 ignored; the previously unfinished test passed unmodified. Controlled runtime/RSS and the larger scaling matrix remain pending. The diagram and PDF retain the original route as a labeled planning snapshot.

The next question is whether the verified allocation reductions improve complete-call runtime and peak memory. The [native allocation study](../f4_f5_row_build_20261011/REPORT.md) measures cumulative allocation traffic, while the [history report](REPORT.md) retains diagnostic CI timings without isolation eligibility or stored spread. Neither establishes a runtime or peak-memory improvement.

The [editable workflow diagram](next-steps.svg) and [one-page PDF](next-steps.pdf) describe the proposed sequence. The PDF's source is [next-steps.typ](next-steps.typ).

## 1. Establish a controlled runtime and memory comparison

- First compare the unoptimized production revision `e957e48c2b8cf1f84a7fc44a1d673103964a789d` with the optimized revision `2e3dffda681a4771be511653180da8354de3bce3`. Use the same frozen driver, prerequisite repairs, compiler, flags, host and inputs in both builds. Keep a later current-main replication separate so unrelated intervening changes cannot be attributed to this optimization.
- Retain the allocation study's eight systems and two seeds, then add the shared matrix-step and solver workloads from [the group registry](../../docs/perf/groebner-groups.json). Measure complete construction, reduction and readback. Record F5 criterion construction separately as a phase; its cost stays inside the complete F5 call.
- Implement the harness, comparison math and isolation envelope natively. Preserve [the repository's isolation contract](../../AGENTS.md#10-performance-work-starts-from-a-clean-recorded-baseline), including an exclusive lock, core reservation, contention/PSI diagnostics and explicit `RAYON_NUM_THREADS=1`. A pinned process alone is insufficient. Freeze the protocol, binaries and budgets before measurement; do not run builds or tests during a timed stage.
- Run at least five paired rounds of baseline A, candidate B and repeat baseline A′, rotate arm order and seed workload shuffling. Keep initialization, warmup and verification outside the timed interval. Time without the allocation-counting instrumentation. Preserve every sample, contended run, failure and timeout.
- Admit runtime verdicts only using [the existing noise gate](../../docs/perf/PERFORMANCE_INDEX.md#the-formula): a 95% bootstrap interval excluding zero and an effect larger than both twice the A/A noise estimate and 2%; otherwise report neutral. Implement the documented formulas in native code and freeze the bootstrap seed and count. Contended timings remain diagnostics.
- Measure memory in separate instrumented runs: peak live allocated bytes over the same complete-call interval, recording the live bytes at entry too. Also retain fresh-process peak RSS with setup costs explicitly included. Cumulative requested bytes, live heap and RSS are distinct measures; do not subtract process high-water marks to invent a phase peak.

Gate: identical results and accounting for each algorithm's baseline/candidate pair, an eligible timing envelope, raw paired samples and explicit memory intervals. This stage may report neutral results or regressions.

## 2. Profile before selecting the next code change

Use separate, worker-inclusive native profiling runs to split symbolic preprocessing, row packing, F5 criteria, pivot search, elimination and readback. Retain rows, columns, matrix density/fill-in, pruning, word-XOR counts and phase costs alongside the input hashes.

Choose one change from the measured bottleneck. Packed elimination/cache layout or dispatched SIMD is a hypothesis if elimination dominates; signature/criterion indexing is a hypothesis if F5 bookkeeping dominates. Keep portable fallbacks and compare one change at a time. Do not assume that reducing temporary allocations makes either algorithm faster.

Gate: a recorded hotspot, an unchanged correctness contract and a matched baseline/candidate comparison. Instruction profiles must include worker threads; historical incomplete profiler coverage is not a substitute.

## 3. Find the F4/F5 crossover on held-out workloads

Preregister a bounded scaling matrix over variables, equations, degree and density, with unseen seeds and time/memory limits. Keep the shared workload definitions fixed and preserve OOMs and timeouts. Test one worker first, then fixed supported multicore configurations. Report the hardware class actually measured.

For F4 versus F5, compare mathematical results using independent row-space/equivalence checks and verified solver roots. Their pruning and operation counters can differ and must remain separately reported; do not require cross-algorithm equality of performance fingerprints. Baseline/candidate equality within each algorithm remains mandatory.

Gate: reproducible complete-call results and held-out crossover evidence. Consider a workload-based algorithm selector only after that evidence exists. General coefficient fields are a separate workstream; the current changes concern Boolean F4/F5.

## 4. Publish accepted measurements together

Deliver each completed implementation or experiment in a linked follow-on PR, with frozen source/input hashes, native records, correctness checks, controls and failures. Extend graph schemas to retain isolation eligibility, spread and worker coverage for the new native measurements. Keep those series separate from legacy CI diagnostics. Regenerate source data, SVGs and PDFs together; allocation controls remain controls rather than extra production revisions.

No quantitative graph changes in this planning update: the four current study groups, `docs/perf/PERFORMANCE_INDEX.md`, the aggregate history, `docs/index-calculus-scoreboard.html` and `docs/ic/progress-timeline.json` retain their existing values. There is no new measurement or curve fact to plot.

## Independent validation work

Original proposal: investigate `cryptanalysis::jv_isogeny_walk::tests::end_to_end_from_a_non_weak_curve_recovers_and_verifies` with a bounded reproducer and diagnose its stalled work. The [prior validation record](records/validation.json) was partial: 3,904 distinct library tests passed, 110 were ignored and this test timed out. In the later hashing study, all library tests completed and this test passed unmodified. The earlier timeout has not been diagnosed; investigate it independently if it recurs rather than skipping the test or weakening assertions. Run library validation separately from timed benchmarks.

| Deliverable | Status |
| --- | --- |
| Source-linked roadmap, explanatory SVG and inspected PDF | Complete planning deliverable |
| Controlled native runtime and peak-memory campaign | Tracked live allocations measured in hashing study; controlled runtime and fresh-process RSS pending |
| Worker-inclusive profiling and next optimization | Completed for criterion hashing in the linked native study; other hotspots remain follow-ons |
| Held-out scaling, crossover and optional selector | Three native holdouts measured; larger scaling matrix, crossover and selector remain pending |
| Complete library verification and investigate prior timeout | Full verification passed in hashing study; test unmodified, prior timeout diagnosis remains open |

Conductor was unavailable; the original planning reservation is preserved in [the planning scope record](records/next-steps-conductor.log). The original update changed documentation and its diagram/PDF only. Later execution and validation are recorded in the linked hashing study; the prior library-validation record remains historical evidence.
