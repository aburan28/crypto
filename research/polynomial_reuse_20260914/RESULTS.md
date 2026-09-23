# Polynomial reuse: measured results

Classification: **engineering experiment; no end-to-end speedup established**. The affine template is correct and reusable across targets. The current serialized cache path adds enough overhead to erase its stage savings in the observed full-run medians. Full parameterized Boolean bases are not promoted. All cache layers remain opt-in.

**2026-09-23, run-002:** a re-run on current `main` reproduces this verdict. No cache mode pays end to end, and the reduction cache is no longer reached by the default Gröbner engine. See [Run 002](#run-002-current-main-against-the-same-reference).

## Frozen scope and checks

Run `results/run-001` contains all **159 process records**: 63 stage processes and 96 full-DLP attempts. Of the stage processes, 54 completed and 9 parameterized-basis jobs hit the declared 30-second limit. All **1,080 completed target/repetition/variant checks** agree with exhaustive ANF truth and both SAT and F4 solution extraction. This represents 152 distinct field/layout/target systems for baseline/template, with 56 also completed by the parameterized route; repetitions are not independent instances.

All 96 DLP outcomes and shared logical counters match the unchanged reference. Eight distinct n=9 fixtures (two solvers, four seeds) verify in all six modes: 48 verified executions. The eight n=7 fixtures remain incomplete in all modes: 48 incomplete executions. Their default factor base has one point and zero usable columns; none is treated as a solved DLP or silently discarded.

Redis cross-process replay recorded **175 hits**, with **zero cache errors or corrupt values** in the matched suite. This was a dedicated local Redis 7.4.2, not AWS ElastiCache. Separate tests cover corruption/recomputation, layer/key separation, bounded retention, preserved solver node budgets, all-target template equivalence, specialization, and outage fallback: **21 tests passed**, plus an unavailable-Redis DLP comparison that preserved the answer and recorded one connection failure. Feature-enabled release build and no-Redis feature check passed. Existing unrelated compiler warnings remain.

## Stage comparison

The table reports medians of three repetitions. One unit: **milliseconds for field/template/basis setup + encoding/instantiation + F4 solving** across the fixed targets. Candidate/reference ratios below 1 are lower observed stage time. Exhaustive certification, SAT cross-checking and the extra root-matrix counter probe are timed/recorded separately and excluded from this F4-stage column. These are shared-host diagnostic timings without a paired confidence claim.

| n | ell | m | Variant | F4 stage ms | Candidate/reference | Template bytes | Parameter basis bytes | Correct/complete |
|---:|---:|---:|---|---:|---:|---:|---:|---|
| 3 | 1 | 2 | baseline | 0.1013 | 1.000 | 0 | 0 | yes |
| 3 | 1 | 2 | template | 0.0647 | 0.639 | 505 | 0 | yes |
| 3 | 1 | 2 | parameterized | 0.1071 | 1.058 | 505 | 148 | yes |
| 3 | 2 | 2 | baseline | 0.2646 | 1.000 | 0 | 0 | yes |
| 3 | 2 | 2 | template | 0.2233 | 0.844 | 803 | 0 | yes |
| 3 | 2 | 2 | parameterized | 32.9658 | 124.602 | 803 | 1360 | yes |
| 5 | 2 | 2 | baseline | 0.8004 | 1.000 | 0 | 0 | yes |
| 5 | 2 | 2 | template | 0.6717 | 0.839 | 1481 | 0 | yes |
| 5 | 2 | 2 | parameterized | 3838.1820 | 4795.114 | 1481 | 3397 | yes |
| 7 | 3 | 2 | baseline | 2.2725 | 1.000 | 0 | 0 | yes |
| 7 | 3 | 2 | template | 2.0717 | 0.912 | 3226 | 0 | yes |
| 7 | 3 | 2 | parameterized | — | — | — | — | 3 × 30 s timeout |
| 9 | 3 | 2 | baseline | 1.5643 | 1.000 | 0 | 0 | yes |
| 9 | 3 | 2 | template | 1.0716 | 0.685 | 4431 | 0 | yes |
| 9 | 3 | 2 | parameterized | — | — | — | — | 3 × 30 s timeout |
| 3 | 1 | 3 | baseline | 0.3671 | 1.000 | 0 | 0 | yes |
| 3 | 1 | 3 | template | 0.3667 | 0.999 | 988 | 0 | yes |
| 3 | 1 | 3 | parameterized | 1.1009 | 2.999 | 988 | 584 | yes |
| 5 | 2 | 3 | baseline | 66.1723 | 1.000 | 0 | 0 | yes |
| 5 | 2 | 3 | template | 65.3199 | 0.987 | 3982 | 0 | yes |
| 5 | 2 | 3 | parameterized | — | — | — | — | 3 × 30 s timeout |

The reusable template stores roughly **505–4,431 serialized bytes** in these cases. Direct template use has F4-stage candidate/reference ratios of about **0.64–1.00** and SAT-stage ratios of **0.51–0.96**. This experiment retains the template as a Rust object; the Redis-backed implementation also pays key construction, serialization/deserialization, allocation and network costs. That distinction matters in the full-run table below.

The n=5, ell=2, m=2 parameterized basis occupies only 3,397 bytes but takes about **3.84 seconds** to prepare. Its setup-plus-F4-stage cost is roughly **4,795×** the unchanged stage in this batch. Storage capacity is not the limiting factor here: constructing the reusable basis is. The n=7/ell=3/m=2, n=9/ell=3/m=2, and n=5/ell=2/m=3 parameterized jobs each time out on all three repetitions. These findings concern this encoding, variable order and existing Boolean Buchberger engine; they do not rule out other parameterized algorithms.

## Complete pipeline comparison

One unit: **end-to-end elapsed milliseconds summed over the four n=9 seeds for each solver**. Each ratio is the sum divided by the corresponding unchanged-reference sum. Every listed execution verifies the expected scalar and point identity. n=7 incomplete attempts remain in raw files and summary but are excluded from this equal-verified-workload timing table. Timing includes cache setup/access, encoding, relation work, linear algebra and verification inside `ic`; process orchestration is separately recorded.

| Solver | Mode | Total elapsed ms | Candidate/reference | Verified |
|---|---|---:|---:|---:|
| groebner | reference | 24.614 | 1.000 | 4/4 |
| groebner | off | 24.232 | 0.984 | 4/4 |
| groebner | preprocess-local | 25.723 | 1.045 | 4/4 |
| groebner | both-local | 32.901 | 1.337 | 4/4 |
| groebner | both-redis-cold | 46.544 | 1.891 | 4/4 |
| groebner | both-redis-warm | 37.134 | 1.509 | 4/4 |
| sat | reference | 10.612 | 1.000 | 4/4 |
| sat | off | 12.564 | 1.184 | 4/4 |
| sat | preprocess-local | 12.516 | 1.179 | 4/4 |
| sat | both-local | 11.840 | 1.116 | 4/4 |
| sat | both-redis-cold | 13.477 | 1.270 | 4/4 |
| sat | both-redis-warm | 12.608 | 1.188 | 4/4 |

Redis warm mode runs the identical workload in a new process. Its hits demonstrate cross-worker transport/reuse, not additional independent relations. Cold and replay results are both retained; cache population is never treated as free. The uncached candidate is also measured to distinguish cache effects from unrelated build/runtime noise.

## Boundaries and decision

For a fixed workload, if solver/preprocessing fraction f has fraction h eliminated, the ideal speedup is at most 1/(1-fh), before cache overhead and precomputation. Exact replay adds no matrix rank. The experiment has no calibrated conversion from ANF XORs, F4 elimination XORs, SAT conflicts and group/linear-algebra operations into one common unit. Therefore **total calibrated operations, S = total operations / sqrt(N), ratios to measured rho and a derived attack floor, and scaling exponent remain null for every variant**. No generic-group improvement, asymptotic claim, or full-DLP speedup follows.

Next engineering step: retain the decoded preprocessing object per worker to avoid repeatedly deserializing a small reusable template, then repeat the same matched suite with larger verified workloads. For parameterized algebra, test target-independent chain-prefix reductions or limited-degree reusable consequences before attempting a full basis over all target bits. The current full-parameter basis route failed its promotion criterion.

## Run 002: current main against the same reference

Re-run on 2026-09-23 with `run.py`, `compare.py` and `contract.json` **unchanged**. It is kept beside run-001 rather than replacing it. The setup:

- **Candidate:** `main` at `6db529ea`, built with `--features redis-cache`. The run includes the decoded-value local layer (PR #531) and the field-closed Gröbner engine (PR #626).
- **Reference:** a fresh build of `01fa842`. Its binary hash differs from run-001's (`3ea7c7ec…` against `a0b9f9e4…`), but its outputs are identical in all 16 reference runs: status, result and every logical counter. It is the same reference.
- **Host:** Redis 7.0.15; 4 cores; `RAYON_NUM_THREADS=1` as the script sets.
- **Timing:** one sample per cell. The load average was 4.18 at the start (the tail of the two builds) and 1.05 at the end.

### What `compare.py` says, and why

**Admission fails**, with `all_dlp_outcomes_match_reference: false`. Everything else is healthy:
- all 159 processes are recorded;
- every stage answer matches exhaustive truth;
- there are 0 Redis errors;
- there are 16 warm new-process hits;
- the same 9 parameterized jobs time out as in run-001.

The failing rows are **20 of 96: every n = 9 Gröbner solve, in every mode, including `off`**.
- **The results are identical and verified.** Only the path differs: seed 101 takes 3 trials and 28 F4 reductions where the reference takes 2 and 24.
- **The cause is the solver, not the cache.** #612 made `InheritedF4` the default engine, and it splits on the smallest free variable. At `01fa842` the default was `MatrixF4 { max_degree: 3 }`. So the candidate's solver is no longer the reference's.
- **#626 is not the cause.** `main` at `4428148f`, just before #626, gives the same counters on all four seeds.

The admission test assumes the candidate solver equals the reference solver. That no longer holds, so the reference column is a correctness anchor (same answers, all verified) and no longer a counter-level match.

What that test protects is that **the cache must not change behaviour**. With the reference behind, the form of that invariant that still applies is: every cache mode must match `off` in the same binary. It holds in **64 of 64** mode runs, in status, result and every logical counter.

### Complete pipeline comparison

The unit is the same as above: end-to-end elapsed milliseconds summed over the four n = 9 seeds. **The ratio to `off` is the primary reading**, because `off` is the same binary with the cache disabled. The ratio to the reference also carries every solver change since `01fa842`.

| Solver | Mode | Total elapsed ms | / off | / reference | Verified | run-001 / off |
|---|---|---:|---:|---:|---:|---:|
| groebner | reference | 30.122 | 2.722 | 1.000 | 4/4 | 1.016 |
| groebner | off | 11.067 | 1.000 | 0.367 | 4/4 | 1.000 |
| groebner | preprocess-local | 12.082 | 1.092 | 0.401 | 4/4 | 1.062 |
| groebner | both-local | 12.279 | 1.110 | 0.408 | 4/4 | 1.358 |
| groebner | both-redis-cold | 14.906 | 1.347 | 0.495 | 4/4 | 1.921 |
| groebner | both-redis-warm | 14.458 | 1.306 | 0.480 | 4/4 | 1.532 |
| sat | reference | 14.030 | 1.293 | 1.000 | 4/4 | 0.845 |
| sat | off | 10.847 | 1.000 | 0.773 | 4/4 | 1.000 |
| sat | preprocess-local | 11.870 | 1.094 | 0.846 | 4/4 | 0.996 |
| sat | both-local | 11.872 | 1.094 | 0.846 | 4/4 | 0.942 |
| sat | both-redis-cold | 14.204 | 1.309 | 1.012 | 4/4 | 1.073 |
| sat | both-redis-warm | 13.017 | 1.200 | 0.928 | 4/4 | 1.004 |

These single-sample differences are backed by the cache's own counters, not left to noise:
- **Preprocessing layer.** It records 1.03–1.11 ms of overhead over the four solves, 8.6–9.4% of elapsed time. That accounts for the measured 1.09×.
- **Why it cannot pay.** Each `ic run` process makes 2.5–3.5 lookups: one cold miss (key construction, admission encoding and hashing) and 1.5–2.5 hits. A template that builds in microseconds, reused once or twice more, cannot repay a cold miss.
- **Redis modes.** They record 3.0–3.7 ms (21–25%) of overhead.
- **The #531 fix.** Holding decoded values cut the local Gröbner penalty (`both-local` from 1.358 to 1.110 of `off`) but did not reverse it.

**The exact-reduction layer is no longer reached on the default engine.** It records **0 lookups** in run-002:
- `InheritedF4` reduces through `reduce_inherited`, which does not consult the memo (`koblitz_groebner.rs`, `solve_rec`). So `IC_REDUCTION_CACHE` has no effect on the default Gröbner path, and the `both-*` modes now differ from `preprocess-local` only by noise.
- This is also why warm new-process Redis hits fell from 175 to 16: the reduction layer's hits are gone.
- Its absence is not a loss. In run-001 that layer recorded 6 local hits in 109 lookups, and 22–48% of elapsed time as overhead.

The Gröbner `off / reference` ratio of 0.367 comes from the engine change, not from caching. The two sides also take different trial paths, so this is not the same work done faster. It is recorded here and is not claimed.

### Stage comparison, run-002

The same checks passed. Cold-F4 candidate/baseline for the **template** is 0.75–1.05 (run-001: 0.64–1.00). The template still saves at stage level, where it is held as a plain Rust object with no cache machinery.

The **parameterized** route measures 22.6× and 776× at n = 3 and 5 with ell = 2 (run-001: 124.6× and 4,795×). The chain criterion made the engine faster since then, and the route is still nowhere near promotion. Its full-basis cache was removed in #626. The same three parameterized cells time out on every repetition.

### Decision after run-002

**Unchanged.** No cache mode improves end-to-end elapsed time at the frozen sizes. Every mode is slower than `off` in the same binary, and the cache's own overhead counters account for the difference. The preprocessing cache stays opt-in and should stay off by default.

What would change this verdict is a workload with many lookups per process: long relation searches at a degree where the template build is a real fraction of each trial. The frozen suite's degrees, 7 and 9, cannot reach that regime.

## Evidence

- [Frozen contract](contract.json)
- [Machine-readable comparison](results/run-001/summary.json)
- [All process statuses and commands](results/run-001/processes.json)
- [Source and binary hashes](results/run-001/hashes.json)
- [Test outcomes](validation/tests.json)
- [Outage comparison](validation/outage.json)
- Run 002: [comparison](results/run-002/summary.json), [process statuses](results/run-002/processes.json), [hashes](results/run-002/hashes.json)

`results/run-001/*.stdout` and `*.stderr` preserve every raw run, including timed-out basis jobs. `validation/dependencies.lock.txt` preserves the dependency resolution. Benchmarks were executed locally; no AWS resources were provisioned.
