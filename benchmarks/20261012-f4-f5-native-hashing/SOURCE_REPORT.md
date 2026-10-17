# Boolean F5 hashing: measured follow-up

The accepted native comparison records **3.71% fewer instruction references for complete F5 calls** and **45.14% fewer for standalone criterion construction** on eleven frozen Boolean systems. F4's control is neutral. Every input/output digest, rank, pruning report, criterion membership digest and word-XOR counter matched. No primary or paired-control workload increased allocation requests, cumulative requested bytes or tracked peak live bytes.

| Arm | Baseline median Ir / 11-system suite | Candidate median Ir / suite | Paired instruction reduction | 95% bootstrap interval | Native verdict |
| --- | ---: | ---: | ---: | --- | --- |
| Complete F4 control | 98,506,024 | 98,471,390 | 0.035% | 0.029–0.042% | Neutral: below the 2% gate |
| Complete F5 | 95,971,529 | 92,406,706 | 3.715% | 3.714–3.726% | Fewer instructions |
| F5 criterion | 6,128,134 | 3,362,140 | 45.136% | 45.121–45.155% | Fewer instructions |

These are Callgrind instruction counts from an instrumented native driver, including allocator accounting, with merged worker coverage. They are not elapsed-time measurements or hardware retired-instruction counters. Primary profiles contain three calls per system and arm; verdicts use five paired A/B/A′ process rounds, each covering one call on every system. The native comparison applies the repository's log-ratio, A/A-noise, 95% bootstrap and 2% effect gate. Every workload verdict is retained; none crossed the instruction-regression gate.

## Implementation and workloads

The [baseline profile](records/baseline-profile.json) attributes 8,509,068 of 18,228,487 criterion self instructions (**46.68%**) to standard hash functions across 33 primary criterion calls. It was captured before production edits. F5 now uses the existing mixed-mask `MaskHasher` for column maps and pruned-monomial sets, and `FxHasher` for dense temporary pivot indices. High-bit mixing avoids clustering sparse masks with many zero low bits. Pruning rules, generator/signature order, row reduction, output form and size limits are unchanged.

Production baseline: `baa624218809b821400cfe33155c6bacb542e9a8`. Common driver/protocol: `1bfb7acadbe6a8dea034ef5707d7dbdcd2611ded`. Candidate production change: `0505fb45ae98b0fe4cb3f86fb17305b67f01d676`. F4 is the control for this iteration; [the earlier row-building study](../f4_f5_row_build_20261011/REPORT.md) remains unchanged.

Eight systems retain the earlier quadratic fixtures and two seeds. Three fresh-seed holdouts add 12-variable degree-three/four cases and a case supported on bits 52–63 with six multiplier variables. Definitions, intervals, budgets and acceptance conditions are frozen in [PROTOCOL.md](PROTOCOL.md). This is separate from the historical n=17/n=23 matrix/solver workload and general coefficient-field algorithms.

## Memory and reports

The allocator records live requested bytes from process initialization, live bytes at call entry and the absolute interval peak. Native preparation also derives peak additional live bytes above that entry, so a smaller retained warm criterion is not attributed to F4's own heap needs. Criterion additional peaks decrease by 96–192 bytes per call, matching smaller stored set headers; request counts are unchanged. Absolute peaks, entry values and cumulative requested bytes remain in the records. Process RSS is unmeasured.

The [combined three-page PDF](reports/run-20261012T0215Z/report.pdf) contains report text and the three SVG groups, with editable Typst:

- [Complete-call instruction graph](reports/run-20261012T0215Z/complete-instructions.svg) and [PDF](reports/run-20261012T0215Z/complete-instructions.pdf)
- [Criterion instruction graph](reports/run-20261012T0215Z/criterion-instructions.svg) and [PDF](reports/run-20261012T0215Z/criterion-instructions.pdf)
- [Peak live allocation graph](reports/run-20261012T0215Z/peak-live.svg) and [PDF](reports/run-20261012T0215Z/peak-live.pdf)

Instruction-history points are measured production revisions; controls are not extra revisions. Memory points are native medians of three primary calls. The performance index links these groups separately from legacy CI diagnostics. Historical quantitative groups, the allocation index, aggregate index, IC/rho scoreboard and progress timeline retain their values: this study has a new eleven-system/profiler contract and supplies no ECDLP measurement or curve fact.

## Reproduction and retained attempts

The [native driver](../../examples/f4_f5_row_build_bench.rs) executes and verifies calls. The [native comparison](../../examples/f4_f5_hash_compare.rs) checks full frozen coverage, successful exits, exact accounting, shared inputs, request dumps and instruction self-cost totals. [comparison.json](records/comparison.json) preserves raw-source BLAKE3 hashes, measured values, native medians, confidence intervals and verdicts. All raw JSONL and Callgrind dumps remain in `records/`; every measurement process completed within its budget.

Callgrind 3.24's incremental `summary` wrapped after zero-stats requests. The parser uses positive file `totals` and independently sums self costs. Original wrapped fields and rejected inspection remain unchanged. A macro compilation error and initial Clippy findings are retained. Final Clippy cleanup reproduced the comparison byte for byte without changing the bootstrap contract.

Build the locked root Cargo examples. Extract the recorded Valgrind 3.24 Debian package outside the repository and set `VALGRIND_LIB` to its `usr/libexec/valgrind`. [run-profiles.sh](run-profiles.sh) is thin command orchestration; binaries use equally long `a/bench` and `b/bench` paths, and existing run directories are refused. Read-only comparison and rendering:

```sh
cargo run --locked --release --example f4_f5_hash_compare -- research/f4_f5_hash_20261012/records
/tmp/f4-f5-report-venv/bin/python research/f4_f5_hash_20261012/render.py --out /tmp/f5-hash-render
```

Rust performs execution and result math; Python lays out recorded values and renders Typst. Dependencies are pinned in `scripts/perf/groebner-render-requirements.txt`. Compiler/host manifests, source/binary hashes, commands and statuses are retained. The first baseline process's exact start timestamp was not captured; its profile-file timestamps are observations, not an invented start time. The report ID uses their UTC minute.

## Validation and remaining work

All **18 focused F5 tests** passed, including an independent high-bit embedding check of criterion membership/accounting and canonical full-versus-pruned row spaces. The native comparison's **five tests** passed. Both changed examples passed direct Clippy checks; the comparison check denies warnings. Native profiling passed all eleven systems, primary repetitions, baseline A/A and five rotated A/B/A′ rounds.

The broader release-library check **passed: 3,906 tests, zero failures and 110 ignored**, completing in 847.81 seconds within the 20-minute budget (exit 0). Its exact catalog and named completion sets are retained. The previously unfinished isogeny test passed unmodified; older timeout records remain historical evidence rather than being erased. The run used three test threads and a network-enabled sandbox; the standalone measurement driver still fixes one Rayon worker. This result establishes successful verification in this run, not a diagnosed fix for the prior timeout.

Controlled runtime, fresh-process RSS, a larger preregistered scaling matrix and an algorithm selector remain pending. This host lacks `/proc/pressure/cpu`, so controlled timing eligibility is unearned. Every page of the four study PDFs and the refreshed roadmap PDF was inspected. Study SVGs, Typst sources and PDFs reproduced byte for byte on an independent render; the native comparison reproduced byte for byte too.

Conductor is unavailable; attempted reservations and scope expansions are retained. The accepted instruction result does not close pending timing, scaling or independent library-validation work.
