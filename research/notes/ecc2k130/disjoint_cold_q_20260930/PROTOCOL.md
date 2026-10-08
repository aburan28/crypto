# Disjoint-Q cold confirmation for the compact W64 index calculus

Status: preregistered before generating any Q for this panel or measuring an arm. This is a six-cell, native CPU confirmation of the current compact producer against the matched strong rho reference. It tests target variation missing from the same-Q cold panel; it does not test an isogeny-descended factor base, a new PDP solver, or the ECC2K-130 challenge instance.

## Hypothesis and fixed boundary

The current W64 compact S3 producer reaches full rank and recovers all targets in the [same-Q cold panel](../compact_ir_cold_gap_20260930/RESULT.md), but its complete-process CPU exceeds matched signed-Frobenius rho in every tested cell. The hypothesis is that the smallest observed batch gap, n37/L1024, remains greater than one on new public Q batches when variation across targets is included. A disjoint-Q result can overturn that hypothesis, but one favorable block cannot: the decision uses all preregistered blocks and costs.

Both arms use the exact twenty-file source in [`compact_s3_prefilter_20260930/FROZEN.json`](../compact_s3_prefilter_20260930/FROZEN.json), whose SHA-256 is `3e9f67cc2cd6de5a8458badb3525983d561d9c3449c05118819fb424c5093b2b`. The freeze pins the compact S3 source at `702a0a05709bc14bc10bafbb08edb2a2f5f86794e970d84c317da3bc65cdcf38`, rho v3 at `98b6e8a27d821ebfc7ae716410184dffe89aadf27dfe36bcac8ce834c39cf04c`, and Cargo.lock at `7d671f48c2da93f133d98802e80f858d1d9ea3b86996f7037f758990e1566627`. Materialize it with the pinned `compact_frozen_source_replay_20260929/materialize.py`; reject any file or binary mismatch. The child-process measurement helper is `compact_ir_ledger_20260930/run_panel.py`, SHA-256 `bb57f12a4f57b6921674af026a5993c661751bdcbc47d0984b6f32af64b014b8`. Archive the full source manifest, toolchain, compiler flags, and binary hashes with the outcome.

| Cell | Useful K | Q per block | S3 prefilter | Blocks | Seed for block b |
|:--|--:|--:|:--|--:|--:|
| n37/L1 | 7 | 1 | off | 20 | 2026093091000 + b |
| n37/L1024 | 42 | 1,024 | off | 5 | 2026093091100 + b |
| n41/L1 | 85 | 1 | off | 5 | 2026093091200 + b |
| n41/L1024 | 255 | 1,024 | blocked | 5 | 2026093091300 + b |
| n53/L1 | 220 | 1 | off | 5 | 2026093091400 + b |
| n53/L1024 | 440 | 1,024 | blocked | 5 | 2026093091500 + b |

The Q generator is the pinned rho v3 with `KIC_RHO_GENERATE_ONLY=1`, `KIC_RHO_BATCH_CORPUS=compact-disjoint-cold-n{n}-L{L}-b{b:02d}-20260930-v1`, and command `rho n 0 signed_frobenius L seed`. Each of the 45 blocks has a distinct seed and a separate point-only file. The fixture scalar is kept in a separate verifier-only label file and never enters an IC or rho command. Within a block, IC-A, rho, and IC-B use the **same** Q file. Across blocks, cells, and prior published point-only files, Q must be disjoint at a given n. This fixes 15,390 target points and 45 paired blocks before measurement.

The prior inventory is exactly the 31 committed `research/notes/ecc2k130/**/*.points.jsonl` files present before this panel: 19,468 rows. For the sorted list of records `{path,sha256,rows}`, serialized as canonical JSON with sorted keys and no extra spaces, its SHA-256 is `30a0b11482bad46c6efe994278001941a2e59c84cf10d366ec433b785637c9cd`. The preparation script must assert this inventory and reject any Q collision, malformed row, wrong curve, wrong subgroup, inconsistent generator, or `[d]G != Q`. A failure invalidates the corpus; it must be recorded and a new protocol frozen before a replacement seed is used. Commit the point and label files plus hashes in `FROZEN.json` before any timed arm. No outcome-dependent selection or resampling is allowed.

## Cold execution and independent replay

Run three fresh processes per block, rotating the arm order `ic_a,rho,ic_b` by `b mod 3`. The compact arms use `construct:n:0:K`, rank seed 7, S3 W64, the fixed prefilter, and one worker. Rho uses 32 walks, DP bits 4, signed Frobenius and normal-basis canonicalization. The complete process charges startup, curve and basis setup, factor-base/index construction, full-rank collection, linear algebra, all target descents, group checks, and output. Neither arm receives fixture scalars or a prebuilt index.

Run on Linux x86-64 with one pinned CPU and a reserved-core isolation monitor. Record `wait4` child user+system CPU, wall, max RSS and observed peak RSS, exact command/environment, raw stdout/stderr, base/rank/target traces, source/input/binary hashes, host CPU details and isolation samples. Each arm has a 900-second timeout and 5-GiB address-space/RSS limit; each cell job has a 90-minute cap. Stop a cell on its first failed child and preserve every completed child and the failure. Do not omit a block or selectively rerun a cell under this protocol.

An independent Python verifier must validate the frozen input hashes, every generated `[d]G=Q`, every compact base point, four-point relation, rank row, representative log, target scalar, and every rho scalar against the point-only Q file. It must reject altered raw bytes, swapped Q files, incomplete ranks, missing targets, changed commands, source mismatch and label leakage. Archive a second-host replay receipt before treating any timing as eligible.

## Analysis and decisions

For each block, calculate `IC_CPU = sqrt(CPU_ic_a * CPU_ic_b)`, `R = IC_CPU / CPU_rho`, and drift `D = CPU_ic_b / CPU_ic_a`. Report all child costs and the median and two-sided 95% log-ratio t interval of R and D for each of the six cells. The paired interval now reflects deterministic fresh-target variation as well as process noise. A cell is timing eligible only if all fixed arms succeed, all scalars and ranks independently replay, the D median lies in [0.9,1.1], the D interval includes one, and the isolation monitor reports zero contention. Any missing gate censors the entire cell; retain its raw evidence and do not impute a speedup.

An eligible R interval wholly below one is a native-CPU crossover **at that n/L/K and host**. It triggers a separately frozen second-host confirmation and a common all-phase operation-equivalent `S=operations/sqrt(r)` calibration before a method-level speedup claim. An interval wholly above one is a quantitative no-go for this fixed compact policy at that cell; report the necessary reduction `1 - 1/R` using the median R. An interval straddling one is unresolved and is not extended adaptively. If any cell is censored, report the cause and keep its ratio unset. The independent instruction ledger remains a separate metric, not a conversion to CPU or generic group work.

Update the canonical scoreboard and boundary ledger only with eligible full-process evidence. This panel cannot establish n83/n131 transfer, isogeny-descendant PDP benefit, a challenge logarithm, or a generic-group crossover; those require their separately frozen policies and costs. The n131 m10 capacity dispatch remains subject to the independent-review and explicit-label gate in [PR #937](https://github.com/aburan28/crypto/pull/937).
