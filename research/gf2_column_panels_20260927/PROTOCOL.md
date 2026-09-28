# Cache-sized column panels for shared GF(2) Macaulay elimination

## Hypothesis and implementation

On the frozen matrix-F5 `n24_m24_d4` step, the shared kernel clears the full matrix after each pivot block. The earlier two-block, one-word experiment in PR #883 was exact but approximately neutral and slowed smaller cases. This distinct candidate keeps a two-word pivot panel current, snapshots each block's Gray-code combinations and row selection patterns, updates selected pivot rows' trailing suffixes on demand, and replays the panel's deferred trailing updates across cache-sized column ranges. It groups several pivot blocks per full matrix pass. The original per-block path is the same-binary reference with `KIC_GF2_COLUMN_PANEL_WORDS=0`; the candidate is `=2`. Neither arm uses the negative deferred-above mode or F5 echelon output.

## Frozen paired workload

Use `examples/f4_f2_bench.rs`, one repetition of all seven matrix-F5 cells per process, primary `f5_n24_m24_d4`. Freeze seed XOR zero and holdout XORs `badc0de1` and `5eed2026`. Pin one Linux CPU, set `RAYON_NUM_THREADS=1`, warm each arm once per workload, then run five A/A reference pairs and five alternating-order A/B pairs. Record every call, failure, timeout, full output, source and binary hashes, CPU and load. Pair only within one workload and runner. Report exact five-pair bootstrap 95% intervals. `reduce_ms` is the Macaulay elimination interval; `wall_ms` is the complete F5 step inside the call, excluding process launch and fixture construction.

## Correctness and promotion

The shared kernel must match textbook RREF bit for bit on randomized, rank-deficient and structured matrices, and preserve echelon row space. Every paired benchmark case must match row fingerprint, rank, pruning and criterion work. Counted elimination word XORs may differ and are retained. Require a frozen primary elimination median reference/candidate ratio of at least 1.20 with its 95% lower bound above one, a complete-call ratio and lower bound above one, both holdout primary gains outside their A/A ranges, and no smaller complete-call median below its own A/A range. If the one-thread gates pass, run the same four-thread control before considering a default. Otherwise preserve the negative result and leave the original default. These are internal solver-stage measurements, not one-target IC or DLP speedup claims.

## Compact fused replay, frozen after the first result

The first column-panel candidate was exact but much slower. Its per-row replay fetched up to four full-stride table entries separately for each pivot block and each narrow column tile. A second distinct candidate copies each block's lookup entries into a compact column-tile table, then applies the entries for one row and block through the existing fused XOR kernel. The pivot scheduling, deferred patterns, two-word panel and one-megabyte tile budget remain fixed. It still uses `KIC_GF2_COLUMN_PANEL_WORDS=2` against the original `=0` reference. Run the same frozen seed and two holdouts, five A/A and five alternating A/B pairs, pinned to one CPU. Use the original correctness and promotion gates. This is one bounded implementation revision; if it misses a gate, retain the negative result without another timing-driven revision.
