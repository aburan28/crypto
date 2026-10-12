# Contiguous GF(2) rows: rejected on the local screen

The frozen screen in [PROTOCOL.md](PROTOCOL.md) completed on 2026-09-29.
Converting F5's built matrix to a contiguous word arena before reduction
preserved its raw output but did not improve the complete call. The
experimental runtime code was removed; its exact source and test are in
[rejected_candidate.patch](rejected_candidate.patch). This is a matrix-F5
solver-stage diagnostic, not an IC online-time or DLP speedup.

The [raw receipt](screen_2026-09-29T180841Z.json) records 24 successful
processes: one warmup per arm and five alternating pairs per seed. It
contains all seven cases from every process, full phase times, statuses,
row signatures, counts, host data, and source/binary hashes. The local
host was Apple ARM64 with macOS 26.6 and one Rayon thread. A timing
outlier occurred on one holdout flat call; no run was dropped.

## Complete n24 degree-4 call

Times are medians of the five paired calls, in milliseconds. Ratios are
medians of the five corresponding prior/flat full-call ratios. The flat
conversion and old row deallocation are included in reduction time.

| Seed | Arm | Full call (ms) | Build (ms) | Reduction including conversion (ms) | Unpack (ms) | Prior / flat |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| Frozen `0` | prior | 76.453 | 4.356 | 51.688 | 19.381 | — |
| Frozen `0` | flat | 78.146 | 4.401 | 53.020 | 20.125 | **0.978×** |
| Holdout `badc0de1` | prior | 77.122 | 4.507 | 53.495 | 19.247 | — |
| Holdout `badc0de1` | flat | 76.507 | 4.343 | 52.298 | 19.773 | **1.002×** |

Both primary ratios miss the preregistered 1.05 local gate, so no x86
measurement or promotion followed. On all seven cases on both seeds, the
arms matched raw row fingerprints, canonical row-space fingerprints,
rank, output term count, built/pruned row and column counts, criterion
work, and reduction word XORs. The opt-in flat layout was actually used
for the degree-4 n20 and n24 cases; the smaller or reduced-form cases
remained on their established paths. No other case showed a repeatable
substantial gain from the flag.

The experimental code passed `cargo test --offline gf2_elim --lib` (six
tests, including an exact vector-versus-slice comparison),
`cargo test --offline matrix_f5_f2 --lib` (ten tests), and
`cargo build --offline --release --example f4_f2_bench`. The source patch
applies to the frozen reference commit for replay. The requested further
2× complete-call gain remains unproved.
