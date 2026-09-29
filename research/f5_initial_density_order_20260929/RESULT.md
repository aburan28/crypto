# Initial-density F5 row order: local progression gate missed

The local ARM64 screen frozen in [PROTOCOL.md](PROTOCOL.md) completed
on 2026-09-29 from candidate head `ca3559a9e`. All 44 benchmark
processes succeeded: two seeds, one warmup per arm, five prior/prior
A/A pairs and five alternating prior/new pairs per seed. Every case
preserved rank, canonical row space, criterion word operations, and
row-build counts. The candidate changed the allowed raw echelon rows,
term counts and reduction work. The unit tests passed, including the
small, sparse, dense and partial-word row-space test, and the release
benchmark built successfully.

The [complete compressed receipt](screen_2026-09-29T190556Z.json.gz)
has SHA-256 `1b8e8b62edce90db81fcb7efd052a22dd6b72fe00229c12af6f60ab38f355065`.
The uncompressed JSON was 396,823 bytes with SHA-256
`853edf8010e3e94f7cf70b6b61871d05a82818791b7d2970788d64ec13e749d9`.
It retains all process outputs, phase costs, source and binary hashes,
row signatures, counted work, statuses, and host load. No failure,
timeout or OOM was omitted.

## Frozen n24 complete call

The paired ratio is the median of five prior/new ratios; values above
one favor the sort. The A/A range is five prior/prior ratios.
Marginal time medians describe the same five paired runs but are not
used to calculate the paired ratio.

| Seed | Prior/new complete call | A/A range | Prior / sorted terms | Prior / sorted reduction word ops |
| --- | ---: | ---: | ---: | ---: |
| `0` | **1.006×** | 0.537–1.026 | 13,734,979 / 13,247,993 | 100,213,183 / 100,032,160 |
| `badc0de1` | 0.987× | 0.937–1.553 | 13,722,549 / 13,182,994 | 100,342,346 / 100,151,146 |

On frozen `0`, marginal complete-call medians were 81.618 ms prior
and 79.068 ms sorted, reduction 57.289 and 55.674 ms, and unpacking
19.599 and 18.727 ms. On `badc0de1`, complete-call medians were
79.481 and 78.190 ms, reduction 54.606 and 55.213 ms, and unpacking
19.481 and 18.690 ms. The paired ratios are more informative than
these marginal medians under the observed load. Across the six smaller
cases, paired medians ranged from 0.955× to 1.022×; no large gain
appeared.

The machine was macOS 26.6 on Apple ARM64, with one Rayon thread.
Its load average was 18.49 at start and 14.02 at end on 14 logical
CPUs. The wide A/A ranges preclude a small timing claim. The
deterministic changes are also modest: 3.5–3.9% fewer returned terms
and about 0.2% fewer reduction word operations. Both primary paired
medians miss the preset 1.05 local progression gate. This screen is
therefore stopped without an eligible x86 run. The opt-in runtime path
has been removed; the exact tested [source patch](rejected_candidate.patch)
and [screen script](screen.py) remain for replay. The further 2×
complete-call goal remains open. This is a matrix-F5 solver-stage
diagnostic, not an IC online-time or DLP speedup.
