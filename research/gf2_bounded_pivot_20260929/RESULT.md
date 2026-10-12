# Eight-candidate pivot search: output sparsity lost, local gate missed

The local ARM64 screen in [PROTOCOL.md](PROTOCOL.md) completed on
2026-09-29 from candidate head `541b7087d`. All 44 benchmark
processes succeeded: two seeds, one warmup per arm, five prior/prior
A/A pairs and five alternating prior/new pairs per seed. Every F5
case preserved rank, canonical row space, criterion word operations,
and row-build counts. The candidate changed the allowed echelon rows,
output terms and counted reduction work. The shared-eliminator unit
suite passed, including the bounded-pivot sparse, dense and
partial-word rank and row-space checks, and the release benchmark
built successfully.

The [complete compressed receipt](screen_2026-09-29T192242Z.json.gz)
has SHA-256 `d1d0553cfaa0bbef3665028986f1507f30cebd5bb77fed1626c7fa02cd94dfd1`.
The uncompressed JSON was 392,350 bytes with SHA-256
`5d355ef50240789235de88cbbe18409fd774d533eab5367db2d9bf7166496b39`.
It retains every process output and status, phase costs, source and
binary hashes, counted work and host load. No failure, timeout or
OOM row was omitted.

## Frozen n24 complete call

The paired ratio is the median of five prior/new complete-call
ratios; values above one favor bounded pivots. The A/A range is five
prior/prior ratios. Marginal medians describe the same paired calls,
but do not calculate the paired ratios.

| Seed | Prior/new complete call | A/A range | Prior / bounded terms | Prior / bounded reduction word ops |
| --- | ---: | ---: | ---: | ---: |
| `0` | 1.051× | 0.771–1.168 | 13,734,979 / 13,547,720 | 100,213,183 / 100,048,539 |
| `badc0de1` | **0.905×** | 0.991–1.042 | 13,722,549 / 13,447,441 | 100,342,346 / 100,169,394 |

On `badc0de1`, marginal complete-call medians were 80.375 ms prior
and 87.295 ms bounded, reduction 54.408 and 62.596 ms, and
unpacking 19.472 and 19.272 ms. The frozen `0` A/A range was wide
and its marginal time was distorted by contention; its apparent
1.051× paired median cannot establish a gain. All six smaller cases
on both seeds regressed, with paired medians from 0.898× to 0.971×.
The eight-candidate limit retained only 1.4–2.0% of output-term
savings and about 0.2% of reduction-word savings, while increasing
the pivot-search cost.

The host was macOS 26.6 on Apple ARM64, one Rayon thread, with load
averages 18.12 at start and 16.19 at end on 14 logical CPUs. This
is a local rejection screen, not an eligible x86 speed claim. The
holdout primary and every smaller case miss the preset progression
gate. The opt-in path is removed. The exact tested
[source patch](rejected_candidate.patch) and [screen script](screen.py)
remain for replay. A follow-on may probe a fixed number of rows
spread across the strip without scanning it all. The further 2×
complete-call goal remains open. This is a matrix-F5 solver-stage
diagnostic, not an IC online-time or DLP speedup.
