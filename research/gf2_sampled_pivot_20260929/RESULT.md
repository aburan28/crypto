# Sampled current-weight pivots: sparse output, slower complete call

The local ARM64 screen in [PROTOCOL.md](PROTOCOL.md) completed on
2026-09-29 from candidate head `7931a52f2`. All 44 benchmark
processes succeeded: two seeds, one warmup per arm, five prior/prior
A/A pairs and five alternating prior/new pairs per seed. Every F5
case preserved rank, canonical row space, criterion word operations,
and row-build counts. The candidate changed the allowed echelon
rows, output terms and counted reduction work. The shared-eliminator
unit suite passed, including sampled-pivot sparse, dense and
partial-word rank and row-space checks, and the release benchmark
built successfully.

The [complete compressed receipt](screen_2026-09-29T191454Z.json.gz)
has SHA-256 `7f074a37b981d6306544b40893cf8df12762030be828772570a16092e402006b`.
The uncompressed JSON was 392,426 bytes with SHA-256
`6f8de2af2a2745fc1875701e4c660ae5fbf68d8bdd304a577907e1e59ca629a6`.
It retains every process output and status, phase costs, source and
binary hashes, counted work and host load. No failure, timeout or
OOM row was omitted.

## Frozen n24 complete call

The paired ratio is the median of five prior/new complete-call
ratios; values above one favor sampled pivots. The A/A range is five
prior/prior ratios. Marginal medians describe the same five paired
runs but do not calculate the paired ratios.

| Seed | Prior/new complete call | A/A range | Prior / sampled terms | Prior / sampled reduction word ops |
| --- | ---: | ---: | ---: | ---: |
| `0` | **0.840×** | 0.875–1.410 | 13,734,979 / 7,546,352 | 100,213,183 / 92,683,711 |
| `badc0de1` | 0.812× | 0.932–1.008 | 13,722,549 / 7,510,747 | 100,342,346 / 92,329,994 |

On frozen `0`, marginal complete-call medians were 76.936 ms prior
and 93.591 ms sampled, reduction 53.581 and 76.740 ms, and
unpacking 19.301 and 12.238 ms. On `badc0de1`, complete-call
medians were 76.350 and 94.258 ms, reduction 52.802 and 77.376 ms,
and unpacking 18.935 and 12.154 ms. All six smaller cases on both
seeds regressed: their paired complete-call medians range from
0.665× to 0.899×. The candidate recovers much of the output
sparsity seen with exact live weights and saves about 7–8% in
counted reduction XORs, but pivot selection raises reduction time
by more than the unpack saving.

The host was macOS 26.6 on Apple ARM64, one Rayon thread, with
load averages 12.41 at start and 14.04 at end on 14 logical CPUs.
The frozen A/A range is wide, so this is a local rejection screen,
not an eligible x86 performance claim. The large loss on both seeds
and every smaller case misses the preset 1.05 local progression gate.
The opt-in path is removed. The exact tested [source patch](rejected_candidate.patch)
and [screen script](screen.py) remain for replay. A follow-on may
bound how many equal-leading-column rows are scored during pivot
search. The further 2× complete-call goal remains open. This is a
matrix-F5 solver-stage diagnostic, not an IC online-time or DLP
speedup.
