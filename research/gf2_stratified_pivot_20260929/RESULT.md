# Spread-out pivot probes: denser F5 output and slower complete call

The local ARM64 screen in [PROTOCOL.md](PROTOCOL.md) completed on
2026-09-29 from candidate head `02c7f06ea`. All 44 benchmark
processes succeeded: two seeds, one warmup per arm, five prior/prior
A/A pairs and five alternating prior/new pairs per seed. Every F5
case preserved rank, canonical row space, criterion word operations,
and row-build counts. The candidate changed the allowed echelon
rows, output terms and counted reduction work. The shared-eliminator
unit suite passed, including stratified-pivot sparse, dense and
partial-word rank and row-space checks, and the release benchmark
built successfully.

The [complete compressed receipt](screen_2026-09-29T193123Z.json.gz)
has SHA-256 `7384636d232b313bcc8fe218c3354711198f4e2a54ed4f522ae0c3824ddc6d6d`.
The uncompressed JSON was 394,518 bytes with SHA-256
`0b282d6ef064f7f412794a491a4a42c1213540046f3a39f6ef0bcdc5fe65e05e`.
It retains every process output and status, phase costs, source and
binary hashes, counted work and host load. No failure, timeout or
OOM row was omitted.

## Frozen n24 complete call

The paired ratio is the median of five prior/new complete-call
ratios; values above one favor the 32 probes. The A/A range is five
prior/prior ratios. Marginal medians describe the same paired calls
but do not calculate the paired ratios.

| Seed | Prior/new complete call | A/A range | Prior / probed terms | Prior / probed reduction word ops |
| --- | ---: | ---: | ---: | ---: |
| `0` | **0.887×** | 0.981–1.444 | 13,734,979 / 16,300,644 | 100,213,183 / 100,330,564 |
| `badc0de1` | 0.863× | 0.921–1.026 | 13,722,549 / 15,980,091 | 100,342,346 / 100,362,691 |

On frozen `0`, marginal complete-call medians were 77.136 ms prior
and 86.996 ms probed, reduction 53.173 and 59.870 ms, and
unpacking 19.331 and 22.441 ms. On `badc0de1`, complete-call
medians were 80.088 and 92.907 ms, reduction 55.381 and
65.577 ms, and unpacking 19.813 and 22.974 ms. All six smaller
cases on both seeds regressed, with paired medians 0.906–0.958×.
The spread-out probes increased output terms by roughly 16–19%; the
particular row choices made the echelon basis denser despite seeking
locally light sampled rows.

The host was macOS 26.6 on Apple ARM64, one Rayon thread, with load
averages 32.89 at start and 24.82 at end on 14 logical CPUs. Timing
is highly contended, so this is a local rejection screen, not an
eligible x86 speed claim. The deterministic output growth, both
primary paired medians, and all smaller-case regressions miss the
preset progression gate. The opt-in path is removed. The exact
tested [source patch](rejected_candidate.patch) and
[screen script](screen.py) remain for replay. The further 2×
complete-call goal remains open. This is a matrix-F5 solver-stage
diagnostic, not an IC online-time or DLP speedup.
