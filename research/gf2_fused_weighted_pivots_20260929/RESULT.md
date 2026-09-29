# Fused weighted-pivot discovery: rejected

The frozen Apple ARM64 screen in [PROTOCOL.md](PROTOCOL.md) completed on
2026-09-29. Moving next-pivot discovery into the strip update preserved
the previous weighted pivot result exactly and shortened some calls, but
the full matrix-F5 call remained slower than the accepted fast arm on both
seeds. The experiment fails its 0.95 local stop gate and was not sent to
the eligible x86 runner. The production eliminator was restored. The exact
experimental source is [rejected_candidate.patch](rejected_candidate.patch).

The [raw receipt](screen_2026-09-29T175730Z.json) contains all 24 process
records, seven cases per process, host data, source and binary digests,
full benchmark output, phases, row signatures, counts and status. Each seed
had one warmup per arm and three rotating paired triads. Every process
succeeded. The host was macOS 26.6 on Apple ARM64, with one Rayon thread.
It had visible timing outliers, so small local differences are not treated
as gains.

## Full n24 call, in milliseconds

Times are medians of three paired calls. Ratios are medians of the three
corresponding per-pair ratios. `fast` uses the accepted packed elimination
and direct scalar unpack; `weighted` adds the previous full-scan pivot
choice; `fused` makes the same choice during strip clearing.

| Seed | Arm | Complete call | Reduction | Unpack | Output terms | Reduce word ops | Fast / fused | Weighted / fused |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Frozen `0` | fast | 91.926 | 60.533 | 25.914 | 13,734,979 | 100,213,183 | — | — |
| Frozen `0` | weighted | 108.484 | 88.526 | 14.183 | 7,216,068 | 91,715,296 | — | — |
| Frozen `0` | fused | 104.935 | 85.919 | 13.285 | 7,216,068 | 91,715,296 | **0.936×** | 1.070× |
| Holdout `badc0de1` | fast | 86.857 | 58.448 | 23.846 | 13,722,549 | 100,342,346 | — | — |
| Holdout `badc0de1` | weighted | 111.857 | 92.442 | 14.009 | 7,180,533 | 90,691,531 | — | — |
| Holdout `badc0de1` | fused | 100.283 | 81.354 | 13.136 | 7,180,533 | 90,691,531 | **0.866×** | 1.071× |

Across all seven cases on both seeds, `weighted` and `fused` matched raw
row fingerprints, output term counts, rank, canonical row-space
fingerprints, criterion and build counts, and reduction word operations.
The `fast` arm also matched rank, row space, criterion work and build
counts. The weighted arms changed echelon rows and terms on degree-4 n20
and n24 cases, as allowed by the protocol. For the primary frozen case,
the fused search recovered about 7% relative to the weighted full-scan
implementation, but still lost to `fast`; the holdout showed the same
pattern. Several smaller cases regressed much more.

The experimental code passed `cargo test --offline gf2_elim --lib` (six
tests, including bit-for-bit weighted/fused output and counted-work checks)
and `cargo build --offline --release --example f4_f2_bench`. The patch
retains that test for replay. No eligible-host 2× result, IC online-time
result or DLP speedup is claimed. The requested further 2× full-call goal
remains open.
