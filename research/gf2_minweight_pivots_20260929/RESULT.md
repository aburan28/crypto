# Current-weight pivot screen: rejected

The local screen specified in [PROTOCOL.md](PROTOCOL.md) completed on
2026-09-29. The result does not pass its stop gate: choosing the lightest
available row for every pivot substantially slows the complete matrix-F5
call. The experimental eliminator was restored; its exact source change is
preserved in [rejected_candidate.patch](rejected_candidate.patch). This is a
solver-stage diagnostic, not an IC online-time or DLP speedup.

The raw receipt is
[screen_2026-09-29T174623Z.json](screen_2026-09-29T174623Z.json). It records
all 24 benchmark processes, seven cases per process, individual phase times,
row signatures, counted work, source and binary hashes, and host details.
Each seed had one warmup per arm and three rotating paired triads. There were
no process failures, timeouts, or OOMs. The host was Apple ARM64 running
macOS 26.6, with one Rayon thread. These timings are a local rejection
screen; they are not an eligible Linux x86-64 speed claim.

## Frozen and holdout n24 result

Values below are medians of the three paired runs, in milliseconds. Ratios
are medians of the three corresponding per-pair complete-call ratios; values
below one mean the candidate is slower. The `fast` arm is the accepted
packed-kernel mode with ordinary scalar unpack, and `unpack` enables direct
scalar unpack. `weighted` adds current-weight pivot choice to `unpack`.

| Seed | Arm | Complete call (ms) | Reduction (ms) | Unpack (ms) | Terms | Reduce word ops | Fast / weighted | Unpack / weighted |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Frozen `0` | fast | 94.847 | 64.425 | 24.703 | 13,734,979 | 100,213,183 | — | — |
| Frozen `0` | unpack | 89.092 | 61.513 | 22.189 | 13,734,979 | 100,213,183 | — | — |
| Frozen `0` | weighted | 119.016 | 99.761 | 13.866 | 7,216,068 | 91,715,296 | 0.806× | 0.788× |
| Holdout `badc0de1` | fast | 92.831 | 62.923 | 24.764 | 13,722,549 | 100,342,346 | — | — |
| Holdout `badc0de1` | unpack | 91.297 | 64.303 | 21.408 | 13,722,549 | 100,342,346 | — | — |
| Holdout `badc0de1` | weighted | 166.362 | 146.639 | 14.593 | 7,180,533 | 90,691,531 | 0.552× | 0.541× |

The candidate produced roughly half as many output terms and reduced
counted word XORs by about 9–10%, but its reduction time rose enough to
overwhelm the unpack saving. The full candidate scans every remaining strip
row at each pivot and recomputes a row popcount after each nontrivial table
update. The timing establishes that this implementation is slower; the
receipt does not isolate a precise cost for either of those two operations.
The local timing had visible outliers, so no small speed difference would be
credible here. The regression is large and appears on both seeds.

For every one of the seven benchmark cases on both seeds, the `fast`,
`unpack`, and `weighted` arms agreed on rank, canonical row-space
fingerprint, criterion work, and F5 build counts. The weighted arm changed
raw echelon rows and term counts on the degree-4 n20 and n24 cases; this is
allowed by the protocol because echelon output is not unique. All seven
cases had `unpack / weighted` paired medians below 1.00 on both seeds. The
frozen n24 medians were 0.788× and 0.541×, below the 0.95 local stop gate.

The experimental code passed `cargo test --offline gf2_elim --lib` (six
tests, including the new minimum-weight rank and row-space test) and the
release benchmark built with `cargo build --offline --release --example
f4_f2_bench`. The patch retains the test for replay. No x86 confirmation was
run after the local stop gate failed. The prior x86 fast-path reference of
125.89 ms belongs to a different host and is not combined with these local
ratios. The requested further 2× complete-call improvement remains unmet.
