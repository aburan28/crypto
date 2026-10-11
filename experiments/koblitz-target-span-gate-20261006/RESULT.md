# Target-independent rank-prefix shortcut fails this frozen gate

The preregistered [target-span protocol](PROTOCOL.md) was committed as
`927751505` and opened in PR #1462 before the analyzer read the frozen
rank/target transcripts. The native Rust analyzer checked all **24/24**
successful held-out runs from the prior base-size panel. Those runs comprise
only **four distinct point/base cells**: six deterministic repeats of each
cell, not 24 independent targets. Each transcript still independently
replays its base logs, rank rows, target relation and scalar in the parent
evidence. The new analyzer additionally checked each modular row against
the archived solved logs and reconstructed the same target scalar from its
four point labels.

| Curve | Orbit columns `K` | Actual usable base points | Admission prefix `floor(0.8K)` | First target-span rank in all six repeats | Rank probes in the frozen stream | Rows/probes saved by an exact prefix stop |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| n41 source Koblitz | 64 | 5,248 | 51 | **64** | 2,993,628 | **0 / 0** |
| n41 source Koblitz | 85 | 6,970 | 68 | **85** | 2,100,243 | **0 / 0** |
| n53 source Koblitz | 160 | 16,960 | 128 | **160** | 32,315,380 | **0 / 0** |
| n53 source Koblitz | 220 | 23,320 | 176 | **220** | 18,413,613 | **0 / 0** |

The n41 runs share public point
`[1071506060992,898053054019]` and recovered scalar `449905522971`;
the n53 runs share `[7960849849661793,7443722527872608]` and scalar
`17385600002971`. Within each cell, all six repeats had the same point,
base digest, scalar and probe count. The target row was outside the span of
every proper prefix and entered the span only with the last independent
rank equation. This is an exact statement about those frozen row streams,
not a probability estimate over unseen points or a lower bound for all
relation policies.

**Decision:** reject a new producer that simply truncates these
target-independent rank streams. The prespecified 80% gate fails in every
cell, and even a stop at the first exact span saves nothing. A target-guided
rank policy would be a different algorithm: its rank continuation depends
on Q and must be charged to the one-target **online** interval. The next
source-curve producer test should hold actual base size and Q fixed while
reducing probes per independent rank row or index construction cost, then
replay full rank, target scalar and same-Q rho under complete cold and online
accounting. Separately, no source-curve result transfers to the degree-263
ECC2K-130 descendant until an equal-useful-base native/transported/pullback
m≥3 PDP comparison has verified natural-target yield and complete costs.
This diagnostic makes no CPU speedup claim and adds no dashboard ratio.

## Reproduction and custody

From the repository root, verify the parent's frozen inputs and run the
checked-in analyzer:

```sh
shasum -a 256 -c experiments/koblitz-base-size-cold-panel-20261004/HOLDOUT_SHA256SUMS --quiet
cargo test --offline --locked --example koblitz_target_span_gate
cargo run --offline --locked --example koblitz_target_span_gate -- "$PWD" experiments/koblitz-target-span-gate-20261006/RESULT.json
shasum -a 256 -c experiments/koblitz-target-span-gate-20261006/RESULT_SHA256SUMS --quiet
```

The input manifest check passed; the analyzer's positive proper-prefix and
missing-pivot control passed; the complete 24-run analysis exited 0 and
reported `admit_new_producer=false`. The source worktree lacked the
repository-ignored `Cargo.lock`, so `cargo generate-lockfile --offline`
created the exact lock used here; its copy is
[`analyzer-Cargo.lock`](analyzer-Cargo.lock). This is analyzer dependency
custody, not a new producer build or timing baseline. The result and source
hashes are in [`RESULT_SHA256SUMS`](RESULT_SHA256SUMS). The environment was
Darwin 25.6.0 arm64, Homebrew `rustc` and `cargo` 1.93.1. No host-isolation
receipt was claimed or needed for this algebraic result.
