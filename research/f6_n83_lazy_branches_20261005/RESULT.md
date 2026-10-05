# n83 F6 wide search: lazy children pass the bounded block gate

This is an exploratory **three-summand planted-block diagnostic**, not a
complete F6-IC speedup. The [preregistered root protocol](PROTOCOL.md)
and [bounded witness amendment](AMENDMENT_1.md) use the pinned K0 n83
curve, standard dimension-12 source subspace, and source indices
`[0,2,4]`. The source has 4,057 curve points; cofactor-four projection
has 4,054 usable subgroup points. No ordinary-query relation was
attempted or counted here.

The search previously copied and substituted the entire reduced
polynomial frame for all **2,029 distinct abscissa choices** at the root.
The protocol called these about 4,057 alternatives, conflating source
points with their sign-paired abscissae; the source has 4,057 points but
only 2,029 distinct abscissae. This correction is recorded here without
rewriting the preregistration.
The candidate holds one shared parent and a small choice record per
child; it substitutes the frame only when that child is visited. The
equations, choice order, node budget, final group verification, and
source-to-subgroup arithmetic are unchanged. The existing small-curve
and n83 wide-backend tests all passed: [eight of eight](wide_tests_candidate.log).

Five alternating pairs ran the **same frozen probe binary source and
input**, with node budget one. Every process exited zero, returned
`exhausted` with no witness, and reported two nodes, two reductions,
166 maximum rows and 15,259 maximum columns. All individual JSON
results, empty stderr files and [statuses](status.tsv) are retained.

| Measure | Eager baseline median | Lazy candidate median | Baseline / candidate |
| --- | ---: | ---: | ---: |
| Query interval | 550.225 ms | 99.900 ms | 5.508× |
| Peak process RSS | 768.393 MB | 20.578 MB | 37.340× |
| Setup plus query, supplementary | 640.404 ms | 176.818 ms | 3.622× |

The five **paired query ratios** were 5.397, 4.085, 4.818, 5.508 and
6.467; median 5.397×, observed range 4.085–6.467×. They are noisy
local observations, not an isolated-host confidence interval. The
query interval includes valid-coordinate construction and the bounded
solver search; setup is reported separately. Process launch is outside
both intervals. Peak RSS comes from `getrusage` after each process run.
The protocol's ≥1.5× query and ≥3× memory gates both passed. These
ratios apply only to the root budget and cannot be transferred to a
complete n83 decomposition or index-calculus call.

The separately registered [128-node planted search](witness_probe.jsonl),
with one Rayon thread and a 120-second cap, **found exactly `[0,2,4]`**
in four nodes and five reductions. Its source sum and cofactor-four
subgroup projection independently replayed; the four-torsion bridge
was checked in the n83 unit test. The query interval was 84.835 ms,
peak RSS 19,972,096 bytes, exit zero, and stderr empty. This makes the
n83 algebraic block executable, but the easy planted indices are a
correctness control. They give no natural relation-yield estimate.

Reproduction used `cargo test --offline --release --lib
cryptanalysis::wide_groebner::tests::`, `cargo build --offline --release
--example f6_n83_wide_block_probe`, `sh
research/f6_n83_lazy_branches_20261005/run.sh`, then `cargo build
--offline --release --example f6_n83_wide_lazy_witness_probe` and
`RAYON_NUM_THREADS=1 gtimeout -k 5s 120s
/Volumes/SSD990/crypto/worktrees/f6-n83-geometric-20261004/target/release/examples/f6_n83_wide_lazy_witness_probe`.
The Cargo commands set `CARGO_TARGET_DIR` to
`/Volumes/SSD990/crypto/worktrees/f6-n83-geometric-20261004/target`
and `TMPDIR` to
`/Volumes/SSD990/crypto/worktrees/f6-n83-batch-size-20261004/target/tmp`.
The paired script uses frozen copies in `target/bench_bins`.
The baseline commit is `9b969fe26130b9eaa2548395bccfc19a0c5ab57a`.
SHA-256: baseline probe binary
`17678d80e0bb221c43b3008e28c73557e79b0b25fec682f3fac6067abb52f16e`,
candidate root probe binary
`16bfcd5408a8c736269b9a2d1359375f8bd3eb2880c194bf11a0333dfd7db22f`,
witness probe binary
`8d79eebbb93634a7a4faa3f58ef3cb1827797d079abe22ef35aedb1dd3cc7cd2`,
solver source `023b8b4147d053dbb5cc6c0b28f89ce02a5465f98791fb9ddb71448398970b0e`,
root-probe source
`fde4a969957de325ab0803adba429382165b119dd97b70ecc0c0d0763166659c`,
witness-probe source
`109c24b6cfcf42b20b016190df04e743299409ea7fa21ed0bc7832551ad78e66`,
and curve/factor-base source
`1c98b492fafe80ae3a9dbea5400e7fceb114dec7b9a548e4d98f33a691f05b5d`.
The protocol and amendment SHA-256 digests are
`7a24d6c7e0488b5f5e7590d7b26ac6e2e54dcc27b44d17b28456f6d8393b52e6`
and `6576807da9e014f51f1cbfeb8605bb96481d46b305176032854b9ecba1f8fa11`.
Hardware: physical Apple M4 Pro, 14 cores (10 performance, four
efficiency), arm64, macOS 26.6, Rust 1.93.1 release profile. The host
was **not isolated**; no CPU-wall speedup claim is promoted under the
repository isolation gate.

The next algorithmic requirement is an eight-or-more-summand ordinary
query solver that returns a verified subgroup relation while charging
all failed attempts. Only then can a complete F6/IC comparison use
candidate IDs, one-target online time and a same-point rho reference.
Those measurements and speedups remain **unknown**.
