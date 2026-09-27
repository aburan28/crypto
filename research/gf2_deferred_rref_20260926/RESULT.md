# Deferred above-row GF(2) RREF: first paired run

The first candidate is **rejected on performance**. `KIC_GF2_DEFER_ABOVE=1` enables it; the unchanged block order remains the default. No F5, IC, or DLP speedup is claimed. The original local host was heavily loaded (one-minute load averages observed at 177, 195, and 63), so the paired comparison ran later on one pinned Ubuntu CI runner.

## Change

The forward pass still chooses the same pivots and clears all rows below each block. It records block starts and pivot columns. The reverse pass clears earlier pivot rows over bounded column ranges, visiting ranges right to left and blocks in reverse order within each range. This order preserves each block's selector bits and ensures its pivot rows have already been reduced by later blocks. `KIC_GF2_DEFER_ABOVE=0` or an unset variable uses the previous RREF path.

## Evidence

| Arm | Correctness | F5 n24/m24/d4 elimination time | Full F5 time | DLP online speedup |
| --- | --- | ---: | ---: | ---: |
| Existing default | Textbook RREF control passed | 161.5 ms | 272.9 ms | unknown |
| Tiled reverse pass | Same exact RREF and rank across the kernel suite and F5 cells | 301.4 ms | 411.4 ms | unknown |

- Base revision: `a71eac61ebc2b28529be2d367a76e0fcc0e56c37`.
- Base `gf2_elim.rs` SHA-256: `5693b918fa8b003448ef340e919e8e80d1635a7c2c3447e4ca2dde98fcdfd258`.
- Command: `TMPDIR=/Volumes/SSD990/llm/tmp rustup run stable rustc --edition=2021 --test src/cryptanalysis/gf2_elim.rs --extern rayon=/Volumes/SSD990/crypto/target/release/deps/librayon-f66c8d45e47d2ddd.rlib --extern rand=/Volumes/SSD990/crypto/target/release/deps/librand-a5b1dd1afa19929b.rlib -L dependency=/Volumes/SSD990/crypto/target/release/deps -o /Volumes/SSD990/llm/tmp/gf2_deferred_tests -C opt-level=1`, then `gf2_deferred_tests --test-threads=1 --quiet`.
- Result: 4 passed, 0 failed. The tests compare both modes bit for bit with textbook RREF on random, structured and rank-deficient matrices; they force two- and three-word reverse tiles across every table-count, SIMD and parallel-threshold setting. A fourth test verifies that one-, two-, three-word and automatic tile widths have the same exact output **and counted word work**. The row-echelon control also passed.
- Compiler: Rust 1.98.0 stable; host: Darwin 25.6 arm64. The shell's Rust 1.93.1 was incompatible with cached dependencies. `/private/tmp` was full, so compiler scratch used `/Volumes/SSD990/llm/tmp`.

The first implementation visited ranges left to right; the two RREF comparison tests failed because an earlier range cleared a block's selector word before later ranges read it. That failure was retained as a design control and resolved by visiting ranges right to left.

The local release builds under Rust 1.93 and 1.98 and a module-test build were stopped under high load without producing benchmark evidence. CI run [36291052511](https://github.com/aburan28/crypto/actions/runs/36291052511) at PR head `4b2866ef` passed the shared-kernel and full matrix-F5 tests, then ran five A/A and five A/B pairs with one Rayon thread and CPU affinity 0. Its complete receipt is [runs/36291052511/ci-result.json](runs/36291052511/ci-result.json), SHA-256 `cadb5884cf23ec4c1d2d3957f62be32aa7648969edd2b1f646b47905bec84b74`. The runner was Linux x86-64 on Azure, Rust 1.98.1, four logical CPUs. Load averaged 1.90 at start and 1.83 at end; a shared virtual runner can still add timing noise.

For the primary n24/m24/d4 cell, the paired median reference/candidate elimination ratio was **0.529** (bootstrap 95% interval 0.515–0.566), outside the A/A range 0.959–1.003. The full F5-call ratio was **0.659** (interval 0.644–0.700), outside the A/A range 0.975–1.009. The candidate reduced counted elimination work from 385,616,841 to 336,771,897 word operations but took longer. All seven cells kept identical fingerprints, ranks, pruning counts and criterion counts. Smaller cells also showed no compelling gain; n20/m20/d4 elimination was 42.7 ms reference versus 52.0 ms candidate.

The first reverse pass still scanned earlier rows separately for every pivot block. The next bounded iteration groups prepared blocks and applies them during one row visit per group. Its timing is unmeasured; this negative result remains the reference for that design decision. No one-target online IC comparison or rho boundary changed, so no scoreboard value is updated.
