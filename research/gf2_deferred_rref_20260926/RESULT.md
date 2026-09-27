# Deferred above-row GF(2) RREF: correctness control

The candidate is **not promoted**. `KIC_GF2_DEFER_ABOVE=1` enables it; the unchanged block order remains the default. No F5, IC, or DLP speedup is claimed. The host was shared and heavily loaded during this attempt (one-minute load averages observed at 177, 195, and 63), so the predeclared A/A and paired A/B wall comparison was not run.

## Change

The forward pass still chooses the same pivots and clears all rows below each block. It records block starts and pivot columns. The reverse pass clears earlier pivot rows over bounded column ranges, visiting ranges right to left and blocks in reverse order within each range. This order preserves each block's selector bits and ensures its pivot rows have already been reduced by later blocks. `KIC_GF2_DEFER_ABOVE=0` or an unset variable uses the previous RREF path.

## Evidence

| Arm | Correctness | F5 n24/m24/d4 elimination time | Full F5 time | DLP online speedup |
| --- | --- | ---: | ---: | ---: |
| Existing default | Textbook RREF control passed | unknown | unknown | unknown |
| Tiled reverse pass | Same exact RREF and rank across the kernel suite | unknown | unknown | unknown |

- Base revision: `a71eac61ebc2b28529be2d367a76e0fcc0e56c37`.
- Base `gf2_elim.rs` SHA-256: `5693b918fa8b003448ef340e919e8e80d1635a7c2c3447e4ca2dde98fcdfd258`.
- Command: `TMPDIR=/Volumes/SSD990/llm/tmp rustup run stable rustc --edition=2021 --test src/cryptanalysis/gf2_elim.rs --extern rayon=/Volumes/SSD990/crypto/target/release/deps/librayon-f66c8d45e47d2ddd.rlib --extern rand=/Volumes/SSD990/crypto/target/release/deps/librand-a5b1dd1afa19929b.rlib -L dependency=/Volumes/SSD990/crypto/target/release/deps -o /Volumes/SSD990/llm/tmp/gf2_deferred_tests -C opt-level=1`, then `gf2_deferred_tests --test-threads=1 --quiet`.
- Result: 4 passed, 0 failed. The tests compare both modes bit for bit with textbook RREF on random, structured and rank-deficient matrices; they force two- and three-word reverse tiles across every table-count, SIMD and parallel-threshold setting. A fourth test verifies that one-, two-, three-word and automatic tile widths have the same exact output **and counted word work**. The row-echelon control also passed.
- Compiler: Rust 1.98.0 stable; host: Darwin 25.6 arm64. The shell's Rust 1.93.1 was incompatible with cached dependencies. `/private/tmp` was full, so compiler scratch used `/Volumes/SSD990/llm/tmp`.

The first implementation visited ranges left to right; the two RREF comparison tests failed because an earlier range cleared a block's selector word before later ranges read it. That failure was retained as a design control and resolved by visiting ranges right to left.

The full F5 benchmark build, module tests and paired timing remain outstanding. A release build under Rust 1.93 was stopped when the host load exceeded 200. A second build under Rust 1.98 was stopped when the load rose from about 63 to 96; neither produced a benchmark binary. The offline `cargo test --lib cryptanalysis::matrix_f5_f2::tests` attempt reached compilation of the large `crypto` crate but was stopped after the load rose to about 105 without producing a test binary. Keep this path opt-in until the protocol's exact F5 output and timing gates pass on an uncontended host.
