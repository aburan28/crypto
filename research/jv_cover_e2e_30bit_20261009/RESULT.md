# ISO-1 whole-method replay on prime subgroups above 30 bits

Date: 2026-10-09 PDT. Requested size gate: **prime subgroup order** greater than 30 bits. Four fresh seed-5 instances were run at `p = 251, 503, 1009, 1511`, each with eight challenge-construction moves. The source is the historical, buildable §18.4 revision `a2c32c4f9b5153f00a5c2e5b2e579bd7ff89a04a`, the same method used for the earlier 32-run table. The raw [JSON](seed5.json) and [log](seed5.log) are retained.

| p | prime subgroup order ℓ | log₂ ℓ | path 2/odd | walk + transport + model, F_p muls | route, F_p muls | total, F_p muls | end-to-end / pooled rho | solved, correct, challenge relation |
|---:|---:|---:|:--|---:|---:|---:|---:|:--|
| 251 | 62,514,729,468,773 | 45.83 | 3/1 | 102,589,542 | 30,635,906,363 | 30,738,495,905 | 8.6305 | yes / yes / yes |
| 503 | 4,049,001,322,975,577 | 51.85 | 4/0 | 107,968 | 41,295,428,476 | 41,295,536,444 | 1.4407 | yes / yes / yes |
| 1009 | 263,807,419,858,934,933 | 57.87 | 1/0 | 159,450 | 13,525,395,634 | 13,525,555,084 | 0.05846 | yes / yes / yes |
| 1511 | 2,975,272,818,899,256,217 | 61.37 | 0/2 | 32,930,550 | 8,463,167,908 | 8,496,098,458 | 0.01093 | yes / yes / yes |

The runner constructs a weak-curve instance with `#E = 4ℓ`, transports its generator and target through eight steps to a non-weak challenge curve, walks to a weak curve, transports the points, changes models, and solves on the cover. `route_correct` compares the recovered scalar with the planted answer; the challenge relation check verifies `[d]G = Q` on the challenge curve when the route is correct. The challenge construction is excluded from the charged total, as specified in §18. The comparison denominator is the **pooled rho reference** `rho_ref = 1.3609` from the earlier campaign, rather than a newly run matched rho trial for each instance. These are constructed reachable classes; the result does not measure the probability of reach from an arbitrary curve.

## Reproduction and checks

From the historical revision:

```sh
cargo build --release --locked --example jv_isogeny_walk
target/release/examples/jv_isogeny_walk --e2e --sizes 251,503,1009,1511 --seeds 5 --e2e-moves 8 --json seed5.json > seed5.log 2>&1
```

The build exited 0 on Apple M4 Pro (`arm64`, Darwin 25.6.0), Rust/Cargo 1.93.1. The run exited 0. Its elapsed seconds by row were 125.583, 111.449, 52.085, and 57.252, respectively, on a contended host. The JSON audit checked four rows, `bits > 30`, non-weak challenge, walk found, solved, correct, challenge relation, stage `ok`, and exact equality between each total and its four charged phases. PARI/GP independently returned `isprime(ℓ) = 1` for all four orders. These controls check the recorded subgroup gate and internal consistency; the curve-order assignment itself uses the code's two-point randomized point counter.

Source file SHA-256: `jv_isogeny_walk.rs` (module) `c6822fb1754315aa8626d3a2adf6381ce5fb9fb7ab428164f12d4776cb68cb73`; `jv_cover.rs` `933334abbf38c8b1db5ada51689b2e0608f3fb89bc8a3c3169574c9fedf89457`; `jv_sieve.rs` `2d249195e343cd3cd7af84e713d3d493d03f10732a9121206614297db54eaf3d`; example driver `51638a7eaa2ee9da270ac2247f205ef25a819e032a2aa415a9aadfece050a6e7`. The executable SHA-256 was `bb22e6f3f4f8acb2482ce0d34ff1588096a81546ee8ec2afb18c9ccaac59cd44`; the JSON SHA-256 was `a8369b2e1ca69fe6b06f09036772b0729fe0bf6dc1bea29c69b1bf41c5710715`, and the log SHA-256 was `5cc201305e79cf1e246d013fd0aa92462529102d83528810653065f1cc684f65`.

The current `main` branch has subsequent census and walk changes; this evidence records the stated historical revision. No new source code is introduced by this evidence commit.

## Release test gate

`cargo test --release --locked --lib` on the measured revision compiled and ran 3,768 tests: **3,648 passed, 11 failed, 109 ignored**. Two `jv_cover` tests failed in that contended run but each passed when rerun alone with `--nocapture`: `replayed_f4_agrees_with_the_full_solver_and_the_oracle` (58.16 s) and `the_solver_loses_no_rational_solution_at_small_p` (18.30 s). The direct `jv_isogeny_walk::tests::end_to_end_from_a_non_weak_curve_recovers_and_verifies` test passed in the full run. Nine `pollard_collab` socket tests failed with `Operation not permitted` under the filesystem/network sandbox; a release rerun of the whole `pollard_collab` module with local socket access passed **37/37** in 0.66 s. The aggregate test command was not rerun after these focused checks, so its recorded exit remains failure.
