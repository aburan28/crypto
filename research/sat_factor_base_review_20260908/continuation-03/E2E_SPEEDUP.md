**End-to-end SAT index-calculus speedup — degree 19**

This continuation measures the complete toy index-calculus algorithm on

\[
K_1:y^2+xy=x^3+x^2+1
\quad\text{over}\quad
\mathbb F_2[z]/(z^{19}+z^5+z^2+z+1),
\]

using the 152-point full-orbit factor base selected in continuation 02. Every admitted run includes curve setup, exact factor-base materialization, repeated native-XOR SAT relation collection, modular relation-matrix solving, recovery of an unknown planted scalar, and independent verification by recomputing the public target point. No direct-relation shortcut occurred.

**Matched arms.** Both arms use the same factor-base points, target scalar, sampler seed, target order, three-summand symmetric `S_4` equations, exact domain trie, native XOR rows, trace constraint, 10,000-conflict per-target cutoff, 1,024-model cap, 64-trial cap, and two-relation fixed-surplus ceiling.

The Frobenius control retains eight ordinary Frobenius relation columns and collects ten relations before solving. The optimized arm identifies `P` with `-P`, leaving four signed Frobenius columns, and attempts the modular solve after the fifth relation. It stops early only when the recovered scalar reproduces the public target; otherwise it continues collecting relations. Thus rank-aware stopping cannot turn an underdetermined or incorrect solve into success.

| Run seed | Unknown scalar | Execution order | Optimized | Frobenius control | Paired speedup |
|---:|---:|---|---:|---:|---:|
| 1 | 4,242 | optimized, control | 85.516961 s | 291.720155 s | 3.411255x |
| 2 | 9,001 | optimized, control | 141.672549 s | 276.528004 s | 1.951881x |
| 3 | 17,001 | optimized, control | 74.646770 s | 194.972048 s | 2.611929x |
| 4 | 30,001 | control, optimized | 126.847856 s | 154.822686 s | 1.220538x |

All four optimized runs and all four controls recovered and verified their unknown scalars. The optimized arm won every pair, including the reverse-order pair.

| Aggregate | Optimized | Frobenius control | Ratio |
|---|---:|---:|---:|
| End-to-end time | 428.684136 s | 918.042893 s | **2.141537x speedup** |
| Target trials | 94 | 184 | 1.957447x fewer |
| SAT conflicts | 865,057 | 1,669,784 | 1.930259x fewer |
| Accepted relations | 20 | 40 | 2x fewer |

The minimum paired speedup is 1.220538x, the median paired speedup is 2.281905x, and the geometric mean is 2.146447x. External `/usr/bin/time` wall measurements agree with the internal end-to-end monotonic timers within two percent for every run.

More than 99.99% of measured algorithm time is relation collection. Factor-base construction takes about one to two milliseconds per process and the dense modular solve takes tens of microseconds. The speedup is therefore an end-to-end reduction of the dominant stage rather than a faster negligible subroutine.

**Conflict-cutoff selection.** Relation collection needs enough easy decompositions; it need not prove every rejected target UNSAT. On seed 1, every cutoff arm completed and verified the same scalar:

| Per-target conflicts | End-to-end time | Trials | SAT UNKNOWN skips | SAT conflicts |
|---:|---:|---:|---:|---:|
| 5,000 | 113.848404 s | 45 | 40 | 216,473 |
| **10,000** | **85.516961 s** | 21 | 16 | 183,775 |
| 20,000 | 151.286810 s | 13 | 8 | 216,927 |

The 10,000-conflict operating point was fixed before the four-pair aggregate. A smaller cap spends too many setup cycles on skipped targets; a larger cap overinvests in hard targets. `UNKNOWN` is counted as an inconclusive skipped attempt and never as a refutation or relation.

**Implementation changes.** `KoblitzIcOptions` now exposes the SAT decomposition controls, signed-column switch, and verified rank-aware stopping. `KoblitzIcReport` records SAT calls, models, conflicts, collection time, linear-algebra time, solve attempts, and the selected column mode. The end-to-end runner emits JSON only after verifying the recovered scalar and rejects direct-relation shortcuts or invalid SAT models.

`verify_e2e_speedup.py` rejects any pair with mismatched seeds, targets, factor bases, conflict budgets, stopping modes, relation counts, direct shortcuts, invalid models, or failed scalar verification. It checks internal versus external time, computes the aggregate only from eight completing runs, and independently verifies that 10,000 conflicts is the fastest completing seed-1 cutoff.

Reproduction commands:

```sh
cp research/sat_factor_base_review_20260908/continuation-01/source_snapshots/Cargo.lock Cargo.lock
cargo build --release --example koblitz_selected_e2e --locked
target/release/examples/koblitz_selected_e2e optimized 1 4242 10000
target/release/examples/koblitz_selected_e2e frobenius 1 4242 10000
target/release/examples/koblitz_selected_e2e optimized 2 9001 10000
target/release/examples/koblitz_selected_e2e frobenius 2 9001 10000
target/release/examples/koblitz_selected_e2e optimized 3 17001 10000
target/release/examples/koblitz_selected_e2e frobenius 3 17001 10000
target/release/examples/koblitz_selected_e2e frobenius 4 30001 10000
target/release/examples/koblitz_selected_e2e optimized 4 30001 10000
python3 research/sat_factor_base_review_20260908/continuation-03/verify_e2e_speedup.py
```

**Evidence boundary.** This is a measured end-to-end constant-factor speedup for complete, verified toy degree-19 index-calculus runs. It does not establish an exponent improvement, a global factor-base optimum, a speedup over Pollard rho, or a threat to deployed binary curves. The conflict cutoff and paired timing result are specific to the recorded implementation and host.
