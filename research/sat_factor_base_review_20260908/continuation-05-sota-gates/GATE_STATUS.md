# Koblitz index-calculus SOTA gate status

Current through Stage 53, 2026-09-12. The machine-readable current audit is
`stage-53-current-gate-audit-20260912/audit.json`; its latest hosted evidence is
sealed in `stage-52-n53-lazy-result-20260912/verification.json`.

The campaign has target-independent algebraic factor bases; matched native-XOR,
WDSat, CryptoMiniSat, direct-MITM, GGMP, and signed-Frobenius-rho controls; a
balanced 160-instance PDP panel through `n=59`; public unknown-scalar end-to-end
runs at degrees 23, 31, and 41; and a same-target known-answer construction,
rank, solve, and rho comparison at `n=53`. It has not passed all seven gates and
does not establish a new state of the art.

| Gate | Status | Current evidence | Remaining requirement |
|:--|:--|:--|:--|
| 1. Charge every stage and resource | **Partial overall** | The corrected Phase-B matrix charges 21,038.596137 core-seconds. Stage 39 charges fixed algebraic construction, relations, LA, five unknown-scalar descents, rho, clean build, wall/core, and tree RSS. Stage 40 charges every one/four/four/one cell. Stage 52 verifies the latest clean build, 1,024/4,096 direct arms, and exact same-target rho: 455.588094 core-seconds, 149.479894 sequential wall-seconds, and 1,718,161,408 bytes maximum sampled tree RSS. | Licensed Magma process resources are absent. Preinstalled OS/toolchain acquisition remains an explicit exclusion. |
| 2. WDSat, CryptoMiniSat, Magma F4, MITM, GGMP | **Partial** | Native XOR SAT, WDSat, CryptoMiniSat, and direct MITM ran on the exact 160-input packet. Standard `n=31`/`n=41` and GGMP `n=31` are represented. The Stage-32 successor removes the two original WDSat buffer errors without rewriting Stage 26. | Execute all 160 Stage-22 Magma inputs on a licensed host under the frozen one-thread/no-retry contract, seal the return before truth scoring, and report F4 resources. |
| 3. Single-core, core-seconds, memory, conflicts, wall | **Partial because Magma is absent** | Every executed Phase-B arm retains wall, core-seconds, peak RSS, conflicts or operations, tree memory, and workflow wall. Stages 46–52 add metered `n=53` direct/rho RSS, CPU, wall, support/query counts, inversion/wave counts, and watchdogs. | Supply the same fields for licensed Magma F4. The original local `n=53` pointwise autolab recorded a 16 GiB cap without enforcing or measuring it; Stage 42 supersedes that resource receipt. |
| 4. Scale through `n=31`, `n=41`, and a larger PDP regime | **Satisfied for finite execution coverage** | Phase B covers `n=31`, GGMP `n=31`, `n=41`, and PDP-only `n=59`. Unknown-scalar end-to-end controls cover `n=31` and `n=41`. Stage 42 constructs, ranks, solves, and verifies an exact same-target known-answer `n=53` instance. | The evidence is finite and toy-sized; it is not an asymptotic scaling law. The `n=59` arm remains PDP-only. |
| 5. Unknown scalar with no constructed factor-base logs | **Satisfied for finite degrees 23, 31, and 41; n=53 known-answer only** | Stage 39 defines the `n=41` factor base from a fixed algebraic rule with zero discovery target samples or scalar labels, derives every factor-base log from relations, and solves five hash-derived targets whose scalars were not constructed. | Stage 42 uses one explicit public `n=53` scalar retained only for validation. An unknown-scalar `n=53` run and unaffiliated replay remain open. |
| 6. Full cost against automorphism-optimized Pollard rho | **Finite online crossover; amortized/full-cost crossover false** | Stage 39's `n=41` online descent is 3.502068 times faster than rho, while amortized IC is 5.025825 times slower and fresh build plus science is 110.791382 times rho. Stage 40 reduces four-core amortized wall to 3.605523 times rho while spending 1.286948 times the one-core CPU. Stage 52 verifies the preferred 4,096-cursor direct on the exact same `n=53` target: it is 5.166323 times slower than rho, 1.073641 times faster in wall and uses 0.956319 times the CPU of 1,024 cursors in the same hosted run; fresh build plus preferred direct is 22.668683 times rho. | The `n=41` crossover assumes the factor-base log database exists. No amortized or whole-process cell crosses. Parallel latency reductions spend additional CPU. All improvements are finite constants, not exponent changes. |
| 7. Independent external reproduction and novelty review | **Missing** | Issue [#97](https://github.com/aburan28/crypto/issues/97) contains current source pins, the Stage-22 Magma packet, Phase-B and Stage-35 artifacts, exact verifier commands, primary-source links, and the `CONCUR` / `QUALIFIED` / `BREAKS` format. Outreach is also open at [mtrimoska/EC-Index-Calculus-Benchmarks#1](https://github.com/mtrimoska/EC-Index-Calculus-Benchmarks/issues/1). | An unaffiliated reviewer must return a sealed reproduction and source-pinned novelty/correctness assessment. Project-authored CI and replays do not satisfy independence. |

## Current Phase-B matrix

The corrected 480-row SAT view has 141 true positives, 120 true negatives, 219
inconclusive outcomes, no false classifications, and no solver errors. Standard
`n=31` and `n=41` resolve all 40 targets under all three SAT backends. The
`n=31` GGMP cell yields 2 native, 1 WDSat, and 18 CryptoMiniSat true positives;
the rest are inconclusive. All three SAT backends are inconclusive on the 40
`n=59` targets under their caps. Direct MITM classifies all 160 inputs: 80 true
positives and 80 true negatives.

| Backend | Core-seconds | Summed process wall | Peak RSS | Conflicts or operations |
|:--|--:|--:|--:|--:|
| Native XOR SAT | 958.931674 | 1,278.114043 s | 140,333,056 B | 9,233,732 conflicts |
| WDSat, corrected view | 7,404.749336 | 9,524.938771 s | 24,252,416 B | 7,599,249 conflicts |
| CryptoMiniSat | 6,237.273520 | 8,137.581155 s | 122,703,872 B | 19,480,743 conflicts |
| Direct MITM | 759.710741 inclusive outer | 766.538655 s one-CPU elapsed | 85,479,424 B tree | 4,804,252 additions; 4,768,800 pair entries |

The original Stage-26 cells remain immutable. Stage 32 adds 426.718011 charged
core-seconds to turn the two `n=59` WDSat buffer assertions into clean
timeout-inconclusive terminals. The classifications do not change.

## Current n=41 and n=53 end-to-end results

Stage 39 fixes the `n=41` factor base algebraically as the two-torsion
saturation of the Frobenius union generated by masks `[1,2,4,8,16,32]`.
Construction uses zero target samples or discrete-log labels: 2,380 abscissae,
4,759 rational points, 60 signed orbits, and 29 projected columns.

It derives and certifies all 29 logs from 35 relations over 114,688 probes.
Five domain-separated public hash targets solve in 0.259641 seconds versus
0.909281 seconds for signed-Frobenius rho. Precomputation is 4.310245 seconds,
so the online ratio crosses while the amortized ratio remains 5.025825. The
four-core ABBA control reduces median precomputation to 2.010618 seconds and
complete five-target IC wall to 2.479844 seconds, while increasing median
process CPU by 28.7 percent.

Stage 42 uses a point-defined `n=53` base with 9,964 points and 94 orbit
columns, selected without scalar labels. Direct and rho receive the same public
point `Q=(2565091273463387,5885236316843894)` and both recover validation scalar
`476811900269`. Stages 46–52 retain the same 189 verified four-summand relation hashes while
reusing query scratch, increasing the parallel batch width, and loading pair
labels only after exact hits. In the latest same-run comparison, 4,096 cursors
reduce wall from 28.352866 to 26.408147 seconds and CPU from 63.847327 to
61.058406 core-seconds versus 1,024 cursors. Packed signed-Frobenius rho takes
5.111594 seconds, so preferred direct remains 5.166323 times slower; fresh build
plus direct remains 22.668683 times slower. Peak direct RSS falls to
1,392,087,040 bytes and the retained support allocation falls by 134,217,728
bytes.

The narrow supported conclusion is unchanged: this is strong internal
engineering and finite public toy-research evidence. The known SAT-based
point-decomposition and Frobenius-invariant-factor-base literature remain the
prior-art baseline. The finite `n=41` online descent is faster after
precomputation; rho remains faster under amortized and full available
accounting, and the `n=53` direct arm loses on the exact same target. Licensed
Magma and unaffiliated reproduction/novelty review remain open. This is not a
new Koblitz index-calculus SOTA result.
