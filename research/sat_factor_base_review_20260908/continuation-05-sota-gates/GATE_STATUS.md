# Koblitz index-calculus SOTA gate status

Current through Stage 59, 2026-09-12. The machine-readable current audit is
`stage-59-current-gate-audit-20260912/audit.json`; its latest hosted evidence is
sealed in the Stage 56 PCLMUL and Stage 58 unknown-scalar results.

The campaign has target-independent algebraic factor bases; matched native-XOR,
WDSat, CryptoMiniSat, direct-MITM, GGMP, and signed-Frobenius-rho controls; a
balanced 160-instance PDP panel through `n=59`; public unknown-scalar end-to-end
runs at degrees 23, 31, 41, and 53; and same-target known-answer and
scalar-blind construction, rank, solve, and rho comparisons at `n=53`. It has not passed all seven gates and
does not establish a new state of the art.

| Gate | Status | Current evidence | Remaining requirement |
|:--|:--|:--|:--|
| 1. Charge every stage and resource | **Partial overall** | The corrected Phase-B matrix charges 21,038.596137 core-seconds. Stage 39 charges fixed algebraic construction, relations, LA, five unknown-scalar descents, rho, clean build, wall/core, and tree RSS. Stage 40 charges every one/four/four/one cell. Stage 56 charges clean build plus portable/PCLMUL/rho at 428.232650 core-seconds, 140.410989 sequential wall-seconds, and 1,709,281,280 bytes tree RSS. Stage 58 charges clean build plus unknown-scalar direct/rho at 368.127516 core-seconds, 117.138948 sequential wall-seconds, and 1,708,228,608 bytes tree RSS. | Licensed Magma process resources are absent. Preinstalled OS/toolchain acquisition remains an explicit exclusion. |
| 2. WDSat, CryptoMiniSat, Magma F4, MITM, GGMP | **Partial** | Native XOR SAT, WDSat, CryptoMiniSat, and direct MITM ran on the exact 160-input packet. Standard `n=31`/`n=41` and GGMP `n=31` are represented. The Stage-32 successor removes the two original WDSat buffer errors without rewriting Stage 26. | Execute all 160 Stage-22 Magma inputs on a licensed host under the frozen one-thread/no-retry contract, seal the return before truth scoring, and report F4 resources. |
| 3. Single-core, core-seconds, memory, conflicts, wall | **Partial because Magma is absent** | Every executed Phase-B arm retains wall, core-seconds, peak RSS, conflicts or operations, tree memory, and workflow wall. Stages 46–58 add metered `n=53` direct/rho RSS, CPU, wall, support/query counts, inversion/wave counts, target derivation, and watchdogs. | Supply the same fields for licensed Magma F4. The original local `n=53` pointwise autolab recorded a 16 GiB cap without enforcing or measuring it; Stage 42 supersedes that resource receipt. |
| 4. Scale through `n=31`, `n=41`, and a larger PDP regime | **Satisfied for finite execution coverage** | Phase B covers `n=31`, GGMP `n=31`, `n=41`, and PDP-only `n=59`. Unknown-scalar end-to-end controls cover `n=31` and `n=41`. Stage 42 constructs, ranks, solves, and verifies an exact same-target known-answer `n=53` instance. | The evidence is finite and toy-sized; it is not an asymptotic scaling law. The `n=59` arm remains PDP-only. |
| 5. Unknown scalar with no constructed factor-base logs | **Satisfied for finite degrees 23, 31, 41, and 53** | Stage 58 derives `n=53` target `Q` from public hash seed 53001 without constructing or supplying its scalar, derives all 94 factor-base logs from 165 relations, and has direct IC and rho independently recover the same `d` with `[d]G=Q`. | Repeat on independent n=53 public seeds and obtain unaffiliated replay; these are strengthening steps rather than prerequisites for the finite gate-5 execution. |
| 6. Full cost against automorphism-optimized Pollard rho | **Finite online crossover; amortized/full-cost crossover false** | Stage 39's `n=41` online descent is 3.502068 times faster than rho, while amortized IC is 5.025825 times slower and fresh build plus science is 110.791382 times rho. Stage 40 reduces four-core amortized wall to 3.605523 times rho while spending 1.286948 times the one-core CPU. Stage 56 verifies PCLMUL on the same binary and target: it is 1.418695 times faster in wall, uses 0.614352 times portable CPU, remains 3.688473 times slower than rho, and remains 21.182086 times rho with fresh build. Stage 58 unknown-scalar direct remains 4.428199 times rho and 25.722944 times with fresh build. | The `n=41` crossover assumes the factor-base log database exists. No amortized or whole-process cell crosses. Parallel latency reductions spend additional CPU. All improvements are finite constants, not exponent changes. |
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
`476811900269`. Stages 46–56 retain the same 189 known-answer relation hashes while reusing
query scratch, widening batches, loading pair labels only after exact hits, and
using x86 carry-less multiplication. In the same-binary Stage 56 comparison,
PCLMUL reduces wall from 26.773640 to 18.872020 seconds and CPU from 61.623592
to 37.858583 core-seconds. It remains 3.688473 times slower than rho and
21.182086 times slower with fresh build. Stage 58 then derives an unrelated
public hash target without constructing its scalar, derives all factor-base
logs from 165 relations, and has direct IC and rho independently recover
`7892094459170` with `[d]G=Q`. Unknown-scalar direct remains 4.428199 times rho
and 25.722944 times rho with fresh build.

The narrow supported conclusion is unchanged: this is strong internal
engineering and finite public toy-research evidence. The known SAT-based
point-decomposition and Frobenius-invariant-factor-base literature remain the
prior-art baseline. The finite `n=41` online descent is faster after
precomputation; rho remains faster under amortized and full available
accounting, and the `n=53` direct arm loses on the exact same target. Licensed
Magma and unaffiliated reproduction/novelty review remain open. This is not a
new Koblitz index-calculus SOTA result.
