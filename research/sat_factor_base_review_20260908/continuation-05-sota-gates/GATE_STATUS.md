# Koblitz index-calculus SOTA gate status

Current through Stage 73, 2026-09-12. The machine-readable current audit is
`stage-73-current-gate-audit-20260912/audit.json`; its latest hosted evidence is
sealed in the Stage 68 fused and Stage 72 unknown-scalar/single-core results.

The campaign has target-independent algebraic factor bases; matched native-XOR,
WDSat, CryptoMiniSat, direct-MITM, GGMP, and signed-Frobenius-rho controls; a
balanced 160-instance PDP panel through `n=59`; public unknown-scalar end-to-end
runs at degrees 23, 31, 41, and 53; and same-target known-answer and
scalar-blind construction, rank, solve, and rho comparisons at `n=53`. It has not passed all seven gates and
does not establish a new state of the art.

| Gate | Status | Current evidence | Remaining requirement |
|:--|:--|:--|:--|
| 1. Charge every stage and resource | **Partial overall** | The corrected Phase-B matrix charges 21,038.596137 core-seconds. Stage 68 charges clean build plus separate/fused direct and rho at 346.758590 core-seconds, 111.312247 sequential wall-seconds, and 1,708,961,792 bytes tree RSS. Stage 72 retains optimized unknown-scalar at 344.394514 core-seconds / 118.116243 sequential wall-seconds / 1,718,992,896 bytes tree RSS and corrected single-core at 253.299090 core-seconds / 89.877466 sequential wall-seconds / 1,715,879,936 bytes tree RSS. | Licensed Magma process resources are absent. Preinstalled OS/toolchain acquisition remains an explicit exclusion. |
| 2. WDSat, CryptoMiniSat, Magma F4, MITM, GGMP | **Partial** | Native XOR SAT, WDSat, CryptoMiniSat, and direct MITM ran on the exact 160-input packet. Standard `n=31`/`n=41` and GGMP `n=31` are represented. The Stage-32 successor removes the two original WDSat buffer errors without rewriting Stage 26. | Execute all 160 Stage-22 Magma inputs on a licensed host under the frozen one-thread/no-retry contract, seal the return before truth scoring, and report F4 resources. |
| 3. Single-core, core-seconds, memory, conflicts, wall | **Partial because Magma is absent** | Every executed Phase-B arm retains wall, core-seconds, peak RSS, conflicts or operations, tree memory, and workflow wall. Stage 68 retains the preferred four-core fused process. Stage 72 retains an inherited Linux CPU-0 affinity receipt: direct 10.945707 wall / 8.576129 core-seconds / 864,251,904 B RSS versus rho 5.490656 wall / 4.334277 core-seconds. | Supply the same fields for licensed Magma F4. |
| 4. Scale through `n=31`, `n=41`, and a larger PDP regime | **Satisfied for finite execution coverage** | Phase B covers `n=31`, GGMP `n=31`, `n=41`, and PDP-only `n=59`. Unknown-scalar end-to-end controls cover `n=31` and `n=41`. Stage 42 constructs, ranks, solves, and verifies an exact same-target known-answer `n=53` instance. | The evidence is finite and toy-sized; it is not an asymptotic scaling law. The `n=59` arm remains PDP-only. |
| 5. Unknown scalar with no constructed factor-base logs | **Satisfied for finite degrees 23, 31, 41, and 53** | Stage 72 derives `n=53` target `Q` from public hash seed 53001 without constructing or supplying its scalar, derives all 94 factor-base logs from 95 verified relations, and has direct IC and rho independently recover `d=7892094459170` with `[d]G=Q`. | Repeat on independent n=53 public seeds and obtain unaffiliated replay; these are strengthening steps rather than prerequisites for the finite gate-5 execution. |
| 6. Full cost against automorphism-optimized Pollard rho | **Finite online crossover; amortized/full-cost crossover false** | Stage 39's `n=41` online descent is 3.502068 times faster than rho, while amortized IC is 5.025825 times slower. Stage 68 fused `n=53` direct is 1.182872 times rho wall and 19.474746 times rho with fresh build. Stage 72 optimized unknown-scalar direct is 1.467611 times rho and 25.359925 times with fresh build. Corrected one-CPU direct is 1.993515 times rho wall and 1.978676 times rho core. | The `n=41` crossover assumes the factor-base log database exists. No amortized or whole-process cell crosses. All improvements are finite constants, not exponent changes. |
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

Stage 68 retains the current point-defined `n=53` base with 9,964 points and 94
orbit columns, selected without scalar labels or target-subgroup enumeration.
Packed witnesses reduce the support table to 738,197,504 bytes. Rank-aware
collection reaches augmented rank 95 with exactly 95 relations, and ordered,
pipelined support expansion plus fixed/fused PCLMUL arithmetic preserves every
relation hash and group check. On the same public point
`Q=(2565091273463387,5885236316843894)`, fused direct takes 6.039837 wall
seconds / 11.269810 core-seconds / 1,051,987,968 B RSS versus 5.106080 seconds
for signed-Frobenius rho: a 1.182872 loss. Fresh build plus direct is 19.474746
times rho.

Stage 72 then retains the optimized unrelated public hash target. It constructs
no scalar, derives all factor-base logs from 95 relations, and has direct IC and
rho independently recover `7892094459170` with `[d]G=Q`. Unknown-scalar direct
takes 6.571510 seconds and remains 1.467611 times rho; fresh build plus direct
is 25.359925 times rho. The corrected Linux single-core receipt pins parent and
child to CPU 0: direct takes 10.945707 wall / 8.576129 core-seconds versus rho
5.490656 wall / 4.334277 core-seconds, a 1.993515 wall and 1.978676 CPU loss.

The narrow supported conclusion is unchanged: this is strong internal
engineering and finite public toy-research evidence. The known SAT-based
point-decomposition and Frobenius-invariant-factor-base literature remain the
prior-art baseline. The finite `n=41` online descent is faster after
precomputation; rho remains faster under amortized and full available
accounting, and the `n=53` direct arm loses on the exact same target. Licensed
Magma and unaffiliated reproduction/novelty review remain open. This is not a
new Koblitz index-calculus SOTA result.
