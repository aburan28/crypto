# Koblitz index-calculus SOTA gate status

Current through Stage 36, 2026-09-12. The machine-readable audit is
`stage-36-current-gate-audit-20260912/audit.json`.

The campaign has target-independent algebraic factor bases; matched native-XOR,
WDSat, CryptoMiniSat, direct-MITM, GGMP, and signed-Frobenius-rho controls; a
balanced 160-instance PDP panel through `n=59`; and public unknown-scalar
end-to-end runs at degrees 23, 31, and 41. It has not passed all seven gates and
does not establish a new state of the art.

| Gate | Status | Current evidence | Remaining requirement |
|:--|:--|:--|:--|
| 1. Charge every stage and resource | **Partial overall** | The corrected Phase-B matrix charges 21,038.596137 core-seconds across acquisition, builds, four one-CPU SAT cells, direct MITM, and the WDSat correction. Stage 35 separately charges fresh Cargo acquisition, a sampled four-job clean build, algebraic factor-base discovery/materialization, 14,336 relation probes, sparse linear algebra, five walked descents, rho, and the one-CPU outer envelope: 308.810720 core-seconds, 105.117902 sequential wall-seconds, and 1,642,135,552 bytes maximum sampled tree RSS. | Licensed Magma process resources are absent. Stage 35 excludes acquisition of the preinstalled operating system and Rust toolchain. |
| 2. WDSat, CryptoMiniSat, Magma F4, MITM, GGMP | **Partial** | Native XOR SAT, WDSat, CryptoMiniSat, and direct MITM ran on the exact 160-input packet. Standard `n=31`/`n=41` and GGMP `n=31` are represented. The Stage-32 successor removes the two original WDSat buffer errors without rewriting Stage 26. | Execute all 160 Stage-22 Magma inputs on a licensed host under the frozen one-thread/no-retry contract, seal the return before truth scoring, and report F4 resources. |
| 3. Single-core, core-seconds, memory, conflicts, wall | **Partial because Magma is absent** | Every executed Phase-B arm retains process wall, total core-seconds, peak RSS, conflict or operation counters, outer process-tree memory, and workflow wall. Stages 33 and 35 add true inherited-one-CPU end-to-end `n=41` receipts. | Supply the same fields for licensed Magma F4. |
| 4. Scale through `n=31`, `n=41`, and a larger PDP regime | **Satisfied for finite execution coverage** | Phase B covers standard `n=31`, GGMP `n=31`, standard `n=41`, and standard `n=59`. End-to-end public unknown-scalar controls now cover `n=31` and `n=41`. | The evidence is finite and toy-sized; it is not an asymptotic scaling law. The `n=59` arm remains PDP-only. |
| 5. Unknown scalar with no constructed factor-base logs | **Satisfied for finite degrees 23, 31, and 41** | Each degree has five public targets whose expected scalars were not constructed. The degree-31 and degree-41 workflows derived factor-base logs from verified relations and group-certified every column. Stage 35 keeps the Stage-33 algebraic search and exact targets while changing only individual descent. | Unaffiliated replay remains part of gate 7. No deployed-size conclusion follows. |
| 6. Full cost against automorphism-optimized Pollard rho | **Finite online crossover; amortized and full-cost crossover false** | Stage 35 uses signed-Frobenius classes (`A=2n`). After the log database exists, the algebraic walked descent is 2.048767 times faster than rho and improves the Stage-33 online ratio by 4.608397 times. Discovery, relations, logs, and descent amortized over five targets remain 28.857106 times slower than rho. Fresh build plus science consumes 105.117902 seconds, 208.529527 times the same five rho runs. | The online crossover assumes a precomputed database. Preinstalled OS/toolchain acquisition is excluded. The result is finite and changes a constant, not an exponent. |
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

## Current n=41 end-to-end result

Stage 35 compares the algebraic Frobenius union generated by masks
`[1,2,4,8,16,32]` with its two-torsion saturation on 1,024 separately sampled
subgroup points. It does not enumerate the subgroup or receive the five future
targets. It selects the saturated base: 2,380 abscissae, 4,759 rational points,
60 signed orbits, and 29 projected columns.

Collection retains 39 verified relations from 14,336 probes and group
multiplication certifies all 29 logs. Collection uses three summands; individual
descent walks probes and asks for two. Five domain-separated public hash targets
solve in 2,417,575 walked probes and 0.246046 seconds total, versus 0.504091
seconds for signed-Frobenius rho. The scientific outer envelope is 15.831464
seconds on one inherited CPU, 15.802432 core-seconds, and 345,108,480 bytes
sampled process-tree peak.

The narrow supported conclusion is unchanged: this is strong internal
engineering and finite public toy-research evidence. The known SAT-based
point-decomposition and Frobenius-invariant-factor-base literature remain the
prior-art baseline. The finite online descent is faster after precomputation;
rho remains much faster under amortized and full available accounting. Licensed
Magma and unaffiliated reproduction/novelty review remain open. This is not a
new Koblitz index-calculus SOTA result.
