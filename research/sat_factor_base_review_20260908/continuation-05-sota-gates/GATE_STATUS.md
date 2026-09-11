# Koblitz index-calculus SOTA gate status

The campaign now has algebraically selected, target-independent factor bases;
matched native-XOR, WDSat, CryptoMiniSat, direct MITM, GGMP, and rho controls;
a balanced 160-instance Phase-B panel through `n=59`; five public targets whose
scalars were never constructed; compact terminal custody; and a project-internal
replay of every retained mathematical witness. It has not established a new
state of the art.

| Gate | Status | Current evidence | Remaining requirement |
|:--|:--|:--|:--|
| 1. Charge every stage and resource | **Partial overall** | Stage 20 seals Phase A, staging, tool builds, the 160-instance Phase-B run, scoring, failed attempts, and non-additivity rules. Its successful path charged 17,864.202394 core-seconds and 19,682.452562 summed receipt wall-seconds. Stage 21 separately charges natural relation-yield discovery. Stage 23 charges a fresh build, two public discoveries, target generation, five IC/rho pairs, and an inclusive outer driver: 669.627462 child core-seconds and 671.154203 enclosing core-seconds. Stage 25 repeats the production path under one-CPU Linux affinity and records 739.0632245 seconds of true single-core elapsed time. | Add the licensed Magma return. A literal full-cost closure also needs prior source/dependency acquisition, toolchain and system-dependency installation, licensing costs, the same affinity measurement for the full Stage-20 backend matrix, and aggregate simultaneous process-tree memory; existing frozen values cannot be relabeled to provide them. |
| 2. WDSat, CryptoMiniSat, Magma F4, MITM, GGMP | **Partial** | Native XOR, WDSat, CryptoMiniSat, direct MITM, and standard/GGMP cells have bounded, target-matched evidence. Stage 22 packages all 160 exact blinded Magma inputs and validators in a checksum-pinned licensed-host packet. | No licensed Stage-22 Magma process has run. Execute the packet once under its one-thread, watchdog, no-retry contract; seal the return before opening truth labels; then score it. Public-calculator basis terminals do not supply the missing process CPU/RSS or all-cell witness evidence. |
| 3. Single-core, core-seconds, memory, conflicts, wall | **Partial; affinity-bound degree-23 control complete** | Executed local processes retain wall, user/system CPU, total core-seconds, peak process RSS, conflicts when exposed, and outer receipts. Stage 23 records 28,422,672 conflicts, 511.091332 IC core-seconds, 0.282682 rho core-seconds, and 1,194,147,840-byte maximum process RSS. Stage 25 restricts the complete five-target Linux run to CPU 0 and measures 739.0632245 seconds of true one-CPU outer elapsed time, 725.634561 total outer core-seconds, and 1,339,957,248-byte peak process RSS. | The historical `single_core_seconds` field remains an aggregate-CPU alias. Affinity-bound time for the full Stage-20 `n=31`/`n=41`/`n=59` backend matrix, simultaneous aggregate process-tree RSS, and licensed Magma process receipts remain absent. |
| 4. Scale through n=31, n=41, and a larger PDP regime | **Satisfied for finite execution coverage; no scaling law** | Stage 20 ran 40 balanced targets in each of `n=31` standard, `n=31` GGMP, `n=41` standard, and `n=59` standard. Across 480 native/WDSat/CMS outcomes it recorded 261 correct, 219 inconclusive, and zero false positives or false negatives. Standard `n=31` and `n=41` were fully resolved; GGMP `n=31` and `n=59` were heavily censored by caps. | Complete the same-input Magma arm. A complexity exponent or crossover claim needs a separately frozen same-family replicated design; these four cells do not justify one. |
| 5. Unknown scalar with no constructed factor-base logs | **Satisfied for the finite degree-23 panel** | Stage 23 generated five targets by public hash-to-affine decoding and cofactor projection without constructing target scalars. Both public discoveries ran before target generation without subgroup enumeration, relation yield, solver timing, or factor-base log labels. Five IC/rho pairs completed. Stage 24 independently rebuilt the field, curves, target stream, 4,281-point factor base, 93 projected columns, all 437 attempt targets, all 252 witnesses/rows, five modular scalars, and all ten `[d]G=Q` equations. | The narrow retained-witness replay is complete. SAT proof/conflict trajectories and the hidden step-by-step rho trajectory were not retained, so the broader mathematical-payload flag remains false. Larger regimes and unaffiliated replay remain open. |
| 6. Automorphism-optimized Pollard rho | **Satisfied as a finite same-target control; crossover open** | On the Stage-23 targets, IC used 511.091332 core-seconds and rho used 0.282682. The online ratio is 1808.008051 in rho's favor; charging both public discoveries raises it to 1810.596058. The Stage-25 one-CPU Linux repeat used 589.051224 IC core-seconds and 0.354414 rho core-seconds, ratios 1662.042764 and 1664.443027 in rho's favor. Rho was faster in all five pairs. | No larger-regime crossover is established. Any such claim needs a newly frozen same-target design with discovery, setup, relation collection, linear algebra, validation, and whole-process accounting. |
| 7. Independent external reproduction and novelty review | **Open** | Project-authored GitHub runs reproduce the workflows on fresh Linux hosts. Issue 97 contains the exact review request, Stage-22 licensed packet, source pins, and `CONCUR` / `QUALIFIED` / `BREAKS` format. The Stage-23 bundle is portable and structurally verified, and Stage 24 supplies one project-internal retained-witness replay with a separately locked official BLAKE3 helper. | A reviewer unaffiliated with the project must return a sealed reproduction and source-pinned novelty/correctness assessment. Project-authored CI, local independent code, and an empty review template do not satisfy this gate. |

## Current finite results

Stage 20's balanced panel contains 20 decomposable and 20 nondecomposable
public-synthetic targets in each of four cells. Standard `n=31` and `n=41`
completed under all three SAT backends. In `n=31` GGMP, CryptoMiniSat resolved
17/40, native XOR 2/40, and WDSat 1/40; the rest were inconclusive. At `n=59`,
WDSat resolved 1/40 and all other SAT outcomes were inconclusive under their
caps. The absence of false classifications establishes bounded correctness, not
an advantage over direct decomposition or rho.

Stage 21 measured target-independent natural relation yield for the selected
`K_0 / F_(2^23)` GGMP divisor base `[0,2]`: 163 of 256 natural targets
admitted direct pair decompositions, a rate of 0.63671875 with Wilson 95%
interval `[0.576185047549112, 0.6932099917094386]`. The planted control was
64/64 and the exact-miss control 0/64. This remains a finite bridge measurement.

Stage 23 strengthens the unknown-scalar boundary. Its five target systems have
94 columns but ranks `[43,60,40,68,41]`; nullities affect factor-base logs while
the target coordinate is invariant. Early target-coordinate recovery reduced
the preserved diagnostic's 470 relations, 832 trials, 54,514,950 conflicts,
and 1,023.366665 IC core-seconds to 252 relations, 437 trials, 28,422,672
conflicts, and 511.091332 IC core-seconds. These are reductions of 46.38%,
47.48%, 47.86%, and 50.06%.

Stage 25 supplies the previously missing literal single-core definition for
that degree-23 panel. Linux restricted the parent and fresh child probes to CPU
0 before and after the run. The inclusive Stage-23 outer wall was 739.0632245
seconds, with 725.634561 outer core-seconds and 98.1830% utilization of the
available core. This is a host-specific finite control; it does not backfill
single-core time for the Stage-20 backend matrix.

The compact successor-03 evidence archive is pinned by SHA-256
`7886a767b0fe0f3372bd034429e3616632c40ea1781de835602fa08b72e9dd44`.
Its trusted structural verifier reconstructs 14 tasks and five rows without
executing bundled code. Stage 24 then independently validates every retained
positive mathematical witness. Neither step reconstructs missing SAT proof
traces or hidden rho states, and neither is an external review.

The narrow supported conclusion remains: this is a strong internal engineering
and finite toy-research improvement. Automorphism-optimized rho remains far
faster on every matched target. The licensed Magma, literal single-core/parallel
memory, larger-crossover, and unaffiliated reproduction/novelty gates remain
open. This is not a new Koblitz index-calculus SOTA result.
