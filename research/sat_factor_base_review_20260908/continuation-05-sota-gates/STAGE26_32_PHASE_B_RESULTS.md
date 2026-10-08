# Phase B one-CPU matrix and WDSat capacity correction

Date: 2026-09-11

## Claim

This is a complete internal execution of the truth-free 160-instance Phase B packet on native XOR SAT, WDSat, CryptoMiniSat, and matched direct meet in the middle. It includes standard factor-base cells at `n=31` and `n=41`, the GGMP cell at `n=31`, and a larger `n=59` PDP cell. Two `n=59` WDSat buffer assertions in the immutable Stage 26 run were rerun in Stage 32 with the exact required capacity and ended as clean timeouts.

The result is a strong internal engineering and public synthetic toy-research improvement. It is **not** a Koblitz index-calculus state of the art result. Licensed Magma F4, unified end-to-end scaling at `n=41` and above, fully comprehensive cost capture, and unaffiliated reproduction and novelty review remain open.

## Algebraic factor-base boundary

The factor bases and target predicates are generated from public algebraic curve parameters and coordinate-subspace or divisor recipes. The target subgroup is not enumerated to discover the factor base, target discrete logarithms are absent from the solver packet, and factor-base logarithm labels are not constructed for discovery.

The degree-23 unknown-scalar control selected the divisor-coordinate recipe `[0,2]`, producing 4,096 abscissae, 4,281 rational points, 95 signed orbits, and 93 projected signed orbits. The degree-31 public-target control selected the two-torsion saturation of divisor `[0,1,5]`, producing 1,986 abscissae, 3,971 points, 66 signed orbits, and 32 projected columns. In the latter control, all 32 factor-base logs were derived from collected relations and certified by group multiplication.

## Same-instance solver matrix

The 160 inputs contain 20 decomposable and 20 nondecomposable targets in each of four cells. The SAT classifications below are taken only after authenticated truth scoring. No false positive or false negative occurred.

| Cell | Native XOR SAT | WDSat | CryptoMiniSat | Direct MITM |
|---|---:|---:|---:|---:|
| `n31-l5-m3-standard-a1-f0` | 20 TP, 20 TN | 20 TP, 20 TN | 20 TP, 20 TN | 20 TP, 20 TN |
| `n41-l5-m3-standard-a1-f0` | 20 TP, 20 TN | 20 TP, 20 TN | 20 TP, 20 TN | 20 TP, 20 TN |
| `n31-l5-m3-ggmp-a0-f0` | 2 TP, 38 inconclusive | 1 TP, 39 inconclusive | 18 TP, 22 inconclusive | 20 TP, 20 TN |
| `n59-l9-m3-standard-a1-f0` | 40 inconclusive | 40 timeout-inconclusive | 40 timeout-inconclusive | 20 TP, 20 TN |

Across the three SAT backends, the matrix produced 141 true positives, 120 true negatives, and 219 inconclusive outcomes. Direct MITM produced 80 true positives and 80 true negatives. Stage 32 changed two WDSat rows from `solver_error` to `timeout_inconclusive`. Both rows were already inconclusive under truth scoring, so the 141/120/219 totals are unchanged. The composed score has zero solver-error rows.

## SAT resources

The process-wall column sums individual solver-process wall time. Peak RSS is the maximum individual solver-process RSS. Conflict totals include only solver rows that expose a conflict counter.

| Backend | Total core-seconds | Summed process wall (s) | Peak RSS (bytes) | Conflicts | Rows reporting conflicts |
|---|---:|---:|---:|---:|---:|
| Native XOR SAT | 958.931674 | 1,278.114043 | 140,333,056 | 9,233,732 | 160 |
| WDSat, corrected panel | 7,404.749336 | 9,524.938771 | 24,252,416 | 7,599,249 | 81 |
| CryptoMiniSat | 6,237.273520 | 8,137.581155 | 122,703,872 | 19,480,743 | 98 |

Each Stage 26 cell ran under singleton CPU affinity. Its outer envelope includes packet verification, setup, source regeneration checks, and all three solvers.

| Cell | Single-core elapsed (s) | Outer core-seconds | Sampled process-tree peak RSS (bytes) |
|---|---:|---:|---:|
| `n31` standard | 343.101152 | 342.765848 | 69,894,144 |
| `n41` standard | 253.859726 | 253.559622 | 71,319,552 |
| `n31` GGMP | 8,238.613331 | 8,237.029054 | 161,984,512 |
| `n59` original immutable cell | 10,582.727936 | 10,580.835720 | 188,895,232 |
| Stage 32 two-row correction | 424.290467 | 424.217075 | 72,146,944 |
| `n59` original plus correction | 11,007.018404 | 11,005.052795 | 188,895,232 |

The Stage 26 workflow wall time was 11,088 seconds. The Stage 32 workflow wall time was 457 seconds, giving 11,545 seconds when the two workflow durations are summed. The four Stage 26 cells overlapped over a 10,750-second parallel span.

The Stage 32 increment charges 1.04 core-seconds for the metered WDSat source clone, 1.460936 core-seconds for the clean build, and 424.217075 outer core-seconds for packet/source verification and both attempts: 426.718011 core-seconds in total. Its two WDSat processes each reached the 120-second watchdog. Required ANF buffer capacities were 32,588 and 32,804; the fresh binary was compiled at the exact maximum, 32,804.

## Direct meet in the middle

The direct method solved all 160 same-instance inputs with 80 true positives and 80 true negatives. It consumed 163.312640 child-process core-seconds, 759.710741 inclusive outer core-seconds, 766.538655 summed one-CPU elapsed seconds, and 322 seconds of workflow wall time. Its maximum sampled process-tree RSS was 85,479,424 bytes. It performed 4,804,252 group additions and materialized 4,768,800 pair-table entries. SAT conflicts do not apply to this construction.

## Charged scope and omitted scope

The composed Phase B matrix charge is 21,038.596137 core-seconds. It includes exact solver acquisition and builds, four one-CPU Stage 26 cell envelopes, matched direct-MITM outer envelopes, and the Stage 32 correction acquisition, build, and one-CPU outer envelope. The largest recorded individual-process or sampled-process-tree memory value in that composition is 1,526,329,344 bytes, from the Rust build.

This number does not absorb the separate unknown-scalar controls below. Their curves and workloads differ, so adding them to the 160-instance PDP matrix would not form a meaningful same-instance comparison. Workflow checkout, artifact transfer, Stage 32's post-clone Git checkout, and aggregate parallel-compiler memory are not fully metered. Magma has not run.

## Unknown-scalar end-to-end controls and rho

The degree-23 single-CPU control completed five public synthetic unknown-scalar targets. Factor-base logs and target scalars were not known by construction. It collected 252 verified relations in 437 attempts across 437 SAT calls and charged 28,422,672 conflicts. The outer envelope consumed 725.634561 core-seconds, 739.063225 single-core seconds, and 1,339,957,248 bytes peak RSS. The algorithmic IC work consumed 589.051224 core-seconds, compared with 0.354414 core-seconds for signed-Frobenius/negation Pollard rho on the same five targets: ratios of 1,662.04 online and 1,664.44 with setup charged.

The degree-31 control completed five public hash-derived targets without constructed target scalars or factor-base log labels. The successful factor-base search, build, collection, sparse linear algebra, descent, and verification path consumed 6,214.831656 core-seconds and 1,958.164505 sequential wall-seconds. Including the bounded failed predecessor raises this to 8,259.024390 core-seconds and 2,858.175432 wall-seconds. It collected 11,312 relations in 16,384 probes, retained 11,256 unique verified rows, solved a 32-column system, certified all 32 factor-base logs, and completed all five descents. The online IC/rho wall ratio was 13.0256; the amortized full path/rho ratio was 844.4785. Pollard rho remains decisively faster in both finite controls.

## Gate status

| Gate | Status | Evidence still missing |
|---|---|---|
| 1. Charge discovery through parallel resources | Partial | One unified same-instance end-to-end experiment at `n=41` and above; some orchestration and aggregate compiler costs |
| 2. WDSat, CryptoMiniSat, Magma F4, direct MITM, GGMP | Partial | Licensed Magma F4 on all 160 exact inputs |
| 3. Single-core, core-seconds, memory, conflicts, wall | Partial | Magma resources and the unmetered workflow portions listed above |
| 4. Scale through `n=31`, `n=41`, and larger | PDP complete through `n=59`; end-to-end IC incomplete | End-to-end IC at `n=41` and a larger regime |
| 5. Unknown-scalar runs | Complete at degrees 23 and 31 | Unknown-scalar end-to-end runs at `n=41` and above |
| 6. Full cost versus automorphism rho | Complete at degrees 23 and 31; broader gate partial | Unified `n=41` and larger comparisons |
| 7. External reproduction and novelty review | Missing | Unaffiliated execution and source-pinned novelty verdict |

## Custody

- Stage 26 workflow: `34632018379`, execution commit `03968a5da2a511abb723652529de34dc60e10942`
- Stage 26 terminal archive SHA-256: `7465181f8c53d057829e883ff0356dd36ce6194edc57777a03085b4703158a20`
- Stage 32 workflow: `34649735091`, execution commit `8ad01b0c8eb7bd90cb83976a080cea5cd518f0fb`
- Stage 32 GitHub artifact: `10283915853`, digest `sha256:9eacd0472395b475129e652658ba4917843f60b9b2d17386903f244e61f9843c`
- Stage 32 committed archive SHA-256: `04a0ebbc099b08bff526df6b0f9c428980c1e21b8e578d2314ede6652dc317e6`
- Composed score SHA-256: `ae7cc2f59a04e153f59ad1de65704c233f6fa8960976dd95593abdd8df9d1897`
- Composed score inventory SHA-256: `a86226326b54e697e4e530b80d218e60d049c37dde2ff9f296b918a4671ea5ef`

The composition verifier extracts both archives, verifies their seals, compares the two corrected source systems and targets byte for byte against the immutable Stage 26 tasks, checks the exact capacity requirements and timeout terminals, recomputes all corrected rows and resource totals, and requires every SOTA gate flag to remain false.
