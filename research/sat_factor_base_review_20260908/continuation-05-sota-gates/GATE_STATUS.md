# Koblitz index-calculus SOTA gate status

The campaign has advanced from a scalar-labelled degree-19 factor base to algebraically defined factor bases, matched PDP exports, two complete scalar-blind toy runs, and a same-target signed-Frobenius rho control. It has not established a new state of the art.

| Gate | Status | Current evidence | Remaining requirement |
|:--|:--|:--|:--|
| 1. Charge every stage and resource | Partial | Stage 1 charges instance generation, per-cell WDSat builds, solver processes and retained internal timers. Stage 2 charges curve/subgroup setup, factor-base construction, relation collection and linear algebra end to end. | Isolate native PDP CPU/RSS from direct MITM, repeat full end-to-end runs beyond degree 9, and charge any parallel speculative work. |
| 2. WDSat, CryptoMiniSat, Magma F4, MITM, GGMP | Partial | WDSat and CryptoMiniSat consume matched exports from one source ANF; MITM and the GGMP degree-31 construction ran. | Magma is unavailable. The standard and GGMP factor-base cells use different targets, so a target-matched causal base comparison remains. |
| 3. Single-core, core-seconds, memory, conflicts, wall | Partial | Fresh-process CPU/RSS/wall receipts exist for external solvers and every complete degree-9 attack. Conflict counts and internal stage timers are retained. | Split native PDP and MITM into separate metered processes and repeat enough samples for distributions. |
| 4. Scale through n=31, n=41 and larger PDP | Passed for the planted-PDP sub-gate | Source-validated SAT models at n=31 and n=41; n=59 ran to bounded native/WDSat/CMS outcomes and direct MITM found the planted witness. | Random-target SAT/UNSAT panels, replications, and a fitted scaling law. n=67 currently exceeds the n<64 field-bitmask implementation. |
| 5. Unknown scalar with no constructed factor-base logs | Passed at toy degree 9 | Two complete runs recovered secrets 53 and 101 from five relations in five trials each. The GGMP predicate is target-independent, no base logs are constructed, no direct shortcut occurred, and the recovered scalars reproduce the targets. | Larger degrees and independent curves. |
| 6. Automorphism-optimized Pollard rho | Passed as a toy implementation control | The signed-Frobenius quotient walk recovered both degree-9 targets; it used fewer core-seconds than index calculus in both pairs. | Independent audit, repeated distributions, larger degrees, and the full crossover comparison. |
| 7. External reproduction and novelty review | Open | An internal independent review is being obtained from a separate agent session. | Reproduction by an external person or host plus an external, source-pinned novelty review. Internal review cannot satisfy this gate. |

The degree-59 result is especially constraining: direct MITM found the planted decomposition in 3.13 seconds, while native SAT stopped at 100,000 conflicts and both WDSat and CryptoMiniSat reached 120-second watchdogs. Those capped solver outcomes are inconclusive, but they do not support a SAT advantage.

The narrow supported conclusion remains: this is a strong internal engineering and toy-research improvement. It is not a new Koblitz index-calculus SOTA result.
