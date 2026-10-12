# Stage 192: selected parallel F4 and direct-MITM boundary

## Result

The selected twelve-worker Stage 190 path remains far slower than same-binary
direct meet-in-the-middle decomposition on the identical public
`n=59, ell=9, m=3` true-negative.

The fixed fresh-process order was direct MITM, selected parallel F4, direct
MITM. The direct reference is the arithmetic mean of both bracketing controls.

| arm | wall seconds | total core-seconds | single-core seconds | peak RSS |
|:---|---:|---:|---:|---:|
| direct MITM before | 0.584046 | 0.302350 | 0.302350 | 41,222,144 B |
| direct MITM after | 0.325632 | 0.310483 | 0.310483 | 40,910,848 B |
| **direct arithmetic-mean reference** | **0.454839** | **0.306417** | **0.306417** | **41,066,496 B** |
| **selected parallel F4** | **16.765098** | **181.297986** | **null** | **4,366,270,464 B** |

The selected parallel F4/direct ratios are:

- wall: `36.859390`;
- total core-seconds: `591.671747`;
- single-core: `null`; and
- peak RSS: `106.321963`.

Parallelism reduces F4 latency relative to Stage 191's strict single-core row,
but raises total CPU and memory substantially. It does not approach direct MITM.

## Correctness and accounting

Both direct controls authenticate the same source and repeat 483 rational
factor-base points, 116,635 pair entries, and 117,369 group additions. The F4
arm uses twelve Rayon workers, batch 512, dense pair selection, the selected
unset build-serial policy, and parallel `BlockTables` elimination. It visits all
512 masks, skips 270 non-rational masks, completes all 242 systems, finds zero
roots, and returns exhaustive `UNSAT`.

Both approaches use the algebraic `span_F2(1,z,...,z^8)` factor base without
target-subgroup enumeration or known log labels. Conflicts remain `null` or
unavailable because neither method exposes a SAT-conflict counter.

The first composition correctly replayed all Stage 192 measurements but used
Stage 190 rather than Stage 191 as the inherited cumulative base. That result
is retained under `development/superseded-cumulative/`; the additive correction
changes only inherited totals. Final native replay passes `24/24` with result
SHA-256 `e1849a1b2ae97bdc3633025eaadcb18438ff2a49448820a936a890db24904a48`.

Stage 192 contributes a measured lower bound of 10 components,
`196.127459` wall-seconds, `961.170802` total core-seconds, and
`5,648,957,440` bytes peak RSS. The cumulative measured campaign lower bound is
582 components, `23,910.746621` wall-seconds, `61,559.696526` core-seconds,
and `6,310,576,128` bytes maximum RSS.

Complete campaign cost remains `null`; development and final-write work is not
all outer-metered. This is a parallel decomposition-stage comparison, not
natural relation yield, relation collection, linear algebra, a recovered DLP,
full rho, external reproduction, novelty, or Koblitz index-calculus SOTA.
