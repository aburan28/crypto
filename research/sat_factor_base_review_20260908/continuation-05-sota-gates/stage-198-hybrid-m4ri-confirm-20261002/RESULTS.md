# Stage 198: 4,096-row hybrid M4RI confirmation

## Decision

`CONFIRMED_FOR_DEFAULT_CHANGE`. The candidate clears the frozen three-pair
median wall-and-total-core gates. It is not yet the repository default; a
separate selected-source commit and unset replay remain mandatory.

| pair | order | wall ratio | total-core ratio | RSS ratio |
|---:|:---|---:|---:|---:|
| 1 | current, hybrid | 0.754489 | 0.884893 | 0.600994 |
| 2 | hybrid, current | 0.956306 | 0.859775 | 0.630062 |
| 3 | current, hybrid | 0.908433 | 0.879279 | 0.770164 |
| **median** | frozen interleaving | **0.908433** | **0.879279** | **0.630062** |

Median wall falls 9.16 percent, total CPU falls 12.07 percent, and RSS falls
36.99 percent. All three wall and CPU pairs individually favor the hybrid.

## Correctness and work

Every control routes zero matrices through full M4RI and performs
147,794,583,858 XORs. Every candidate routes exactly 481 matrices through full
M4RI and performs 102,707,985,015 XORs. Profiling counters are zero.

All six runs authenticate the same public source and equation fingerprint,
visit all 512 fixed-X1 masks, skip 270 non-rational masks, complete all 242
rational systems, find zero roots, and return exhaustive `UNSAT`. The
algebraic factor base enumerates neither the target subgroup nor known
discrete-log labels.

The confirmation inherits the exact Stage 197 screen result at SHA-256
`84b92c9d0e2b028095b48553ce33d7f07077d1590259dfd831422acc9e8de241`.

## Verification and accounting

The Rust verifier authenticates all commands and modes, threshold 4,096, exact
0/481 routing, zero profile counters, source/equation identities, six terminal
records, pairwise ratios, median decision, the Stage 197 parent, and unchanged
default behavior. All 12 charged metrics files have one authenticated receipt.
Final replay passes `17/17`; result SHA-256 is
`38dd39b13cbcf26ac0fd6e5f43f27bf6e733488d349cc24c197c4842047c5198`.

Stage 198 contributes a measured lower bound of 12 components,
`249.727289` wall-seconds, `1,362.341551` total core-seconds, and
`5,192,499,200` bytes peak RSS. The cumulative measured campaign lower bound
is 682 components, `25,322.254375` wall-seconds, `67,684.185543`
core-seconds, and `6,310,576,128` bytes maximum RSS.

Complete campaign cost remains `null`. This is a confirmed one-target F4
engineering improvement, not relation-yield, unknown-scalar, full-rho,
independent-review, novelty, or Koblitz index-calculus SOTA evidence.
