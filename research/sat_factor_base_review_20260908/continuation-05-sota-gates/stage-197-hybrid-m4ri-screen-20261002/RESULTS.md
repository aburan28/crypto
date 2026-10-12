# Stage 197: 4,096-row hybrid M4RI screen

## Decision

`CONTINUE_TO_CONFIRMATION`. The unprofiled hybrid clears the frozen
strict-below-`0.97` wall-and-total-core screen gates. It is not yet selected
as the repository default.

| arm | wall seconds | total core-seconds | peak RSS |
|:---|---:|---:|---:|
| current five-column BlockTables | 38.813683 | 190.758951 | 4,065,132,544 B |
| full M4RI only at 4,096 or more rows | 29.680804 | 168.655318 | 3,379,003,392 B |
| **candidate / current** | **0.764699** | **0.884128** | **0.831216** |

The control wall time contains visible host descheduling, but the candidate
also reduces the more stable total CPU by 11.59 percent. The frozen rule uses
both fields and therefore admits confirmation.

## Work and correctness

The control routes zero matrices through full M4RI. The candidate routes
exactly the 481 matrices identified by Stage 196 and keeps every smaller matrix
on current BlockTables.

| counter | current | hybrid |
|:---|---:|---:|
| row-equivalent logical XORs | 319,313,687,585 | 318,703,372,596 |
| actually performed XORs | 147,794,583,858 | 102,707,985,015 |
| full-M4RI matrices | 0 | 481 |
| full-M4RI blocks | 0 | 351,164 |

Both arms authenticate the same source and equation fingerprint, visit all 512
fixed-X1 masks, skip 270 non-rational masks, complete all 242 rational systems,
find zero roots, and return exhaustive `UNSAT`. The algebraic factor base
enumerates neither the target subgroup nor known discrete-log labels.
Profiling counters are zero in both runs.

The threshold support remains opt-in. Unset behavior and the selected runtime
are unchanged. The exact implementation is preserved in
`candidate-hybrid-m4ri-4096.patch`.

## Verification and accounting

The Rust verifier authenticates commands, explicit modes, threshold 4,096,
zero profiling, exact full-M4RI routing, source/equation identities, terminal
correctness, ratios, decision, candidate patch, and unchanged default. All 12
charged metrics files have exactly one authenticated receipt. Final replay
passes `19/19`; result SHA-256 is
`84b92c9d0e2b028095b48553ce33d7f07077d1590259dfd831422acc9e8de241`.

Stage 197 contributes a measured lower bound of 12 components,
`265.834501` wall-seconds, `1,101.564255` total core-seconds, and
`5,246,697,472` bytes peak RSS. The cumulative measured campaign lower bound
is 670 components, `25,072.527086` wall-seconds, `66,321.843992`
core-seconds, and `6,310,576,128` bytes maximum RSS.

Complete campaign cost remains `null`. A separate three-pair confirmation in
the already frozen interleaving is required before any default change. This
one-target screen changes no relation-yield, unknown-scalar, full-rho,
independent-review, novelty, or Koblitz index-calculus SOTA gate.
