# n37 shared-rank index calculus against strong signed-Frobenius rho

Status: preregistered before planning or solving the new public target. This
follows the target-blind rank and 16-point correctness gates in PRs #1302 and
#1309. It tests one previously unseen target; the 16-point control is not a
speed baseline.

## Question and fixed inputs

On `icv1-f2m37-tm534059-32aad96b`, compare the verified online solve of one
public point by the #1309 shared-rank method and `rho.signed_frobenius_strong`.
The point is `ecbench`'s `hash_to_subgroup_v1` target at seed
`202610041751`, index 0, on Koblitz `a=0,n=37`. No producer receives a
known scalar. The target must be outside the n37 signed-Frobenius orbit
inventory frozen in #1309; if it is not, record the collision and stop this
gate without replacing the point.

The IC arm uses `compact-orbit-scan:columns=42,raw_x_cap=1000000`, the exact
counted signed-Frobenius m3 table, target-blind uniform `[a]G` rank queries,
rank seed `202610031137`, at most 1,000,000 rank trials, incremental Gaussian
elimination to all 42 columns, independent point checks of every column log,
and at most 64 target residual attempts. Attempt zero tests Q; later attempts
use `Q+[a]G` with a run-seeded RNG. The rank transcript must reproduce the
#1309 base digest
`8460ac4c28515db701c3897a03b4ce0f28abf7cd56fd759ad095f98436a76dcf`.
The rho reference uses one target, 32 lockstep lanes, 8 distinguished-point
bits and a step cap factor of 2000; its table is discarded after each run.

The attached `ecbench.spec/v1` fixes the above algorithms, one public
workload, one warm-up and five measured rounds with alternating arm order.
It includes a repeated IC control arm for A/A noise. Each arm gets the
same algorithm seed in a round. Freeze the spec and implementation in the PR
before running the session. Do not change the point, seeds, budgets or
accounting after observing outcomes. Preserve every timeout and failure.

## Measurement and decision

The primary metric is verified one-target online wall time, from the first
target-dependent operation after rank/base/table setup through `[d]G=Q`.
IC records exclusive target query, PDP (all hits and misses), witness check,
descent and scalar replay phases; rho records walk, collision and recovery.
The five IC times must sum exactly to the online interval. Pair arms on the
same Q and resource envelope. Set an online speedup only when both answers
verify, both runs earn L2 isolation and an independent other-host replay of
the exact session is available. One target cannot establish a population
claim; five repetitions measure run noise, not target variance.

Report cold factor-base construction, table, relation search, matrix solve,
column checks and target work separately, together with counted group-addition
equivalents divided by `sqrt(r)`. The derived generic floor is
`sqrt(pi/(2*74))` in S for signed-Frobenius rho. Native field arithmetic,
hashing, allocation and modular combination remain unpriced until calibrated;
any S ratio with those missing is a lower-bound diagnostic, not a speedup.
Do not infer an ECC2K-130 crossover from n37.

Correctness passes only if the fixed public Q has a verified log from both
methods in every measured round, the rank/base identities match, phase sums
and group counts replay, and `ecbench verify --replay-all` passes. Otherwise
record the failing row and stop the comparison. The performance hypothesis is
that IC online time is lower than rho outside the IC A/A interval on an L2
host; L0 local runs cannot establish it. An other-class audit receipt is
required before building an admissible `vs_rho` claim. If this gate remains
bounded or L0, the decision is explicitly inconclusive and the next work is
native-work calibration and an isolated Linux session, not a speed claim.
