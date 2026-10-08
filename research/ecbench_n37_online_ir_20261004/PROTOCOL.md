# Fresh n37 K8/K16 target-only instruction gate

Status: preregistered before a target plan, solve outcome, or profile for this
panel. Commit and push this protocol and `SPEC.json` before planning or
measuring. Preserve any later amendment or failed run in a separate commit.

## Question and exact workload

[The preceding cold census](../ecbench_n37_k8_k16_20261004/RESULT.md)
selected K8 for the complete method-solve instruction route: K16/K8
Callgrind Ir was 1.1903 [1.1584, 1.2218] on 16 public points. It did not
measure a target-only instruction interval. K16 averaged far fewer target
attempts and counted target group additions, so cold and online base choices
may differ. Test that distinction on 16 **new** public one-target workloads
of `icv1-f2m37-tm534059-32aad96b`, exact EC1 representation
`EC1N37Ce0h0c51aa4aa7c3`, subgroup order 230603167. The same strong
signed-Frobenius rho implementation solves each point.

`SPEC.json` freezes target law `hash_to_subgroup_v1`, seed 202610044181,
16 indices, four arms, rank seed 202610042032, algorithm seed schedule,
five measured rounds, one warmup, and a 120-second per-child cap. K8 and
K16 keep the prior factor-base and solver parameters: nominal eight and
16 folded columns, with 592 and 1,184 previously enumerated usable
points. The new native build must re-enumerate and record the actual base
sizes and exact IC1 candidate IDs; neither nominal size nor a historical
ID substitutes for that inventory. `ic-k16-control` repeats the complete
K16 method. Before measurement, commit the exact `ecbench plan` and
verify every workload ID and point differs from the preceding 16-point
panel and the older eight-point factor-base sweep. If not, amend the seed
in a committed protocol change before any solve.

The one-target boundary is after reusable base, pair table, full-rank
relation collection, verified column logs, and strong-rho jump table are
ready. The IC interval starts at Q subgroup validation, charges all direct
and shifted residual attempts, exact m3 PDP lookups, full-point witness
checks, log combination, and scalar replay, and ends when `[d]G = Q`
verifies. The five exclusive online wall phases are `target_query`,
`target_PDP`, `target_relation_check`, `target_descent`, and
`target_recovery_check`. Rho's interval starts at the first
target-dependent walk and ends after independent recovery verification.
The reusable preparation, online interval, and post-interval work remain
separate; failures and timed-out attempts stay charged to their target.

## Two measured units and a correctness gate

First, run the unchanged one-target ecbench accounting contract into a
never-reused session directory. Preserve all 384 warmup and measured
execution records, failures, host capsule, code and plan hashes, exact
candidate manifests, phase counters, full local `--replay-all` audit and
independent Linux replay. The counted group-addition-equivalent ledger
leaves field, hash, allocation, and modular work unpriced; its quotient
does not bound true IC/rho cost. A macOS session earns L0 and its wall
times remain exploratory.

Second, freeze the exact first measured round's 64 child inputs in a
committed `JOBS.json` **before** profiling. On one Ubuntu 24.04 x86-64
release binary under Valgrind 3.22.0, set both
`ECBENCH_CALLGRIND_SOLVE=1` and `ECBENCH_CALLGRIND_TARGET=1`, then profile
each frozen child using `--cache-sim=no --branch-sim=no`. The new
`ecbench callgrind-online-ir` parser must find exactly one ordered pair
of target interval markers inside the outer solve markers, reject missing
or duplicate markers, and return separate pre-target, target, and
post-target Ir whose sum equals complete solve Ir. The IC and rho markers
must bracket the same online events as their wall clocks. Preserve raw
numbered Callgrind parts, inputs, child outputs, binary and lockfile hashes,
CPU feature dispatch, Valgrind logs, and per-file SHA-256 receipt. An
error, timeout, missing part, or scalar mismatch remains visible and
invalidates a base selection.

Callgrind Ir is a **simulated user-space instruction** unit. Report
target-only and complete-solve Ir separately against the same-point rho
in a four-arm table. Do not call either an isolated native wall result,
a PMU retired-instruction count, or a primary online speedup.

## Frozen decision

For each unit use the ratio of sums over the 16 point pairs and a
20,000-sample target-block bootstrap with seed 202610044184. Report
per-point counts, 95% percentile intervals, K8/rho and K16/rho, and the
duplicate K16 A/A maximum paired relative deviation in both units.
Require 64/64 verified profiles, independent replay of every measured
native execution, and target-only A/A maximum deviation at most 2%.
If the target-only K8/K16 Ir interval is wholly above 1.10, prioritize
K16 for an **isolated n37 online wall** gate. If wholly below 0.90,
prioritize K8. Otherwise retain both without a target-only choice.
The cold Ir decision from the prior panel cannot override this gate,
and this diagnostic cannot override a future isolated online wall
measurement. Carry both bases into n41/n53 until those sizes are measured;
an n37 result alone does not establish transfer.

This panel makes no wall-speedup, n41/n53, n83, or ECC2K-130 claim.
The next primary gate requires an auditable host-level CPU-isolation
receipt and same-point online IC/rho wall intervals. Until that exists,
leave `online_speedup` unknown even if target-only Ir favors one base.
