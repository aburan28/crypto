# Untouched n37 K8/K16 one-target confirmation

Status: preregistered; no target plan, solve outcome, or new profile has been
examined for this panel. Push this protocol and `SPEC.json` before planning or
measuring. Keep any subsequent failure or protocol amendment in a new commit.

## Question and frozen inputs

The prior [six-base sweep](../ecbench_n37_rank_columns_20261004/RESULT.md)
selected K16 by the smallest *incomplete* cold counted cost, but its K16/K8
ratio was 0.990 [0.927, 1.057] on eight targets. A separate
[whole-solve instruction census](../ecbench_callgrind_solve_20261004/RESULT.md)
priced K16 and K42, leaving K8 unprofiled. Test whether K8 or K16 is the
better implementation-cost base on 16 **new** public one-target workloads of
`icv1-f2m37-tm534059-32aad96b`, exact EC1 representation
`EC1N37Ce0h0c51aa4aa7c3`, subgroup order 230603167. The same strong
signed-Frobenius rho implementation solves each point as the reference.

`SPEC.json` fixes the public target law, seed 202610044081, 16 indices,
curve, method parameters, algorithm seed schedule, five measured rounds,
one warmup, and 120-second per-child limit. The K8 and K16 factor-base
constructions and rank seed are unchanged from the prior sweep: 592 and
1,184 actual subgroup-usable points before folding, with eight and 16 matrix
columns respectively. Their existing candidate IDs are
`IC1N37Ckb0fb592PDP3mitmfrobeniuscountedRCsampleLAgaussTDpdpISO0h0eb4fc2e0d54`
and
`IC1N37Ckb0fb1184PDP3mitmfrobeniuscountedRCsampleLAgaussTDpdpISO0h44f5af6dc772`.
The `ic-k16-control` arm duplicates the complete K16 method as an A/A check.
Before measurement, `ecbench plan --spec SPEC.json --json` must confirm all
16 workload IDs and target points differ from the earlier eight-target panel;
otherwise freeze a new seed in a committed amendment before running.

## Two cost boundaries and verification

First run the native `ecbench` session into a never-reused directory. Preserve
all 384 warmup and measured records, including failure/timeout rows, host
capsule, binary and plan hashes, exact candidate manifests, scalar replay,
`ecbench verify --replay-all`, and the independent Linux full replay. Pair
each K8, K16, control, and rho run on the same public point and round. The
counted cold group-addition-equivalent ledger is an **incomplete diagnostic**:
field operations, hashing, allocation, and modular combination remain
unpriced. The session's macOS wall measurements, if made there, are L0 and
exploratory, not an online speed result. Preserve the cold preparation and
target-dependent stage counters separately.

Second, profile exactly the first measured round's 64 child inputs from that
sealed session, one per target and arm, under one release binary on Ubuntu
24.04 x86-64 with Valgrind 3.22.0. Use `ECBENCH_CALLGRIND_SOLVE=1`,
`--tool=callgrind --cache-sim=no --branch-sim=no`, and the existing
`ecbench callgrind-ir` boundary. Freeze the exact sequence numbers, workload
IDs, method IDs, algorithm seeds, and recovered scalars as `JOBS.json` before
the first profile. The counted boundary is all user-space instructions within
`methods::solve`, including factor-base/rank preparation, failed target
attempts, descent, modular recovery and scalar replay. It excludes input
reconstruction and the external independent auditor. Preserve every raw
profile part, child input/output, Valgrind log, and SHA-256 receipt; commit the
raw archive or a durable content-addressed copy plus analyzer. Check every
recovered scalar against the independently replayed session row. An error,
timeout, missing profile part, or scalar mismatch stays in the table and
invalidates a selection claim. Callgrind Ir is a simulated instruction unit,
not native wall time or a hardware-counter value.

## Predeclared decision

The primary base-choice diagnostic is K16/K8 ratio of sums in Callgrind Ir
over the 16 point pairs, with a 20,000-sample target-block bootstrap using
seed 202610044084 and a
bootstrap 95% percentile interval. Report all per-target counts, the K16
duplicate A/A maximum relative deviation, and K8/rho and K16/rho in the same
unit. Require all 64 verified profiles and an A/A maximum deviation at most
2%. Select K16 for the next n41/n53 base-size gate only if the K16/K8 interval
is wholly below 0.95; select K8 only if it is wholly above 1.05. Otherwise
retain **both** bases for that gate. The counted cold `ecbench` ratio is
reported alongside this decision but cannot override the complete
method-solve instruction result. If the A/A or verification gate fails, do
not select either base from this run; preserve the failure and diagnose it.

This panel makes no native wall-speedup claim. In particular, it cannot
replace the primary physically isolated, same-point, one-target online
IC-versus-rho measurement. No n41, n53, n83, or ECC2K-130 transfer follows
from an n37 base choice alone. The next larger-field round must carry both
bases if this decision is inconclusive.

## Build-only amendment before profiling

[CI run 37191621092](https://github.com/aburan28/crypto/actions/runs/37191621092)
stopped at `cargo test --locked` because this repository does not track
`Cargo.lock`; it produced no profile. The corrected Linux workflow builds
without `--locked` and archives its generated lockfile and SHA-256 with the
binary provenance. This amendment also writes the bootstrap seed explicitly
before any Callgrind result is observed. The frozen targets, jobs, method
parameters, cost boundary, and selection thresholds above are unchanged.
