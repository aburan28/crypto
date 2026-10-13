# Claim round: the root-finding oracle on a Frobenius-orbit base, so that a `vs_rho` claim can be built

Registered 2026-10-13, after the ladder of `../` and before any session
of this directory ran.

## Why

Every `vs_rho` claim build in the ladder was refused: the IC1 candidate
identity names a factor base by its sign-and-Frobenius orbits, and the
standard subspace base (kept for comparability with the F6-IC ladder) has
more columns than orbits.  A base whose columns *are* the signed Frobenius
orbits is `koblitz-orbit`: the Frobenius-invariant subspace cut out by
factors of the Frobenius characteristic polynomial, each orbit folded onto
one column.  This round reruns the root-finder on that base so that the
claim machinery can name the candidate, and builds a claim for every
workload.

## Sizes and bases

Only sizes where an invariant subspace of useful dimension exists (the
dimensions are sums of the degrees of the factors of `x^n − 1`):

| `n` | `icv1` slug | base | points | columns |
|--:|:--|:--|--:|--:|
| 23 | `icv1-f2m23-t5197-69e76b73` | `koblitz-orbit:divisor=1` | 2,025 | 45 |
| 31 | `icv1-f2m31-tm90707-c95f16f5` | `koblitz-orbit:divisor=2;3` | 869 | 15 |

(`ecbench fb` on the binary that runs the sessions; `n = 13, 19, 37` have
no invariant subspace between dimension 1 and `n − 1` worth using, and
`E_0` at `n = 17` has no usable subgroup in the constructor.)

## Arms

`rho-strong` (reference, `rho.signed_frobenius_strong`, lanes 32,
`dp_bits` 4, `step_cap_factor` 2000), `ic-rootfind` (candidate,
`ic.pipeline`, the base above, `rootfind-mitm:m=3`), `ic-mitm` (baseline,
same base, `mitm:m=3,negation_folded=1`), `rho-strong-aa` (control).
Eight public one-target workloads a size, three measured rounds, `warmup = 1`,
`order = alternate`, `isolation_required = L2`, `max_trials = 200000`.

## What is read

The audit (`verify --replay-all`), the paired comparisons, and one
`ecbench claim build --ic ic-rootfind --rho rho-strong` per workload (first
measured round), checked by `ecbench claim check`.  Without a replay
receipt from a host of another class, the checker fails a claim on exactly
`independent_validation` and the two certificates; that is reported as the
claims' status, not hidden.  Prediction: every claim's speed comparison
shows index calculus slower than the strong rho (the ladder's reading),
and every claim is otherwise well formed.

## Amendment 1 (before the rerun)

The first run of this round (kept in `superseded-pre-merge/`) used a
binary built before `origin/main` was merged into this branch.  Its 16
claims built, but every one lacked the `operation_accounting` and
`record_class` stage fields: that binary's `ecbench claim build` did not
write them.  Main's commit `2bcca4148` ("Emit truthful claim classification
and native counter accounting") makes the builder write both.  A claim is
built only by the binary that measured its session, so the round is rerun
whole, unchanged in specs, seeds and arms, with the merged binary, and the
claims rebuilt from the new sessions.  The first run's readings
(`online_speedup` 0.097–0.392 at `n = 23`, 0.007–0.020 at `n = 31`; audits
`OK`) were seen before this amendment; no decision rule changes.
