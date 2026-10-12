# Frozen n37 shared-rank factor-base column sweep

Status: preregistered, no outcomes measured. This round varies the number of
signed-Frobenius factor-base columns in the native `ic.shared_rank` pipeline.
It follows the one-target same-point comparison in PR #1313 and the independent
replay in PR #1316. The source implementation is the unmodified Rust tree at
`694ac8706d6b5d7ee4083b9fe195f9a25efeecd6`; this protocol and its
`SPEC.json` must be committed and pushed before `ecbench plan` constructs the
new public points or a measured run begins. No curve, seed, cap, or decision
rule may be changed after seeing an outcome. Amendments require a new
protocol and separate runs; preserve every original row.

## Question and boundary

On the registered Koblitz `a=0,n=37` curve
`icv1-f2m37-tm534059-32aad96b` (subgroup order `r=230603167`, cofactor
`596`), does a smaller folded factor base reduce the **cold, counted** cost
of a verified single-target shared-rank solve? The generic signed-Frobenius
floor is `S_floor=sqrt(pi/(2*74))=0.14569480906717377`, where
`S=total charged group-addition equivalents/sqrt(r)`. The reference is the
same-target 32-lane strong signed-Frobenius rho arm with eight distinguished
bits and step-cap factor 2000. It discards its collision table after each
target. The unmodified 42-column IC arm is the baseline; an identical
42-column control diagnoses local A/A wall noise.

The six IC configurations request `K={4,8,12,16,24,42}` columns from
`compact-orbit-scan`, each with raw-x scan cap 1,000,000, the same fixed
target-blind rank seed `202610042032`, at most 100,000 rank trials, and at
most 512 target residual attempts. A successful base is expected to have
`74K` usable physical points before orbit folding; this is a **prediction**,
not an `fb` identifier or measured inventory. The actual point count,
point-set digest, effective column count, rank, and implementation source
hashes must be recorded for every variant. The names `ic-k4` through
`ic-k42` are local arm labels, with `candidate_id: null` until those exact
records have been constructed and checked.

`SPEC.json` fixes eight independent public-point workloads from
`hash_to_subgroup_v1` at target seed `202610042033`. Each row is a **separate
one-target run**; no IC rank table or rho collision work is shared across
points. There is one warm-up and five measured rounds per point, alternating
arm order at seed `202610042034`. The resource limit is 120 seconds per arm.
The IC and rho arms receive the same point and run seed in each round. The
same rank seed across `K` makes the target-blind query streams comparable;
the point generator and scalar are never given to a solver. If a target
collides with a prior archived target orbit, retain the row and disclose it;
do not replace the point.

## Accounting and decision

Use the native `ecbench` record, including failures. Charge base construction,
folded-table construction, all successful and failed rank queries, relation
verification, final matrix algebra, every base-column scalar check, all
target residual attempts, witness checking, combination, and final full-point
scalar replay. Report these stages, actual factor-base size, folded columns,
rank trials/hits, target attempts, peak memory, and each candidate's cold
counted `S` beside same-target rho. Keep the online five-phase IC interval
and the rho online interval in raw records, but **do not promote or compare
wall-time speedups from the local macOS L0 host**. Field arithmetic, hash and
allocation, and modular combination remain unpriced in the pinned-only GAE
calibration, so counted `S` and IC/rho ratios are **lower-bound diagnostics**;
the fully priced end-to-end cost and speedup remain null.

Before choosing a winning `K`, require all eight target logs and all 40
measured runs of that arm to verify, full rank and exact phase-sum invariants,
and an exact independent replay. Preserve timeouts, exhausted searches,
wrong answers, and OOM rows. Among eligible arms, select the one minimizing
the sum of cold counted GAE over the 40 measured runs, with the smaller `K`
breaking an exact tie. Compare it with the 42-column baseline by paired
workload/round ratios and a workload-resampled 95% interval. A reduction of
at least 20% in this **counted lower bound** with an interval wholly below
one is an engineering lead, not an attack-speed result. Otherwise stop this
column sweep and prioritize a different oracle or a complete native-cost
calibration. A failure to verify is a failure row, never a win. Do not infer
an n41/n53/n83/n131 crossover from n37.

Regenerate the plan, run, audit, and table with the native Rust binary:

```sh
cargo build --release --locked --bin ecbench
target/release/ecbench plan --spec research/ecbench_n37_rank_columns_20261004/SPEC.json
target/release/ecbench run --wait --allow-busy --spec research/ecbench_n37_rank_columns_20261004/SPEC.json --out research/ecbench_n37_rank_columns_20261004/sessions/mac_arm64_l0_01
target/release/ecbench verify --dir research/ecbench_n37_rank_columns_20261004/sessions/mac_arm64_l0_01 --replay-all --out research/ecbench_n37_rank_columns_20261004/AUDIT.json
target/release/ecbench table --dir research/ecbench_n37_rank_columns_20261004/sessions/mac_arm64_l0_01
```

`--allow-busy` may lower the isolation grade but cannot change deterministic
operation counts. Never overwrite a session. The final PR must contain the
spec, plan, host, all raw records, failure rows, replay receipt and hash,
candidate manifests/IDs when exact bases exist, comparisons, decision,
source/input hashes, and a scoped canonical scoreboard update. A separate
physical-host isolation campaign is required before any CPU timing claim.
