# n37 ecbench whole-solve instruction census

Status: preregistered; no Callgrind outcome has been observed. This is a
follow-on to the independently replayed [n37 factor-base column sweep](../ecbench_n37_rank_columns_20261004/RESULT.md),
whose pinned group-addition-equivalent count left field arithmetic, hashing,
allocation, and modular combination unpriced. Freeze this protocol,
`JOBS.json`, the profiling code, and the Linux workflow in a pushed PR before
the first profile. Preserve failures and amend the protocol in a new commit
without altering earlier output.

## Frozen question and input

Does the n37 K16 factor base still beat K42 when the **entire native
`methods::solve` implementation** is counted in one software instruction
unit? Profile all eight independent public single-target workloads of the
archived session's first measured round, exactly once for each of `ic-k16`,
`ic-k42`, `ic-k42-control`, and `rho-strong`. `JOBS.json` freezes all 32
sequence numbers, points through workload IDs, common per-target algorithm
seeds, and archived recovered logarithms. The existing session spec and its
plan, not knowledge of a planted scalar, rebuild each child input. The
archived session is `ECBS1h9a1d373bbe54`; its `records.jsonl` SHA-256 is
`82c8b31640074b2f829d88e546bf487ee4b598edf4ef26f779984d625ed7de3c`.
All 32 archived rows are verified and the session's 320 measured runs passed
independent Linux replay. Profile a single frozen release binary for every arm
on one Ubuntu 24.04 x86-64 runner with Valgrind 3.22.0 and
`--tool=callgrind --cache-sim=no --branch-sim=no`.

`ecbench profile-input --spec ... --seq N` exports the same explicit-curve
child input as the measured runner. Set `ECBENCH_CALLGRIND_SOLVE=1` on
`ecbench exec`. A Callgrind dump-and-reset immediately before
`methods::solve` excludes workload reconstruction, input parsing, method
resolution, and process startup. The second marker follows the solve even on
an error; count all numbered parts after the first marker through the second,
including any internal phase dumps. `ecbench callgrind-ir --prefix ...` refuses
missing, duplicated, reordered, or incomplete markers. Keep each raw part,
Valgrind stderr, child JSON, input JSON, and parsed Ir receipt as a CI artifact;
commit a compact result and the artifact hash/URL after validation.

The profile must reproduce each archived workload ID, method ID, recovered
logarithm, and successful result. Check the recovered scalar against the
archived independently replayed point record; a failure remains a row and
cannot enter a ratio. Retain toolchain versions, binary SHA-256, Valgrind
version, CPU model/features, environment, target IDs, raw failure outputs, and
all profile-part hashes. The eight `ic-k42`/`ic-k42-control` pairs are an A/A
control: if any pair differs by more than 2% in Ir, keep the results but make
the cross-arm comparison exploratory until the cause is identified. Require
all 32 profiles and all scalar checks before a column decision.

## Decision boundary and claim limit

Use one table with all four arms in the **Callgrind Ir** unit, both raw per
target and divided by `sqrt(r)` for the n37 subgroup. Report paired K16/K42,
K16/rho, and K42/control instruction ratios, including the eight-target range
and a workload-resampled 95% interval. A K16/K42 interval wholly below 1
would retain K16 as an implementation-cost lead; otherwise the counted-GAE
selection is not stable and the next factor-base sweep should stop. A ratio to
rho is a same-binary instruction diagnostic, not an online wall-time speedup.
Do not compare this unit to the generic group-addition floor or use it to
claim an n41/n53/n131 crossover. Callgrind can select a different CPU feature
path than native execution; record the selected path or mark it unknown, and
keep CPU performance claims gated on physical-host isolation and matched
native wall measurements. The frozen Mac L0 wall times remain exploratory.

This census prices every user-space instruction in `methods::solve`, including
the adapter's repeated factor-base inventory and reporting, so it is a
conservative **implementation** comparison. It does not price the external
runner's independent audit or hardware transfer, nor replace the primary
one-target online interval. The next end-to-end gate must reconcile those
boundaries and report online wall time under isolation.
