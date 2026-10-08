# F6-IC adaptive pair-sum closure pilot

Registered before code or measurement. The completed n17 F6-IC panel
showed a 1.291–1.583× exploratory complete online advantage over inherited
F4, short of the user's 2× goal. Its exact one-fixed-source geometric
closure repeatedly computes two point additions for each eligible second
base point. This pilot tests an exact pair-sum index built only after enough
such work has already accumulated within one decomposition attempt.

Use `icv1-f2m17-tm101-00378d4e`, the frozen 62-point usable standard
dimension-6 base with 29 folded columns, three summands, degree three,
8,192 Boolean nodes, at most 32 target trials, one Rayon thread, and the
imported certified log state in
`research/f6_ic_geometric_closure_20261003/eight_target/`. Freeze its
previously unseen public targets T1 and T7, their archived prepared inputs,
and algorithm seed `20261004039`. T1 is the low-work control; T7 had the
largest observed geometric work and 11 target PDP attempts. Do not select
different targets after seeing the pilot result. Compare the inherited F4
worker, unchanged F6-IC, and an opt-in `f6_ic_pair` worker on each same
public point. Preserve the original inputs and create exact new candidate
inputs differing only in the declared solver code. New IC1 candidate and
run IDs must bind the source, flags, and prepared-state hashes before a
speed claim.

The candidate pair index contains **every unordered pair, including equal
indices**, of usable factor-base points, keyed by their exact full-group
sum. Build it with `FastCurve::add_many`; charge all pair additions and
batch work inside target PDP. At a one-fixed-source branch, compute
`target − first_point`, look up its pair sums, test both pair orientations
against the two remaining partial source-code assignments, and replay any
witness in general group arithmetic. A missing key or exhausted bucket is
an exact branch refutation. Keep the general-arithmetic path and prior
batched path as fallbacks. Build the index lazily only when cumulative
one-fixed-closure logical additions in that attempt reach
`B × (B + 1)`, where B is the actual usable base size; do not build it on
a branch with no matching first point. The index is rebuilt per PDP call,
so its full cost is charged to that call. No target or base cache is added.

Before timing, use native tests to compare every branch decision and
independently replay witnesses on small exact bases, including duplicate
pair sums and both pair orientations. Run each frozen target with two
fresh-process repetitions per arm in alternating order, preserving every
failure and the five exclusive online phases. All scalar answers must
replay on the public target. Report per-target paired online wall ratios,
attempts, reductions, pair-index builds, geometric additions, RSS and an
identical F6-IC A/A control. The Mac host is unisolated; wall ratios are
exploratory, and no CPU speedup claim is promoted from this pilot.

Pilot success requires exact verified results, at least a 25% decrease in
geometric additions on T7, no more than a 5% increase on T1, and an
exploratory complete online F4/`f6_ic_pair` ratio of at least 2× on T7.
A failed pilot is a retained result, not authority to change the frozen
threshold or targets. If it passes, preregister and run all eight earlier
public targets, then seek an auditable isolated-host confirmation before
claiming 2× F6. The end-to-end IC-versus-rho objective still requires a
same-point single-target solve and rho reference under accepted isolation;
this pilot does not substitute for it.
