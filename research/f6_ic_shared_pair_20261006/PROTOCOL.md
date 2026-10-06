# F6-IC shared exact pair index: bounded follow-up

Registered before implementation and timing. Depends on the exact pair
lookup and negative pilot in PR #1457. On the frozen n17 T7 target, its
attempt-local index was constructed ten times and counted 62,295 geometric
point additions versus 79,635 in original F6-IC. Nine repeated builds cost
`9 × 1,953 = 17,577` counted additions. The hypothesis is that keeping a
curve-and-factor-base keyed exact pair index across PDP attempts within a
single target solve will cut T7 additions to at most 44,718, without
changing branch decisions, attempts, scalar, or five-phase accounting.

Keep the n17 `icv1-f2m17-tm101-00378d4e` curve, 62 actual usable base
points, 29 folded columns, dimension-6 standard subspace, three summands,
degree three, 8,192 Boolean nodes, 32 maximum target trials, seed
`20261004039`, one Rayon thread, and certified prepared logs from
`research/f6_ic_geometric_closure_20261003/eight_target/`. Freeze the same
public T1 and T7 inputs and workloads as PR #1457. This is a mechanistic
follow-up on previously seen targets, not independent validation of a
general speedup. Compare inherited F4, original F6-IC, attempt-local
`f6_ic_pair`, and new `f6_ic_shared_pair` from one rebuilt worker binary.
All arms get fresh candidate IDs and measured run IDs.

The shared index contains all unordered usable-point pairs, including equal
points, keyed by exact packed full-group sum. Cache it only on the immutable
factor-base object's derived state, keyed by its curve and ordered points.
Cloning or changing the base must invalidate the derived cache. The first
one-fixed branch to cross the **unchanged** `B × (B + 1)` prior-addition
threshold builds the index. Build work is charged to that first target PDP
attempt; later attempts reuse it with zero new build count. The pair lookup
in later attempts begins at their first eligible one-fixed branch; the
threshold applies only to the first build. This clarifies the registered
reuse policy before any timed run. The pair lookup
must check both partial-code orientations and replay every proposed witness
in general group arithmetic. No point outside the exact usable base may
be returned. Keep the original attempt-local and original F6 arms available
as controls. The shared cache is process-local and has no cross-target
reuse in the primary one-target measurement.

Before timing, extend native tests to verify identical branch decisions and
witness replay across the original, attempt-local, and shared paths for
small bases, including duplicate pair sums and source-code assignments.
Verify cache reuse and invalidation on cloned or edited bases. Build an
offline release worker and freeze binary/source/input/candidate hashes.
Run T1 and T7 with two fresh-process repetitions per arm, rotating arm
order; retain timeouts, failures, stdout, stderr, exit status, exact phase
sum, peak RSS via the native host time tool when available, and all hashes.
No Sage job is involved. Derive one IC measurement row per run and report
same-point paired online wall ratios and additions. The local Mac host is
unisolated, so any wall ratios remain exploratory.

Pilot success requires all verified scalars and exact phase sums, no more
than one pair-index build per T7 target solve, at least 40% fewer counted
geometric additions on T7 than original F6-IC, no more than 5% increase
on T1, and **both** T7 repetitions at F4/shared ≥2× and F6/shared ≥1.1×
exploratory online wall ratios. If any gate fails, retain the result and
do not expand to eight targets. If it passes, preregister all eight prior
public targets and seek auditable isolated-host confirmation before making
a controlled speedup claim. End-to-end IC-versus-rho on a previously unseen
single public target is a separate gate; this experiment cannot establish
an IC attack win or n83 transfer.
