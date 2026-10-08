# F6-IC v3: batch the exact point additions, preregistered

V2's new public target needed 44,449 packed point additions and 22,223 exact
residual lookups. It completed in roughly 83 ms versus 119 ms for inherited
F4 on this unisolated Mac, short of the exploratory 2x target. Version 3
keeps the same candidate logic and 256-point bound, but changes its inner
kernel. For each fixed first point, select second base points consistent with
the partial assignment. Use `FastCurve::add_many` to compute all
`first + second` points with one batch inversion. Negate that result vector
and use a second `add_many` call for `target - (first + second)` with one
batch inversion. Preserve the exact packed point lookup and general group
replay of each proposed witness. Count every logical point addition and
residual lookup; additionally count nonempty batch calls. Construction of
scratch vectors and maps is target-dependent PDP cost. Where packed
arithmetic is unavailable, retain v1's general path. A mismatched witness
falls back to general enumeration.

Before the pilot, compare batched, scalar packed and general decisions on
exhaustive n9 positive and negative branches with both partial-bit values;
run the prepared n17 scalar control. Freeze the new source/binary and one
previously unseen public target from fixture seed `20261003022`, algorithm
seed `20261003037`, the same certified n17 state, 62 actual usable points,
29 folded columns, `m=3`, degree 3, 8192 nodes, one worker, one query per
trial, at most 32 queries, 600-second per-process cap and five exclusive
online phases. Commit candidate/workload identities and input hashes before
the first timing. Run exactly six fresh processes in order F4/F6, F6/F4,
F4/F6. Keep every failure and timeout.

All six calls must recover and independently verify the same scalar, keep
all target-dependent attempts and close all five phase costs. The batch and
scalar decisions must agree; batch counts must be positive when the gate
fires. The exploratory engineering target is at least 2x lower **complete
online** F6 time in each paired run, not just fewer F4 reductions. The Mac
cannot promote a controlled CPU wall-time claim. That needs the isolated
benchmark service, and a cross-method rho comparison additionally needs a
sealed `ecbench` session and independent replay. If the full call remains
below the 2x target, preserve the receipt and stop this F6 branch.
