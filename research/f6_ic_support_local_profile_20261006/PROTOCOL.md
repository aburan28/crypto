# Attribute F6 support-local Macaulay build time

Registered before instrumentation or timing. The exact ranked-column
pilot in PR #1464 matched the full n17 T7 one-target solve but changed
neither its Macaulay build nor complete online wall appreciably. The
support-local builder has three major activities: produce parity-cancelled
rows, collect exact columns, and pack those rows. This diagnostic times
those activities separately on the inherited F6-IC path. They are
nested timers within the existing F4 build timer, **not** extra exclusive
online phases. Instrumentation is opt-in, charged to the online call, and
does not change rows, columns, solver decisions, or the default path.

Freeze one new IC1 candidate for the prepared n17
`icv1-f2m17-tm101-00378d4e` curve, 62 actual usable standard-base
points, 29 folded columns, archived public targets T1 and T7, imported
certified logs, three summands, degree three, 8,192 node budget, 32
target trials, default algorithm environment and one Rayon thread.
Use the original F6-IC algorithm with all earlier packing, streaming,
and bitmap pilots disabled. The only new flag enables profiling. Freeze
source, worker binary, candidate and input hashes before timing. Run
two fresh-process repetitions per target, preserving all raw statuses,
failures, timestamps, five exclusive online phases, scalar replay,
attempts, reductions, matrix dimensions and word operations.

Before timing, test that the timed builder returns byte-for-byte equal
columns and packed rows to the uninstrumented builder, including a
full-support delegation and a zero-row case. Run the existing F6 closure
control. Require each complete target to return the archived scalar and
the five exclusive phases to sum exactly to online wall. Retain the
support-local builder call count, delegated count, row/column/packing
nanoseconds, and row/column totals. Compare their sum against the
existing nested `F4Profile.build_ns` and target PDP; label the residual
unattributed. A missing component is unknown, not zero.

Decision rule: if one component is at least half the F4 build timer on
both T7 repetitions, prioritize it. If the timed components together
cover less than 80% of F4 build in either T7 repetition, inspect the
remaining basis construction before implementing another row/column
micro-optimization. Otherwise prioritize the largest component while
preserving uncertainty. This is a diagnostic on seen public targets,
not a speedup claim. The Mac host is unisolated; n83 ordinary relation
yield and one-target IC-versus-rho remain unmeasured by this protocol.
