# Compact inherited F6 bases once they prove one

Registered before implementation or timing. The inherited-basis F6
profile in PR #1466 attributes 91.7% of T7 F4 build time to child
specialisation. A basis that has derived constant one now retains its
full column layout and row storage; later child specialisations clone
that layout even though the solver reads only the refutation. This
opt-in candidate represents a refuted basis as the single column
`1` and single row `[1]`, drops its obsolete layout history, and lets
subsequent children clone the compact state. The old representation is
the same-binary baseline. No accepted solution or oracle decision may
change.

Freeze two exact IC1 candidates differing only in `compact_refuted`.
Use the prepared n17 Koblitz curve, 62 actual usable base points, 29
folded columns, archived public T1/T7, imported certified logs, three
summands, degree three, 8,192 node budget, 32 target trials, one
Rayon thread and default algorithm environment. Disable earlier
profiling and experimental optimization flags. Freeze code, binary,
candidate and input hashes before timing. Run two paired fresh-process
repetitions per target in ABBA order, preserving raw failures and
timestamps. Time the complete one-target online interval with its five
exclusive phases and scalar replay.

Before timing, compare the compact and historical refuted basis through
at least two further specialisations and exact decisive output; run the
prepared F6 geometric-closure control. A successful measured pair must
complete, independently replay the same archived scalar, match target,
attempt count and per-attempt outcomes, and sum its five phases exactly
to online wall. F4 operation counts may change; retain them with basis
counts, oracle calls/refutations, and memory data if available.

Retention gate: both T7 paired complete-online ratios baseline/candidate
at least 1.10, no T1 paired ratio below 0.95, and all correctness checks
pass. Otherwise keep this control default-off. The unisolated Mac makes
CPU ratios exploratory; an isolated-host receipt is needed for a
promoted speedup. This pilot does not establish n83 ordinary relation
yield or IC-versus-rho online speedup.
