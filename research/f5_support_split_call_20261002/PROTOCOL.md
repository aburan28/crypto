# Exact support-separated F5 complete-call candidate

## Source and invariant

Build on the exact support-split diagnostic at
`research/f5_support_split_20261002/RESULT.md`. Add an explicit output
form that uses the full packed n=24 degree-4 matrix. It partitions
degree-4 rows by whether they have a degree-4 term supported inside
variables 0–19, checks the inner and outer projected ranks, then checks
each lower highest-degree row group on columns of its own degree. Return
the unchanged original rows only when **every** group has full row rank.
If any group is dimensionally impossible, rank deficient, empty in an
unexpected way, or otherwise cannot certify the whole matrix, run exact
full-matrix echelon reduction and return its rank and row space. Report
the actual route and all counted XORs. Keep default reduced F5 output
and existing selective echelon and selected-column forms unchanged.

## Correctness and local screen

Before timing, release tests must cover an independent mixed-degree
matrix, a dependent matrix that takes the fallback, zero/partial-word
rows, and unchanged existing GF(2) and matrix-F5 tests. Use the seven
`examples/f4_f2_bench.rs` F5 cases and seed XORs `0`, `badc0de1`,
`5eed2026`, `f5c02a28`. Compare selective echelon (form 2) with the
support-split form in separate processes of one release binary, one
Rayon thread. Per seed run one warmup per arm, five reference/reference
A/A pairs and five alternating reference/candidate pairs. Preserve all
calls, failures, timeouts, source/binary hashes, host, output fields,
route, exact rank, canonical row space, counted work, exclusive phases
and complete-call timings. The primary may return different raw rows
and term count, but must match rank, canonical row space, F5 criterion,
built/pruned counts and columns. Smaller cases must match every output
and route field of selective echelon except the requested form label.

The Apple ARM64 local screen is nonpromoting. Advance to the isolated
Linux gate if all cases and seeds are exact, the support route certifies
the primary on all four seeds, and its counted reduction work is at most
55% of the corresponding selective-echelon reference on all four.
Report the local paired complete-call ratios and A/A spread whether or
not they look favorable; do not select seeds based on timing.

## Final gate

The unchanged requested further-2× claim requires physical isolated
Linux x86-64 one-thread complete-call **median and exact five-pair
bootstrap 95% lower bound both above 2.00×** versus selective echelon
on each of the four seeds. All smaller-case medians must avoid regression
beyond their A/A lower bound, and the separate two-thread controls must
avoid regression beyond their A/A lower bound. Every run must return
verified rank and canonical row space with the selected route reported.
Timeouts, failures, contention and unavailable hardware remain explicit
receipts and cannot be counted as a pass. This is a solver-stage result,
not an IC online or matched-rho speedup.
