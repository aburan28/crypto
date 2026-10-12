# n83 F6 wide search: delayed child assignments

Registered before implementation or measurement on 2026-10-05. The
`wide_groebner` search currently builds one full substituted polynomial
frame for every valid next-point coordinate at each split. The first
node of the pinned n83, dimension-12, three-summand block has about
4,057 such alternatives. The hypothesis is that storing a shared
parent frame and the coordinate choice, then substituting only when a
child is visited, reduces peak memory and root-search wall time without
changing which assignments are searched or their outcomes.

Freeze the K0 curve `icv1-f2m83-tm6151469093347-debefd74`, the standard
dimension-12 source subspace, planted source indices `[0,2,4]`, the
two-word field table, default wide-search options and release Rust
profile from the immediately preceding n83 block admission. The
baseline source is commit `9b969fe26130b9eaa2548395bccfc19a0c5ab57a`;
its frozen native probe binary has SHA-256
`17678d80e0bb221c43b3008e28c73557e79b0b25fec682f3fac6067abb52f16e`.
The previous root receipt is `research/f6_n83_wide_block_20261005/root_probe.jsonl`.
Use the identical existing example and input for the candidate, with a
node budget of one. The only algorithmic change is child scheduling and
materialization. Preserve the old binary; never rebuild it in place.

Correctness gate: all `wide_groebner` tests pass, including exhaustive
small-curve decompositions, field-table checks and the n83 planted
equation/cofactor checks. The candidate probe must return the same
`exhausted` status with no witness, equal node and reduction counts,
and must never call a budget stop a proved miss. Any found witness in a
later diagnostic must replay exactly in the curve group.

Performance gate: run five paired baseline/candidate probes in
alternating order on this physical arm64 Mac with a 120-second timeout
per probe, recording each JSON line, exit status and stderr. Use the
median of query intervals and peak RSS. A successful local optimization
has at least 1.5× lower query time and at least 3× lower peak RSS,
without a correctness change. Preserve regressions and all failures.
These timings are **exploratory** because this host lacks the required
CPU-isolation receipt. The measured slice is a planted block root, not a
complete F6 decomposition, ordinary relation, full IC pipeline, or
same-target rho comparison. Those speedups remain unknown regardless of
this gate.

Stop if the candidate test or probe exceeds 120 seconds, or if peak RSS
exceeds 8 GiB. Record source, binary and input hashes, exact build/run
commands, compiler/OS/CPU, all outcomes, and the result decision in this
PR. Keep target-independent setup separate from query time.
