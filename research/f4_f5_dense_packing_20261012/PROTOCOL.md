# Bounded direct indexing for Boolean F5 packing

Requested: enhance the existing F4/F5 implementation, keep improvements in a
PR, and add plot groups, graphs and PDFs. Extend open PR #1723. Production
baseline is `0677b1eb0eb4275f5c3366a48689d5620c0bf0b8`; this is the previous
index/readback/fusion candidate, not its earlier baseline. No compounded
performance claim across unrelated main updates.

Hypothesis fixed before production edits: F5's shared `pack_rows` still hashes
each monomial, whereas F4 already has bounded `u16` direct column lookup. Reuse
that selection rule in packing. Select the representation once per matrix,
with separate specialized packing loops so high-bit fallback rows do not
dispatch on the representation for every term. Dense storage must fit within
the minimum bytes of the existing hash entries; reject sparse/high-bit ranges,
empty column sets, overflow and more than 65,535 columns before allocating it.
Retain the hash path and an explicit `KIC_F4_HASH_COLUMN_INDEX=1` control.
Preserve row and column order, cancellation, caps, pruning, rank and word-XOR
accounting. Column-collection bitmaps were considered from the profile, but
are not part of this experiment.

Freeze `examples/f4_f5_row_build_bench.rs` unchanged. All eleven existing
systems remain: `(n,m,D)=(8,8,3),(8,8,4),(10,10,4),(10,10,5)`, each at
`0x2545f4914f6cdd1d` and `0x9e3779b97f4a7c15`; `(12,12,3),(12,12,4)`
and the six-equation shifted-support degree-four system at
`0xd1b54a32d192ed03`. Shifted polynomial coordinates have 64 variables and
support bits `[52,54,57,60,62,63]`; there are six multiplier variables.
The previously frozen holdouts are now known regression fixtures.

Build unmodified code and save the common driver before edits. Record source,
binary, compiler, profiler-package hashes and host details. Use the same
compiler/flags and equal-length baseline/candidate paths, one Rayon worker,
Callgrind 3.24, client-request intervals and merged worker dumps. Complete F4
and F5 calls include construction, packing, elimination and readback; the
standalone criterion is a separate diagnostic. Input construction, warmup,
digesting, verification and serialization are excluded. Report Ir and native
allocator requests/cumulative bytes/live-at-entry/absolute peak; additional
peak is derived above that interval's entry. These are not RSS measurements.
CPU-pressure telemetry is unavailable: no elapsed-time speedup claim.

Fresh campaign: three primary repetitions per variant, three baseline A/A
repetitions, five interleaved A/B/A′ rounds with rotating process order. Eighteen
processes, 792 audited calls, 180-second budget per process. Do not overwrite
runs; retain every failure, timeout, partial dump and negative result. Extra
candidate hash-only control: all eleven systems, one repetition, 33 calls,
outside the acceptance campaign. Native Rust checks exact input/output,
criterion membership, ranks, reports, accounting and F4/F5 reduced outputs.
Use positive Callgrind `totals`, with independent self-cost reconciliation;
retain wrapped `summary` fields unchanged.

Acceptance fixed for this F5-specific change: F5 aggregate fewer instructions,
F4 aggregate not more instructions, no more-instruction verdict in any of the
22 complete-call cases, and no allocation-request, cumulative-byte, absolute
peak or additional-peak regression in any primary or paired call, including
criterion diagnostics. Median paired log-ratios, A/A MAD noise floor, 2%
minimum threshold, 2,000 bootstrap resamples, seed `0x243f6a8885a308d3`, 95%
percentile interval; preserve the existing native implementation and all old
acceptance fields. Add a separately named F5-with-complete-call-guards field.
Standalone criterion instruction verdicts remain fully reported diagnostics;
its unchanged source does not imply identical instruction costs. This scope
is declared before measurements, not changed after observing phase results.

Correctness tests independently reconstruct packed words from column positions,
cover duplicates, arbitrary column order, constants, holes, 63/64 boundaries,
padding, bit 63, sparse fallback and the 65,535-column sentinel limit. Check
serial and Rayon packing without changing the process-wide environment in
parallel unit tests. Existing F4/F5 row-space and solver tests remain intact.

Verification: focused Boolean Gröbner, matrix F5 and PQ Gröbner tests, native
comparator tests and Clippy, touched Cargo example tests, and the mandatory
`cargo test --locked --release --lib -j2 -- --test-threads=3` with a 1,200-second
budget and ordinary library Rayon defaults. Do not weaken tests; retain actual
exit codes, failures, ignored names and unfinished coverage.

Delivery: source-linked Markdown, reproducible native records and hashes,
instruction-history, per-case and allocator graph groups, inspected PDFs,
editable rendering sources, copies in `benchmarks/`, performance-index entry,
commits pushed to the existing ready PR. Check affected canonical graphs and
document frozen ones left unchanged. No CI monitoring or merge in this session.
Conductor is unavailable; its failed scope reservation is retained in records.
