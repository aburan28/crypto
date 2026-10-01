# Discovery source evolution

`probe_01` preserves the first complete implementation, including the serial native
zero reduction and inlined original-system check. It completed all 24 discovery
fixtures and 10,944 observations. The 64-point native arm improved six of nine
groups; no arm passed a dramatic group. This source and result remain sealed.

Before `probe_02` timing, two hot-loop implementation changes were made:

- The native eight-vector zero test uses an explicit balanced unsigned-minimum
  tree on AArch64 and a balanced comparison/OR tree on x86_64. Its truth condition
  is unchanged.
- The complete original-system check is a cold, non-inlined function. It still
  checks every original equation for every projected hit and retains all counters.

The input grids, arm roster, block widths, dispatcher, assignment order, caps,
control solvers, verification predicates and thresholds are unchanged. These edits
were made on discovery data only. They are tested as a combined engineering
candidate; neither edit receives an isolated causal speedup claim.

A full run must bind the selected sealed discovery manifest and identical timed
Rust sources. No holdout timing has informed these edits.

## October 1 continuation

Both original probes completed and remain immutable. Probe 2 regressed: no arm
passed a dramatic group, and the 64-point native arm passed five incremental
groups. Passing the whole mutable result to an outlined helper may inhibit scalar
replacement of loop counters; this is a source-based hypothesis, not an isolated
causal performance finding.

The next source keeps counter updates inline and outlines only the read-only
original-equation predicate. It retains the balanced reduction. Point-budget
subtraction uses the established invariant `points <= cap`: each complete block
is admitted before its points are added. No budget or equation check is removed.

Runtime-gated AArch64 three-input XOR kernels are added for both full 32-bit and
projected 16-bit syndromes. The installed Rust intrinsic implementation and live
feature detector were checked; this local machine reports SHA3 support. Portable
fallbacks preserve the corresponding two-XOR policy. Full-word controls distinguish
the hardware operation from narrowing; hardware availability is recorded per run.

The checkout was advanced to current main before further work. Current AGENTS.md
section 10 requires reserved/pinned CPU isolation and an A/A noise measurement.
New timing uses `isolated_run.py` on Linux and the frozen repository isolation tool.
It refuses unsupported platforms rather than falling back to unisolated timings.
The earlier unisolated probes are historical discovery diagnostics and cannot
discharge this new qualification gate. No new qualified result or holdout result
has been obtained yet.
