# Direct baseline-to-final combined component gate

Registered after the separate #1399→#1416 and #1416→raw-key panels,
before this direct comparison. Multiplying their ratios would join
different contended sessions. Freeze the #1399 full-base binary
`f6_pmull_full_candidate` (SHA-256
`8923553ad79d2647383fd8a929b12fd4ddbfe0f880d95b9d641de295fc6e541b`)
and the raw-key final binary `f6_xmap_full_candidate` (SHA-256
`75009df1bfe412f00390b1d7e707ba5d4cb64d196a143f7888e6bfe235a15d95`).
Both use the same registered n83 K0 curve, dimension-12 cofactor-projected
base of 4,054 usable points, 8,219,485 pairs, public T001, and exact
four-summand miss query.

Run baseline, final, final, baseline, with no builds between arms and a
120-second cap per process. Preserve raw stdout/stderr/exit status,
index-build and exact-query intervals, peak RSS, representative counts,
and correctness. The stage gate is descriptive: report the two-run
median ratio and the two adjacent baseline/final query ratios. A direct
twofold query result requires both adjacent ratios to exceed 2.0 and
identical exact outcomes. Otherwise state that the observed panel did
not establish 2×. This is still an unisolated four-summand component
diagnostic; it cannot establish a complete F6, F4/F5, or IC speedup.
