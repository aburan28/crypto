# Inherited F4: eager rewrite without BMI2, through shift runs (2026-10-16)

Before this change the eager rewrite of kept rows (`KIC_F4_EAGER`, #1671)
ran only on BMI2 hosts; every aarch64 host fell back to the lazy rewrite
with the bit scatter. A step-plan segment maps its source bits to its
destination bits in increasing order, so it splits into a few runs moved by
one common shift; `pdep(pext(x, src), dst)` is then a handful of
mask-and-shift operations (`ShiftRun`, `apply_step_runs` in
`src/cryptanalysis/inherited_f4.rs`). Hosts without BMI2 now take the eager
rewrite through the runs. `KIC_F4_EAGER=runs` forces the runs on a BMI2 host;
`KIC_F4_EAGER=0` and `KIC_F4_REWRITE=scatter` are unchanged controls. The test
`rewrite_plan_matches_the_bit_scatter` checks the runs and the BMI2 step
against the bit scatter on random composed maps, on every architecture.

Run report: [`reports/run-20261016T2335Z/report.pdf`](reports/run-20261016T2335Z/report.pdf).

**Result** (x86 proxy; one binary from `1750ec759`; 9 cells, 5 rounds × 3
repetitions; L0):

- 80/80 targets decided identically in every arm; stage units identical
  across rounds; the runs charge exactly the BMI2 path's units;
- against the bit scatter (the old non-BMI2 path): median CPU 1.22–1.70×
  lower on six of nine cells on both engines, 1.27× / 1.07× on K1/2^19, level
  on K1/2^17 d6 and d7;
- against `pdep ∘ pext`: 2–28% slower, so BMI2 hosts keep the BMI2 path.

**Stage units.** Both eager paths charge more word XORs than the lazy ones
(every kept row is rewritten and charged) while using less CPU. A non-BMI2
host now charges the eager count, which is the count x86 hosts have charged
since #1671: the two architectures charge the same units again. Committed
records made on aarch64 under the lazy path stay evidence at their commit.

**Not measured.** No aarch64 host was available; the ratio there is
unmeasured (task: Arm64 measurements via CI runners).
