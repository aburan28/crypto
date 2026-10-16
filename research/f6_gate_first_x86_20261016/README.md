# Gate-first F6-IC made the default: x86 confirmation (2026-10-16)

`research/lit2026_w5_lane_batched_f4_20261011` measured the gate-first order
(`KIC_F6_GATE_FIRST`) on aarch64 (lazy rewrite) on four cells and left it off,
a default-on candidate "after a check on more cells and on x86 with the eager
rewrite". This study is that check, and the change it supports: gate-first is
now on by default in `koblitz_groebner::gate_first`, and
`KIC_F6_GATE_FIRST=0` restores the specialise-then-ask order as a same-binary
control.

Run report: [`reports/run-20261016T2232Z/report.pdf`](reports/run-20261016T2232Z/report.pdf)
(source `report.typ`, records in `raw/`, `table.tsv` from `summarize.py`,
graphs from `plot.py`; arms from `scripts/run_arms.sh`).

**Result** (one binary from `f25d5e520`, 9 cells, 5 rounds × 3 repetitions,
4 vCPU Xeon with BMI2, L0):

- every one of the 80 targets decided identically in every arm and round;
  stage units identical across rounds;
- F6-IC with the pair index: 1.336–1.733× fewer word XORs on all nine cells,
  median process CPU 1.56–2.12× lower, worst round pairing ≥ 1.25;
- inherited F4: unchanged stage units (no node oracle), CPU within noise;
- the skipped-node counts on the four W5 cells (447 / 432 / 840 / 405) equal
  the aarch64 study's.

**Not changed here.** The RREF root (`KIC_F4_INHERIT_ROOT=rref`), the
study's other default-on candidate, is mixed on this host: 1.10–1.23× fewer
F6-IC stage units on the four larger cells, 0.89–0.98× on four smaller ones.
It stays off.

**Effect on ecbench.** `pdp3-koblitz` with `engine=f6-ic` runs the defaults
(it refuses `KIC_*` overrides) and charges this word count as its unpriced
solver unit. Its records from this commit on show the lower count; verdicts,
relations and group additions are unchanged. The committed sessions that use
the arm (`research/f4_f5_f6_cross_family_20261011/sessions/`) are not among
the sessions CI re-audits; as with the eager-rewrite default (#1671), they
stay evidence at their own commit and are not edited.

**Limits.** L0, one host, CPU time only; m = 3, `max_degree` 3,
standard-subspace bases. No IC-versus-rho statement.
