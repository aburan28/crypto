# Confirmed production F5 versus inherited F4 on one new public target

**Admission: PASS for the preregistered F5-to-inherited-F4 development gate.**
All nine processes exited zero, reported `complete`, recovered scalar `852`,
and replayed it against the same public target `[6259,101971]`. Each IC row's
five exclusive online phases sum exactly to its reported interval. The fixed
seed was `20261003018`; no scalar was constructed or supplied to any arm.
The one-target workload is `Wf7e18c28f828`. The binary, algorithm inputs and
run order were frozen in the branch before the first measured process.

| Repetition | F5 online ms | Inherited F4 online ms | Paired F5 / F4 | Same-point rho online ms | Rho / inherited F4 |
| --- | ---: | ---: | ---: | ---: | ---: |
| R1 | 150213.866 | 1117.529 | 134.42× | 0.079834 | 0.00007144× |
| R2 | 151753.339 | 668.441 | 227.03× | 0.073416 | 0.00010983× |
| R3 | 184015.748 | 729.949 | 252.09× | 0.082292 | 0.00011274× |
| Median of paired ratios | | | **227.03×** | | **0.00010983×** |

The F5/F4 paired-ratio range is 134.42–252.09×. The median paired
`rho_online_ns / F5_online_ns` is 0.000000484×. Inherited F4 is about
9,105× slower than rho by the median paired inverse ratio on this small
curve. These are one-target online intervals from first target-dependent
query/walk through independent scalar replay; imported certified base logs,
base construction, symbolic template and fixture generation are outside
both online intervals. The rho reference used one signed-Frobenius walk on
the same target, with 39 iterations and 38 walk additions in every run.

The candidate identities are
`IC1N17Ckb1fb62PDP3f5RCsampleLAgaussTDpdpISO0h31be729d8daa` and
`IC1N17Ckb1fb62PDP3f4RCsampleLAgaussTDpdpISO0ha58786952741`.
The factor base had 63 geometric points, 62 usable projected points and 29
sign/Frobenius-folded columns, all with certified precomputed logarithms.
The preparation's original ordinary relation yield and rank receipts remain
in `prepared-f5-runtime-v1`; no ordinary relations were recollected in this
online panel. The complete candidate records, workload record, inputs,
per-run raw stdout/stderr, process statuses and phase-level one-row-per-run
tables are adjacent to this result. Peak RSS was unavailable and is null in
the tables; no memory comparison is claimed.

All 17 target attempts in each IC run used the identical seeded queries:
16 `proved_unsat`, then one witness. Target PDP took at least 99.986% of
every inherited-F4 online interval and 99.997% of every F5 interval. On R1,
F5 recorded 42,549 splits and 110,964 reductions across those attempts;
inherited F4 recorded 1,120 and 2,675. The solver selection changes the
inherited-basis path **and** the automatic split rule, so this result
measures the complete solver-family choice, not the effect of a single
matrix kernel or heuristic in isolation.

The result confirms a large opportunity within this F5-using IC candidate,
but it does not improve the previously measured n53 `PDP4root` path, which
does not call F4 or F5. The n17 rho comparison shows that even the much
faster inherited-F4 target PDP is far from competitive on this curve. No
ECC2K-130 security or complexity claim follows. A larger fresh-target panel
would be needed to estimate natural attempt-count variation and promotion
beyond this development gate.
