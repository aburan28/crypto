#set page(margin: 17mm)
#set text(size: 9pt)

= Rootfinding relation collection: lane, B, and coverage panels

Registered prime subgroup `r=562427`; five frozen local ecbench sessions; release binary SHA-256 `b24171038803d4873d63dbdf5b5c95d5958b0549a7c40266f559fb3331405105`. Every session is L0 on a shared macOS host. See `REPORT.md` and each session directory for complete receipts.

*Execution integrity.* 1,088/1,088 planned executions verified; 816/816 measured replays were identical. All targets were public. The 4/8/16 serial and batch lane arms matched in 72/72 paired counted outcomes and relation counters. Their counted GAE ratios were 1. Native field-operation counters were absent from the sealed records, so the batch-inversion mechanism is supported by the focused source test rather than a measured inversion count in these sessions.

*B sweep.* B48 was selected on train at mean 7,957.42 cold counted GAE. On the independent holdout, B56 averaged 8,151.71 and B48 averaged 8,328.29; the paired B56/B48 ratio was 0.979, with a 95% target-cluster interval of 0.919–1.038. The train winner did not yield a resolved holdout advantage.

*Coverage.* Exported B48 bases at strides 1, 7 and 23 had 2,286, 2,301 and 2,298 unique signed pair-table keys. For an independent uniform nonidentity subgroup point, exact hit probability is `2K/(r-1)`. All observed/uniform target-cluster intervals covered 1 in train and holdout. The holdout stride-23/stride-1 paired cold counted ratio was 0.963 (0.841–1.114). Target-pinned stopping produces unequal probe budgets and the walk probes are correlated.

#image("coverage-efficiency.svg", width: 100%)

*Equal-prefix follow-up.* On eight fresh public targets frozen before the run, each B48 base consumed the same first 8,192 surviving batch-8 walk points. At the predeclared 4,096-probe comparison, stride 23 minus stride 1 yielded 2.750 independent rows per target (paired bootstrap interval 1.000–4.625); rank per lower-bound cold GAE ratio was 1.095 (1.032–1.176). This bounded early-prefix effect was smaller than in the exploratory panel and attenuated by 8,192 probes. See `FIXED_PREFIX_FRESH.pdf` and the raw receipt. The equal-prefix component result does not price the full pipeline.

*Local component diagnostic.* Inversion 79.51 ns; coordinate round trip 4.22 ns; jump choice 1.28 ns; table lookup 3.70 ns; guard insert 3.80 ns per 200,000-sample microloop. Restart-pool draws averaged 5,143.88 ns across 256 draws. These numbers are not isolated whole-collector timings.

*Outstanding gates.* Serialize native field-operation counts; price native costs in a common unit; run the frozen L2 panel with A/A noise and independent-host replay. A matched rho comparison for batch-8 was 7.463 in lower-bound cold GAE (5.817–9.626), with unpriced work on both arms. No larger-field or qualified wall conclusion is drawn here.
