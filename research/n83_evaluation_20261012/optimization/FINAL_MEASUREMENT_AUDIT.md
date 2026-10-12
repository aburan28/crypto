# Final fairness and evidence audit

Date: 12 October 2026. Scope: examination of completed native statistics, stage manifests, source/executable bindings and independent correctness receipt. No timed solve or heavy test suite was run by this audit.

## Outcome qualification

All requested stages completed: 18 kernel solves, 108 policy solves, 30 held-out rho solves, and one IC process containing five held-out target solves. All were producer/parent verified. The standalone C++ receipt reports 156 rho answers, five IC answers and five four-point IC witnesses independently replayed, zero invalid rows and status `PASS`. There are 157 worker processes; 161 is the scalar-answer count.

The selected configuration is `fold-L64-DP6-E0`. The held-out primary default/selected full-process geometric ratio is **1.7566708916385678**, source 95% interval **[1.540157442393286, 2.033353415469057]**. Descriptive arithmetic means are **3.9152050208 s default** and **2.152615742733333 s selected**. This is a diagnostic improvement on the registered held-out workload and host, not global rho optimality, another subgroup/hardware result, or an admitted ledger promotion.

Prepared IC/selected-rho online geometric ratio is **0.43101140623628703**, interval **[0.13798766377195482, 1.776387799724782]**. Parity is inside the interval, so no online winner is established. IC is faster in four point comparisons, but the first target costs **13034.487993 ms**. All five targets remain present. Arithmetic means are **3008.6820282 ms IC** and **1950.3966396 ms selected rho**, descriptive ratio **1.5426000881631132**, which reverses the geometric ordering. Preserve both estimands with their meanings; do not select an estimator or discard an outlier after observing the result.

IC preparation is **1149.151699878 s** and its complete cohort process wall is **1164.439319887 s**, from the retained base with one shared setup. This is not five cold single-target replicas. The old **218.034569** fresh-process ratio remains historical and unchanged; new rho times are not joined to old IC timings.

## Selection, holdout and failure gates

Three training points and five held-out points are distinct by signed-Frobenius orbit from the exposed Q, retained base and each other. The law is sampling without replacement conditional on those exclusions, not unconditional IID targets. Training and holdout seeds are disjoint. Every training candidate receives the same nine cases; no fastest seed or per-target policy is selected.

Kernel completion is 18/18, with all nine paired outcomes/trajectories/charges matching. Policy completion is 108/108; all 12 policies are eligible. Selection uses the prospectively declared arithmetic mean complete-process criterion and 2% default preference. Online is secondary for the cold-selected policy. The staged design tests two kernels at one policy, then 12 policies on the selected kernel; it does not cover all 24 kernel-policy combinations or their interactions.

Chronological freeze order is consistent: kernel manifest `1791782879409`, completion `1791782958364`, selection `1791782958366`; policy manifest `1791782978301`, completion `1791783318577`, selection `1791783318584`; holdout manifest `1791783367468`, completion `1791783459429`; IC manifest `1791783502195` (epoch milliseconds). Policy and holdout manifests bind their earlier selection hashes, which match the current serialized selections.

All 30 held-out rho draws and all five IC rows are verified, uncensored and present. Native statistics have no input issues, wrong answers, missing blocks or rejected soundness. The complete-panel ratios therefore use actual measured costs rather than failure penalties. The separate capped expected-cost statistic remains secondary. Native source withholds full-panel primary/online aggregates for incomplete or invalid blocks; wrong answers reject soundness. No silent exclusion affected this result.

Rho's cap is 60 seconds per complete process; IC's is 1800 seconds for the shared cohort, with finite index exhaustion and no per-target hard probe cap. Their unequal caps remain disclosed. Uncensored completion permits observed complete-cost comparison in this panel, not equal success-probability or capped-cost claims outside it. Operation-equivalent S remains unknown because probes and rho steps have different costs.

## Native statistical semantics

Each rho target cost is the arithmetic mean of its three seed draws. The native primary is the geometric mean of the five paired target ratios. Bootstrap resamples paired target blocks, preserving each target mean, with 20,000 percentile draws: seed `202610128399` for rho and `202610128400` for IC/rho. Direction labels and inverted intervals are consistent. There are five target clusters, not 15 independent targets. The bootstrap explicitly approximates exchangeability conditional on the exclusion law, fixed selection and one host.

The IC comparison uses one immutable preparation at rank seed 29 and one extraction sample per target. It measures neither setup variation nor replicated IC timing on each point. Rho and IC run sequentially, with a long IC preparation; target resampling cannot remove temporal host drift. IC assumes prevalidated fixtures while rho charges Q membership/subgroup validation online. Retain that conservative asymmetry and producer phase residuals; use complete intervals instead of summing nested timers.

## Source and receipt bindings checked

Mechanical SHA-256 checks match current runner sources to embedded compiled-source hashes: main `40fe6db5d5a085d57fef8db9608367df3a2bf66900e8ad1056f56a14f1dda408`, monitor `582caadc0c58170c159177014b4e153dc2abbe97245d2e3c9c0c09c350cb6158`, statistics `5bbc8a90eb8cc64ad4bb3b56bf6144dfe067459f3e3cdc4951d0120d184e40cd`.

Worker source `fb14072f4fa5d08e18d8b79130b6d58a220adae0f76369cc9be949691c83eaa0` and the strong-rho, WideGf2 and existing IC field source hashes also match the worker's compiled-source receipts. Both backends use actual runtime hardware CLMUL; Fold records tail-fold reduction and bit-spread squaring. Sampled child start receipts bind executable, point and configuration bytes plus resource limits. The compiled git revision is null, so exact source hashes are the operative binding.

Frozen executable hashes match: rho `21827dcc645fd4591b39dd913619f047d267a98d541854bd7c83ee7143d8bbcb`, runner `b1b94b24fea42a7ac02e2ce21dfd913e7ba4a2b66a536f48d65bb3cfe8044a27`, immutable IC `983f28400b7852fc6628ba97d3ea83b2ab9e32d9944658c7eea831ee89a4abf2`. Machine protocol `bdbed17abe2c9073089d52498d481e42fe67df6b6ae1953467b20125c1faa9db` binds human protocol `012cbae4dabe692ceebd98b012e220f0bd77e6199a59fb120e15ae04cf02e3fa`, which also matches current bytes.

Native statistics input hashes match copied holdout manifest/records and IC manifest/observations. Kernel/policy manifest and record hashes match their selection files. Selection hashes match downstream stage references. These checks found no source or receipt mismatch; they do not replace a full reproducible build.

Authoritative statistics SHA-256: `6250c42f5a7f6644db5c7c39bd5c712155b2667d081309b5e99202f64f8918f7`. Independent arithmetic receipt SHA-256: `6640a9f69e8ef66f3e31a3eab70c18925bdc8140a623aac3559c01fdeb57e294`.

## Correctness and admission limits

The C++ receipt reproduces fixtures/exclusions, checks the retained 99,600-point base/labels, replays all scalar answers and five decomposition witnesses. It does not reconstruct the full 600-row rank/factor-log transcript, reproduce timing/operations, or supply independent-host admission. Correct final answers do not fill those gates.

All four reservation receipts exit zero and record zero observed contended samples on the virtualized AMD EPYC 9V74 host. They do not prove physical-host exclusivity, fixed frequency, absent short interference or sealed ecbench admission. The common new process monitor uses a 5 ms exit poll; outer wall includes observer/scheduling latency and differs from producer online time.

Recorded validation passed four worker tests, 20 native-runner tests and 3947 library tests, with 112 library tests ignored. This audit read logs and did not rerun tests. Initial failed attempts remain preserved.

## Closeout recommendations

Lead the report with completed bounded rho tuning and its diagnostic scope. Present IC's geometric interval, arithmetic means and all five targets together, retaining the 13.034-second case. Separate shared IC preparation from cold single-target costs and keep historical ratios in their original epoch/units. Update pending scaffolds, retain missing full-rank transcript/common S/ecbench/host gates, and leave canonical performance graphs and ledger unpromoted.

Sources: [protocol](PROTOCOL.md), [native statistics](records/panel/statistics.json), [kernel selection](records/panel/kernel/selection.json), [policy selection](records/panel/policy/selection.json), [IC observations](records/panel/ic/observations.json), [independent receipt](verification/final_receipt.json), [runner](native_runner/README.md), [worker](../../../examples/n83_rho_optimization_eval.rs). This audit changes only this file and the optimization README. Conductor pre-edit check returned `command not found`; no heavy job was started.
