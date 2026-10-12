# Fair rho optimization on the exact legacy N83 workload

**Completed diagnostic: 156 rho solves and five prepared IC target solves; all 161 answers and five IC witnesses independently replayed.** This follow-up strengthens rho under the user's requested bounded search. It does not establish a globally optimal implementation or promote a ledger result. The exact model is `icv1-f2m83-t6151469093347-9eba6cd5`, the legacy a=1 approximately 52.928-bit prime subgroup.

Training selected **`fold-L64-DP6-E0`**: existing IC field primitives through a thin rho adapter, 64 lockstep lanes on one CPU, distinguished-point bits 6, and fruitless escape disabled. The [configured worker](../../../examples/n83_rho_optimization_eval.rs) retains signed-Frobenius strong rho, normal-basis canonicalization and batched inversion. Both backends use runtime hardware carry-less multiplication; selection includes preparation and complete process costs.

On five disjoint held-out targets, with three fresh rho seeds per target and arm, the preregistered **default/selected fresh-process geometric mean ratio is 1.756671**, with 95% target-block percentile-bootstrap interval **[1.540157, 2.033353]**. Descriptive mean process costs are **3.915205 s default** and **2.152616 s selected**. These arithmetic means have a different ratio from the paired geometric estimator. The result is conditional on this frozen search, target law and virtualized host.

The secondary prepared IC/selected-rho online geometric ratio is **0.431011**, interval **[0.137988, 1.776388]**; the interval includes parity. IC is faster on four of five targets, while the first target takes **13.034488 s**. With all five targets retained, descriptive arithmetic mean IC time is **3008.682028 ms**, versus **1950.396640 ms** for selected rho, giving the opposite ordering: arithmetic cost ratio **1.542600**. These estimands answer different questions. The experiment does not establish a prepared-online winner, and no estimator or target was selected after observing that difference.

IC measured **1149.151700 s of reusable preparation** and **1164.439320 s for the complete five-target process** starting from the retained K600 base. It used the immutable earlier IC binary, one rank seed 29 and one shared index/log setup. This is one prepared five-target cohort, not five cold single-target replicas. The original five-pair **218.034569** cold ratio remains unchanged historical evidence against the then-frozen default strong rho; old IC timings are never divided by newly optimized rho timings.

| Requirement | Completed evidence |
| --- | --- |
| Freeze target law and selection rules | [PROTOCOL.md](PROTOCOL.md); [machine protocol](records/protocol.json) |
| Compare arithmetic backends | 18 kernel solves; all nine paired verified outcomes, walk counters and reported charges match; no full step trace exported; [selection](records/panel/kernel/selection.json) |
| Evaluate declared policy grid | 108 solves across 12 policies, all completed; [selection](records/panel/policy/selection.json) |
| Freeze choice before held-out timing | Stage manifests bind selection hashes and earlier selection timestamps |
| Default/selected held-out rho | 30 completed verified solves; [native statistics](records/panel/statistics.json) |
| Same-target prepared IC | Five completed verified target rows; [IC observations](records/panel/ic/observations.json) |
| Independent arithmetic replay | 161 scalar replays and five witnesses; [final receipt](verification/final_receipt.json), zero invalid rows |
| Isolation observations | Four stage receipts, all exit 0 and zero observed contended samples; [records](records/) |
| Verification tests | Before timing: four worker,20 native-runner and3947 library tests passed. Integrated delivery:3992 library and four worker tests passed,0 failed,112 library tests ignored; test-only upstream corrections and preserved failures are in [logs](checks/) |
| Statistical/accounting review | [FINAL_MEASUREMENT_AUDIT.md](FINAL_MEASUREMENT_AUDIT.md) |

Native statistics SHA-256: `6250c42f5a7f6644db5c7c39bd5c712155b2667d081309b5e99202f64f8918f7`. Independent final-receipt SHA-256: `6640a9f69e8ef66f3e31a3eab70c18925bdc8140a623aac3559c01fdeb57e294`. [Reports](reports/) contain the source-linked rendering work; layout and packaging completion are tracked separately from this measured outcome.

No timed arm failed, reached its cap, or was excluded. Rho's 60-second process cap and IC's 1800-second cohort cap remain unequal; IC scans its finite index once per target and has no hard per-target probe cap. The five-target bootstrap is conditional on the orbit-exclusion target law and approximates exchangeability; three seeds do not create 15 independent targets. One IC setup does not measure preparation variation.

The full 600-row rank/factor-log transcript, calibrated common operation unit, applicable sealed legacy-wide ecbench adapter and independent-host/isolation promotion gates remain unresolved. Zero observed contention does not prove physical-host exclusivity. Canonical ledger/graphs remain unchanged because this is diagnostic evidence. See [MEASUREMENT_REVIEW.md](MEASUREMENT_REVIEW.md) and [OPTIMIZATION_SCOPE_AUDIT.md](OPTIMIZATION_SCOPE_AUDIT.md).

Original delivery PR 1719 had merged before this follow-up. This work uses branch `codex/n83-rho-optimization-20261012`; the parent study owns commit/PR delivery, integrity manifests and artifact packaging. Frozen original evidence is preserved.

Separate post-merge delivery checks pass: release library 3992 passed / 0 failed / 112 ignored in 402.32 s, and the configured example 4 passed / 0 failed in 12.19 s. These checks do not replace the pre-measurement results or measured-source identities.
