# Measurement review of the bounded rho optimization panel

Date: 12 October 2026. Status: design audit and protocol recommendations; no new timed run was performed for this review.

The requested extension strengthens rho before comparing it with IC, so that an optimized IC implementation is not credited for beating an unnecessarily weak reference. The approved design is a bounded kernel and policy search, followed by disjoint validation. It does not certify a globally optimal rho implementation, and it does not reopen or rewrite the completed five-pair study.

The exact legacy model is `icv1-f2m83-t6151469093347-9eba6cd5`, with polynomial low terms `[0,2,4,7]`, curve `y²+xy=x³+x²+1`, subgroup order `8569786107849059`, and cofactor `1128547018`. This is the approximately 52.928-bit subgroup, distinct from the a=0 degree-83 primary confidence gate. Preserve the full curve identity, generator, public target hashes and the metadata erratum in every new record.

## Approved design and its interpretation

| Stage | Frozen cases | Purpose | Selection or reporting |
| --- | --- | --- | --- |
| Kernel | Wide and Fold at 32 lanes, DP4, escape off; three fresh training points × three seeds per point: 18 solves | Compare arithmetic implementations at identical policy and random streams | Select by capped arithmetic mean fresh-process wall; default preference within 2%; online secondary |
| Policy | Selected kernel; lanes 16/32/64 × DP4/6 × escape off/on; the same three training points × the same three seeds: 108 solves | Tune the complete rho configuration, including collision behavior | Same declared selection rule; freeze one choice before holdout |
| Rho holdout | Default Wide/32/DP4/off and the selected configuration; five new points × three new seeds: 30 solves | Measure whether the selected configuration generalizes and compare its complete costs | All outcomes retained; independent-target clustered inference |
| IC prepared holdout | Frozen stock IC binary, retained K600 base, one rank seed 29, one immutable preparation, five holdout targets | Compare each verified target interval with the selected rho policy | Secondary prepared workload; preparation separate; not five cold IC replicas |

Rho has one pinned CPU (CPU 3), 12 GiB address-space limit and a 60-second fresh-process cap. The rho default and the selected policy receive the same frozen caps. IC has the same placement/memory envelope and an 1800-second whole-process cap. Its frozen extraction scans at most one full cycle of `index.states`, but has no hard per-target probe or wall cap. These cap definitions are unequal and remain explicit. The CLI's third argument is a rank seed, not a probe cap.

The two-stage search evaluates only one backend's policy grid after the kernel stage. A kernel-policy interaction could favor a discarded backend at different parameters. Describe the winner as “selected by the registered staged search” rather than “best of all 24 backend-policy combinations.” Searching rho in response to this request while retaining IC's pre-existing K600 policy is an explicit asymmetry in research effort. It is appropriate for strengthening this reference, but supplies no claim that either algorithm has exhausted its optimization space.

## Selection must choose a policy, not a lucky walk

Freeze the three training points, three seed pairs per point, configuration list, execution ordering, score and tie rules before timing. Every policy receives the same cases. The seeds are samples used to evaluate the policy; none is selected as the policy's preferred fast seed. The old `202610061`–`202610065` observations are already exposed and must not become holdout seeds.

For an eligible configuration c, define the training score as the arithmetic mean of its nine whole-process wall costs, with each target weighted equally and its three seeds weighted equally. For a 60-second cap, the scored run is `min(T_process,60s)` under the actual declared timeout boundary. A capped run contributes 60 seconds to this operational score; its uncensored solve time and correctness remain unknown. Retain measured kill/drain duration separately when the outer process lasts longer than the nominal deadline. Do not silently discard the cost beyond the deadline from the raw record.

A wrong answer, invalid relation or failed replay invalidates the configuration; it is not a faster observation with a verification warning. An execution error or OOM remains an outcome and must not benefit a configuration by being assigned its short time-to-error. State prospectively whether those outcomes disqualify the configuration or receive the common full-cap penalty in the training score. The safer default is disqualification for invalid arithmetic/producer errors, and cap penalty plus an explicit completion-rate field for operational budget exhaustion. Preserve all requested but unexecuted cells with their actual reason if an invalid configuration stops early.

Define the default preference precisely: use the incumbent unless the lowest eligible score is at least 2% below the incumbent score. The kernel incumbent is Wide at the default policy; the policy incumbent is the selected backend at 32/DP4/escape off. Record a deterministic ordering for exact ties among nondefault configurations. Report the score, completion count, memory and all nine costs, rather than only the selected row.

Two percent is a conservative selection preference, not a confidence test or proof of improvement. Nine stochastic solves on three points can select a lucky policy. The disjoint holdout is what assesses that selection. Kernel comparisons can be much more tightly paired when both implementations produce identical trajectories, counters and answer on the same seed. Changing lane count, distinguished-point rule or cycle escape can change trajectories/RNG consumption, even with the same seed; common seeds do not guarantee common walk states.

The primary selected policy minimizes fresh-process wall, including reusable rho preparation. Its held-out online time is a secondary measurement of that cold-selected policy. Do not call it the optimized online rho or select another online winner after looking at IC. If online optimality becomes a separate requested question, preregister and freeze a distinct online selection rule and winner before new holdouts.

## Freeze and protect the holdout

Generate five independent uniform nonzero subgroup fixtures under a declared law. Keep validation scalars outside solver input. Independently check each point's membership and subgroup, and reject signed-Frobenius-orbit overlap with the training points and exposed old point where feasible. Freeze three held-out jump/walk seed pairs per point, disjoint from training and old evaluated streams, before selecting the policy; reveal their results only after the selection manifest is sealed.

Run default and selected rho in balanced randomized or alternating order within each target/seed block. No builds, tests, rendering or heavy audits run concurrently with a timed process. Retain the native reservation/monitor receipt, code/compiler/kernel dispatch, exact input and executable hashes, affinity, caps, exit/signal, raw output, wall/CPU/RSS, steps, ideal-step model, distinguished-point entries, fruitless cycles/escapes and replay outcome. A fast per-step field kernel is not sufficient when its walk loses useful collision progress.

The default rho heldout arm is a reference, not a second opportunity to optimize after validation. If the selected policy underperforms it, report that result and retain both rows. Do not choose whichever arm makes a more favorable IC headline on the holdout. A later adoption decision can prefer the default, but a new tuned policy would need fresh independent validation.

The holdout has five independent target blocks, not 15 independent targets. For target j, average the three verified, uncensored rho intervals under each fixed policy. Compare selected/default target means; resample complete target blocks in a native cluster bootstrap. The IC comparison likewise pairs its one target interval with that target's three-seed mean rho interval. Define whether the cross-target statistic is an arithmetic ratio of means or a geometric mean of the five target-level ratios before execution; do not switch after inspection. Preserve all five target-level ratios, quantiles when estimable and the actual uncertainty scope. Five blocks can give a broad, unstable interval and do not justify population or hardware generalization.

## IC preparation and accounting boundaries

The frozen stock IC path loads the retained point base, builds its index, acquires rank, solves logs once, and then loops over the supplied targets. The index/base/log state is immutable across those target queries; sharing it is an explicit prepared-workload assumption. Each target row retains its complete target-dependent interval through scalar replay. Do not divide preparation by five and call that a cold single-target cost. The one whole-process IC row is the cost of one setup and a five-target cohort from a retained base file.

One rank seed 29 makes preparation reproducible, but does not sample the setup-cost distribution. Report that single measured setup separately from the old five independently prepared children. It starts with an already selected/constructed base; original representative selection/construction and OS-cold caching remain outside the boundary. Reusing OS page cache is not the same as reusing factor logs.

The frozen binary's target extraction is deterministic from each public Q. There is one IC extraction sample per held-out target, versus three rho seeds per target. This is declared asymmetric replication, not missing IC runs or 15 paired IC observations. It is defensible for this bounded validation because the logarithm table and query policy are fixed, but timing drift remains possible. Record target order and target start/end producer timestamps where already available; identify the long IC preparation and rho measurement ordering instead of presenting a temporally adjacent pair that did not occur.

IC stock assumes prevalidated fixture inputs while strong rho charges Q membership/subgroup validation online. Retain those rho costs and disclose this conservative asymmetry; no retrospective subtraction from the frozen online estimator. Keep independent external replay cost separately if it lies outside both producer intervals. Producer phase residuals remain visible; the complete online interval is the metric, and nested phase times are not added again.

The stock IC binary emits completed target rows but no target-start or rank-complete event before the target loop. The supervisor therefore cannot enforce an exact new per-target 60-second online deadline without changing the worker. `extract128` scans each index state at most once and returns `None` after exhaustion; it has no hard per-target probe cap. Preserve the 1800-second outer deadline as its only process time cap. On a process kill, any unreported targets remain UNKNOWN; do not assign them zero online cost, a fabricated 60-second time, or UNSAT. The common-unit operation cost S remains null because a probe and a rho step are incomparable.

## Unequal caps and incomplete holdouts

The rho training score uses a common cap within rho and is valid as a bounded operational selection score. It is not an estimator of unconditional mean time to a successful solve. In the IC/rho holdout, rho's per-process timeout and IC's cohort-process timeout plus finite index exhaustion prohibit a common capped-time penalty comparison.

If all 15 selected-rho runs and all five IC targets complete and verify without censoring, compare their observed complete online costs with the declared target-level estimator and interval. Continue to disclose the caps, because absence of failures in a small panel does not establish equal success probability in general.

If any selected rho seed or any IC target fails, caps, is killed or cannot be verified, the full-panel head-to-head winner is unknown. Keep the target block, all consumed wall/CPU/memory, status and available phases. A ratio computed from the completed subset is explicitly conditional on both arms completing and does not establish a full-workload win. Do not replace a missing rho seed by the mean of the other two and then present a three-seed target mean. Do not discard a failed IC point and bootstrap the remaining points as if they were the original holdout population.

Report complete counts and completion-at-budget by arm alongside any conditional cost summaries. A right-censored time is not a mathematical no-decomposition proof. A failed independent check invalidates a solved row even if the producer exited zero.

## Historical observations and updated claim wording

The existing five-pair geometric means remain frozen: IC/strong online `1.382330`, 95% source interval `[1.035774,1.882300]`; fresh-process retained-input wall `218.034569`, interval `[172.358384,285.573727]`. They refer to one exposed point, five seed rounds, the then-current Wide/32/DP4/escape-off strong policy, and that study's outer timing boundary. They are not newly tuned results or five independent target samples.

Do not combine those old IC timings with new optimized-rho timings, or replace the old denominator. New hardware conditions, execution epochs, seed populations, target populations and estimator definitions remain separate. The old fresh-process ratio includes one complete IC setup per child; the new five-target IC cohort shares one setup. Their totals cannot be compared as equivalent cold workloads.

Suitable new wording after successful execution is: “Rho was tuned over the registered two-stage kernel/policy design on three training points. The selected configuration was frozen before five independent held-out points. On those points, the report compares its verified prepared-online intervals with the fixed K600 IC implementation and reports preparation separately.” Name the actual selected policy, host, all counts, interval, default comparison and caps. Use “fastest tested under the registered selection criterion” only for the training selection, and state what holdout established separately.

Remaining admission gates include exact candidate/workload manifests, appropriate host tier, independent-host replay, full factor-log/rank evidence and calibrated common costs where required. The bounded panel can answer the user's baseline-fairness question without pretending that these gates or global optimization are complete.

## Budget and completion tracking

The declared rho count is 156 solves: 18 kernel, 108 policy and 30 heldout. At the prior approximately five-second complete-process average, that is roughly 13 minutes before supervisor overhead, validation and builds; it is a planning estimate, not a measurement. With 60-second per-run caps, the maximum scheduled rho timeout allowance is 9,360 seconds (156 minutes), plus kill/drain and launch overhead. A separate total campaign ceiling, if imposed, must be prospective; retain every not-attempted cell and do not declare a complete-grid winner from unequal partial exposure.

IC's 1800-second process cap is separate. It does not create five 1800-second cold replicas or five individual online caps. Complete the exact authorized count or report each unfulfilled obligation; do not silently shrink the grid or holdout when observed times exceed the nominal estimate.

Sources checked: the original [PROTOCOL.md](../PROTOCOL.md), [README.md](../README.md), frozen `runner/protocol.json`, `records/panel/statistics.json`, the source kernel audit, `examples/koblitz_orbit_dlp_fast.rs`, `examples/n83_strong_point_eval.rs`, and `src/cryptanalysis/koblitz_strong_rho.rs`. This review adds no new mathematical or timing finding and changes no canonical performance graph. The root study report will incorporate its protocol limits and newly recorded results.

Conductor pre-edit check was attempted for this path and returned `command not found`. This file is the scoped audit handoff requested by the parent task; no measured worker, frozen receipt, other report or repository graph was edited by this review.
