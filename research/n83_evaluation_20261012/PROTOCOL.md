# Degree-83 a=1/K600 validation against strong rho

Frozen 12 October 2026 before new timed solves. Source baseline: `a8194d616ed8c04489ec46177d369f1b1670360a`.

## Requested question

Evaluate the retained degree-83 IC online lead. This means the legacy a=1 curve with a 53-bit prime subgroup, not the a=0 degree-83/81-bit ECC2K confidence gate. Run the same frozen public target through native IC and the current strong signed-Frobenius rho; retain reusable preparation and cold process costs separately.

## Exact fixture and reference

- Curve `icv1-f2m83-t6151469093347-cdcc5432`: `y²+xy=x³+x²+1`, modulus `z^83+z^7+z^4+z²+1`.
- Subgroup order `8569786107849059`, cofactor `1128547018`, generator `(6614244651398184860292015,8868260410029249391539529)`.
- Public target `(2655898345537166180532960,3881526307991673265292615)` supplied as a point; no scalar passed to either solver.
- Retained base: `experiments/koblitz-single-target-n83-20261006/base_n83_K600.jsonl`, 600 signed-Frobenius columns/99,600 points. Input hashes are recorded before launch.
- Reference: existing `StrongRhoG<WideGf2,u64>` through the new exact-public-point example, default 32 lockstep lanes, DP bits 4, step cap factor 2000, fruitless escape disabled. Normal-basis and G-only jumps are reusable and excluded from online time. Q validation, target-dependent walk starts, walk, collision recovery and separate generic scalar replay are online.
- IC: baseline `koblitz_orbit_dlp_fast`, same retained base/point, rank threads 1 and target threads 1. No new selector, oracle or collector is being claimed. Archive its raw five stages and any interval residual explicitly.

## Execution and stop rules

One frozen point per child; fresh process and empty DP/index/log state. Common enforced envelope: one pinned CPU (CPU 3), 12 GiB virtual-address limit, process-group timeout 1800 seconds per arm. Record actual user/system CPU, monotonic process wall, Linux wait4 RSS in bytes, raw stdout/stderr, semantic status and source/binary/input hashes. Do not accept exit 0 as IC success without one verified target, no failed target, full rank and correct target/scalar.

First independent archive/fixture validation and narrow-vs-wide native tests. Then one excluded warm-up per method, five strong/strong A/A pairs on identical within-pair seeds, and five interleaved strong/IC pairs. Strong seeds 202610061..202610065; IC rank seeds 21..25. Alternate order by round. No builds, tests or independent heavy jobs during a timed process. All failures remain rows; a failed prerequisite or cap is reported rather than changing the protocol after results.

This is five timing repeats of one point, not five independent DLP targets. Analyze within-point paired ratios and noise only; no population claim. Preserve online intervals, reusable setup, complete-process CPU/wall and memory as separate metrics. No calibrated conversion exists between IC root probes and rho steps; common operation-equivalent S remains null.

## Admission and scope

Use the native `isolated_bench` reservation/monitor for the diagnostic executions. It does not by itself prove the full current ecbench L2/physical-host requirements. The current sealed ecbench path lacks this exact legacy wide-curve/IC adapter: it selects the distinct registered wide modulus and its IC plugins stop at one-word fields. Record that gap; do not substitute the a=0 curve or a different generator/field.

New timings are diagnostic unless all applicable isolation, identity, common-cost and independent-host gates are met. A faster verified result is not automatically an admitted `vs_rho` promotion. Correctness must be independently replayed with standalone native arithmetic. The historical packed-reference ratios are not evidence against the required strong baseline.

Conductor CLI invocation was attempted and returned command-not-found; this checkout has no tracked `.conductor/WORKFLOW.md`. Work is isolated on feature branch `codex/n83-strong-evaluation-20261012`; no existing sealed evidence is edited. Repository-delivery PR will include code, protocol, raw receipts, independent check, report and graphs. No large cloud rental or production target is involved.
