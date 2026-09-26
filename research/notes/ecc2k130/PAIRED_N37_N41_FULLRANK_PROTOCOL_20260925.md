# Paired cold full-rank control for n=37 and n=41

Status: six clean main-derived same-Q pairs independently replayed and archived; final PR-head CI 36110446397 revalidated source hashes and replay on current-main merge ref 910a20ee1cba0e9608b2d8d4e379b34d95d2166b. This is a negative control, not a crossover or shared-log batch claim.

## Question and fixed arms

Test whether a point-defined, signed-Frobenius factor base with exact pair support can recover a public hash-derived target at n=37 and n=41 at lower **cold, fully charged** cost than packed signed-Frobenius Pollard rho on the same target. The published n=37 1024-fixture sweep is useful historical evidence, but its two arms derive different seeded scalar targets. Here both fixtures receive the identical `hash:SEED` and must emit identical Q coordinates and recovered d.

Pinned source ref: `3bec69ca636cebc71f828a3035c49b8258e13f5f` (the parent of this branch). Producer Git blob IDs: rank `13e289839fa823db418aa67cddf343066784a51a`, rho `b08cdeb7e2e47ce44e00f22c943272a345b379a9`. Run from a clean checkout of this ref or a documented reconstruction with all source overlays hashed.

| Fixed setting | Value |
| --- | --- |
| Curves | `K_0/F_{2^37}`, then `K_0/F_{2^41}` after n37 pilot passes |
| Targets | public hash seeds `202609250037` (pilot), `202609250137`, `202609250237` (holdouts); n41 analogues `202609250041`, `202609250141`, `202609250241` |
| Direct | `koblitz_rank_fixture N 0 1 2 13737 signed_expanded independent pair_pair_16 1 hash:SEED` |
| Rho | `koblitz_rho_fixture N 0 signed_frobenius 1 packed 13737 hash:SEED` |
| Env | `KIC_RANK_SURPLUS=0`, `KIC_INCREMENTAL_RANK_CROSSCHECK=1`, `RAYON_NUM_THREADS=1`; default reference relation validation stays on |
| Repetitions | 1 pilot plus 2 held-out target seeds, each arm a fresh process, alternate arm order on holdouts |
| Stop | Save every failure. Stop n41 if n37 rank/replay fails, process exceeds 120 s or 2 GiB RSS. A stop is not a favorable censor. |

No scalar label is supplied to either solver. The hash-to-curve target is public and identical in both. The direct arm must reach augmented matrix rank `K+1`, solve, and emit a receipt for every admitted relation. Rho must recover and independently verify the same scalar; report `walk_steps/ideal_steps` and `A=2n` health. The paired verifier must recompute group sums for every direct relation, its orbit-column row and modular rank/solution, all representative logs and `[d]G=Q`, independently of the Rust assertions.

## Cost and decision

Primary comparison: each fresh process's monotonic wall time, user+system CPU time, and peak RSS. Wall includes curve construction, base/support setup, target hash-to-curve, relation search and failures, matrix solve, in-process reference validation, and JSON output. The independent Python replay audits correctness after each producer process exits and is not charged to either arm. Keep producer stage timers as diagnostics only. Record source/executable/input/output SHA-256, host/toolchain, exact commands and env, exit status, stderr and complete JSONL stdout (compressed if needed). Report rho/IC ratios only when both arms finish and replay on the **same Q**. A single pilot cannot justify a crossover. A favorable performance claim requires all three seeds at a rung, 95% paired interval excluding parity, no missing completions, and full cost and operation accounting; otherwise record a control observation, not a win.

This fixture does not amortize factor-log linear algebra across targets: its batch mode runs a full-rank solve for each target, even though it reuses the support table. A true shared-log batch requires a first full-rank solve followed by one or more target relations per new Q, with setup charged once and all targets replayed. Do not label a batch of repeated full-rank solves as the shared-log method. The compact-orbit producer is a separate chain; this PR measures the existing point-defined baseline so a later shared-log comparison has a matched control.

## Historical boundary

The committed [2026-09-12 n37 sweep](../../sat_factor_base_review_20260908/autolab/evidence/20260912-koblitz-vs-rho-no-crossover/README.md) reported 14.778 ms/target for its best direct mode versus 12.822 ms/target for rho over 1024 independently generated targets, a 0.868 rho/IC ratio. Its aggregate receipts did not preserve per-relation replay in Git. That figure is a cross-cohort historical control, not this paired result.
## Measured panel and promotion gate

The six clean-source pairs and all raw/replay receipts are committed under [paired_fullrank_clean_evidence_20260925](paired_fullrank_clean_evidence_20260925/clean_archive_manifest.json). On the measured PR merge ref `a045a0a80d3a8af00fedc85c0342f4e7d026aacd`, median rho/IC cold wall ratios were 0.670 at n37 and 0.0634 at n41; all 820 relations were independently replayed. This is a negative point-defined control. The local source-overlay panel remains labeled provisional. Final PR-head CI 36110446397 revalidated every archived source hash and all six replays on main-derived merge ref 910a20ee1cba0e9608b2d8d4e379b34d95d2166b (main parent 76869fd7e4979f68b601c22df8c61c21884bdd0d).
