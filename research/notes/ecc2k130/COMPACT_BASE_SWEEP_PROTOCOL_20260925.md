# Compact four-sum base-size and root-occupancy sweep at n=37/41

Status: preregistered protocol. Measurements and decision belong in this PR after the frozen inputs and producer source are committed. This is a public synthetic, bounded experiment, not an ECC2K-130 attack claim.

## Question and boundary

Test whether reducing the deterministic signed-Frobenius orbit base reveals a useful nonsaturated four-sum regime: measurable target success, rank, root occupancy, and failed-query work at a cold cost that could plausibly beat the same-Q signed-Frobenius rho reference after shared-base amortization. The existing n=53 R=220 base gives 512/512 hits but F^4/q about 14,053, so that hit rate says little about sparse coverage. PR #739 gives only the necessary upper bound p <= min(1,F^4/q), not a random-sum forecast. PR #737 is the clean same-Q full-rank control for a different producer: median rho/direct wall ratios 0.670 at n37 and 0.0634 at n41. It is a comparison boundary, not a timing baseline for the compact extractor.

The source boundary is merged main after PR #737, `ab4a25ea46fcbef8b5f0f3d77cfe66196912df06`. This PR may add only query telemetry and a sweep driver/verifier; any algorithm change invalidates a direct comparison and requires a new source hash and protocol. The factor base is the first R distinct, complete signed-Frobenius orbits from the producer's point-defined deterministic x scan. All sizes at a rung form nested bases. The exact q, R, F and eta arguments are frozen below. Eta fractions have denominator 1,000,000 and are chosen to force exactly the specified R by the existing balanced-orbits loop.

| n | q | R columns | F=2nR points | eta numerator/1,000,000 | F^4/q necessary hit bound |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 37 | 230,603,167 | 1 | 74 | 146 | 0.1300354 |
| 37 | 230,603,167 | 2 | 148 | 1,317 | 2.0805665 (vacuous above 1) |
| 37 | 230,603,167 | 3 | 222 | 5,125 | 10.5328677 (vacuous above 1) |
| 41 | 549,756,390,943 | 4 | 328 | 7 | 0.0210535 |
| 41 | 549,756,390,943 | 8 | 656 | 71 | 0.3368566 |
| 41 | 549,756,390,943 | 12 | 984 | 255 | 1.7053366 (vacuous above 1) |

Tuple permutations, collisions, invalid lifts, and the current root index retaining one witness per root can only lower actual success below this loose upper bound. Do not interpret the bound as a Poisson prediction or an observed probability.

## Frozen inputs and staged execution

For each n, the training scalar stream is the first 4,096 distinct values from `1 + SHA256(b"ECC2K-COMPACT-BASE-SWEEP-20260925-v1/" + str(n).encode() + b"/" + str(counter).encode()) mod (q-1)` for counter=0,1,..., preserving first appearance. This stream is identical across the three bases of a rung. The first 64 targets are the nonclaim calibration stage; after reviewing its resource receipt, run the first 512 targets per arm if every cold process stays below 8 GiB RSS and 1,800 seconds. If a 512-target arm has fewer than R+4 hits and the measured per-failure time projects the full 4,096 targets within 1,800 seconds and 8 GiB, extend only that arm to 4,096 in a new cold process. No early stop on favorable hits. A timeout/OOM is retained as a censored failure, and later arms at the same rung continue when feasible. The final panel identifies exact N per arm and gives binomial 95% Wilson intervals; zero hits remains an upper interval, never p=0.

Three point-only holdouts at each rung are the already frozen PR #737 public hash-derived Q values for seeds `202609250037`, `202609250137`, `202609250237`, and the analogous n41 seeds `202609250041`, `202609250141`, `202609250241`. Their exact coordinates and independent validator scalars come from the committed clean archive in PR #737; none appears in the training scalar schedule. Feed coordinates only to the compact point input. The rho reference uses `koblitz_rho_fixture N 0 signed_frobenius 1 packed 13737 hash:SEED`, one cold process per same Q, with seed 13737 and `RAYON_NUM_THREADS=1`. Reuse the existing clean rho archive only for Q/scalar provenance; rerun rho locally to match the sweep host and commit its raw receipt. If training reaches rank R, recover a holdout log only from the independent base-log solution and a holdout relation, then check [d]G=Q and the PR #737 scalar. A failed holdout relation counts as an unsuccessful full-DLP attempt; do not select a successful Q after seeing outcomes.

Exact compact command per size: `koblitz_s5_sat_instance N 0 ETA_NUM 1000000 natural 1 2000 1 internal`, with `KIC_ALGEBRA_ENCODING=orbit_factorized`, `KIC_ORBIT_LAZY_RELATIVE_SUPPORT=1`, `KIC_ORBIT_BRANCH_ORDER=pair_then_pair`, `KIC_ORBIT_REP_ENCODING=one_hot`, `KIC_ORBIT_BATCH_ONLY=1`, and either `KIC_ORBIT_TARGET_SCALARS` (training) or `KIC_ORBIT_TARGET_POINTS_JSONL` (holdouts). Each arm is a fresh process with its own base scan and root index. The driver freezes input files before launching processes, records source/executable/input/output SHA-256, exact command/env, exit and stderr, monotonic wall, user+system CPU, and peak RSS. No scalar labels reach the holdout producer.

## Measurement, replay, and decision

Record R, F, q, point-set hash, x scan count, regular/exceptional states, two-root candidate entries, unique indexed roots and occupancy `unique/(2*regular)`, cold base/scan/index/target-loop timing, CPU and RSS. Retain every target's hit/miss, query time, S3 calls, partner roots, index hits, and group-lift attempts. An unsuccessful affine query should execute exactly `2*n*regular_states` S3 calls; independently assert that invariant. Keep the full 4,096-scalar schedule even when only a prefix is run and preserve both stages without overwriting. The independent Python replay must reconstruct all admitted relations, verify S3 identities and group lifts, reconstruct and solve the modular rank system, check every recovered base log, and verify every recovered holdout scalar. Save raw failures and any aborted process.

Cold rho/compact wall and calibrated operation counts are reported only for equal completed Q workloads; failed compact Q arms have null speed ratio. Shared-base amortization may be displayed as a separate counterfactual with its measured base/setup and query costs, explicit L, and no victory claim unless full rank and all L holdouts complete. A favorable performance claim needs all three preregistered Q values at both rungs, correct full-rank logs, a 95% paired interval excluding parity, and a fully charged operation comparison. Otherwise classify the result as an occupancy/rank boundary diagnostic and choose whether a changed pair-index architecture is justified. Do not promote an n131 feasibility claim from these toy rungs.
