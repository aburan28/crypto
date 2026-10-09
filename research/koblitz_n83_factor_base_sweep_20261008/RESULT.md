# Degree-83 factor-base sweep: observed evidence

Started 2026-10-08. Public known-answer research. The requested minimum **complete cold index-calculus runtime remains unresolved**. The finite design, native construction panel, arithmetic replay, bounded two-summand diagnostics and S3 storage are complete. Three capped complete-pipeline checks of the separate 53-bit diagnostic arm returned `UNKNOWN_budget`. The primary 81-bit arm has no complete cold measurement.

## Requirement status

| Requested work | Evidence and limit |
| --- | --- |
| Review repository prior art | Complete source/literature inventory in PRIOR_ART.md, with inspected main revision pinned |
| Every combination | All 45,360,000 tuples in the declared finite grid receive one disposition and a reproducible ordinal; choices outside that grid remain unsearched |
| Splitting, symmetry, WDSat, Gray, FES, double large primes, Frobenius | Included in the design; compatibility/capacity audit is complete. Signed-Frobenius construction and Gray-prefix policies are executed. Most solver combinations require adapters and have no full runtime evidence |
| Empirical factor-base panel | 54 constructions complete, containing 42 distinct point sets; exact scan counts, all points and labels retained |
| Arithmetic replay | PASS: 54 bases, 2,748,960 point records and 16,560 representatives checked with generic multi-limb arithmetic |
| Relation-stage experiments | 1,728 fixed public two-summand probes, 87,966,720 exact complement lookups, zero relations; other arities and complete runtimes remain separate |
| Exact support screen | 30 base-size/arity cases derived from the retained point counts and subgroup orders; uniform-target Markov ceilings only, with no solver or fixed-fixture yield inference |
| S3 storage | PASS: 54 compressed objects uploaded, downloaded and byte-hash matched; content-addressed panel receipts uploaded |
| Best total runtime | Unresolved: full relation collection, rank, linear algebra and individual-log phases have not completed the comparison gate |

## Exact revisions and instances

Repository inspection used `c70c32d486a3ac7531fe27f7193d9f09caa58344`. That main snapshot has existing duplicate-definition and scope errors; the initial exporter test build failed before compiling the example. The captured repeated main library build in `verification/main-lib-build.log` was interrupted and must not be described as a completed test run.

The isolated library baseline is `8ab924b935923df9faac25915ed7d9849974de0b`, a common ancestor of the inspected main and the recorded shared branch. The library is unchanged. Export producer revision: `52bb2c4b0477a506ef433fab64f9afa45faf2854`. Bounded relation-stage implementation revision: `41fba579a` (full hash retained in its eventual native receipt). Every panel binds the exact exporter source BLAKE3 as well as the revision. The feature branch is `codex/koblitz-n83-factor-base-sweep-20261008`.

The primary `a=0` arm has subgroup order `2417851639230796216685689` and cofactor 4. The diagnostic `a=1` arm has subgroup order `8569786107849059` and cofactor `1128547018`. Both use modulus `z^83+z^45+z^2+z+1`; their generators and registered model identifiers are explicit in every object header. The 81-bit and 53-bit subgroups stay separate.

## Construction observations

The producer completed all 54 planned constructions in **166.577547042 seconds** of descriptive process elapsed time. It retained **2,748,960 point records**, **16,560 orbit-column records**, and **83,328,324 compressed bytes**. These totals include repeated points across nested or duplicate bases. They are artifact counts, not a count of distinct points in the union.

| Curve arm | Policy | Export variants | Distinct point sets |
| --- | --- | ---: | ---: |
| Primary a=0 | Sequential public x | 9 | 9 |
| Primary a=0 | Hash-derived public x | 9 | 9 |
| Primary a=0 | Gray-prefix public x | 9 | 9 |
| Diagnostic a=1 | Sequential public x | 9 | 3 |
| Diagnostic a=1 | Hash-derived public x | 9 | 9 |
| Diagnostic a=1 | Gray-prefix public x | 9 | 3 |
| Total | | 54 | 42 |

For the diagnostic sequential and Gray-prefix policies, the three seed labels at each K produced exactly the same point set. The sorted coordinate hash detects this independently of header metadata, row order and representative choice. Twelve exported variants are redundant as sets. All variants remain retained because ordering and coefficient representations can matter to a future solver. These labels cannot be treated as independent base draws. K=64, 256 and 600 contain 10,624, 42,496 and 99,600 signed-orbit points per complete base.

Full comparison data is in `pilot-01/manifest.json` and `pilot-01/duplicate-point-sets.json`. No relation yield, matrix rank or complete IC runtime is inferred from the construction observations. Construction timings exclude some preparation and publication work and have L0 status. Other workloads and verification tests overlapped on the shared Mac.

## Replay, relation diagnostic and S3 evidence

`pilot-01/replay.json` is PASS for all 54 bases. It checks every point, candidate index, projection, subgroup relation, Frobenius phase, signed coefficient, closure and distinctness with generic arithmetic against the producer's wide-word kernel. Replay took 1,189.33 seconds of process wall on this host. The compressed-object hashes and panel manifest hash are bound in the receipt.

`pilot-01/probes.json` records 54 completed base rows and 32 fixed public targets per curve arm. The exact two-summand complement oracle made **87,966,720 lookups across 1,728 probes** and found zero relations. Every failed scan covered the corresponding complete base. The matched targets and repeated bases make these dependent observations; no independent success-rate interval or higher-arity relation-yield conclusion follows. Probe execution took 225.24 seconds of process wall. Validation scalars are in a separate public sidecar and were not inputs to that point-only oracle.

`pilot-01/upload-receipt.json` is PASS with 54 object receipts whose downloaded BLAKE3 hashes matched. The uploader exited successfully after also sending the panel manifest, replay receipt, upload receipt, duplicate-set audit and relation-probe sidecars. Its process wall was 141.11 seconds. The exact destination is `s3://crypto-autoresearcher/factor-bases/icv1/etc/koblitz_n83_factor_base_sweep_20261008/`; panel metadata is under `panels/3e46b5b87250e67f18284e03dcbd587d3b5ca502afd5b240bb95037a22bc02b5/`. The construction/replay/probe/upload stages plus the three cold caps consumed approximately **3,530.848 seconds** of active pilot execution, below the authorized 3,600-second limit. Construction uses its internal elapsed timer and the other stages use `/usr/bin/time` process wall, so the sum is a budget audit rather than an admitted performance comparison. Idle time between resumed work sessions is excluded.

## Capped cold diagnostics

The integrated `cold` subcommand accepts the retained, replayed S3-backed base and public fixture zero, rechecks every factor-base point and label before work, and runs the archived compact S3 four-summand index/rank/target pipeline for the 53-bit `a=1` subgroup. The original ordered-pair K=256 run reached its 600-second cap after 606.68 seconds including adapter overhead. The ordered-pair K=64 run reached its cap after 601.32 seconds. The unordered-pair K=64 run reached its cap after 600.59 seconds. All three `cap.json` files say `UNKNOWN_budget` and all three `cold-run.jsonl` files are empty; none supplies a completed cold total, rank receipt, column-log verification or individual-log result. Exact binary SHA-256 digests, source commits, run directories and stage clocks are in `pilot-01/cold-diagnostics.json`, generated by `summarize_cold.py` from retained logs.

An optional unordered-pair index reduces duplicate summand-pair states under Frobenius canonicalization. The test checks that every canonical root key in a small ordered index exists in the unordered index. This is an exploratory implementation choice outside the frozen v1 grid. Its K=64 index completed in **0.771463875 seconds**, with 172,640 regular states and 339,776 root-table entries. Its remaining 600-second solver window did not produce a completed rank/target record. The ordered runs predate phase checkpoints, so there is no matched ordered index timing in this panel; no index or total-runtime speedup is claimed.

## Post-pilot exact support screen

`SUPPORT_MOMENTS.md` proves an exact first-moment upper bound for full smooth unordered m-summand relations with repetition, using each verified base's distinct point count and exact subgroup order. `pilot-01/support-moments.json` retains 30 exact integer/fraction cases for both curve arms, K=64/256/600 and m=2..6. This is a uniform-target mathematical ceiling, not a rate estimate for the fixed public fixtures or a solver benchmark. At primary a=0, K=600, the m=4 ceiling is **1.695987e-6** and the m=5 ceiling is **0.03378543**; at m=6 the bound is vacuous. Thus the current small primary bases offer little full-smooth coverage at m=4/5 under the uniform-target model, while m=6 still needs a verified feasible solver. Splitting, Gray/FES and symmetry do not enlarge the fixed full-smooth sumset; double-large-prime partials require their own graph/rank model. None of these bounds identifies a total-runtime winner.

`RANK_QUERY_SCREEN.md` and `pilot-01/rank-query-screen.json` apply that exact multiset count to a separate one-row rank-query model. For uniform whole-subgroup query targets and at most one full-smooth row per query, no independence assumption is needed to show `Pr(rank >= K) <= min(1, q min(M,r)/(rK))` after `q >= K` queries. At primary a=0 K=600, the upper bound cannot reach one half until at least **176,888,106** four-summand queries or **8,880** five-summand queries. The six-summand bound is vacuous and reduces to `q >= K`. These are necessary query counts only: they do not measure solver time, imply row independence or cover biased queries, multirow solvers and large-prime partials. The current cold driver has not demonstrated uniform full-width sampling on the primary arm.

The generic Koblitz IC source now draws uniform full-width additive scalars and full-width nonzero coefficients when the subgroup order exceeds 64 bits, while preserving the previous at-most-64-bit random stream. The factor-base-log precompute skips its single-word rank tracker for wide moduli and uses the existing BigUint solve. This removes low-limb sampling and rank-gate errors for the 81-bit primary order. It does not turn the retained S3 bases into a completed primary cold run, establish solver throughput, or change any factor-base ranking. No timed pilot work was added.

## Validation and preserved failures

The release profile uses optimization level 3, 256 codegen units and four Cargo build jobs. Toolchain and host facts are in `host.json`; the resolved dependency lockfile is retained in `verification/Cargo.lock`.

| Check | Result | Evidence |
| --- | --- | --- |
| `cargo test --release --lib` on isolated baseline | 2,235 passed, 94 ignored, zero failed | `verification/lib-tests-final.log` |
| Touched export/probe example release tests | 7 passed, zero failed | `verification/export-probe-tests-final.log` |
| Integrated export/probe/cold example release tests | 12 passed, zero failed | `verification/export-cold-tests-final.log` |
| Relevant boundary autolab Python suite | 16 passed | `verification/boundary-tests.log` |
| Exact support screen small-group controls | 2 passed | `verification/support-moments-tests.log` |
| Post-screen release library recheck | 2,235 passed, 94 ignored | `verification/lib-tests-rank-screen-escalated.log` |
| Post-screen touched example recheck | 12 passed | `verification/example-tests-rank-screen.log` |
| Rank-query small-group, threshold and receipt checks | 3 passed | `verification/rank-query-tests.log`, `verification/rank-query-replay.log` |
| Post-screen support and boundary Python rechecks | 2 and 16 passed | `verification/support-moments-tests-rank-screen.log`, `verification/boundary-tests-rank-screen.log` |
| Wide-order sampler and rank-gate release library recheck | 2,236 passed, 94 ignored | `verification/wide-sampler-lib-tests.log` |
| Post-source touched example recheck | 12 passed | `verification/wide-sampler-example-tests.log` |
| Post-source study and boundary Python suites | 5 and 16 passed | `verification/wide-sampler-python-tests.log`, `verification/wide-sampler-boundary-tests.log` |
| Native finite-grid audit | All 45,360,000 dispositions accounted for | `verification/design-audit.log`, `design.json` |

Earlier build and test outcomes remain retained. Initial dependency resolution failed under restricted networking; `--locked` could not be used before this older baseline resolved a lockfile. Some build attempts were interrupted during baseline recovery. The first complete isolated library run had two loopback-network permission failures. The next run passed those tests and failed the pre-existing randomized SQIsign wrong-message assertion. Its focused replay passed; the subsequent full suite passed. The SQIsign source is unchanged, and the earlier failure remains visible rather than being relabeled as a pass.

The post-rank-screen sandboxed full suite again reached the same two loopback bind denials (2,233 passes, 2 failures; `verification/lib-tests-rank-screen.log`). Rerunning the required full suite with loopback access passed all 2,235 tests (`verification/lib-tests-rank-screen-escalated.log`). This recheck changed no Rust source or factor-base object.

Example tests exercise finite coverage and ordinal bounds, the same complete set under binary/Gray enumeration, field encoding rejection, separate subgroup moduli, generic arithmetic replay, corrupted compressed bytes, semantic coefficient corruption with recomputed byte hashes, and point-only relation discovery, negative coverage and budget censoring.

Generic-versus-wide arithmetic replay uses a separate process on the same host and repository. It checks all point and coefficient records, public candidate indices, subgroup membership, source projection, eigenvalue labels, closure and distinctness. External independent-host replay remains pending. This construction panel supplies no boundary-ledger promotion or canonical performance finding.

## Remaining comparison work

The full primary driver needs wide field/point and complete subgroup-modulus integration. An available BigUint elimination routine does not by itself remove single-word sampling or rank gates in coupled drivers. Large-prime graph adapters, WDSat sealed ANF capacities, FES full-system verification, cofactor-aware subspace lifting and domain-preserving symmetry remain explicit gates in `protocol.json`.

Select a base only after complete, matched one-target cold runs, natural rank gain, verified column logarithms and an individual logarithm, all failed work charged, disjoint holdouts, A/A controls, at least five paired rounds, admitted uncertainty and independent validation. The local L0 observations do not establish that minimum. Literal coverage of every possible point set and parameter value is outside this finite experiment.
