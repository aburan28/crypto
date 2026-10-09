# Degree-83 factor-base sweep: observed evidence

Started 2026-10-08. Public known-answer research. The requested minimum **complete cold index-calculus runtime remains unresolved**. The finite design, native construction panel, arithmetic replay, bounded two-summand diagnostics and S3 storage are complete. Three capped complete-pipeline checks of the separate 53-bit diagnostic arm returned `UNKNOWN_budget`. The primary 81-bit arm has a replay-bound cold runner with frozen pilot K=64 zero-trial preflight and one exact m=2 enumeration trial, plus capped m=3 and K=256 attempts, but no complete cold measurement.

## Requirement status

| Requested work | Evidence and limit |
| --- | --- |
| Review repository prior art | Complete source/literature inventory in PRIOR_ART.md, with inspected main revision pinned |
| Every combination | All 45,360,000 tuples in the declared finite grid receive one disposition and a reproducible ordinal; choices outside that grid remain unsearched |
| Splitting, symmetry, WDSat, Gray, FES, double large primes, Frobenius | Included in the design; all finite tuple dispositions are audited. Signed-Frobenius construction and Gray-prefix policies are executed. N83 solver capacity and most backend adapters remain unvalidated, with no full runtime evidence |
| Empirical factor-base panel | 54 constructions complete, containing 42 distinct point sets; exact scan counts, all points and labels retained |
| Arithmetic replay | PASS: 54 bases, 2,748,960 point records and 16,560 representatives checked with generic multi-limb arithmetic |
| Relation-stage experiments | 1,728 fixed public two-summand probes, 87,966,720 exact complement lookups, zero relations; other arities and complete runtimes remain separate |
| Exact support screen | 30 base-size/arity cases derived from the retained point counts and subgroup orders; uniform-target Markov ceilings only, with no solver or fixed-fixture yield inference |
| S3 storage | PASS: 54 compressed objects uploaded, downloaded and byte-hash matched; content-addressed panel receipts uploaded |
| Primary base to generic solver adapter | PASS: replay-bound K=64 S3 object reconstructed as generic ordinary and signed Frobenius orbits; solver stage unexecuted |
| Primary cold runner | K=64, m=2 exact-enumeration preflight PASS; one trial returned `UNKNOWN_trial_cap` with zero relations. K=64 m=3 and K=256 m=2 attempts reached `UNKNOWN_budget` caps. All retain null total runtime and winner |
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

A subsequent primary K=64 adapter check added **5.222421667 seconds** of process wall, bringing the active pilot audit to approximately **3,536.070 seconds of 3,600**. Its budget and semantic receipts are `pilot-01/primary-adapter-a0-hash-64-budget.json` and `pilot-01/primary-adapter-a0-hash-64.json`. The check loaded the retained S3-backed compressed object through the verified local panel, reconstructed 10,624 points in 128 ordinary Frobenius and 64 signed orbits, and passed the pinned-curve, label, coefficient, membership, closure, uniqueness and point-set-hash checks. This is an adapter validation result only; `solver_stage_executed=false`, `total_index_calculus_runtime_ms=null` and `selected_best_total_runtime=null` are recorded in the receipt.

## Capped cold diagnostics

The integrated `cold` subcommand accepts the retained, replayed S3-backed base and public fixture zero, rechecks every factor-base point and label before work, and runs the archived compact S3 four-summand index/rank/target pipeline for the 53-bit `a=1` subgroup. The original ordered-pair K=256 run reached its 600-second cap after 606.68 seconds including adapter overhead. The ordered-pair K=64 run reached its cap after 601.32 seconds. The unordered-pair K=64 run reached its cap after 600.59 seconds. All three `cap.json` files say `UNKNOWN_budget` and all three `cold-run.jsonl` files are empty; none supplies a completed cold total, rank receipt, column-log verification or individual-log result. Exact binary SHA-256 digests, source commits, run directories and stage clocks are in `pilot-01/cold-diagnostics.json`, generated by `summarize_cold.py` from retained logs.

An optional unordered-pair index reduces duplicate summand-pair states under Frobenius canonicalization. The test checks that every canonical root key in a small ordered index exists in the unordered index. This is an exploratory implementation choice outside the frozen v1 grid. Its K=64 index completed in **0.771463875 seconds**, with 172,640 regular states and 339,776 root-table entries. Its remaining 600-second solver window did not produce a completed rank/target record. The ordered runs predate phase checkpoints, so there is no matched ordered index timing in this panel; no index or total-runtime speedup is claimed.

## Post-pilot exact support screen

`SUPPORT_MOMENTS.md` proves an exact first-moment upper bound for full smooth unordered m-summand relations with repetition, using each verified base's distinct point count and exact subgroup order. `pilot-01/support-moments.json` retains 30 exact integer/fraction cases for both curve arms, K=64/256/600 and m=2..6. This is a uniform-target mathematical ceiling, not a rate estimate for the fixed public fixtures or a solver benchmark. At primary a=0, K=600, the m=4 ceiling is **1.695987e-6** and the m=5 ceiling is **0.03378543**; at m=6 the bound is vacuous. Thus the current small primary bases offer little full-smooth coverage at m=4/5 under the uniform-target model, while m=6 still needs a verified feasible solver. Splitting, Gray/FES and symmetry do not enlarge the fixed full-smooth sumset; double-large-prime partials require their own graph/rank model. None of these bounds identifies a total-runtime winner.

`RANK_QUERY_SCREEN.md` and `pilot-01/rank-query-screen.json` apply that exact multiset count to a separate one-row rank-query model. For uniform whole-subgroup query targets and at most one full-smooth row per query, no independence assumption is needed to show `Pr(rank >= K) <= min(1, q min(M,r)/(rK))` after `q >= K` queries. At primary a=0 K=600, the upper bound cannot reach one half until at least **176,888,106** four-summand queries or **8,880** five-summand queries. The six-summand bound is vacuous and reduces to `q >= K`. These are necessary query counts only: they do not measure solver time, imply row independence or cover biased queries, multirow solvers and large-prime partials. The current cold driver has not demonstrated uniform full-width sampling on the primary arm.

The generic Koblitz IC source now draws uniform full-width additive scalars and full-width nonzero coefficients when the subgroup order exceeds 64 bits, while preserving the previous at-most-64-bit random stream. The factor-base-log precompute skips its single-word rank tracker for wide moduli and uses the existing BigUint solve. This removes low-limb sampling and rank-gate errors for the 81-bit primary order. It does not turn the retained S3 bases into a completed primary cold run, establish solver throughput, or change any factor-base ranking. That source-only change added no timed pilot work.

The primary adapter now imports the replayed `a=0` K=64 hash-derived base into the generic `FrobeniusFactorBase` representation without using its single-word field builder. The release example tests exercise a small imported base through the generic solver's zero-trial preflight and reject corrupted construction metadata and coefficients. The retained K=64 process check validates the real object and orbit representation, but does not call relation search. Therefore this step changes neither the runtime ranking nor the primary cold-run status.

The bounded primary cold runner next used that same K=64 base and public fixture zero. Two preliminary receipts bind adapter-source BLAKE3 `039222d7f9abe00f25b64db5f6d5199b5d7b48210b8afd97e4df2c5e1f9c51f1`, which precedes the SAT fail-closed guard and is retained only as exploratory evidence. Their zero-trial preflight finished in **3.96936575 seconds** of outer process wall (`PREFLIGHT_ONLY`), and an m=2 exact-enumeration run with one relation trial finished in **5.45737475 seconds** (`UNKNOWN_trial_cap`, zero relations). Their configurations, flushed events, reports and separate outer budget receipts are in `pilot-01/primary-cold-a0-hash-64-m2-enumerate-{preflight,one-trial}-01/` and sibling `-budget.json` files. Both had an inner 20-second wall cap and retained a null total runtime and winner. These 9.4267405 seconds brought the earlier active pilot audit to approximately 3,545.497162167 seconds.

After building the guarded pilot source at commit `7b956daaf` (adapter BLAKE3 `876576c9a3d46333ab1f9d7d7a63ffb249e11e9f2b1ccfce262080b047d29193`), the same K=64 m=2 preflight finished in **2.322702084 seconds** of outer process wall (`PREFLIGHT_ONLY`). One m=2 exact-enumeration trial finished in **3.699593292 seconds** (`UNKNOWN_trial_cap`, zero relations); the report records 1.794292 seconds inside relation collection. An m=3 exact-enumeration trial reached its inner 15-second cap after **15.035790625 seconds** of outer process wall (`UNKNOWN_budget`); flushed events show base ready and relation collection started, with no completed trial or summary. A K=256 m=2 preflight reached its 11-second cap before the base-ready event (**11.061315917 seconds** outer wall). A K=256 m=2 one-trial attempt reached the base-ready event with 42,496 points and 256 signed orbits after 12.268 seconds, then hit its 14-second cap during relation collection (**14.057679125 seconds** outer wall); it has no completed trial or summary. The frozen pilot directories are `pilot-01/primary-cold-a0-hash-{64,256}-m{2,3}-enumerate-*-final/`, each with a sibling budget receipt. Another workspace was compiling on the shared host, so these frozen pilot wall observations are functional and budget evidence only. The five frozen pilot cases charged **46.177081043 seconds** and bring the approximate active-run audit to **3,591.674243210 of 3,600 seconds**, leaving **8.325756790 seconds**. `pilot-01/primary-cold-budget-audit.json` accounts for all seven primary cold receipts. No stage duration here is a matched throughput or complete-runtime estimate.

A source audit found the original generic SAT S4 finite-coordinate-domain trie accepted `u64` codes and extracted only `raw_bits().first()` from F2^83 points. This lost the high 19 coordinate bits. The source now stores exact `BigUint` codes in both SAT domain paths, with a degree-83 trie regression that distinguishes coordinates sharing their low word. The primary cold CLI still rejects `sat-m3` before it creates a run directory; its full N83 S4 encoding, model lifting and resource capacity have not been independently replayed on a retained base. The declared SAT/WDSat grid remains an unexecuted recipe. This restriction does not affect the exact m=2 enumeration receipts above.

A source audit of cofactor-class admission found a distinct wide-field correctness gate: the general-field orbit walk hashed points with `pack_point`, which is documented as injective only through degree 62. The walk now uses exact multi-limb `point_key` identities, including its reference implementation; degree 63 also routes around the packed fast admission path. This matters when future N83 bases include nontrivial cofactor classes; all 54 retained bases were cofactor-projected into their respective subgroups, so the correction does not change their stored point sets or establish a runtime ranking. The regression constructs two valid degree-83 points with the same packed key but different Frobenius orbits.

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
| Post-adapter release library suite | 2,236 passed, 94 ignored, zero failed | `verification/primary-adapter-lib-tests.log` |
| Replay-bound primary adapter example tests | 13 passed | `verification/primary-adapter-example-tests.log` |
| Post-adapter study and boundary Python suites | 5 and 16 passed | `verification/primary-adapter-python-tests.log`, `verification/primary-adapter-boundary-tests.log` |
| K=64 retained-object primary adapter check | PASS; 5.222421667 seconds process wall | `verification/primary-adapter-64.log`, `pilot-01/primary-adapter-a0-hash-64-budget.json` |
| Wide-key degree-83 on-curve collision regression | 1 passed | `verification/wide-key-focused.log` |
| Post-wide-key release library suite | 2,237 passed, 94 ignored, zero failed | `verification/wide-key-lib-tests.log` |
| Post-wide-key touched example suite | 13 passed, zero failed | `verification/wide-key-example-tests.log` |
| Post-wide-key study and boundary Python suites | 5 and 16 passed | `verification/wide-key-python-tests.log`, `verification/wide-key-boundary-tests.log` |
| Primary cold runner release example suite | 14 passed, zero failed | `verification/primary-cold-example-tests.log` |
| Post-runner release library suite | 2,237 passed, 94 ignored, zero failed | `verification/primary-cold-lib-tests.log` |
| Post-runner study and boundary Python suites | 5 and 16 passed | `verification/primary-cold-study-python.log`, `verification/primary-cold-boundary-python.log` |
| Wide SAT coordinate-domain regression | 2 focused tests passed, including degree-83 high-bit membership | `verification/wide-sat-domain-focused.log` |
| Post-wide-SAT release library suite | 2,238 passed, 94 ignored, zero failed | `verification/wide-sat-lib-tests.log` |
| Post-wide-SAT touched example and Python suites | 14, 5 and 16 passed | `verification/wide-sat-example-tests.log`, `verification/wide-sat-study-python.log`, `verification/wide-sat-boundary-python.log` |
| Pinned solver availability and core-size audit | PASS: four main-only backend modules and 910 native S4 core variables checked | `verification/solver-gates-source-check.log`, `SOLVER_GATES.md` |
| Post-audit release library and touched example suites | 2,238 passed with 94 ignored; 14 passed | `verification/solver-gates-lib-tests.log`, `verification/solver-gates-example-tests.log` |
| Post-audit study and boundary Python suites | 5 and 16 passed | `verification/solver-gates-study-python.log`, `verification/solver-gates-boundary-python.log` |
| Current S4 encoder auxiliary bound | PASS; at least 593,364 SAT variables and 2,349,149 AND-definition clauses at `l=83` | `verification/s4-auxiliary-bound-check.log`, `verify_s4_aux_bound.py` |
| Wide batch-add and retained-primary witness-order checks | 1 focused library and 1 focused example test passed | `verification/wide-batch-focused.log`, `verification/wide-enum-focused.log` |
| Post-wide-enumerator release library, touched example and Python suites | 2,239 passed with 94 ignored; 15, 5 and 16 passed | `verification/wide-enum-lib-full.log`, `verification/wide-enum-example-full.log`, `verification/wide-enum-study-python.log`, `verification/wide-enum-boundary-python.log` |
| Factored S4 source bound and small-system equivalence | PASS static bound; 2 focused release tests passed | `verification/s4-factored-bound-check.log`, `verification/s4-factored-focused.log` |
| Factored S4 generic decomposition integration | PASS small explicit-orbit group-lifting comparison and fail-closed unsupported modes; N83 capacity still untested | `verification/s4-factored-integration-focused.log` |
| Post-integration release suites | 2,242 library tests passed (94 ignored); two edited example targets passed (zero embedded tests), N83 example 15 passed; study and boundary Python 5 and 16 passed | `verification/s4-factored-integration-lib-full.log`, `verification/s4-factored-integration-examples-full.log`, `verification/s4-factored-integration-n83-example-full.log`, `verification/s4-factored-integration-study-python.log`, `verification/s4-factored-integration-boundary-python.log` |
| Post-factoring release library, touched example and Python suites | 2,241 passed with 94 ignored; 15, 5 and 16 passed | `verification/s4-factored-lib-full.log`, `verification/s4-factored-example-full.log`, `verification/s4-factored-study-python.log`, `verification/s4-factored-boundary-python.log` |
| Factored S4 capacity guard preparation | Shared finite-domain model builder and cgroup-guarded, no-network construction-only worker; synthetic success, timeout, memory-kill and malformed-receipt outcomes retained. No N83 model construction was run. | `verification/capacity-guard-smoke/`, `sat_capacity_supervisor.py` |
| Linux capacity worker cross-build and startup | PASS cross-build with native Redis TLS disabled only for this isolated worker; guarded startup reached the deliberately empty manifest and reported `PRODUCER_FAILURE_worker_exit`. This is a malformed-input check, not a retained-base construction. | `verification/sat-capacity-linux-build.log`, `verification/sat-capacity-linux-build-command.txt`, `verification/capacity-worker-startup/` |
| Post-capacity-builder release and Python suites | 2,242 library tests passed with 94 ignored; N83 example 15 passed; study and boundary Python 5 and 16 passed | `verification/sat-capacity-lib-full.log`, `verification/sat-capacity-example-full.log`, `verification/sat-capacity-study-python.log`, `verification/sat-capacity-boundary-python.log` |
| Native finite-grid audit | All 45,360,000 dispositions accounted for | `verification/design-audit.log`, `design.json` |

Earlier build and test outcomes remain retained. Initial dependency resolution failed under restricted networking; `--locked` could not be used before this older baseline resolved a lockfile. Some build attempts were interrupted during baseline recovery. The first complete isolated library run had two loopback-network permission failures. The next run passed those tests and failed the pre-existing randomized SQIsign wrong-message assertion. Its focused replay passed; the subsequent full suite passed. The SQIsign source is unchanged, and the earlier failure remains visible rather than being relabeled as a pass.

The post-rank-screen sandboxed full suite again reached the same two loopback bind denials (2,233 passes, 2 failures; `verification/lib-tests-rank-screen.log`). Rerunning the required full suite with loopback access passed all 2,235 tests (`verification/lib-tests-rank-screen-escalated.log`). This recheck changed no Rust source or factor-base object.

Example tests exercise finite coverage and ordinal bounds, the same complete set under binary/Gray enumeration, field encoding rejection, separate subgroup moduli, generic arithmetic replay, corrupted compressed bytes, semantic coefficient corruption with recomputed byte hashes, and point-only relation discovery, negative coverage and budget censoring.

Generic-versus-wide arithmetic replay uses a separate process on the same host and repository. It checks all point and coefficient records, public candidate indices, subgroup membership, source projection, eigenvalue labels, closure and distinctness. External independent-host replay remains pending. This construction panel supplies no boundary-ledger promotion or canonical performance finding.

## Remaining comparison work

The imported primary base still needs an end-to-end solver run with an audited wide-field relation oracle, full subgroup-modulus handling, natural rank, verified column logarithms and target extraction. An available BigUint elimination routine does not by itself establish those coupled phases. Large-prime graph adapters, WDSat sealed ANF capacities, FES full-system verification, cofactor-aware subspace lifting and domain-preserving symmetry remain explicit gates in `protocol.json`.

`SOLVER_GATES.md` pins the immediate source and capacity blockers. At `n=l=83`, the original native S4 encoder needs at least **593,364 SAT variables and 2,349,149 AND-definition clauses** from its unreduced x/e correspondence and default ordering. An experimental factored encoder now has a source-derived ceiling of **158,032 variables and 469,881 AND-definition clauses**, before finite-base domain constraints and ordering clauses. Its small-system equivalence and explicit-orbit group-lifting checks pass, but no N83 model construction, capacity or lifting receipt exists. The exact wide coordinate trie has a focused regression. The WDSat, FES and double-large-prime modules in the inspected main snapshot are absent from this isolated branch and cannot directly represent the primary N83 system as implemented there. A further enumeration cap would add another censored observation without addressing those backend gates.

After the pilot, the exact enumerator gained a `u128` batched-add path for degree 83. Focused generic-arithmetic and retained-primary-base witness-order tests pass. These tests establish correctness for the checked cases; no cold workload has been rerun on the new source, so all earlier capped timings remain pinned to their original implementation. The exponential enumeration count and the missing SAT, WDSat, FES and large-prime adapters still prevent a total-runtime factor-base ranking.

Select a base only after complete, matched one-target cold runs, natural rank gain, verified column logarithms and an individual logarithm, all failed work charged, disjoint holdouts, A/A controls, at least five paired rounds, admitted uncertainty and independent validation. The local L0 observations do not establish that minimum. Literal coverage of every possible point set and parameter value is outside this finite experiment.
