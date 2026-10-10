# N83 solver availability and capacity gates

This is a source audit, not a timed solver result. The backend-availability comparison distinguishes the inspected main snapshot `c70c32d486a3ac7531fe27f7193d9f09caa58344` from the isolated study branch at `001acc47870e935a17f4d8155c2cb934ace58284`. The latter builds the retained primary base and exact enumeration runner; subsequent study-branch code adds a wide-word enumeration path. The inspected main snapshot contained additional backends absent from that pinned study revision. Reproduce this historical distinction with `git ls-tree -r --name-only REV -- src/cryptanalysis` at each pinned revision; `verification/solver-gates-source-check.log` records the corresponding blob identities and source-width assertions. The current integrated-source status is recorded at the end of this document. A path in the prior-art inventory is not proof that an N83 run can execute it.

| Backend | Pinned source fact | N83 consequence |
| --- | --- | --- |
| Exact enumeration | `koblitz_index_calculus.rs::enumerate_decompose` now uses the existing `u128` field path through degree 127, with a batched affine add for ordinary pairs. The retained-base test checks first-witness agreement with a separate generic recursion; the cold runner retains a hard wall cap. | Executable reference oracle with `|F|^(m-1)` search cost per target. The earlier capped receipts predate this code path and do not measure its runtime or establish a complete primary runtime. |
| Incremental relation rank | `koblitz_relation_solver.rs::RelationSolver` selects the native-word echelon for small subgroups and exact BigUint echelon for the 81-bit primary subgroup. The index-calculus driver feeds each accepted relation to that solver, so independent/dependent counts and a target-column pivot have the same meaning at both widths. `WideRankTracker` gates factor-base-log precomputation on full coefficient rank and its report counts dense solve attempts. High-limb planted systems and coefficient-rank prefixes are checked against dense elimination. | A wide target is reported only when its column is algebraically determined; rank-deficient log collections skip the dense solve and return no purported table. These are correctness and accounting paths. The retained N83 base has produced no natural rank or column-log verification receipt. |
| Original native SAT S4 | `binary_semaev_s4.rs::S4System` uses `3l` x variables and `6l-3` symmetric-function variables. `semaev_sat.rs::encode_semaev_s4_with` adds `2l` lex-order auxiliaries by default and one auxiliary per distinct non-linear ANF monomial. The degree-83 coordinate-domain trie now retains all bits. | At `n=l=83`, this encoder needs **at least 593,364 SAT variables and 2,349,149 AND-definition clauses**, before e-space auxiliaries, parity machinery and three finite-coordinate-domain tries. |
| Experimental factored SAT S4 | `encode_semaev_s4_factored_with` shares the three pair-product coefficient arrays and builds the e-space equations without materializing the cubic x-ANF. It accepts native XOR only. `SatDecompositionOptions::factored_s4` now selects it for ambient-basis three-summand unions; other shapes return inconclusive. Complete small x-assignment equivalence, a positive n=15 scalar S4 round trip and a small explicit-orbit group-lifting comparison pass. | At `n=l=83`, source structure gives **at most 158,032 SAT variables and 469,881 AND-definition clauses**, before finite-domain constraints and ordering clauses. No N83 construction, memory, lifting or cold-run receipt exists. This encoding is outside the frozen v1 grid and the primary CLI still gates `sat-m3`. |
| Balanced-S5 SAT example | `koblitz_s5_sat_instance.rs` contains four-summand native-XOR, orbit-factorized and pair-domain encodings, but admits only n=7,11,13,17,19,23,37,41,53, with `u64` point/label/rank paths and `pair_sum_trie` capped at n<=53. | Requires full-width field, point, domain and rank transfer with original group/model replay. It is not an executable N83 SAT arm on this branch. |
| Wide compact-orbit four-sum runner | `koblitz_orbit_dlp_fast.rs` has a u128 S3 index at n=83 but uses `KoblitzCurve::new`'s polynomial basis, a different point-defined-base JSONL schema and u64 subgroup order/rank arithmetic. Historical n=83 `a=1` online clocks exclude reusable preparation. | Pinned study objects use a different polynomial basis; the primary a=0 order is 81 bits. See `S5_CAPACITY.md` for the source-derived `83K²` index-state count and guarded-transfer sequence. No stored N83 base has run through this producer. |
| WDSat, inspected main only | `wdsat_oracle.rs::model_to_u64` returns `None` for a model longer than 64 bits. Its ANF comes from the single-word `F2BoolPoly` path. | An ambient three-summand model has 249 x bits; the pinned adapter cannot lift it. A split/reduced representation needs a sealed complete-system verifier and capacity receipt. |
| FES and Gray, inspected main only | `mq_fes.rs` requires quadratic forms. Its Möbius routine admits at most 24 variables and 64 equations; Gray variants have their own small limits. `mq_monica.rs::monica_search` refuses more than 40 variables or 128 equations. | The full ambient S4 model has cubic x/e correspondence and at least 910 core variables. These kernels can only serve a justified split/filter with verification against the complete original equations. |
| Double large primes | Inspected main's `ecbench_large_prime.rs` imports the field modulus, points, subgroup order and cofactor through `u64` corpus fields. The isolated branch has a separate `koblitz_large_prime.rs` exact BigUint eliminator with signed-Frobenius residual keys, source coefficients, group replay, explicit caps and a full-rank matrix bridge. | The bridge accepts verified subgroup-point partials and verifies solved orbit and target logs. It does not discover partials or run the retained N83 base. A source-pinned producer, cofactor-policy decision, measured graph/rank behavior and matched cold run are still required. |

The retained signed-Frobenius bases at K=64, 256 and 600 have 5,312, 21,248 and 49,800 distinct x coordinates respectively (`83K`, from the full-orbit and point-distinctness replay checks; the K=64 and K=256 cold imports also check this invariant). Each SAT S4 instance constrains three 83-bit coordinate blocks to one such finite set. These counts describe input size, not solver time or memory use. The current finite design still addresses every declared tuple, while its backend dispositions remain recipes until their adapters and original-system checks are executable.

The original S4 bound follows from its unreduced x/e correspondence. Write each summand coordinate as `X_i(z) = sum_j x_(i,j) z^j`, with three disjoint blocks of `l` Boolean variables. The three pair products in `sigma2 = X0 X1 + X0 X2 + X1 X2` contribute `3l²` distinct quadratic x-monomials. `sigma3 = X0 X1 X2` contributes `l³` distinct cubic x-monomials. No field reduction occurs before these correspondence rows, and disjoint variable blocks prevent cancellation. `encode_semaev_s4_with` allocates one auxiliary for each such monomial and emits three AND-definition clauses per quadratic or four per cubic. With default ordering, its source-derived lower bounds are `11l-3 + 3l² + l³` SAT variables and `9l² + 4l³` AND-definition clauses. For `l=83` these are **593,364 variables** and **2,349,149 clauses**. This bound is specific to the original encoder.

The factored correspondence uses `3l²` pair-product AND gates, `3(2l-1)` shared pair-coefficient variables, and `l(2l-1)` further AND gates for the third product. The e-space S4 equations are quadratic, so they need at most `C(6l-3,2)` distinct e-monomial auxiliaries. Including the same `11l-3` core and default ordering variables gives at most **158,032 variables** at `l=83`; three clauses per two-input AND give at most **469,881 AND-definition clauses**. Lex-order clauses, finite-base domain constraints and their auxiliaries are additional. `verification/s4-factored-bound-check.log` pins the source assertions and arithmetic; focused equivalence tests are in `verification/s4-factored-focused.log`. Neither bound is a measured N83 memory or runtime result.

The next runnable gate is a bounded, source-pinned construction of the factored N83 S4 model with the retained K=64 coordinate domain. `build_union_s4_encoding` now shares the exact domain and trace construction with the solver, and `sat_capacity_supervisor.py` prepares a Linux worker in a no-network Docker container with a hard cgroup memory ceiling, zero swap and a wall deadline. Its success, wall-cap and memory-cap paths passed synthetic guard checks (`verification/capacity-guard-smoke/`); no retained N83 model was constructed in those checks. The real gate still needs an authorized run and a retained-base construction receipt. Any subsequent solve needs independent scalar S4 and group checks of lifted models, with `UNKNOWN` on budget exhaustion. Only then can this experimental encoding enter a versioned sweep and a matched one-target cold comparison. WDSat and FES need wide adapters and original-system validation; the large-prime arm needs a source-pinned partial producer and measured graph/rank behavior. Spending more time on the existing exact enumerator alone would not rank factor bases.

The local Docker server exposed 7,529,156,608 bytes of VM memory at the guard check, less than the 48 GiB macOS host. The static x86-64 Linux worker ran inside an arm64 container under emulation. Its future construction wall time is therefore only a resource-capacity diagnostic; it cannot be used as a matched macOS index-calculus timing or a factor-base ranking.

## Current-main applicability audit (2026-10-09)

The refreshed `origin/main` source snapshot is
`09750b8108b2cb1e99bac24069618461b3d7936e`; the isolated study branch
was `e276b4f2dc3da51bd454c5a490ff0cde60f0ae97`. This is a source audit,
not an N83 solver measurement. The main-branch
`koblitz_factor_base_search.rs::measure_solve_cost` explicitly accepts only
`FactorBaseDomain::LinearSubspace`: its Weil-restricted solver prices the
entire span, not an explicit orbit subset. The main-branch
`RESEARCH_FACTOR_BASE_SOLVE_COST.md` §6 reports a 22.41× difference in
Gröbner word-XOR cost between trials- and solve-cost-selected bases across
72 complete, verified runs on K_1/2^15 at m=2. That is a solver-stage
counter on a different domain and cannot be imported as a cold-runtime ratio
or as a score for the stored N83 bases.

The WDSat exporter calls `build_decomposition_system_reusing` and converts
its single-word `F2BoolPoly` to ANF; its model lift stops above 64 bits.
`mq_fes.rs` accepts quadratic systems with at most 24 variables and 64
equations for its Möbius path, and refuses chained cubic m>=3 Semaev systems.
`semaev_higher.rs` handles prime-field short-Weierstrass curves, not the
binary Koblitz curve. `large_prime_filter.rs` stores modular coefficients in
`u64`, and `koblitz_sparse_la.rs::modulus_supported` permits at most 63-bit
orders; the primary order is 81 bits. These modules supply design prior art,
but no directly executable N83 explicit-orbit higher-arity/large-prime/sparse
pipeline. The study-local BigUint row adapter supports checked algebra and
small fixtures; it still lacks a retained N83 partial producer.

A dry `git merge-tree --write-tree --name-only --messages` of those two pinned
revisions identified content conflicts in `Cargo.toml`,
`src/bin/ic/experiment.rs`, `src/cryptanalysis/binary_semaev_s4.rs`,
`src/cryptanalysis/koblitz_index_calculus.rs`, `src/cryptanalysis/mod.rs` and
`src/cryptanalysis/semaev_sat.rs`. It did not edit the worktree. A deliberate
integration must resolve and retest those contracts before importing a
current-main backend. The next runtime-relevant implementation gate is a
full-width producer for m>=5 or verified partial relations, with a complete
finite-base model and independent group replay. A guarded K=64 four-sum
probe remains useful for capacity and correctness, but the one-row uniform
rank screen in `RANK_QUERY_SCREEN.md` rules out treating it as a likely
full-rank comparison within a short pilot. This is a necessary-condition
screen for that oracle model, not a lower bound on every IC algorithm.

### Refreshed upstream and semantic merge gate (2026-10-09)

At refreshed `origin/main` `416dffe827704cb8b9ee5e61daa4bb3b46cce763`,
the Git blobs for `wdsat_oracle.rs`, `mq_fes.rs`, `mq_monica.rs`,
`ecbench_large_prime.rs`, `koblitz_index_calculus.rs`,
`koblitz_factor_base_search.rs`, and `koblitz_sparse_la.rs` are byte-identical
to the `09750b8` snapshot audited above. The width and system-shape gates in
that audit therefore still apply. The current main driver has WDSat,
MQ-FES, Crossbred and batched relation-attempt branches; their presence on
main does not make them executable on this isolated study branch or establish
N83 throughput.

An actual no-commit merge into study-branch `f1a154026` confirmed six content
conflicts: `Cargo.toml`, `src/bin/ic/experiment.rs`,
`src/cryptanalysis/binary_semaev_s4.rs`,
`src/cryptanalysis/koblitz_index_calculus.rs`, `src/cryptanalysis/mod.rs`, and
`src/cryptanalysis/semaev_sat.rs`. The main driver changed by 16,005 insertions
and 2,837 deletions relative to the merge base, while the study driver changed
by 854 insertions and 85 deletions. Git aligned a moved decomposition block
twice, producing duplicate `groebner_decompose` and SAT definitions in the
unresolved file. Applying the study driver's patch over the main file left 33
of 51 hunks rejected. These are source-integration diagnostics, not timings or
solver verdicts. The merge was aborted and the study branch was restored clean.

The integration must use the current main driver as the structural base and
port the study changes by contract: exact multi-limb point identities and
subgroup coefficients; `RelationSolver`/`WideRankTracker` rank and column-log
gates; the guarded five- and six-summand `ChainedS3` strategy; then the
S3-bound `primary-chain-cold` adapter and its source/input supervisor. Main's
coordinate-domain trie already uses `u128` codes, so its 83-bit paths should
be retained and checked against distinct coordinates that share their low
64 bits. WDSat, MQ-FES, Crossbred and batched attempts must remain reachable
after adding `ChainedS3`. A merged binary needs a fresh source hash, Linux
build and startup receipt, release library and touched-example tests, the
study and boundary Python suites, and independent replay before any N83
backend or runtime claim. This source integration uses no part of the
exhausted one-hour N83 experimental allowance.

### Source integration status (2026-10-10)

The isolated integration branch merges upstream `845fd323015e28d32b73fa25e551f830d5826e6c`
with the study code using the upstream Koblitz driver as the structural base.
The six content conflicts were resolved by contract. The merged library
compiles, the 83-bit coordinate trie retains high bits, and the main-branch
WDSat, FES, Crossbred and batched dispatch remain in the source. Their width
and system-shape limits above still apply; an available enum variant is not an
N83 solver result.

The generic driver now feeds full-width coefficients to `RelationSolver`, and
its chained-S3 branch requires explicit nonzero variable, domain, model and
conflict caps before search. Its zero-cap path retains an inconclusive report.
The one-shot factor-base-log path samples the full subgroup and checks every
relation in the original group before exact coefficient-rank gating and dense
solve. One-shot descent likewise samples full-width `a,b`, retains decimal
coefficients in its witness and verifies the resulting logarithm in the group.
The older distributed work-unit collector, streamed log solver and observed
descent ledger still encode scalars as `u64`; they reject a subgroup wider than
64 bits. Sparse LA still rejects the 81-bit primary modulus. A deterministic
N83 one-orbit known-answer test exercises the split precompute/descent
arithmetic, while the retained S3 bases have not supplied a natural relation
or complete single-target cold runtime. The source merge does not promote a
factor-base winner or reset the exhausted local-pilot clock.
