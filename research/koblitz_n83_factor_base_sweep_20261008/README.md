# Factor-base sweep for degree-83 Koblitz index calculus

Date: 2026-10-08. The primary objective is the lowest fully charged index-calculus runtime for one public known-answer target. Construction, orbit/index preparation, failed relation attempts, lifting, partial-relation processing, linear algebra, individual logarithm, verification and artifact I/O all belong in the cold cost. Reusable online-only cost is a separate secondary measurement. A candidate is a factor base together with its oracle, solver and resource configuration; a base cannot be selected by cardinality or relation yield alone.

The exact primary model is `icv1-f2m83-tm6151469093347-debefd74`, with equation `y^2 + xy = x^3 + 1`, polynomial-basis modulus `z^83 + z^45 + z^2 + z + 1`, subgroup order `2417851639230796216685689`, and cofactor 4. Its generator is pinned in the native exporter. The diagnostic model is `icv1-f2m83-t6151469093347-cdcc5432`, with `a=1`, the same modulus, subgroup order `8569786107849059`, and cofactor `1128547018`. These are different workloads: an 81-bit subgroup and a 53-bit subgroup. The degree alone does not identify a workload. The slugs were checked against the registry at `c70c32d486a3ac7531fe27f7193d9f09caa58344`; exact representations, generators and subgroup parameters are retained in every exported header.

## Requested scope and completion conditions

| Requirement | Deliverable or completion gate |
| --- | --- |
| Review prior art | Source and evidence inventory in PRIOR_ART.md, pinned to the inspected revision |
| Every combination | Complete ordinal-addressable finite Cartesian grid in design.json; literal universal coverage remains outside this finite design |
| Splitting, symmetry, WDSat, Gray codes, FES, double large primes, Frobenius | Explicit axes, compatibility dispositions, capacity checks and adapter gaps; a named option is never evidence of execution |
| Store factor bases in S3 | Full point sets, subgroup metadata, source points, orbit labels, content hashes, cross-backend replay and verified upload/download receipts |
| Rigorous empirical pilot | Public-data construction panel, exact counted candidates and point counts, preserved failures and budget, separate replay process |
| Best total runtime | Complete matched pipeline measurements with independently verified answers, disjoint holdouts and admitted uncertainty; remains unresolved until that evidence exists |

## The finite search

The native `plan` command enumerates the complete compatibility audit of **45,360,000 tuples**. The design has 1,050 conditional base specifications and 43,200 solver configurations per base. No tuple is silently dropped. `case ORDINAL` reconstructs any tuple in `0..45360000` without writing millions of duplicate JSON rows. Disposition counts must sum to the full population; the executable's test checks that invariant.

Base families are sequential public x, hash-derived public x, Gray-prefix public x, polynomial subspaces, random linear subspaces, affine subspaces, rational 2-torsion coordinates, rational 4-torsion coordinates, and Frobenius-invariant module blocks (named `frobenius_divisor` in the grid). Orbit-column counts are 32, 64, 128, 256, 600 and 900. Subspace dimensions are 4, 8, 12, 16, 20, 28, 41 and 82. Invariant-block dimensions are 1, 82 and 83; 82 is a module dimension, not a divisor of the field degree or a subfield degree. Each applicable family crosses the frozen three public seeds and closure policies `none`, `sign` and `signed_frobenius`. Invariant-block bases have a deterministic seed and signed-Frobenius closure.

The solver axes cross summands 2 through 6; split bits 0, 4, 8 and 12; no symmetry, summand ordering, signed Frobenius, 2-torsion and 4-torsion; compact S3 four-sum, MITM, F4, SAT-CDCL, WDSat, Gray FES, Moebius FES and Monica FES; binary or Gray enumeration; zero, one or two large primes; 1, 4 or 12 threads; and dense single-word, sparse single-word or BigUint modular linear algebra. In solver names, S3 denotes the third Semaev polynomial; artifact storage uses Amazon S3.

This is exhaustive over those declared choices. All other point sets, bases of subspaces, affine shifts, field representations, solver versions, split heuristics, restart schedules, numeric values and seeds remain outside it. A finite pilot cannot establish a universal optimum. Adding an axis or level requires a new version of the frozen design, with its own hash and coverage audit.

## Structural and implementation gates

Direct modular computation gives `ord_83(2)=82`. The nontrivial factor of `X^83-1` has degree 82 and is irreducible over GF(2). Since 83 is odd, the Frobenius module is semisimple; invariant linear-subspace dimensions are exactly 0, 1, 82 and 83. The small dimension-one points are rational torsion and disappear under the primary cofactor projection. This obstruction applies to invariant linear subspaces. It does not obstruct a small union of point orbits or the closure of a general subspace. The computation is pinned by a native test and is a derived structural statement, not a timed experiment.

The current large-prime adapter represents fields, points and subgroup orders with single words, so the degree-83 arm needs an adapter before combining it with zero/one/two-large-prime experiments. The exact double-large-prime model is two residual point/orbit columns, not integer factorization of a point coordinate. Every partial relation must retain exact modular coefficients, signs and Frobenius phases; cycles must be group-checked before filtering into the main matrix. Count duplicate partials, singleton vertices, cycle closure, cycle length, resulting matrix density and rank gain; charge graph construction, filtering and unsuccessful partials.

FES kernels are quadratic Boolean-system solvers with explicit variable/equation capacity limits. Degree-83 systems must be reduced or split into a represented filter and a complete original-system verifier. Higher-arity chained systems are not automatically quadratic. Gray code is an enumeration method inside FES and subspace construction, so treating Gray and FES as independent multiplicative improvements is unjustified. On the same complete subspace, Gray and binary enumeration must produce the same set; the native test checks that control. Truncated first-K prefixes can produce different bases and are recorded as different construction policies.

WDSat needs a frozen ANF, derived static allocation capacities, a sealed binary/configuration hash, model verification and explicit SAT/UNSAT/UNKNOWN outcomes. Splitting must preserve branch coverage; exhausting a conflict/node/time budget is UNKNOWN. Summand-order symmetry is applicable only to exchangeable domains; Frobenius symmetry requires the actual domain to be closed under it. Torsion symmetrization changes the coordinate domain and must be checked separately, including excluded and exceptional points. The exact primary subgroup order exceeds u64; single-word modular linear algebra is inadmissible there. A BigUint option in the design is an unexecuted recipe until the complete backend integration is verified.

## Pilot and storage

The authorized local pilot is capped at 3,600 seconds. Its construction panel has 54 actual bases: two curve arms, three public-x policies, three seeds and K=64, 256 or 600 signed-Frobenius columns. Each complete base has exactly `166*K` distinct subgroup points. The exporter retains each source point before cofactor clearing and every expanded point with column, phase, sign and exact coefficient. Hash-derived public x sampling uses a domain-separated BLAKE3 stream and no known base logarithms. Sequential and Gray policies scan public polynomial words; the post-projection coordinates are explicitly not asserted to lie in the seed subspace.

Gzip objects are named by the BLAKE3 hash of their exact compressed bytes, with a second hash for the uncompressed JSONL. A third hash covers sorted point coordinates, independent of metadata, row order and orbit representative choice, to detect duplicate bases. Accepted candidate indices bind source points to the exact public construction policy. Full objects live under `s3://crypto-autoresearcher/factor-bases/icv1/etc/koblitz_n83_factor_base_sweep_20261008/a0/objects/` and the corresponding `a1/objects/`. Manifests and replay receipts live under a content-addressed `panels/` prefix. The uploader accepts only the declared destination, requires a matching successful replay receipt, rechecks local hashes, then downloads every uploaded object and verifies its bytes. Local full objects are excluded from Git; manifests, validation receipts and research source are committed.

The replay uses the generic multi-limb `BinaryCurve` implementation against the producer's `Gf2_128` arithmetic. It checks source and representative membership, projection, subgroup order, Frobenius eigenvalue, every expanded point and coefficient, distinctness, closure and exact counts. This is a separate process and arithmetic implementation on the same host and repository. Independent-host validation remains an additional obligation. Construction elapsed times on macOS have L0 status and are descriptive only; they cannot rank runtime candidates.

The bounded `probes` command follows successful replay. It generates 32 domain-separated public fixtures per curve arm and gives its exact two-summand complement oracle point inputs only. Validation scalars stay in `probe-validation.json`. Each base uses the same arm's public corpus, with complete failed lookups counted and UNKNOWN retained on a cap. Generic arithmetic verifies every returned relation. The oracle supports repeated summands; its tests cross-check a positive and an exhaustive negative case with generic arithmetic and check censoring. These rows share a dictionary across the 32 stage probes and cannot supply independent cold one-target totals. Shared corpora, nested bases and duplicate sets also make rows dependent observations. Higher-arity, WDSat, FES, splitting and large-prime execution coverage remains separate.

`compact_cold.rs` pins the diagnostic degree-83 `a=1` curve to the exact exported modulus, generator, 53-bit subgroup and cofactor. It adapts the repository's archived compact four-sum wide driver to a retained S3-backed base and fixture zero of the public corpus. The adapter regenerates the hash-derived base from the frozen seed and checks every point and coefficient against the replayed object before rank work starts. Its outer process caps rank and target work, retains a censored result on expiry, and checks each solved column logarithm with generic group arithmetic. This complete diagnostic pipeline keeps the 81-bit primary arm separate; its older rank logic still uses single-word scalars.

Ordered-pair K=256 and K=64, plus exploratory unordered-pair K=64, each reached a 600-second cap before a complete summary; all three `cap.json` files retain `UNKNOWN_budget`. The `unordered` pair-index mode uses summand exchange and Frobenius canonicalization to avoid redundant pair states. A test checks that every canonical root key in a small ordered index remains in the unordered index. It is an added implementation experiment outside the frozen 45,360,000-tuple v1 grid. Its K=64 index completed in 0.771 seconds, but neither that stage nor the key-set test establishes a total-runtime gain. The three cold cases and approximate 3,530.848-second pilot budget audit are in `pilot-01/cold-diagnostics.json`.

`SUPPORT_MOMENTS.md` and `pilot-01/support-moments.json` add a post-pilot exact support-count screen for every retained K and m=2..6. Its uniform-target Markov ceilings concern the existence of full smooth decompositions, not the fixed public fixtures or solver runtime. At the primary a=0 K=600 base size, the m=4 ceiling is about 1.70e-6 and m=5 is about 0.0338; m=6 becomes vacuous. This directs future full-smooth work toward arity/size regimes with non-negligible potential coverage, while the actual cost of the needed solver and large-prime graph variants remains unmeasured.

`RANK_QUERY_SCREEN.md` and `pilot-01/rank-query-screen.json` derive a further necessary condition for a uniform-target oracle returning at most one full-smooth row per query. Even allowing dependent queries, the primary K=600 four-summand base requires at least 176,888,106 queries before the first-moment bound can permit a 50% chance of 600 rank rows; five summands require at least 8,880. These are exact query-count gates under the stated model, not predicted completion times. They do not apply to biased sampling, multirow solvers or large-prime partials. The current capped driver has not established uniform full-width primary sampling.

The generic Koblitz IC driver now draws full-width scalars for subgroup orders above 64 bits. Its additive randomizer ranges uniformly over the entire subgroup, so `[a]G+[b]Q` is uniform for fixed nonzero `b`; its nonzero coefficient uses the full order. The existing at-most-64-bit random stream is unchanged. The factor-base-log precompute also disables its single-word rank gate for a wide modulus and falls back to BigUint elimination. This fixes two primary-arm adapter prerequisites, but the retained S3 bases, higher-arity solvers and complete cold pipeline are not yet integrated or measured on that path.

## Subsequent total-runtime experiment

Freeze the exact public target corpora, resource envelope and native backend versions before the first relation measurement. Use disjoint tuning and holdout targets, with validation-only known-answer scalars kept out of solver input. Keep one result per independent one-target workload. Preserve every phase cost, operation unit, cap, failed attempt, OOM and timeout. Baseline and candidate receive the same target and budget; randomize/interleave their order with A/A controls. Require at least five paired rounds and a paired 95% interval outside the measured noise floor for a runtime improvement. Run the admitted measurements on an L2-capable host; local macOS L0 timings cannot establish that comparison.

Screen by verified natural relation yield, rank gain per charged unit, matrix density, partial-graph behavior and memory. All screening metrics are secondary and censored failures remain visible. Fully charge survivors through linear algebra and one individual logarithm, verify the answer in the group and cross-check the validation sidecar. Confirm finalists on disjoint holdouts and an independent host. Preserve a Pareto set for total charged work, total wall, memory and success probability; select the requested runtime winner only among complete, comparable, verified runs. Leave the winner null when an arm has missing phases or unverified completion. Rho is a matched reference, not a substitute objective or an online-only denominator.

## Reproduce

```sh
export CARGO_PROFILE_RELEASE_CODEGEN_UNITS=256
cargo test --release --lib
cargo test --release --example koblitz_n83_factor_base_export
cargo build --release --example koblitz_n83_factor_base_export
target/release/examples/koblitz_n83_factor_base_export plan research/koblitz_n83_factor_base_sweep_20261008/design.json
target/release/examples/koblitz_n83_factor_base_export case 45359999
target/release/examples/koblitz_n83_factor_base_export pilot research/koblitz_n83_factor_base_sweep_20261008/pilot-01 3600
target/release/examples/koblitz_n83_factor_base_export replay research/koblitz_n83_factor_base_sweep_20261008/pilot-01
python3 research/koblitz_n83_factor_base_sweep_20261008/support_moments.py
python3 -m unittest discover -s research/koblitz_n83_factor_base_sweep_20261008 -p test_support_moments.py -v
python3 research/koblitz_n83_factor_base_sweep_20261008/rank_query_screen.py
python3 -m unittest discover -s research/koblitz_n83_factor_base_sweep_20261008 -p test_rank_query_screen.py -v
target/release/examples/koblitz_n83_factor_base_export probes research/koblitz_n83_factor_base_sweep_20261008/pilot-01 600
target/release/examples/koblitz_n83_factor_base_export upload research/koblitz_n83_factor_base_sweep_20261008/pilot-01
target/release/examples/koblitz_n83_factor_base_export cold research/koblitz_n83_factor_base_sweep_20261008/pilot-01 64 600
# Exploratory pair-index mode, with a fresh immutable attempt directory:
target/release/examples/koblitz_n83_factor_base_export cold research/koblitz_n83_factor_base_sweep_20261008/pilot-01 64 600 unordered
python3 research/koblitz_n83_factor_base_sweep_20261008/summarize_cold.py
```

The recorded release profile retains optimization level 3 and uses 256 codegen units, with four Cargo build jobs. `verification/Cargo.lock` preserves the resolved dependency versions. Artifact-writing commands require new paths and refuse to overwrite completed evidence. The native tests also check invalid encodings, ordinal bounds, separate subgroup sizes, corrupted compressed artifacts and semantic coefficient corruption with recomputed byte hashes. See RESULT.md for actual commands, revisions, outcomes and remaining requirements. The study diagram is `design.svg`; its editable source is `design.mmd`. No canonical performance graph changes unless a new admitted performance finding is established; this pilot does not supply one.

For exact replay of the pinned panel, use a new directory because `replay` creates an immutable receipt. Restore both S3 arm object directories there. The committed manifest names and hashes every required object, and `replay` rejects wrong compressed or uncompressed bytes:

```sh
n83_replay_dir=$(mktemp -d /private/tmp/n83-replay.XXXXXX)
cp research/koblitz_n83_factor_base_sweep_20261008/pilot-01/manifest.json "$n83_replay_dir/manifest.json"
aws s3 sync s3://crypto-autoresearcher/factor-bases/icv1/etc/koblitz_n83_factor_base_sweep_20261008/a0/objects/ "$n83_replay_dir/objects/"
aws s3 sync s3://crypto-autoresearcher/factor-bases/icv1/etc/koblitz_n83_factor_base_sweep_20261008/a1/objects/ "$n83_replay_dir/objects/"
target/release/examples/koblitz_n83_factor_base_export replay "$n83_replay_dir"
```
