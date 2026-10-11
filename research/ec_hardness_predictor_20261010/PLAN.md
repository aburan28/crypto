# Geometric hardness predictor for ECDLP across isogeny classes: end-to-end plan

Frozen 2026-10-10. Status: **plan; nothing run.** Result class when run:
**level-1 learned-correlation study with level-2 exact certificates, not an
ECDLP speedup.** `S`, end-to-end cost and speedup are **unset** for every
artifact this plan produces and stay unset whatever a model predicts. A
model prediction is a lead; a ledger row moves only through a frozen
paired `ecbench` challenge (`AGENTS.md`, `docs/ic/boundary_targets.json`).

This document is the hand-off to the executing agent. It is meant to be
followed without this conversation: every stage names the producer, the
file it writes, the check that gates it, and the number it is budgeted at.

Companion repositories:

| repo | role | path |
|:--|:--|:--|
| `crypto` (this branch `feat/elliptic-curve-complexity-ml-4b4228`) | Rust producers: class enumeration, exact invariants, counted labels | `/Volumes/SSD990/crypto/elliptic-curve-complexity-ml-4b4228` |
| `ml-cryptanalysis` (`github.com/aburan28/ml-cryptanalysis`) | Python learner: dataset join, splits, nulls, models, reports | `/Volumes/SSD990-2/ml-cryptanalysis` |

---

## 1. Question, hypotheses, and what a result would be

**Question.** From the algebraic data of a curve and its position in its
isogeny class, namely the endomorphism order `End(E)` and its conductor,
the Frobenius order `Z[π]` and the conductor evaluations `v_ℓ(f_π)`, the
rational and extension-field torsion structure, the automorphism group, and
the local isogeny neighbourhood, can a learned model predict the algebraic
index-calculus labels (relation yield `γ`, first-fall degree `d_ff`,
refutation/solving degree `D*`, semi-regular degree, Macaulay rank profile,
solver work, whole-pipeline `S_IC`) better than the size-only baseline, on
classes and field sizes it never saw?

**The theorem that fixes the prior** (`research/isogeny_class_difficulty_20261008/README.md`):
a separable `ℓ`-isogeny with `ℓ ≠ r` is a group isomorphism on `E[r]`, so
every generic quantity is constant across the `F_q`-isogeny class and only
three things can vary: transport cost, `|Aut|` (only at `j ∈ {0, 1728}`), and
**model-dependent non-generic structure**, which is exactly what index
calculus reads off the equation. The labels this plan predicts are all of
the third kind, so variation inside a class is *possible* and is the thing
to measure; variation of a class invariant (`#E`, `r`, `t`, `D_π`,
embedding degree) is a bug.

Hypotheses, registered now:

- **H0 (class invariance).** For a fixed factor-base family and dimension,
  each label is constant across a class up to the A/A interval, and a
  model's held-out skill above the size baseline is explained entirely by
  `(log q, log r, ℓ_fb, m, |Aut|)`. Expected on prime fields
  (`PROTOCOL-I4`, H4-null) and on prime-degree binary fields; known at
  `n = 17, 19` for Macaulay rank and first-fall degree of isogenous
  neighbours (`ml-cryptanalysis/docs/PLAN.md`, measured obstructions).
- **H1 (model-dependent effect).** Some label varies with features that are
  not class invariants: volcano level `v_ℓ(f_End)`, `|Aut|`, the model's
  abscissa distribution relative to the factor base, subfield-definedness
  of `j` or coefficients over extension fields, the GHS magic number over
  composite-degree binary fields (known positive control), torsion over
  `F_{q^k}` for small `k`.
- **H2 (neighbourhood signal).** Information in the `ℓ`-isogeny
  neighbourhood (neighbour multiset, ball growth, crater length, spectral
  data) predicts a label beyond the node's own invariants. The Cayley-graph
  argument of §9.3 says this is impossible for class-invariant labels and
  possible only through model-dependent node features; the GNN ablation
  tests exactly that.

What a result is: a held-out R² above the size baseline that survives the
isogeny-class split, the field-family split, the permutation null, the
structure-blind control label, and Benjamini–Hochberg across every label
tested in the round, *and* replicates on a fresh seed over a field size the
search never saw. Anything else is a "nothing here" with a stated minimum
detectable effect (§10.4). Both outcomes are deliverables.

---

## 2. What exists today (read these before building anything)

### 2.1 `crypto` (Rust; counted, certified)

| capability | where | state | relevance |
|:--|:--|:--|:--|
| `F_p`-class walker: BFS over `Φ_ℓ mod p`, `ℓ ≤ 61`, kernel-certified `IW1` edges, `curves.yaml` / `isogeny_routes.json` / `walk.json`, class audits, trait detectors, `plan`/`traits`/`collect` sharding over `taskq`, million-curve streaming store | `src/bin/isogeny_walk.rs`, `src/cryptanalysis/isogeny_walk/{walk,record,store,million,traits,queue}.rs` | done; 1,000,000-curve P-256 grid replayed (`research/p256_isogeny_million_20261006/`) | class enumeration engine; `ClassInfo` already carries `D_π`, its trial factors, twist order, embedding degree, per-`ℓ` kind |
| Isogeny algorithm library: Vélu/√élu/Kohel, BMSS, `Φ_ℓ` by Hecke/CRT, SEA with isogeny cycles, Kohel End-conductor (`path/endo.rs`), class group & relation lattice (`path/relation.rs`), CSIDH, quaternion/KLPT/Brandt for supersingular, `F_{p²}`, `F_{p⁴}`, `GF(2ⁿ)`, `GF(3ⁿ)`; CLI prints one JSON per run with `PASS/FAIL/INDETERMINATE` | `isogeny_algos/` (`docs/USAGE.md`, `docs/SURVEY.md`) | done, 109 tests | point counting at any size, End conductor, class groups, extension-field arithmetic for the `F_{p^n}` walker |
| Prime-field index calculus (`S₃`, abscissa factor base, Gaussian elimination mod `n`) | `src/cryptanalysis/ec_index_calculus.rs`, `ec_index_calculus_j0.rs` | toy (< 40 bits) | `γ`, `S_IC` on `F_p` |
| Binary/Koblitz IC pipeline with oracles `enumerate`, `sat`, `groebner`, `f4-f2`, `crossbred-f2`, `wdsat`; subspace bases; `ic run --degree n --solver …` | `src/bin/ic.rs`, `docs/ic/README.md`, `src/cryptanalysis/{semaev_*,koblitz_groebner,sat,crossbred}.rs` | done to `n = 23` sweeps | `γ`, solver work, Macaulay stats on `F_{2^n}` |
| First-fall-degree harness (binary Semaev `S₃` after Weil descent; Macaulay rank at every degree vs generic Hilbert series) | `src/cryptanalysis/ffd_harness.rs` | done (binary only) | `d_ff` label |
| Refutation-degree harness (`D*` of a non-decomposable PDP instance, factor base a subspace `V`) | `src/cryptanalysis/pc_degree_harness.rs`, `pc_degree_avg.rs` | done (binary only) | `D*` label |
| `F_p` and tower F4, sparse Macaulay, Gaudry cubic (`F_{p³}` symmetrised `S₄`) | `f4_fp.rs`, `f4_fp_tower.rs`, `sparse_macaulay.rs`, `gaudry_cubic.rs` | done | building blocks for the missing `F_{p^n}` degree harness |
| `ecbench`: counted unit `S`, interleaved sessions, A/A, `ic_yield` view, `run_solver` table (`semi_regular_degree`, `solving_degree_max`, `macaulay_{rows,columns,degree,rank}`, SAT stats) | `src/bin/ecbench.rs`, `src/cryptanalysis/ecbench/`, `docs/ecbench/schema.sql` L390 | done; no algebraic oracle on prime curves | the label ledger of record |
| Yield sweep protocol (bases × oracles × solvers on identical targets) | `research/ecbench_yield_sweep_20261004/` | run at `n ≤ 23`, `p ≤ 2^20` | template for the labelling spec |
| Preregistered class-difficulty protocols I-1..I-5 | `research/isogeny_class_difficulty_20261008/` | frozen, not run | the controls this plan inherits; I-4 and I-5 are the labels' null |
| Large-prime-degree edges (rational-torsion window, class-group route) | `research/large_prime_isogeny_degree_20261008/` | protocol | vertical edges with `ℓ > 61` |
| Walker engineering: precomputed `Φ_ℓ` tables, BMSS kernels | `research/isogeny_walker_engineering_20261008/` | protocol K-1, K-2 pending | lifts the walker's `ℓ ≤ 61` cap and `ℓ²` step cost |
| Distributed execution: `taskq` (Redis/k8s, pinned commits), `isolab` (Linux hub/worker, NUMA, counters) | `taskq/`, `isolab/` | done | sharded generation and labelling |
| GPU kernels: batched `F_p` Macaulay reduction, Semaev pair sweep | `gpu/macaulay/`, `gpu/semaev/` | **never run on a GPU** | optional accelerators for labelling at the top tier |
| Curve identity: ICV1 slugs, typed curve links, registry | `docs/curve-identities.md`, `docs/curves/registry.json`, `docs/curves/ic/curve-links/` | done | every node gets an identity |

### 2.2 `ml-cryptanalysis` (Python; the learner)

| capability | where | state |
|:--|:--|:--|
| Toy-tier curve table: families `random/j0/j1728/cm_h1`, `p` 16–21 bits, rho `aut`/`plain` arms, 24 seeds | `csd/dataset.py`, `scripts/generate_dataset.py`, `data/curves.jsonl` (1,608 rows), `data/curves_holdout.jsonl` (300), `data/enumeration.jsonl` (6,032 isomorphism classes over three fields) | done |
| Invariants: `D_K, f, h_K`, certified `f_End` (division-polynomial Frobenius test, `MAX_TORSION = 64`, Hilbert `h=1` certificate), volcano height/level for `ℓ ∈ {2,3,5,7}`, exact rational `ℓ`-isogeny counts, observed `ℓ = 2, 3` neighbourhoods and crater lengths (Vélu, `CRATER_CAP = 60`), `Cl(End)` structure, `ord_ℓ`, 2/3-torsion, embedding degree, twist data | `csd/invariants.py`, `csd/isogeny.py`, `csd/endo.py`, `csd/classgroup.py`, `csd/volcano.py` | done; `O(p)` vectorised point counting caps `p` near `2^21` |
| Splits by isogeny class `(p, |t|)` and by field; label permutation within prime; leak audit | `csd/splits.py`, `csd/experiment.py` | done |
| Models: ridge, 2-layer MLP; exact-rule miner; symbolic predicate search with permutation band; autoencoder novelty; anomaly replication | `csd/models.py`, `csd/mine.py`, `csd/symbolic.py`, `csd/anomaly.py` | done |
| Results: planted `√|Aut|/2` law recovered; Kohel vertical rules and crater = `ord_ℓ` recovered exactly on 6,032-class enumeration; rho residual null; 0 replicated anomalies | `experiments/{planted-aut,isogeny-structure,anomaly-search}/report.md` | done |

---

## 3. Gaps (what is missing for the stated goal)

Numbered so the phases in §11 can cite them.

- **G1. No algebraic labels on the ML side.** `ml-cryptanalysis` measures
  rho only. `γ`, `d_ff`, `D*`, Macaulay/solver statistics and `S_IC` exist
  only inside `crypto` and only for binary fields (degree harnesses) or
  prime fields (yield, no degree). There is no join between the two.
- **G2. No degree-of-regularity instrument over `F_{p^n}`.** `ffd_harness`
  and `pc_degree_harness` are binary-Semaev-after-Weil-descent only. The
  natural algebraic PDP over odd characteristic is `E/F_{p^n}`, `n ∈ {2,3,4,5}`,
  with the factor base `x ∈ F_p` (Gaudry/Diem/Joux–Vitse) and the
  Weil-descended `S_{n+1}` system over `F_p`. `f4_fp_tower.rs` and
  `gaudry_cubic.rs` are the building blocks; the harness is absent.
- **G3. No prime-field degree label at all.** Over `F_p` with the abscissa
  base `{x < B}` the decomposition is root-finding plus a range condition,
  not a polynomial system; the honest labels there are `γ`, its product-law
  ceiling, `S_IC`, and the oracle constant. The plan records this rather
  than inventing a degree.
- **G4. The walker is `F_p` only.** `isogeny_walk` walks prime fields. Classes
  over `F_{p^n}` (needed for G2) and `F_{2^n}` (where the degree harnesses
  live; a `Φ_ℓ mod 2` neighbour path exists in `isogeny_algos::binary` and
  in `binary_isogeny`, per `PROTOCOL-I3`) have no class-enumeration driver
  with certified edges, identities and sharding.
- **G5. `ℓ ≤ 61` cap and `ℓ²` per-step cost** (`PROTOCOL-K1`); vertical
  edges at primes `ℓ | f_π` with `ℓ > 61` unreachable except through the
  rational-torsion window (C4) which is a protocol, not code.
- **G6. End conductor at scale.** `csd/endo.py` certifies `f_End` by division
  polynomials up to prime powers `≤ 64` at toy `p`. `isogeny_algos`
  `path/endo.rs` (Kohel) exists but is not exposed through the walker's
  node record; the walker records `D_π` and `End = Z[π]` only when `D_π` is
  fundamental. Per-node `f_End`, `D_End`, `Cl(End)` are needed for every
  node of every class.
- **G7. Torsion structure is thin.** Only `E(F_p)[2]` rank and `3`-torsion
  presence. Needed: the group structure `Z/n₁ × Z/n₂` of `E(F_q)`, the
  `ℓ`-primary parts, `E(F_{q^k})[ℓ]` for `k ≤ 6` via the trace recurrence,
  the Frobenius eigenvalue orders `r(λ), r(μ)` mod `ℓ` for `ℓ` up to a
  window bound (`PROTOCOL large_prime_isogeny_degree`, E1), and the twist's
  structure.
- **G8. Neighbourhood features are local and tiny.** `ℓ ∈ {2,3}` degree-1
  neighbourhoods and crater length. Needed: radius-`R` balls per `ℓ`,
  neighbour level multisets, ball growth, cycle counts, `Φ_ℓ` root
  multiplicities (ramified edges), spectral statistics of the crater
  component, and the typed multigraph itself as model input.
- **G9. Model-dependent features that IC actually reads are absent** from
  the vocabulary by design (they were "hidden" for the planted-law test):
  `a, b, j`, `|Aut|`, the abscissa distribution of the chosen factor base on
  this model, `c₄, c₆`, the isomorphism-class representative. For
  predicting `γ` they are *the* candidates; the vocabulary needs two
  modes (§5, `hidden` vs `model` block) and the leak audit needs to be
  per-label.
- **G10. Scale.** Python `O(p)` point counting and enumeration cap the
  learner at `p ≈ 2^21`; the Rust side has SEA and `Φ_ℓ` at 256 bits. There
  is no bulk path "Rust emits invariants → Python learns" and no storage
  format for `10^5`–`10^6` nodes with edges.
- **G11. No geometric model.** Ridge and MLP only; no graph model, no typed
  edges, no class pooling, no expressivity diagnostics.
- **G12. No ceilings beside the labels** in the ML table: the product-law
  yield ceiling, the Bardet–Faugère–Salvy semi-regular degree, the generic
  Macaulay rank, and the rho birthday bound must be columns so skill is
  measured against theory, not against the mean.
- **G13. No power analysis.** A null result without a minimum detectable
  effect is not a result.
- **G14. Operational.** `/Volumes/SSD990` is at 100 % (16 GB free); all
  generated data must go to `/Volumes/SSD990-2`. Conductor's control plane
  timed out at plan time (`conductor check` TLS handshake); the executing
  agent must re-run the claim. The GPU kernels have never run on a GPU and
  there is no NVIDIA device on the local M4 Pro (14 cores, 48 GB).

---

## 4. Labels (what is predicted), with producers and ceilings

Every label row carries: curve identity (ICV1 slug + model hash), field
type, `(ℓ_fb, m, family)` of the factor base, target index, seed, status
(`ok | timeout | exhausted | unverified | skipped`), producer binary and
commit, and the ceiling/baseline columns of §4.4. Unit for cost is the
`ecbench` group operation; solver work is in its own unit and marked
unpriced, exactly as `ecbench_yield_sweep` does.

### 4.1 Field types and which labels apply

| field | IC family | `γ` | `d_ff` | `D*` | Macaulay/F4/SAT stats | `S_IC` | rho `S` (control) |
|:--|:--|:--:|:--:|:--:|:--:|:--:|:--:|
| `F_p` | abscissa base `{x < B}` or `hash_to_subgroup` base; `S₃` (`m = 2`), `S₄` (`m = 3`, `semaev_higher`) | yes | **no** (G3) | no | no (root-finding oracle) | yes | yes |
| `F_{p^n}`, `n ∈ {2,3,4,5}` | base `x ∈ F_p` (Gaudry–Diem); Weil-descended `S_{n+1}` over `F_p`; symmetrised variant | yes | yes (new, G2) | yes (new) | yes (`f4_fp_tower`, `sparse_macaulay`) | yes | yes |
| `F_{2^n}`, `n` prime and composite | subspace base `V`, `dim ℓ_fb`; `S₃`/`S₄` after Weil descent | yes | yes (`ffd_harness`) | yes (`pc_degree_harness`) | yes (`ecbench run_solver`) | yes | yes |

### 4.2 Definitions

- **`γ` relation yield** = verified relations / trials over a fixed trial
  count on a fixed `(family, ℓ_fb, m)`; recorded with `trials`,
  `relations`, `lookups`, `lift_failures` (`ic_yield` columns). Also
  `γ_class` = `γ` per class-of-points after `Aut` folding (I-5 Q4).
- **`d_ff` first-fall degree** = smallest Macaulay degree `D` at which the
  rank deficit exceeds the generic deficit (`ffd_harness` definition),
  with the full rank-vs-degree profile `rank(D)`, `D = 2..D_max` kept.
- **`D*` refutation degree** = smallest degree at which `1` enters the
  Macaulay row space of a *non-decomposable* instance (`pc_degree_harness`),
  averaged over `k` non-decomposable targets (`pc_degree_avg`), with the
  per-target vector kept.
- **Semi-regular degree** `D_reg` (Bardet–Faugère–Salvy) from `(n_vars,
  n_eqs, degrees)` — analytic, not measured; recorded as a baseline column
  and, when the solver reports `semi_regular_degree`, both are kept.
- **Solving degree** `D_solve` = `solving_degree_max` from the F4/Gröbner
  run; **solver work** = `solver_ops` in `solver_op_unit`, `sat_conflicts`,
  `sat_decisions`, `macaulay_rank/rows/columns`.
- **`S_IC`** whole pipeline in group operations over `√r`, every phase
  charged, `Q = [d]P` verified on the original curve after transport back
  when the arm is a neighbour (I-4 instrument).
- **Rho `S`** = structure-blind control label (`rho.negation` on `F_p`,
  `rho.signed_frobenius_strong` on Koblitz); any model that predicts it from
  structure beyond `r` and `|Aut|` is wrong.

### 4.3 Labels at larger sizes

Above the labelling tier (§6.3) no label is measured; nodes carry features
only. They are used for (i) extrapolation tests of a model trained below,
scored only where a cheap label (`γ` at reduced trial count) can still be
produced on a sampled subset, and (ii) unsupervised structure checks. No
prediction on an unlabelled registered curve (P-256, secp256k1, ECC2K-130)
is reported as a hardness statement; it is a model output with its
calibration interval and the note that it is out of the training support.

### 4.4 Ceiling and baseline columns (G12)

| column | formula | source |
|:--|:--|:--|
| `gamma_ceiling` | product law for `m`-sums of a `|F|`-point set in a group of order `r` (`RESEARCH_ECC2K130_DECOMPOSITION.md` §5) | analytic |
| `D_reg_bfs` | Bardet–Faugère–Salvy semi-regular degree for the system shape | analytic |
| `rank_generic(D)` | `cols(D) − max(H_gen(D), 0)` (`ffd_harness` doc) | analytic |
| `S_rho_expected` | `√(π r / 2) · √(2/|Aut|)` normalised | analytic |
| `S_floor` | `√(π / 2A)` generic floor as `ecbench` records it | `ecbench` |

Skill is always reported as improvement over a model that sees only these
columns plus `(log q, log r, ℓ_fb, m)`.

---

## 5. Feature vocabulary

Three blocks by leak status. Every feature names its producer and whether
it is a class invariant (constant across the `F_q`-class) — the model is
also trained with class invariants removed, so H2 is tested cleanly.

### 5.1 `class` block (class invariants; constant across the class)

| feature | definition | producer |
|:--|:--|:--|
| `log2_q, char, ext_deg n` | field | walker/CLI |
| `t, trace_ratio, log2_N, log2_r, cofactor, N_factorization` | Frobenius trace, group order, largest prime, cofactor factored | SEA (`isogeny_algos count`) |
| `D_pi, D_K, f_pi, f_pi_factors, v_ell(f_pi)` for `ℓ ≤ L` | `t² − 4q = f_π² D_K`; conductor evaluations | `ClassInfo` + factoring (trial to `2^20`, then ECM/`factorint` via PARI where needed; `composite_unfactured` kept as status) |
| `h_K, Cl(O_K)` structure (`exponent, cyclic, 2-rank`) | class group of the maximal order | `csd/classgroup.py` (toy), `isogeny_algos path/relation.rs` (large) |
| `h(Z[pi])`, `Cl(Z[π])` structure, `ord_ℓ` of primes above `ℓ` in `Cl(Z[π])` | the Frobenius order's class group: crater lengths at the floor | same |
| `kron(D_K, ℓ)`, `kron(D_pi, ℓ)` | split/inert/ramified per `ℓ` | arithmetic |
| `n_isog_exact(ℓ)` at the floor and surface | `0, 1, 2, ℓ+1` from Frobenius eigenlines | `csd/isogeny.py`, walker `ells` |
| `embedding_degree` (cap 24), `twist_order`, `twist_r_bits`, `twist_factors` | pairing/twist | walker |
| `E(F_q)` structure `(n₁, n₂)` with `n₂ | n₁`, `ℓ`-primary ranks | group structure (same for all curves in the class? **no**: the group structure is *not* a class invariant; it is in 5.2) | — |
| `frob_eigen_order(ℓ)`: `r(λ), r(μ), r_min(ℓ)` for `ℓ ≤ L_eig` | multiplicative orders of Frobenius eigenvalues mod `ℓ`; `r_min(ℓ)` = least `k` with `ℓ | #E(F_{q^k})` by the `s_k` recurrence | new (`eigenvalue_orders.py` is a prototype) |
| `log2_#E(F_{q^k})`, `k ≤ 6`, with largest prime factor bits | extension orders | recurrence |
| `volcano_shape(ℓ)` | height `v_ℓ(f_π)`, crater length, whether `ℓ` splits/ramifies/is inert | arithmetic |
| `anomalous, MOV_small, smooth_cofactor` flags | audits | walker class audits |

### 5.2 `node` block (varies across the class, allowed for all labels)

| feature | definition | producer |
|:--|:--|:--|
| `f_End, D_End = f_End² D_K, log2_f_End, v_ℓ(f_End)` for `ℓ ≤ L` | endomorphism order; **certified** (`divpoly` to `MAX_TORSION`, Kohel `path/endo.rs`, Hilbert `h=1`), with `f_End_method` and `certified` flag | `csd/endo.py` (toy), `isogeny_algos` Kohel (all sizes), walker when `D_π` fundamental |
| `End_maximal`, `level(ℓ) = v_ℓ(f_End)`, `depth(ℓ) = v_ℓ(f_π) − v_ℓ(f_End)` | position in each `ℓ`-volcano | derived |
| `h_End, Cl(End)` structure, `ord_ℓ` in `Cl(End)` | class group of the node's own order | `csd/classgroup.py` / `relation.rs` |
| `E(F_q) ≅ Z/n₁ × Z/n₂`, `ℓ`-primary type `(a_ℓ, b_ℓ)` for `ℓ ≤ L`, `tors_rank(ℓ)` | rational torsion structure (changes along vertical edges) | Weil pairing / random points + division (`isogeny_algos`), `csd/divpoly.py` at toy size |
| `E(F_{q^k})[ℓ]` rank for `k ≤ 6`, `ℓ ≤ L` | extension torsion | eigenvalue orders + structure |
| `aut_order ∈ {2,4,6}` (char > 3), `{2,4,6,12,24}` in char 2, 3; `j ∈ {0,1728}`; `twist_class` (which twist of `j`) | automorphisms | model |
| `n_isog_observed(ℓ)`, `n_h(ℓ), n_u(ℓ), n_d(ℓ)` | degree-1 neighbourhood by direction, `ℓ ≤ L_walk` | walker edges |
| `phi_root_mult(ℓ)` | multiplicity pattern of `Φ_ℓ(X, j)` roots (double roots = ramified edge) | walker |
| `crater_len(ℓ)` observed, `on_crater(ℓ)` | horizontal cycle | walker |
| `ball_size(ℓ, R)`, `R ∈ {1,2,3}`; `ball_level_hist(ℓ, R)` | growth and level multiset of the `ℓ`-ball | walker |
| `cycle_count(ℓ, len ≤ 6)` through the node, in the typed multigraph | local cycle structure | graph pass |
| `spectral(ℓ)`: top-3 eigenvalues of the normalised adjacency of the connected component restricted to the node's level, `λ₂` gap | spectral summary | graph pass (networkx/scipy) |
| `dist_to_surface(ℓ)`, `dist_to_floor(ℓ)` | hops | BFS |
| `path_from_root` encoding (`r.d3.d3.h5`) and depth | walker route | walker |

### 5.3 `model` block (depends on the equation; the IC candidates; G9)

| feature | definition | why |
|:--|:--|:--|
| `a, b` reduced to `min(v, q − v)`, `coefficient_bits`, `a_minus_3_model`, `qr_prefix_64` | the walker's trait detectors | model structure IC can see |
| `j`, `j_bits`, `j_in_subfield(d)` for `d | n`, `coeffs_in_subfield(d)` | subfield-definedness (extension fields) | I-3 screens |
| `c4, c6`, discriminant class mod squares | isomorphism data | model |
| `fb_abscissa_stats`: for the registered factor base on this model, the count of base points, their `x`-distribution moments, fraction of base points with rational `2`-torsion image | the factor base *as realised on this model* | `γ` is a function of this set |
| `ghs_m_min`, `ghs_m_by_factorisation` (binary composite degree) | GHS magic number | positive control |
| `cover_status`, `jv_cover` genus | cover existence | I-3 |
| `subfield_coeffs` (binary) | | |

### 5.4 `hidden` block (never shown; used for the planted-law positive controls)

For each experiment the plan names what is hidden. Default hidden set for
`γ`-prediction: nothing (model block allowed). For the planted controls
(§8.4): `aut_order, j, a, b` hidden, exactly as `csd/invariants.HIDDEN`.

### 5.5 Feature production contract

- Rust emits exact invariants per node as one JSON object per line
  (`features.jsonl`), schema `ec-hardness/features/v1` (§6.6), with every
  unknown `null` and a `status` per block (`exact | certified | observed |
  uncertified | unsupported`). Python never recomputes an exact invariant
  it was given; it computes derived and graph features only.
- Every feature has a unit test against an independent computation at toy
  size (`csd` already does this for `f_End`, volcano rules, crater = `ord_ℓ`).

---

## 6. Data generation design

### 6.1 Field types and sizes (tiers)

| tier | field | size | classes | nodes per class | labels | purpose |
|:--|:--|:--|--:|--:|:--|:--|
| T0 | `F_p` | 16–24 bits | 300 | full class (`h(D_π)·Σ` levels; cap 5,000) | all `F_p` labels, rho control, 8 targets × 2,000 trials | main supervised set |
| T0b | `F_{p^n}`, `n ∈ {2,3}` | `p^n` 18–30 bits | 150 | full class (cap 5,000) | `γ, d_ff, D*, Macaulay, S_IC`, rho | the degree labels |
| T0c | `F_{2^n}`, `n ∈ {11,13,15,17,19,21}` incl. composite `{12,15,21}` | | all ordinary classes at `n ≤ 15` (complete), sampled at `17–21` | `γ, d_ff, D*`, solver stats, GHS `m` | binary degree labels + GHS positive control |
| T1 | `F_p` | 24–40 bits | 300 | sampled 2,000 per class (BFS from root + random crater jumps by `Cl` action) | `γ` (500 trials), `S_IC`, rho | size extrapolation |
| T1b | `F_{p^n}` | 30–44 bits | 100 | 1,000 per class | `γ, d_ff` (D_max capped), `D*` on a 100-node subsample | degree extrapolation |
| T2 | `F_p` | 40–64 bits | 100 | 1,000 per class | `γ` on 50-node subsample, otherwise features only | feature-only extrapolation |
| T3 | registered curves: P-256 million grid (exists), secp256k1, ECC2K-130 neighbours (`Φ_ℓ mod 2`) | | existing stores | none | unlabeled structure; out-of-support predictions reported as such |

Totals: roughly `3·10^5` labelled `F_p` nodes, `10^5` labelled `F_{p^n}`
nodes, `10^5` binary nodes, `5·10^5` feature-only nodes, plus the existing
`10^6` P-256 grid. The "million curves with full graphs" requirement is met
at T0–T1 by full/large class enumerations with every edge certified, and at
T3 by the existing grid.

### 6.2 Families per tier (controlled distribution; every node tagged)

1. `random`: uniform `(a, b)`, rejection to cofactor `≤ 4` (I-1 rule) and
   `r ≥ 2^{11}` bits at toy; one class per `(p, t)` kept.
2. `j0`, `j1728`: `D_K = −3, −4` classes (the `|Aut|` control); all sextic /
   quartic twists.
3. `cm_h1`: class-number-one `D` with ordinary reduction (no extra `Aut`).
4. `cm_small_h`: `D_K ∈ {−15, −20, −23, −24, −31, −35, −39, −40, −47, …}` via
   Hilbert class polynomials mod `p` (`isogeny_algos` has `Φ_ℓ`; Hilbert
   polynomials for `|D| ≤ 2000` to be generated by CM or read from a table).
5. `volcano_tall`: `p, t` chosen so `v_ℓ(f_π) ≥ 3` for `ℓ ∈ {2, 3, 5}`
   (`csd/volcano.py` does this at toy size); gives vertical depth.
6. `large_prime_conductor`: `ℓ | f_π` with `61 < ℓ < 2^{12}` (I-5 input);
   one vertical edge per class through the rational-torsion window.
7. `anomalous_near`, `low_embedding` (embedding degree `≤ 6`): known-weak
   controls, so the model sees class-level weakness and the splits can
   hold it out.
8. Binary: random ordinary `y² + xy = x³ + a x² + b` over `F_{2^n}`, Koblitz
   `a ∈ {0,1}` with `b = 1`, composite-`n` classes for GHS.
9. `F_{p^n}`: random `E/F_{p^n}` with `j ∉ F_p` and the subfield-curve
   control `j ∈ F_p` (positive control for subfield features).

Sampling is deterministic from one `seed` per tier; primes are drawn as
`csd/dataset.default_primes` does (balanced residues mod 12 so both
automorphism groups are possible), extended to the tier's bit range.

### 6.3 Class enumeration algorithm (per root curve)

1. **Count** `#E(F_q)` (SEA via `isogeny-algos count`, certified ICV1), derive
   `t, D_π, f_π, D_K`; factor `f_π` (trial to `2^{20}`, PARI `factorint` for
   the cofactor at T1+; status kept).
2. **Degrees to walk**: `L_walk = {ℓ prime ≤ 61 : (D_π/ℓ) ≠ −1 or ℓ | f_π}`
   (walker `ells`), plus every `ℓ | f_π` with `ℓ > 61` through the
   rational-torsion window (`r_min(ℓ) ≤ 6`) — new code (G5).
3. **BFS** from the root over `Φ_ℓ(X, j)` roots, each edge certified by its
   kernel polynomial (`IW1`), dedup by canonical model; record direction
   (`h/u/d`) by comparing `v_ℓ(f_End)` of the endpoints; stop at the class
   cap with the frontier recorded (`PROTOCOL-I3` budget rule).
4. **Crater jumps** at T1+: when the class exceeds the cap, sample nodes by
   acting with random `Cl(Z[π])` elements (ideal → kernel via
   `path/relation.rs` and `path/endo.rs`), so the sample is not a BFS ball
   around the root; tag `discovery_route ∈ {bfs, class_group, torsion_window, vertical}`.
5. **Per node**: certified `f_End` (Kohel), torsion structure, `Aut`, model
   traits, factor-base realisation stats; write `features.jsonl`.
6. **Edges**: `edges.jsonl` with `(src_slug, dst_slug, ℓ, direction, kernel_hash, phi_root_mult)`.
7. **Verify** (`isogeny_walk verify` semantics): replay identities, orders and
   every edge certificate from the stored routes; a class whose replay
   fails is dropped with its failure kept.

For `F_{p^n}` (G4): same algorithm with `isogeny_algos` `fp2.rs`/`ext.rs`
arithmetic and `Φ_ℓ mod p` evaluated at `j ∈ F_{p^n}`; the End conductor by
Kohel over `F_{p^n}`. For `F_{2^n}`: the `Φ_ℓ mod 2` neighbour path
(`isogeny_algos::binary`, `binary_isogeny`), char-2 Vélu/Kohel kernels for
certificates, complete enumeration at `n ≤ 15`.

### 6.4 Identity and storage

- Every node: ICV1 slug (`docs/curve-identities.md`), canonical model hash,
  `class_id = (field, t)`, `discovery_route`, `path_from_root`.
- Layout under `/Volumes/SSD990-2/ec-hardness/` (G14):

```
ec-hardness/
  manifests/<tier>/<class_id>.json      # config, seed, commit, hashes of every shard, status
  classes/<tier>/<class_id>/nodes.jsonl # one object per node: identity + 5.1/5.2/5.3 blocks
  classes/<tier>/<class_id>/edges.jsonl
  classes/<tier>/<class_id>/walk.json   # ClassInfo, counts, frontier
  labels/<tier>/<class_id>/<label>.jsonl# one object per (node, target, arm)
  ecbench/<tier>/sessions/…             # raw ecbench sessions + audits (unit of record)
  tables/<tier>.parquet                 # joined, one row per (node, label arm) — produced by Python
  graphs/<tier>/<class_id>.npz          # CSR typed adjacency for the GNN
```

- JSONL is the interchange; Parquet (`pyarrow`) the analysis format; SQLite
  via `ecbench` for the counted labels. Shards are hashed (SHA-256) into the
  manifest; `collect` refuses missing, duplicated or edited shards (reuse
  `isogeny_walk collect` semantics).
- Size: node record ≈ 2 KB, edge ≈ 150 B, label row ≈ 400 B. `10^6` nodes ≈
  2 GB; `5·10^6` edges ≈ 0.75 GB; `3·10^6` label rows ≈ 1.2 GB; ecbench
  sessions dominate (≈ 20–50 GB at T0 with per-run records). Budget 100 GB.

### 6.5 Compute budget (order of magnitude; measure in Phase 1 and replace)

| stage | per unit | units | total |
|:--|:--|--:|--:|
| SEA count at `≤ 64` bits | ms | `10^3` roots | minutes |
| `Φ_ℓ mod p` build, `ℓ ≤ 61` | once per `(ℓ, p)`: `≈ ℓ⁵` mults, ~1 s at `ℓ = 61` | 18 primes × 900 fields | ~5 core-h (eliminated by K-1 tables) |
| one certified edge, `ℓ ≤ 61`, `p ≤ 2^{64}` | 1–50 ms | `5·10^6` | ~20–70 core-h |
| Kohel `f_End` per node | 1–100 ms | `10^6` | ~10 core-h |
| `γ` label: 8 targets × 2,000 trials of `S₃` root-finding at `p ≤ 2^{24}` | ~5 ms/target | `3·10^5` nodes | ~4 core-h; T1 at 500 trials ≈ 2 core-h |
| `d_ff` + `D*` over `F_{p^n}` / `F_{2^n}`, Macaulay to `D ≤ 5`, `n ≤ 19` | 1–60 s | `2·10^5` | 100–3,000 core-h (the dominant cost; subsample `D*` to 100 nodes/class if over budget) |
| rho control, 24 seeds | ~50 ms | `3·10^5` | ~100 core-h |
| GNN training (CPU/MPS) | | | hours |

Local M4 Pro (14 cores): T0 and T0c fit in days. T0b/T1b degree labels and
T1/T2 enumeration go through `taskq`/`isolab` workers; confirm the fleet
before Phase 3 (G14). GPU Macaulay kernels are optional and gated on a
first-ever GPU run with the Python oracle check (`gpu/macaulay/README.md`).

### 6.6 Schemas (to be written as JSON Schema files in `docs/ec-hardness/schema/`)

- `features/v1`: `{slug, model_hash, class_id, tier, field:{char, n, modulus}, model:{a,b,…}, class:{…5.1}, node:{…5.2}, model_feats:{…5.3}, status:{block→status}, producer:{bin, commit, args}}`
- `edges/v1`: `{src, dst, ell, direction, kernel_hash, phi_root_mult, certificate_ref}`
- `labels/v1`: `{slug, label, family, ell_fb, m, target_index, seed, arm, value, aux:{trials, relations, rank_profile, …}, ceiling:{…4.4}, unit, status, producer, ecbench_record_id}`

---

## 7. Labelling pipeline

### 7.1 `F_p` (`γ`, `S_IC`, rho)

- Instrument: native `ic` pipeline with the abscissa base at the size's
  registered dimension (16/32/64 points at 16/18/20 bits in the yield sweep;
  register dimensions for 22–40 bits by the same rule) and
  `hash_to_subgroup_v1` base as a second family; oracle `subtract` and
  `mitm:negation_folded=1`; `S₃` and `S₄` (`m = 2, 3`).
- Per node: 8 hidden targets, 2,000 trials (T0) / 500 (T1), seed
  `20261020 + tier`; `γ` with `gamma_ceiling`; `S_IC` whole pipeline;
  matched `rho.negation` arm, A/A arm first, in one interleaved `ecbench`
  session per class (the I-4 layout: root, two small-`ℓ` neighbours, one
  large-`ℓ`/vertical neighbour, isomorphic-model control, then every
  remaining node as a `γ`-only arm).
- Required new `ecbench` capability: a `γ`-only arm type that records the
  relations phase and skips linear algebra, so `10^5` arms per tier are
  affordable; and an **isomorphic-model control** arm (`u⁴a, u⁶b`) on 5 % of
  nodes — a `γ` difference between isomorphic models of the same curve is
  the floor of what counts as "model-dependent".

### 7.2 `F_{p^n}` (`γ`, `d_ff`, `D*`, Macaulay, F4, `S_IC`) — new instrument (G2)

- `src/cryptanalysis/ffd_harness_fp.rs`: build `S_{n+1}` on `E/F_{p^n}`,
  restrict `x_i ∈ F_p`, Weil-descend to `n` equations in `n·(n)` … variables
  over `F_p` (follow `gaudry_cubic.rs` for `n = 3`), symmetrise by `S_{n}`
  action (`symmetrized_semaev.rs`), build the Macaulay matrix at
  `D = 2..D_max` with `sparse_macaulay.rs`, rank over `F_p`, report
  `rank(D)`, `rank_generic(D)`, `d_ff`.
- `pc_degree_harness_fp.rs`: the same on non-decomposable targets (verified
  non-decomposable by `enumerate` at toy size), `D*` per target, `k = 8`.
- F4 run via `f4_fp_tower.rs` with `solving_degree_max`, op counts into the
  `run_solver` block; the `ecbench` `descent-algebraic` arm must accept a
  prime-power field. `ecbench` presently runs no algebraic oracle on prime
  curves (`ecbench_yield_sweep/PROTOCOL.md`); extend `methods.rs`.
- Cross-check: `enumerate` oracle agrees with F4 on every decomposition
  verdict (I-4 Q3) on T0b.

### 7.3 `F_{2^n}` (existing)

- `ic run --degree n --solver {enumerate, sat, f4-f2, crossbred-f2}` on the
  class's nodes with subspace base dimension `ℓ_fb` fixed per `n`
  (`ecbench_yield_sweep` values: 6 and 8), `m = 2, 3`; `ffd_harness` and
  `pc_degree_harness` per node; GHS `m(b)` from `audit_curve`.
- Reference `rho.signed_frobenius_strong` on Koblitz, `rho.negation` else.

### 7.4 Rules inherited from I-4/I-5

Never change `ℓ_fb`, `m`, family or trial count between arms of a class;
never report `γ` without its ceiling; keep timeouts and failures as rows;
verify `Q = [d]P` on the original curve for every transported arm; `sat`
and `enumerate` must agree.

---

## 8. Dataset assembly, splits, leak audit, controls

### 8.1 Join (Python, `csd/assemble.py`)

Read `nodes.jsonl`, `edges.jsonl`, `labels/*.jsonl`; compute derived and
graph features (5.2 ball/spectral/cycle features with `networkx`/`scipy`);
write `tables/<tier>.parquet` with one row per `(node, label, arm)` and
`graphs/<class_id>.npz`. Assert every exact invariant that `csd` can
recompute at toy size agrees with the Rust value (`tests/test_bridge.py`).

### 8.2 Splits

- **Isogeny-class split**: key `(field, |t|)`; twists stay together.
- **Field-family split**: hold out whole fields; the largest size is test.
- **Size-extrapolation split**: train `≤ 24` bits, test `24–40`.
- **Route split**: train on `bfs` nodes, test on `class_group`-sampled nodes
  of the same classes (does the model generalise away from the root ball?).
- **Vertical split**: train on surface nodes, test on floor nodes and vice
  versa (does `v_ℓ(f_End)` carry anything?).

### 8.3 Leak audit, per label

For each label, assert the hidden block is absent; fit with and without each
block (`class`, `node`, `model`, `graph`); report the drop. For `γ`, assert
that no feature is a deterministic function of the label's own trial
outcomes (e.g. `fb_abscissa_stats` is computed from the base, not from the
relations).

### 8.4 Nulls and positive controls (every round)

- Label permutation within `(field, class)` for regressors and inside the
  symbolic search.
- Structure-blind control label: rho `S`.
- **Planted controls** (must be recovered before any blind result is read):
  (a) the `√|Aut|/2` rho law with `aut_order, j, a, b` hidden (exists);
  (b) `γ` with `|Aut|` folding on `j ∈ {0,1728}` nodes (`γ_class` differs
  from `γ` by exactly `|Aut|/2`; hide `aut_order`); (c) GHS: on composite
  binary classes the node with minimal `m(b)` has a descent-aligned `γ`
  outside the A/A interval (I-4 Q4); (d) `F_{p^n}` subfield control: `j ∈
  F_p` nodes have a lower `d_ff` than `j ∉ F_p` nodes of the same class
  size, if and only if they do in the measurement (pre-registered as a
  *check*, not a prediction).
- Oracle ceiling: when a label has a known exact law (planted), report
  recovered fraction of the law's own R².

---

## 9. Models

### 9.1 Baselines (must be beaten; reported always)

1. Ceiling-only: predict from `(log q, log r, ℓ_fb, m, |Aut|)` + §4.4 columns.
2. Ridge on the full tabular vocabulary (`csd/models.fit_ridge`).
3. Gradient-boosted trees (LightGBM or sklearn `HistGradientBoosting`) with
   monotone constraints on size columns; SHAP for attribution.
4. MLP (`csd/models.fit_mlp`).
5. Exact-rule miner and symbolic predicate search (`csd/mine.py`,
   `csd/symbolic.py`) — readable conjectures with permutation band.

### 9.2 Geometric model

- Input: the class's typed multigraph. Node features = 5.2 + 5.3 (+ 5.1
  broadcast); edge type = `(ℓ, direction)`; edges undirected for `h`,
  directed for `u/d`.
- Architecture: relational message passing (R-GCN / typed GINE) with one
  weight set per `(ℓ, direction)`, 3–4 layers (radius matches 5.2 balls),
  residual, LayerNorm; node head for node labels; class pooling (sum + max)
  for class labels; multi-task heads for `γ`, `d_ff`, `D*`, `S_IC`, with
  the rho head as the control (its loss must stay at the size baseline).
- Framework: PyTorch Geometric on CPU/MPS; graphs from `graphs/*.npz`.
- Uncertainty: deep ensembles (5 seeds) + split-conformal intervals on the
  held-out classes; calibration plots per tier.
- Training: class-level batching, early stopping on held-out classes,
  weight decay, seeds fixed; log every run (`experiments/hardness/runs/*.json`).

### 9.3 Expressivity diagnostic (pre-registered)

The horizontal `ℓ`-graph at a fixed level is a Cayley graph of `Cl(O)` (I-1
reasoning; crater = `ord_ℓ` is a unit test). Cayley graphs are vertex
transitive, so a message-passing network with *structural features only*
assigns the same embedding to every node of a crater; any node-level skill
must come from node features. Test: train the GNN with node features
removed (constant input) — its node-level skill must be zero on horizontal
labels; train with `model` block only versus `node` block only — this
separates H1 from H2. Report the 1-WL colour-refinement partition sizes per
class as a column.

### 9.4 Interpretation

SHAP on the GBDT; integrated gradients on the GNN; symbolic regression on
the top-k features; every discovered rule re-checked for zero violations on
the complete binary enumerations and the toy `F_p` enumeration
(`data/enumeration.jsonl` style), as `csd/mine.py` does.

---

## 10. Evaluation protocol and decision rules (registered)

### 10.1 Metrics

R² and MAE in the label's natural scale (`log γ`, integer degrees treated as
ordinal with exact-match rate, `log S`), per tier, per split, with bootstrap
95 % intervals over classes (not rows). Skill = metric minus the
ceiling-only baseline's metric.

### 10.2 Pass/fail

- **Q1 (controls).** All planted controls recovered at ≥ 80 % of their
  oracle R²; the rho control head shows no skill above baseline (interval
  contains 0). Fail → fix the pipeline; read nothing else.
- **Q2 (H0).** For each label, skill over the ceiling-only baseline on the
  isogeny-class split and on the field-family split, BH-adjusted across
  labels × splits. Interval contains 0 on both → **H0 holds at this tier**:
  boundary; the label is class-invariant up to the A/A spread and the
  model is a size model.
- **Q3 (H1).** Skill > 0 on both splits, survives the leak audit, the
  isomorphic-model control (the effect must be absent between isomorphic
  models), and the fresh-seed replication at an unseen field size with the
  same sign → **reproducible model-dependent effect**; name the features
  carrying it; this earns its own `isogeny_class_difficulty` protocol at the
  next size. It is not an advance: no exponent moved.
- **Q4 (H2).** GNN skill minus tabular skill with identical node features,
  BH-adjusted; the structural-only GNN must be at zero (§9.3). Positive and
  surviving → neighbourhood information matters; otherwise the
  neighbourhood is redundant given the node's invariants.
- **Q5 (extrapolation).** Size-extrapolation split skill reported; no
  claim at a size the training never covered beyond the interval.

### 10.3 Inadmissible moves

Changing `ℓ_fb`/`m`/family between arms; reading a solver-stage number as
a yield; a hardness statement about P-256, secp256k1 or ECC2K-130 from a
model output; dropping failed/timed-out rows; promoting any model
prediction to a ledger without a frozen paired `ecbench` challenge.

### 10.4 Power (G13)

Before T0 labelling, simulate the detectable effect: with `C` classes of
`k` nodes and per-node `γ` standard error `σ/√trials`, compute the smallest
class-internal ratio detectable at 80 % power under the class split; fix
`trials` and `k` so that a 5 % yield effect is detectable at T0 and 10 % at
T1. Record the numbers in the dataset card; report "nothing here" relative
to them.

---

## 11. Phased work plan (for the executing agent)

Each phase lists tasks, files, acceptance gate, estimate. Commit per phase
on `feat/elliptic-curve-complexity-ml-4b4228` (crypto) and a branch
`hardness-predictor` (ml-cryptanalysis); open PRs ready for review (not
draft); update titles as scope grows (`AGENTS.md`).

### Phase 0 — contracts and skeleton (1–2 days)

- `docs/ec-hardness/schema/{features,edges,labels}.v1.schema.json`; a
  `schema_check` test in both repos.
- `ml-cryptanalysis`: `csd/bridge.py` (read Rust JSONL → feature dicts),
  `csd/assemble.py` (join + parquet + graph npz), `csd/graphfeat.py` (balls,
  cycles, spectral), `csd/ceilings.py` (§4.4 formulas with unit tests
  against hand values), `csd/power.py` (§10.4), `pyproject` adds `pyarrow,
  scipy, networkx, scikit-learn, torch-geometric` as optional extras.
- Re-run `pytest` in `ml-cryptanalysis` (baseline green) and
  `cargo test --release --lib` in `crypto` before touching anything.
- Gate: schema tests pass; `assemble` round-trips the existing
  `data/curves.jsonl` into parquet with identical values.

### Phase 1 — Rust feature producer (3–5 days)

- `src/bin/ec_hardness_features.rs` (or `isogeny_walk features`): per node of
  a stored walk, emit `features/v1` using: `ClassInfo` (class block), Kohel
  End conductor from `isogeny_algos` (add a dependency path or a thin FFI/CLI
  call), torsion structure `(n₁, n₂)` and `ℓ`-primary parts, eigenvalue
  orders `r(λ), r(μ)` for `ℓ ≤ 1000`, `#E(F_{q^k})` for `k ≤ 6`, `Aut`,
  trait detectors, `Φ_ℓ` root multiplicities, factor-base realisation stats
  for the registered bases. (G6, G7, G9)
- `isogeny_walk walk`: record `direction` on every edge (requires `f_End`
  of both endpoints) and `discovery_route`. (G8)
- Tests: at toy `p` every exact invariant equals `csd`'s value on the 1,608
  existing curves (`tests/test_bridge.py` drives the binary).
- Gate: 100 % agreement on exact invariants; `null` only where `status`
  says uncertified.

### Phase 2 — class generation at T0 and T0c, local (3–5 days, mostly wall time)

- `scripts/gen_classes.py` (ml side, orchestration only) or a Rust
  `ec_hardness gen` subcommand: sample primes/families (§6.2), count, walk
  full classes (cap 5,000), verify, write manifests to
  `/Volumes/SSD990-2/ec-hardness/`.
- Binary classes at `n ≤ 15` complete via the `Φ_ℓ mod 2` path; Koblitz and
  composite-`n` roots registered by slug.
- Gate: `verify` replays every class; node/edge counts match
  `h(D_π)`-based expectations where the class is complete (Kohel structure
  theorem counts); manifests hashed.

### Phase 3 — labelling at T0/T0c and the power run (5–10 days, wall time dominated)

- `ecbench` additions: `γ`-only arm type; isomorphic-model control arm;
  session generator `scripts/spec_from_class.py` producing one spec per
  class in the I-4 layout. (§7.1)
- Run T0 `F_p` `γ`/`S_IC`/rho locally (14 cores); run T0c binary `γ`,
  `d_ff`, `D*` locally with `D_max = 5` and 100-node `D*` subsample per
  class when the budget is exceeded.
- `csd/power.py` run on the first 30 classes; fix `trials`, `k`.
- Gate: every session audited (`ecbench` audit JSON), A/A interval
  recorded, I-4 Q2 (ceiling) and Q3 (oracle agreement) pass on every arm.

### Phase 4 — `F_{p^n}` instruments and T0b (5–8 days)

- `ffd_harness_fp.rs`, `pc_degree_harness_fp.rs`, F4 arm over prime-power
  fields in `ecbench::methods`; `F_{p^n}` class walker (`isogeny_algos`
  `fp2/ext` arithmetic, Kohel over `F_{p^n}`). (G2, G4)
- Tests: at `n = 2`, `p ≈ 2^{10}`, `d_ff` and `D*` agree with a brute-force
  Macaulay rank in Python (`sympy`/`numpy` mod `p`) on 20 instances; the
  subfield control (`j ∈ F_p`) behaves as measured, not as assumed.
- Generate and label T0b.
- Gate: cross-check passes; label rows carry `D_reg_bfs` and
  `rank_generic(D)`.

### Phase 5 — learning and reports at T0 (3–5 days)

- `scripts/run_hardness.py`: assemble → splits → baselines → GBDT → MLP →
  GNN (PyG) → nulls → planted controls → §9.3 diagnostic → BH → report
  `experiments/hardness/report.md` with every table of §10, the dataset
  card and the model card.
- Gate: Q1 passes before Q2–Q5 are read; report states H0/H1/H2 verdicts
  per label per tier with intervals and the minimum detectable effect.

### Phase 6 — scale: T1/T1b/T2 on the fleet, K-1 tables, optional GPU (2–4 weeks wall)

- `Φ_ℓ` tables for `ℓ ≤ 199` (K-1) and BMSS kernels (K-2) to lift the
  walker cap; rational-torsion window for vertical `ℓ > 61` (C4). (G5)
- `taskq`/`isolab` plans: one spec per class, pinned commit, shards
  collected with the refusal rules; results to `SSD990-2`.
- Size-extrapolation and route splits; refit; re-report. Out-of-support
  predictions on T3 reported only with the §4.3 caveat.
- GPU Macaulay kernels only after the first GPU run passes the Python
  oracle check in `gpu/macaulay/README.md`.

---

## 12. Risks and mitigations

| risk | mitigation |
|:--|:--|
| H0 holds everywhere (most likely on `F_p`) and the project "finds nothing" | the deliverable is the bounded null with minimum detectable effect (§10.4) and the planted controls proving the pipeline could have seen an effect; this closes I-4/I-5 at toy size with data |
| degree labels dominate compute | subsample `D*`; cap `D_max`; keep rank profile to `D_max` and mark truncated; GBDT on `d_ff` first |
| class caps bias toward BFS balls around the root | `class_group` sampling route and the route split (§8.2) |
| leakage through class invariants (every node in a class shares `r`) | class split is mandatory; skill reported over the ceiling baseline that already knows `r` |
| `f_End` uncertified at some nodes | `status` column; models trained with and without uncertified rows |
| walker cap `ℓ ≤ 61` hides vertical large-`ℓ` structure | Phase 6 K-1/C4; until then `large_prime_conductor` family is limited and labelled as such |
| disk | `SSD990-2` only; manifests hashed; raw ecbench sessions compressible |
| Conductor unreachable | re-run `conductor check` before each phase; report scope; do not edit outside `research/ec_hardness_predictor_20261010/`, `docs/ec-hardness/`, `src/bin/ec_hardness_*`, `src/cryptanalysis/{ffd_harness_fp,pc_degree_harness_fp}.rs`, `src/cryptanalysis/isogeny_walk/`, `src/cryptanalysis/ecbench/methods.rs` without reporting |

---

## 13. Reporting (per `AGENTS.md` and `.agents/skills/report-evidence`)

Each phase ends with a status table: requirement → `verified complete |
implemented but unverified | partial | blocked | not attempted`, with the
exact command, commit, inputs, environment, wall time and the gap. Failed
runs, timeouts and regressions stay in the record. Interpretation is
labelled separately from measurement. No "breakthrough" language; no
statement about deployed curves.
