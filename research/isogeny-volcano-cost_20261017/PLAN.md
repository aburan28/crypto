# Do isogenous curves differ in ECDLP cost? Experiment plan

Status: **draft plan (v0), not frozen.** `PROTOCOL.md` will freeze it after
the tooling in §8 exists and the dry run in §9 has passed. No result below
is a finding. The pre-flight numbers in §7 only size the design.

## 1. Question and the claims each outcome licenses

On curves linked by ℓ-isogenies (one volcano, same `F_p`, same `#E`), does
the cost of solving one ECDLP instance depend on the curve? That covers
rho, BSGS and index calculus, separately for the volcano **level**
(crater … floor) and for the individual **curve**. Each candidate cause is
measured on its own.

| outcome | claim it licenses (fixed wording) |
|:--|:--|
| level contrast's 90% CI inside ±M | "level changes *method* cost by less than ±M on these classes and sizes" |
| level contrast significant after correction, positive control passed | "level changes *method* cost by d% [CI]", with its cause named in §3 |
| curve heterogeneity significant, model heterogeneity not | "isogenous curves differ by τ% [CI] beyond coordinate effects" |
| model heterogeneity also significant | "the effect is a coordinate (model) effect, not an isogeny effect" |
| positive control fails | "inconclusive: the instrument could not see a known effect" |

Not in scope: any security statement, and anything above the sizes measured.
"No difference" is never written; only an equivalence bound or a measured
effect is.

## 2. What is already known (do not re-ask)

- `runs/run-20261017T0415Z-instance` (`research/isogeny-volcano-rho`): the
  generator and k do not change rho's cost (5 arms × 1,000 solves; every
  |z| < 0.8).
- `runs/run-20261017T0111Z`: under the pull-back walk, walks are identical
  on every curve.
- Curves A–E: plain rho shows no level effect. The j = 0 crater fold gives
  1.89 ± 0.09× and j = 1728 gives 1.32 ± 0.06× fewer operations.
- Earlier crater/floor studies (`research/*conductor_gap*`,
  cryptanalysis `experiments/volcano-ic`): every same-method ratio
  replicated at 1. IC is 0.998 [0.970, 1.026]. The only real effects are
  automorphism or Frobenius folds.

  Weaknesses this plan fixes:
  - **Few curves:** 1–6 fixed pairs, with targets counted as if they were
    curves.
  - **Wide intervals:** ±30–150% at 99%, so no equivalence was shown.
  - **Controls:** fixed walk seeds, an A/A arm that duplicated its
    partner, and no isomorphic-model control for IC.
  - **Untested costs:** transport never priced.

## 3. Sources of variation and how each is isolated

| # | source | design element |
|:-:|:--|:--|
| 1 | rho walk randomness | many paired targets per curve; distributions compared, not single solves; the pull-back walk as a zero-variance cross-check on a subset |
| 2 | instance (P, k) | carried by isogeny: transported generator `G_c = φ(G_root)`, `planted_shared` k, so `Q_c = φ(Q_root)` on every curve |
| 3 | curve model (a, b within the isomorphism class) | **M = 4 isomorphic models of every curve** (`(u⁴a, u⁶b)`, generator mapped). Model is a random effect nested in curve. The rho hash and the IC small-abscissa base both depend on it. |
| 4 | group structure | record each curve's rational ℓ-torsion (`kernels`) and `E(F_p)` structure; test cost against structure as well as level |
| 5 | automorphisms (j = 0, 1728) | those classes are the **positive control**, analysed apart from the pooled contrast |
| 6 | per-operation arithmetic | operation counts (ecbench GAE/S) are the primary metric; wall time per operation only at isolation L2 (§6, family F5) |
| 7 | implementation artifacts | one binary, one thread, the same parameters on every curve; restarts and cycle escapes logged; `ecbench verify --replay-all` |
| 8 | host noise | counted units are host-independent; timing runs at L2 on an isolated runner with shuffled order and an A/A arm |
| 9 | transport cost | Vélu-chain operations per carried point counted and reported as a fraction of solve cost |
| 10 | size | primary 24-bit `p`; secondary 20 and 28 bits (rho also at 32 and 40) |

Every effect is located in a nested model,
`cost ~ level + (1 | class / curve / model)`, with target and walk seed as
the residual.

## 4. Population

- **Classes (K = 8).** Six ordinary volcanoes plus two positive-control
  classes, all from `volcano-rho --bits B --depth 3 --seed s`:
  - ℓ = 2: 3 classes (ramified, split and inert crater as available);
  - ℓ = 3: 3 classes;
  - one j = 0 class (`--j0`) and one j = 1728 class (`--j1728`).

  Seeds are fixed in the protocol, and every curve the walk finds is used
  (whole volcanoes, about 30 curves per ℓ = 2 class).
- **Null population.** Isomorphic models of the same curves. They have the
  same group, generator and instance, and differ only in coordinates, which
  is exactly what an isogeny effect must beat.
- **Sizes.** 24-bit `p` for every family. 20 and 28 bits for F1 and F3 on
  2 classes, for a size trend. Rho alone at 32 and 40 bits on 2 classes.

## 5. Arms (one ecbench session per class and size; same spec template)

| arm | method | role |
|:--|:--|:--|
| `rho-neg` | `rho.negation` | reference (matched rho) |
| `rho-neg-aa` | `rho.negation` with a different `measurement.seed` | A/A noise floor (independent seeds this time) |
| `bsgs-neg` | `bsgs.negation` | exact-invariance check |
| `ic-mitm` | `ic.pipeline` `prime-abscissa:size=32`, `mitm:negation_folded=1` | IC primary |
| `ic-mitm-33` | the same with `size=33` | **planted IC positive control** (yield should rise by ≈ (33/32)² = +6.3%) |
| `ic-isogeny` | `ic.pipeline` `prime-isogeny:t=…` | a model-free factor base, to separate a curve effect from a model effect |
| `rho-cubic` / `rho-quartic` | j = 0 / j = 1728 classes only | positive control (√3 and √2) |

Target law: `planted_shared` with the transported generator via
`prime_explicit`. The same `target_seed` on every curve gives paired
targets.

## 6. Metrics, test families and decision rules

`m = 4` corrected tests (F1 and F3 × {level contrast, curve heterogeneity}),
α = 0.05, Bonferroni α/m = 0.0125, equivalence margin **M = ±5%** (90% CI).

- **F1 rho operations (S).**
  - Statistic: the floor/crater ratio of mean S, pooled over classes.
    Inference uses a cluster bootstrap over classes, then curves, then
    targets.
  - Level test: a permutation of level labels within class.
  - Heterogeneity: Cochran Q over curve means against its permutation
    distribution, compared with the same Q over models.
- **F2 BSGS operations.** Deterministic. Prediction: for each target, S is
  identical on every curve and model. One mismatch falsifies the transport
  or the instrument. No statistics.
- **F3 IC relation yield** `Y = Σ relations / Σ trials`, per curve and per
  level. The same F1 tests are applied to Y. IC end-to-end S per curve is
  reported with the AGENTS.md §5a band: unpriced counters at 0, 0.1 and 1
  GAE on every arm, the rho reference included. The gap to the matched rho
  is reported as distance remaining.
- **F4 group structure (secondary, exploratory).** Cost against
  rational-torsion class within each level. This is reported only as a
  hypothesis for a later protocol.
- **F5 time per operation (only at L2).** ns per GAE per curve, with the
  A/A arm as noise floor. A difference counts only if it lies outside the
  A/A band. **Not run on this VM, which is L0** (§10).

**Power check (skill §3).** Permutation tests use 10,000 permutations, so
the smallest attainable p is 2/10,001 = 2·10⁻⁴ < α/m = 0.0125, and
`N + 1 ≥ 4m/α = 320` holds. Pass.

**Sample sizes from the pre-flight variances (§7):**
- **F1:** sd/mean of S is 0.46. A ±2.5% half-width on a level ratio (90%)
  needs about 2·(1.645·0.46/0.025)² ≈ 1,830 solves per level, pooled.
  - Plan: 64 targets × all curves, plus a crater top-up session (the
    crater has 1–2 curves) to 512 solves per crater curve.
  - That gives ≥ 3,000 solves per level over the six ordinary classes.
- **F3:** one `ic-mitm` run finds 20 relations. ±2.5% on Y per level (95%)
  needs about (1.96/0.025)² ≈ 6,150 relations per level per class, i.e.
  ≈ 310 runs per level.
  - Per curve: ≥ 1,000 relations (50 runs), ±6.2% per curve.
  - The planted +6.3% control must come out with its CI excluding 0.
- **Models:** M = 4 per curve on F1 and F3 (16 targets each per extra
  model).

## 7. Pre-flight (scratch session `ECBS1h003daf885073`, L0, `preflight/`)

One 24-bit ℓ = 2 depth-3 volcano: 30 curves plus 2 isomorphic controls,
2 paired targets each, 512 runs, all verified, 1 min 42 s.

| arm | S mean | S sd | wall ms/run (mean, max) |
|:--|--:|--:|:--|
| rho-neg | 1.165 | 0.541 | 10.8, 14.9 |
| bsgs-neg | 0.664 | 0.097 | 9.3, 14.6 |
| ic-mitm | 61.2 | 22.2 | 52.1, 109.1 |
| ic-rootfind | 3,823 | 1,404 | 455, 988 (dropped: 9× the cost of mitm, same question) |

- **BSGS:** for each target, all 32 curves gave the identical S (1
  distinct value per target). The transport and pairing work as designed.
- **IC yield:** 0.95–2.39 ×10⁻⁴ per curve. That 2.5× range comes from only
  40 relations per curve (Poisson ±31%), which is why F3 needs about 50
  runs per curve.

**Runtime estimate.**
- **Per 24-bit class:** 30 curves × 4 models × (64 rho + 16 BSGS + 50 IC
  runs, plus the IC positive control) ≈ 120 × (0.7 + 0.15 + 2.6 + 2.6) s ≈
  12 min. Add the crater top-up, about 1 min.
- **All 24-bit classes:** 8 classes ≈ 1.8 h single-core.
- **Size arms:** 28-bit IC costs about 16× per run, ≈ 3 h for 2 classes at
  1 model; 20-bit is negligible.
- **Rho at 32 and 40 bits:** ≈ 1 h.
- **Total:** ≈ 6 h single-core at L0, about 2 h at 4-way parallel sessions.
  Fixed costs (volcano search and walk, generator transport) are seconds.

## 8. Tooling to build before freezing

1. `volcano-rho --ecbench-spec FILE [--models M] [--targets T]`: write the
   `prime_explicit` spec for a walked volcano directly, with transported
   generators and M models per curve. Also write a `curves.jsonl` sidecar
   with level, kernels, group structure, model u and parent edge, so the
   spec is generated natively, not by Python.
2. `volcano-rho` records each curve's `E(F_p)` structure (n₁, n₂) and its
   rational 2- and 3-torsion.
3. A native analysis binary (`volcano-stats`, Rust): it reads ecbench
   `records.jsonl` and the sidecar and computes the F1–F3 statistics
   (cluster bootstrap, within-class permutations, Cochran Q, TOST) into
   `analysis.json`. Python only draws (`plot_*.py`, the volcano heat maps).
4. Register the class root curves in `docs/curves/registry.json` before
   citing them in prose (AGENTS.md §11).

## 9. Dry run

Run the full pipeline on one class at N/8 into scratch, with every family,
the positive controls and the analysis binary. It must show:
- BSGS exactly equal across curves and models;
- the `size=33` IC control detected with the right sign;
- the j = 0 `rho-cubic` control detected;
- the analysis reproducing `ecbench compare` on a two-arm subset.

Then freeze `PROTOCOL.md` with its hash, seeds and N, and disclose
anything the dry run showed about the classes.

## 10. Host and timing

This container is a shared 4-vCPU VM, so ecbench runs at **L0**.

- **Counted units** (S, relations, trials) are host-independent and fully
  admissible from here.
- **Wall time and ns/op (F5)** need L2: an isolab runner
  (`ecbench-independent-runner`), which also provides the independent
  replay a claim needs.
- Until then, F5 is reported as descriptive only, or not at all.

## 11. Does not count

- Changing N, the margin, the metric, the class seeds or the arms after
  seeing results.
- Dropping a curve, a model or a failed run (`incomplete` rows stay).
- Selecting a "fast" curve and reporting it without re-measuring on fresh
  targets and fresh classes.
- Pooling sizes into a scaling claim.

Any such change is a new, labelled protocol version, and both versions are
reported.

## 12. Deliverables

Sessions, audits with `--replay-all` receipts, `analysis.json`, a README
with the §5 report table in the skill's fixed wording, the volcano heat
maps per class (operations and yield panels), forest plots of the level
contrasts against the ±5% margin, a PDF, and a `benchmarks/` row.
