# Factor-base designs for prime and binary curves: register, Lean proofs, experiments

**Date:** 2026-10-15. **Branch:** `research/factor-base-designs-lean`.
**Kind of work:** a design round (24 factor-base strategies), the Lean
proofs of every statement in it that can be formalized, and the first
preregistered experiments (protocol frozen in [`PREREGISTRATION.md`](PREREGISTRATION.md);
rows in `runs/`, results in §6). Nothing here prices a pipeline or
reports `S`, and nothing here is a speedup.

**The rule this round follows.** A design is only *closed* when a Lean theorem
disposes of it. A design whose closing argument cannot be stated as a checked
theorem stays **open**, with an experiment. The first draft of this round
"killed on paper" eight designs. Under this rule, three are closed by Lean
theorems (K7, K8, K5). Two are closed in the mean only (K1 at set level, K2),
with an open remainder. Three are reopened (K3 rests on an unformalized
literature theorem, and K4 and K6 have no checkable statement). K1's paper
argument was wrong: the combination it dismissed had never been tested.

## 1. The objective: degree slope, not degree

The repository's yield conservation (`AR:KN-FIND-007`) makes the yield per
trial a function of `|V|` and `m` only. So "more yield" means a bigger base,
and a design is described by one number: how fast the solving degree grows
with base entropy, `σ = ΔD / Δ log₂|V|` per summand.

**Break-even.** The cost model per relation is
`cost_ω(ℓ, D) = C(mℓ, D)^ω · 2^{n − mℓ}`, at `m = 4`, `n = 131`.
`BreakEven.lean` checks, in the kernel, that:
- at slope 1 (the measured `m = 4` fit, `RESEARCH_DREG_MEASUREMENT.md`) and
  at slope 1/2, every step up in `ℓ ∈ [8, 32]` raises the cost, for both
  `ω = 2` and `ω = 3`;
- at slope 0, every step lowers the cost;
- at slope 1/4, a step lowers the cost at `ω = 2` and raises it at `ω = 3`.

So under the measured slopes the optimum is the smallest base: there is no
middle ground. An interior optimum needs `σ ≲ 1/4`, and that is the
preregistered **candidate** threshold. These are statements about the model,
not about real systems.

## 2. Design register

Status vocabulary:
- **open (E…)**: live, with an experiment;
- **measured**: preregistered arm run in this PR, with the outcome in §6;
- **closed (Lean)**: a Lean theorem disposes of it, with any bridge named;
- **closed-in-mean (Lean)**: the theorem covers the mean over targets, and
  the remainder is open;
- **reopened**: the earlier paper kill did not survive.

### Prime fields: new membership presentations

| # | design | status | Lean | experiment |
|:--|:--|:--|:--|:--|
| 1 | Exponent-interval bases `{h·g^k : k < 2^t}` via bit products; no smoothness of `p ± 1` needed | **measured: no** (§5: D saturated at 8; widest matrices, first to hit memory) | `exp_interval_identity`, `bitsVal_lt`, `iterate_sq` (the presentation cuts out exactly that set) | A1z, A1g |
| 2 | Torus exponent intervals via a Lucas ladder | open | `lucas_double`, `lucas_double_succ` (ladder steps) | E2 |
| 3 | Elliptic exponent intervals `{x(kQ′)}` on an auxiliary curve | open | — (x-only ladder correctness not formalized) | E3 |
| 4 | Parallel (sumset) towers `x = u + v` | **measured: no** (§5: D = A0 + 1 at equal `t`; wider matrices) | — | A4 |
| 5 | Radix-3 towers on auxiliary curves | open | `three_maximizes_bits_per_degree` (`q³ < 3^q` for all `q ≥ 2`, `q ≠ 3`); `p…_minus/plus` (field-side odd parts are tiny on the NIST primes, so the 3-power torsion must come from `E′`) | E5 |
| 6 | P-521 torus tower | open | `p521_plus` (`p + 1 = 2^521`) | E6 |
| 7 | Fan-in roots | open | — | E7 |
| 8 | Lopsided depths per summand | open | — | E8 |
| 9 | Weighted-homogeneous towers | open | — | E9 |
| 10 | Split `S₃` chain + towers | open | — | E10 |

### Changing the equation

| # | design | status | Lean | experiment |
|:--|:--|:--|:--|:--|
| 11 | Norm-divisibility chain over a Kummer fibre (prime) | **measured: open, replication owed** (§5: at `m = 3`, D 4 vs 6 and ≈ 4× narrower matrices, same growth per bit) | `chain_sound` (soundness at each root), `pencil_norm_factor` (the `m = 2` norm factorization), `norm_chain_excess` (one more equation than unknown, every `m`) | A11 |
| 12 | Binary norm chain as an auxiliary-variable solver input (credits `research/nagao_relations/` open item 1) | open | `chain_sound` applies verbatim | E12 |
| 13 | Algebraic large prime (partial divisibility) | open | — | E13 |

### Learned and searched

| # | design | status | Lean | experiment |
|:--|:--|:--|:--|:--|
| 14 | Learned coordinates in `L(3O)`, `L(4O)` (ml-cryptanalysis CSD) | open | — | E14 |
| 15 | Family-slope selection (score rules by `σ`) | open | `BreakEven` (the threshold) | E15 |
| 16 | Arity portfolio | open (engineering) | — | E16 |

### The eight earlier "paper kills", re-examined

| # | design | status now | what Lean proves | what stays open |
|:--|:--|:--|:--|:--|
| K1 | Gauss-period window in `s = x + 1/x` (type-II ONB) | **reopened** | `window_eq`: in characteristic 2, `span(c₁…c_ℓ) = span(c₁, c₁², …, c₁^ℓ)`, so the *set* is a scaled `C1` window. `dickson_mul`: `cᵢcⱼ = c_{i+j} + c_{i−j}`. `typeII_*`: the bases exist at `n = 11, 23, 113, 131, 233, 239` | The paper kill also claimed the `s`-coordinate base is already in the repository. It is not: `koblitz_symmetrised.rs` uses the `u`-frame with Frobenius-stable `V`. A window **in `s`** is untested: E-K1 |
| K2 | Frobenius-rotated slots | **closed-in-mean** | `patterns_neutral`, `sum_count`: every pattern over a domain of the same size has the same total yield, so `k` patterns cost `k` solves for `k` times the mean yield | A single fixed target; a merged solve cheaper than `k` separate ones |
| K3 | Optimizing index sets in `ℤ/(2n+1)` | **reopened** | — | Its kill rests on the linear Vosper theorem (Bachoc–Serra–Zémor 2017), which is not formalized here. Formalizing it, or measuring, is owed |
| K4 | MinRank / big-field low-rank modelling | **reopened** | — | No checkable statement disposes of it: E-K4 |
| K5 | Restricted intermediates in `S₃` chains | **closed (Lean)** | `restrict_le`: restricting the domain never increases the yield at any target | — |
| K6 | Shared-trunk towers | **reopened** | — | The paper kill was a cost model, not a theorem: E-K6 |
| K7 | Multiplicative towers in `F_{2^n}^*` | **closed (Lean, with bridge)** | `k7_power_step_trivial`: for every `(n, B_n)` in {(113, 3391), (131, 263), (163, 150287), (233, 1399), (239, 479), (283, 9623), (409, 10⁶), (571, 5711)}, every power step of degree `2 ≤ d < B_n` has fibres of size ≤ 1. `bounds_attained`: `B_n` divides `2^n − 1` (except 409, where `10⁶` is a search bound) | Bridge: `F_{2^n}^*` is cyclic of order `2^n − 1` (Mathlib has it; not imported, see §4) |
| K8 | Doubling/halving fibres on the target curve | **closed (Lean, with bridge)** | `doubling_fibre_subsingleton`: in `ℤ/r` with `r` odd, every fibre of `[2^t]` has at most one point | Bridge: a group of prime order is cyclic (Mathlib has it; not imported) |

## 3. Proof map

All declarations are in `lean/FactorBaseDesigns/`, namespace
`FactorBaseDesigns`. `#print axioms` output is in
[`proof_logs/axioms.txt`](proof_logs/axioms.txt). Every theorem uses only
the standard axioms `propext`, `Classical.choice` and `Quot.sound`, or no
axioms at all. There is no `sorry`, no custom axiom and no `native_decide`.
Every `decide` runs in the kernel (`decide +kernel`).

| statement | declaration | file | status |
|:--|:--|:--|:--|
| `k·a ≡ k·b (mod N)`, `gcd(k, N) = 1`, `a, b < N` ⇒ `a = b` | `mul_mod_injective`, `fibre_subsingleton` | `CyclicMul.lean` | proved |
| K8 in the `ℤ/r` model | `doubling_fibre_subsingleton` | `CyclicMul.lean` | proved; bridge to `E(F_p)` named |
| gcd certificate ⇒ every `d` in range is coprime to `N` | `coprime_of_prod`, `prodR_eq`, `dvd_pr` | `NumberFacts.lean` | proved |
| K7 for 8 binary degrees | `k7_coprime`, `k7_power_step_trivial`, `cert*`, `bounds_attained` | `NumberFacts.lean` | proved; bridge to `F_{2^n}^*` named |
| small parts of `p ± 1` (cofactor coprime to all primes ≤ 31) | `p256_minus` … `p521_minus`, `p521_plus` | `NumberFacts.lean` | proved |
| type-II ONB conditions | `typeII_11`, `_23`, `_113`, `_131`, `_233`, `_239` | `NumberFacts.lean` | proved |
| degree 3 maximizes bits per unit degree | `three_maximizes_bits_per_degree` | `NumberFacts.lean` | proved for all `q` |
| Dickson product identity | `dickson_mul` | `Dickson.lean` | proved (any commutative ring) |
| K1 set level: Gauss-period window = power window | `window_eq` | `Dickson.lean` | proved (characteristic 2; additive span) |
| total yield = domain size; K2; K5 | `sum_count`, `patterns_neutral`, `restrict_le` | `Yield.lean` | proved |
| #1 presentation correct | `exp_interval_identity`, `bitsVal_lt`, `iterate_sq` | `Encodings.lean` | proved |
| #2 ladder steps | `lucas_double`, `lucas_double_succ` | `Encodings.lean` | proved |
| #11 soundness, `m = 2` factorization, counting | `chain_sound`, `pencil_norm_factor`, `norm_chain_excess` | `Encodings.lean` | proved; completeness checked per instance by the harness |
| break-even instances of the cost model | `BreakEven.slope_*`, `binom_sanity` | `BreakEven.lean` | proved (instances of the model) |

## 4. Toolchain and reproduction

- Lean `v4.34.1` (`lean/lean-toolchain`), with **no dependencies** (Lean
  core only, as in `research/ic_tree_split_structured_v_20261011/lean/`).
  Mathlib was tried first. Its cache did not fit the container's disk
  allowance beside the toolchain, so the two bridges above, which Mathlib
  would close, are named rather than imported. Porting them is the one
  formalization step owed.
- `cd lean && lake build` (about 100 s, dominated by the kernel gcd
  certificates). Then `lake env lean AxiomAudit.lean`, and
  `lake env leanchecker --fresh FactorBaseDesigns`. The outputs are committed
  in `proof_logs/`.
- CI: `.github/workflows/factor-base-designs-lean.yml` runs all three and
  fails on any `sorry` or `native_decide` in the sources.

## 5. Results (p₁ = 786433, run-20261015T1657Z)

Rows: `runs/*.jsonl`. Analysis: `runs/analysis.txt`. Run report:
[`reports/run-20261015T1657Z/REPORT.pdf`](reports/run-20261015T1657Z/REPORT.pdf).
Review: [`reports/review-20261015T2030Z/REVIEW.pdf`](reports/review-20261015T2030Z/REVIEW.pdf).
Harness: `examples/fb_designs_slope.rs`. 140 planted systems: 116 completed,
20 stopped at the 300 s budget and 4 aborted for memory. Every planted
solution checked.

- **Instrument check passed.** A0 rebuilt with the pilot's seeds matches all
  20 committed pilot systems on D and on every structural field
  (`runs/instrument-check.txt`).
- **Preregistered verdicts.** *Candidate*: every series except A0, N and A4
  at `m = 3`, which read *no middle ground* (σ = 0.5). The four points below
  bear on what those labels mean. They change none of them.
  1. At `m = 2` the incumbent and the null both read *candidate*
     (σ = 0.19), so the threshold, which was derived at `m = 4`, does not
     discriminate there.
  2. N equals A0 in every cell despite different curves, sets and targets.
     At these sizes D reflects a presentation's shape, not the curve.
  3. A1z and A1g sit at D = 8, their input degree, while their matrices and
     times grow fastest. σ = 0 is saturation, not cheapness.
  4. A11 at `m = 3` is the one arm below the incumbent on both D (4 vs 6) and
     matrix width (≈ 4× narrower, ≈ 9× faster at `t = 4`). Its growth per bit
     is the same as A0's, though, so it is a constant-factor shift, not a
     flatter slope. It rests on 3 cells, one prime, and no null of its own
     shape.
- **Predictions (§3).** A1z = A0 ± 1 failed (+2 to +4). A1g = A1z held
  (dropping smoothness costs nothing in D). Neither A4 prediction held
  (observed: A0 at equal `t`, plus 1).
- **Owed before any further statement** (each by dated amendment): §5
  replication of A11 at `m = 3` on p₂; an A11-shape null; the `t = 5` cell
  completed; an `m = 4` series.

## 6. What this round does not establish

- No design here is shown to reduce any cost. The Lean theorems show that
  presentations are correct, that some designs are empty or neutral in the
  mean, and what the break-even means *within the cost model*. The question
  that matters, `σ` on real systems, is the experiments' job.
- "Never seen before" means not in `docs/ic/factor-base-encyclopedia/` or
  `research/nagao_relations/`, plus a literature check from memory. The
  closest prior work: Petit–Kosters–Messeng 2016 (towers), Nagao and
  Joux–Vitse (Riemann–Roch decompositions), FPPR 2012, Petit–Quisquater
  2012, Huang–Kosters–Yeo 2015, Bachoc–Serra–Zémor 2017. A novelty audit is
  owed before any design is described as new outside this repository.
