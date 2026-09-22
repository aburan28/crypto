# mq-fes Semaev oracle in the full index-calculus pipeline (2026-09-22)

Question: how much of `ic run --solver mq-fes` (Koblitz index calculus with
the quadratic `m = 2` Semaev oracle on Gray FES) can be removed without
changing what the oracle answers, measured end to end, cold, from curve
construction to a verified logarithm?

## Boundaries (existing, not chosen for this round)

- **Floor.** The free-oracle floor of
  `research/notes/ecc2k130/RESEARCH_WDSAT_IC_UNIFICATION.md`: targets and
  linear algebra charged, the decomposition decision free. No oracle
  change moves it, so every step here is at best **engineering**.
- **Reference (same decomposition problem, same instances).** The
  pipeline's two non-algebraic `m = 2` oracles, run as arms of the same
  matched round: `enumerate` (walk the factor base, test `R − P`) and
  `pair-table` (meet in the middle over a precomputed pair-sum table).
- **Reference (same DLP).** Pollard rho on the same subgroup. Its
  operation count is not converted into this pipeline's units, so `S`
  and the rho ratio are **null** for every row (see "Accounting").

## What changed (`src/cryptanalysis/mq_fes_semaev.rs`)

| step | commit | class | what moved |
|:--|:--|:--|:--|
| oracle profile (points walked, build/walk/lift ns) | `0236d06` | accounting | counters only; this binary is the baseline |
| packed template + swap-symmetric walk + lift on the fly | `248440d` | engineering | build 2.1 ms → 0.05 ms per call; `2^{2ℓ}` → `≈2^{2ℓ−1}` points |
| linear split (bilinear in `(a, s)`) | `89e07ac` | engineering | oracle `2·4^ℓ` → `2^ℓ(ℓ²/2 + 5ℓ/2 + 1)` word ops |
| 8 lanes, branch-free, AVX2 dispatch | `13223c1` | engineering | 490 µs → 74 µs per system at `ℓ = 11` |

The linear split is the substantive step. With `x₁ = Σaᵢvᵢ`, `x₂ = Σbᵢvᵢ`
and `s = a ⊕ b`, the two-summand Semaev polynomial becomes
`s²x_R² + (a² + as)x_R + a⁴ + a²s² + b`. Squaring is `F₂`-linear, so
the only quadratic Boolean monomials are `aᵢsⱼ`. For each `a` in Gray
order the system is `n` linear equations in `ℓ` unknowns, updated by XOR
deltas and eliminated. Unordered roots are visited once (`a ≤ b`). This
puts the algebraic oracle in the same `2^ℓ` class as `enumerate`.
`tests::linear_split_matches_symmetric_walk_on_semaev_systems` checks the
bilinearity on every curve used here.

## Matched full-DLP round (`results/round_02`)

Arms: the baseline and the two candidates are the binaries built at the
revisions above; the references use the candidate binary with
`--solver enumerate` / `--solver pair-table`. Everything else is
identical: legacy factor base index 0, `m = 2`, relation batch 1, target
seeds 1–8 (train) and 1001–1008 (independent holdout), three
repetitions interleaved across arms. 960 runs, every one `complete` and
verified (`[k]G = Q` against the planted `k`, zero inconsistent
relations, zero verification failures), none dropped.

Unit: cold end-to-end wall seconds of the whole `ic run` (curve,
factor base, template, every trial, lifting, filtering, linear algebra,
verification). Ratios are baseline ÷ arm: the paired geometric mean
over seeds of per-seed medians, with a 95% t interval.

| instance | runs verified | baseline s (8 seeds) | symmetric | **linear** | ref enumerate | ref pair-table | linear vs pair-table | trials base/linear | oracle µs/call base → sym → lin |
|:--|:-:|--:|--:|--:|--:|--:|--:|--:|--:|
| n=13 a=0 holdout | 24/24 | 0.117 | 1.55 [1.44,1.67] | 1.60 [1.47,1.74] | 1.80 [1.65,1.97] | 0.04 [0.04,0.05] | 36.6 | 92/92 | 562 → 117 → 93 |
| n=13 a=0 train | 24/24 | 0.137 | 1.76 [1.64,1.88] | 1.82 [1.70,1.96] | 2.08 [1.92,2.26] | 0.05 [0.05,0.06] | 35.7 | 126/126 | 564 → 99 → 75 |
| n=17 a=1 holdout | 24/24 | 0.124 | 4.40 [3.41,5.68] | 4.53 [3.40,6.04] | 0.75 [0.68,0.82] | 5.32 [3.97,7.15] | **0.85** | 158/172 | 699 → 85 → 76 |
| n=17 a=1 train | 24/24 | 0.136 | 4.74 [3.73,6.03] | 4.89 [3.73,6.41] | 0.79 [0.71,0.89] | 5.70 [4.28,7.57] | **0.86** | 171/181 | 713 → 83 → 72 |
| n=23 a=0 holdout | 24/24 | 1.931 | 4.49 [3.96,5.09] | **13.05 [11.51,14.81]** | 0.27 [0.26,0.29] | 3.06 [2.50,3.74] | 4.27 | 548/564 | 3388 → 617 → 141 |
| n=23 a=0 train | 24/24 | 2.326 | 5.04 [4.66,5.44] | **14.51 [13.30,15.83]** | 0.27 [0.26,0.28] | 3.70 [3.29,4.17] | 3.92 | 662/671 | 3393 → 608 → 135 |
| n=23 a=1 holdout | 24/24 | 2.105 | 4.97 [4.72,5.23] | **13.95 [11.99,16.24]** | 0.26 [0.23,0.29] | 2.95 [2.28,3.83] | 4.73 | 572/589 | 3554 → 620 → 141 |
| n=23 a=1 train | 24/24 | 2.127 | 5.33 [4.97,5.72] | **14.76 [13.41,16.24]** | 0.30 [0.28,0.34] | 3.22 [2.68,3.88] | 4.58 | 581/554 | 3535 → 599 → 143 |

Reading it:

- At `n = 23` the linear split is 13–15× faster than the baseline and
  3.9–4.7× faster than `pair-table`, the fastest pre-existing oracle,
  with every interval clear of 1 on holdouts as well as training seeds.
- At `n = 17` it is **slower than `pair-table`** (0.85×): the table's
  cold build is cheap there and its lookups are faster than a `2^8` split.
- At `n = 13` (`ℓ = 12`, 13 equations, many roots per target) every
  oracle stops early and process overhead dominates; `enumerate` is
  fastest and `pair-table` pays a 2.6 s table build.
- The oracles decide the same targets but may return a different
  decomposition, so the relation matrix fills differently. Trial counts
  move by −5% to +9%, and that cost is inside the end-to-end ratios.

## Oracle scaling, stage diagnostic (`results/oracle_scaling_01`)

`examples/mq_fes_oracle_scaling.rs`: both oracles on the same 12
targets per curve, every curve with `n ≥ 2ℓ` for `ℓ ∈ {3, 4, 5, 8, 11,
12, 14}` (overdetermined, so almost every target is refuted and pays the
full enumeration). The oracles agree on every target. The linear split's
one-off template (42 µs to 37 ms) is timed separately. This is a stage
diagnostic of the decomposition oracle, not a speed of the method.

| arm | fitted slope of log₂(word ops/call) on ℓ | prediction | fitted wall slope |
|:--|--:|:--|--:|
| baseline | 2.000 | `2·4^ℓ` → 2 | 1.03 |
| linear split | 1.316; **1.014** with its polynomial divided out | `2^ℓ(ℓ²/2 + 5ℓ/2 + 1)` → 1 | 0.83 |

Wall slopes sit below the word-op slopes because both oracles also pay a
polynomial per-call build (the baseline rebuilds its Weil system
symbolically every call), which dominates until `ℓ ≈ 11`. Per-call wall
ratio at `ℓ ≥ 11`: 31× (`n = 23`), 85× and 115× (`ℓ = 12`, `n = 39, 45`),
141× (`ℓ = 14`, `n = 43`).

`results/ell_sweep_01` is an earlier sweep that mixed in `n < 2ℓ` curves.
There early exits on underdetermined systems confound the fit (baseline
slope 1.6). It is kept as evidence, not used.

## Accounting

- The frozen WDSat regression suite
  (`research/index_calculus_baseline_20260914/regression`) drives an
  external `wdsat_solver` with the WDSat command/output protocol and
  conflict counter. This oracle is in-process Gray FES and linear
  algebra with no conflict counter, so that suite cannot run it. This
  round instead uses the matched full-DLP protocol above, under the
  parent contract's rules: identical inputs, independent holdouts,
  retained raw runs, verified answers, paired intervals.
- **Speedup in operations: not established.** `word_ops` covers only the
  oracle's walk/solve stage (two per Gray point; each delta and masked
  reduction XOR of the split). Template build, lifting, filtering and
  the linear algebra have no conversion into that unit. `S`, the rho
  ratio and the floor ratio stay null. The end-to-end figures above are
  paired wall time on one shared 4-vCPU host: a runtime claim with its
  interval, secondary by the repository's rule.
- **Class of the round: engineering.** Oracle cost fell, the free-oracle
  floor did not move, and the question each target is asked is
  unchanged. The linear split changes the oracle's exponent in `ℓ`, not
  the method's position against the floor.

## Falsification target for the next round (stated now)

Keep the linear split as the default `mq-fes` oracle only while, on a
fresh matched round at `n = 23`, holdout seeds give linear vs
`pair-table` ≥ 1.5 with the 95% interval above 1, every run verified,
and zero disagreements with `mq_fes_decompose_reference` over at least
160 targets per curve. Abandon the claim that it is the fastest `m = 2`
oracle here if any `n ≥ 23` instance puts `pair-table` ahead.
Inadmissible: changing the factor base, dropping template build from the
cold cost, favourable seeds, or comparing warm against cold.

## Reproduce

```bash
python3 run.py --arm baseline=BIN@0236d06:mq-fes:0236d06 \
  --arm mqfes-symmetric=BIN@248440d:mq-fes:248440d \
  --arm mqfes-linear=BIN@13223c1:mq-fes:13223c1 \
  --arm ref-enumerate=BIN@13223c1:enumerate:13223c1 \
  --arm ref-pair-table=BIN@13223c1:pair-table:13223c1 --output results/round_NN
python3 summarize.py results/round_NN
cargo run --release --example mq_fes_oracle_scaling > results/oracle_scaling_NN/raw.jsonl
python3 fit_scaling.py results/oracle_scaling_NN
```

`BIN@rev` is `target/release/ic` built at that revision; `manifest.json`
in each round records the binary hashes.
