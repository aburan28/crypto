# Unifying Koblitz index calculus with Trimoska WDSat

**Modules:** `src/cryptanalysis/wdsat_oracle.rs`,
`DecompositionStrategy::Wdsat` in `koblitz_index_calculus.rs`,
`ic run --solver wdsat --wdsat-binary PATH`
**Related:** [`RESEARCH_TRIMOSKA_BENCHMARKS.md`](RESEARCH_TRIMOSKA_BENCHMARKS.md),
[`RESEARCH_SAT_SEMAEV.md`](../index-calculus/RESEARCH_SAT_SEMAEV.md),
[`RESEARCH_ECC2K130_DECOMPOSITION.md`](RESEARCH_ECC2K130_DECOMPOSITION.md),
[`research/index_calculus_baseline_20260914/`](../../index_calculus_baseline_20260914/)

## Boundary (stated before measuring)

- **Reference.** Pollard rho on the same Koblitz subgroup, automorphism
  discount `√(2n)`, same operation accounting as the rest of the
  repository. On ECC2K-130 that is `2^60.8090` (`S ≈ 0.0774`).
- **Floor.** Free-oracle floor for an `m`-summand subspace factor base:
  charge targets and linear algebra, charge nothing for the
  yes/no decomposition decision. Formula as in
  `RESEARCH_ECC2K130_DECOMPOSITION.md` §5.3. The floor moves only with
  `B` and `m`; a better SAT engine cannot cross it.

**Falsification target.** This thread is a *unification* thread, not an
attack-improvement thread. It succeeds if (i) the Semaev / Weil system
used by native CDCL is emitted in Trimoska ANF, (ii) an external WDSat
binary returns a model that lifts to a verified point decomposition on
every planted prime-degree toy below, and (iii) the lifted witness
agrees with native SAT after sorting. It is abandoned if ANF emission
and WDSat disagree on any planted instance after independent group
verification, or if the adapter silently treats a capacity warning as a
refutation.

Inadmissible: claiming an end-to-end ECDLP speedup, moving the free-oracle
floor, or treating the ONB / Hamming-weight CryptoMiniSat path in
`ic fixed` as interchangeable with this Semaev ANF path without a matched
encoding.

## What was unified

Two previously separate stacks now share one decomposition API:

| Stack | Encoding | Solver | Factor base |
|---|---|---|---|
| Native IC (`--solver sat`) | Semaev Weil restriction → XOR/CNF | in-crate CDCL | Frobenius-invariant subspace |
| Trimoska / WDSat regression | Semaev Weil restriction → ANF | external `wdsat_solver` | same |
| **This PR (`--solver wdsat`)** | **same Semaev system → Trimoska ANF** | **external WDSat** | **same** |

The ECC2K-130 field `F_2^{131}` is a **prime-degree** extension, so it
sits in the same family as the Trimoska corpus (`n ∈ {15,17,19,…}`). The
unified oracle is the Semaev subspace question on that family. The
durable `ic fixed` workflow for ECC2K-130 still uses the ONB circuit +
CryptoMiniSat path for Hamming-weight bases; that is a different
encoding and is not relabelled here.

## Measured agreement (engineering)

Planted two-summand decomposition on `K_1 / F_2^7` (prime degree),
identical subspace factor base, independent group lift:

| variant | class | planted lift | agrees with native SAT |
|---|---|---|---|
| native `--solver sat` | reference | yes | — |
| `--solver wdsat` (WDSat `61c6ff3`) | **engineering** | yes | yes |
| `--solver mq-fes` (ALMASTY Möbius + Gray early-exit) | **engineering** | yes | yes |

Ratio to the free-oracle floor is unchanged: every oracle answers the
same algebraic question. Class is **engineering** by §3 of `AGENTS.md`.

The `mq-fes` backend ports the ALMASTY
[mq](https://gitlab.lip6.fr/almasty/mq) solvers (public domain):

| backend | when used | measured vs Möbius / baseline |
|---|---|---|
| Incremental Gray (libfes FFS, `L=4` + `Fl[0..4]` register block) | `find_one` / Semaev lift (early exit); default full enum | **~337×** wall faster than full Möbius on planted early root; **~37–50×** vs prior O(n)/step Gray on full `n=18` enum |
| `u32` L=4 twin (`m≤32`) | opt-in | ties or loses to `u64` at fit sizes (~0.9–1.0×); tables already L1-resident |
| Parallel outer specialisation (rayon, 4 outer bits) | opt-in (`gray_ffs_parallel_outer`) | within noise of serial `L=4` at `n=20` on a 4-core host (~0.96–1.34×); specialisation tax dominates below that |
| Möbius transform (`moebius.c`) | `find_all` for `n ≤ 24` | reference for all-roots |
| Monica hybrid (`monica.c`) | `n > 24` (range extension) | does **not** beat Möbius inside `n ≤ 24` (release wall on `n=14,m=32` was ~0.22×); calibrated cost model agrees |
| AVX2 Gray (`avx2_8x32` ideas → 4×u64 lanes, ± batch) | opt-in only (`mq_fes_avx2`) | **does not** beat packed-u64 scalar `L=4` on single-system Semaev (~0.42× per-step; ~0.6–0.8× batch on `n=16,m=24`); correct vs Möbius/L=4 |
| Hardcoded scalar `L=8` / `L=4` batch probe | opt-in / experimental | within noise of or **slower** than plain `L=4` at the sizes that fit (I-cache; well-predicted zero-checks beat rewind) |

Falsification for the “faster than Möbius” claim: a release run of
`gray_early_exit_beats_moebius_find_one_wall` must keep ratio `≥ 1.5` on
the fixed dense quadratic with a Gray-index-2000 planted root. Falsification
for the Gray speedup itself: `gray_ffs_beats_on_step_full_enum_wall` must
keep FFS/`L=4` ≥ 1.5× the prior O(n)-per-step Gray on full `n=18` enum.
Parallel outer must agree with serial (`parallel_outer_agrees_with_serial`);
auto-select stays off while `n=20` walls oscillate around 1×. AVX2 stays
opt-in. Monica is kept as a capacity extension, not as an in-cap speedup.
Cubic chained (`m ≥ 3`) Semaev systems are refused; those stay on SAT /
WDSat.

### mq-fes in the full pipeline: packed, bilinear oracle (2026-09-22)

Evidence, raw runs and reproduction: [`research/mq_fes_ic_pipeline_20260922/`](../../mq_fes_ic_pipeline_20260922/README.md).
Boundaries are the ones above (free-oracle floor; rho), plus the
pipeline's own `m = 2` oracles `enumerate` and `pair-table` as measured
references on the same instances.

Pricing the oracle inside `ic run` showed the earlier Gray work had
optimised the smaller cost. At `n = 23` a call spent 2.1 ms rebuilding
the Weil-restricted `S₃` system symbolically and 1.3 ms walking all
`2^{2ℓ}` assignments. The new default oracle
(`src/cryptanalysis/mq_fes_semaev.rs`) changes three things:

- **Packed template.** The system is affine in the target's bits, so it
  is stored once per factor base as bit-sliced words, and a call costs
  at most `n` table XORs.
- **Bilinear split.** With `s = a ⊕ b`, the only quadratic monomials are
  `aᵢsⱼ`. So each Gray-enumerated `a` leaves a linear system in `s`,
  which is eliminated in eight branch-free lanes (AVX2 when available).
  The walk drops from `2·4^ℓ` word ops to `2^ℓ(ℓ²/2 + 5ℓ/2 + 1)`.
- **Lift on the fly.** Roots are lifted as they are reached, with no
  64-root cap.

The previous oracle is kept as `mq_fes_decompose_reference` and
cross-checked on every target of a sweep.

| variant (matched full DLP, cold) | n=17 holdout | n=23 a=0 holdout | n=23 a=1 holdout | class |
|:--|--:|--:|--:|:--|
| mq-fes reference (baseline) | 1.00 | 1.00 | 1.00 | reference |
| + packed template, swap-symmetric walk | 4.40 [3.41, 5.68] | 4.49 [3.96, 5.09] | 4.97 [4.72, 5.23] | engineering |
| + bilinear split, 8 lanes | 4.53 [3.40, 6.04] | 13.05 [11.51, 14.81] | 13.95 [11.99, 16.24] | engineering |
| ref `enumerate` | 0.75 [0.68, 0.82] | 0.27 [0.26, 0.29] | 0.26 [0.23, 0.29] | reference |
| ref `pair-table` | 5.32 [3.97, 7.15] | 3.06 [2.50, 3.74] | 2.95 [2.28, 3.83] | reference |

Cells are baseline ÷ variant, cold end-to-end wall time, as a paired
geometric mean with a 95% interval over 8 holdout seeds (training seeds
agree; `n = 13` is in the round README). All 960 runs verified. `S`, the
rho ratio and the floor ratio are null: only the oracle stage has an
operation counter, and the other phases have no conversion into it.

The frozen WDSat suite cannot run an in-process oracle, as the README
explains. The stage fit over seven `n ≥ 2ℓ` sizes gives word-op slopes
of 2.000 for the baseline and 1.014 for the split once its polynomial
factor is divided out, matching the predictions.

At `n = 23` the split beats `pair-table`, the fastest pre-existing
oracle, by 3.9–4.7×. At `n = 17` it loses to it (0.85×). The round is
**engineering**: the oracle's exponent in `ℓ` fell and the floor did
not move. The next round's falsification target is stated in the round
README.

```text
WDSAT_BINARY=/path/to/wdsat_solver cargo test --lib \
    wdsat_agrees_with_native_sat_on_prime_degree -- --ignored
```

Build the pinned WDSat binary with the baseline pilot:

```text
python3 research/index_calculus_baseline_20260914/pilot/build_pilot.py \
  --cache-dir /tmp/ec-baseline-cache --build-dir /tmp/wdsat-ic-build
```

## ECC2K-130

Nothing in this round moves the ECC2K-130 attack cost. The decomposition
oracle floor and the rho reference in
`RESEARCH_ECC2K130_DECOMPOSITION.md` still apply. WDSat's static
`config.h` allocation and the `m · 2^{131}` product of the free-oracle
argument both remain binding at full size. What changed is that the IC
pipeline can hand the *same* Semaev instance to Trimoska's solver on the
prime-degree ladder that leads there, and that the quadratic `mq-fes`
path can answer `find_one` without paying a full Möbius transform when a
root appears early in Gray order.

## Verdict

**Engineering unification landed; no advance against the floor.** Native
SAT and WDSat agree on the planted prime-degree toy. Incremental Gray
beats Möbius on `find_one` wall time (engineering, floor ratio flat);
`Fl[0..4]` register blocking and a dedicated no-`Vec` `find_one` walk
are further micro-optimisations of that path.
Monica extends past the Möbius `n ≤ 24` table rather than beating it
inside the cap. An AVX2 4×u64 port of libfes `avx2_8x32` ideas is
correct but **slower** than packed-u64 scalar `L=4` on single-system
instances, so it stays opt-in. Parallel outer specialisation and
hardcoded `L=8` / batch probe likewise stay opt-in.

In the full pipeline, the packed bilinear-split oracle makes
`ic run --solver mq-fes` 13–15× faster end to end at `n = 23` than the
previous oracle, and 3.9–4.7× faster than `pair-table` there. It is
slower than `pair-table` at `n = 17`. The oracle's cost falls from
`4^ℓ` to `2^ℓ·ℓ²`, the free-oracle floor is untouched, and `S` is
unmeasured. Full-size ECC2K-130 index calculus remains above rho for
every oracle this repository has priced.
