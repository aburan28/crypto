# Unifying Koblitz index calculus with Trimoska WDSat

**Modules:** `src/cryptanalysis/wdsat_oracle.rs`,
`DecompositionStrategy::Wdsat` in `koblitz_index_calculus.rs`,
`ic run --solver wdsat --wdsat-binary PATH`
**Related:** [`RESEARCH_TRIMOSKA_BENCHMARKS.md`](RESEARCH_TRIMOSKA_BENCHMARKS.md),
[`RESEARCH_SAT_SEMAEV.md`](RESEARCH_SAT_SEMAEV.md),
[`RESEARCH_ECC2K130_DECOMPOSITION.md`](RESEARCH_ECC2K130_DECOMPOSITION.md),
[`research/index_calculus_baseline_20260914/`](research/index_calculus_baseline_20260914/)

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

| backend | when used | measured vs Möbius (`n=18`, planted early Gray root) |
|---|---|---|
| Incremental Gray (libfes FFS, `L=4` unroll) | `find_one` / Semaev lift (early exit) | **~273×** wall faster than full Möbius; **~33×** vs prior O(n)/step Gray on full `n=18` enum |
| Möbius transform (`moebius.c`) | `find_all` for `n ≤ 24` | reference for all-roots |
| Monica hybrid (`monica.c`) | `n > 24` (range extension) | does **not** beat Möbius inside `n ≤ 24` (release wall on `n=14,m=32` was ~0.22×); calibrated cost model agrees |

Falsification for the “faster than Möbius” claim: a release run of
`gray_early_exit_beats_moebius_find_one_wall` must keep ratio `≥ 1.5` on
the fixed dense quadratic with a Gray-index-2000 planted root. Falsification
for the Gray speedup itself: `gray_ffs_beats_on_step_full_enum_wall` must
keep FFS/`L=4` ≥ 1.5× the prior O(n)-per-step Gray on full `n=18` enum.
Monica is kept as a capacity extension, not as an in-cap speedup. Cubic
chained (`m ≥ 3`) Semaev systems are refused; those stay on SAT / WDSat.

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
beats Möbius on `find_one` wall time (engineering, floor ratio flat).
Monica extends past the Möbius `n ≤ 24` table rather than beating it
inside the cap. Full-size ECC2K-130 index calculus remains above rho for
every oracle this repository has priced.
