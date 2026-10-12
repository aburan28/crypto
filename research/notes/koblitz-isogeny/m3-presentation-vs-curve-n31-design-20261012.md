# m = 3 curve-vs-presentation at n = 31 (design, 2026-10-12)

This is an additive amendment to `presentation-vs-curve-design-20261003.md`. The model,
the estimator and the thresholds are unchanged:
- **H0:** the 95 % upper bound on ICC_curve is < 0.05 on every metric.
- **H1:** any cost metric has ICC ≥ 0.10 with its lower bound > 0.
- Otherwise: **inconclusive**.

It follows the n = 19 m = 3 arm (`m3-presentation-vs-curve-design-20261011.md`, which
was inconclusive). This arm switches to the `PinX3` solver and uses 5× the probes, to
tighten the yield bound.

## Arm

| item | value |
|---|---|
| class | n = 31, a₂ = 0, from the walk (the same class as the m = 2 n = 31 arm) |
| sample | 32 curves (the sweep's seeded hash, K_a included) |
| subspaces | geometric V_s, s = 1..4 |
| decomposition | m = 3, l = 10, `M3Solver::PinX3` |
| probes | 2000 full-group probes per cell |
| order | cells shuffled with a fixed seed |
| command | `KOBLITZ_PRES_M3=pin KOBLITZ_PRES_GEOMETRIC=1 KOBLITZ_PRES_FFD_TARGETS=0 koblitz_presentation_sweep 31 0 10 2000 32 1,2,3,4` |
| output | `experiments/koblitz_presentation/pres_31_0_l10_geom_m3pin.jsonl` |
| analysis | `examples/koblitz_presentation_icc.rs` |

## Metrics

- **yield per |F|³.**
- **ln µs per call.**
- **work per call:** (s, p) roots per call. PinX3 records them in `sym_roots`.
- **ln projected solve.**

First fall is not measured.

## Expected precision

The pilot full solves at n = 31 gave about 0.09 relations per probe, so a cell has
about 180 relations.

## Budget

About 128 × 2000 × 7 ms ≈ 0.5 core-h.

## Results (2026-10-12)

The run completed all 128 cells in 872 s on 4 cores.
- **Soundness:** every cell matched BSGS, with 0 inconsistent relations.
- **Yield:** 20661 relations from 256000 probes; median 160 relations per cell.
- **Cost:** median 16.6 ms per call.

Analysis output: `experiments/koblitz_presentation/icc_n31_m3_geom.json`.

### Pre-registered metrics

| metric | ICC_curve | 95 % CI | subspace share |
|---|---|---|---|
| yield per \|F\|³ | 0.000 | [0, 0.116] | 0.000 |
| ln µs per call | 0.048 | [0, 0.216] | 0.011 |
| **work per call** ((s, p) roots) | **0.228** | **[0.032, 0.403]** | 0.009 |
| ln projected solve | 0.067 | [0, 0.219] | 0.017 |

**Verdict by the pre-registered rule: H1, on one cost metric (work per call).** Its ICC
is ≥ 0.10 and its lower bound is above 0. The other three metrics are consistent with
no curve effect, but they do not clear the 0.05 H0 bound.

### Mechanism (post hoc, not pre-registered)

`KOBLITZ_ICC_DIAG=1`, output in `icc_n31_m3_geom_diag.json`:

| diagnostic | ICC_curve | 95 % CI |
|---|---|---|
| factor-base size \|F\| | 0.239 | [0.025, 0.437] |
| PinX3 pins per call | 0.230 | [0.035, 0.406] |
| roots per pin | 0.000 | [0, 0.138] |
| work per call ÷ \|F\| | 0.086 | [0, 0.297] |

- PinX3 tries every x₃ in the factor base, so its work per call is proportional to
  |F|.
- **The curve effect on work is the curve effect on |F|.** Per pinned x₃, the linear
  algebra shows none (roots per pin: 0.000).
- Yield normalised by |F|³ shows none either.

|F| is curve-dependent in this class, consistently across the four V. The same
diagnostic on the earlier data:

| data | ICC(\|F\|) |
|---|---|
| n = 19, m = 3, a₂ = 1 | 0.000 [0, 0.000] |
| n = 23, m = 2, a₂ = 1 | 0.011 [0, 0.085] |
| n = 19, m = 2, a₂ = 1 | 0.000 [0, 0.016] |
| n = 31, m = 2, l = 16, a₂ = 0 | 0.087 [0, 0.202] |

The a₂ = 1 classes have cofactor 2. This class has cofactor 1492 = 4·373, and
factor-base entries are folded by the cofactor projection. A plausible cause is that
the cofactor-torsion structure differs across the isogeny class; the sample spans
volcano depths 0 and 1. That is consistent with the cofactor-balance yield spread in
`yield-spread-results-20261003.md`. **It is not tested here.**

### Reading

At n = 31, the curve changes the size of the factor base its x-coordinates produce
in V, and so the work an attacker spends per probe and the relations they must
collect. It does not change the per-element decomposition algebra, or the yield per
|F|³.

Scope: m = 3, l = 10, geometric V, this class (a₂ = 0, large cofactor) and this
solver. The effect is a factor-base-size effect, not a hardness difference in the
decomposition problem.

**Next test:** ICC of |F| against volcano depth and cofactor-torsion type, and the same
arm on an n = 31 class with a small cofactor, if one exists.
