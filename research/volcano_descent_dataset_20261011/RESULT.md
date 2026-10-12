# Result: tall isogeny volcanoes, descent algorithm, dataset, per-layer ECDLP cost

**Date.** 2026-10-11. **Class.** Accounting and measurement at toy sizes. **Not claimed:** any
ECDLP speedup, anything about a deployed curve, anything about binary fields beyond one step.
Method and file guide: [`README.md`](README.md). Every number below is printed by
`analyze.py` into [`RESULT_TABLES.md`](RESULT_TABLES.md) from the files in `data/`.

## 1. Answer

1. **The descent algorithm works and is cheap.** From any curve it reaches a floor curve, and
   reads off the endomorphism-ring conductor, in a few milliseconds even at depth 20.
2. **Tall volcanoes are easy to build and explore exactly.** 102 volcanoes, 565,026 curves
   explored, 204,922 stored. Deepest in the dataset: 15 levels. Stress runs: 20 levels
   (1,048,576 curves, 71 s) and, for ell = 3, 12 levels (531,441 curves, 52 s).
3. **Every label is independently checked.** All 102 volcanoes pass all five checks
   (section 3). 204,922 of 204,922 measured group structures match the ideal-quotient
   prediction. 2,448 of 2,448 blind descents reach a correct floor curve.
4. **ECDLP cost with a generic solver does not depend on the layer.** Pohlig-Hellman plus
   rho, with or without the negation map, shows no layer difference larger than its bootstrap
   intervals, about +-8% (section 5). The one measured exception is the extra units at
   `j = 0` and `j = 1728`, which exist only at the crater.
5. **Observables do not reveal the level.** Where every curve in a volcano has a cyclic group,
   no cheap feature beats a constant guess. The group structure reveals the level only when it
   varies with level. Section 6.

## 2. Dataset

| | |
|---|---|
| volcanoes | 102 (ell = 2, 3, 5, 7, 11, 13) |
| curves explored exactly / stored | 565,026 / 204,922 |
| depth range | 3 to 15 levels (stress: 20) |
| field size | 40 to 68 bits |
| crater conductors c > 1 | yes (c = 2, 3, 5) |
| level-dependent group structure | 80 of 102 volcanoes |
| size on disk | about 11 MB (gzip CSV) |

Columns and their roles (labels, observables, graph-derived leaks) are in the README.

## 3. Verification

| Check | Result |
|---|---|
| V1 level sizes equal class numbers `h(D_K (c ell^i)^2)` | 102 of 102 |
| V2 vertex degree is ell + 1 (non-floor) or 1 (floor) | 102 of 102 |
| V3 measured group structure equals `O_d/(pi - 1)` prediction | 102 of 102 volcanoes, 204,922 of 204,922 curves |
| V4 point count equals `p + 1 - t` (PARI SEA, sampled) | 102 of 102 |
| V5 blind descent matches BFS truth | 102 of 102 volcanoes, 2,448 of 2,448 descents |
| unit tests | 11 of 11 |

Two bugs were found by the tests and fixed before the dataset was built, and are covered by
regression tests:

* **Degenerate crater.** At `j = 0` or `j = 1728` the extra automorphisms merge several
  kernels onto one j-invariant. A walk keyed on distinct roots saw the crater as a floor and
  stopped. Descent now counts kernels with multiplicity.
* **Chinese remaindering.** `pari.chinese` takes two arguments. Replaced by a plain CRT.

The first dataset build used `ellgroup`, which intermittently hung on some curves and drew from
PARI's unseeded generator. Both group structure and point sampling now use a seeded generator
and point orders against the known factorisation of N. The first build was discarded and the
full campaign rerun from scratch (219 s).

## 4. ECC2K-130, one real step

`ecc2k130_level1.py` builds `Phi_263(X, 1)` over `F_(2^131)` for the Koblitz curve
`y^2 + xy = x^3 + 1`. It has degree 264, took 61 s to build and 0.36 s to factor. It splits
completely: 264 linear factors, `j = 1` with multiplicity 2 (the two horizontal 263-isogenies)
and 262 further distinct roots. `h(O_263) = 262`, so those are exactly the level-1 curves of
the 263-volcano. This confirms that `pi_q` acts as a scalar on `E[263]`.

**Not done.** (a) The second step. In characteristic 2 PARI's `polmodular` accepts only a
prime-subfield evaluation point, so a level-1 curve cannot be expanded; its single rational
neighbour is predicted by theory, not computed. (b) The P-volcano, `P = 146505763881528721`,
needs an isogeny of degree about 2^57 and is out of reach. (c) A first attempt printed degree-131
factors. PARI had factored over `F_2` instead of `F_(2^131)`; the cause was found and the script
now uses `factorff` over the full field.

**Correction to an earlier statement in this session.** I said factoring `Phi_263` over
`F_(2^131)` takes well under a second. Factoring does; building the polynomial takes about a
minute.

## 5. ECDLP cost per layer

Solver: pure-Python Pollard rho with distinguished points and Pohlig-Hellman, counted unit one
affine group operation. This is not the ecbench unit. Largest prime factor of the point order
`r` is 28 to 32 bits. Costs are `ops_rho / sqrt(pi r / 2)`, so volcanoes of different `r` are
comparable. Bootstrap intervals resample solves and ignore clustering by vertex.

**E1, generic layers.** 6 volcanoes (ell = 2, depth 12, `D_K` not -3 or -4), 152 vertices,
4,864 solves, all correct.

| solver | crater | middle | floor |
|---|---|---|---|
| plain rho | 0.859 [0.787, 0.934] | 0.914 [0.895, 0.934] | 0.956 [0.887, 1.03] |
| negation map | 0.624 [0.566, 0.684] | 0.675 [0.660, 0.690] | 0.712 [0.662, 0.764] |
| plain / negation | 1.38 | 1.35 | 1.34 |

The slope of cost against level is within +-0.011 per level in every volcano. The ell-power part
of the group is removed by Pohlig-Hellman: at `m = 12` the naive `sqrt(ord G)` would be
8.2e8 operations and the measured total is 6.2e4. I read this as no layer dependence for a
generic solver at this resolution. The crater and floor intervals for negation (0.566 to 0.684 and
0.662 to 0.764) overlap, so I do not claim a crater-to-floor difference.

**E2, automorphisms.** 4 volcanoes (`D_K` = -3 and -4, with `m` = 0, 6, 12), 96 vertices, 8,352
solves, all correct. Only the crater has `j = 0` (mu_6) or `j = 1728` (mu_4).

| `D_K` | layer | negation / full orbit | plain / full orbit |
|---|---|---|---|
| -3 | crater | 1.48 [1.33, 1.64] (ideal 1.73) | 2.01 [1.82, 2.24] (ideal 2.45) |
| -3 | not crater | 1.00 [0.96, 1.05] | 1.38 [1.33, 1.44] |
| -4 | crater | 1.31 [1.20, 1.44] (ideal 1.41) | 1.84 [1.68, 2.02] (ideal 2.00) |
| -4 | not crater | 1.01 [0.97, 1.05] | 1.33 [1.27, 1.38] |

In all four volcanoes the crater lower bound exceeds 1. The solver reaches about 80 to 85% of the
ideal gain, so these are measurements of this implementation, not of the theoretical limit.
Control: for non-crater layers the `neg` and `aut` solvers run identical code, so the true ratio
is 1. The middle row covers 1; the floor row, 0.875 [0.79, 0.97], does not. That means the
intervals are optimistic for layers with few vertices, and the crater intervals should be read
that way too.

**Superseded run.** The first E1 run (`data/superseded/`) used a solver with two scalar
multiplications per trail start and no fruitless-cycle escape. It showed plain/negation of 1.02
instead of 1.4 and was discarded as a solver defect. A first E2 run was stopped for the same
reason before producing output.

**Interpretation (separate from the measurements).** For a generic solver the only
layer-dependent inputs are the group structure and the units of the endomorphism ring. The
group structure is neutralised by Pohlig-Hellman. The units matter only at the crater and only
for `D_K` = -3, -4. Whether a non-generic attack distinguishes layers is untested here.

## 6. What can be predicted from observables

Target: round(10 * height / k), 11 classes, held-out volcanoes, gradient boosting. Chance is
the constant-prediction balanced accuracy.

| volcanoes | features | balanced accuracy | chance |
|---|---|---|---|
| all curves cyclic (22) | j, a, b statistics | 0.111 | 0.111 |
| all curves cyclic | + degrees for other ell | 0.110 | 0.111 |
| all curves cyclic | + group structure | 0.110 | 0.111 |
| all curves cyclic | + own-ell degree (graph-derived) | 0.279 | 0.111 |
| group varies with level (80) | j, a, b statistics | 0.130 | 0.091 |
| group varies with level | + degrees for other ell | 0.146 | 0.091 |
| group varies with level | + group structure | 0.576 | 0.091 |

Two cautions for anyone building on this. First, the small lift from raw j, a, b is an artefact
of CM construction: 80 of the 134 stored craters have a tiny signed j-invariant (for example
8000), against 196 of 204,788 other rows. Second, `m` in `volcanoes.csv.gz` is a requested lower
bound; use the observed `n1` to decide whether a volcano is cyclic.

## 7. Requirement status

| Requirement | Status | Evidence or gap |
|---|---|---|
| Build the descent algorithm | verified complete for prime fields | tests, 102 volcanoes, 2,448 descents |
| Same for binary fields | partial | one step only (section 4); second step blocked by PARI |
| Try it on tall volcanoes | verified complete | depth 15 in the dataset, 20 in the stress run |
| Find curves with a lot of depth | verified complete | construction from `p = (t^2 - D_K f^2) / 4`; depth limited by `p >= \|D_K\| f^2 / 4` |
| Build a dataset | verified complete | `data/`, schema in README |
| Use ML for endomorphism-ring properties | partial | dataset and a baseline only; no model for the user's question was trained beyond section 6 |
| ECDLP hardness across layers | partial | generic and automorphism solvers at `r <= 2^32`; no non-generic attack, no larger sizes |
| ecbench integration | not attempted | counted unit differs; do not mix |
| Binary-field tall volcanoes (Koblitz) | not attempted | needs a bivariate `Phi_ell` mod 2 |
| `cargo test --release --lib` | not run | no Rust file touched |
| Conductor reservation | blocked | control plane unreachable; warn-through path |
