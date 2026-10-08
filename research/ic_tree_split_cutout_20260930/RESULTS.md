# Is the `m = 4` tree-split set cut out at low degree? Results

This was written after the run. [PREREGISTRATION.md](PREREGISTRATION.md) is unchanged since
its registration commit `3f8beb7a`. The readout is
[runs/registered/readout.txt](runs/registered/readout.txt); the kernel recheck is [runs/recheck.txt](runs/recheck.txt).

## Run

| | |
|:--|:--|
| code | commit `3b37b1a1` (`runs/registered/commit`). Its one change after registration is a `mkdir -p` in `run.sh`, a harness fix with no protocol change. |
| kernel | `gf2_span`, sha256 `f21e9bbb…7b7b` (`runs/registered/sha256`) |
| wall | 2026-09-30 01:18Z to 02:04Z, not timed and not pinned (§6 of the registration) |
| cells | all 48 exited 0 inside their limits; no arm hit the 90,000-column cap before its `D_s`; no censoring |

**Lint-only kernel change.** A clippy/rustfmt-only rewrite of `examples/gf2_span.rs`
(`0f093dcb`) was pushed while the run was in progress. It uses iterator loops with the same
arithmetic. The run's binary was built before it and was not replaced until the run had
ended: the new binary's mtime is 02:04:15, and the last cell exited at 02:03:57.

After the run, the rebuilt kernel (sha256 `fc470ed0…3045`) reproduced four registered cells
exactly, excluding the logged `secs`:

- `T2-n17a1-l4-s20260930`
- `T2-n19a1-l8-s20261001`
- `Z-n13a1-l4-s20260930`
- `Z-n15a0-l5-s20260930`

## Registered verdict: CLOSED, for both `T2` and `Z`

**At every cell, the curve's set has the same sampled cut-out degree as a random set of
its size.**

- The gap is 0 in all 40 `T2` cells and all 8 `Z` cells, against every control.
- `D_s` climbs with `ℓ` at the random sets' own rate.
- No excess closure was seen at `D_s(R)`: the false-zero fraction is 0 in every cell.

### `T2` (primary): `D_s(curve)` = `D_s(R)` in every cell

| `ℓ` | `|T2|` (range over 8 series) | `n = 17` (4 series) | `n = 19` (4 series) |
|--:|:--|:--|:--|
| 4 | 25–100 | 2 2 2 2 | 2 2 2 2 |
| 5 | 169–439 | 3 3 3 3 | 3 3 3 3 |
| 6 | 827–1,222 | 4 4 4 4 | 3 4 3 3 |
| 7 | 2,788–4,581 | 5 5 5 5 | 4 4 4 4 |
| 8 | 13,456–16,628 | 6 6 6 6 | 5 5 5 5 |
| slope | | 1.0 in every series | 0.7 in every series |

The median slope is 0.85, against the ALIVE threshold of `≤ 0.15`.

- **At `ℓ = 8`, `T2` is already 10–12% of `GF(2)^17`.** At that size its degree-5 and
  degree-6 behaviour is that of a random set.
- **The slope is not structural.** It is the counting threshold, set by `|T2| ≈ 2^{2ℓ−2}`
  against `M(D)`, which is why it is lower at `n = 19`.

### `Z` (secondary): `D_s(curve)` = `D_s(ZR)` = `D_s(R)` in every cell

| series | `ℓ = 4`: `|Z|`, `D_s` | `ℓ = 5`: `|Z|`, `D_s` |
|:--|:--|:--|
| `n = 13`, `a = 0` | 325, 2 | 18,915, 5 |
| `n = 13`, `a = 1` | 3,081, 4 | 23,220, 5 |
| `n = 15`, `a = 0` | 1,225, 3 | 19,306, 4 |
| `n = 15`, `a = 1` | 666, 3 | 14,365, 4 |

`Sym²` adds no low-degree structure: the `ZR` control is also indistinguishable from `R`.
The curve adds none either.

## The one structure found: a linear equation from the 2-torsion

Two cells show a positive deficit, and it has a known cause:

- **`T2` at `n = 17`, `a = 1`, `ℓ = 4`, seed 20260930:** deficit 1 at `D = 1`. About half
  the samples (1,005 of 2,000) are false zeros there, because the set lies in an affine
  hyperplane.
- **`Z` at `n = 13`, `a = 1`, `ℓ = 4`:** deficits 1, 26 and 326 at `D` = 1, 2, 3.
  - For `N = 26` these are exactly 1, `N` and `C(N, 2) + 1`: the multilinear multiples of
    one linear form, and nothing else.
  - The false-zero fraction is 989/2000 until `D = 4`, where it falls to 0, the same as
    the controls.

**Cause.** On `K₁` with `n` odd, `#E = 2p`, so `2E` is the index-2 subgroup. A point is in
`2E` exactly when `Tr(x) = Tr(a) = 1`.

- In both cells, every factor-base `x` has trace 0; checked directly for these seeds.
  This happens with probability `2^{−ℓ}`.
- So every `P + Q` lies in `2E`, and `Tr(x(P + Q)) = 1` on all of `T2`.
- In `Z`, `Tr(e + f) = 0` then gives the one linear equation.
- A factor base that mixes both cosets removes it: the same curve at `ℓ = 5` has no deficit.

This is the ordinary halving/trace constraint. It lowers no cut-out degree in any cell.

## Predictions (the prior)

The registration's prior was closure (§1): "`T2` is the image of `2ℓ` bits under a map of
algebraic degree well above 1 in `n` bits, and such images usually look random at low
degree." **That holds.**

## What this changes

- **The tree split is a table method at these sizes.** Neither `T2` nor `Z = Sym²(T2)` has
  a low-degree description. The affine-linearity of `S₃` in `(e + f, e·f)` is real, but
  all the degree sits in `Z`, and `Z` needs the degree of a random set of its size.
  - The survey's §3.5 bound for table decompositions (`r^{j/(2j−1)}`) applies to this
    route unchanged.
  - The reopening condition of
    [DECISION-20260929-stop.md](../ic_candidate_tournament_20260915/campaign_20260916/DECISION-20260929-stop.md)
    ("a decomposition that is not a table") stays unmet for the tree split.
- **Not tested, and still open:**
  - Other models of `Z ∩ A_R`: linearisation with auxiliary variables, where
    `(e, f) ↦ (E₁, E₂)` is kept as equations rather than eliminated. That is the chained
    system the `m = 4` audits measured.
  - Nagao / Riemann–Roch decompositions.
  - Frobenius-orbit or summand-permutation coordinates.

## An unregistered observation, reported and unresolved

- **The smoke.** In the pre-registration timing run (`Z` at `n = 11`, `ℓ = 5`, seed 1, curve
  arm only, §5 of the registration), 10 of 2,000 samples were false zeros at `D = 5`.
- **Why it is unresolved.** No control was run for that cell, and `n = 11` is below the
  registered sizes. The registered `n = 13` and `n = 15` cells at `ℓ = 5` show no excess.
  So it is either a small-field effect or a random-set effect at that size; it is not
  resolved here.
- **What it cannot be.** An excess closure means the set is **harder** to cut out than
  random. It cannot be a lever in either case.

## Scope

- **Curves and sizes.** Two Koblitz curves, `n` 13–19, `ℓ` 4–8 for `T2` and `ℓ` 4–5 for `Z`.
  One random subspace per seed: two seeds for `T2`, one for `Z`.
- **Metric.** A sampled cut-out degree, capped at `D = 6` and 90,000 monomials.
  - 2,000 samples cannot see a degree-`D` closure excess below about 0.15% of the space
    (§3 of the registration).
  - The rank and deficit profile, which has no such limit, is also random-like in every
    cell except the two trace cells above.
- **Nothing beyond these sizes is claimed.** In particular, this is no statement about
  `n ≈ 83` or `n = 131`.
