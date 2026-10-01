# Symmetry-lever slope test: results

This was written after the run. [PREREGISTRATION.md](PREREGISTRATION.md) is unchanged
since its registration commit `ea9f9592`. The readout is
[runs/readout.txt](runs/readout.txt).

## Run

| | |
|:--|:--|
| binary | `sym_degree_ladder-9dcdadd5`, sha256 `193d61a8…1161` (`runs/binary.sha256`) |
| wall | 2026-09-29 17:38:41Z to 18:28:43Z; pinned to CPU 3 under the benchmark lock |
| cells | all 14 exited 0 inside their limits, 4 unsatisfiable draws each; no censoring |

## Registered verdict: constant lever

**Torsion symmetrisation does not flatten the refutation degree.** `D` rises by about 1
per unit `ℓ` on every curve: `s̄ = 1.033`, far above the 0.35 threshold. Both `ℓ = 6`
cells read `≥ 8` on every draw.

| `ℓ` | unknowns | `icv1-f2m13-t181-515ee569` | `icv1-f2m17-tm101-00378d4e` | `icv1-f2m19-tm797-9c54981b` |
|--:|--:|:--|:--|:--|
| 2 | 4 | 4 4 4 4 | 4 4 4 4 | 4 4 4 4 |
| 3 | 7 | 5 5 5 5 | 4 5 4 5 | 5 4 4 4 |
| 4 | 10 | 6 6 6 6 | 6 6 5 6 | 5 5 5 5 |
| 5 | 13 | 7 7 7 7 | 7 7 7 7 | 6 7 7 7 |
| 6 | 16 | — | ≥8 ×4 | ≥8 ×4 |
| **`s_n`** | | **1.000** | **1.100** | **1.000** |

The system degree is 4 in every cell. On the refuting draws, the Macaulay degree needed
is about `ℓ + 2`.

## Predictions

1. **Constant lever: holds.**
2. **`D` non-decreasing in `n` at fixed `ℓ`: fails.** At `ℓ = 3` and `ℓ = 4` the median
   is lower at larger `n`: 5, 4, 4 and 6, 6, 5. At fixed `ℓ`, more field equations over the
   same unknowns refute at the same or a lower degree. The prediction was the wrong way
   round, and the prediction was wrong, not the data.
3. **Satisfiable fraction rising with `ℓ`: not supported.** Only one draw in the whole run
   was satisfiable, at `K0n13l4`. The cells stop at four unsatisfiable draws, and at these
   `(n, ℓ)` the expected root count is below one, so the fraction could not be estimated.

## Robustness, not registered

- **The `ℓ = 2` point sits at a ceiling.** There the system has 4 unknowns, and a Boolean
  refutation degree cannot exceed the unknown count. The point therefore sits at its
  maximum and could flatten or steepen the fit.
- **Without it, the verdict is unchanged.** Fitting `ℓ = 3–5` only gives slopes of 1.0,
  1.5 and 1.5.

## Against the chained system

The committed chained `m = 3` ladder ([dreg_ell_grid_20260925](../dreg_ell_grid_20260925/RESULTS.md))
reads 5, 6, 6, `≥ 7` at `ℓ = 2`–5. The two are not directly comparable:

- the fields, draws and bases differ;
- the chained system has `n + 3ℓ` unknowns, while the symmetrised one has `3ℓ − 2`.

**The qualitative reading.** Symmetrisation cuts the unknown count, a large constant, and
starts one degree lower at small `ℓ`. It grows at least as fast in `ℓ`, so it is a
constant-factor lever, not a slope lever.

## What this decides

- **What was gated.** Per §1, an `m = 5` audit was worth building only alongside a lever
  on the degree slope. This is the one lever the library can test, torsion symmetrisation,
  and at `m = 3` it is a constant.
- **Not worth building.** Deriving symmetrised `S₅`/`S₆` so that the lever can run at
  `m = 4`/`5` is therefore not worth its cost on this evidence. **`m = 5` is not started.**
- **What reopens it.** Any of the following, taken from survey §3.2 and the first audit's
  §11:
  - a symmetry that is not tested here: the symmetric-group action on summands, or
    Frobenius-orbit coordinates;
  - a published degree bound for symmetrised or Weil-descended systems;
  - a solver family other than Macaulay/F4 degree growth.
- **Scope.**
  - One lever, `m = 3`, at `n ∈ {13, 17, 19}` and `ℓ ≤ 6`;
  - the solving degree, not word operations or time;
  - no asymptotic statement.
