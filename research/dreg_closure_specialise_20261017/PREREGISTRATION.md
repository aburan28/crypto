# Pre-registration: does the closure's rise follow the largest summand dimension or the total?

Registered 2026-10-17, before any specialised system below was measured. It
is the cheap separating test the closure-driver study
(`research/dreg_closure_driver_20261016/RESULTS.md`) could not run directly.

## Why

Two explanations survive the closure-at-scale and closure-driver studies
(56 draws). The MutantXL closure stops refuting at degree 4:

- **H_max:** once some summand's subspace dimension reaches **7**;
- **H_sum:** once the summands' dimensions sum to **21 or more**.

At m = 3, ℓ = 7 both predict the observed rise. Native cells that
separate them (m = 4, ℓ = 6) are out of reach in this container: at the
needed surplus, unsatisfiable draws are essentially absent.

**Specialisation separates them on the draws that already rise.** In the
chained system, summand `i` is `Σ_j c_{iℓ+j} v_j` over the shared subspace
basis `v_0 … v_{ℓ−1}`, so variables `iℓ … iℓ+ℓ−1` are its coordinates.
Fixing coordinate `c_{i,ℓ−1}` confines that summand to 2^{ℓ−1} points:

- to the 6-dimensional subspace when the value is 0;
- to its coset when the value is 1.

Every specialisation of an unsatisfiable system is unsatisfiable. So each
branch is a well-posed refutation problem in which one summand has
dimension 6.

## Draws (frozen)

All 8 unsatisfiable draws of the closure-at-scale study's m = 3, ℓ = 7
cells:

- `research/dreg_closure_scale_20261016/draws/m3-19-7.jsonl`, draws 0, 3, 5
  and 7 (N = 40);
- `research/dreg_closure_scale_20261016/draws/m3-21-7.jsonl`, draws 0, 2, 3
  and 6 (N = 42).

On every one of them the unspecialised degree-4 closure saturated without
refuting.

## Arms (per draw; each a separate process)

Variable indices are for m = 3, ℓ = 7.

| arm | `--fix` | summand dimensions | Σ | max | H_sum | H_max |
|---|---|---|--:|--:|---|---|
| A0 | `6=0` | 6, 7, 7 | 20 | 7 | refutes at 4 | does not |
| A1 | `6=1` | 6 (coset), 7, 7 | 20 | 7 | refutes at 4 | does not |
| B | `6=0,13=0,20=0` | 6, 6, 6 | 18 | 6 | refutes at 4 | refutes at 4 |
| C (null) | `21=0` | 7, 7, 7 (one intermediate coordinate fixed) | 21 | 7 | does not | does not |

- **B is the positive control.** If fixing one coordinate per summand does
  not reproduce the native ℓ = 6 behaviour, refuting at 4, then
  specialisation does not emulate a smaller ℓ, and arm A is uninformative.
- **C is the null control.** It removes one unknown without changing any
  summand. If C refutes at 4, the degree-4 failure is fragile to any
  fixed variable, and arm A is uninformative.

## Instrument

- `examples/dreg_closure --fix …` substitutes the constants into the
  equations, cancels equal terms, and then runs the unchanged dense closure
  scan at degrees 3 and 4.
  - The command is
    `prlimit --as=15000000000 timeout 7200 dreg_closure --to 4 --max-cols 130000 --draw K --fix …`.
- **Pre-registration check.** Fixing a variable that does not occur (60)
  on the (16, 6) draw 0 reproduced the recorded outcome and both recorded
  ranks exactly.
- **XL identity check.** It is skipped under `--fix`, because the committed
  XL outcome is for the unspecialised system.

## Verdict rules

Per run:

- **refutes at 4:** the closure refutes at degree 3 or 4;
- **does not:** degree 4 is saturated without `1`;
- **unmeasured:** the run was killed.

**Controls valid** when B refutes at 4 on ≥ 7 of 8 draws and C does not
refute on ≥ 7 of 8.

If the controls are valid, arm A has 16 runs (A0 and A1 × 8 draws):

- **H_sum favoured** if ≥ 14 of the 16 refute at 4;
- **H_max favoured** if ≤ 2 of the 16 refute;
- **mixed** otherwise. A0 and A1 are also reported separately.

If the controls are not valid, the study is reported as **uninformative
for H_sum vs H_max**, with every run shown.

## Predictions

- B refutes at 4: high confidence.
- C does not: high confidence.
- A: no prior between H_sum and H_max (50/50).

## Scope and class

- **Class:** stage diagnostic. No `S`, no rho ratio, no speed claim.
- **Container:** 4-core, 16 GB, not isolated (`AGENTS.md` §10).
- **The result chooses an extrapolation rule** for the closure degree. It
  measures no ECC2K-130 cost. The descent route's verdict against rho is
  degree-independent.
- **Not a hybrid algorithm.** The branches here are probes of the
  mechanism, not a cost proposal. Fixing coordinates to bring ℓ down to 6
  at the priced ℓ ≈ 33 would cost `2^{m(33−6)}` guesses.
- **Conductor:** no task was taken, because the CLI is not installed in
  this container.
