# Disclosed-point F4/F5 standard-subspace d6 pilot — result

Status: **dispatch diagnostic complete; no complete recovery; not a family
qualification**. Comparison kind: `factor-base-policy`.
`promotion_eligible=false`. Competitive `S` / rho ratios / scoreboard rows
remain unset.

## Boundaries (unchanged)

| Boundary | Value |
| --- | --- |
| Encoder floor | `m·ℓ+(m−2)·n` with `m=3`, `ℓ=6` ≤ 64 (35–49 on the five cells) |
| Ambient-orbit contrast | Registered v2 F4/F5 cells need `4n` vars and are statically encoder-`unsupported` ([STATIC-FEASIBILITY](../generic-backend-qualification-v2/STATIC-FEASIBILITY.md)) |
| Reference | Optimized pair-table incumbent not re-timed here; wiring only |

## Measured series

| Item | Binding |
| --- | --- |
| Worker | commit `765c3c5f19032bd852163805f257c56babef2040`, rustc 1.94.1 |
| Build | `build-identity.json` — worker SHA-256 `ed8e5a2f7aa5924e81587b9b0ed5bb11935166e226db640b720bc28842faf798` |
| Points | Five inventory-control targets; inventory JSON SHA-256 `be3053b4b637311e8255f294807fde412510fbda2464e678a6a16fdea397d403` |
| Budget | `max_trials=1`, `batch_trials=1`, child wall 180s, `RAYON_NUM_THREADS=1`, no `KIC_*` |
| Arms | `f4` and `f5`, dense relation LA, `standard_subspace` dimension 6 |
| Manifest | [`runs/manifest.json`](runs/manifest.json) |

Retained operational failures (not solver mathematics):

- `runs-env-error/`: exclusive-phase worker rejects undeclared `KIC_F5_AVX512_UNPACK`.
- `runs-max64-timeout/`: `max_trials=64` hit the 300s cap on n17a1 with empty stdout (no receipt).

## Table (one unit: verified complete one-target IC solve under the frozen budget)

| Cell | Solver | Collector strategy | Engine | `unsupported` | PDP outcome | Relations | Solutions | Wall s | Class |
| --- | --- | --- | --- | --- | --- | --- | --- | ---: | --- |
| n17a1 | f4 | Groebner | MatrixF4 d≤3 | false | incomplete | 0 | 0 | 8.54 | engineering / diagnostic |
| n17a1 | f5 | Groebner | MatrixF5 d≤3 | false | incomplete | 0 | 0 | 5.43 | engineering / diagnostic |
| n19a0 | f4 | Groebner | MatrixF4 d≤3 | false | incomplete | 0 | 0 | 15.36 | engineering / diagnostic |
| n19a0 | f5 | Groebner | MatrixF5 d≤3 | false | incomplete | 0 | 0 | 10.15 | engineering / diagnostic |
| n23a0 | f4 | Groebner | MatrixF4 d≤3 | false | incomplete | 0 | 0 | 23.80 | engineering / diagnostic |
| n23a0 | f5 | Groebner | MatrixF5 d≤3 | false | incomplete | 0 | 0 | 13.01 | engineering / diagnostic |
| n23a1 | f4 | Groebner | MatrixF4 d≤3 | false | incomplete | 0 | 0 | 23.45 | engineering / diagnostic |
| n23a1 | f5 | Groebner | MatrixF5 d≤3 | false | incomplete | 0 | 0 | 14.21 | engineering / diagnostic |
| n31a0 | f4 | Groebner | MatrixF4 d≤3 | false | incomplete | 0 | 0 | 128.77 | engineering / diagnostic |
| n31a0 | f5 | Groebner | MatrixF5 d≤3 | false | incomplete | 0 | 0 | 42.06 | engineering / diagnostic |

Every job exited with report `status=incomplete`, exit code 2, and
`unsupported_false=1` / `unsupported_true=0` on the single ordinary query.
Node-budget exhaustion (`exhausted: true` on the n17a1-f4 probe) is incomplete,
not UNSAT and not encoder rejection.

## Decision

1. **Hypothesis on dispatch: held.** Dimension-6 `standard_subspace` clears the
   Semaev encoder cap on all five disclosed cells; both `f4` and `f5` enter the
   named MatrixF4/MatrixF5 engines with `unsupported: false`.
2. **Hypothesis on bounded recovery under `max_trials=1`: failed.** Zero
   accepted relations and zero verified solutions on every cell. A larger trial
   budget is a separate amendment; the max-64 series timed out without receipts
   on n17a1 within 300s, so recovery cost under the v2-style 65536 budget is
   not established here.
3. **Not a qualification.** Do not rank F4/F5 against the incumbent or rho from
   this panel. Any later competitive registration needs a new seed, full v2
   exposure exclusions, `generic_solver_feasibility.py --require-pass`, and the
   `factor-base-policy` label.

## Next command

[Actions 36580669479](https://github.com/aburan28/crypto/actions/runs/36580669479)
is [operationally censored](../generic-backend-qualification-v2/RESULT.md); never
redispatch seed `2026092902`. If a
fresh F4/F5 complete-solve campaign is justified, register a new protocol with
feasible layouts and withhold promotion until verified full-panel recovery.
