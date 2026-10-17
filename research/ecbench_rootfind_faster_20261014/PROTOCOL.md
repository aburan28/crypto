# Protocol: ecbench rootfind actual-speed loop (2026-10-14)

Extends PR #1774 / `research/ecbench_rootfind_parity_20261013/`.
Classification: **engineering**. Synthetic public targets only.
No key-recovery, asymptotic sub-rho, or deployed-curve claims.
Multi-target amortization stays secondary (claim hygiene).

## Goal

1. Close **complete-S** gap vs F4 under apples-to-apples costing
   (price the ops F4 left out, and/or drop complete rootfind S).
2. Surface honest **wall** (median solve_wall_ms) beside S in `ecbench table`.
3. Stretch: match/beat rho-strong on the same panel under the same rules,
   or document the hard stop.

## Starting point (parity claim-n13-v2 @ `3f4a87aae0`)

| arm | mean S | lower bound |
|---|---:|---|
| rho-strong | 14.359 | yes |
| ic-c2-frob-counted-m2 | 15.552 | yes |
| ic-f4-m3-d5 | 16.749 | yes |
| ic-negfold-m2-d5 | 23.220 | **no** |

## Dead ends (do not repeat)

orbit divisor subtract; window_div; window@d4.

## Conductor

Task T-2514 on worktree `wt-ecbench-rootfind-improve`.
Evidence: `research/ecbench_rootfind_faster_20261014/`.
Private `CARGO_TARGET_DIR=…/target-faster`.
