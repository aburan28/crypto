# Protocol: ecbench rootfind vs-rho chase (2026-10-14)

Extends PR #1774 / `research/ecbench_rootfind_faster_20261014/`.
Classification: **engineering**. Synthetic public targets only.
No key-recovery, asymptotic sub-rho, or deployed-curve claims.
Multi-target amortization stays secondary (claim hygiene).

## Goal

Close or document the cold complete-S gap vs `rho-strong` on the
same n=13 public panel under the same costing rules (XOR-priced F4
world). Stretch: lower complete IC wall further.

## Starting point (faster claim-n13-v3 @ `f7a4ba497f`)

| arm | mean S | lower bound | med wall ms |
|---|---:|---|---:|
| rho-strong | 14.822 | yes | 0.085 |
| ic-c2-frob-counted-m2 | 15.648 | yes | 0.169 |
| ic-negfold-m2-d5 | 22.972 | **no** | 0.262 |
| ic-f4-m3-d5 | 421.369 | **no** | 154 |

Complete IC / rho ≈ 1.55×. Negfold GAE: oracle_setup ≈ 305 (pair
table), relations ≈ 700 (walk + trials).

## Dead ends (do not repeat)

orbit divisor subtract; window_div alone on full m=3 without fold;
d=3 standard subspace for cold S; comparing LB compact-frob to
complete IC as a speedup.

## Conductor

Task T-2575 on worktree `wt-ecbench-rootfind-improve` (Conductor
daemon unreachable at session start; check retried).
Private `CARGO_TARGET_DIR=…/target-faster`.
