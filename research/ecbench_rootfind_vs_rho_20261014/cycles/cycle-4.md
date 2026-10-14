# Cycle 4 — d=6 negfold + window-negfold

Session: `sessions/cycle4` (`ECBS1h7331f9437767`), 36/36 verified,
audit OK, 8 replays identical.

## Hypothesis

(A) denser `d=6` cuts trials enough to offset larger pair-table setup.
(B) `rootfind-window-negfold` m=3 `window_div=2|3` cuts relation search
vs complete negfold m=2 d=5.

## Outcome: **dead end** (cold S)

| arm | mean S | LB? | med wall ms | S / rho |
|---|---:|---|---:|---:|
| rho-strong | 14.335 | yes | 0.079 | 1.000 |
| ic-c2-frob-counted-m2 | 15.690 | yes | 0.168 | 1.094 |
| **ic-negfold-m2-d5** | **24.289** | **no** | 0.274 | 1.694 |
| ic-winneg-m3-d5-div3 | 28.480 | no | 0.366 | 1.987 |
| ic-winneg-m3-d5-div2 | 30.638 | no | 0.257 | 2.137 |
| ic-negfold-m2-d6 | 39.087 | no | 0.347 | 2.727 |

d=6 and both window-negfold arms are worse than d=5 negfold on cold S.
Secondary amS@64 for d=6 ≈ 13.5 (beats rho cold) — not a cold claim.
