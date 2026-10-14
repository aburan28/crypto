# Cycle 5 — d=4/d=7 negfold + plain rootfind

Session: `sessions/cycle5` (`ECBS1h4b924adb208a`), 36/36 verified,
audit OK, 8 replays identical.

## Hypothesis

(A) d=4 or d=7 moves the setup/hit-rate trade-off favorably.
(B) plain `rootfind` (full PairSumTable) at d=5 beats negfold despite
larger setup.

## Outcome: **dead end** (cold S); d=5 negfold still best complete

| arm | mean S | LB? | med wall ms | S / rho |
|---|---:|---|---:|---:|
| rho-strong | 14.145 | yes | 0.084 | 1.000 |
| ic-c2-frob-counted-m2 | 15.456 | yes | 0.169 | 1.093 |
| **ic-negfold-m2-d5** | **21.803** | **no** | **0.229** | **1.541** |
| ic-negfold-m2-d4 | 26.595 | no | 0.265 | 1.880 |
| ic-rootfind-m2-d5 | 27.523 | no | 0.383 | 1.946 |
| ic-negfold-m2-d7 | 108.913 | no | 1.347 | 7.700 |

Plain rootfind loses to negfold (~1.26× worse). d=7 is dominated by
pair-table setup. Hard stop vs rho unchanged at ~1.54× for best
complete IC.
