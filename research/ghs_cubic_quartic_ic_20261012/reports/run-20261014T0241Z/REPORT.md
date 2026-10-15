# Run report 20261014T0241Z

Recorded with `rustc 1.93.1 (01f6ddf75 2026-02-11) (Homebrew)` on macOS (Apple silicon),
`rustc --edition 2021 -O quartic_ic.rs`. Command:

```bash
for q in 5 7 11; do ./quartic_ic $q 10 rho; ./quartic_ic $q 10 full; ./quartic_ic $q 10 2lp; ./quartic_ic $q 1 --check-lin-equiv; done
```

Results: [../../results/run-20261014T0241Z.jsonl](../../results/run-20261014T0241Z.jsonl).
Plot: S_mean.png / S_mean.svg. The prior smoke run `run-20261014T0200Z` was reproduced
exactly (identical `k_planted`, `ops`, `S`) before this run was recorded.

90/90 trials `correct:true` (10 targets per variant per field, all three fields
`q ∈ {5, 7, 11}`); IC recovery is `recovery:ic` on every trial (no
`transfer_pushforward_dlog` fallback). `--check-lin-equiv` (`pic_jacobian`) passes
at `q = 5, 7, 11`.

| q | variant | trials | mean S (ops/√ℓ) | mean rho steps |
|---|---------|--------|-----------------|----------------|
| 5 | rho | 10 | — (steps only) | 9.0 |
| 5 | full / 2lp | 10 each | 650 896.03 | — |
| 7 | rho | 10 | — (steps only) | 11.5 |
| 7 | full / 2lp | 10 each | 2 751 662.59 | — |
| 11 | rho | 10 | — (steps only) | 30.7 |
| 11 | full / 2lp | 10 each | 37 920 187.89 | — |

Classification: toy planted-scalar diagnostics. The rho arm records step counts,
not priced operations, so no IC-vs-rho `S` ratio or "below rho" claim is made
from this vehicle. No deployed-curve or sub-ρ claims.
