# Cycle 3 — claim panel (parity check)

Session: `sessions/claim-n13` (`ECBS1h5834884bc326`), 256/256 verified, audit OK.

## Mean cold S (8 public × 3 rounds)

| arm | mean S | S / strong-rho | lower bound |
|---|---:|---:|---|
| rho-strong | 14.414 | 1.000 | yes |
| ic-f4-m3-d5 | 16.480 | 1.143 | yes |
| **ic-negfold-m2-d5** | **24.462** | 1.697 | no |
| ic-negfold-m2-d4 | 25.992 | 1.803 | no |
| ic-rootfind-m2-d4 | 27.802 | 1.929 | no |
| ic-rootfind-m2-d5 | 30.182 | 2.094 | no |

Paired ops:

- `ic-negfold-m2-d5 / ic-rootfind-m2-d5` = **0.810** [0.797, 0.822] (real set-up win)
- `ic-negfold-m2-d5 / ic-f4-m3-d5` = 1.484 [1.349, 1.628] (**bounded**; not a speedup)

## Parity

**Not achieved.** Best rootfind cold S 24.462 > F4 IC 16.480 (gap ≈ 1.48×).

Secondary amS@64: negfold-m2-d5 ≈ 17.53 vs F4 ≈ 16.34 — near but not under;
labelled secondary only.
