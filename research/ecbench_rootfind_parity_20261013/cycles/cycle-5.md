# Cycle 5 — compact-orbit counted frob (columns 2/3)

Session: `sessions/cycle5` (`ECBS1h3600bd37dbee`), 42/42 verified, audit OK.

## Hypothesis

`mitm-frobenius-counted:m=2` on `compact-orbit-scan:columns=2` prices
table canons/frobenius/lookups and still lands admitted cold S ≤ F4
m=3/d=5 on the same panel.

## Mean cold S (2 public × 2 rounds)

| arm | mean S | lower bound |
|---|---:|---|
| ic-c2-frob-m2 | 14.055 | yes |
| **ic-c2-frob-counted-m2** | **15.271** | yes |
| ic-f4-m3-d5 | 15.950 | yes |
| ic-c3-frob-counted-m2 | 18.344 | yes |
| ic-c3-frob-counted-m3 | 18.382 | yes |
| ic-negfold-m2-d5 | 24.134 | no |

(rho-strong mean skewed by one L0 outlier S≈108; ratios vs rho not used.)

## Decision

**Provisional parity vs F4 on admitted (lower-bound) S:** counted c2
15.271 ≤ F4 15.950. Both arms still leave work unpriced (F4: word XORs;
compact: inversions + FB-construction natives). Complete-S negfold remains
above F4. Promote to claim panel (cycle-6 / claim-n13-v2).
