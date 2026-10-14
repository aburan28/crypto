# Cycle 1 — pin inherited-F4 word-XOR op_unit (apples-to-apples)

**Change:** `is_word_xor_op_unit` matches
`word XORs (elimination, specialisation and linear elimination only)`
so pinned `ns_per_word_xor` prices F4 into admitted S (`priced_by=pinned`).
Previously only the bare unit `word XORs` matched, so F4 stayed a permanent
lower bound beside a calibrated ratio.

Session: `sessions/cycle1` (`ECBS1haab050283f6e`), 24/24 verified, audit OK,
8/8 replays identical. Binary SHA-256
`27fd2653b706b6e82e8da1c39a8c26c8f686e3f6e4831953b85eb9798610b76c`.

## Mean cold S (smoke, 2 public × 2 rounds)

| arm | mean S | lower bound | med wall ms |
|---|---:|---|---:|
| rho-strong | 13.961 | yes | 0.073 |
| ic-c2-frob-counted-m2 | 15.513 | yes | 0.284 |
| **ic-negfold-m2-d5** | **22.228** | **no** | **0.244** |
| **ic-f4-m3-d5** | **540.708** | **no** | 224.530 |

## Outcome

**Win — complete-S parity under apples-to-apples costing:**
complete negfold 22.228 ≤ complete F4 540.708 (ratio ≈ 0.041).
Wall (L0 descriptive): negfold / F4 ≈ 0.0011.

Rho still wins cold S. Compact counted remains best LB IC arm (~15.5).
Classification: engineering.
