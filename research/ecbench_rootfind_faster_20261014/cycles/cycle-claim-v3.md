# Claim panel v3 — XOR-priced F4 vs complete negfold

Session: `sessions/claim-n13-v3` (`ECBS1h1a87cdf5adc1`), 224/224 verified
(7 arms × 8 public × 3 rounds + controls), audit OK, 12 replays.

## Mean cold S

| arm | mean S | 95% CI | lower bound | med wall ms |
|---|---:|---|---|---:|
| rho-strong | 14.822 | [14.174, 16.324] | yes | 0.085 |
| ic-c2-frob-counted-m2 | 15.648 | [15.461, 15.831] | yes | 0.169 |
| **ic-negfold-m2-d5** | **22.972** | [21.265, 24.816] | **no** | **0.262** |
| ic-negfold-m2-d4 | 25.729 | [22.163, 29.289] | no | 0.327 |
| **ic-f4-m3-d5** | **421.369** | [356.623, 494.193] | **no** | 154.013 |

Paired ops (complete/complete):

- `ic-negfold-m2-d5 / ic-f4-m3-d5` = **0.0545** [0.047, 0.064]
- wall median ratio ≈ **0.0017**

Vs rho (diagnostic; rho still has unpriced natives):

- `ic-negfold-m2-d5 / rho-strong` = 1.550 [1.381, 1.710]
- `ic-c2-frob-counted-m2 / rho-strong` = 1.056 [0.960, 1.109] (LB)

## Parity status

| question | status |
|---|---|
| complete-S ≤ F4 (apples-to-apples) | **Y** — negfold 22.97 ≤ F4 421.37 |
| cold S ≤ rho-strong | **N** — hard stop; complete IC ~1.55×; LB compact ~1.06× |
| wall vs F4 | **Y** (L0 descriptive) — ~600× lower |

d=4 does not beat d=5 on the claim panel (25.7 vs 23.0).
