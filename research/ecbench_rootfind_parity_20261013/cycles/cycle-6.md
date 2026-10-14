# Cycle 6 — claim panel v2 (compact-orbit counted frob)

Session: `sessions/claim-n13-v2` (`ECBS1h2ef2bb46ad6c`), 256/256 verified
(8 arms × 8 public × 3 rounds + controls), audit OK, 12 replays.

## Mean cold S

| arm | mean S | 95% CI | lower bound |
|---|---:|---|---|
| rho-strong | 14.359 | [14.021, 14.878] | yes |
| **ic-c2-frob-m2** | **14.336** | [14.217, 14.457] | yes |
| **ic-c2-frob-counted-m2** | **15.552** | [15.433, 15.673] | yes |
| ic-f4-m3-d5 | 16.749 | [16.150, 17.299] | yes |
| ic-c3-frob-counted-m2 | 18.643 | [18.519, 18.765] | yes |
| ic-negfold-m2-d5 | 23.220 | [21.447, 25.416] | **no** |

Paired ops (diagnostic; unpriced work on both sides):

- `ic-c2-frob-counted-m2 / ic-f4-m3-d5` = **0.928** [0.901, 0.960]
- `ic-c2-frob-counted-m2 / ic-negfold-m2-d5` = **0.670** [0.612, 0.726]

## Parity

**Achieved under harness admitted cold S** (both arms lower-bound):
counted compact-orbit frob 15.552 ≤ F4 16.749.

**Not** a complete-S vs complete-S win: negfold (complete) remains 23.220.
F4 still omits word XORs; compact omits inversions + FB-construction
natives. Secondary pinned-XOR sensitivity on F4 raises its mean to
~421 (`tables/f4-sensitivity-claim-n13-v2.txt`) — would widen the gap
further in rootfind’s favour, but is not the admitted figure.

Amortized secondary: counted-c2 amS@64 ≈ 11.16 vs F4 ≈ 16.56.
