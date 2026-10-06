# P-256 low-delta cyclic factor-base tradeoff, round 31: result

Date run: 2026-10-06

A two-delta known-log factor base does not reach Pollard-rho parity.  Reducing
the exact hidden update state from round 30's 139,592 classes to two makes the
local collision constant attractive, but simultaneously confines signed
17-sums to too few arithmetic runs.

At the largest rare-edge count still compatible with the hidden-state parity
budget, `R=4778`, the deliberately generous support upper bound covers only
`2^-9.131045` of the group.  It therefore requires at least 560.684 target
retries and costs at least **560.683 times rho**.  The support bound first
reaches the group at `R=6934`, after the hidden-state ratio has risen to
1.016460; that cycle is non-primitive.  The best primitive cycle is `R=6935`
at an optimistic lower bound of **1.016468 times rho**.

That last number is not an achieved attack cost.  It gives the candidate full
group coverage even though only an upper bound on support reached the group.
The actual support may be much smaller, as it is in every complete toy cell.
No P-256 usable-relation probability is proved and no full-depth relation was
attempted.

## Complete fixed-size sweep

The sweep evaluated all 65,729 integer choices
`1 <= R <= floor(131458/2)`.  For `p=(B-R)/B` and `q=R/B`, it credited the
optimistic two-delta collision factor `1/sqrt(p^2+q^2)` and the frozen
round-26 local-oracle constant.  It then charged at least `n/U(R)` target
trials, where

```text
U(R) = min(n, (2R)^17 * (34B+1)).
```

`U` overcounts ordered terms, signs, duplicate columns and offset collisions;
it is a support upper bound, never a measured yield.

| row | R | gcd(B-R,R) | hidden / rho | support upper bits | retry lower bound | complete lower bound / rho |
|:--|--:|--:|--:|--:|--:|--:|
| one rare edge | 1 | 1 | 0.964344 | 39.091706 | 1.9765e65 | 1.9061e65 |
| last primitive parity-side landmark | 4,777 | 1 | 0.999990 | 246.863821 | 562.683 | 562.677 |
| largest hidden-only parity row | 4,778 | 2 | 0.999997 | 246.868955 | 560.684 | 560.683 |
| first support saturation / global minimum | 6,934 | 2 | 1.016460 | 256.000000 | 1 | 1.016460 |
| **primitive-cycle minimum** | **6,935** | **1** | **1.016468** | **256.000000** | **1** | **1.016468** |
| quarter-width landmark | 32,864 | 2 | 1.219796 | 256.000000 | 1 | 1.219796 |
| half-width landmark | 65,729 | 65,729 | 1.363778 | 256.000000 | 1 | 1.363778 |

The ordered sweep digest is
`fe3ba23e8d7fa2fc74fb34325d70974593f44400bbea839778b89b07a62570c8`.
Concentrating all rare mass in one second delta maximizes the credited
collision probability, so splitting that same mass across more delta classes
cannot improve this hidden-state boundary.

## Exact native coefficient replay

For every registered landmark the binary solved

```text
(B-R) * 1 + R * b = 0 mod n
```

and distributed the `R` exceptional edges with the frozen mechanical word.
Every valid row contains 131,458 distinct, nonidentity P-256-order scalar
coefficients, unique up to sign.  Replaying each adjacent coefficient equation
is exact group replay because `[x]G + [d]G = [x+d]G` in the prime-order P-256
subgroup.  The three primitive landmarks replayed 394,374 edges with zero
failures.  Rows with non-coprime edge counts decomposed into shorter cycles and
were reported invalid rather than silently repaired.

| candidate | R | primitive | sign-unique | edges replayed | failures | coefficient bytes |
|:--|--:|:--:|:--:|--:|--:|--:|
| `LD2R1h89f683aa7341` | 1 | yes | yes | 131,458 | 0 | 4,206,656 |
| `LD2R4777h7b483d6354c5` | 4,777 | yes | yes | 131,458 | 0 | 4,206,656 |
| `LD2R6935h980917981827` | 6,935 | yes | yes | 131,458 | 0 | 4,206,656 |

The failed family is not assigned an `FB1` identity or emitted as a native
factor-base dump.  Doing so would promote a screened coefficient construction
despite its failed cost, probability, degree and end-to-end gates.  Its curve
identity remains the registered ICV1 slug
`icv1-fp256-t89188191154553853111372247798585809583-f188c491`.

## Complete toy references

Ten exact cells enumerate 5,752,512 signed distinct-column sums.  They replay
196 cycle edges and have zero false positives, false negatives, edge failures
or bound violations.  Their ordered digest is
`e0a23406351e37e0d8a679097958a741e563d48f61226aacc1439e606b76b4a5`.

The registered bound is intentionally loose.  Even in the seven cells where
it saturates at the whole toy group, exact support occupies only 0.252% to
26.070% of the bound; the three unsaturated cells occupy 2.866% to 21.401%.
This validates the obstruction use of the bound but supplies no P-256
coverage estimate.

## Gates and decision

Exact cycle construction, coefficient replay, sign uniqueness, fixed signed
S17 and one known-log class pass.  The time gate fails at 1.016468 times rho
even under full optimistic support credit.  P-256 usable-relation probability
is unproved; structured residual degree remains unknown; relation collection,
per-row cost, materialized storage and a non-generic end-to-end algorithm are
not established.  The candidate is not promoted.

Reject the two-delta cyclic family for rho parity.  This closes the proposed
constant-delta-state repair, not arbitrary factor bases or arbitrary ECDLP
algorithms.  A successor inside fixed signed S17 would need to break the
arithmetic-run support versus delta-entropy tradeoff, rather than merely encode
the same delta state more compactly.

## Reproduction

```bash
cargo test --bin p256_low_delta_factor_base
cargo run --release --bin p256_low_delta_factor_base -- \
  --round24 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_negation_cross_colour_round24_20261006/negation-result.json \
  --round29 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_memoryless_closure_round29_20261006/memoryless-closure-result.json \
  --round30 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_hidden_state_round30_20261006/hidden-state-result.json \
  --out research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_low_delta_round31_20261006/low-delta-result.json
```

Canonical result artifact: `low-delta-result.json`, 19,312 bytes, SHA-256
`9b9ee16c5a24e868d1b9304f08aa80ca39ac23b1b4586ef730534e739172c435`.
