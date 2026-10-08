# P-256 root-orbit factor-base screen, round 23: result

Date run: 2026-10-06

No screened root-orbit factor base passes.  Sixteen new FB1 bases provide the
intended exact global coordinate symmetry, but the complete `[1..16]P_i`
census finds zero cross-column elliptic-group relations in 34,099,472 images.
The symmetry therefore does not transport factor-base logarithms: every arm
retains one independent logarithm per column.

The smallest arm, `FB1h8fc8b5fd8529`, has the lowest formal cross-colour
boundary at 491.494 times rho.  That is only a 0.41% width effect relative to
round 22's smallest geometric base; it is not structural compression and it
cannot be credited as a better relation mechanism.  The bases admit
quadratic-generator membership chains, but their structured residual and
unsplit S17 degrees of regularity remain unknown.  Generator degree two is
not solving degree.

## Exact screen

The least primitive root of the P-256 base field under the frozen complete
factorisation of `p-1` is 6.  The implementation verified exact subgroup
orders 260,270 and 272,425, then derived eight deterministic nonsingular PGL2
maps per order.  All 4,261,560 source abscissae closed under both
`u -> omega*u` and the induced fractional-linear recurrence.  Quadratic
lifting retained 2,131,217 P-256 columns.

Every base has a canonical `ecbench.factor_base/v1` manifest, a full point-set
digest using `ecbench.factor_base_dump/v1-wide` point-key semantics, and an
`FB1h<12 hex>` identity.  The first 4,096 sorted columns of every arm were
independently regenerated with BigUint field arithmetic: 65,536 point replays,
zero failures.  Point keys were sorted and unique, every retained point was on
P-256, and every subgroup/order/closure check passed.

| candidate | FB1 | columns | exact images | rank / K | S_cross | / rho |
|:--|:--|--:|--:|:--|--:|--:|
| d260270 arm 0 | `FB1h0a763532252c` | 130,158 | 2,082,528 | 0 / 130,158 | 639.632 | 492.025 |
| d260270 arm 1 | `FB1h8688c3fdd265` | 130,404 | 2,086,464 | 0 / 130,404 | 640.236 | 492.489 |
| d260270 arm 2 | `FB1h073bace3e392` | 130,287 | 2,084,592 | 0 / 130,287 | 639.949 | 492.269 |
| d260270 arm 3 | `FB1hca1721dbc049` | 130,469 | 2,087,504 | 0 / 130,469 | 640.396 | 492.612 |
| d260270 arm 4 | `FB1h5d205e2512e5` | 130,393 | 2,086,288 | 0 / 130,393 | 640.209 | 492.469 |
| d260270 arm 5 | `FB1h8fc8b5fd8529` | 129,877 | 2,078,032 | 0 / 129,877 | 638.942 | **491.494** |
| d260270 arm 6 | `FB1hc67754b68b17` | 130,400 | 2,086,400 | 0 / 130,400 | 640.226 | 492.482 |
| d260270 arm 7 | `FB1h4e1d6910749b` | 130,049 | 2,080,784 | 0 / 130,049 | 639.365 | 491.819 |
| d272425 arm 0 | `FB1h437e407515c5` | 136,259 | 2,180,144 | 0 / 136,259 | 654.444 | 503.418 |
| d272425 arm 1 | `FB1h3b9c039e801a` | 136,094 | 2,177,504 | 0 / 136,094 | 654.048 | 503.113 |
| d272425 arm 2 | `FB1heb5cb94a083d` | 135,679 | 2,170,864 | 0 / 135,679 | 653.050 | 502.346 |
| d272425 arm 3 | `FB1hd59ed6581112` | 135,910 | 2,174,560 | 0 / 135,910 | 653.605 | 502.773 |
| d272425 arm 4 | `FB1h74b1984b0856` | 136,529 | 2,184,464 | 0 / 136,529 | 655.091 | 503.916 |
| d272425 arm 5 | `FB1hd3e7d589cfa2` | 136,155 | 2,178,480 | 0 / 136,155 | 654.194 | 503.226 |
| d272425 arm 6 | `FB1h9e3c6143abb6` | 136,254 | 2,180,064 | 0 / 136,254 | 654.432 | 503.409 |
| d272425 arm 7 | `FB1hcbee29703684` | 136,300 | 2,180,800 | 0 / 136,300 | 654.542 | 503.494 |

There were zero cross-column candidates, zero reported group relations, zero
replay failures, and zero inconsistent cycles.  The round-22 exhaustive
detector controls remain the frozen zero-FP/zero-FN reference.

## Degree result and obstruction

For matrix `(a,b;c,e)`, membership away from the rejected pole is

```text
((e*x-b)/(a-c*x))^d = 1.
```

A binary addition chain expresses this with equations of generator degree at
most two.  That is a real improvement in the maximum input degree relative to
writing the degree-`d` equation densely, but it does **not** establish degree
of regularity at most 5.  No structured residual or unsplit S17 Gröbner run was
completed for this family, so both fields remain `null` in the artifact.  The
PGL2 recurrence permutes factor-base abscissae, not elliptic-curve points by a
group homomorphism.  Since P-256 has generic `j`, this nonidentity degree-one
coordinate action is not a P-256 curve automorphism.  The exact rank-zero
census confirms that the coordinate orbit yields no small-multiplier log
transport.

This is the dominant obstruction: a low-degree membership representation is
local algebraic structure, while relation collection needs a global
elliptic-log quotient or a non-generic target mechanism.  This family supplies
neither.

## Counted work and end-to-end projection

Candidate construction counted 4,261,560 orbit multiplications, 12,784,680
transform multiplications, 12,784,152 batch-inverse multiplications, 528 batch
inversions, 153,416,160 square-root multiplications, 1,090,959,360 square-root
squarings, and 12,784,680 closure multiplications.  Transport counted
31,968,255 group additions and 172,115,152 normalization field
multiplications.  Peak logical memory for one transport arm was 100,485,344
bytes; disk traffic was zero.

For the best-width arm, the optimistic independent two-colour lower boundary
is:

- `S_cross = 638.942`, or 491.494 times rho;
- fixed 138,031-row collection `log2 = 137.363` operations;
- `log2 = 120.289` operations per usable row;
- sparse-linear-algebra lower bound `log2 = 39.188` operations;
- materialized two-list storage `log2 = 142.494` bytes.

All promotion gates fail.  In particular, residual degree is unmeasured,
collection exceeds `2^120`, per-row cost exceeds `2^103`, storage exceeds
`2^50`, and the complete boundary is nowhere near rho.  No full-depth
unplanted P-256 relation was attempted.

## Failed iteration and reproduction

The first release attempt stopped before producing an artifact because its
independent checker compared fixed-width field bytes with minimally encoded
BigUint bytes.  The retained source left-pads the independent encoding to 32
bytes and the complete rerun has zero mismatches.  `INVALIDATED.md` records
that stopped run; it receives no experimental credit.

```bash
cargo test --bin p256_root_orbit_factor_base
cargo run --release --bin p256_root_orbit_factor_base -- \
  --round22 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_factor_base_symmetry_round22_20261006/symmetry-result.json \
  --out research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_root_orbit_fb_round23_20261006/root-orbit-result.json
```

Canonical artifact: `root-orbit-result.json`, 77,729 bytes, SHA-256
`f315ef150c123824c4d8954f55e62120579fd666ddf42051bda5351a798bac07`.
