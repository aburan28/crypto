# P-256 indexed S3-image transfer round 8: result

Date run: 2026-10-04

**PASS on both frozen targets at the actual P-256 modulus.**  The indexed
quadratic image and an independent signed-point subtraction oracle returned
identical ordered factor-base column pairs, with zero false negatives and zero
false positives.  The complete selected factor base rebuilt and verified as
`FB1h2f8621cda105`: 131,458 x-coordinate columns and 262,916 signed points.

The planted target returned exactly `(0,1)` and `(1,0)`.  The preregistered
hash-derived public target returned no pair in this factor base, consistently
in both algorithms.  Every emitted planted pair also passed a separate
four-sign group-addition check against `Q` or `-Q`.

## Boundary table

The structural boundary is the rejected 131,458-squared ordered-column scan.
The indexed path performs one quadratic solve for each left column and one
index lookup for each returned root.  Pair trials, field multiplications and
group operations are deliberately separate units; no conversion between them
is implied.

| target | rejected column pairs | quadratic solves | returned roots / lookups | algebra pairs | reference pairs | FN / FP | algebra field muls | reference group ops | exact | class |
|:--|---:|---:|---:|---:|---:|:--|---:|---:|:--:|:--|
| planted-positive | 17,281,205,764 | 131,458 | 262,916 / 262,916 | 2 | 2 | 0 / 0 | 90,443,104 | 262,916 | yes | engineering / transfer validation |
| hash-public | 17,281,205,764 | 131,458 | 262,916 / 262,916 | 0 | 0 | 0 / 0 | 90,443,104 | 262,916 | yes | engineering / transfer validation |

In frontier counts, the quadratic-solve count is exactly 131,458 times smaller
than the rejected pair count.  Charging a solve and each of its two root
lookups as one abstract visit gives 394,374 visits, 43,819.33 times fewer than
the pair frontier.  These are structural count ratios, not an end-to-end
speedup or a conversion of field multiplications to group operations.

Each nondegenerate solve charged exactly 688 P-256 field multiplications: 11
for coefficient construction, 2 for the discriminant, 289 for the fixed
`p = 3 mod 4` square-root exponentiation and check, 384 for Fermat inversion,
and 2 for root construction.  All 262,916 quadratics across the two targets
returned two field roots.  There were no linear or universal degeneracies;
almost all roots simply missed the factor-base x index.

## Frozen identities and targets

- curve: `icv1-fp256-t89188191154553853111372247798585809583-f188c491`;
- factor-base spec:
  `dickson-torus:depth=18,root_exponent=0x2b6fdc73dc04e7667129`;
- FB1 SHA-256:
  `2f8621cda105a51f8703dabe38dec87f3715ab855770acd7b2af0a14906f9e42`;
- point-set SHA-256:
  `70ab21cbda76205ff43ab6bf85c68607be9f6b89df6399e6fa98e3489b888cf1`;
- terminal:
  `0x5b17195299a3158b93389ad04c776fff2a8bb23ca8659b5b0c2b75f9009b65e5`;
- planted target:
  `x=0x202f1d492dfdaf6d2b8ca40585f2ac8a812284d03f62fa6d8f37ac91cda07b74`,
  `y=0x6447d09272db8edf7e0b767016c00e8e1fef5c0570168dffd6d7c60a147e6dc3`;
- public hash scalar:
  `0x191dff478ed845a34715294b813d9e3cb29d1e524ea2a0a223794eebe4fef3d3`;
- public hash target:
  `x=0x1d5ed1307fcc472716e7f0e485072d5712f5340c1f526e1f2e56d73b29cc0269`,
  `y=0xb1db008d254fa19cf203ba7a126ca279b0d4b501fe72b1579de93249a346f9e6`.

The planted and public pair-set digests are respectively
`6101ec16015a127b9102cbc34d949c2a38335e11642cfd15c2f992197efa7722` and
the SHA-256 empty-string digest
`e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855`.

## Reproduction and evidence

```bash
cargo test --bin p256_s3_image_transfer
cargo clippy --bin p256_s3_image_transfer -- -D warnings \
  -A clippy::mismatched_bit_width_type
cargo run --release --bin p256_s3_image_transfer -- \
  --out research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_s3_image_transfer_round8_20261004/transfer-result.json
```

`transfer-result.json` is 3,407 bytes with SHA-256
`c2a5178aea46223fb1e6f82a65df4d273c68aeacdf5c664a827bec50b144a6e7`.
An immediate full rebuild and experiment replay was byte-identical.

The native unit tests compare the derived coefficients with direct `S3`
evaluation at the P-256 modulus and recover constructed quadratic roots plus a
linear case.  The main run independently checks the algebraic pair set against
public-point subtraction over every signed factor-base row.

## Decision and degree consequence

Accept the indexed quadratic image on the exact selected P-256 factor base.
For this final known-target two-summand gate, generic Gröbner/F4 is removed
rather than assigned a lower degree of regularity: the gate is solved directly
as 131,458 univariate quadratics.  Round 7's retained toy verification systems
had solving degree 3, but that number is not a measured degree for a full P-256
summation-polynomial system.

This round does **not** solve `S18`, propagate unknown intermediate
x-coordinates, collect a length-17 relation, recover a logarithm, or establish
an end-to-end improvement over Pollard rho.  The next credible boundary is a
preregistered residual-depth-2 join: propagate two indexed final-S3 images
through one unknown intermediate while checking whether the frontier remains
linear/indexed or becomes quadratic again.
