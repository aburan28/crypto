# P-256 complementary-cover quotient, round 302: result

## Verdict

The complementary elliptic quotient of the registered split degree-two cover
does **not** provide a second P-256 order-`n` relation row.

For

```text
E:  y^2 = x^3 + a*x + b
H:  v^2 = u^6 + a*u^2 + b
```

the exact maps

```text
pi_1(u,v) = (u^2,v)                         in E
pi_Q(u,v) = (u^2,u*v)                       in Q: Y^2=X^4+aX^2+bX
pi_2(u,v) = (b/u^2,b*v/u^3)                 in E': V^2=U^3+aU^2+b^2
```

replay wherever defined over a complete 994-point small-field control (with
the two `u=0` exceptions counted explicitly) and on 64 deterministic P-256
cover points.  The complementary curve is nonsingular.
However, a deterministic point `R in E'(F_p)` satisfies `[n]R != O` under two
independent scalar algorithms.  The Hasse interval for `E'` lies below `2n`,
so `n` cannot divide `#E'(F_p)`.  There is therefore no rational P-256
order-`n` subgroup in the complementary quotient and no second row to credit.

This leaves one order-`n` row per cover event.  The optimistic free-perfect-
oracle boundary for the minimum-width base remains **13.920747 times rho**;
the registered 17-term base `FB1h2f8621cda105` remains **394.425280 times
rho**.  This is a narrow negative result for this split cover, not a theorem
about every cover, Jacobian, or non-homomorphic relation mechanism.  No
unplanted full-depth P-256 relation was attempted because the promotion gates
failed.

## Requirement-to-evidence status

| requested direction / gate | status | Round-302 evidence or gap |
|:--|:--:|:--|
| pursue many independent rows per event | refuted for this construction | the complementary quotient has no rational order-`n` subgroup |
| exact quotient construction | verified | three displayed maps; nonsingular `E'`; complete toy and deterministic P-256 replays |
| exact group replay | verified | independent LSB/MSB computations of `[n]R` agree on the same nonidentity point |
| zero false positives / negatives on checked instances | verified | zero map or equation failures over all 994 toy cover points and 64 P-256 controls |
| two independent P-256 rows per event | failed | exactly one primary order-`n` row survives |
| structured residual degree `<=5` | not established | no promoted second-row relation system; degree is unset |
| complete collection below `2^120` | failed / unset | no surviving collector or recovery pipeline |
| usable relation below `2^103` | failed / unset | no second-row relation oracle |
| storage below `2^50` | failed / unset | no promotable construction |
| complete cost at or below rho | failed | lower boundary remains 13.920747 times rho |
| actual unplanted P-256 relation | correctly not attempted | second-row, degree, cost, and storage gates fail |

## Boundary table

One unit is the result's optimistic ratio to Pollard rho.  Values other than
rho are lower boundaries with a free perfect relation oracle, not complete
measured attacks.

| variant | independent order-`n` rows / cover event | optimistic ratio to rho | correctness / scope |
|:--|--:|--:|:--|
| Pollard rho | n/a | **1.000000** | complete reference |
| minimum-width independent-log `FB1hc72514a2a8d3` | 1 | **13.920747** | free-oracle lower boundary |
| split cover plus complementary quotient, same base | 1 | **13.920747** | exact transport refutation; no second-row credit |
| registered 17-term `FB1h2f8621cda105` | 1 | **394.425280** | unchanged free-oracle comparison |

Class: **negative**.  Neither boundary nor total cost fell.

## Exact controls

### Complete toy field

The preregistered cell is `p=1019`, `a=2`, `b=3`.

| quantity | exact value |
|:--|--:|
| `#E(F_p)` | 1,032 |
| `#E'(F_p)` | 984 |
| affine points on `H` | 994 |
| primary / quartic replays | 994 / 994 |
| complementary replays | 992 |
| explicitly counted `u=0` exceptions | 2 |
| replay failures | **0** |
| primary fibres | 2 of size 1; 496 of size 2 |
| complementary fibres | 496 of size 2 |

The complete control used 16,100 field additions, 21,062 multiplications,
15,502 squarings, 1,984 inversions, 3,055 Legendre tests, and 497 square-root
exponentiations.  The full exponentiation accounting is preserved in the
result JSON.

### Deterministic P-256 map control

Scanning `u=1,...,142` produced 64 accepted cover points.  Every `pi_1`,
`pi_Q`, and `pi_2` equation replayed, with zero exceptional points and zero
failures.  The canonical encoding of all coordinates has SHA-256

```text
9a20adaf3e12b5ea3b616e158f0202940f706d0b82c5abce15c9f127ce86849c
```

This control used 796 field additions, 1,116 multiplications, 1,052
squarings, 128 inversions, 142 Legendre tests, and 64 square-root
exponentiations.

## Order-`n` transport certificate

The first admissible complementary point is selected at `U=0`:

```text
R = (0,
  41058363725152142129326129780047268409114441015993725554835256314039467401291)
```

Both scalar algorithms return

```text
[n]R = (
  76074407470426473796363728253042208915998852166885294523260263405079725826771,
  84895104842189184257970919770094616861457938270131543429056297146488700524075)
       != O.
```

The independently replayed scalar paths used respectively 166 additions plus
256 doublings and 166 additions plus 255 doublings.  Their coordinates agree
exactly.

The certified coarse Hasse upper bound is

```text
115792089210356248762697446949407573530766708149132191122460380523730634276864
```

and is strictly smaller than `2n`.  If `n` divided `#E'(F_p)`, then the only
positive multiple of `n` in the Hasse interval would force `#E'(F_p)=n`, and
every rational point would satisfy `[n]P=O`.  The displayed witness contradicts
that consequence.  Thus `n` does not divide `#E'(F_p)` and the complementary
quotient cannot carry a nonzero rational P-256 order-`n` image.

## Promotion decision and open boundary

Only dependency integrity and exact replay passed.  The two-row, inverse-
recovery, degree, complete-cost, per-relation, and storage gates failed or are
unset.  Discarded probabilistic branches were not counted as exhaustive.

The narrowest supported finding is:

> The complementary elliptic quotient of
> `H:v^2=u^6+a*u^2+b` cannot supply a second rational P-256 order-`n`
> relation row.

The remaining open obligation is a different P-256-specific cover or a
non-homomorphic mechanism whose extra collapsed equation demonstrably remains
in the order-`n` quotient, followed by executable recovery and complete
below-rho cost.  Nothing here excludes such a construction.

## Isolation, preserved failure, and artifact integrity

The first canonical run exited before mathematics because the preregistered
Round-300 hash contained a transcription error.  That receipt is preserved.
The protocol and runner were corrected in a committed pre-run change; the two
subsequent isolated executions succeeded and emitted byte-identical result and
assessment files.

| run | exit | wall | user | peak RSS | contention |
|:--|--:|--:|--:|--:|:--|
| canonical-v1, bad dependency hash | 1 | 0.005535 s | 0.004937 s | 12,160 KiB | none |
| canonical-v2 | 0 | 0.121349 s | 0.120577 s | 12,160 KiB | none |
| independent-v2 | 0 | 0.082303 s | 0.077612 s | 12,160 KiB | none |

Artifact SHA-256 values:

| artifact | bytes | SHA-256 |
|:--|--:|:--|
| `complementary-quotient-result.json` | 11,099 | `a9a12d823cbca6df965d28f8372217b08e65c04326e13e408bba9dac9f1819c1` |
| `transfer-assessment.json` | 4,163 | `7ec39e854d60599657f4c58fc3ee72736816f00eaed95bde5f1cfc567cc5a097` |
| `isolation.jsonl` | 6,753 | `7c6dd0b63b5b9d5b96537295d213bc53798023e96a54259719baeb795f5f00fc` |
| `PROTOCOL.md` | 6,807 | `e2b48da63b06c9e744693050a7c295cee524763b04270e5d4f3fb1e38205c1b2` |

The result's semantic-evidence hash is
`674ee12ea9a4eb2c9d8afc2499b5797563565ed124b23688c8bf5ac92680239d`.
The transfer assessment's is
`2dcf6c82aad7225c024d9fbac45fe2a0bb98d76761954a948749b913f6fa6aa8`.

The transfer workflow was available, but its referenced
`references/methodology.md` and `assets/assessment-template.json` resources
were unavailable in this environment.  The typed graph, obligations,
controls, scope, and weakest open obligation are therefore emitted directly
and that resource limitation is recorded in the assessment JSON.

## Reproduction

```bash
cargo test --bin p256_complementary_quotient
cargo clippy --bin p256_complementary_quotient -- -D warnings
cargo build --release --bin p256_complementary_quotient --bin isolated_bench
target/release/isolated_bench run --wait --cpus 4 \
  --label p256-complementary-quotient-round302-canonical-v2 \
  --out research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_complementary_quotient_round302_20261008/isolation.jsonl -- \
  target/release/p256_complementary_quotient \
  --registry docs/curves/registry.json \
  --covers docs/curves/covers.json \
  --round300 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_cover_fiber_round300_20261007/cover-fiber-result.json \
  --round301 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_algebraic_escapes_round301_20261007/algebraic-escape-result.json \
  --out research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_complementary_quotient_round302_20261008/complementary-quotient-result.json \
  --assessment research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_complementary_quotient_round302_20261008/transfer-assessment.json
```
