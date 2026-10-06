# P-256 isogenous-model Dickson factor-base screen, round 28: result

Date run: 2026-10-06

The staged screen found a slightly denser geometric factor base, not a lower
degree system or a cheaper attack.  Isogeny-walk node 854 lifts 131,531 of the
selected Dickson fibre's 262,144 abscissae, 73 more columns than registered
P-256.  Its modeled probability that a random target has at least one signed
distinct-column S17 representation rises from 96.400018% to 96.511732%.
Because every column remains an independent logarithm, the optimistic
cross-colour boundary moves in the wrong direction, from 494.475 to 494.612
times rho.

No full-depth unplanted relation was attempted.

## Certificate and selection ladder

The input is the existing SHA-256-pinned degree-11 walk certificate.  The
binary replays its shape and all order witnesses before using any curve:

| certificate check | result |
|:--|--:|
| continuous degree-11 edges | 4,096 |
| curve models / distinct j-invariants | 4,097 / 4,097 |
| exact-order witnesses replayed | 4,097 |
| continuity / witness failures | 0 / 0 |
| models with j=0 or 1728 | 0 |

The 262,144 selected-fibre abscissae have SHA-256
`80498d55a2e31d9958f93fd3ae06ebaf45fbe9b89c1aeeb521c2f7f5a6448caa`.
The registered selection ladder then evaluated:

| stage | models | abscissae per model | retained | ordered rows SHA-256 |
|:--|--:|--:|--:|:--|
| depth 8 | 4,097 | 256 | 64 | `e67657f7cb23eec39bb50c7d79dc7b01140b3623a09ec4277d70f79b67159eca` |
| depth 12 | 64 | 4,096 | 8 | `2c9e16bd12b26b13643b3d9daa7e4b669b48a3266deafa` |
| depth 18 | 8 selected plus root | 262,144 | best node 854 | exact rows in JSON |

The root had 131 of 256 lifts at depth 8 and did not enter the selected 64;
it was nevertheless evaluated at full depth as the frozen comparison.  The
discarded prefix branches are not counted as exhaustive depth-18 failures.

## Full-depth candidates

The unit in the last column is the exact independent-log two-colour lower
boundary divided by `1.3*sqrt(n)`.  A larger factor base improves the S17
existence probability but increases the number of independent log classes.

| node | columns | modeled S17 success | boundary / rho |
|--:|--:|--:|--:|
| 0, registered P-256 | 131,458 | 96.400018% | 494.475 |
| **854** | **131,531** | **96.511732%** | **494.612** |
| 1,186 | 131,422 | 96.343992% | 494.407 |
| 2,586 | 131,090 | 95.797655% | 493.782 |
| 2,753 | 131,287 | 96.128338% | 494.153 |
| 2,859 | 131,170 | 95.934254% | 493.933 |
| 2,864 | 131,271 | 96.102193% | 494.123 |
| 3,072 | 131,169 | 95.932566% | 493.931 |
| 3,399 | 131,169 | 95.932566% | 493.931 |

For node 854, the fixed 138,031-row relation lower bound is
`2^137.363459` operations, or `2^120.288827` per usable row.  The optimistic
materialised two-list peak is `2^142.502917` bytes.  Its sparse Wiedemann plus
Berlekamp--Massey lower bound is only `2^39.207017` operations and is
negligible beside collection.

The smallest full-depth base in this selected set is node 2,586.  Its
493.782-times-rho boundary is just the expected reduction from fewer columns;
it also has worse S17 existence probability.  It is not a solver improvement.

## Native factor base

The densest full-depth candidate was materialised in the repository's ICV1,
EC1, FB1 and wide-dump formats:

```text
curve       icv1-fp256-t89188191154553853111372247798585809583-8ae9296a
EC1         EC1P256Cfph72dde3c076e4
factor base FB1hb229f878e201
columns     131,531
points      263,062 signed
points SHA  a0e86d0c297298a9aebbaf1a2e54a0e105e52e662709e7d7ea31fd14edeb9fd5
dump bytes  66,167,787
dump SHA    3070c34daa6495cad9e0bc590d53eccc47b7843f2319b32066027e5ea7763321
```

Every square root and curve equation replays, every point is nonidentity,
all columns are unique up to sign, and the wide dump rebuilds the same FB1
identity.  False positives, false negatives, zero-right-hand-side events and
point-lift failures are all zero.  The measured logical native-build peak is
82,945,003 bytes; the failed storage gate concerns the projected relation
collector, not this factor-base file.

The deterministic screen charged 3,932,416 Legendre tests,
1,002,766,080 Legendre squarings, 503,349,248 Legendre multiplications,
7,864,832 curve-right-hand-side multiplications, 131,531 square-root
exponentiations and 4,097 exact-order witness scalar multiplications.

## Degree and transport result

All 4,097 models have generic j-invariant, so their base-field curve
automorphisms remain identity and negation.  More directly, changing `a,b`
does not remove the universal degree-four support of the short-Weierstrass S3
polynomial, including

```text
x1^2*x3^2 - 2*x1*x2*x3^2 + x2^2*x3^2 + x1^2*x2^2.
```

The coefficient changes can affect lower terms, but no completed structured
degree-of-regularity measurement at most five was obtained.  Pulling the
membership equations back through a rational degree-11 isogeny adds map and
denominator-clearing constraints; it is not a lower-degree representation.

An isogeny transports the DLP between curves, but it does not create exact
relations among these geometric factor-base columns.  No independent-log
quotient reduction is credited.  Thus the degree, relation-collection,
per-row, storage, complete-cost and non-generic-transport gates all fail.

## Decision

Reject the screened isogenous-model route.  Node 854 is a reproducible native
factor-base inventory improvement of 73 columns and 0.112 percentage points
of modeled S17 existence probability.  That local gain slightly worsens the
independent-log boundary, leaves the universal input degree unchanged, and
does not approach rho end to end.

## Reproduction

```bash
cargo test --bin p256_isogenous_dickson_factor_base
cargo run --release --bin p256_isogenous_dickson_factor_base -- \
  --certificate research/p256_isogeny_walk_20261004/certificate.jsonl.gz \
  --round27 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_affine_invariant_round27_20261006/affine-invariant-result.json \
  --out research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_isogenous_dickson_round28_20261006/isogenous-dickson-result.json \
  --dump target/p256-isogenous-dickson-round28.factor-base.json
```

Canonical result artifact: `isogenous-dickson-result.json`, 587,117 bytes,
SHA-256
`950d4b7f4997e146ee30fc7c61174070595fa5d344b9a30cbc0ca578ae7e60f4`.
