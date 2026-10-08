# P-256 batched quadratic inversion round 17: result

Date run: 2026-10-05

**PASS: exact batched inversion cuts the measured fixed-atom oracle by more
than 2x in every frozen P-256 cell.**  The candidate replaces independent
quadratic denominator inversions with one Montgomery batch inversion per image
join or target query.  Coefficients, discriminants, square roots, roots, image
sets, and target lookups are unchanged.

Across the prefix atom and three hash-frozen atoms:

- image-build multiplications fall from 104,576 to 52,010, a 2.01069x
  reduction;
- warm-query multiplications fall from 88,064 to 39,679, a 2.21941x
  reduction;
- cold one-target multiplications fall from 192,640 to 91,689, a 2.10102x
  reduction; and
- two-target work after one build falls from 280,704 to 131,368, a 2.13678x
  reduction.

All four planted targets still hit exactly once.  All four hash-public targets
still miss.  Scalar and batched paths have identical half-image and hit
digests, and fresh exhaustive signed P-256 group references agree with all
eight classifications.  There are zero false negatives, false positives,
linear degeneracies, universal degeneracies, or intermediate-image errors.

The local maximum algebraic degree remains 2.  This is not the degree of
regularity of an unsplit P-256 `S17` or `S18` ideal.

## Boundary table

Every sample has the same counts; the two target rows differ only in their
positive/negative outcome and hit digest.

| sample | target | reference / scalar / batch | scalar muls | batch muls | reduction | solves / roots / lookups | inversion muls | FN / FP | degree | class |
|:--|:--|:--|---:|---:|---:|:--|:--|:--|---:|:--|
| prefix | planted | yes / yes / yes | 192,640 | **91,689** | **2.101x** | 128 / 256 / 256 | 107,520 → 6,569 cold | 0 / 0 | 2 | stage advance |
| prefix | public | no / no / no | 192,640 | **91,689** | **2.101x** | 128 / 256 / 256 | 107,520 → 6,569 cold | 0 / 0 | 2 | stage advance |
| hash-0 | planted | yes / yes / yes | 192,640 | **91,689** | **2.101x** | 128 / 256 / 256 | 107,520 → 6,569 cold | 0 / 0 | 2 | stage advance |
| hash-0 | public | no / no / no | 192,640 | **91,689** | **2.101x** | 128 / 256 / 256 | 107,520 → 6,569 cold | 0 / 0 | 2 | stage advance |
| hash-1 | planted | yes / yes / yes | 192,640 | **91,689** | **2.101x** | 128 / 256 / 256 | 107,520 → 6,569 cold | 0 / 0 | 2 | stage advance |
| hash-1 | public | no / no / no | 192,640 | **91,689** | **2.101x** | 128 / 256 / 256 | 107,520 → 6,569 cold | 0 / 0 | 2 | stage advance |
| hash-2 | planted | yes / yes / yes | 192,640 | **91,689** | **2.101x** | 128 / 256 / 256 | 107,520 → 6,569 cold | 0 / 0 | 2 | stage advance |
| hash-2 | public | no / no / no | 192,640 | **91,689** | **2.101x** | 128 / 256 / 256 | 107,520 → 6,569 cold | 0 / 0 | 2 | stage advance |

The solve count remains 152 for the image build and 128 per target.  This round
reduces work inside those solves; it does not relabel them as fewer solves.

## Where the reduction comes from

The scalar build charges 58,368 multiplication-equivalents to denominator
inversion.  The candidate charges 5,802:

- eight two-leaf joins remain scalar batch-size-one fallbacks;
- four four-leaf joins use batch size 4; and
- two eight-leaf joins use batch size 64.

Each target query batches all 128 denominators.  Its inversion charge falls
from 49,152 multiplications to 767: one 384-multiplication exponentiation plus
`3*128 - 1 = 383` prefix and reverse products.  The other 38,912 warm-query
multiplications are unchanged, giving 39,679 total.

| phase | scalar total | batch total | reduction | scalar inversion | batch inversion |
|:--|---:|---:|---:|---:|---:|
| image build | 104,576 | 52,010 | 2.010690x | 58,368 | 5,802 |
| warm target | 88,064 | 39,679 | 2.219411x | 49,152 | 767 |
| cold target | 192,640 | 91,689 | 2.101015x | 107,520 | 6,569 |
| two targets, one build | 280,704 | 131,368 | 2.136776x | 156,672 | 7,336 |

The accounting charges every prefix product, the product inversion, reverse
recovery product, and root-construction multiplication.  It is not an
unpriced “one inversion” abstraction.

## Exact state and factor-base identity

Round 16's 36-byte packets and 8,192-byte transient reconstruction boundary
are unchanged.  Every packet decoded and re-encoded exactly through:

- curve:
  `icv1-fp256-t89188191154553853111372247798585809583-f188c491`;
- FB1: `FB1h2f8621cda105`;
- FB1 SHA-256:
  `2f8621cda105a51f8703dabe38dec87f3715ab855770acd7b2af0a14906f9e42`;
- point-set SHA-256:
  `70ab21cbda76205ff43ab6bf85c68607be9f6b89df6399e6fa98e3489b888cf1`.

The modeled `m=17` relation success remains 96.400018%.  Batch inversion
changes neither factor-base coverage nor relation yield.

## End-to-end effect

This round establishes a **2.101x deterministic operation-count reduction for
the cold fixed-atom oracle**.  It does not establish a 2.101x end-to-end index
calculus speedup.

If a future complete collector measures fraction `f` of its charged work in
this oracle, the corresponding operation-count ceiling is:

```text
E2E speedup = 1 / ((1 - f) + f / 2.101015)
```

`f` remains null because the repository does not yet contain a full P-256
length-17 relation collector.  Branch-cell enumeration, relation acquisition,
matrix construction, linear algebra, and logarithm recovery are unmeasured and
must not be assigned zero cost.  Accordingly the end-to-end `S` and rho ratio
remain null.

## Reproduction and evidence

```bash
cargo test --bin p256_s3_image_transfer
cargo clippy --bin p256_s3_image_transfer -- -D warnings \
  -A clippy::mismatched_bit_width_type
cargo run --release --bin p256_s3_image_transfer -- \
  --batch-round16 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_s17_packed_atom_round16_20261005/packed-result.json \
  --out research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_s17_batch_inversion_round17_20261005/batch-result.json
```

`batch-result.json` is 25,857 bytes with SHA-256
`57cf49ffbab65356b4ef69848e7824f0377bc058ddd57fe0bf6a8f7b77089ce7`.
The input round-16 result matched its required SHA-256
`af2874ef5f0bdabec75b9bf707067232d0d6e3985f0cce4eba6834c2b5790e44`.

## Decision and next experiment

Accept batched inversion as the default candidate for reconstructed fixed-atom
image builds and target queries.  It passes both preregistered 2x arithmetic
gates without changing state, degree, roots, or classifications.

The next admissible end-to-end step is a charged outer relation scan that uses
the batched oracle and measures what fraction of total work it occupies.  That
scan must retain branch enumeration, negative trials, relation verification,
matrix-row emission, and setup costs.  Until it exists, no P-256 relation,
logarithm, end-to-end speedup, or rho win is claimed.
