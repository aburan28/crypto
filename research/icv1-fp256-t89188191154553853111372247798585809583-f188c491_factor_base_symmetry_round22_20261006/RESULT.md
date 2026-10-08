# P-256 globally structured factor-base screen, round 22: result

Date run: 2026-10-06

No screened factor base passes.  The exact screen processed 20,990,320
small-multiplier images from ten actual P-256 geometric bases and found zero
cross-column equalities.  Every geometric quotient therefore retains all of
its independent factor-base logarithms.  Seven deliberately scalar-structured
controls do collapse to one known logarithm class, but none is a closed finite
orbit and even the ideal target-collision lower bound is 1.3638 times rho
before state maintenance, verification, or replay.  This is a reproducible
negative result and an accounting correction, not an attack improvement.

## One boundary table

The unit is P-256 group-addition equivalents divided by `sqrt(n)`.  `S_cross`
is the exact independent two-colour lower boundary
`sqrt(pi*K/p_disjoint)` for collecting `K` useful cross-stream collisions.
`log2 collect` is the fixed 138,031-row projection; `log2/row` divides that
work by 138,031; `log2 bytes` is the optimistic materialized two-list storage.
The sparse-Wiedemann plus Berlekamp--Massey lower bound is 2^39.20--2^39.21
operations for the geometric rows and is included in the complete projected
lower bound; it is negligible beside 2^137.36 relation collection.  These are
lower boundaries, not measured end-to-end costs.

| variant | kind | columns | exact rank / K | degree evidence | S_cross | / rho | log2 collect | log2 / row | log2 bytes | decision |
|:--|:--|--:|:--|:--|--:|--:|--:|--:|--:|:--|
| Pollard rho | reference | -- | -- | -- | 1.300 | 1.000 | -- | -- | -- | reference |
| Dickson coset 0, `FB1h6b7699a7d171` | geometric | 131,343 | 0 / 131,343 | residual max 4; unsplit unknown | 642.536 | 494.258 | 137.363 | 120.289 | 142.502 | reject |
| Dickson coset 1, `FB1h085b2065fa5b` | geometric | 131,034 | 0 / 131,034 | residual max 4; unsplit unknown | 641.780 | 493.677 | 137.363 | 120.289 | 142.500 | reject |
| Dickson coset 2, `FB1h6837e846b36c` | geometric | 131,136 | 0 / 131,136 | residual max 4; unsplit unknown | 642.030 | 493.869 | 137.363 | 120.289 | 142.501 | reject |
| Dickson coset 3, `FB1h4d26bec094de` | geometric | 130,966 | 0 / 130,966 | residual max 4; unsplit unknown | 641.614 | 493.549 | 137.363 | 120.289 | 142.500 | reject |
| selected Dickson coset 4, `FB1h2f8621cda105` | geometric | 131,458 | 0 / 131,458 | residual max 4; unsplit unknown | 642.817 | 494.475 | 137.363 | 120.289 | 142.503 | reject |
| Dickson coset 5, `FB1hf3f83ccb1b6b` | geometric | 131,201 | 0 / 131,201 | residual max 4; unsplit unknown | 642.189 | 493.991 | 137.363 | 120.289 | 142.501 | reject |
| Dickson coset 6, `FB1h142ab37c877c` | geometric | 131,123 | 0 / 131,123 | residual max 4; unsplit unknown | 641.998 | 493.845 | 137.363 | 120.289 | 142.501 | reject |
| Dickson coset 7, `FB1h0a4989de61a1` | geometric | 130,955 | 0 / 130,955 | residual max 4; unsplit unknown | 641.587 | 493.528 | 137.363 | 120.289 | 142.500 | reject |
| terminal-zero Dickson, `FB1h2ea06bef7f7a` | geometric | 131,239 | 0 / 131,239 | residual max 4; unsplit unknown | 642.282 | 494.063 | 137.363 | 120.289 | 142.501 | reject |
| affine bitbox, `FB1ha8fa8dae90f9` | geometric | 131,440 | 0 / 131,440 | measured max 6; unsplit unknown | 642.773 | 494.441 | 137.363 | 120.289 | 142.502 | reject degree and cost |
| additive interval, `FB1h98da82cd1af8` | known-log encoding | 131,458 | 131,457 / 1 | none | 1.773 | 1.364 | 137.363* | 120.289* | 134.000* | relabelling; not closed |
| geometric scalar 2, `FB1hc86f774b602e` | known-log encoding | 131,458 | 131,457 / 1 | none | 1.773 | 1.364 | 137.363* | 120.289* | 134.000* | relabelling; not closed |
| geometric scalar 3, `FB1h06333173229b` | known-log encoding | 131,458 | 131,457 / 1 | none | 1.773 | 1.364 | 137.363* | 120.289* | 134.000* | relabelling; not closed |
| geometric scalar 5, `FB1h3b326d663cee` | known-log encoding | 131,458 | 131,457 / 1 | none | 1.773 | 1.364 | 137.363* | 120.289* | 134.000* | relabelling; not closed |
| geometric scalar 7, `FB1he64604080984` | known-log encoding | 131,458 | 131,457 / 1 | none | 1.773 | 1.364 | 137.363* | 120.289* | 134.000* | relabelling; not closed |
| geometric scalar 11, `FB1h4822890a7d37` | known-log encoding | 131,458 | 131,457 / 1 | none | 1.773 | 1.364 | 137.363* | 120.289* | 134.000* | relabelling; not closed |
| geometric scalar 65,537, `FB1h0ad5f09fb8de` | known-log encoding | 131,458 | 131,457 / 1 | none | 1.773 | 1.364 | 137.363* | 120.289* | 134.000* | relabelling; not closed |

`*` is a hypothetical demand for 138,031 generic target collisions and a
materialized two-list.  These controls already know every base logarithm, so
they do not require a relation matrix; their actual remaining problem is one
target decomposition, whose lower boundary is `S=1.773`.  Counting the known
labels as an index-calculus gain would simply rename the DLP.

The smallest geometric base has the smallest quotient boundary, 493.528 times
rho, but this 0.19% change from the selected base is only the change in column
count.  It is not transport, and its smaller representation domain cannot be
credited as a better decomposition probability.

## Exact census and structural checks

Every geometric dump was rebuilt in native Rust and matched its frozen FB1
identity.  Each complete arm sorted `[a]P_i` for `1 <= a <= 16` by full affine
P-256 abscissa.  The ten arms used 19,678,425 group additions and 105,974,976
normalization field multiplications, with 96,753,088 peak logical bytes and
zero disk traffic.  They emitted zero cross-column candidates, zero relations,
zero replay failures, and zero inconsistent cycles; every column is isolated.
The round-21 positive orbit and exhaustive prefix controls are the frozen
detector reference, with zero false positives and false negatives.

All seven scalar controls materialized and verified 131,458 actual P-256
columns in the repository's wide factor-base encoding.  Their FB1 names in the
table are derived from the full canonical preimages.  No arm contains zero or
a signed duplicate.  Their adjacent construction equations have rank 131,457,
but the next point after the prefix lies outside every base, so none is a
finite invariant factor base.  The additive interval has the stronger exact
obstruction: every signed 17-term scalar sum lies in a set of fraction at most
`2^-233.908` of the P-256 group.

The exact P-256 j-invariant is
`0x1198954424ebb0f8479de43131caece8ee0a9b13a558c21e0b2f74e3fcd36aa3`,
neither 0 nor 1728.  Hence the base-field degree-one automorphisms are only
identity and negation, and negation is already folded.  The 258-bit Frobenius
discriminant has trial factors 3 and 5 through `2^20` and a 255-bit residue.
That is context only: it does not determine the full endomorphism-ring
conductor or rule out higher-degree CM maps.

## Accounting correction and decision

Round 21 used the one-stream birthday constant `sqrt(pi/2)` for its ideal
memoryless quotient row.  A valid S17 relation is a collision between an
eight-term stream and a target-shifted nine-term stream.  The first
cross-colour collision has constant `sqrt(pi)`, so the selected base's
optimistic memoryless lower boundary is 494.475 times rho, not 349.646 times.
No algorithm improved; this row is classified as accounting.

Every promotion gate fails.  The geometric arms have no rank compression,
the affine arm also exceeds the degree gate, the unsplit S17 regularity remains
unknown, fixed-row collection is `2^137.363`, the optimistic per-row cost is
`2^120.289`, and materialized storage is about `2^142.5` bytes.  The scalar
controls have no algebraic degree credit and fail rho parity even at `K=1`.
No unplanted full-depth P-256 relation was attempted.

## Reproduction and artifacts

```bash
cargo test --bin p256_factor_base_symmetry
cargo run --release --bin p256_factor_base_symmetry -- \
  --round21 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_fb_log_transport_round21_20261006/transport-result.json \
  --round2 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_dickson_terminal_degree_round2_20261004/factor-base-result.json \
  --round1 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_affine_bitbox_degree_round1_20261004/factor-base-result.json \
  --out research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_factor_base_symmetry_round22_20261006/symmetry-result.json
```

| artifact | bytes | SHA-256 | status |
|:--|--:|:--|:--|
| `symmetry-result.json` | 46,395 | `3116c678d257794040c5519c85a8037177e12381e5f948a47cf1335d0cafde16` | canonical |
| `symmetry-result-preprojection-fix.json` | 39,523 | `8539546a9a3498944f9c6a6aa37cabe8848467b67abab65f434894b5f69287eb` | superseded: row projection double-counted rows |
| `symmetry-result-pre-la-fix.json` | 43,662 | `435b7c381689fa3c28665841ca20c91b10e305721e83ea1aad4a8a548148a45d` | superseded: sparse-LA term absent |

After removing only projection blocks, the first draft and canonical artifact
are byte-identical.  The two superseded files are retained as failed reporting
iterations and must not be cited for projections.  A fresh release replay of
the final code was byte-identical to the canonical artifact, including SHA-256
`3116c678d257794040c5519c85a8037177e12381e5f948a47cf1335d0cafde16`.
