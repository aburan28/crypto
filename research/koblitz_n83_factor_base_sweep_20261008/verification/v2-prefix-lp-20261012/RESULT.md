# Exact K=1,182 to K=2,048 residual-pool prefix

This is a retained-object identity check for the primary Koblitz F2^83 arm,
not a relation search or a runtime comparison. The two already-published
`public_x_hash`, seed `2026100801` objects were downloaded from their S3 URIs
in [receipt.json](receipt.json); the current S3 length, ETag and modification
metadata are retained in [the small-object head](s3-head-0.json) and
[the large-object head](s3-head-1.json). Their full compressed BLAKE3 hashes were
recomputed with the repository's `blake3` Rust dependency and matched the
content-addressed manifest names. The verifier additionally checks frozen
compressed SHA-256 hashes and sizes, reads both gzip streams through their
checksums, compares the construction headers and every shared point row,
and checks each row's orbit-column label. The independent generic point and
Frobenius replays and S3 round trips for the individual objects remain in
[`two-hour-extension-20261010`](../two-hour-extension-20261010/REPORT.md)
and [`remaining-budget-k2048-20261010`](../remaining-budget-k2048-20261010/RESULT.md).

The 1,182 representatives and all 196,212 point records are an **exact
ordered prefix** of the 2,048-column object. Its suffix has 866 signed
Frobenius orbits and 143,756 point records. The byte-prefix result plus the
individual objects' full distinct-orbit replays makes this a concrete,
disjoint residual-point pool for a K=1,182 small base. The residual point
rows have SHA-256
`0d008d5a81da8ca34ecd9e8606d0ac3e58c6e325dc7329d69b3f568f224c4bf0`.

For five summands, the exact counts of unordered point multisets with
respectively zero, one and two residual points are
`2,423,643,343,511,761,316,264,808`,
`8,878,309,426,598,156,108,089,140`, and
`13,009,394,065,076,279,872,900,344`. Dividing each by the number of
nonidentity subgroup points gives nominal first moments of about 1.0024,
3.6720 and 5.3806. These counts allow repetitions and include duplicate
group sums. They neither predict solver yield nor graph completion, rank,
linear algebra or complete index-calculus time.

The source's exact orbit-domain size formula and the previously built
K=2,048 all-affine capacity receipt give the following **unordered,
five-summand model preflights**. All rows have 4,980 XOR rows. The total
clause calculation subtracts the K=2,048 domain clauses from its measured
1,140,557-clause model to obtain the common 248,087-clause chain, then adds
each mixed domain's exact clause count. These are construction sizes, not
solver time or memory measurements for the mixed arms.

| Domain positions | Variables | Domain clauses | Total clauses |
| --- | ---: | ---: | ---: |
| Five K=1,182 small | 102,818 | 520,090 | 768,177 |
| Four K=1,182 small, one 866-orbit residual | 102,186 | 492,914 | 741,001 |
| Three K=1,182 small, two 866-orbit residual | 101,554 | 465,738 | 713,825 |
| Five K=2,048 full | 111,478 | 892,470 | 1,140,557 |

Reproduce the byte-prefix and support check after placing the two verified
compressed S3 objects locally:

```text
python3 research/koblitz_n83_factor_base_sweep_20261008/v2_prefix_lp_verify.py \
  research/koblitz_n83_factor_base_sweep_20261008/verification/two-hour-extension-20261010/v2-object/manifest.json SMALL_OBJECT.jsonl.gz \
  research/koblitz_n83_factor_base_sweep_20261008/verification/remaining-budget-k2048-20261010/manifest.json LARGE_OBJECT.jsonl.gz NEW_RECEIPT.json
```

The verifier's three small-fixture controls pass and reject one changed
shared point row and one changed representative. Its source SHA-256 is
`3340b0050cf974547d596873fd117860aca71d47421a1001510cc9977551d30d`;
a second run produced a byte-identical receipt. No new N83 SAT trial was
run for this check; the prior 10,800-second local experimental charge is
unchanged. The next empirical gate must build the suffix as an independently
replayed residual domain, then measure unsuccessful partial searches,
one- and two-residual yield, graph completion and rank under a source-pinned
wall and memory cap. The already-planned matched full-smooth K=1,182 versus
K=2,048 five-summand gate should remain separate, so neither component
result is mistaken for a total-runtime winner.
