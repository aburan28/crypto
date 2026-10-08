# P-256 Dickson branched quadratic-S3 round 4: result

Date run: 2026-10-04

Branching one or two terminal levels does **not** lower solving degree: every
one of the 126 component systems still reached degree 4.  The strict degree
hypothesis is falsified.

It does produce a substantial engineering gain.  After two branch levels,
the largest productive matrix fell from 714 to 196 columns and total field
operations over all 16 components per target fell by about 90%.  This total
includes every positive and negative component; no witness branch was selected
with oracle knowledge.

All 126 components completed and matched exhaustive signed-point addition,
including 104 negative components.  There were no timeouts, degree
truncations, or extension-field false positives in these cells.

| terminal | branch depth | components/target | targets | positive | negative | correct | complete | max degree | median degree | max columns | total field ops | ops / unbranched | classification |
|---:|---:|---:|---:|---:|---:|:--:|:--:|---:|---:|---:|---:|---:|:--|
| 0 | 0 | 1 | 3 | 3 | 0 | yes | yes | 4 | 4 | 714 | 55,498,248 | 1.000 | reference |
| 0 | 1 | 4 | 3 | 4 | 8 | yes | yes | 4 | 4 | 360 | 11,492,555 | 0.207 | engineering |
| 0 | 2 | 16 | 3 | 4 | 44 | yes | yes | 4 | 4 | 196 | 5,719,800 | 0.103 | engineering |
| 369 | 0 | 1 | 3 | 3 | 0 | yes | yes | 4 | 4 | 714 | 52,441,153 | 1.000 | reference |
| 369 | 1 | 4 | 3 | 4 | 8 | yes | yes | 4 | 4 | 360 | 10,731,575 | 0.205 | engineering |
| 369 | 2 | 16 | 3 | 4 | 44 | yes | yes | 4 | 4 | 196 | 5,400,695 | 0.103 | engineering |

The branch boundaries were:

- terminal 0: depth 1 `{48,1103}`, depth 2 `{230,240,911,921}`;
- terminal 369: depth 1 `{109,1042}`, depth 2 `{294,548,603,857}`.

## Reproduction and evidence

```bash
cargo run --release --bin p256_dickson_branched_s3 -- \
  --out research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_dickson_branched_s3_round4_20261004/degree-result.json
```

`degree-result.json` is 83,619 bytes with SHA-256
`0ea857544e28c7a640ba7fcc5a11a5010c9475fd844bdd3dce20128483bce5a8`.

The observed cost curve justifies one more bounded rung at branch depths 3
and 4: those components retain only two or one Dickson variables per summand.
That rung must keep charging every component and retain the same correctness
gate.  Until it is measured, the best verified solving degree remains 4 and
the best P-256 factor-base inventory remains `FB1h2f8621cda105`.

No P-256 Gröbner degree or ECDLP speedup is claimed.
