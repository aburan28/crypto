# P-256 Dickson deep-branch round 5: result

Date run: 2026-10-04

**PASS on the frozen small-prime boundary.**  Fixing three or four terminal
Dickson levels lowered the maximum certified component solving degree from 4
to 3 for both tested fibres.  All 1,920 component systems completed and agreed
with exhaustive signed-point addition, including 1,900 negative components.

This is the first measured regularity reduction in the sequence.  It comes
from decomposing the factor base into shallow triangular components, not from
changing the terminal constant or merely rewriting quartics as quadratics.

| terminal | branch depth | components/target | targets | positive | negative | correct | complete | max degree | median degree | max columns | total field ops | ops / unbranched | classification |
|---:|---:|---:|---:|---:|---:|:--:|:--:|---:|---:|---:|---:|---:|:--|
| 0 | 3 | 64 | 3 | 6 | 186 | yes | yes | 3 | 3 | 44 | 877,767 | 0.0158 | degree advance |
| 0 | 4 | 256 | 3 | 6 | 762 | yes | yes | 3 | 3 | 24 | 605,746 | 0.0109 | degree advance |
| 369 | 3 | 64 | 3 | 4 | 188 | yes | yes | 3 | 3 | 44 | 855,143 | 0.0163 | degree advance |
| 369 | 4 | 256 | 3 | 4 | 764 | yes | yes | 3 | 3 | 24 | 605,814 | 0.0116 | degree advance |

The fixed unbranched references were degree 4 with 55,498,248 field
operations for terminal zero and 52,441,153 for terminal 369.  Every deep
branch row therefore has degree ratio `3/4 = 0.75`.  The operation ratios are
stage diagnostics on these small instances; they are not an ECDLP speedup.

## What transfers—and what does not

The result establishes a mechanism: when the residual Dickson chain is at
most two levels deep, the direct quadratic `S3` components decide at degree 3
under grevlex.  It does **not** measure the P-256 system.

Applying the same residual-depth rule to the depth-18 P-256 factor base would
fix about 16 levels.  That creates `2^16` components per summand and
`4^16 = 2^32` component pairs per two-summand target before any longer-relation
effects.  The lower component degree must therefore be balanced against an
exponential branch frontier.  No practical P-256 attack or exponent
improvement follows from this round alone.

The best verified P-256 inventory remains round 2's nonzero Dickson coset
`FB1h2f8621cda105`, with 131,458 columns and 96.400018207451% modeled
`m=17` relation success.  Its actual P-256 Gröbner degree remains unmeasured.

## Reproduction and evidence

```bash
cargo run --release --bin p256_dickson_branched_s3 -- \
  --branch-depths 3,4 \
  --out research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_dickson_deep_branch_round5_20261004/degree-result.json
```

`degree-result.json` is 1,232,501 bytes with SHA-256
`0934de288fbfccd5de656ee978c9e3641bdde36a937a2579626658eb5a4ef870`.

The next credible experiment is a scaling table over depths 5–8 that fixes
enough tail levels to leave residual depth 1, 2, and 3, with branch count and
total operations as the boundary.  That would test whether degree is governed
by residual depth and quantify the branch/degree crossover.  It must precede
any extrapolation to depth 18.
