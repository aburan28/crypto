# C1800 striped S3 target-state schedule

The N53 C1800 striped target schedule preserves scalar recovery and exact
factor-base geometry, but increases S3 decomposition work on all three
checked targets that reach its reordered region. The geometric ratio of the
existing prefix schedule's S3 calls to the striped schedule's calls is
**0.870606** across the four diagnostic points. The prefix schedule remains
the selected C1800 target path.

The implementation changes only target-state order. It visits the first 256
embedded C1200 states in their original point-seeded order, then interleaves
eight distant ranges of the rest of the C1200 prefix while querying the full
C1800 root table. The complete C1800 state traversal remains the fallback.
Two release schedule tests verify that the first 256 states are unchanged and
that the stripe and fallback visits are permutations of their declared state
sets. The T205 correctness control found exactly the prefix build's first
relation witness and used the same 128 S3 calls.

The source-bound candidate is
`IC1N53Ckb1fb190800PDP4rootRCsampleLAgaussTDdirectISO0h57c66f00aeab`.
Its immutable manifest records the exact schedule, 190,800 distinct
subgroup-usable factor-base points, 1,800 relation columns, 1,200 prefix
columns, batch-eight target inversion, source digest, and shared library tree.
The reference is the existing C1800 prefix candidate
`IC1N53Ckb1fb190800PDP4rootRCsampleLAgaussTDdirectISO0h1ab23390e700`.
Both use relation seed `20260928`, 14 preparation workers, and one target
worker. Every checked process reached full rank, recovered the public-point
scalar with native group verification, and passed checked Sage scalar replay.

| Frozen public point | Prefix S3 calls | Stripe S3 calls | Prefix / stripe |
| --- | ---: | ---: | ---: |
| T205 | 128 | 128 | 1.000× |
| T209 | 560 | 672 | 0.833× |
| T211 | 624 | 864 | 0.722× |
| T212 | 672 | 704 | 0.955× |

These points were chosen from an earlier eight-point panel because T205
finishes within the unchanged lead and the other three reach beyond it.
Their call counts test the schedule mechanism; they are not a fresh-target
performance panel. Process wall times are excluded from the comparison.
`CALL_DIAGNOSTIC.json` binds the source commit, release binary, candidate
manifest, previous baseline ledger, workload IDs, report hashes, first
witnesses, and scalar certificates. `local/sage_replay.json` binds the four
independent scalar checks to the checked Sage runtime receipt.
`local/raw.tar.zst` contains each striped report and process stream; its
SHA-256 is
`221e6a6cc80f8519633cf425cf63f62b0a02a0c86ca53deb23418102789c7078`.

The record can be regenerated while the local executable and raw directories
exist:

```sh
python3 experiments/koblitz-s3-c1800-stripe-screen-20261010/record.py \
  --stripe-binary /private/tmp/s3-crosscurve-build/release/examples/s3_partial_prepared_fullrank_prefetch_projective_c1800stripe_targetspan \
  --baseline-dir /private/tmp/s3-batch32-source-matched-screen-20261010 \
  --stripe-dir /private/tmp
```

`record.py` rejects an existing record; use a clean copy of this experiment
directory to regenerate it. The archived prefix and stripe raw reports
preserve the comparison after temporary directories are cleared.
