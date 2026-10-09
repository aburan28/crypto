# N53 S3 factor base expanded to 552 folded columns

Expanding the signed-Frobenius factor base from 396 to 552 columns increased
its usable points from 41,976 to 58,512 and cut target S3 calls by an
exploratory 1.905× geometric ratio. All 32 paired target solves verified.
Every 552-column target row was spanned at prepared rank 551, so no
target-charged fallback relations were needed. The larger root index doubled
its slot payload to 1,073,741,824 bytes. Complete online time improved by
an exploratory 1.392× geometric ratio. Checked [Sage replay](local/sage_replay.json),
with its [runtime receipt](local/sage_runtime_info.json), independently
verified all 32 scalars, both factor bases (100,488 point labels), and 1,470
distinct relation witnesses across two relation systems.

## Exact candidates and workload

The [freeze](FREEZE.json) fixes curve `EC1N53Ckb1h738a572b91f2` over
`GF(2^53)` with modulus `x^53+x^6+x^2+x+1`, curve `y²+xy=x³+1`, subgroup
order 21044858204113, 14 relation and target workers, and public points
T025–T032. Both arms use the same S3 solver binary, source, target points,
relation seed, and two runs per target in alternating ABBA and BAAB order.
The factor-base size and corresponding requested pre-target rank are the
controlled variables. The [candidate manifests](identities/candidates/)
record their exact constructions and algorithm settings.

| Arm | Candidate ID | Usable points | Folded columns | Prepared rank | Index slots |
| --- | --- | ---: | ---: | ---: | ---: |
| 396 | `IC1N53Ckb1fb41976PDP4rootRCsampleLAgaussTDdirectISO0h1e2914fbd8c2` | 41,976 | 396 | 396 | 33,554,432 |
| 552 | `IC1N53Ckb1fb58512PDP4rootRCsampleLAgaussTDdirectISO0h6ee2f33a9114` | 58,512 | 552 | 551 | 67,108,864 |

The [552-column base dump](local/n53_base_c552.json.gz) fixes all 58,512
encoded points, labels, and representatives. Its raw SHA-256 is
`7b6decbbaf34ec391371670f776f82b801872b27e3db09ccad666c2828bc9632`.
Its enumerated-set digest is
`35209ec8dca60f5f6e5735f9b4548159a0932f206214587be352edbd5f2435e4`.
The natural low-`x` scan ends at `x=1099` for this base. The 396-column
base dump and its checked Sage replay are retained in the earlier
[base-size experiment](../koblitz-s3-base-size-20261009/README.md).

## Paired one-target result

The timer begins after reusable factor-base, index, and rank preparation,
then charges target query generation, all target S3 work, relation checking,
descent, and scalar replay. The [raw ledger](local/raw_paired32.jsonl.gz)
retains every run, phase cost, operation count, preparation witness, and
status. The [analysis](local/analysis.json) checks candidate IDs, frozen
points, scalar replay, and exclusive phase sums.

| Measure | 396 columns | 552 columns |
| --- | ---: | ---: |
| Native verified runs | 16/16 | 16/16 |
| Median target S3 calls | 6,260 | 2,770 |
| Median charged target PDP | 2.260 ms | 1.699 ms |
| Median target descent | 0.116 ms | 0.254 ms |
| Median online time | 2.495 ms | 2.073 ms |
| Median reusable index build, supplementary | 3,202 ms | 7,807 ms |
| Median reusable rank collection, supplementary | 1,340 ms | 1,634 ms |
| Verified preparation relations / novel rows | 588 / 396 | 882 / 551 |
| A/A p95 absolute log online ratio | 2.875 | 2.432 |

The paired exploratory online geometric ratio `396 / 552` is **1.392×**,
with target bootstrap 95% interval `[0.587, 3.579]`; four of eight targets
favor the larger base and three exceed 2×. The A/A values exceed the 10%
host-noise gate, and this host has no exclusive-core isolation receipt, so
controlled speedup remains unknown. The larger base improves relation yield
but its index and target descent costs limit the complete online gain. The
next implementation experiment batches target-state inversions across S3
queries while preserving the exact relation order.

The baseline release binary SHA-256 is
`8217647e609942bb1b2876cd1cd05dde4d6a7dd2b6f4d3713421607d2d3a3920`.
The plan SHA-256 is
`41cfd1266845d4feff29ada0a10813da6145ce96319aba1c2ee2239d05944a2b`,
and the uncompressed run ledger SHA-256 is
`d8e0581a89c06de8723250d7f27ce84fe63bc1cfa8cb6d5b19eab29d8e7c1655`.

To reproduce, build the release
`s3_partial_prepared_parallel_targetspan` example using
`CARGO_TARGET_DIR=/private/tmp/s3-prepared-target-20261009 cargo build
--release --offline --example s3_partial_prepared_parallel_targetspan`.
Run `run_paired.py --output <new-directory> --binaries
/private/tmp/s3-prepared-target-20261009/release/examples`, then
`analyze_paired.py --ledger <new-directory>/raw.jsonl --output
<new-directory>/analysis.json`. For independent replay, save
`/Volumes/SSD990/cryptanalysis/sage --runtime-info` before launching
`/Volumes/SSD990/cryptanalysis/sage -python replay_sage.py` with a new
receipt path.
