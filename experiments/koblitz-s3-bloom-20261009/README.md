# Blocked Bloom filter before S3 root lookup

A 16 MiB blocked Bloom filter rejects absent canonical roots before probing
the 536.9 MB N53 root table. It never supplies a witness: every positive
continues through the exact `u64` key lookup. Six release tests passed,
including absent-key and synthetic false-positive checks. One-worker N41 and
N53 controls matched the existing table's prepared relations, target
relation, S3-call count, and recovered scalar exactly. The frozen paired
panel completed 32/32 verified native solves. Checked [Sage replay](local/sage_replay.json),
with its [runtime receipt](local/sage_runtime_info.json), independently
verified all scalars, 41,976 factor-base point labels, and 588 distinct
relation witnesses.

## Exact candidates and clock

The [freeze](FREEZE.json) fixes curve `EC1N53Ckb1h738a572b91f2` over
`GF(2^53)` with modulus `x^53+x^6+x^2+x+1`, curve `y²+xy=x³+1`, subgroup
order 21044858204113, and public targets T025–T032. The natural low-`x`
factor base has 41,976 subgroup-usable points and 396 signed-Frobenius
columns. A target-independent preparation batch requests rank 395 and reaches
rank 396. Each target starts 14 rotating workers and charges all completed
work through cooperative cancellation and worker drain. The online interval
charges target query generation, all S3 attempts, relation checking,
descent, and scalar replay. Factor-base, root-index, Bloom-filter, and rank
preparation are reusable setup outside this interval.

| Arm | Candidate ID | Root-table payload | Filter |
| --- | --- | ---: | ---: |
| Existing | `IC1N53Ckb1fb41976PDP4rootRCsampleLAgaussTDdirectISO0h1e2914fbd8c2` | 536,870,912 bytes | none |
| Bloom | `IC1N53Ckb1fb41976PDP4rootRCsampleLAgaussTDdirectISO0h899a5cd48de6` | 536,870,912 bytes | 16,777,216 bytes |

The filter has 64 independent shards, each with 4,096 blocks of 64 bytes.
Five bits in one block are selected by a separate 64-bit mixer; the final
root-table lookup still compares the full key and retains the first inserted
witness. The [candidate manifests](identities/candidates/) fix this policy,
source hashes, factor base, and run parameters.

## Paired result

Each target was solved twice per arm in alternating ABBA and BAAB order. The
[raw ledger](local/raw_paired32.jsonl.gz) keeps every run and phase cost; the
[analysis](local/analysis.json) checks candidate IDs, frozen points, scalar
replay, and online phase sums.

| Measure | Existing | Bloom |
| --- | ---: | ---: |
| Native verified runs | 16/16 | 16/16 |
| Median online time | 3.375 ms | 5.872 ms |
| Median charged target PDP | 2.362 ms | 5.619 ms |
| Median target S3 calls | 6,365 | 4,112 |
| Median reusable setup, supplementary | 5,158 ms | 4,993 ms |
| A/A p95 absolute log online ratio | 3.893 | 2.700 |

The exploratory existing/Bloom online geometric ratio was **0.925×**,
with target bootstrap 95% interval `[0.477, 1.741]`. Five of eight points
had a ratio above 1×, and two exceeded 2×. The A/A noise exceeds the 10%
host-noise gate, and this host has no exclusive-core isolation receipt;
controlled speedup remains unknown. The filter preserved the algorithmic
trace but did not improve the complete target solve in this panel. The next
experiment increases the factor base within the following fixed root-table
allocation step to raise relation yield.

The Bloom source SHA-256 is
`4123d9db02db3945ffffcef82d36723a0d8879473181331da2431c17a8de91d5`;
its release binary SHA-256 is
`e3ec5f728b7d896201f0989d6155b738fae42e2d757be7f531e7e81671f523bf`.
The baseline release binary SHA-256 is
`8217647e609942bb1b2876cd1cd05dde4d6a7dd2b6f4d3713421607d2d3a3920`.
The plan SHA-256 is
`357cd3cd4f7a79c248e4bb022c84cfdad5872ac7d18080ecabb479e8e4ab9e87`.
The uncompressed ledger SHA-256 is
`a1bad2a83da5e14945187344f88245b6072e47d39cb357741d9e8809ed56bb56`.

To rebuild the two release examples, use
`CARGO_TARGET_DIR=/private/tmp/s3-prepared-target-20261009 cargo build
--release --offline --example s3_partial_prepared_parallel_targetspan
--example s3_partial_prepared_bloom_targetspan`. Then run
`run_paired.py --output <new-directory> --binaries
/private/tmp/s3-prepared-target-20261009/release/examples`, followed by
`analyze_paired.py --ledger <new-directory>/raw.jsonl --output
<new-directory>/analysis.json`. For independent replay, save
`/Volumes/SSD990/cryptanalysis/sage --runtime-info` before launching
`/Volumes/SSD990/cryptanalysis/sage -python replay_sage.py` with a new
receipt path.
