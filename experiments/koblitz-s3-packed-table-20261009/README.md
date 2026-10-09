# Packed inline S3 root table at 396 factor-base columns

An inline `(u64 key, u32 state)` slot reduces the N53 root-table payload by
134,217,728 bytes (25%) while preserving exact canonical-root lookup and
first-inserted witnesses. Six release tests passed. One-worker N41 and N53
controls matched the existing table's prepared relations, target relation,
S3-call count, and recovered scalar exactly. The frozen eight-point N53 panel
completed 32/32 verified native solves. Checked [Sage replay](local/sage_replay.json),
with its [runtime receipt](local/sage_runtime_info.json), independently verified
all recovered scalars, 41,976 factor-base point labels, and 588 distinct
relation witnesses.

## Candidate and workload

The [freeze](FREEZE.json) fixes curve `EC1N53Ckb1h738a572b91f2` over
`GF(2^53)` with modulus `x^53+x^6+x^2+x+1`, curve
`y²+xy=x³+1`, subgroup order 21044858204113, and the parent public points
T025–T032. The low-`x` factor base contains 41,976 subgroup-usable points and
396 signed-Frobenius columns. A target-independent preparation batch requests
rank 395 and reaches rank 396. Each online target uses 14 rotating workers,
first verified relation, cooperative cancellation, and worker drain. Reusable
factor-base, index, and rank preparation are outside the one-target online
interval. The timer charges query generation, target S3 attempts, relation
checking, descent, and scalar replay.

| Arm | Candidate ID | Slot layout | Payload, 33,554,432 slots |
| --- | --- | --- | ---: |
| Existing | `IC1N53Ckb1fb41976PDP4rootRCsampleLAgaussTDdirectISO0h1e2914fbd8c2` | `u64` key and `u64` value | 536,870,912 bytes |
| Packed | `IC1N53Ckb1fb41976PDP4rootRCsampleLAgaussTDdirectISO0h401faac6ebc9` | inline `u64` key and `u32` value | 402,653,184 bytes |

The 32-bit state uses ten bits per representative index and six bits per
Frobenius shift. The slot is `repr(C, packed(4))`; copied keys can be
unaligned, and Rust emits portable reads. The exact-key and synthetic
collision tests cover duplicate and miss behavior. The [candidate
manifests](identities/candidates/) fix the full method and implementation.

## Paired result

Each public point was solved twice per arm, in alternating ABBA and BAAB
order. The [raw ledger](local/raw_paired32.jsonl.gz) retains every run,
including failures, phase times, operations, and witnesses. The
[analysis](local/analysis.json) validates identities, target pairing, scalar
replay, and phase sums.

| Measure | Existing | Packed |
| --- | ---: | ---: |
| Native verified runs | 16/16 | 16/16 |
| Median online time | 3.593 ms | 8.285 ms |
| Median charged target PDP | 3.252 ms | 7.359 ms |
| Median target S3 calls | 7,469 | 6,256 |
| Median reusable setup, supplementary | 4,174 ms | 4,360 ms |
| A/A p95 absolute log online ratio | 3.806 | 2.937 |

The exploratory existing/packed online geometric ratio was **0.717×**,
with a target bootstrap 95% interval `[0.447, 1.230]`. Only two of eight
points had a ratio above 1×. The A/A noise is far above the 10% host-noise
gate and the host has no exclusive-core isolation receipt, so the controlled
speedup remains unknown. The measured layout saves table memory but did not
improve the complete target solve on this panel. A membership filter is the
next candidate to test against absent-key root-table probes.

The packed source SHA-256 is
`f22909478353494823b18413bdc0b05178cf0965fe91bf4406d27416bc05b4af`;
its release binary SHA-256 is
`0376f465efbabd8f335645c490fa3a3dc6b0edb4a534503fbe0c0901220dab7e`.
The baseline binary SHA-256 is
`8217647e609942bb1b2876cd1cd05dde4d6a7dd2b6f4d3713421607d2d3a3920`.
The plan SHA-256 is
`aa8fa43da04d83dfadcfd5ed0cc2ba8fd84392aa1de9f5ef88fbf6715f89355a`,
and the uncompressed ledger SHA-256 is
`a8aa468e09aa8b30c7d77d414aeb3d9e8835d5361a7c10a5f3cff7fe383c64c7`.

To rebuild both release examples from the checked source tree, use
`CARGO_TARGET_DIR=/private/tmp/s3-prepared-target-20261009 cargo build
--release --offline --example s3_partial_prepared_parallel_targetspan
--example s3_partial_prepared_packed_table_targetspan`. Then run
`run_paired.py --output <new-directory> --binaries
/private/tmp/s3-prepared-target-20261009/release/examples` and
`analyze_paired.py --ledger <new-directory>/raw.jsonl --output
<new-directory>/analysis.json`. For independent replay, save
`/Volumes/SSD990/cryptanalysis/sage --runtime-info` before launching
`/Volumes/SSD990/cryptanalysis/sage -python replay_sage.py` with a new
receipt path.
