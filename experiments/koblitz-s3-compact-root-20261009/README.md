# Compact full-rank S3 root index

The [new full-rank solver](../../examples/s3_partial_prepared_fullrank_compact_targetspan.rs)
stores each S3 root-table slot in **13 bytes instead of 16**. Each key remains
a 64-bit canonical root. For the frozen N41 C800 and N53 C1200 bases, a
40-bit value holds two 12-bit representative indices and two 8-bit relative
or Frobenius shifts. The Rust packed slot copies its fields before use, and a
compile-time size assertion fixes the 13-byte layout. The candidate rejects
base sizes at or above 4,096 columns so labels cannot be truncated.

| Full-rank base | N41 C800 | N53 C1200 |
| --- | ---: | ---: |
| Root-table slots | 134,217,728 | 536,870,912 |
| Original 16-byte slot payload | 2,147,483,648 B | 8,589,934,592 B |
| Compact 13-byte slot payload | 1,744,830,464 B | 6,979,321,856 B |
| **Structural payload saved** | **402,653,184 B** | **1,610,612,736 B** |

The [N41 candidate](identities/candidates/IC1N41Ckb1fb65600PDP4rootRCsampleLAgaussTDdirectISO0heee01c30874c.json)
and [N53 candidate](identities/candidates/IC1N53Ckb1fb127200PDP4rootRCsampleLAgaussTDdirectISO0ha9eb2dd51ea9.json)
retain their respective exact factor bases, full-rank preparation policies,
target query, and direct recovery. Their manifests record the new source SHA,
40-bit value encoding, 13-byte slot layout, and table payload. The source is
separate from the frozen 16-byte candidate so its prior ledgers remain
reproducible.

## Correctness controls

Six release example tests passed, including deterministic serial/parallel
root-table construction, both-curve S3 group roots, and target-span linear
algebra. On previously examined T025 controls, each compact run reached
full rank, recovered and verified the fixture scalar, and used zero online
fallback attempts. The [control verifier](verify_controls.py) checked that
**every non-timing report field** matched its original 16-byte counterpart
after excluding the compact slot metadata. That includes the complete
prepared relation witness list, first-witness root choice, target relation,
factor-base digest, rank, S3-call count, and recovered scalar on both curves.
The [receipt](local/controls.json) retains raw and binary hashes. The
[N41 reference](local/reference_n41_t025.json.gz),
[N41 compact](local/compact_n41_t025.json.gz), and
[N53 compact](local/compact_n53_t025.json.gz) control reports are retained;
the N53 reference is the T025 expanded row of the existing
[full-rank ledger](../koblitz-s3-target-conjugate-20261009/local/raw_c1200_fullrank_paired32.jsonl.gz).

The saved bytes are the exact root-table payload difference. The T025
controls served correctness and trace equivalence. A host-isolated receipt
is required to promote any CPU wall-time ratio.

## Peak process memory control

The [RSS freeze](RSS_FREEZE.json) fixed the source, binary, point, workers,
and 180-second limit before the compact measurement. Its
[verified compact receipt](local/peak_rss_n53_t025.json) and
[raw report](local/raw_peak_rss_n53_t025.json.gz) measured **17,988,763,648
bytes** of peak child-process RSS through macOS `getrusage` on N53 T025. A
fresh run of the frozen 16-byte reference measured **19,425,591,296 bytes**,
with [receipt](local/reference_rss_repeat.json) and
[raw report](local/reference_rss_repeat_raw.json.gz). Their observed
difference is **1,436,827,648 bytes**. The earlier 16-byte
[receipt](../koblitz-s3-target-conjugate-20261009/local/memory_peak_c1200_t025.json)
measured 26,739,802,112 bytes for the same binary and point. This reference
spread leaves the reproducible **1,610,612,736-byte table payload saving**
as the firm memory result; additional controlled repetitions are needed for
a stable peak-RSS estimate. The two RSS controls are memory and correctness
observations, not a paired online timing panel.
