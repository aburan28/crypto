# Dense full-rank S3 root table

The [dense solver](../../examples/s3_partial_prepared_fullrank_dense_targetspan.rs)
halves the root-table slot count on both frozen full-rank bases while
preserving every verified S3 relation and recovered scalar. Each shard now
allocates `next_power_of_two(input_roots + 1)` slots. At most one distinct
key can be inserted per input root, so every shard retains an empty slot and
unsuccessful linear probes terminate. Equal roots still choose the first
state/root witness in the original enumeration order. Slots retain the
13-byte key/value encoding of the [compact candidate](../koblitz-s3-compact-root-20261009/README.md).

| Exact root-table payload | N41 C800 | N53 C1200 |
| --- | ---: | ---: |
| Original 16-byte layout | 2,147,483,648 B | 8,589,934,592 B |
| 13-byte compact layout | 1,744,830,464 B | 6,979,321,856 B |
| **13-byte dense layout** | **872,415,232 B** | **3,489,660,928 B** |
| **Saved vs original** | **1,275,068,416 B** | **5,100,273,664 B** |
| Distinct roots / dense slots | 26,186,260 / 67,108,864 | 76,301,710 / 268,435,456 |

The new [N41 candidate](identities/candidates/IC1N41Ckb1fb65600PDP4rootRCsampleLAgaussTDdirectISO0hab6cf6576605.json)
and [N53 candidate](identities/candidates/IC1N53Ckb1fb127200PDP4rootRCsampleLAgaussTDdirectISO0hebadb4a47304.json)
retain the exact factor bases, full-rank preparation, target query, and
direct scalar recovery. Their canonical manifests record the source hash,
13-byte slot, 40-bit value, exact slot count, and capacity policy.

## Verified controls

Seven release example tests passed, including a near-full shard with 15
distinct input roots in 16 slots and a terminating missing-key probe. On
N41 and N53 T025, each dense run reached full rank, recovered and verified
the fixture scalar, and used zero online fallback. The
[control verifier](verify_controls.py) compared every non-timing report field
with the frozen 13-byte candidate, excluding only the capacity policy and
slot counts. The complete prepared relation witnesses, target relation,
first-witness choice, S3-call counts, ranks, and scalars matched exactly.
The [receipt](local/controls.json) retains the source, binary, point, and
report hashes; [N41](local/dense_n41_t025.json.gz) and
[N53](local/dense_n53_t025.json.gz) raw reports are preserved.

The payload savings are structural and exact. Peak process memory and paired
one-target timing need separate frozen measurements; T025 was a correctness
control. The local CPU host has no exclusive-partition receipt for a
controlled wall-time claim.

## Peak process memory control

The [RSS freeze](RSS_FREEZE.json) fixed the dense source, binary, point,
workers, and resource cap before the run. The [verified receipt](local/peak_rss_n53_t025.json)
and [raw report](local/raw_peak_rss_n53_t025.json.gz) measured **18,581,438,464
bytes** of peak child-process RSS on N53 T025. The fresh 16-byte reference
[receipt](../koblitz-s3-compact-root-20261009/local/reference_rss_repeat.json)
measured **19,425,591,296 bytes**, an observed difference of **844,152,832
bytes**. The 13-byte compact candidate's separate receipt measured
17,988,763,648 bytes, and an older 16-byte reference run measured
26,739,802,112 bytes. These single-run process peaks vary enough that the
larger exact table-payload saving cannot yet be promoted to a stable
peak-RSS reduction. The RSS runs are memory and correctness controls; they
do not form a paired online timing panel.
