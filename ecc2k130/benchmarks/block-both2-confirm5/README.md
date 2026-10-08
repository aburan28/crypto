# Five-pair confirmation of both2 on block-v3

Status: completed from the preregistered protocol at commit
`d9f575c6bd7ed5748d867df996dfa3fb8dc3f16f`.

The preceding bounded A/B at commit
`8fc284f0d1c334309afab66118c0c2200c15d5ab` measured the
`PACKED_SQUARE_TABLE=1 PACKED_INV_POLY=2` candidate above its block-v3
control in all three long pairs. Its paired ratios were 1.0078054433,
1.0053694525, and 1.0077395917. That clears the follow-up's 1.005 screen
gate narrowly, but its 1.0077395917 median remains below the older roofline
experiment's 1.040 primary-success threshold. This protocol resolves whether
the roughly 0.7% result is stable enough to select as engineering.

## Measured result

The fixed five-pair run passed its promotion gate:

| Pair | Order | Control B/s | Both2 B/s | Ratio |
|---:|---|---:|---:|---:|
| 1 | control, both2 | 5.069008 | 5.095344 | 1.005195 |
| 2 | both2, control | 5.069321 | 5.103107 | 1.006665 |
| 3 | control, both2 | 5.064024 | **5.097573** | 1.006625 |
| 4 | both2, control | 5.069412 | 5.107611 | 1.007535 |
| 5 | control, both2 | 5.067662 | **5.100950** | 1.006569 |

The paired median is **1.0066249686** and every pair is at least 1.005. The
matched control/control ratios were 1.000052, 1.001571, 0.999952, 1.000047
and 0.998972, with median **1.0000465898**; all satisfy the frozen A/A noise
limits. Both arms replayed 300/300 reports with zero drops and produced the
same sorted multiset of 1,480,278 unique table-v3 records, SHA-256
`2fc6f84d4b1969240418cf6f8e779df9538b1d8f01eb6ff4c6a00be7324fae01`.
An independent spread replay matched 300/300 per arm with 299 nonzero trails
per arm. The first 300 persisted records remain zero-step starts, so the
spread check is the evidence that nonzero trails were exercised post hoc.

The selected outcome is `PROMOTE_BLOCK_V3_BOTH2_ENGINEERING`. The older
roofline experiment's 1.040 primary threshold is still not met, so this is
`Partial` under that older protocol. It is a roughly 0.66% same-walk kernel
gain, not an end-to-end ECDLP result. Compact evidence is in
[`result.json`](result.json), [`independent-audit.json`](independent-audit.json)
and [`artifact-manifest.json`](artifact-manifest.json); [`raw-files.sha256`](raw-files.sha256)
binds every extracted file. The 68 MB raw archive
is hash-pinned and retrievable through the Modal identifiers in the manifest;
it is not committed.

## Frozen arms and work

Both arms use the exact v3 point map and the selected B16/T512 block schedule:

```
TABLE_SPLIT_FORWARD=1 TABLE_BATCH_HINTS=1 TABLE_BLOCK_HINTS=1
TABLE_HINT_QUEUE=512 CYCLE_FAST2=1 BATCH=16 THREADS=512 MINBLOCKS=1
```

The control fixes `PACKED_SQUARE_TABLE=0 PACKED_INV_POLY=0`. The candidate
uses `PACKED_SQUARE_TABLE=1 PACKED_INV_POLY=2`. Runtime geometry is fixed at
96,256 threads x 16 slots = 1,540,096 live walks. A measured sample is 1,024
steps x 64 launches = 100,931,731,456 complete updates. Rates are only the
single final `finished:` rate from each raw log.

Before timing, each arm must:

1. build from the same clean, recorded commit and source manifest;
2. identify every common and arm-specific runtime marker;
3. run 96 steps x 6 launches at DP weight 48;
4. replay 300/300 producer reports with zero drops, `MISMATCH`, or `OVERFLOW`;
5. produce the same complete header-aware, sorted table-v3 record multiset.

The binaries, build logs, replay logs, complete DP corpora, source manifest,
and binary manifest are retained. A failed preflight suppresses timing.

## Frozen timing order

Two warmups, control then candidate, are excluded. Ten measured pairs follow.
Every sample uses 64 launches. A/B and A/A panel order alternates by round so
neither panel always occupies the earlier thermal position:

| Round | first pair | second pair |
|---:|---|---|
| 1 | A/B: control, candidate | A/A: aa1, aa2 |
| 2 | A/A: aa2, aa1 | A/B: candidate, control |
| 3 | A/B: control, candidate | A/A: aa1, aa2 |
| 4 | A/A: aa2, aa1 | A/B: candidate, control |
| 5 | A/B: control, candidate | A/A: aa1, aa2 |

`aa1` and `aa2` are labels for the byte-identical control executable. The A/A
ratio is always aa2/aa1, independent of execution order. The A/B ratio is
always candidate/control.

## Frozen decision

All correctness and evidence gates above are mandatory.

- **Promote** the both2 arithmetic knobs into the block-v3 engineering preset
  only when all five A/B ratios are at least 1.005, their paired median is at
  least 1.005, every A/A ratio lies in [0.995, 1.005], and the A/A paired
  median lies in [0.998, 1.002].
- **Reject** and keep the arithmetic knobs off when a correctness/corpus gate
  fails or the five-pair A/B median is at most 1.0.
- Otherwise **retain optional**. A noisy A/A panel cannot qualify a marginal
  A/B result.

The older roofline classification is reported separately and cannot be
weakened: paired median at least 1.040 plus every A/B repetition above 1.0 is
its primary performance condition. A smaller positive ratio is `Partial` in
that protocol.

This is a same-walk GPU kernel benchmark. It does not measure collision time,
natural relation yield, extraction, total ECDLP cost, or a cryptanalytic break.

## Reproduction

The native helpers self-test before either build. The measured job was
dispatched only after this protocol and its harness were committed:

```sh
modal run --detach modal_job.py \
  --job benchmarks/block-both2-confirm5/gpujob.sh \
  --out /tmp/ecc2k-block-both2-confirm5 --gpu RTX-PRO-6000
```
