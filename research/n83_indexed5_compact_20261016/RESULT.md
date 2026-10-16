# K1152 compact exact-root table: half the table, higher cold cost locally

The opt-in 2^27-slot exact-root table reduced K1152 table allocation from
**1,280 MiB to 640 MiB** and measured peak RSS from **1.945 GiB to
1.316 GiB** on the local macOS host. All four fresh processes reconstructed
the same 191,232-point factor base and 110,146,029 exact root entries,
attained rank 1152, and passed independent four-point relation and scalar
replay. The compact table's mean producer cold time was **130.648 s**, versus
**109.619 s** for the default layout. Both A/A spreads exceeded the frozen
10% gate, so this panel remains a local diagnostic and does not advance the
compact mode to a Tier 1 speed comparison.

| Arm | Cold process, s | Setup, s | Relation/rank, s | Peak RSS, GiB | Table, MiB |
| --- | ---: | ---: | ---: | ---: | ---: |
| Default A | 103.753 | 45.156 | 58.325 | 1.945 | 1,280 |
| Compact A | 138.846 | 49.702 | 88.853 | 1.316 | 640 |
| Compact B | 122.450 | 52.287 | 69.902 | 1.316 | 640 |
| Default B | 115.485 | 49.072 | 66.105 | 1.945 | 1,280 |

The default A/A spread was 10.70%; the compact A/A spread was 12.55%.
The compact mean's 21.029-second extra cold time includes about 3.880 s
more setup and 17.163 s more relation/rank work. This follows the expected
cost of denser linear probing, but the phase attribution is an inference from
the producer intervals. The largest observed absolute A/A difference was
16.396 s. The predeclared gate requires both spreads at most 10% and a
positive compact advantage above that difference; this panel does not meet
either requirement. The 20-value original-source Tier 1 K-grid remains the
factor-base selection experiment.

The one-point public input SHA-256 was
`2c068d2f659588edcdbaa685fc73bb146e2bfee3be123d377d537658f8021c6d`.
The source SHA-256 was
`e1b143ec8ba80425a1102b65c3e7c0b67a4b4e9a7e0bc116bcaf13e76c5651f9`,
and the macOS release binary SHA-256 was
`ee1338cb57e7b7f8711e297777ef381e361b5437de506ff1c1e5b760ab84e005`.
The 21 focused Rust example tests passed. The exact compact flag, occupancy
guard, per-arm command, source/binary hashes, base dumps, recovered answers,
RSS telemetry and replay are retained in the archive.

The first local launch was confined to the sandbox. The macOS `time -l`
wrapper returned `sysctl kern.clockrate: Operation not permitted` after the
first producer completed, so its original receipt said
`PRODUCER_FAILURE` even though stdout and the target file were written.
This is a telemetry failure at the wrapper. Its raw directory is preserved
alongside the complete rerun; the runner now records
`INFRA_TELEMETRY_FAILURE` for that exact condition. The four-arm rerun used
the same source and public point with working RSS telemetry.

The complete raw archive, including the exact binary and both the failed
wrapper attempt and four-arm rerun, is
`/Volumes/SSD990-2/n83-archive-20261016/raw/n83-indexed5-compact-local-20261016.tar.zst`,
SHA-256
`a8bd3b2a535e7b5c002c38029c7610913858dc204b6623722d3e65fe39a11366`.
`zstd -t` and a full 47-entry tar listing passed. The versioned checker can
replay the extracted rerun:

```sh
python3 research/n83_indexed5_compact_20261016/verify.py \
  /path/to/n83-indexed5-compact-local-escalated-20261016
```

The analysis JSON SHA-256 is
`fc42679066dd10bfee1dd00c55f287393c5ec5ff260c84766f6b5647e87cb19f`.
