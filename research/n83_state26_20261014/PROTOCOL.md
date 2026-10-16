# N83 26-byte state factor-base capacity gate

The optional `KIC_STATE26=1` layout stores two 83-bit normal-basis roots in
22 bytes and packs `(left,right,relative)` into four bytes. It preserves the
existing symmetric-S3 state traversal and exact five-byte root index. This
gate tests whether the smaller state array admits more usable factor-base
points under the same 14-GiB child RSS guard.

## Frozen instance and implementation

- Curve: `y^2+xy=x^3+x^2+1` over `F_2^83`, polynomial low terms
  `[0,2,4,7]`; subgroup order `8569786107849059`, cofactor `1128547018`,
  generator `(6614244651398184860292015,8868260410029249391539529)`.
- Public point input: `research/n83_indexed5_20261013/primary_target.jsonl`,
  SHA-256 `21c8dfc7d488909e5c4a42c839fcc01bdefa9821c11ba001626a1d15697240e8`.
  The known-scalar sidecar has SHA-256
  `651eb7d0ef08a1655036746311b19b6989e340ca4752b6d675606055e59c12a1`
  and is used only by the independent checker.
- Producer source: `examples/koblitz_orbit_dlp_fast.rs`, SHA-256
  `0d78b826504daa6f9d05579286241ce4bfa011ea3978b9c6a91d0cfad0000058`.
  `Cargo.lock` SHA-256 is
  `ac5cd7a3459be776180afe1eb2975c76b1a88f84bd1e99aa20475927ec260a90`.
- All arms use symmetric S3, exact five-byte index with the 75-percent
  request, no Bloom filter or sparse canonicalization, rank seed 29, four
  rank threads, one target thread, and a fresh empty process. The only
  controlled changes are K and the 26- versus 28-byte state encoding.

| Order | Arm | K | Predicted full-orbit usable points | State bound | State plus slot payload GiB | Child cap |
| ---: | --- | ---: | ---: | ---: | ---: | ---: |
| 1 | 26-byte | 2940 | 488,040 | 358,710,870 | 13.686 | 2400 s |
| 2 | 26-byte | 2920 | 484,720 | 353,847,060 | 13.568 | 2100 s |
| 3 | 26-byte | 2900 | 481,400 | 349,016,450 | 13.451 | 1800 s |
| 4 | 26-byte | 2543 | 422,138 | 268,375,505 | 11.498 | 900 s |
| 5 | 28-byte | 2543 | 422,138 | 268,375,505 | 12.00 | 900 s |

The payload calculation uses `stateBound(K)=83*K*(K-1)/2+42*K`,
26 or 28 bytes per state, and `2^30` five-byte slots. It excludes allocator
and other process memory. The K2543 symmetric receipt measured about
0.297 GiB above the two arrays. All arms retain sampled RSS, `wait4`
resource use, cgroup counters, wall caps, failures, base dumps, and
target-independent setup. A child above 14 GiB is stopped and recorded as
`STOPPED_RESOURCE_GUARD`; it is not a performance row.

The same-point strong-rho arm uses the archived `fold-L64-DP6-E0` config
SHA-256 `cf76ebb5372acfc1355b88374f875980be66cdffba1a82cb6e4c2f28bcfa723b`,
seed `2026101301`, one OS thread, and a 120-second cap. The independent
verifier replays every completed base, subgroup, four-point relation,
scalar, and five exclusive online phases. The K2543 pair checks base hash,
exact root-key count, and recovered scalar across both state encodings.
The source test checks every inserted canonical key and first witness for
six actual N83 representatives in both encodings, plus 2,047 probes.

This is one instance-scoped Tier 1 capacity execution per arm. Its timing
ratios remain exploratory until five alternating same-point pairs pass the
tiered noise gates. Actual usable point counts and canonical `IC1` manifests
are materialized only for completed, independently replayed arms. The base
selection for a controlled online comparison will use a newly generated
public target after this capacity gate.
