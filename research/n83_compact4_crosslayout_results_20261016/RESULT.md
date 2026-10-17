# N83 compact Bloom cross-layout result

The new `compact4` root Bloom profile halves the K3596/root4 Bloom allocation
from 1 GiB to 512 MiB while retaining exact root-table verification. Both
K3596/compact4 repetitions solved and independently replayed all 8,000 public
targets. Their A-minus-B target-online mean was -0.290 ms, with the frozen
99% contiguous-block interval [-0.323, -0.257] ms, inside the preregistered
±1 ms repeatability gate. All six process arms passed that gate, replay, the
15.25 GiB RSS guard, and the zero-throttling and zero-OOM gates.

K3450/root11/wide3 had the lowest descriptive mirrored target-online mean,
27.672 ms, followed by K3500/root4/standard at 27.987 ms and
K3596/root4/compact4 at 28.651 ms. The K3450-minus-K3500 mean was -0.315 ms
with a 99% interval [-1.523, +0.903] ms. K3450-minus-K3596 was -0.978 ms
with a 99% interval [-2.114, +0.120] ms. The frozen selector therefore
records `BOUNDED_UNKNOWN`; neither contrast establishes an online winner at
its 99% threshold.

| Complete pipeline | Usable base points | Folded columns | Bloom | Mirrored online mean | A-minus-B 99% interval | Setup, A/B | Peak RSS, A/B |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| K3450/root11/wide3 | 572,700 | 3,450 | 1 GiB | 27.672 ms | [-0.006, +0.049] ms | 547.880 / 550.487 s | 11.503 / 11.504 GiB |
| K3500/root4/standard | 581,000 | 3,500 | 1 GiB | 27.987 ms | [-0.108, -0.047] ms | 652.648 / 573.207 s | 14.636 / 14.656 GiB |
| K3596/root4/compact4 | 596,936 | 3,596 | 512 MiB | 28.651 ms | [-0.323, -0.257] ms | 609.954 / 610.575 s | 14.459 / 14.459 GiB |

Six fresh processes ran A/B/C/C/B/A on one GitHub-hosted Ubuntu 24.04
allocation exposing four AMD EPYC 7763 CPUs and 16,373,452 KiB of host
memory. Each process built an empty base and index, reached full rank, and
solved the same 8,000 previously unseen public subgroup points. Independent
replay reconstructed all six factor bases and all 48,000 four-point target
relations and scalar answers. Each arm recorded zero new cgroup CPU-throttled
periods and OOM kills. Major faults were 89/12 for K3500, 16/0 for K3596,
and 0/0 for K3450. The complete process wall times, respectively, were
924.122/846.090 s, 889.596/890.599 s, and 817.047/819.552 s. These
process intervals include reusable setup and 8,000 targets; they are
secondary batch measurements. The online intervals above start after
target-independent setup and include target query, point decomposition,
relation check, descent, and producer recovery check.

The curve is `EC1N83Ckb1a1haacb69ee9f37` over polynomial-basis `GF(2^83)`
with defining-polynomial bitmask `9671406556917033397649557`, model
`y² + xy = x³ + x² + 1`, prime subgroup order `8569786107849059`, and
cofactor `1128547018`. The public-point SHA-256 is
`6c5558613d4f1115dc7bdb09b74d569645c0e85328ddad6555a2484317803a08`;
the verifier-only answer sidecar SHA-256 is
`82decefca05b868e2e9753de3f4892e9bcccfab97ed2c8f6cfbdfd871e13d847`.
The workload ID is `7646bf6fc1e2`. The three executed-snapshot candidate
manifests are [K3450/root11/wide3](candidate_K3450_root11_wide3.json),
[K3500/root4/standard](candidate_K3500_root4_standard.json), and
[K3596/root4/compact4](candidate_K3596_root4_compact4.json).

## Frozen evidence

[GitHub Actions run 38848691554](https://github.com/aburan28/crypto/actions/runs/38848691554)
executed commit `4a791be1d7fc53e4d3c5a7072d5a0bafa885de2a`, IC source SHA-256
`f78382ada0ee4eabfa6820d805d0e51598dbed82154b43ea309bcd489305e2aa`,
and producer binary SHA-256
`0365c555ed95121fc6bdcb76227d9223c9f0df4edd08a144ad0f3bc0656513bc`.
Artifact `12028240172` has ZIP SHA-256
`cff3b20efe5189d0ea952cbb0dc2f2e856c71fc7547ac4698ca582241b3f9c2e`,
matching GitHub's published digest. Cloud and separate local replay produced
the same `analysis.json` SHA-256,
`3e3079efb0b10759b6bc1f43af0ee52081eaf772e49d6b6fe104991ce5ac7921`.
The durable raw ZIP, extracted receipts, and local replay are under
`/Volumes/SSD990-2/n83-archive-20261016/raw/n83-compact4-crosslayout-38848691554`.
The compressed six-row measurements and workload, their decoded digests,
every failure/status field, and the frozen selector are retained in this
directory.

The next comparison should measure one-target cold setup and online work
across a smaller-K factor-base grid, because the current 548–653 s setup
intervals dominate one-target cold cost. A separate larger fresh-target
paired panel can resolve the sub-millisecond online gap between K3450 and
K3500. Both comparisons need their own frozen inputs and resource envelopes;
the five-target same-point Tier 1 IC/rho gate remains the promotion step for
a selected online pipeline.
