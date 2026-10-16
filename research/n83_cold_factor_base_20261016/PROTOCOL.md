# N83 one-target cold factor-base geometry screen

This screen measures the complete cold process interval for one fresh public
point as factor-base size falls below the online-focused K3450 region. Every
arm starts with an empty factor base, exact root index, and factor-log system;
the producer receives the same point and never receives its known scalar.
The target and four untouched confirmation points were frozen before this
runner was timed. The primary public-point SHA-256 is
`2c068d2f659588edcdbaa685fc73bb146e2bfee3be123d377d537658f8021c6d`;
the answer-sidecar SHA-256 is
`b71c5a8bfdae603931fa6c03fc6f6faa3b261288cf787ac682bf8a3e3afe262d`.
The four reserve-point and answer SHA-256 values are respectively
`6ab08db41c5e86c96162f538a0c45d53331c2f9f823a66eea9bf0b458dd3554d`
and `b21d3a1d54c701255a079a3f1d88940f46563a2f857303b4e3bb98ac40483f39`.
`prior_manifest.json` fixes the hashes of 23 earlier public-point files and
the generator also checks the concurrent wide3 point domain.

The exact curve is `EC1N83Ckb1a1haacb69ee9f37`, polynomial-basis
`GF(2^83)` with modulus bitmask `9671406556917033397649557`, model
`y²+xy=x³+x²+1`, prime subgroup order `8569786107849059`, cofactor
`1128547018`, and rank seed 29. The executed IC source SHA-256 is
`f78382ada0ee4eabfa6820d805d0e51598dbed82154b43ea309bcd489305e2aa`.
The [8,000-point cross-layout result](../n83_compact4_crosslayout_results_20261016/RESULT.md)
admitted this cold screen with selection digest
`75a8f9f7638839d1d3d37e91576a490597f97fb0dfa5cab8b7526d38176415be`.

| Order | Pipeline | Distinct usable base points | Root index | Child cap |
| ---: | --- | ---: | --- | ---: |
| 1 | K2543/root22/off A | 422,138 | 22-byte states, dynamic five-byte load-75 slots | 900 s |
| 2 | K256/root22/off | 42,496 | same | 900 s |
| 3 | K512/root22/off | 84,992 | same | 900 s |
| 4 | K768/root22/off | 127,488 | same | 900 s |
| 5 | K1024/root22/off | 169,984 | same | 900 s |
| 6 | K1536/root22/off | 254,976 | same | 900 s |
| 7 | K2048/root22/off | 339,968 | same | 900 s |
| 8 | K3450/root11/wide3 | 572,700 | 11-byte states, fixed five-byte slots, 1 GiB Bloom | 900 s |
| 9 | K2543/root22/off B | 422,138 | repeated base and index control | 900 s |

The usable count is exactly `166*K`, before sign and Frobenius folding; the
effective relation-matrix column count is `K`. Root22 requests
`ceil(4*S/3)` exact root entries for `S = 83*K*(K-1)/2 + 42*K` retained
states. Its constructor doubles this request and rounds the slot count to a
power of two. The K3450 root11 control retains the prior fixed 2^30-slot
five-byte table and wide3 Bloom profile. All arms use four rank threads,
one target thread, symmetric S3 states and sparse canonicalization. A
900-second cap is above the previously measured 548–550-second K3450 cold
setup interval. Preserve a timeout as a bounded row, with its measured
work and RSS, rather than assigning it a zero cost.

Run the nine arms in the listed order on one Ubuntu 24.04 allocation with
at least four visible CPUs and 15.5 GiB host memory. Require observed peak
RSS at most 15.25 GiB, zero new cgroup throttle, OOM and OOM-kill events,
and at least 2 GiB free disk at admission. Record exact source and binary
hashes, visible affinity, cgroup limits, major faults, per-arm wall and CPU
time, and all failed attempts. After the IC arms, run the frozen strong-rho
binary once on the same public point with seed `2026101601` and a 120-second
cap. Independently replay every complete factor base, four-point relation,
IC scalar, and rho scalar.

The cold comparison clock is the producer's internal `process_total` from
the first cold setup operation through the verified one-target result. Retain
the outer `wait4` wall separately, together with exclusive `setup_total`,
`rank_stage`, `targets_total`, and the five target-online phases. The
one-target online interval begins after reusable setup and includes query,
point decomposition, relation check, descent and producer recovery check.
The known-scalar fixture construction, process launch and independent replay
are separately recorded. The control A/B cold-process difference must be at
most 10% of its mean before using a descriptive ranking. This is a coarse
cloud-instance screen: rank complete resource-admitted pipelines by cold
process time on the single common point and carry the two lowest into a
separate, already frozen four-point confirmation. Keep the selected
factor-base ID null until that confirmation resolves variability. Report
rho's same-point online interval beside each IC online interval; a promoted
Tier 1 speedup still requires alternating same-point pairs and the tier's
calibration gate. Operation-normalized `S` requires a complete operation
counter and conversion boundary in the confirmation run.
