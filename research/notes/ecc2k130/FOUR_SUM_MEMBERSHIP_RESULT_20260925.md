# Complete four-sum support versus compact extraction: fresh n=37/41 streams

**Decision (Gate 1): support/representation is the next lever.** On three
preregistered, disjoint toy-point streams, the complete elliptic-curve
four-sum oracle and the unchanged compact S3 extractor agree on every target.
There is no observed root-index/extractor miss on these streams. At n=37 R3,
164/512 targets have a four-point sum in the frozen base, compared with the
exact **whole-group necessary counting ceiling** 0.450829; the older training
stream's extractor yield was 175/512. The ceiling counts potential unordered
multisets before collisions and is not the measured support probability. The
fresh 164/512 result is 0.3203125. Replacing a root index cannot recover the
348 targets proved absent from this particular base. Agreement on finite
streams is not a proof that the extractor is globally complete.

| Frozen arm | F | Exact unordered ceiling | Oracle+/extractor+ | Oracle+/extractor− | Both− | Extractor+/oracle− | Distinct four-multiset witnesses |
|---|---:|---:|---:|---:|---:|---:|---:|
| n37 R3, 512 targets | 222 | 0.450829 | 164 | 0 | 348 | 0 | 223 |
| n41 R8, 512 targets | 656 | 0.014164 | 10 | 0 | 502 | 0 | 11 |
| n41 R12, 128 targets | 984 | 0.071490 | 3 | 0 | 125 | 0 | 3 |

The deterministic SHA-256 stream is disjoint from all 4096 #747 training
scalars and three independent point holdouts at each n. If those hash-derived
scalars are *modeled* as independent uniform subgroup draws, descriptive 95%
Wilson intervals for the oracle membership fractions are 0.2814–0.3619,
0.0106–0.0356 and 0.0080–0.0666 respectively. These are model-based
sampling descriptions, not coverage proofs; in particular 10/512 at n41 R8
may exceed its 0.014164 whole-group ceiling by finite-sample fluctuation.
The n41 R12 targets are literally the first 128 of the n41 stream used at R8,
so the two n41 rows are not statistically independent samples.

For each arm the Rust oracle materialized every unordered pair `i≤j`, retaining
doubling, opposite pairs, repeated indices and infinity. It queried `T−S`
for **every unique exact group-law pair sum** `S`, then enumerated every pair
bucket match and deduplicated sorted four-index multisets. Every four-sum has
a partition into two unordered pairs, so a completed miss scan is exhaustive
for the frozen point set. The independent Python replay used separate affine
formulas and a polynomial-Euclid inverse cross-checked against the archived
exponentiation inverse on 128 field elements per n. It rebuilt every pair sum,
checked 136 random/exceptional additions per n, re-added all 237 reported
witnesses, independently scanned 16 SHA-selected misses per arm, re-lifted
every compact extractor relation and checked target scalar/group membership.
All controls passed. The root-x collision flag is only an endpoint-x
multiplicity proxy; no missed positive case was available to classify causally.

| Arm | Pair entries | Unique pair sums | Pair collisions | Infinity pairs | Raw matched pair-partition products | Duplicate partition products |
|---|---:|---:|---:|---:|---:|---:|
| n37 R3 | 24,753 | 23,977 | 776 | 111 | 1,314 | 1,091 |
| n41 R8 | 215,496 | 213,201 | 2,295 | 328 | 66 | 55 |
| n41 R12 | 484,620 | 481,177 | 3,443 | 492 | 18 | 15 |

The negative-symmetry closure contributes one infinity pair per `F/2` at
minimum. All 237 positive witnesses admit all three finite balanced pair
partitions; 12 n37 witnesses repeat an index. The 60/2/0 witness root-x
collision-proxy flags at n37/R8/R12 did not cause a missed target. These
figures explain why a counting ceiling can substantially overstate actual
support and why tuple-partition multiplicity must not be mistaken for distinct
target yield.

All costs below are **stage-only** whole-process receipts on local macOS. They
include each producer's own setup and every frozen target but no relation
matrix, scalar recovery, cofactor accounting, or matched rho run. They do not
establish an index-calculus speed advantage.

| Arm | Oracle wall / peak RSS | Compact extractor wall / peak RSS | Exact completion |
|---|---:|---:|---|
| n37 R3 | 0.868 s / 10.6 MB | 15.860 s / 15.7 MB | 512/512 each |
| n41 R8 | 5.823 s / 35.3 MB | 203.127 s / 18.3 MB | 512/512 each |
| n41 R12 | 2.451 s / 117.8 MB | 134.858 s / 17.8 MB | 128/128 each |

The n37 extractor's first 15.755 s receipt overlapped a separate local
7.56 s input-generation task (10:19:42–10:19:50 UTC); it was preserved and
rerun cleanly above. An n41 R12 input-plumbing error made the first extractor
read all 512 n41 points despite the preregistered count 128. Its complete
518.152 s receipt and 512-target raw output are archived as an **invalid
attempt for this comparison**; none of its extra 384 outcomes enters the
128-target inference. The literal corrected 128-line file is SHA-256
`1b154bdd8aa9dabdb37d2dd5a7bfee69a1fa8284ac5db678ef73eaf2c2f297a5`
and is byte-for-byte the first 128 lines of the frozen 512-line n41 stream.
Both R12 producers were rerun against that same file after the amendment
commit. No valid arm hit the preregistered 900 s/2 GiB sampled-RSS stop.

The [pre-outcome protocol](four_sum_membership_20260925/PROTOCOL.md) was first
committed in draft PR #757 at `b9b4020c437864c46b0394f8f7fdc3e059cdff7c`.
The source and target stream were frozen at `12d2e489de77dc633f02df47de2debf865a3f301`,
format-only oracle correction at `089052161aaf8a3c37fc7346418596ac1e650552`,
and literal R12 prefix correction before its valid rerun at
`e50a9d3d26caca3855eaf472cf98af1eafd157a4`. The final oracle Rust
source SHA-256 is `6185018d17e9541245456ffbfec2df36630fe5a4d7e47e8f638ff59f072bdd8a`;
the unchanged #747 extractor source is
`c2bc8b05087df69bef9593363e9d7c112e843ef16da122da50eb29ab22115f09`.
The complete compressed [raw archive](four_sum_membership_20260925/evidence/archive_manifest.json)
contains all six valid runs, all three superseded/invalid receipts, three
independent validations and source snapshots. The archive-manifest SHA-256 is
`c2615eb2adcb1a0a737cd8f896a33d791b00db23dce4383320225729d3ec3ea7`;
the [SHA256SUMS ledger](four_sum_membership_20260925/evidence/SHA256SUMS)
hash is `eee6dec54fccb590091abffff954d42e349becae94e68288a73729c5548bd840`.
`replay_all.py` checks every archive entry, source/input manifest, witness,
complete pair table and selected miss.

The next bounded direction is to change **support**, not merely tune the root
index: test factor-base representations or a higher summand count under a
frozen exact group-law support gate before building another solver. A separate
proposed gate will examine Frobenius-rotated, pairwise disjoint noninvariant
normal-basis subspaces for m=5/6, first proving n131 structural and
cofactor-4/torsion-coset accounting, then measuring tiny-rung support. It is
not an S6/S7 solver commitment. The [n131 materialized four-sum bound](../../sat_factor_base_review_20260908/autolab_orbit_extract_20260924/COMPACT_ORBIT_N131_UNORDERED_ADDENDUM_20260925.md)
still rules out the current explicit-root design at useful coverage/memory;
none of these toy counts changes that bound or demonstrates a full DLP.
