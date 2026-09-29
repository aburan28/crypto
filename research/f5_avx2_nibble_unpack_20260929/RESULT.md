# AVX2 nibble compaction: exact output, complete-call regression

The paired experiment in [PROTOCOL.md](PROTOCOL.md) completed in
[CI run 36613705628](https://github.com/aburan28/crypto/actions/runs/36613705628)
from PR head `2f8a4a8b85fe498b7a4998bdc20a1d7e4062e793` (the Actions
checkout reported synthetic merge SHA `fc101107a3c9995f791f3407760236a80e0da70b`).
All 88 processes succeeded. The x86 release-mode F5 tests passed. Every
case and seed matched raw row fingerprints, canonical row-space
fingerprints, rank, output term count, row and column counts, criterion
operations, and reduction word operations. The route checks confirmed
that only the new arm selected AVX2 nibble unpacking.

The runner was an AMD EPYC 9V74, Linux x86-64, Rust 1.98.1, with AVX2
and BMI2, one pinned visible CPU and one Rayon thread. The same binary
ran both arms; this host differs from the EPYC 7763 used for the prior
accepted-path estimate, so only the within-run paired ratios below are
used for this decision. The [complete compressed receipt](runs/36613705628-t1.json.gz)
has SHA-256 `c61621ef7e9465477665223101a7fff2cacaa2ce27e9086b71a6d6e087d52548`;
the uncompressed JSON has SHA-256
`a931cc4c9560b1097c5e8866f8feb3ca89109f0e6ff5eadb2d4a9972e3d4f652`
and 856,353 bytes. It retains every call, host load, affinity, source
and binary hashes, exact outputs, and phase timings. No failure,
timeout, or OOM row was dropped.

## Frozen and holdout complete-call comparison

Each ratio is the median of five paired prior/new ratios. Intervals
are exact 3,125-resample bootstrap 95% intervals over five pairs.
`wall_ms` charges criterion, row build, elimination, and full output
unpacking. A/A is five prior/prior pairs on the same workload. Values
above one favor AVX2 nibble compaction.

| Seed | Complete call prior/new (95% interval) | A/A range | Unpack prior/new (95% interval) |
| --- | ---: | ---: | ---: |
| Frozen `0` | **0.912×** (0.906–0.925) | 0.994–1.002 | 0.826× (0.823–0.845) |
| Holdout `badc0de1` | 0.912× (0.910–0.920) | 0.980–1.013 | 0.827× (0.822–0.836) |
| Holdout `5eed2026` | 0.923× (0.892–0.931) | 0.988–1.045 | 0.831× (0.790–0.837) |
| Holdout `f5c02a28` | 0.915× (0.891–0.926) | 0.903–1.096 | 0.828× (0.795–0.835) |

On the frozen primary, the marginal complete-call medians were
129.195 ms prior and 141.465 ms new; reduction was 62.959 and
63.351 ms, and unpacking was 57.222 and 69.328 ms. All six smaller
cases regressed on every seed; their complete-call paired medians
range from 0.886× to 0.964×. Their individual A/A ranges and all
phase samples are in the receipt.

The frozen complete-call ratio and its lower bound miss both the
requested further 2× target and the preregistered 1.05× incremental
gate. The option also violates the smaller-case guard. The AVX2 path
is therefore rejected and removed from production. Its exact tested
[source patch](rejected_candidate.patch),
[workflow](WORKFLOW.yml), and [paired script](paired_f5.py) remain for
replay. The separate clippy diagnostic on this experimental source
was two unnecessary casts in the shuffle-index table; it has no
bearing on the measured algorithm and is absent with the rejected
source. This is a matrix-F5 solver-stage result, not an IC online-time
or DLP speedup. The further 2× complete-call goal remains open.
