# n83 F6 raw x-key lookup: retained after small-base replay

The [preregistered full-base gate](PROTOCOL.md) and the
[small-base adjudication amendment](AMENDMENT_1.md) both pass. This
branch is stacked on #1416. `F6SignedPairIndex` now hashes affine
`u128` x coordinates directly and stores the point at infinity in a
separate slot. This removes the 32-byte `PointXKey` enum from the
4,108,723-entry lookup without reserving a field value as a sentinel.
The insertion and query order, exact sign check, portable fallback,
and full-group witness replay remain unchanged.

The frozen K0 curve is
`icv1-f2m83-tm6151469093347-debefd74`. The standard dimension-12
cofactor-projected base has **4,054 actual distinct subgroup-usable
points**, 2,027 sign-folded columns, 8,219,485 unordered pairs, and
4,108,723 signed-sum representatives. Both arms queried the same
public T001 point
`(355fb5df7a905f16921eb,5900a390f42d290f1bbe)` and returned
the same exact `no_witness` result. The separate full-base planted
`[0,2,4,6]` control returned the same witness from portable and PMULL
queries and replayed the full group sum.

| Full-base order | Arm | Index build (ms) | Exact query (ms) | Peak RSS (B) | Outcome |
| ---: | --- | ---: | ---: | ---: | --- |
| 1 | #1416 baseline | 7,544.584 | 2,634.295 | 677,314,560 | exact miss |
| 2 | raw x map | 5,103.850 | 2,295.292 | 543,653,888 | exact miss |
| 3 | raw x map | 4,479.343 | 2,117.776 | 543,391,744 | exact miss |
| 4 | #1416 baseline | 5,662.234 | 2,418.608 | 678,051,840 | exact miss |

All six first-panel processes exited zero. The two-run medians are
**6,603.409 versus 4,791.597 ms build** and **2,526.451 versus
2,206.534 ms query**, baseline versus candidate. The query is 12.66%
lower (1.145× exploratory ratio); adjacent query ratios are 1.148×
and 1.142×. Maximum process RSS is 19.82% lower. The candidate build
also improves. Individual outputs, timings and empty stderr files are
preserved beside the [first-panel status table](status.tsv) and
[runner](run.sh). There was no build between arms.

The first small-base process had an anomalous dimension-8 query median:
1.640 ms baseline versus 15.913 ms candidate, while dimension 10
improved. Because this could have overturned retention on a smaller
case, the amendment froze a second A/B/B/A replay before it ran. The
four dimension-8 **process medians** were 4.100, 3.306, 1.843 and
4.283 ms in that order. The two-process median is **4.191 ms baseline
versus 2.574 ms candidate**; the candidate is 38.6% lower and passes
the amendment's no-more-than-20%-regression gate. The original anomaly
is retained in [`baseline_small.jsonl`](baseline_small.jsonl) and
[`candidate_small.jsonl`](candidate_small.jsonl), not overwritten.
The complete replay is in the four `small_replay_*.jsonl` files and
[`small_replay_status.tsv`](small_replay_status.tsv). All replayed
dimension-8 and dimension-10 outputs and representative counts match.
This replay does not make the unisolated timing a controlled speedup.

The [planted control](planted.jsonl) exited zero. All nine focused
geometry tests passed on physical ARM64, covering ordinary and
exceptional additions, identity handling, exact stored-coordinate
PMULL keys, and exhaustive small-base four-sum membership over 41
target scalars. The first SSD-backed test build failed during LLVM
output with **no space left on device**; that raw
[`geometry_tests_no_space.log`](geometry_tests_no_space.log) is preserved
as an infrastructure failure, not a test failure. The same source and
offline release flags were retried with `CARGO_TARGET_DIR` under
`/private/tmp`; [`geometry_tests_retry.log`](geometry_tests_retry.log)
records 9/9 passes. The original release [build log](build.log) and
the [SHA-256 receipt](SHA256SUMS) preserve provenance.

Baseline #1416 full and small binary SHA-256 values are
`4b5a292164aec2e4584dec1f2470f694f78b71f510eabf47782e3f33afe093d7`
and `76400560b6481206d5e6934c7ce7fa7faf4ebc51bbe03f7b2e65f6715cff1ab8`.
Candidate full, small and planted binary hashes are
`75009df1bfe412f00390b1d7e707ba5d4cb64d196a143f7888e6bfe235a15d95`,
`8d2bed511af624b35e5bd48d16b9b42c1920ddaea6ef6163e0c01867a2b6a0e4`,
and `4dda20b3331d799036c48bc8362f78ff773c3b78aadfb48c090e48a4ab7987b6`.
Final F6 geometry source SHA-256:
`feb569917be8b526eee743e988bc01850c3bba32a312f272bd50e42313c3df25`.
The small-replay runner and amendment hashes are
`7c4f8ef43a0a037e380aaf7782c7445471b275d48d96a8ba13c01521a399fd5b`
and `4c974adf9614f7827952a2071f5c8cacfa1636f3ba5d7467601d1231d1a147b5`.

## Direct #1399-to-final gate

The [second amendment](AMENDMENT_2.md) froze a direct matched comparison
so the separate #1399→#1416 and #1416→final ratios would not be
multiplied across sessions. The unchanged #1399 full binary and final
full binary ran baseline, final, final, baseline on the same public T001.
All four processes exited zero, produced 4,108,723 signed sums, and
returned the same exact miss. The raw [direct panel](direct_status.tsv),
its four `direct_*.jsonl` files, empty stderr files and
[runner](direct_panel.sh) are preserved.

| Direct order | Arm | Index build (ms) | Exact query (ms) | Peak RSS (B) |
| ---: | --- | ---: | ---: | ---: |
| 1 | #1399 baseline | 15,123.246 | 4,263.074 | 1,214,693,376 |
| 2 | final raw x map | 4,314.180 | 2,187.778 | 544,587,776 |
| 3 | final raw x map | 7,616.908 | 3,428.531 | 542,752,768 |
| 4 | #1399 baseline | 14,148.535 | 4,431.044 | 1,268,269,056 |

The direct medians are **4,347.059 versus 2,808.154 ms exact query**
(**1.548×** exploratory baseline/final), and **14,635.890 versus
5,965.544 ms index build** (2.453× exploratory). The adjacent query
ratios are **1.949× and 1.292×**; neither clears the preregistered
requirement that both exceed 2.0. Maximum RSS fell 57.1%. This panel
therefore **does not establish a twofold query speedup**. Its wide
adjacent spread is another reason not to promote a controlled CPU
ratio. The baseline full binary hash was
`8923553ad79d2647383fd8a929b12fd4ddbfe0f880d95b9d641de295fc6e541b`;
the final hash is recorded above. The direct runner hash is
`bfff0f2efc3ada3d2d61afbc5ad06ef2362d3e127f1ab59484f6c65a3bbda046`.

Host: physical Apple M4 Pro, arm64 macOS 26.6, Rust 1.93.1. It had no
auditable exclusive CPU partition and other work ran during the
panels. All wall ratios are **exploratory four-summand component
diagnostics**, not controlled complete-call speedups. Index building
is reusable target-independent preparation; exact query is a
target-dependent component. The K0 base's uniform-target four-summand
coverage ceiling is `4.662e-12`, so these exact misses do not estimate
a useful natural relation yield. There is no complete n83 F6
decomposition, F4/F5/F6 whole-call comparison, recovered target
scalar, paired one-target rho run, or IC online speedup. Those remain
unknown.

Decision: retain the raw-key map after the amended smaller-case gate.
The direct twofold query gate was not met.
The next algorithmic task remains a higher-arity ordinary relation
solver; further four-summand lookup constants alone cannot establish
an F6 or IC crossover.
