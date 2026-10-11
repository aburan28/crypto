# Exact n=83 field-reduction experiment for the F6 pair kernel

The [initial protocol](PROTOCOL.md) tested a fixed-width reducer for the
exact K_0/K_1 polynomial `z^83 + z^45 + z^2 + z + 1`. The reducer used
three `u128` shift-XOR folds in place of the generic 38-bit chunk path;
all other field polynomials kept the existing path. The candidate passed
the independent bitwise-reduction test, word multiplication and squaring
tests, and all four n=83 F6 pair-closure tests. All benchmark arms used the
same actual K_0 subgroup bases, public T001, Rust 1.93.1 release build, and
an unisolated arm64 macOS 26.6 host. Every result was an exact no-witness
query with the expected pair count.

| Actual usable points | Build baseline / candidate | Query baseline / candidate |
| ---: | ---: | ---: |
| 258 | 18.878 / 11.659 ms (three-run medians) | 16.636 / 11.151 ms (three-run medians) |
| 1,048 | 294.968 / 209.785 ms (three-run medians) | 282.573 / 194.443 ms (three-run medians) |
| 4,054, first run | 5.684 / 7.398 s | 4.510 / 11.414 s |

The first full-base candidate run failed the registered no-more-than-10%
regression gate. Its 8,219,485-pair index and exact query completed, but
the host showed memory pressure. The first [small baseline](baseline_small.jsonl),
[small candidate](candidate_small.jsonl), [full baseline](baseline_full.jsonl),
and [full candidate](candidate_full.jsonl) are retained without substitution.

The separately [registered paired replay](PAIRED_REPLAY_PROTOCOL.md) froze
two different release binaries and ran them in baseline, candidate,
candidate, baseline, baseline, candidate order. Every process exited 0,
returned the same exact no-witness result, and reported about 2.49 GiB
peak RSS. The six [status rows](replay_status.tsv), individual
`replay_<number>_<arm>.jsonl` files, empty stderr files, and the
[runner](replay.sh) preserve the full sequence.

| Full-base statistic, initial plus three replay runs per arm | Baseline | Candidate | Baseline / candidate |
| --- | ---: | ---: | ---: |
| Median index build | 6.026 s | 6.063 s | 0.99× |
| Median exact T001 query | 5.207 s | 4.323 s | 1.20× |
| Median peak process RSS | 2,678,423,552 bytes | 2,678,095,872 bytes | essentially equal |

The replay query median improved, but the full-base build median did not
meet the separately registered 10% improvement threshold. **The candidate
was reverted.** Both the initial failed gate and the follow-up failed gate
remain visible. The wide spread, including a 11.414 s candidate query in
the first run and 2.885 s in the last replay, bars a controlled wall-time
speedup claim. The index build is target-independent while the query is
target-dependent; a verified complete one-target IC timing remains unknown.

The baseline `f2m.rs` SHA-256 was
`15d8f5f357e863e0f16d49f04e88bbbd3269a8ddbb930079f41dff8ef2c0f197`;
the candidate SHA-256 was
`a782b47cf54feeb3e988f2b3c69c99d8fe485fbd09f4fd320d1291f49a557420`.
The complete rejected [source patch](reducer_candidate.patch.gz) can be
decompressed with `gunzip -c`; its uncompressed SHA-256 is
`fdffda412b1302ea3181d983efe212b465149dcaf09793e42b40bf0c4d369f9f`.
The initial and replay protocol SHA-256 values are
`6b6cccb377528fe3f5e14e438ce987aea8645dab62b1f24875e7afadcc1ed385`
and `092111315ec8391473cf8fea6f3c5fbfb9fdc600659d030dccf19bdeb1f1de7c`.
The baseline and candidate replay executable SHA-256 values are
`625950ba3d7dd4ed5b92ba9b1e77aebf0730a8df7c67b2889e5f3866a9e9f713`
and `bc359477c72931fd28f95bd4f796dbccd96c2b60b9af2d0803c393f35935adf3`.
