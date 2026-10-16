# Compact refuted F6 bases did not improve target time

**Decision: reject for promotion; keep default off.** Replacing a
refuted inherited basis with its single constant-one column and row
preserved every verified result, but both T7 paired complete-online
ratios missed the preregistered 1.10× gate: 0.956× and 0.980×. T1
ratios were 0.906× and 0.945×, also below its 0.95 floor. The candidate
reduced the T7 cumulative F4 column count from 900,702 to 831,162,
showing that refuted bases were compacted, but F4 word operations and
oracle decisions did not change. Any saved clone work failed to pay for
the compacting work on these targets.

The [protocol](PROTOCOL.md) was committed as `d89bd2258` before code
or timing. The opt-in implementation, two-descendant correctness test,
worker configuration control and paired scripts were committed as
`43d92e849`. The [candidate and input freeze](FREEZE.tsv),
[source hashes](SOURCE_SHA256SUMS), [binary identity](BUILD_IDENTITY.txt)
and [manifests](candidates/) were committed as `5409d4c42` before
timing. The [raw runs](runs/), [measurement rows](measurements.jsonl),
[paired rows](pairs.jsonl), and [derivation check](DERIVATION_CHECK.json)
are retained. This candidate adds the executed `inherited_f4.rs` source
hash explicitly to its immutable manifest.

Baseline candidate:
`IC1N17Ckb1fb62PDP3f6RCsampleLAgaussTDpdpISO0h9628dfd41b76`.
Compact candidate:
`IC1N17Ckb1fb62PDP3f6RCsampleLAgaussTDpdpISO0h3cc84cd498bd`.
They fix the prepared n17 Koblitz curve, 62 actual usable factor-base
points, 29 folded columns, archived T1/T7 public points, imported
certified logs, three summands, degree three, one Rayon thread, and
default algorithm environment. They differ only in the
`compact_refuted` flag. All eight fresh-process target solves completed
and independently replayed the archived scalar. Paired attempts,
per-attempt outcomes, oracle calls/refutations/witnesses, reductions,
and geometric additions matched. The five exclusive online phases
summed exactly to each complete online wall interval.

| Target | Rep | Baseline online ms | Compact online ms | Baseline/compact | Baseline F4 build ms | Compact F4 build ms |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| T1 | 1 | 3.194 | 3.526 | 0.906× | 2.570 | 2.918 |
| T1 | 2 | 4.418 | 4.674 | 0.945× | 3.401 | 3.706 |
| T7 | 1 | 158.196 | 165.540 | 0.956× | 130.375 | 136.423 |
| T7 | 2 | 156.077 | 159.252 | 0.980× | 128.911 | 131.448 |

Both T7 arms had 691 F4 basis reads, 17,167,040 F4 word operations,
1,013 oracle calls and 321 geometric refutations. The compact arm's
69,540 fewer cumulative columns came from refuted states; it did not
remove work on the surviving bases. The focused test proved that a
compact refuted child retained the exact constant-one decisive result
after two more specialisations. The prepared F6 closure control and
default-false `effective_config` serialization control passed.

These CPU ratios are exploratory: the Mac host was unisolated, its
physical CPU model was unavailable inside the sandbox, and peak RSS was
not measured. There was no paired one-target rho solve or n83 ordinary
relation. The measured hot path remains specialisation of non-refuted
bases, especially row rewrite and reduction. No 2× F6 gain or
end-to-end IC speedup is established by this pilot.

Rebuild derived rows and checks with
`sh research/f6_ic_compact_refuted_20261006/derive.sh`.
