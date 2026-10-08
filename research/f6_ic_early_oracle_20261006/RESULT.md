# F6 early geometric oracle did not save basis work

**Decision: reject for promotion; keep default off.** Calling the exact
F6 node oracle before inherited-basis specialisation preserved every
oracle decision and verified answer, but did not lower F4 basis-read or
word-operation counts on the archived T1/T7 targets. Both T7 paired
complete-online ratios missed the preregistered 1.10× retention gate:
0.993× and 1.051×. The 321 T7 geometric refutations therefore did not
translate into material skipped basis work on these runs.

The [protocol](PROTOCOL.md) was committed as `0f0ffea34` before code or
timing. The opt-in implementation, exact decision-trace and
configuration tests, and paired scripts were committed as `746b93c18`.
The [two IC1 manifests](candidates/), [input and candidate freeze](FREEZE.tsv),
[source hashes](SOURCE_SHA256SUMS) and [binary identity](BUILD_IDENTITY.txt)
were committed as `ab44b3203` before timing. The [raw runs](runs/),
[measurement rows](measurements.jsonl), [paired rows](pairs.jsonl) and
[derivation check](DERIVATION_CHECK.json) are retained.

The baseline candidate is
`IC1N17Ckb1fb62PDP3f6RCsampleLAgaussTDpdpISO0he8f5dccdff53`;
the early-oracle candidate is
`IC1N17Ckb1fb62PDP3f6RCsampleLAgaussTDpdpISO0hd4a2678ff351`.
Both fix the prepared n17 Koblitz curve, 62 actual usable base points,
29 folded columns, archived public T1/T7, imported certified logs,
three summands, degree three, one Rayon thread, and default algorithm
environment. They differ only in the oracle-order flag. The online
interval starts after reusable setup and charges query, all PDP
attempts, relation check, descent and scalar replay. All eight
fresh-process solves completed and independently replayed the archived
scalar. Every paired target had the same attempt outcomes, oracle
calls, geometric refutations and witnesses, reductions and geometric
additions. The five exclusive online phases summed exactly to wall.

| Target | Rep | Baseline online ms | Early online ms | Baseline/early | Baseline F4 build ms | Early F4 build ms |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| T1 | 1 | 2.864 | 2.915 | 0.983× | 2.283 | 2.332 |
| T1 | 2 | 2.719 | 2.770 | 0.982× | 2.236 | 2.285 |
| T7 | 1 | 144.769 | 145.793 | 0.993× | 119.766 | 120.460 |
| T7 | 2 | 147.645 | 140.475 | 1.051× | 121.483 | 116.312 |

T7 had 1,013 oracle calls and 321 refutations in every arm. F4
basis-read count was 691 and word operations 17,167,040 in every T7
arm; T1 likewise held at nine reads and 422,488 word operations.
The small apparent time differences are within this unisolated Mac
experiment's noise, and the optimization did not change the counted
work. The decision-trace unit and prepared F6 geometric-closure
control passed; the default-false `effective_config` control passed.

No CPU speedup is promoted: this host lacked an isolation receipt and
physical CPU model in the sandbox, and peak RSS was not measured. There
was no paired one-target rho solve or n83 ordinary relation. The
remaining measured hot path is `ReducedBasis::specialise_shared` on
surviving nodes, especially row rewrite and reduction. A next pilot
needs to remove work there rather than reordering the oracle.

Rebuild derived rows and checks with
`sh research/f6_ic_early_oracle_20261006/derive.sh`.
