# Where F6 inherited-build time goes

**Decision: profile basis construction and specialisation next.** The
support-local Macaulay builder's three timed components cover only
4.75–4.98% of the existing F4 build timer on the eleven-attempt T7 target.
This triggers the preregistered residual rule. The 118.1–121.7 ms residual
includes root echelon reduction and inherited-basis specialisation; it is
not yet attributed between them. Another row, column, or packing change
has a small ceiling on this workload.

The [protocol](PROTOCOL.md) was committed as `1ee8e8e9f` before
instrumentation. The opt-in timers, exactness test, F6 closure control,
worker serialization control, and run scripts were committed as
`5cacf451a`. The [candidate and input freeze](FREEZE.tsv),
[binary identity](BUILD_IDENTITY.txt), [source hashes](SOURCE_SHA256SUMS),
and [candidate manifest](candidates/f6_ic_profile.json) were committed as
`f0e5f46ea` before timing. The [raw runs](runs/),
[measurement rows](measurements.jsonl), and
[derivation check](DERIVATION_CHECK.json) are retained. All four runs
completed, independently replayed the recovered scalar, and matched
attempts, reductions, geometric additions, matrix counts, and word
operations within each target. The five exclusive online phases sum
exactly to each online wall interval.

The frozen candidate is
`IC1N17Ckb1fb62PDP3f6RCsampleLAgaussTDpdpISO0h9700c05ef9e4`.
It uses the prepared n17 Koblitz curve, 62 actual usable base points,
29 folded columns, archived public T1/T7, imported certified logs,
three summands, degree three, the default algorithm environment, and
one Rayon thread. The online interval starts after reusable setup and
charges query, all target PDP attempts, relation check, descent, and
scalar replay. The support-local counters are nested inside F4 build;
they are not additional online phases.

| Target | Rep | Online ms | F4 build ms | Rows ms | Columns ms | Pack ms | Unattributed build ms | Timed share |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| T1 | 1 | 3.415 | 2.693 | 0.370 | 0.274 | 0.083 | 1.965 | 27.03% |
| T1 | 2 | 2.970 | 2.439 | 0.206 | 0.235 | 0.097 | 1.900 | 22.10% |
| T7 | 1 | 149.098 | 123.958 | 2.174 | 2.795 | 0.919 | 118.071 | 4.75% |
| T7 | 2 | 154.784 | 128.090 | 2.333 | 2.991 | 1.050 | 121.716 | 4.98% |

T1 required one support-local root and nine inherited-basis reads;
T7 required eleven support-local roots and 691 inherited-basis reads.
No call delegated to the full-support builder. Each T7 run generated
4,636 root rows and observed 24,200 root columns across its attempts.
The T7 F4 build timer was 82.75–83.14% of complete online wall.
These are two-run observed ranges, not confidence intervals.

The uninstrumented and instrumented builders returned byte-for-byte
equal columns and packed rows in the focused mixed-support unit, including
zero-row and full-support delegation cases. The prepared F6 closure
control passed in its required isolated test process. No algorithmic
speedup is claimed: profiling itself is charged and only the profiled
candidate was timed here. The macOS host was unisolated, the physical CPU
model was unavailable inside the sandbox, and peak RSS was not measured.
No paired one-target rho solve or n83 ordinary relation was obtained.

Rebuild the rows and check with
`sh research/f6_ic_support_local_profile_20261006/derive.sh`.
