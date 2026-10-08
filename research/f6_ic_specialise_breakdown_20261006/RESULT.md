# F6 child specialisation: pivot reduction dominates the timed path

The preregistered diagnostic found pivot insertion/reduction at
**60.8–60.9%** of `ReducedBasis::specialise_shared` on the eleven-attempt
prepared n17 T7 target. Displaced-row rewriting took 23.6%. The same
ordering held on the one-attempt T1 control. The ≥50% pivot-reduction
gate passed, so the next F6 stage candidate should change the insertion
algorithm or its row representation and compare complete calls. This is
cost attribution, not a measured speedup.

The exact candidate was
`IC1N17Ckb1fb62PDP3f6RCsampleLAgaussTDpdpISO0hc0e41ab92a2a`, on
curve `EC1N17Ckb1hbbe2b5b6b1e6`: 62 actual usable factor-base points,
29 folded columns, prepared certified logs, degree-three F6-IC, and the
archived public T1/T7 workloads. The [protocol](PROTOCOL.md) was committed
as `96d09e1a0` before the [opt-in profiler](../../src/cryptanalysis/inherited_f4.rs)
(`ca4076f87`). The [source, binary, inputs and manifests](FREEZE.tsv) were
frozen as `36d3b1ece` before any panel run. The release binary SHA-256 was
`9eba0550b481c6ff97b3e28aa6e077ec4ece034e78e434ad472713ede1269312`.

| Target | Rep | Online ms, instrumented | Specialise ms | Pivot reduction ms / share | Rewrite ms / share | Layout + bookkeeping ms | Counted F4 word ops | Verified |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| T1 | 1 | 2.847 | 1.440 | 0.881 / 61.2% | 0.320 / 22.2% | 0.158 | 422,488 | yes |
| T1 | 2 | 3.243 | 1.513 | 0.898 / 59.4% | 0.336 / 22.2% | 0.194 | 422,488 | yes |
| T7 | 1 | 197.658 | 142.189 | 86.491 / 60.8% | 33.570 / 23.6% | 14.994 | 17,167,040 | yes |
| T7 | 2 | 180.853 | 133.460 | 81.243 / 60.9% | 31.478 / 23.6% | 13.441 | 17,167,040 | yes |

The nested timer starts inside `specialise_shared`; its pivot-reduction
bucket includes lazy materialisation of hit pivots, the word XORs, pivot
search and row insertion. Completion and closure were each below 0.06 ms
on T7. Timer calls, displaced-row sorting and miscellaneous control work
account for 5.0–5.4% of its total as an explicit residual. The profiler
observed 2,470 basis specialisations and 114,197 displaced rows on each
T7 run. All four runs returned the archived scalar (T1 4785, T7 2391),
with unchanged attempt counts and counted F4 word operations. Five
exclusive online phases summed exactly to charged wall. The prepared
F6 enabled/disabled test and the inherited row-space test passed.

Profiling adds per-row clock reads and raised T7 complete online wall
above the uninstrumented 139.96–140.19 ms observation in
[#1469](https://github.com/aburan28/crypto/pull/1469). The percentages
are therefore exploratory attribution under instrumentation; they are
not an A/B performance ratio. The host was unisolated Darwin arm64,
physical CPU model unavailable to the sandbox, and peak RSS unmeasured.
The [four raw runs](runs/), [derived rows](measurements.jsonl),
[derivation check](DERIVATION_CHECK.json), source and raw SHA-256 receipts,
and [runner](run.sh) preserve failures and exact accounting. All four
exited zero; the derivation check reports four verified results, correct
phase sums and exact counted controls.

No n83 ordinary relation, complete n83 F6 call, F4/F5 comparison,
same-target rho interval for this candidate, or 2× complete F6 speedup
was measured. The end-to-end IC objective remains open.
