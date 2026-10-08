# Inherited F6 basis cost is child specialisation

**Decision: inspect and optimize `ReducedBasis::specialise_shared`.** On
the eleven-attempt T7 target, inherited-basis specialisation took
105.4–105.6 ms, or 91.68–91.74% of F4 build, across two repetitions.
Root construction took 9.41–9.51 ms. Within specialisation, the calls
to `specialise_shared` itself took 104.6–104.8 ms; child-system
preparation took 0.41–0.43 ms. This passes the preregistered ≥50%
specialisation gate. The first support-local profile's unattributed
build time was therefore overwhelmingly inherited child work.

The [protocol](PROTOCOL.md) was committed as `4489592b2` before code
or timing. Opt-in instrumentation, exact-result and configuration tests,
and scripts were committed as `fc0780ac8`. The exact
[candidate/input freeze](FREEZE.tsv), [source hashes](SOURCE_SHA256SUMS),
[binary identity](BUILD_IDENTITY.txt) and
[IC1 manifest](candidates/f6_ic_basis_profile.json) were committed as
`e6b181563` before the runs. All [raw runs](runs/) and
[measurement rows](measurements.jsonl) are retained; their
[derivation check](DERIVATION_CHECK.json) passed.

The candidate is
`IC1N17Ckb1fb62PDP3f6RCsampleLAgaussTDpdpISO0hfd2a48433e09`.
It fixes the prepared n17 Koblitz curve, 62 actual usable factor-base
points, 29 folded columns, archived T1/T7 public points, imported
certified logs, three summands, degree three, default algorithm
environment, and one Rayon thread. All four fresh-process target solves
completed and independently replayed the archived scalar (T1 = 4785,
T7 = 2391). Their attempts, reductions, geometric additions, matrix
counts and word operations matched within each target. The five
exclusive online phases summed exactly to online wall.

| Target | Rep | Online ms | F4 build ms | Root ms | Specialise ms | Of which basis ms | Child prep ms | Unassigned ms |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| T1 | 1 | 2.916 | 2.363 | 0.945 | 1.417 | 1.401 | 0.010 | 0.001 |
| T1 | 2 | 2.824 | 2.325 | 0.925 | 1.398 | 1.382 | 0.010 | 0.001 |
| T7 | 1 | 137.395 | 115.168 | 9.514 | 105.581 | 104.829 | 0.414 | 0.073 |
| T7 | 2 | 137.869 | 114.849 | 9.413 | 105.366 | 104.594 | 0.426 | 0.071 |

T1 had one root build, 38 specialisation calls and 35
`specialise_shared` calls. T7 had eleven roots, 2,792 specialisation
calls and 2,470 `specialise_shared` calls. Root time includes the earlier
support-local row, column and packing timers; child preparation and
`specialise_shared` are nested inside total specialisation. The
unassigned F4 build residual is calculated as build minus root minus
total specialisation. The two-run ranges are observations, not
confidence intervals.

The prepared F6 result was unchanged by enabling the timers: recovered
scalar, replay result, F4 basis-read count and word-operation count
matched. The geometric-closure control and default-false
`effective_config` serialization control passed. Profiling overhead is
charged to the online interval; these are diagnostics, not a speedup
comparison. The macOS host was unisolated, physical CPU model could not
be read inside the sandbox, and peak RSS was not measured. There was no
paired one-target rho run or n83 ordinary relation.

Rebuild the rows and checks with
`sh research/f6_ic_inherited_basis_profile_20261006/derive.sh`.
