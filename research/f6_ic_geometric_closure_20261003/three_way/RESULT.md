# Matched F4/F5/F6-IC prepared one-target comparison

The [precommitted protocol](../THREE_WAY_PROTOCOL.md) froze one new public
target `[79391,5777]` on `icv1-f2m17-tm101-00378d4e` (EC1 alias
`EC1N17Ckb1hbbe2b5b6b1e6`) from seed `20261004001`, without constructing
or supplying its scalar. Its workload ID is `efc282fb1df2`. All three arms used
the same certified log table, 62 actual usable factor-base points, 29 folded
columns, two target-dependent attempts (one proved unsatisfiable, then one
witness), one CPU worker, and the same 32-trial/600-second limits. The nine
processes all exited zero, recovered `32917`, and passed independent
general-group scalar replay. Every five-phase online ledger summed exactly
to its reported interval. The [freeze receipt](freeze.txt) binds the common
binary, source, candidate manifests, fixture and three inputs; the
[keyed run rows](measurements.jsonl) retain every attempt, phase and status.

| Repetition | F4 online ms | F5 online ms | F6-IC online ms | F4 / F6 | F5 / F6 |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 31.777 | 6216.168 | 21.483 | 1.479× | 289.347× |
| 2 | 32.134 | 6295.054 | 20.906 | 1.537× | 301.113× |
| 3 | 32.196 | 6338.419 | 20.779 | 1.549× | 305.037× |
| Median | **32.134** | **6295.054** | **20.906** | **1.537×** | **301.113×** |

These are complete **prepared one-target online calls**, from the first
target-dependent query through independently checked recovery. Target PDP
consumed 99.83–99.86% of F6, 99.92–99.93% of F4, and over 99.999% of F5.
The measured process repetitions describe timing variation on **one fixed
target**, not natural target-to-target variation or a confidence interval.
The Mac had no auditable host-level CPU partition, so all wall ratios are
**exploratory**, not controlled speedup claims.

For the same two attempts, F4 made 228 Boolean reductions and 95 splits;
F5 made 8,184 reductions and 3,138 splits; F6-IC made 98 reductions and 53
splits. F6 also charged 10,713 logical point additions, 5,349 exact residual
lookups and 170 point-addition batch groups to target PDP. These counts are
diagnostics, not a common cross-solver operation unit. F5 used its production
matrix solver and lowest-free split; F4 and F6 used the inherited-basis and
highest-free route. Thus the ratios compare complete solver-family choices,
not one isolated matrix kernel. Candidate IDs are in the frozen manifests;
these new manifests also correct the relation, final-LA and descent source
digests that the earlier F6 pilot inherited from an old F4 record. The
executed source and binary match the earlier v3 pilot.

**Decision:** F6-IC beats both comparators descriptively on this fixed
target, but the preregistered 2× complete-call F4 gate fails in all three
pairs. The earlier v3 F4/F6 target likewise missed 2×. F6-IC remains an
IC-specific hybrid experiment rather than a promoted replacement. On this
small prepared workload, F5 is a poor comparator and inherited F4 is the
meaningful incumbent. The factor-base/log preparation was already available;
the worker's local setup timings are retained separately in the raw output,
while the imported log table's original construction receipt is in the
prepared-state record. Natural relation yield was not remeasured. Peak
memory is unavailable. No same-target rho arm ran through `ecbench`, so
IC/rho online ratio, operation-normalized `S`,
strong-rho crossover, and any change to the current n53 `PDP4root` path are
**unknown**. No larger-field or asymptotic conclusion follows.
