# Group and point order measurements, 2026-10-17

This is a separately frozen instrumentation phase on the existing exact
46-curve/100-instance corpus. It changes no original fixture, native solver,
rho/BSGS checkpoint or earlier measurement. No new random input is sampled.

PARI/GP 2.17.2 creates a fresh curve object for every vertex and times, in
order, `ellcard(E)`, factorization of its returned cardinality, and
`ellgroup(E)`. Later stages can reuse PARI's cached curve information and
are explicitly not cold independent curve-count runs. It then calls
`ellorder(E,P,[N,factor(N)])` on both P and Q for all 100 fixtures per curve:
9200 actual API calls, including the workload's repeated generator points.
Record both CPU and wall milliseconds. GP's timers have 1 ms resolution;
a returned zero is below the timer resolution, not zero actual cost.
Whole-process resources include all parsing, output, validation and shutdown.

Independently run a native reference point-order algorithm on the same
9200 P/Q inputs. Start with the supplied exact curve cardinality N=4r and
try each prime factor (2 then r), dividing the candidate order only when
the corresponding multiple is the identity. Native wall/process-CPU clocks
and field/group counters cover this reduction algorithm, including its
scalar multiplications. Input curve-membership checks, checking the returned
order, cardinality computation and factorization are outside that interval.
This is a reference implementation, not a fastest-tool comparison.

Count actual field add/subtract calls, multiplications, squares, inversions
(Fermat exponentiation included), scalar calls, group-operation calls and
their nontrivial/identity/inverse-pair branches. Multiplication-including-square
already includes the separately reported square count; never add them again.
Constant field multiplications are counted. Record all counts, not estimates.
These counters describe point-order computation, not the uninstrumented
rho/BSGS solver internals. Do not fill those solver counters from this phase.

Positive correctness controls on the source curve have known orders 1, 2,
r and 2r. The independent uninstrumented native EC implementation validates
every returned prime point order. PARI's actual point-order results provide
an external replay on all 9200 points. Exact curve orders and group invariant
factors must match the original graph; every point order must equal r.

Sources, binary, GP script, fixtures and graph identities are hashed in
ORDER_FROZEN.json before either measured main execution. The build failure
and syntax-probe preflight are preserved separately and not pooled as main
measurements. Main runs are write-once; incomplete or failed evidence remains
retained. All measurements are CPU-only L0 with concurrent BSGS work. No speed,
hardness, isolated timing, or whole-campaign completion claim follows.
