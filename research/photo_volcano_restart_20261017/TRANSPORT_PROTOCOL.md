# Full source-to-vertex transport instrumentation

Use the existing 100 public P/Q pairs on source C0 and each vertex's frozen
first-discovery path, including the exact codomain isomorphisms. This is a
separate reference-instrumentation phase: do not change frozen maps, fixtures,
orders or native solver cases. Freeze source/binary/graph/fixtures first.

For every one of 46 destinations, run all 100 point pairs through the entire
path five times. Time actual reference map evaluations with wall/process-CPU
nanosecond clocks and count every field add/subtract, multiply and inversion,
including inverse exponentiation. Parse/precompute paths and kernels outside
the measured interval. Verify both resulting P/Q images against the original
destination fixtures outside the interval. All point orders were independently
checked in the original phase and the order phase; do not repeat order work
inside the transport interval.

Retain all 23000 case records and 46 exact paths. The empty C0 path has zero
field operations and is valid, not missing. A single reference edge transports
the pair using 174 counted field multiplications; the full-path counter is
checked against that count times the number of steps. All original exact
rational map/dual proofs and subgroup injectivity checks remain applicable.

These are instrumented reference maps, not a fastest-map comparison. Timings
are CPU-only L0 with concurrent BSGS. Five repetitions of a paired target are
not independent scalar samples. Do not add a map timing mean to an unrelated
solver timing and call it a separately measured end-to-end solve. Derived
sums, if reported, must retain that label and boundary. This phase does not
collect missing rho/BSGS internal field counters or complete IC recovery.
