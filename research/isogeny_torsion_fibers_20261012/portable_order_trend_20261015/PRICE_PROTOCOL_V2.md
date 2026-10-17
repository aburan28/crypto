# Complete portable-cost price gate, V2

Status: frozen before a new price-instrumented run. V1's
`S_group=(adds+doubles)/(N sqrt(r))` is a checked **lower bound** on
the complete portable price. The four cells, seven arms, public targets,
caps and memory limit remain those in `BATCH_RUN_SPEC_V1.json`; do not
replace a failed child or tune a conversion after seeing its ranking.

## Accounting contract

Each child must emit a disjoint phase ledger for common validation,
field/curve construction, map and factor-base setup, relation collection,
row recovery, target solve, source/image verification, and output. Record
the following native primitive counts **inside and outside** the group
calls: base-field additions/subtractions, multiplications and modular
inversions; exponentiation and square-root branch steps; modular row
multiply/add/inversion; hash/table probes and inserts with key width;
point comparisons/canonicalizations; bytes allocated, peak resident
memory, bytes read/written and process setup. Do not add a base-field
count to `S_group` if it was already paid inside a counted group call.
Record failed candidates, rho restarts and all cap usage.

The canonical `S` price is the cold sum of group additions/doublings plus
each *disjoint* native count times its pinned addition-equivalent factor,
divided by `N sqrt(r)`. Derive a factor and an interval from preregistered
native calibration kernels on the host class, then freeze them before
ranking the 252 children. Repeat the kernels on the second host. If a
factor's intervals conflict by host, retain host-specific `S` and the
unpriced vector; do not silently call a host-dependent conversion
portable. Keep CPU, wall and memory in their physical units alongside
`S`, never fold them into `S` twice. A conservative interval for a
candidate/control quotient must include calibration and paired-run
uncertainty.

Instrument the current Rust path directly. The historical Python records
are preserved but do not execute in this protocol. Add a native checker
that recomputes phase sums, verifies all scalars and relation witnesses,
and checks that the same public inputs and acceptance caps were used.
Cross-check instrumented and uninstrumented V1 group and Fp counts on a
small fixed public prefix to catch accounting drift; any change in
answers or relation witnesses fails the price run. Confirm counter
coverage with source review and at least one independent operation trace.

## Decision

The price gate passes only if every required native category is counted
or has a defensible bounded interval in one group-addition-equivalent
unit, every child and checker passes, and the conversion is pinned and
rederived on a different host class. A remaining unpriced category keeps
`complete_portable_S` unset. A componentwise candidate loss can be
reported for the recorded axes without a point conversion; it does not
waive the complete-price gate. To propose a frontier change, require the
complete cold `S` interval below matched rho on growing orders, plus
separately admitted isolated wall and A/A evidence.
