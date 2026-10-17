# Native solver arithmetic replay, full frozen corpus

This is a separate instrumentation replay of all 4600 rho cases and all
460 BSGS generator batches (4600 queries), using the same frozen P/Q inputs,
walk seeds, algorithm parameters, table size and single-threaded native calls.
It changes no baseline source, binary, timing, checkpoint, or reported case.
The instrumented source copy and every inserted hook are retained separately.

Counters record executed modular add/subtract, plain modular multiplication,
Montgomery multiply/square/reduction, modular inverse, extended-Euclid division
iterations, and explicit EC coordinate-complement calls. Partition counters
by F_p, scalar modulus r and other modulus. Montgomery reductions are included
inside multiplies/squares as well as conversions: do not add these categories
as if they were independent field operations. Inverse counters do not stand
for unit-cost multiplies; this backend uses extended Euclid, unlike the
separately measured Fermat reference maps/order code. These are source-level
primitive invocations, not retired CPU instructions or every integer/RNG/hash
operation. Full per-step state streams and PMU counters remain uncollected.

Initialize/encode the group before the measured API interval. Reset and enable
thread-local counters immediately before each native solve. For rho, threads=1
executes the native worker inline; no counter aggregation over spawned threads
is claimed. For BSGS, measure table setup, every query and table destruction
separately. Read baseline answer records only after native return, and verify
the public point equation independently. Require identical returned scalar,
group-operation count, iteration count, collision count and table entries to
the corresponding baseline case. Never feed baseline answers into either solver.

Freeze source, hooks, binaries, parameters, fixtures, public seeds and initial
frozen identity before the first replay. Atomic per-case files and complete
job markers support resume. Preserve failures and interrupted attempts; a
rebuilt table has a new setup receipt. A complete marker requires every case
in that job to match the baseline and pass the independent equation check.

Run two rho replay workers while the original BSGS campaign finishes. Start
the two BSGS replay workers only after all 460 baseline BSGS job markers exist
and validate. At most two 4-GiB BSGS tables are resident in either phase, plus
small rho state, under the existing 16-GiB cap. This adds no remote instance.
All instrumented timings are CPU-only L0, with counter overhead and concurrent
work; they are not pooled with original timings or used for speed claims.

Passing modular/group/rho/BSGS tests and one actual paired replay are execution
checks, not full completion. Full coverage is exactly 4600 matched cases per
algorithm. Full IC recovery and unsupported instruction/state metrics remain
separate obligations, even if this replay completes.
