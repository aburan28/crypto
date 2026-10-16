# Disjoint primitive price for the inverse-Vélu fiber panel

Frozen before any run of the instrumented source. This is a follow-on to
`inverse_velu/FOUR_GATE_PROTOCOL_V1.md`, not a revision of the historical
V2 records. The hypothesis is that direct inverse construction costs less
than full-image census after every native primitive is charged. The
same-target source BSGS and signed rho remain the DLP references.

## Boundary and inputs

The generic two-class collision floor is `S=√(π/4)=0.886226925` in
group-addition equivalents divided by `√r`. The comparison reference is
the fully charged same-target signed rho child, with rho A/A as a noise
control. Use the public order-1,201, `p=1237` file
`../portable_order_trend_20261015/order1000/BATCH_1237_1201_V2.public.json`
(SHA-256 `24f2c56a9ba56f91c954f38b7752edf24fee5ed027b49b1ac0747771945099c4`)
and the original inverse runner snapshot
`../portable_order_trend_20261015/inverse_velu/FROZEN_RUNNER_V2.rs`
(SHA-256 `7587fc43ffbfe6042c3915f90d4c7dc164337fce61d0599776f2b607150660c3`).
The separate historical checker and scalar-only input remain the answer
authority. No scalar enters the instrumented solver.

## Instrumentation and price

Copy the frozen runner into a new native Rust example. Add disjoint
integer counters for Fp and Fp² addition, subtraction, multiplication,
squaring and inversion inside versus outside group operations. Count
square-root/Legendre exponent branches and their underlying arithmetic;
polynomial multiplication, reduction, gcd, split and rejected roots;
table comparison, lookup and insertion; modular row arithmetic; heap
allocation count and bytes; startup, input/output bytes, verification,
and every rejected construction or relation attempt. A primitive may
occupy exactly one priced bucket. Group additions and doublings remain
audit counters and are **not** added a second time to field work.

Calibrate each priced bucket against the same native source group-addition
kernel on Darwin and an independent Linux x86-64 host. Pin CPU placement
when supported, retain affinity, host, compiler, binary hash, run order,
frequency/load and repeated A/A samples. Freeze conservative factor
intervals before ranking children. Keep wall and process CPU as separate
observations. Any omitted bucket or conflicting host calibration leaves
portable `S` unset; report the disjoint count vector and host-specific
intervals instead.

## Execution and acceptance

Run the five original arms at prefixes `1,4,16,64`, three rotated
rounds, for 60 cold children. Preserve `8r` relation trials, 32 rho
restarts of at most `20r`, and the 300-second/512-MiB child caps.
Run the original uninstrumented child beside each instrumented child
and compare answer, witness, relation stream, rank and group counters
exactly. The independent Rust checker must replay all 60 instrumented
answers and relation witnesses and recompute every phase count and
total. A mismatch stops pricing and remains in a versioned failed-run
record. Publish cold N1 and whole N64 means and paired uncertainty
intervals against rho and BSGS. Setup reuse is secondary.

The first exploratory extension tests whether raw preimage witnesses
are needed for a projected destination factor base. Its protocol and
inputs must be frozen separately before a run. This price protocol
does not treat a changed base or a projection-only arm as paired
relation-yield evidence for the original five arms.

## Records and stop rule

Keep source/input/binary hashes, exact commands, raw children, failed
attempts, checker receipts, calibration samples and environment in this
directory. A result also needs a source-linked Markdown report, SVG,
PNG and PDF, mirrored into `benchmarks/` with an index row. Stop if
instrumentation changes any V2 answer or witness, a child exceeds a
hard cap, or the independent checker disagrees. A capped or unavailable
host is a censored resource outcome. No frontier or speedup claim follows
from a single host-specific phase ratio.
