# Degree-indexed Boolean support-envelope experiment

This standalone experiment asks whether exact product-schedule reuse can beat
fresh Boolean Macaulay matrix construction when coefficients change. Its inputs
are generated public quadratic systems with at most 12 Boolean variables. It
does not run a curve, key, relation-collection, or production solver path.
These measurements are construction-stage diagnostics; full-IC cost, calibrated
operation cost, and the rho ratio remain null.

## Frozen mathematical contract

The predecessor `research/boolean_support_envelope_20260922` proved exact reuse
across coefficient changes but passed 0/48 performance gates. Its worker,
protocol, and completed-run manifest have the SHA-256 digests recorded in
`protocol.json`. The new binary retains its five arms: direct construction,
verified column-layout reuse, exact product scheduling, exact matrix caching,
and the original eager support envelope. It adds two treatments in that same
binary; all seven arms return the same ordered, compact matrix or the same
explicit cap error as the independent set-parity oracle.

For each ordered generator envelope, compilation retains the predecessor's
groups mapping `(coefficient slot, multiplier)` to output monomial
`source OR multiplier`. The parity of all active slots in a group is evaluated
at each application; cancelled terms and empty rows disappear. The current
generator's actual degree determines eligible multipliers. A zero generator
contributes no rows, a constant generator admits the complete multiplier
range, and support escapes invoke the unchanged direct fallback. Current row
and column caps are applied after cancellation and compaction.

The first new arm precomputes, for every permitted multiplier-degree cutoff,
the indices of eligible plans in their original numeric multiplier order.
It visits only those plans at application. This keeps every row in exactly the
old order, including rows newly needed when a generator's degree drops.
The second arm adds one sorted universe of possible output columns. Each
application marks columns from odd-parity terms only, scans the occupancy
bitmap to create the actual compact columns, and maps row terms into their
current dense indices. No absent column is retained in the returned matrix.
The universe is a routing aid, not a cached completed matrix.

Compilation refuses if any plan-row, group, or retained-capacity bound in the
frozen protocol is exceeded. Both new arms use the same 512-row and 512-column
current-output caps as all controls. Their retained bytes include the extra
indices and column universe. Per-arm output payload and whole-worker peak RSS
are reported separately; whole-worker RSS cannot be assigned to one arm.

## Inputs, costs, and acceptance

The fixture generator and five controls are ported from the hash-pinned
predecessor. Sizes are n=6, 8, 10, 12, matrix degree three; families are
repeat, changing coefficients, degree cycle, and support escape; batches are
1, 4, 16, 64. Each phase has two fixed seeds, 128 cells, seven arms and ten
rotated repetitions. The old discovery and holdout seeds are development
references only. `protocol.json` fixes new discovery and untouched holdout
seeds before performance work. No holdout timing is examined before the timed
source and native verifier are frozen.

Each A/B observation charges a newly compiled arm, every application and
fallback, fresh output construction, exact result validation, and destruction.
Common fixture generation and oracle preparation are outside arm clocks and
inside the worker resource receipt. For every cell and repetition, two calls
to the same matrix-cache control form a separate A/A calibration. The worker
returns output digests and counted work; a separate Rust verifier recreates
the expected ordered matrices using the independent set-parity oracle.

The primary dramatic criterion is unchanged-output completion in every cell
and a paired 95% bootstrap lower bound above both **2.0** and that cell's
97.5-percentile symmetric A/A ratio, against the pointwise fastest of the five
old controls. It must hold for either fixed new arm in **all 32** combinations
of four sizes, batches 16/64, changing-coefficient/degree-cycle families,
and two fresh holdout seeds. Improvement above 1.0 is a separate engineering
diagnostic. Discovery cannot establish the primary claim. A positive full
result requires unchanged-source confirmation on further unused seeds.

For correctness, exhaust the predecessor's 8,192 coefficient, degree and
active-mask cases, plus focused cancellation, constant, zero, degree-drop,
support-escape and cap cases. Every expected in-envelope hit and fallback is
checked. The predecessor's 256-cell data and timings remain immutable. A
failed, contended, timed-out, capped or incomplete new run is retained but
contributes no accepted timing. The n6–n12 test does not imply scaling or a
cryptanalytic result.

## Resource and implementation rule

The new producer, verifier, manifest writer and result renderer are Rust.
Shell only starts builds and runs. The repository's required
`tools/isolated_bench.py` remains the external isolation controller. On Linux
ARM64 a single campaign worker executes the frozen cells in order under one
reserved core and one whole-campaign resource receipt. Each solve still has
its own timer, output digest and A/A pair. The receipt's 10% other-CPU and 5.0
PSI limits are unchanged; its scope is the entire campaign, so a short cell
does not receive an individual resource-certification claim. The campaign
retains every fixture and never pads, divides or retries timed samples.
The host, CPU feature set, compiler, binary/source hashes, A/A noise, wall
time, peak RSS, failures and exact commands accompany the result. Separate
fresh discovery and full campaigns use the same timed Rust source hashes.
