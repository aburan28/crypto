# Stage 200: size-gated parallel full-M4RI row clearing

## Decision

`REJECTED_SCREEN`. The corrected candidate reduces total CPU but regresses
wall time, so it misses the preregistered joint strict-below-`0.98` gate.
Confirmation is prohibited and the Stage 199 serial-clear hybrid remains the
repository default.

| arm | wall seconds | total core-seconds | peak RSS |
|:---|---:|---:|---:|
| selected serial row clearing | 20.914397 | 173.316228 | 3,372,613,632 B |
| size-gated parallel row clearing | 22.144250 | 166.943149 | 3,399,729,152 B |
| **candidate / control** | **1.058804** | **0.963229** | **1.008040** |

Parallel clearing saves 3.68 percent total CPU but costs 5.88 percent wall and
0.80 percent RSS. Against the inherited Stage 192 direct-MITM mean, the
candidate still consumes `48.685879x` wall, `544.824280x` CPU, and
`82.785957x` RSS. The decomposition boundary remains decisively negative.

## Exact work and correctness

Both corrected arms authenticate the same public source and equation
fingerprint, retain algebraic factor base `span_F2(1,z,...,z^8)` without
target-subgroup enumeration or known discrete-log labels, visit all 512 masks,
skip 270 non-rational masks, complete all 242 rational systems, find zero
roots, and return exhaustive `UNSAT`.

Both route exactly 481 matrices and 351,164 pivot blocks through full M4RI,
spend 9,222,665,017 table-preparation word XORs, and perform exactly
102,707,985,015 word XORs. After removing only phase timing, process timing,
and the three new scheduling counters, their complete solver JSON records are
byte-equivalent under canonical JSON ordering.

The candidate sends 75,803 pivot blocks, 611,986,047 target-row visits, and
106,712,818,981 scheduled target-row words through the existing Rayon pool.
The control reports zero for all three counters.

## Preserved failure and correction

The first candidate attempt terminated with return code 101 after nested Rayon
work re-entered a worker's thread-local M4RI scratch while its outer call still
held a `RefCell` borrow. It consumed 22.299680 wall-seconds, 161.557050
core-seconds, and 3,357,638,656 bytes peak RSS. It produced no solver JSON and
is not timing evidence.

The correction takes the scratch value out of the thread-local slot before
entering nested work. A regression runs twelve outer M4RI calls concurrently,
each forcing parallel row clearing, and requires exact rows, pivots, row space,
logical XORs, performed XORs, and table-preparation XORs. A fresh corrected
binary then reran the complete frozen `serial, parallel` pair. The failed
attempt and both source builds are preserved and charged.

The corrected candidate is retained in
`candidate-parallel-m4ri.patch` at SHA-256
`bf4269950abd6c5b2da741e0394739d8b68bde2ff32edcda7bce0f1e9d871a79`;
the runtime source is reverted.

## Verification and accounting

The Rust verifier authenticates both exact commands, the absence of full-M4RI
mode and row-threshold overrides, selected-default routing, every receipt,
exact terminal and work counters, normalized report equality, ratios, the
frozen decision, Stage 199 parent hash, failed-attempt custody, candidate
patch, and runtime reversion. Final replay passes `27/27`; result SHA-256 is
`5a1530a9f89645076b99c1c027996793c107f7eacb1d4907ec1cdedb76756f31`.

Stage 200 contributes a measured lower bound of 26 components,
`1,088.370919` wall-seconds, `4,553.408345` total core-seconds, and
`5,592,170,496` bytes peak RSS. The cumulative measured campaign lower bound
is 722 components, `26,616.264678` wall-seconds, `73,077.015408`
core-seconds, and `6,310,576,128` bytes maximum RSS.

Complete campaign cost remains `null`. This is a rejected one-target F4
scheduling experiment, not relation-yield, unknown-scalar, full-rho,
independent-review, novelty, or Koblitz index-calculus SOTA evidence.
