# Stage 201: exact combination ends in full M4RI

## Decision

`REJECTED_SCREEN`. Exact combination ends reduce performed work and total CPU,
but the wall-time ratio misses the preregistered strict-below-`0.98` gate.
Confirmation is prohibited and the Stage 199 block-wide-end implementation
remains the repository default.

| arm | wall seconds | total core-seconds | peak RSS | performed XORs |
|:---|---:|---:|---:|---:|
| selected block-wide end | 25.383130 | 163.804192 | 3,250,962,432 B | 102,707,985,015 |
| exact combination ends | 26.508356 | 161.673973 | 3,296,854,016 B | 102,515,969,703 |
| **candidate / control** | **1.044330** | **0.986995** | **1.014116** | **0.998130** |

The candidate saves 192,015,312 performed word operations, or 0.186953
percent, and 1.30 percent total CPU. Its bookkeeping and trailing-zero scans
cost 4.43 percent wall and 1.41 percent RSS. Against the inherited Stage 192
direct-MITM mean, it still consumes `58.280711x` wall, `527.628156x` CPU, and
`80.280870x` RSS. The decomposition boundary remains decisively negative.

## Mechanism and exact work

The candidate keeps an exact last non-zero word for every full-M4RI
combination while preserving the control's block-wide row-end metadata. It
therefore changes neither later pivot/table schedules nor logical work.

It trims 5,057,385 table entries and 106,080,228 row lookups. The exact savings
are:

- 7,706,619 table-preparation words; and
- 184,308,693 row-application words.

Those sum to 192,015,312 and replay exactly to the control/candidate performed-
XOR difference. Candidate table-preparation work is 9,214,958,398 versus
9,222,665,017 in the control. Both arms route 481 matrices and 351,164 blocks
through full M4RI and report the same 318,703,372,596 logical XORs.

Both arms authenticate the same public source and equation fingerprint,
retain algebraic factor base `span_F2(1,z,...,z^8)` without target-subgroup
enumeration or known discrete-log labels, visit all 512 masks, skip 270
non-rational masks, complete all 242 rational systems, find zero roots, and
return exhaustive `UNSAT`. After removing only timing, exact-end mechanism,
performed-work, and table-memory fields, the solver records are identical.

The differential test covers varied row ends, partial blocks,
non-consecutive pivots, and word boundaries. It requires exact pivot rows,
rank, canonical row space, logical work, and both savings identities.

The rejected candidate is retained in `candidate-exact-end-m4ri.patch` at
SHA-256
`bc1cd7d4b54ba751a4a3b5201037be6caeca923ba8c1f144a61ca26b69dbdd3f`;
runtime source is reverted.

## Verification and accounting

The Rust verifier authenticates both commands, selected-default routing, every
receipt, exact terminal and work counters, savings identities, normalized
report equality, ratios, the frozen decision, Stage 200 parent hash, candidate
patch, and runtime reversion. Final replay passes `22/22`; result SHA-256 is
`1f267740b9823c73f8bceb687a46f4e91199ce512fd74f208a18a8cbdd0a2813`.

Stage 201 contributes a measured lower bound of 16 components,
`1,016.740619` wall-seconds, `2,349.721545` total core-seconds, and
`4,559,831,040` bytes peak RSS. The cumulative measured campaign lower bound
is 738 components, `27,633.005297` wall-seconds, `75,426.736953`
core-seconds, and `6,310,576,128` bytes maximum RSS.

Complete campaign cost remains `null`. This is a rejected one-target F4
implementation experiment, not relation-yield, unknown-scalar, full-rho,
independent-review, novelty, or Koblitz index-calculus SOTA evidence.
