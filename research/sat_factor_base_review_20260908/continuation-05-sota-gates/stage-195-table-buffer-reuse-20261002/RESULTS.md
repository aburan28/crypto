# Stage 195: same-add BlockTables word-buffer reuse

## Decision

`REJECTED_SCREEN`. The candidate produces a favorable wall sample but does not
reduce total core-seconds, so it fails the preregistered joint
strict-below-`0.98` wall-and-CPU gate. Confirmation is prohibited and the
selected runtime remains unchanged.

| arm | wall seconds | total core-seconds | peak RSS |
|:---|---:|---:|---:|
| current fresh allocation per table run | 17.787954 | 184.000982 | 3,639,541,760 B |
| same-add retired-buffer reuse | 16.902013 | 184.334805 | 4,032,495,616 B |
| **candidate / current** | **0.950194** | **1.001814** | **1.107968** |

Wall falls 4.98 percent in this ordered pair, but CPU rises 0.18 percent and RSS
rises 10.80 percent. The frozen rule does not promote a wall-only sample.

## Mechanism and correctness

The control allocates 728,503 table word buffers. The candidate reuses 26,591
retired buffers and allocates 701,912, reducing allocation count by
3.650088 percent. Reused buffers cover 24,469,504 of 1,700,553,592 scheduled
table words, or 1.438914 percent. The algorithmic table-allocation peak changes
from 53,070,336 to 53,260,128 bytes.

Both arms exactly agree on:

- 5,500 table-add calls, 728,503 table runs, and 1,700,553,592 scheduled table
  words;
- 147,794,583,858 actually performed table/row XORs and 319,313,687,585
  row-equivalent logical XORs;
- all 512 fixed-X1 masks, 270 non-rational skips, and 242 completed rational
  systems;
- source/equation identities, equation and term counts, pair and field-pair
  counters, reducers, matrices, degrees, basis, and zero roots; and
- exhaustive `UNSAT`, with target-subgroup enumeration and known
  discrete-log-label use both false.

Aggregate in-call matrix-build time moves from 45.059463 to 43.468140 seconds,
while aggregate elimination moves from 122.265632 to 134.631584 seconds. These
are sums across concurrent F4 calls; charged process CPU and wall decide the
screen.

The candidate implementation is preserved in
`candidate-reuse-table-buffers.patch` and reverted from runtime source.

## Verification and accounting

The Rust verifier authenticates commands, explicit arm modes, meter receipts,
source/equation identities, algebraic factor-base claims, invariant work,
allocation/reuse counters, ratios, decision, candidate patch, and runtime
reversion. All 16 charged metrics files have exactly one authenticated receipt.
Final replay passes `19/19`; result SHA-256 is
`505b772e029c529d573856f06f906ebd3dbe3d0e5cffba40acf0565435ece628`.

Stage 195 contributes a measured lower bound of 16 components,
`219.075331` wall-seconds, `1,204.681181` total core-seconds, and
`5,302,697,984` bytes peak RSS. The cumulative measured campaign lower bound
is 646 components, `24,599.183898` wall-seconds, `64,219.173358`
core-seconds, and `6,310,576,128` bytes maximum RSS.

Complete campaign cost remains `null`. This rejected one-target allocation
experiment changes no natural relation-yield, unknown-scalar recovery, full
automorphism-aware rho, independent reproduction, novelty, or Koblitz
index-calculus SOTA gate.
