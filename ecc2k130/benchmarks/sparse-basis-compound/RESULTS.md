# Sparse-basis compound route: negative static result

Decision: **do not build the whole-walk GPU candidate under this protocol.**
Both sparse fields and every proposed linear map are exact, and the sparse
reducer is much smaller than the current direct reducer at the source-circuit
level.  The required dense linear maps move the wall to shared-memory lookups.
Transferring the RTX PRO 6000's measured random-`LDS.U8` rate to those lookups
gives a conditional lookup-only projection of **15.283375 B complete scalar
updates/s**, below
the confirmed fused reference of **15.436677 B/s**.  This gives the arithmetic,
state traffic, loop/control work, and all non-table memory accesses zero cost.
The reference is pinned to commit
`e0b0858b179fe28f753353e7a7313cfe23d71884`, headline-result SHA-256
`d31f2d758d2e787abe8f043b1b2ae69d337ebad9a30403afe74d50d4a003dd8c`,
and independent-audit SHA-256
`1b1a46ca52aa14ef86b96fc9b765c9c10e2815c1fe9bcad867422785ab7976d1`.

This is negative static engineering evidence under the frozen transfer model,
not an unconditional shared-memory hardware ceiling.  No GPU was launched, and no
search, distinguished-point collection, collision solve, or key recovery was
performed.

## Exact native checks

The C++17 audit independently constructs each field, an order-263 element, the
isomorphism to the repository beta polynomial basis, the normal-coordinate
map, and `L_j(a) = a + a^(2^j)` for `j=3..10`.  The first attempted protocol
incorrectly placed the order-263 root in a quadratic extension.  The native
run rejected that construction before emitting a result; the correction that
`2^131 mod 263 = 1` was committed before this admitted run.

| check, per candidate | count | result |
|:--|--:|:--|
| raw-product reduction basis directions | 261 | pass |
| deterministic dense raw products | 4,096 | pass |
| conversion basis round trips | 131 | pass |
| dense conversion checks | 8,192 | pass |
| normal / eight `L_j` basis-map checks | 1,179 | pass |
| normal / eight `L_j` dense-map checks | 36,864 | pass |
| field product basis pairs | 17,161 | pass |
| deterministic dense field products | 4,096 | pass |

Both pentanomials pass the complete prime-degree Rabin conditions.  Both basis
conversion matrices and both sparse-to-normal matrices have rank 131.  Every
`L_j` matrix has rank 130, the required kernel dimension for `I + sigma^j` in
prime degree 131.

The matrices are dense even though multiplication reduction is sparse:

| modulus | beta to sparse ones | sparse to beta ones | sparse to normal ones | `L_3` .. `L_10` ones |
|:--|--:|--:|--:|:--|
| `z^131+z^8+z^3+z^2+1` | 8,551 | 8,559 | 8,541 | 1,474 / 2,567 / 3,894 / 5,302 / 6,511 / 7,261 / 7,790 / 8,266 |
| `z^131+z^8+z^5+z^2+1` | 8,515 | 8,619 | 8,493 | 1,486 / 2,602 / 4,107 / 5,573 / 6,765 / 7,512 / 8,023 / 8,350 |

The first pentanomial is weakly preferable for any later non-table synthesis
because seven of its eight `L_j` matrices and the combined two-way basis maps
have fewer one-bits.  The admitted three-bit evaluator has the same runtime
shape for both, so this does not move the present decision.

## Reduction circuit

For either sparse modulus, write

```
R(z) = 1 + z^2 + z^b + z^8,  b in {3,5}
D    = input >> 131
S    = D * R
E    = S >> 131
out  = (input mod z^131) + (S mod z^131) + E * R
```

`deg(D) <= 129`, `deg(S) <= 137`, and `deg(E*R) <= 14`, so exactly two
folds suffice.  Equality on all 261 raw-product basis vectors proves this
linear reducer for every degree-at-most-260 input; the 4,096 dense cases are
an implementation guard.

The exact 32-bit source circuit count is:

| reducer | shifts | OR | XOR | AND | total | sparse / current |
|:--|--:|--:|--:|--:|--:|--:|
| current `reducePolynomial131` direct circuit | 86 | 39 | 53 | 2 | 180 | 1.000 |
| either sparse two-fold circuit | 40 | 16 | 24 | 2 | **82** | **0.4556** |

These are expression-level word operations.  They are exact for the checked
C++ circuits but are not Blackwell SASS counts; ptxas can fuse shifts and ORs
into `SHF` and XOR trees into `LOP3`.  The decision below does not assume any
particular compiled reduction cost: it makes all reduction free.

## Complete batch-16 ledger

The fused sigma schedule performs 85 polynomial products per 16 updates: 30
in the forward prefix, eight in the beta-basis inverse, 31 in the reverse
prefix, and 16 for the final coordinate product.  The 31 reverse products
already include all 16 `inv*W_i` products that yield lambda together with the
15 `inv*D_i` updates; charging a separate lambda stage would count them twice.
The compound route therefore has:

- 5.8125 sparse reductions/update, including the polynomial lambda square;
- zero direct beta reductions inside the unchanged default inverse: its eight
  products use `fromPolynomialProduct131` and remain a common cost;
- 38.125 `CLMAD` lane-instructions/update, unchanged from the fused arithmetic:
  510 for products, 20 for five normal-basis inverse squarings, and 80 for the
  16 polynomial lambda squarings per batch;
- 3.125 dense map applications/update: sparse X to normal, direct `L_j` on X
  and Y, plus sparse-to-beta and beta-to-sparse once each per batch; and
- 990.625 table-evaluator data-ALU operations/update.

At the source-circuit level, the reductions fall from 1,046.25 to 476.625
word operations/update.  The selection/direct-map stage, including the two
amortized batch-boundary maps, falls from the prior composed model's 1,230 to
990.625 under the frozen conservative charge (984.375 for the exact
straight-line evaluator).  The conservatively charged combined static saving
is 809 word operations/update.  That saving is real but does not pay for the
lookup pipe it introduces.

## Shared memory and live storage

Each map has 44 three-bit chunks, eight entries/chunk, and five output words:
7,040 bytes.  Nine hot maps and two batch-boundary maps occupy **77,440
bytes**.  With the driver's 1,024-byte reservation the block needs **78,464
bytes**, below the RTX PRO 6000's 101,376-byte opt-in limit.  It permits one
block per SM.

The aligned layout uses one `LDS.128` for the low four words and one `LDS.32`
for the top word per chunk.  A map therefore costs 88 LDS instructions and
220 scalar-word equivalents.  Its exact straight-line evaluator has 315
data-ALU operations; the frozen admission ledger conservatively retains the
prior direct-sigma charge of 317.  At 3.125 maps/update that is **275 LDS
instructions**, **687.5 scalar-word equivalents**, and 990.625 charged
data-ALU operations per update.  The evaluator needs five input words, five
output words, and one index live.  Adding those 11 words to the fused kernel's
measured 126 registers gives a naive 137-register overlap estimate; replacing
20 named current-reducer temporaries with 10 gives a separate net named-word
estimate of 127 if that reuse transfers.

The distinction matters for geometry.  An SM has 65,536 registers, so a
512-thread block can use at most 128 registers/thread.  The 127 estimate fits;
the 137 estimate does not.  Because 78,464 shared bytes permit only one block,
a 256-thread build exposes eight resident warps and a 384-thread build twelve,
versus the fused reference's sixteen.  Only an actual CUDA compile can resolve
the overlapping lifetimes, spills, and whether a 512-thread block remains
launchable.  The lookup-only rejection assumes the favourable geometry and
therefore does not depend on this uncertainty.

## RTX PRO 6000 roofline

The table uses the most favourable recorded clock, 2.430 GHz, and 188 SMs.
The rate source is `benchmarks/fast-clmad/probe/summary.json` (SHA-256
`57d8541c28648911e3e59a9a6bdd34f114658f5df0b299bab651193293b4dc50`);
register/shared limits use `benchmarks/hardware-limits/result.json` (SHA-256
`31d1c91ca70006086d0454f51be5a5f5866afdbab6c06adcda88026c88a54c1a`).
Every row makes all other work free.  The first row conditionally transfers a
measured rate from random `LDS.U8`; the second shows the exact table layout's
unmeasured 16-lane sensitivity.

| limiting work | rate used | work/update | modelled B/s | / 15.436677 reference | / 16.208511 admission |
|:--|--:|--:|--:|--:|--:|
| three-bit tables, measured random-shared rate | 9.2 lanes/SM-clock | 275 LDS | **15.283375** | **0.9901** | **0.9429** |
| three-bit tables, sensitivity at ideal 16 lanes/SM-clock | 16 | 275 LDS | 26.579782 | 1.7218 | 1.6399 |
| table-evaluator ALU only | 64 lanes/SM-clock | 990.625 ALU | 29.514458 | 1.9120 | 1.8209 |
| unchanged carryless arithmetic only | 2 lanes/SM-clock | 38.125 `CLMAD` | 23.965377 | 1.5525 | 1.4786 |

The frozen gate requires every modelled one-pipe rate to clear the fused
reference by 5%.  The transferred random-shared row fails before charging the
rest of the kernel.
The unchanged carryless count also caps this route below the separate 26 B/s
objective even under its ideal two-lane rate.

The 16-lane table row is sensitivity, not evidence that this access pattern
achieves that rate.  The table keys are data-dependent, and the repository's
measured random-shared stream is 9.2 lanes/SM-clock.  A small, separately
frozen `LDS.128`/broadcast microprobe could test whether this exact eight-entry
layout behaves more like the ideal row.  A whole-walk sparse-basis prototype
is not justified unless that prerequisite establishes at least the frozen 5%
margin after charging its own extraction work.

## Reproduce

From `ecc2k130/`:

```sh
benchmarks/sparse-basis-compound/run.sh
```

The command compiles and runs only the native C++ audit, rewrites the canonical
[`result.json`](result.json), and records hashes in
[`files.sha256`](files.sha256).  It does not invoke CUDA or a remote provider.
