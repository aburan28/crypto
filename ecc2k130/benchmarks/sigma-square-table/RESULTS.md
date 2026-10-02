# Sigma-fused lambda square table: native admission result

Status: **native pass; exact sm_120 CUDA compile and resource audit pending.**
No GPU was launched and there is no throughput result.

## Novelty decision

This exact arm has not been measured before.  The existing
`PACKED_SQUARE_TABLE` path is confined to `ECC_WALK_TABLE`.  Its only retained
GPU promotion evidence is combined with `PACKED_INV_POLY=2` on the different
B16/T512 block-v3 table walk, where the paired median was 1.006625.  The
sigma-fused tuning star forced `PACKED_SQUARE_TABLE=0`; its rejected
`PACKED_ALU_SQUARE=1` arm is a separate logic-spread implementation.

## Implemented one-knob path

`SIGMA_SQUARE_TABLE=1` is default-off and requires the polynomial-state
sigma-fused walk.  It is deliberately incompatible with
`PACKED_SQUARE_TABLE=1` and `PACKED_ALU_SQUARE=1` for this first comparison.

The host fills the existing exact 2,080-word table through
`fillSquareTable131` and stores it in the ABI-compatible `twConsts` field.  The
kernel copies it to 8,320 bytes of dynamic shared memory and synchronizes the
whole block before the partial-block return.  The reverse pass calls the same
compile-time `sigmaLambdaSquare131` wrapper exercised by the native replay.
The control keeps the original `squarePolynomial131` call.

Source shared allocation is 1,792 static bytes for shared-sigma masks plus
8,320 dynamic table bytes, or 10,112 bytes/block and 20,224 bytes for two
blocks.  Target ptxas register, spill, static-shared, and binary/SASS evidence
remain pending; source arithmetic alone does not establish occupancy.

## Native evidence

Both compile-time arms replayed the same 20,134 inputs:

- 131 basis vectors;
- zero, one, and the all-ones canonical edge;
- 20,000 deterministic dense polynomial-basis elements.

Every output matched both `squarePolynomial131` and the independent field
reference.  The two 402,680-byte output streams are byte-identical.  A separate
strict C++ audit reopened both streams, required canonical outputs, recomputed
all 2,080 table words from basis squares, verified that every fixed
`(word,window)` address has bank index equal to the 5-bit entry, and reopened
the source bindings and guards.

The repository's broader four-build `test-packed-network` matrix also passed.
Each build independently reported 2,526 long-division reduction controls,
20,134 table-square controls, 2,133 polynomial-inversion controls, and the
packed multiplication/Frobenius suite.

## Exact source ledger

There is one affected lambda square per complete scalar update:

| per update | control | candidate | delta |
|---|---:|---:|---:|
| carry-less 32-bit spreads | 5 | 2 | -3 |
| dense beta reductions | 1 | 0 | -1 |
| 5-bit window extractions | 0 | 13 | +13 |
| 32-bit shared reads | 0 | 65 | +65 |
| XOR accumulations | 0 | 65 | +65 |

Per candidate block and kernel launch, the cooperative prologue additionally
performs 2,080 global reads, 2,080 shared stores, and one barrier.  These are
source and memory-operation counts, not SASS counts or a predicted speedup.

## Next gate

The focused PR's no-GPU CI compiles both exact B16/T256/minBlocks2 arms for
native `sm_120` with CUDA 13.3.73 and retains their binaries, ptxas logs,
resource reports, SASS, and hashes.  Those artifacts require a separate native
audit before the frozen matched GPU run is admissible.  Until then the flag
remains default-off and the decision remains
`NATIVE_PASS_CUDA_COMPILE_PENDING`.
