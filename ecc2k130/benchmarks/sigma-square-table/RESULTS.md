# Sigma-fused lambda square table: native admission result

Status: **native and exact sm_120 compile pass; independent review and device
execution pending.** No GPU was launched and there is no throughput result.

The protocol was frozen at `18ee85c1075010ce83be60a4a02636e9aa4ac83a`.
Before device dispatch, the implementation was merged through current main
`bc217d318dde444014cde4e81dea02c70b126995`, including the selected GPU-wide
hint evidence and the separate fused shared-scratch static screen.  Those
merges do not enable either feature in the frozen sigma baseline or candidate;
all native evidence below was regenerated afterward.

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

The source-visible masks occupy 1,792 static bytes and the table requests 8,320
dynamic bytes.  Exact CUDA 13.3.73 compilation adds the compiler/runtime's
existing 1,024-byte shared base: `cuobjdump` reports 2,816 compiled shared bytes
for both walk binaries.  The candidate therefore requests 11,136 bytes/block
and 22,272 bytes for two blocks.  That is a resource-feasibility result; the
device occupancy API remains an execution gate.

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

## Exact CUDA 13.3.73 / sm_120 compile

The no-GPU producer built both exact B16/T256/minBlocks2 clients and retained
their binaries, ptxas logs, resource reports, and complete SASS.  A separate
native parser reproduced:

| walk resource | control | candidate |
|---|---:|---:|
| registers/thread | 126 | 126 |
| stack/spill/local bytes | 0 | 0 |
| ptxas user static shared | 1,792 | 1,792 |
| binary shared | 2,816 | 2,816 |
| launch dynamic shared | 0 | 8,320 |

The compiled walk uses 64,512 raw registers for two 256-thread blocks and
22,272 compiled-plus-dynamic shared bytes for two blocks.  The launch bound
compiled without spills; actual two-block residency remains pending on the
named device.

Static walk SASS contains 5,264 control instructions and 5,408 candidate
instructions.  Those totals are code size, not dynamic work.  The candidate
has exactly three fewer `CLMAD.LO`, 65 additional `LDS`, and one additional
barrier, matching the frozen source structure.  All 121 candidate `LDS`
instructions occur after its second initialization barrier; the 65-instruction
delta is therefore in the post-prologue walk body rather than the cooperative
table copy.

The producer is GitHub run `37044967938`, job `110964020023`, artifact
`11243409010`; the immutable ZIP SHA-256 is
`8d4044f1aef862396374c7961009584b3ad95451fdbf95bcbf9844dfd53669d3`.
[`compile-artifact.json`](compile-artifact.json),
[`compile-files.sha256`](compile-files.sha256), and
[`compile-result-core.json`](compile-result-core.json) retain the bindings and
core-client audit.  Independent review correctly blocked GPU dispatch because
that artifact did not contain the six preregistered CUDA probe objects; the
next exact compile artifact must add and bind them before dispatch.

Two fail-closed harness attempts remain recorded.  The first omitted a
Makefile-supplied macro and incorrectly used CUDA-13.3 CLMAD in the generic
CUDA-12.9 syntax job; neither binary was produced.  The second target compile
passed, while the Linux native replay exposed an exact C++ ABI-type mismatch
in its test-only limb array.  The repair changes that array to the reference
API's `unsigned long long` type and does not touch CUDA kernel sources.

## Next gate

The compile artifact requires an independent read-only audit, followed by the
frozen matched device replay, checkpoint, corpus, occupancy, A/A, and A/B
gates.  Until then the flag remains default-off and the decision is
`COMPILE_RESOURCE_PASS_DEVICE_OCCUPANCY_AND_RUNTIME_PENDING`.
