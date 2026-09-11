# Hybrid sparse tensor raw-product probe

This isolated experiment computes the integer convolution used to reconstruct
a **raw 128-bit carryless product**. It is not the production walk multiplier,
does not include the top three field bits or GF(2^131) reduction, and has no
throughput result yet.

The checked-in CUDA source passed on one RTX PRO 6000 Blackwell Server Edition:
16,384 bit-basis pairs, 64 edge pairs and 64 deterministic dense pairs. Both
implementations check every integer convolution coefficient and reconstruct
the raw carryless product against a separate bit-by-bit reference. In total,
4,227,072 coefficients and 33,024 raw products were compared. The test performs
no timing.

| Compiled kernel | Native matrix instructions | Registers | Stack/spills |
|---|---|---:|---:|
| Dense control | 3 × `IMMA.16832.U8.U8` | 28 | 0 |
| Hybrid | 1 × `IMMA.16832.U8.U8` + 1 × `IMMA.SP.16864.U8.U8` | 24 | 0 |

These are actual CUDA 13.3.73 `sm_120` compiler observations. A sparse
instruction need not have the same issue cost as a dense instruction, so
the count change is not a measured speedup.

## Reproduce

From this directory, using an authenticated Modal client:

```sh
python3 check_mapping.py
modal run compile.py
modal run validate.py
```

`compile.py` performs one CPU-only CUDA compile using a pinned CUDA 13.3.1
image, retains full compiler/disassembly output and executable hashes in
Modal Volume, and downloads a verified receipt to `build/compile-result.json`.
It does not execute the binary. `validate.py` loads that exact binary, verifies
its source and binary hashes, disables PTX JIT, and runs the fixed correctness
panel on one RTX PRO 6000. Its receipt is `build/gpu-result.json`.

The `evidence/` directory contains immutable results from the original
validated run. Local driver paths in those historical records describe that
run; the portable drivers write fresh results under ignored `build/`.

## Mapping

Each input byte represents two bits as `a[2*i] + 128*a[2*i+1]`. The 128
output matrix elements hold integer convolution coefficients. Extracting
three parity positions recovers the raw carryless result, with one bit-128
correction when both original inputs are all ones. Every coefficient is
bounded by 1,065,024, so signed 32-bit accumulation is exact in this fixture.

The dense control uses K slices 32–63, 64–95 and 96–127. The middle slice is
kept dense. For the remaining slices, 17 always-active columns pair with
17 always-zero columns; 15 pairs of moving edge columns are complementary.
One global permutation gives exactly two designated A positions in every
four-column group in every row. One sparse K64 operation then replaces two
dense K32 operations. This uses NVIDIA's
[structured sparse MMA format](https://docs.nvidia.com/cuda/parallel-thread-execution/index.html#warp-level-matrix-instructions-mma-sp).

Metadata follows the documented warp layout, which differs from A-value
ownership. With `g = lane >> 2` and `t = lane & 3`, metadata nibble `j` belongs
to row `g + 8*(t & 1)` and group `8*(t >> 1) + j`. Indices within each nibble
are sorted, and the sparse selector is zero. Both matrix instructions are
executed uniformly by all 32 lanes, and all four accumulators per lane are
stored and checked.

The source and packing were independently reviewed against the PTX layouts
and [CUTLASS sparse MMA implementation](https://github.com/NVIDIA/cutlass/blob/main/include/cutlass/arch/mma_sparse_sm80.h).
The original CUDA source SHA256 is
`6836bfef024c2f635e5f25018c19bea842ad0bdd80a2466190c7883aa2e441a9`.

## Scope

Input matrix preparation is on the CPU in this correctness fixture, and raw
parity reconstruction is also checked on the CPU. Neither is charged as GPU
field work. A useful performance experiment must account for packing,
reconstruction, field reduction and the complete walk before claiming a
walk-rate improvement. This probe does not establish 15 billion walk
iterations per second or any other full-field/whole-walk rate.
