# Hybrid sparse tensor raw-product probe

This isolated experiment computes the integer convolution used to reconstruct
a **raw 128-bit carryless product**. It is not the production walk multiplier,
does not include the top three field bits or GF(2^131) reduction. The separate
capacity panel below measures repeated matrix sequences, not complete fields
or walk updates.

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

## Measured matrix-sequence capacity

`bench.cu` reuses the validated packing and matrix helpers. It runs four or
eight independent accumulator chains for 1,024 rounds, retaining every
component. Each CUDA-event interval includes device fragment loads, seed
initialization and output stores. CPU packing, readback and verification are
outside the interval; full CPU verification separates successive launches.

| Mode | Median B raw-core equivalents/s | Range across five measurements |
|---|---:|---:|
| Dense, 4 chains | 35.859551 | 35.640470–36.005424 |
| Hybrid, 4 chains | 52.982523 | 52.746626–53.516437 |
| Dense, 8 chains | 36.205183 | 36.101639–36.282799 |
| Hybrid, 8 chains | 53.819402 | 53.688071–54.316710 |

The ratios of hybrid to dense mode medians are 1.477501 and 1.486511.
Medians of the five individual paired ratios are 1.479964 and 1.489431;
all ten pairs favored the hybrid. Eight warmups are excluded. These are
per-launch rates under the stated verification cadence, not uninterrupted
sustained capacity or universal hardware ceilings.

Each launch uses 192,512 threads, hence 6,016 warps. The exact denominators
are 24,641,536 or 49,283,072 raw-core equivalents. Dense executes three
matrix instructions per equivalent; hybrid executes one dense and one
sparse instruction. Counts use warps, not individual lanes. All
129,368,064 accumulator outputs across the 28 launches were checked against
initial C plus 1,024 times an independent integer convolution. The reviewed
maximum accumulated value is 1,124,139,008, below signed 32-bit overflow.

The current walk schedule would need `15B * 165/32 = 77.34375B` raw cores/s
at 15B complete updates/s, even before its other work. The measured
hybrid-eight median is 69.5847% of that demand. This comparison does not
prove that the target is impossible; it means the measured sequence does
not yet demonstrate the required capacity for replacing all those cores.
No complete-field or whole-walk rate is measured here.

Full source, native-code and runtime audits are retained in
`evidence/capacity-source-review.json`, `capacity-code-review.json` and
`capacity-gpu-review.json`. The executed loops contain exactly the expected
unpredicated MMAs, distinct feedback chains and complete stores, with no
local/shared memory or spills. The original correctness kernels retain
their exact validated instruction encodings.

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

To reproduce the reviewed capacity panel:

```sh
python3 -m unittest test_capacity_gate.py test_assertion_mode.py
modal run compile_capacity.py
modal run run_capacity.py
```

The capacity compiler checks a fresh build against the retained source,
complete native instruction encodings and resources before creating a
binary-bound review under `build/`. A mismatch stops before GPU timing.
`run_capacity.py` requires that review and first reruns the complete
16,512-case correctness fixture. It then executes the fixed paired panel,
checks every output and records exact counts and event intervals. Fresh
results are written to `build/capacity-gpu-result.json`; historical evidence
remains unchanged.

The Python fixtures use assertions for verification. Every entry point and
reference/admission module rejects `-O`, `-OO` and `PYTHONOPTIMIZE` before
importing dependencies or setting up remote work, so disabling assertions
cannot silently produce an accepted result.

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
