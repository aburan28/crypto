# Native squaring candidate (unmeasured)

This is an opt-in engineering experiment, not a measured throughput gain.
The measured RTX PRO 6000 preset still defaults to `CLMAD_SQUARE=0`.
The source baseline is `88de2c239958e5bd348376496b6ea88820269cad`.

## Mechanism and boundary

`spread32p(x)` currently inserts a zero bit between every input bit using five
64-bit shift/OR/mask stages. For a zero-extended 32-bit polynomial `a`,

```
a(X)^2 = sum_i a_i X^(2i),    deg(a^2) <= 62.
```

Thus `clmad.lo.u64(a, a, 0)` computes exactly the same 64-bit value. Cross
terms cancel in characteristic two, and no high-half product is needed.
The candidate replaces this helper only on the device. It affects both
`squarePolynomial131` and `sqr131`, including the latter's uses in inversion
and initialization. Field reduction, inversion chains and walk rules remain
the same. Host compilation retains the shift/mask implementation.

[NVIDIA's PTX specification](https://docs.nvidia.com/cuda/parallel-thread-execution/index.html#integer-arithmetic-instructions-clmad)
defines the carryless product, low-half selection and XOR accumulation.
The existing CLMAD compiler/architecture guards apply; the candidate additionally
requires `PACKED_CLMAD=1` and rejects square flag values outside 0/1.

This does not change the generic-group boundary or the number of field
operations. At batch 16 the walk still spends `5+5/16 = 5.3125` field products
per scalar update. The ratio of field-product work to the unchanged baseline
is **1.0**. Any device gain would be engineering through cheaper squaring,
not an improvement to the discrete-log exponent. One PTX instruction is not
a claim of one SASS instruction, one cycle, or a measured speedup.

The latest retained public-command audit in [SHARED-SIGMA.md](SHARED-SIGMA.md)
provides context, not a matched control for this candidate:

| Variant | Class | Benchmark B updates/s | DP34 B updates/s | Paired throughput ratio | Correctness |
|---|---|---:|---:|---|---|
| Existing shared-sigma public audit | engineering | 14.637530 | 14.106673 | reference only; different run | upstream GPU audit passed |
| Same-source square=0 control | engineering | unmeasured | unmeasured | pending | GPU gate pending |
| Same-source square=1 candidate | engineering | unmeasured | unmeasured | pending | host fallback checked; GPU gate pending |

The 26 B/s objective remains unachieved. No new result is added to the
index-calculus scoreboard because this change contains no GPU measurement.

## Run the validation

From `ecc2k130/`, use the same public audit with the one explicit override:

```sh
make test-clmad-square
make audit-rtx-pro6000 RTX_PRO6000_CLMAD_SQUARE=0
cp build/rtx-pro6000-audit.json build/square-control-audit.json
make audit-rtx-pro6000 RTX_PRO6000_CLMAD_SQUARE=1
cp build/rtx-pro6000-audit.json build/square-candidate-audit.json
```

These two Modal invocations may use different GPU allocations. They establish
separate correctness/audit results, not a paired speedup. For direct builds,
add `PACKED_CLMAD_SQUARE=1 PACKED_CLMAD=1` to the existing matching Make flags;
use `make -B gpu` when changing flags because Make does not track flag values.
For other Modal commands set `ECC_PACKED_CLMAD_SQUARE=1` alongside the existing
CLMAD-enabled configuration. The RTX preset override takes precedence over
the corresponding environment variable for the preset targets.

The candidate bit is propagated through the baked image, rebuild cache,
compile flags, benchmark metadata, runtime markers and device arithmetic audit.
Missing, duplicate, invalid or different markers invalidate samples. The
arithmetic test receives the same flag as the client and must pass before
the public audit proceeds to integration and timing. Its existing independent
polynomial-square and Frobenius checks exercise both affected square functions.

## Frozen acceptance rule

Before changing the default, run both binaries on one physical GPU with the
same toolkit, native target, batch 16, 256-thread blocks, minBlocks 2, 385024
workers, seeds, DP cutoff, steps and launches. Retain source/binary hashes,
GPU identity, compiler output, registers/spills and raw outputs.

1. Both device arithmetic suites and client replay/resume/guard checks pass.
   Compare cross-binary normal-basis checkpoints and sorted full DP records,
   including duplicates. Reject any state difference, mismatch or dropped DP.
2. Exclude warmups; alternate control/candidate for at least three pairs of
   complete benchmarks and three pairs of DP34 collections. Each timed sample
   must complete 201863462912 scalar updates. Count complete updates, including
   synchronization/reseed overhead, with no change to the walk population.
3. Require at least 1% improvement in the ratio of medians for **both** workloads,
   with every pair favoring the candidate. Repeat a borderline result on a
   fresh allocation. A compile failure, correctness failure, collection
   regression or failure to meet the threshold leaves the default off.

The local environment used for this change has no CUDA toolkit, GPU or Modal
runner. Preprocessing checks are not CUDA compilation; host fallback checks
are not execution of the new instruction. Device performance remains unknown.

## GPU-free symbolic validation

`python3 codegen/prove_native_square.py` checks exact Boolean-polynomial
identity for every one of the 2^32 inputs, without enumerating them. It reads
`spread32p` from the actual header, accepts only the specified native branch
(including the zero extension, identical operands, low-half selection and
zero addend), and interprets the fallback's actual shifts and masks.

The checker represents each output bit in algebraic normal form. It models
OR as `a XOR b XOR (a AND b)` rather than silently treating OR as XOR. It
independently expands the full carryless product, cancels cross terms, and
checks both equality of all 64 retained output bits and zero in every
discarded high-half bit. The result has 32 nonzero output bits.

`testclmad.py` runs this check and rejects a corrupted mask, high-half
selection and a changed native operand. All seven tests pass locally. These
are exact source-semantic checks under the documented instruction semantics,
not CUDA compilation or device execution. The existing GPU gate still applies.
