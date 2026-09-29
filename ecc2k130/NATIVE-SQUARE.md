# Native squaring candidate

**Accounting note (2026-09-24).** The instruction path that this note introduced
behind `PACKED_CLMAD_SQUARE=1` is now selected by `ECC_USE_CLMAD_INSN` whenever
`PACKED_CLMAD=1` compiles for sm_80+: the native square was absorbed into the
carryless product build so a fat `ARCHES="75 …"` client can keep the software
product on Turing. `PACKED_CLMAD_SQUARE` remains a baked-image / audit identity
bit and still defaults to 0; it no longer gates the asm. The measured opt-outs
are `PACKED_ALU_SQUARE` / `PACKED_ALU_SQR` (see [TWO-CHAINS.md](TWO-CHAINS.md)
and the fast-clmad jobs). Restoring the old selector would change every recent
CLMAD=1 / CLMAD_SQUARE=0 binary's Frobenius path and is not done here.

The measured RTX PRO 6000 campaign preset still reports `CLMAD_SQUARE=0`.
The source baseline for the original experiment is
`88de2c239958e5bd348376496b6ea88820269cad`.

## Mechanism and boundary

`spread32p(x)` inserts a zero bit between every input bit. For a zero-extended
32-bit polynomial `a`,

```
a(X)^2 = sum_i a_i X^(2i),    deg(a^2) <= 62.
```

Thus `clmad.lo.u64(a, a, 0)` computes exactly the same 64-bit value. Cross
terms cancel in characteristic two, and no high-half product is needed.
On sm_80+ with `PACKED_CLMAD=1`, `spread32p` uses that instruction via
`ECC_USE_CLMAD_INSN`. Host compilation and Turing slices retain the five
shift/OR/mask stages. Field reduction, inversion chains and walk rules remain
the same.

[NVIDIA's PTX specification](https://docs.nvidia.com/cuda/parallel-thread-execution/index.html#integer-arithmetic-instructions-clmad)
defines the carryless product, low-half selection and XOR accumulation.
The existing CLMAD compiler guards apply; Turing soft-falls back rather than
`#error` so mixed-fleet fat binaries stay buildable.

This does not change the generic-group boundary or the number of field
operations. At batch 16 the walk still spends `5+5/16 = 5.3125` field products
per scalar update. The ratio of field-product work to the unchanged baseline
is **1.0**. Device gains from cheaper squaring are engineering, not an
improvement to the discrete-log exponent.

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
selection and a changed native operand. These are exact source-semantic checks
under the documented instruction semantics, not CUDA compilation or device
execution.
