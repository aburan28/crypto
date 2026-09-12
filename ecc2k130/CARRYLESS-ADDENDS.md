# Optional carryless addend fusion

`PACKED_CLMAD_FUSED` selects how the packed field product uses the addend
operand of native carryless multiplication. Its default is `0`, and a
nonzero value requires `PACKED_CLMAD=1` with the existing CUDA 13.3 / sm_80
minimums.

| Value | Product construction |
|---:|---|
| 0 | Existing native product and software top-three-bit correction |
| 1 | Fold the two middle Karatsuba reconstruction words into native addends |
| 2 | Also use native addends for the top-three-bit cross terms |

CLMAD first selects the low or high 64 bits of a carryless product, then
XORs its 64-bit addend into that selected half. The implementation uses this
operation to combine product terms without materializing every intermediate
XOR. [NVIDIA PTX specification](https://docs.nvidia.com/cuda/parallel-thread-execution/index.html#integer-arithmetic-instructions-clmad).

The field type requires canonical 131-bit inputs, with only the low three
bits of the fifth word populated. Mode 2 explicitly masks both top limbs.
Its equivalence covers that field domain; arbitrary noncanonical raw words
are outside the claim. Walk denominators clear their packed jump tags before
entering the multiplier.

For a local native build, pass `PACKED_CLMAD=1 PACKED_CLMAD_FUSED=2` with the
other selected packed settings. Modal uses `ECC_PACKED_CLMAD_FUSED=2`.
The image environment, baked-binary identity, rebuild command and benchmark
metadata all retain the exact mode. Arithmetic and timed audit records must
contain a single matching mode marker before their results are accepted.

Host field arithmetic can be checked with:

```sh
make test-clmad PACKED_CLMAD_FUSED=2
```

Host execution simulates each native carryless primitive in software. It
does not execute GPU instructions. The host field suite passed for modes 1
and 2 with undefined-behavior checks enabled; seven preprocessing tests and
60 report/wrapper tests also passed.

A CPU compilation screen at B16/T256/minBlocks2 retained 122 walk registers,
98 initialization registers, and zero stack/shared/local bytes for all
three variants. Mode 2 reduced the ordinary / paired / normal-basis product
helpers from 171 / 325 / 353 to 109 / 210 / 294 non-NOP instructions. It
also increased the number of native carryless instructions. These are
compiler observations, not throughput measurements.

A separately frozen control-versus-mode-2 GPU comparison is pending. It
requires device arithmetic, full client replay, common-state, checkpoint
and exact-count checks before timing. The RTX preset continues to use mode
0 until a complete matched comparison supports a change.
