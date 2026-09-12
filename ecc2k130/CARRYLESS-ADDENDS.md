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
| 3 | Use native low products with software high carries and a software three-bit product |

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

## Mode-2 GPU screening result

The [controlled GPU screen](benchmarks/carryless-addends/comparison.json)
passed device arithmetic, full client replay, common-state, checkpoint and
exact-count checks on one RTX PRO 6000 Blackwell Server Edition. Every
timed row completed 201,863,462,912 complete scalar updates.

| Screening order | Mode | Complete scalar rate B/s |
|---:|---:|---:|
| 1 | 0, control | 13.363256 |
| 2 | 2, native top terms | 10.280322 |
| 3 | 0, control | 13.280746 |

Mode 2 failed the predefined qualification threshold. No repeated
confirmation or timed DP34 collection was triggered. Two warmups are
excluded from this table. The [retained-result audit](benchmarks/carryless-addends/comparison-review.json)
passed all five timed rows and 20 checkpoint children, including two
expected worker-geometry rejections.

The shorter native code did not improve this workload's throughput. This
screen supplies no hardware-pipe ceiling measurement. Mode 2 is not selected
by the RTX preset, and the implementation PR remains a draft while a hybrid
with fewer native top-term multiplies is evaluated. The 15 B/s target remains
unachieved.

## Mode-3 hybrid result

Mode 3 reduces the ordinary field product to ten native carryless
instructions (seven low, three high), with software operations for the
two-bit high carries and the three-bit top product. Its ordinary / paired /
normal-basis helpers compiled to 134 / 249 / 323 non-NOP instructions at
the same 122 walk registers and zero spills.

The [fixed-geometry GPU screen](benchmarks/carryless-addends/hybrid-comparison.json)
passed both arithmetic suites, full client checks, common-state comparisons
and all 20 checkpoint children. Each of the five timed rows completed
201,863,462,912 scalar updates. Screening measured 13.486482 B/s for the
first control, **12.601867 B/s for mode 3**, and 13.414589 B/s for the final
control. The candidate failed qualification, so no repeated confirmation or
timed DP34 collection followed. The [artifact audit](benchmarks/carryless-addends/hybrid-comparison-review.json)
passed.

Mode 3 remains opt-in and is not selected by the RTX preset. A separate
dependency assessment tested independent native products before combining
their results. Its measured outcome is recorded below.

## Independent-product follow-up

The retained [prototype header](benchmarks/carryless-addends/independent-top-prototype.h)
adds experimental modes 4 and 5. They are not exposed by the public mode
selector. Mode 4 restores independent Karatsuba leaves; mode 5 also computes
all four top-term native products with zero addends before combining them.
PTX and SASS analysis found maximum native multiply dependency depths of
three and one, respectively, versus four in mode 3.

The [GPU comparison](benchmarks/carryless-addends/independent-comparison.json)
passed arithmetic, client, normalized-state and all 36 checkpoint children,
including three expected worker-geometry rejections. Every timed sample
completed 201,863,462,912 scalar updates.

| Screening order | Variant | Complete scalar rate B/s |
|---:|---|---:|
| 1 | Control | 13.384979 |
| 2 | Independent core, mode 4 | 12.411916 |
| 3 | Independent core and top products, mode 5 | 12.232020 |
| 4 | Control | 13.327505 |

Neither candidate qualified. Three warmups are excluded from the table;
no repeated confirmation or timed collection was triggered. The
[artifact audit](benchmarks/carryless-addends/independent-comparison-review.json)
passed all seven timed rows. Reducing dependency depth did not improve
throughput in this screen. The tested native top-term constructions are
therefore retained as unsuccessful experiments, with no RTX preset change.
