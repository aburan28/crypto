# The top three bits of every product, on the carryless unit

`PACKED_TOP_CLMAD=1` moves the 3-bit top-word cross terms of `product131`
from masked shifts onto `clmad`. Off by default. It requires `PACKED_CLMAD=1`.
SASS counted −5.4% ALU / +41% clmad and the pipe model said about +4%. On
one RTX PRO 6000 the documented Make pair ran **15% slower** than the
matched shipping product. The knob stays off.

## Why this product term, measured

`P131` holds 131 bits in five words. The 4x4-word `clmul128` covers 128 of
them with six `clmad`s; the remaining three bits of each operand contribute
`a4*B_lo + A_lo*b4 + a4*b4*x^256`, which the shipping code builds as six sign
masks, twelve masked ANDs and ten funnel-shift accumulations.

Compiled with nvcc 13.3.73 for `sm_120` at the audited RTX PRO 6000 preset,
counted with `cuobjdump` (`./sass_cost.py`, less the one-instruction kernel
overhead):

| routine | ALU | clmad |
|---|---:|---:|
| `product131`, as shipped | 77 | 6 |
| `product131` with the top-word correction removed (wrong, for the count) | **12** | 6 |
| `reducePolynomial131` | 77 | 0 |

So 65 of the product's 77 ALU instructions are the correction for three bits,
and the four-word carry-less product itself is 12 — ptxas already absorbs the
Karatsuba glue, so there is nothing left there for `clmad`'s free addend to
take. At 5.3125 products per update that correction is about 345 of the
2,189.75 instruction visits the benchmark receipts count per update: **15.8%
of the walk, the largest single addressable item in it.**

## What the variant does

`a4*B_lo` is two 3x64-bit products whose low 64 bits each fit one `clmad.lo`,
and `clmad`'s third operand is XORed in for free, so the accumulate into the
product costs nothing:

    c[4..5] = clmad.lo(a4, B0, clmad.lo(b4, A0, c[4..5]))

What does not fit is the two bits each 3x64 product pushes past bit 63; those
are `(a4 * (B0 >> 62)) >> 2`, a 3x2-bit product done as masked shifts. With
`a4*b4` that residue is about 30 ALU.

The device path is inline `clmad.lo.u64` with a register addend. The host
path emulates it through the existing software `clmul64`, so
`make test-packed-network` checks the device formulation bit for bit against
the independent reference: every unit-coefficient pair, 160 dense cases, the
paired product, squaring, inversion and the Frobenius powers.

## Measured in SASS

Same compiler, same preset, `ECC_PACKED_TOP_CLMAD` 0 against 1:

| routine | ALU 0 | ALU 1 | clmad 0 | clmad 1 |
|---|---:|---:|---:|---:|
| `product131` | 77 | **52** | 6 | 10 |
| `mulPolynomial131` | 164 | **128** | 6 | 10 |
| `mulPolynomialPair131` | 315 | **279** | 14 | 20 |
| `mul131` (inverse chain) | 347 | 320 | 6 | 10 |
| `inv131` | 1853 | 1824 | 26 | 30 |

`mulPolynomial131` loses 36, not 25: with the shifted copies gone ptxas fuses
more of the reducer that follows.

Per scalar update at batch 16, over the routines `sass_cost.py --update`
weights: **1,621 -> 1,534 ALU slots (-5.4%)**, 29.1 -> 41.2 clmads (+41%).
Against the 2,189.75 visits per update, about -4.0%.

The walk kernel's registers fall from 104 to 94 (the masks and shifted copies
are gone), static SASS from 4,048 to 3,992. Nothing about occupancy moves.

## What the pipe model predicts

Nsight Compute put the ALU pipe at 87.3% and the FP64 pipe that carries
`clmad` at 51.4% at this preset. Scaling the carryless load by 41.2/29.1 puts
that pipe at 73% at today's rate; if the ALU cut converts to rate, about 4%
more, at 76%. Both stay under the 87% the ALU pipe already sustains, so the
model does not predict the carryless unit binding, and its estimate is the
full ALU saving: **about +4%**.

Two ways that is wrong. A product issues its ten `clmad`s close together, and
a unit at 76% average may stall on that burst. And 30 ALU of residue scheduled
between `clmad`s may cost more than 30 slots. Either shows up only on a card.

## Measuring it

Alternate on one allocation; the two builds differ by one define.

```sh
make bench-rtx-pro6000
make bench-rtx-pro6000 RTX_PRO6000_TOP_CLMAD=1
```

The benchmark prints `packed top clmad: N` and the runner rejects a rate
whose printed identity differs from the requested build, as it does for
every other packed knob.

## Measured on a card

Two successive Modal allocations, both NVIDIA RTX PRO 6000 Blackwell
Server Edition, driver 580.95.05, CUDA 13.3.1 image, same shipping
preset (batch 16, 256 threads, minBlocks 2, 385,024 workers, 1,024
steps, 32 launches, 3 repeats). Each sample completed 201,863,462,912
scalar updates with zero drops. Identity matched the requested build
on every sample: control printed `packed top clmad: 0` at 104
registers; candidate printed `packed top clmad: 1` at 94 registers,
the SASS count. Different GPU UUIDs; this is not an interleaved pair
on one allocation.

Receipts: [benchmarks/top-clmad/summary.json](benchmarks/top-clmad/summary.json),
[control.json](benchmarks/top-clmad/control.json),
[candidate.json](benchmarks/top-clmad/candidate.json).

| variant | median B/s | samples | / control | / 14.637530 | / 23 B floor | / 25 B floor | correctness | class |
|---|---:|---|---:|---:|---:|---:|---|---|
| shipping product (matched control) | **15.115792** | 15.205241, 15.115792, 15.026171 | 1.000 | 1.033 | 0.657 | 0.605 | top clmad 0 | reference |
| `PACKED_TOP_CLMAD=1` | 12.859376 | 12.820670, 12.859376, 12.868661 | **0.851** | 0.879 | 0.559 | 0.514 | top clmad 1 | engineering |

Unit is billions of completed scalar updates per second. The one-add
floor is the 23–25 B/s bound from [THROUGHPUT-30B.md](THROUGHPUT-30B.md).
The 14.637530 column is the public shipping audit
([SHARED-SIGMA.md](SHARED-SIGMA.md)); the ratio that decides the knob
is the matched control from this pair. Every pair was slower. The
predeclared 1% acceptance rule fails. Default remains 0.

The SASS ALU cut happened (104 → 94 registers, identity gate on). The
pipe model's +4% did not. The note already named the two ways that
prediction fails: ten `clmad`s issued close together can stall a unit
that is only 76% on average, and the 30 ALU of residue scheduled
between them can cost more than 30 slots. The card chose one of those.
This is not an advance against the one-add floor; the ratio fell.

DP34 collection was not run. The complete-scalar bench already
falsifies promotion; a second workload cannot loosen that rule.

## What was ruled out on the way

Each of these was costed before this variant was chosen, and each is dead for
a reason that is a number.

- **A general linear map by table lookup**, which
  [THROUGHPUT-30B.md](THROUGHPUT-30B.md) names as lever 2: nibble-indexed
  tables need 33 lookups of 5 words per application, two per update. Every
  lane has its own `j`, so one LDG gathers from about 16 cache lines; at
  14.6 B updates/s that is ~5.6 L1 cycles per SM-clock against one available.
  Shared memory is no better: random rows across 32 banks cost ~3 cycles an
  LDS, and 10.3 LDS per update x 3 is 31 cycles against 29.6 per update.
  The vector4 calibration that made memory look free was 42 coalesced
  instructions, not 330 gathers.
- **The sparse pentanomial basis**, lever 1: it pays only with lever 2,
  because in that basis the map to normal-basis coordinates for the class
  weight is a dense matrix with no butterfly structure.
- **A five-word Frobenius network**: only 2 of the 64 emitted operations
  touch padding words 5-7; the generator already skips the rest.
- **A shifted modulus** `f(x+1)`, reachable by a cheap Sierpinski butterfly:
  45 terms against 19.
- **Hoisting compact-state address arithmetic**: the SASS is one
  `IMAD.WIDE` and one `IADD.64` per access with the per-thread base already
  outside the step loop.
- **`clmad`'s free addend on the base product**: 12 instructions left in the
  whole 4x4-word product; nothing to absorb.

## Getting the compiler without root

`sass_cost.py` needs nvcc 13.3 or newer for `clmad`, and neither the pip
wheels nor the Ubuntu package carry one (`nvidia-cuda-nvcc-cu13` on PyPI is a
1 KB placeholder). NVIDIA's apt repository serves plain `.deb` files, and a
`.deb` is an `ar` archive:

```sh
R=https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2404/x86_64
for p in cuda-nvcc-13-3 libnvvm-13-3 cuda-crt-13-3 cuda-cudart-dev-13-3 \
         cuda-cuobjdump-13-3 cuda-nvdisasm-13-3; do
  f=$(curl -s $R/Packages | awk -v P=$p '/^Package: /{p=$2} /^Filename: /{f=$2} /^$/{if(p==P)print f;p=""}' | sort -V | tail -1)
  curl -sO $R/$f && ar x $(basename $f) && tar -xf data.tar.*
done
./sass_cost.py --nvcc usr/local/cuda-13.3/bin/nvcc --cuobjdump usr/local/cuda-13.3/bin/cuobjdump
```

About 40 MB. Every SASS number above came from that.
