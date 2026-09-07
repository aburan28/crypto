# Optimising the elliptic-curve kernels for Blackwell

Companion to `gpu/ecc/`. This document separates three kinds of claim, and
labels every number accordingly:

- **Measured** — produced in this repository by the commands shown. All the
  measurements below are *static*: instruction counts from the emitted PTX,
  and register/stack/occupancy figures from `ptxas`.
- **Derived** — arithmetic on measured quantities plus published hardware
  parameters.
- **Expected** — a prediction that needs a GPU to confirm.

**No kernel in this directory has been run on a GPU.** There was no NVIDIA
device in the environment where it was written. The arithmetic is verified
exhaustively against a Python oracle on the CPU (`make test`), and the device
path has a self-check (`./bench selftest`) that compares every kernel against
that same host code — run it first on real hardware. Treat every throughput
figure here as a model to be falsified, not a result.

## 1. What the workload actually costs

Both kernels reduce to one primitive: multiplication in F_p for a 256-bit
prime. Everything else is bookkeeping.

| Operation | Modular multiplies | Notes |
|---|---|---|
| Field squaring | 1 | same code path as multiply |
| Field inversion (secp256k1) | 270 | 255 squarings + 15 muls, addition chain |
| Field inversion (generic) | ~316 | 4-bit window over p-2 |
| Jacobian doubling (a = 0) | 7 | 2M + 5S |
| Jacobian + affine addition | 11 | 7M + 4S |
| 256-bit scalar multiplication | ~1050 | 4-bit window: 256 dbl + 64 add + table |
| **Pollard-rho step, batch W** | **6 + 270/W** | 3 for the affine add, 3 for Montgomery's trick |

The last row is the one that matters for cryptanalysis, and it is why the
kernels batch. An unbatched affine rho step costs one inversion — 276
multiplies — while at W = 64 the same step costs about 10. That is a 27x
difference and it dwarfs every microarchitectural concern below.

At W = 8 the step costs 39.8 multiplies; at W = 16, 22.9; at W = 32, 14.4;
at W = 64, 10.2; at W = 128, 8.1. Doubling W past 32 buys progressively
less, which is what makes the memory-pressure trade-off in section 3 worth
having: the arithmetic argument for very large W is weak, so occupancy wins.

## 2. Why 32-bit limbs

Measured, with `clang++ --cuda-gpu-arch=sm_90 -O3` (see `./ptx_stats.sh`):

| Field multiply | PTX instructions |
|---|---|
| secp256k1 special reduction, portable C++ | 472 |
| secp256k1 special reduction, inline-PTX carry chains | **210** |
| Montgomery CIOS, portable C++ | 475 |
| Montgomery CIOS, inline-PTX carry chains | 352 |

The secp256k1 figure was 328 before the reduction and the product rows
were rewritten; see "Where a multiply actually goes" below.

The portable path emits no `mad.lo.cc`/`madc.hi.cc` at all: clang's NVPTX
back end decomposes the `uint64_t` accumulator into `mul.lo`/`mul.hi` pairs
plus explicit carry material (8 `selp` per multiply survive into the final
PTX). Writing the carry chain by hand removes 55% of the instructions. That
is the single largest code-level win available, and it is available on every
architecture, not just Blackwell.

Compile with `-DFP_PTX=1` to switch it on. It is off by default because
the assembly has never been executed: this was written without a GPU. It
does assemble — `ptxas -arch=sm_90/sm_100/sm_120` accepts every block —
but that only rules out syntax and register-constraint errors, not wrong
arithmetic.

Validating it is one command, though, and the reason is worth stating.
Every asm block is guarded `#if FP_PTX && defined(__CUDA_ARCH__)`,
so the host **always** compiles the portable path. `./bench selftest` built
with `-DFP_PTX=1` therefore compares device-assembly results against
host-portable results over full rho walk state — which is exactly the
differential test the assembly needs, and it now covers seven asm blocks
(`mp_add`, `mp_sub`, `mp_mac_row`, `mp_mul_row0`, `mp_mac_row9`,
`mp_add_shift32`, `mp_add_small`).

If that passes, turn `FP_PTX` on and take the 55%. It is by a wide margin
the largest win available in this code, and the only thing standing
between it and the default is one run on any CUDA device. nvcc's NVVM may already generate
better carry code than clang does here; measure both before assuming the
55% transfers.

### Where a multiply actually goes

Splitting a field multiply into its two halves shows where the time goes,
and both halves turned out to have slack:

| | portable | inline PTX, was | inline PTX, now |
|---|---|---|---|
| 512-bit product | 236 | 171 | **137** |
| secp256k1 reduction | 146 | 130 | **77** |
| full multiply | 472 | 328 | **210** |

**The reduction had no carry chains at all.** It barely moved when
`FP_PTX` was switched on (142 to 130) and so grew to 40% of the whole
multiply. Two changes fixed that:

- **Fold 1** (`T = lo + hi*977 + hi<<32`) now uses the same
  multiply-accumulate row primitive as the product, plus a new
  `mp_add_shift32` for the shifted add.
- **Folds 2 and 3** were a multiply interleaved into a carry loop. Because
  `hi2 * (2^32 + 977)` is only three limbs, each fold is really one
  straight carry chain over the low half, which is what `mp_add_small`
  now does.

**The product was zeroing an accumulator it never needed.** Schoolbook
rows were written as `t[0..9] += a * b[i]` over an 18-limb accumulator
zeroed up front — 18 `mov`s, plus two carry-propagation instructions per
row spent adding into limbs that were still zero. But a row can never
carry past its ninth limb: the most it can produce is
`(2^256 - 1) + (2^256 - 1)(2^32 - 1) = (2^256 - 1) * 2^32 < 2^288`. So the
row can *write* its top limb instead of accumulating into it
(`mp_mac_row9`), and the first row, which has nothing to add onto, is
plain `mul.lo` in its low half (`mp_mul_row0`). Each row then defines the
limb the next row's carry-out overwrites, the accumulator needs no
initialisation, and it is 16 limbs rather than 18.

`mp_mac_row` stays as it was for CIOS Montgomery multiplication, where the
accumulator genuinely is nine limbs wide at row entry.

Together: reduction 130 to 77 (−41%), product 171 to 137 (−20%), full
multiply 328 to 210 (−36%), and against the portable path the multiply is
now **55% fewer instructions**.

The portable path costs 0.9% more than it did (468 to 472), all of it in
the reduction, whose restructured folds suit a carry-flag machine rather
than a 64-bit accumulator. The row rewrite left it untouched at 236,
because clang was already eliminating the accumulator zeroing on its own
— which is exactly why that win shows up only where the carry chains are
hand-written. The regression is the right trade, the portable path being
the correctness reference and the fallback rather than the throughput
path, but it is a regression and worth knowing.

### The multiply is carry-bound, not multiply-bound

Worth knowing before optimising the wrong thing. Per squaring on the
special-reduction path, the emitted opcode mix is:

| opcode | count |
|---|---|
| add | 178 |
| shr | 102 |
| and | 100 |
| mul | 46 |
| everything else (`selp`, `cvt`, `setp`) | 18 |

Only 46 of 444 instructions are multiplies. The rest is carry emulation —
the portable path builds a 64-bit accumulator out of 32-bit shifts and
masks. That is why the inline-PTX carry chains are worth 55% while
reducing the multiply count is worth almost nothing.

Two consequences, both measured:

**A dedicated squaring routine is not worth writing.** A square needs only
the 28 off-diagonal products plus 8 diagonal ones, 36 against the
schoolbook's 64, which suggests a 44% saving. It does not materialise:
clang's common-subexpression elimination already collapses `a_i*a_j` and
`a_j*a_i` in `mul(a, a)`, so both forms emit exactly 46 multiplies. A
hand-written `mp_sqr_full` measured 427 instructions against the 439 the
portable path then cost — 2.7% — and on the inline-PTX path it was
**23.5% worse** (405 against that path's 328 at the time),
because `mul(a,a)` there uses the efficient asm carry chains while a
portable squaring routine does not. It was implemented, measured, and
reverted.

**Any further multiply-side work has to be written in the same idiom as
the carry chains it sits next to**, or it loses more on carries than it
gains on multiplies.

### What is left

The 512-bit product is now 137 instructions against a hard floor of 128 —
schoolbook needs 64 `IMAD.WIDE`-equivalent products, two PTX instructions
each — so there is nothing left there short of Karatsuba, which trades 16
of those products for several extra 128-bit adds and would very likely
lose on a carry-bound multiply.

The reduction's remaining 77 splits roughly as fold 1 (27), folds 2 and 3
(28) and the final conditional subtraction (~20). **That conditional
subtraction is the last structural win, and it is not a free one.** The
value entering it is already below 2^256 and congruent to the right
residue, so dropping it and carrying field elements in `[0, 2^256)`
instead of `[0, p)` is arithmetically sound — the reduction accepts any
512-bit input, so unreduced values feed straight back in. What it breaks
is uniqueness of representation: `eq`, `is_zero`, distinguished-point
detection and the negation map's canonical-x comparison all assume the
canonical form. Only `2^32 + 977` of the `2^256` representable values are
non-canonical, so a walk would disagree with the reference roughly once in
`2^224` steps — which is precisely what makes it dangerous to adopt
without a canonicalising step at every decision point. Worth about 10% of
a multiply; not taken.

32-bit limbs are right for every NVIDIA architecture through Blackwell for
a structural reason: the integer datapath is 32-bit, `IMAD.WIDE.U32` is the
widest single-instruction multiply, and the extended-precision carry flag
(`mad.lo.cc` / `madc.hi.cc` / `addc`) is only exposed at 32-bit granularity.
There is no 64-bit integer multiplier to address.

## 3. Occupancy is the binding constraint

Measured with `./ptx_stats.sh`, `k_rho_walk_lowmem<8>`, block size 128:

| Arch | `RHO_MIN_BLOCKS` | Registers | Stack | Resident threads/SM |
|---|---|---|---|---|
| sm_90 | unset | 224 | 792 B | 256 |
| sm_90 | 3 | 168 | 816 B | 384 |
| sm_90 | 4 | 128 | 968 B | 512 |
| sm_100 | unset | 224 | 792 B | 256 |
| sm_100 | 3 | 168 | 816 B | 384 |
| sm_100 | 4 | 128 | 968 B | 512 |
| sm_120 | unset | **255** | 800 B | 256 |
| sm_120 | 3 | 168 | 1072 B | 384 |
| sm_120 | 4 | 128 | 1248 B | 512 |

Two things fall out of this.

**Hopper and datacenter Blackwell are effectively identical.** Across all
ten kernels at three launch-bound settings, `sm_90` and `sm_100` agree on
23 of 30 register allocations and on every stack frame; the seven that
differ do so by at most 6 registers and none of them changes the resulting
occupancy. The register file is 64K 32-bit registers per SM on both and
ptxas makes essentially the same decisions. Nothing about datacenter
Blackwell changes this arithmetic, so a kernel tuned on H100 transfers to
B200 unchanged.

**Consumer Blackwell does not.** On `sm_120`, ptxas left to itself goes to
the 255-register ceiling on *every* rho kernel, and its spill frames are
consistently larger — 1248 B against 968 B at four blocks per SM. Whatever
the reason (a different spill-cost model, or different scheduling
pressure), the practical consequence is direct: **`sm_120` needs the
launch-bounds constraint more than `sm_100` does, not less.** Anyone
developing on an RTX 50-series card and deploying to B200 should tune
`RHO_MIN_BLOCKS` separately for the two.

A caveat on the last column: it assumes 2048 threads per SM, which holds
for `sm_90` and `sm_100`. Consumer parts have capped SM occupancy at 1536
threads since Ada, so on `sm_120` the same register count buys a higher
*fraction* of a smaller maximum. `./bench` prints the real limit at
startup; trust that over any table.

Left free, ptxas spends everything on registers and lands at 256 resident
threads, which is far too few to hide the dependent-issue latency of the
multiply chains. Forcing four blocks per SM doubles residency for a couple
of hundred bytes of stack. **Expected**: `RHO_MIN_BLOCKS=4` is the better
operating point everywhere; verify with `./bench rho` against the default.

Two source-level changes already cut the stack by 4-5x, and both are
measured:

| Kernel | Stack before | Stack after |
|---|---|---|
| `k_rho_walk<8>` | 4664 B | 1616 B |
| `k_rho_walk_lowmem<8>` | 4008 B | 792 B |
| `k_rho_walk_ref` | 3728 B | 536 B |

The two changes:

1. **Never `__forceinline__` a routine containing an inversion or a ladder.**
   Marking `Fp::inv`, `Curve::scalar_mul` and friends `__noinline__`
   (`FP_BIG` in `fp256.cuh`) cut compile time from over ten minutes to 46
   seconds and stopped several hundred inlined field multiplies from
   landing in every call site.

2. **Size the rare path for stack, not for speed.** `double_scalar_mul`
   builds two 16-entry Jacobian tables — 3 KB of frame — and the rho walk
   calls it only when a walk reseeds, a few times per million steps. ptxas
   sizes every thread's frame for the worst path through the kernel, so
   that 3 KB was charged to every resident thread permanently.
   `double_scalar_mul_small` does the same job with a table-free binary
   ladder at half the speed and 3 Jacobian points of frame.

The same reasoning gave `inv_addchain_secp256k1`: the libsecp256k1 addition
chain needs 270 operations against the window method's 316, and 11 live
field elements against a 16-entry table, so it is both faster and 512 bytes
lighter. It is validated against the Python oracle by `make test`.

## 4. Blackwell-specific opportunities

Ordered by expected return. Confirm each device parameter with
`./bench` (it prints SM count, shared memory per SM and per block, and L2
size at startup) rather than trusting a table.

### 4.1 Fit the jump table in shared memory, then make it bigger

`k_rho_walk` stages the r-adding jump table in shared memory. An entry is
17 words (x, y, inf) = 68 bytes, so 2^r_bits entries cost:

| r_bits | Table size |
|---|---|
| 8 | 17 KB |
| 10 | 68 KB |
| 12 | 272 KB |

Datacenter Blackwell (`sm_100`) carries the Hopper-class 228 KB of shared
memory per SM, against 164 KB on Ampere and 100 KB on consumer parts
(`sm_120`), so it can hold an r = 2^11 table resident where Ampere cannot.
Bigger tables matter twice over: the walk is closer to random, and — the
reason it is worth real shared memory — **fruitless cycles scale as 1/R**.
With the negation map on, a 2-cycle forms with probability about 1/(2R) per
step, so going from r_bits = 8 to 11 cuts the escape rate eightfold. The
cycle counter in `bench rho` reports the measured rate.

**The `inf` word is load-bearing.** A 17-word stride is odd, hence
invertible mod 32, so 32 lanes reading 32 different table entries hit 32
different banks. Packing the entry down to 16 words to "save space" makes
the stride even and costs a 16-way bank conflict on the hottest read in the
kernel. This is noted in `kernels.cuh` because it is exactly the kind of
optimisation that gets "fixed" later.

### 4.2 Share one table across a thread-block cluster

Thread block clusters and distributed shared memory arrived with Hopper
(`sm_90`) and Blackwell keeps them. A cluster of 4 blocks can address one
another's shared memory, so a 68 KB table can live once per cluster instead
of once per block. **Expected**: this is what makes r_bits = 11-12 practical
without destroying occupancy, since shared memory per block is what limits
blocks per SM. Not implemented here — `rho_stage_table` copies per block —
and it is the first thing to try if `bench rho` shows occupancy limited by
shared memory rather than registers.

### 4.3 Stage the table with TMA instead of a copy loop

`rho_stage_table` is a strided `for` loop over `blockDim.x`. On `sm_90` and
later, `cp.async.bulk` (the Tensor Memory Accelerator) moves a contiguous
block into shared memory with one instruction and no register traffic.
The table is loaded once per kernel launch, so this only matters for short
launches — but the rho kernel is often launched with small `iters` to keep
DP latency low, and then the staging cost is real. Measure by comparing
`bench rho --iters 16` against `--iters 1024`.

### 4.4 Do not reach for the tensor cores

Blackwell's headline is its fifth-generation tensor cores with FP4/FP6 and
tensor memory. They are the wrong instrument here. Big-integer multiplication
via INT8 tensor-core MMA is a real technique, but it needs the operands
split into 8-bit digits, accumulates in INT32, and pays a lane-permutation
and carry-propagation tax on the way out. For a 256-bit modulus the digit
count is small enough that the setup dominates: the schoolbook needs 64
32x32 products, which the IMAD pipe does directly. The tensor path becomes
interesting at 2048 bits and beyond (RSA, FHE), not at 256.

Concretely: the FP4/FP6 formats and TMEM add nothing to an integer modular
multiply, and the INT8 IMMA path is not obviously faster than IMAD.WIDE at
this size. Spend the effort on section 2 instead.

### 4.5 Distinguished-point output

DP emission uses a global `atomicAdd` on one counter. At a 2^-20 DP rate
with a few million walks, that is a handful of atomics per launch and costs
nothing. If the DP rate is raised (short launches, small search spaces), a
warp-aggregated atomic — one `atomicAdd` per warp using `__ballot_sync` and
a leader lane — removes the contention. Not implemented; it is a five-line
change in `rho_emit_dp` and only worth doing if a profile shows it.

### 4.6 The persistent-kernel question

The walk state lives in global memory in structure-of-arrays layout, so a
warp reading limb *l* of its 32 walks touches 128 contiguous bytes — a
single sector-aligned transaction. Per step, a walk reads and writes x and y:
128 bytes of traffic for roughly 10 modular multiplies, or about 13 bytes
per multiply. **Derived**: at a few hundred SASS instructions per multiply
this is nowhere near the bandwidth roofline, so keeping the state in global
memory is the right call and a persistent kernel holding state in registers
would buy little. The `_lowmem` variant deliberately reads the state a
second time for exactly this reason: the extra load is cheap, the saved
registers are not.

## 5. Reproducing the static measurements

The environment used had no CUDA toolkit, so the numbers above came from
clang's CUDA support plus a standalone `ptxas`:

```bash
pip install nvidia-cuda-nvcc-cu12 nvidia-cuda-runtime-cu12 \
            nvidia-cuda-cccl-cu12 nvidia-curand-cu12
# assemble a CUDA_PATH from the wheels (include/, nvvm/libdevice, bin/ptxas)
clang++ --cuda-device-only --cuda-path=$FAKE --cuda-gpu-arch=sm_90 \
        -Wno-unknown-cuda-version -O3 -std=c++17 -I. \
        -DGPU_ECC_CURVE_HEADER='"curve_secp256k1.h"' -S -o inst.ptx inst.cu
# clang 18 tops out at sm_90; retarget the PTX and let ptxas 12.9 do Blackwell
sed -e 's/^\.version .*/.version 8.7/' -e 's/^\.target sm_90.*/.target sm_100/' \
    inst.ptx > inst_100.ptx
ptxas -arch=sm_100 -O3 -v inst_100.ptx -o /dev/null
```

`ptxas` from the CUDA 12.9 wheel accepts `sm_100`, `sm_103` and `sm_120`,
so register and occupancy figures for Blackwell are obtainable without the
hardware. Instruction *scheduling* is not — for that you need a device.

## 6. First things to do on real hardware

1. `./bench selftest` — nothing else means anything until this passes,
   especially with `-DFP_PTX=1`.
2. `./bench field` — the measured cost of one multiply and one inversion.
   Everything in section 1 is denominated in these.
3. `./bench rho --variant reg` vs `--variant lowmem` vs `--variant ref`,
   sweeping `--w 8 16 32` and `RHO_MIN_BLOCKS`. The model says lowmem at
   W = 32 with four blocks per SM; the model has not been tested.
4. `./bench rho --neg 1 --rbits 8` vs `--rbits 11`, watching the reported
   cycle-escape percentage. This is where the shared-memory capacity of
   `sm_100` should show up as a real advantage over `sm_120`.
