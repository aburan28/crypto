# gpu/ecc2k — GPU kernels for Koblitz curves and ECC2K-95

CUDA kernels for binary-field elliptic curves, aimed at ECC2K-95: the
Certicom Koblitz challenge

```
E_0 : y^2 + x y = x^3 + 1     over  F_{2^97},  f(t) = t^97 + t^6 + 1
```

whose group has order 4·r with r a 95-bit prime. It was solved in 1998 by
Harley and collaborators, so it is a realistic retired target rather than
an attack on anything live. Two toy Koblitz curves (21-bit and 40-bit
subgroups) ship alongside it so the whole pipeline can be run to a real
discrete logarithm in seconds.

Sibling directory `gpu/ecc/` does the same job for secp256k1 over a prime
field. Reading the two together is the point: the same attack, two field
arithmetics, and a very different balance of costs.

**No kernel here has been run on a GPU.** The arithmetic is verified against
a Python oracle on the CPU (`make test`), including an end-to-end discrete
log, and the device path has a self-check (`./bench2k selftest`) to run
first on real hardware.

## Why Koblitz curves get their own kernels

On a Koblitz curve the Frobenius map τ(x, y) = (x², y²) is a group
endomorphism, and it costs two squarings. So from any point P the 2m points

```
{ ± τ^i(P) : 0 ≤ i < m }
```

are all reachable for almost nothing. A rho search that walks on these
*equivalence classes* rather than on points finishes in √(2m) fewer steps:
13.9× for m = 97. That is why the ECC2K challenges fell before their
prime-field counterparts of similar size.

Making the walk descend to classes requires an iteration function that
commutes with the class action. This uses the ECC2K-130 shape:

```
j(P) = ((g(x_P) / 2) mod 8) + 3
P   -> P + τ^j(P)
```

where `g` is the **normal-basis Hamming weight** of x. Squaring x rotates
the normal-basis coordinates, so g is constant on Frobenius orbits, giving

```
f(τP) = τ(f(P))        because j(τP) = j(P)
f(-P) = -f(P)          because x(-P) = x(P)
```

A distinguished point is `g(x) ≤ threshold`, also class-invariant, and it
is reported as the canonical class representative: the smallest x among the
m Frobenius images. Negation needs no handling in the representative
because −(x, y) has the same x.

**The `/2` is load-bearing.** `g` is even for every point of the curve —
checked exactly over all 45,562 classes of the m = 23 subgroup, and by
sampling at m = 97. Drop the halving and `g mod 8` only ever takes the
values {0, 2, 4, 6}: the walk runs on four branches rather than eight, and
collides later. Measured on the 40-bit curve, 64 walks, dp `g ≤ 11`, 160
trials: 187,831 steps for the four-branch rule against 151,005 for this one.
This is the rule Bailey et al. specify for ECC2K-130, and the one
`ecc2k130/include/walk.h` implements.

Arithmetic stays in the polynomial basis, where multiplication is
practical; only `g` needs the normal basis, and that is a fixed linear map
applied through a 4-bit-window table. m = 97 has no optimal normal basis
(2m+1 = 195 is composite, unlike ECC2K-130's 2·131+1 = 263), which is
exactly why the ONB-based bitsliced approach used against ECC2K-130 does
not transfer here.

A step multiplies rather than adds coefficients: P → (1 + τ^j)P, and on the
prime-order subgroup τ acts as multiplication by an integer s with
s² + s + 2 ≡ 0 (mod r). So a walk at aG + bQ moves to
(a(1+s^j), b(1+s^j)). The host tracks that; the device tracks only points.

## Build and test

```bash
make test         # generate headers, build and run three CPU test suites
make bench        # CUDA benchmark binary (needs nvcc)
make bench ARCH=sm_100
```

| Suite | Curve | Subgroup | What it covers |
|---|---|---|---|
| `test_ecc2k95` | F(2^97), a = 0 | 95-bit | the challenge curve |
| `test_k23` | F(2^23), a = 0 | 21-bit | end-to-end solve, 24 trials |
| `test_k41` | F(2^41), a = 0 | 40-bit | end-to-end solve with real DP filtering |

## How this is verified

Every header is host/device dual-compiled, so the CPU harness exercises the
code the kernels run, not a model of it.

1. **Against Python.** Field vectors (add, mul, sqr, inv, τ³), point vectors
   (add, double, scalar mul, τ), 128-step walk traces with the per-step
   Frobenius exponent, and canonical class representatives.
2. **Structural identities the oracle is not needed for.** τ is a field and
   group homomorphism, τ^m is the identity, negation preserves x, and
   `g` is invariant under all m Frobenius images of every test vector.
3. **The identity the whole attack rests on**: τ(G) = s·G, checked both
   against Python and by computing s·G on the device path.
4. **Class equivariance**: f(τP) = τf(P) and f(−P) = −f(P) on sample points.
5. **Three steppers agree.** The batched, low-memory and unbatched
   implementations must produce bit-identical state over 300 iterations and
   the same multiset of distinguished points.
6. **End to end, and quantitatively.** The toy curves are solved through the
   full pipeline — walks, distinguished points, host replay, class
   relation, linear solve — and the recovered logarithm is verified.

That last one is worth showing, because a single solve proves correctness
but says nothing about whether the Frobenius speedup is real. Averaging
over many instances does:

```
[solve] 21-bit Koblitz DLP, 128 walks, 24 trials, dp threshold g <= 23 (rate 1/1)
  mean 251 steps over 24 solves
  sqrt(pi/2 * r/(2m)) = 268  -> measured 0.94x   [class walk]
  sqrt(pi/2 * r/2)    = 1283 -> measured 0.20x   [negation only]

[solve] 40-bit Koblitz DLP, 128 walks, 8 trials, dp threshold g <= 11 (rate 1/463)
  mean 191840 steps over 8 solves
  sqrt(pi/2 * r/(2m)) = 145129 -> measured 1.32x   [class walk]
  sqrt(pi/2 * r/2)    = 929276 -> measured 0.21x   [negation only]
```

The 21-bit run lands on the class-walk prediction. The 40-bit run sits 32%
above it, and that excess is accounted for: 128 walks each need about 463
steps to reach a distinguished point, so 59k steps of the 192k are the DP
tail, and 145k + 59k ≈ 204k.

Rho collision times have a standard deviation about equal to their mean, so
24 and 8 trials resolve these to roughly ±20% and ±35%. They confirm the
√(2m) speedup; they are too noisy to rank iteration functions. The
four-versus-eight-branch comparison quoted earlier used 160 trials against a
fixed collision criterion rather than the full solve pipeline.

## What the field arithmetic actually costs

Measured statically (`./ptx_stats2k.sh`, clang for sm_90):

| Operation | PTX instructions | Notes |
|---|---|---|
| F(2^97) multiply | 405 | 9 carry-less 32×32 products, Karatsuba twice |
| F(2^97) squaring | 79 | 0.20 of a multiply — bit spreading, no multiplier |
| secp256k1 multiply | 472 | for comparison, from `gpu/ecc`; 210 with its inline-PTX carry chains |

**A 97-bit binary field multiply costs about as much as a 256-bit prime
field multiply here** — because this backend multiplies in software. The
32×32 carry-less product is synthesised from sixteen ordinary widening
multiplies with the interleaved-mask trick — split each operand into four
subsets by bit index mod 4, multiply as integers, mask. Within a subset the
per-bit partial sums cannot exceed 8, which fits in the 4-bit gap between
kept bits, so no carry ever crosses into a bit that matters.

**That is an artefact of the emulation, not a property of GPUs.** The
hardware instruction exists: PTX ISA 9.3 defines `clmad.lo.u64` /
`clmad.hi.u64` for `sm_80` and later, and CUDA 13.3 or newer emits it.
`ecc2k130/NATIVE-CARRYLESS.md` measured **+22.4%** on a complete ECC2K-130
walk by switching to it, with the inline asm in
`ecc2k130/include/packed131.h`. This backend has not been ported. Static
counts here (clang, sm_90) put ~83% of a 128-bit carry-less product in the
software `clmul32` emulation and ~55% of a squaring in the generic two-pass
reducer, so both are worth more here than the 22.4% measured there —
`ecc2k130/codegen/gendirectreduce.py` already generates the specialised
reducer.

What Koblitz curves give back:

- **Squaring is a fifth of a multiply**, so τ is nearly free. This is the
  whole basis of the 13.9× class speedup.
- **Inversion is cheap.** Itoh–Tsujii computes a^(2^m−2) with about 7
  multiplications and 96 squarings, roughly 23 multiply-equivalents,
  against 270 for a Fermat inversion modulo a prime.

That last point changes the kernel's shape. Batching inversions across
walks is worth 27× over a prime field; here it is worth 3–4×, and W beyond
16 buys little. Counted at m = 97 with the distinguished-point test
disabled, in multiply-equivalents at 1 squaring = 0.2 multiplies:

| Batch W | Per rho step | Breakdown |
|---|---|---|
| 1 | 30.95 | 9.00 mul + 109.75 sqr |
| 2 | 19.31 | 7.00 mul + 61.56 sqr |
| 4 | 13.51 | 6.00 mul + 37.55 sqr |
| 8 | 10.62 | 5.50 mul + 25.60 sqr |
| 16 | 9.17 | 5.25 mul + 19.62 sqr |
| 32 | 8.44 | 5.12 mul + 16.59 sqr |

(Per step: 2 multiplies for the addition, 1 squaring for λ², 3 multiplies
for Montgomery's trick, 2j squarings for the two Frobenius chains, one class
weight, and one inversion — 97 squarings and ~7 multiplies — split W ways.)

**This table counts arithmetic only, and arithmetic is not what picks W.**
`ptxas` puts `k2k_rho_walk<W>` at 344, 464, 656, 1040 and 1816 bytes of
per-thread local memory for W = 1, 2, 4, 8, 16 (sm_90, block 128,
`R2K_MIN_BLOCKS=4`), and `k2k_rho_walk_lowmem<W>` at 400, 464 and 592 for
W = 4, 8, 16. At 512 resident threads/SM, W = 16 on the batched variant is
930 KB of local footprint per SM. Since the inversion being amortised is
only ~23 multiply-equivalents to begin with, going from W = 8 to 16 buys
1.45 multiplies against nearly doubled spill traffic. Nothing here has run
on a GPU, so the crossover is unmeasured — pick W on hardware, not from this
table, and weigh `lowmem` seriously: it costs 1.46 more multiply-equivalents
per step at W = 8 and less than half the local memory.

**Nothing in the step is computed twice.** A step needs τ^j(P) and g(x_P),
and the obvious arrangement pays for each twice: both phases apply τ^j (and
to both coordinates, though phase A only reads x), and the class weight is
evaluated for the distinguished-point test and again for the next step's j.
Instead phase A walks only the x chain and hands phase B the denominator
x₁ + x₂ it already formed, phase B walks only the y chain, and the weight
is carried in the walk state. At W = 8 that is 12.91 → 10.62
multiply-equivalents and one class-weight evaluation per step instead of
two — and that table is the hottest in the kernel. The `lowmem` variant
takes the weight saving but still re-walks the x chain, 15.52 → 12.08.

Carrying the denominator is a trade, not a free win: it costs one extra
field element per walk, and that is why `k2k_rho_walk<8>` went from 872 to
1040 bytes of local memory while `lowmem<8>` went from 448 to 464.
`ecc2k130` makes both choices under a flag — `PACKED_CACHE_DENOM`, described
in its `DENOMINATOR-CACHE.md`, against the recompute path its `walk.h`
prices at "786 instructions for 1048 bytes of traffic per slot". Which side
wins is a memory-hierarchy question, so the two variants here keep opposite
answers until someone measures on a GPU.

## Occupancy: measured

`./ptx_stats2k.sh`, `k2k_rho_walk_lowmem<8>`, block size 128:

| Arch | `R2K_MIN_BLOCKS` | Registers | Stack | Resident threads/SM |
|---|---|---|---|---|
| sm_90 | unset | 142 | 456 B | 384 |
| sm_100 | unset | 142 | 456 B | 384 |
| sm_120 | unset | 148 | 456 B | 384 |
| sm_90 / sm_100 / sm_120 | 4 | 128 | 464 B | 512 |

(clang 18.1.3 for sm_90, ptxas 12.9.86.)

Two things stand out against the prime-field walk in `gpu/ecc`, which needs
226 registers and 792 bytes of stack for 256 resident threads:

- **The binary-field walk is much lighter.** A field element is four words
  instead of eight, the addition is a handful of XORs, and no Montgomery
  form has to stay live. It reaches 384 threads per SM untuned and 512 with
  the launch bound, so 1.5–2× the occupancy of the prime-field kernel.
- **All three architectures agree.** The prime-field kernel hits the
  255-register ceiling on `sm_120` while using 226 on `sm_100`; nothing like
  that happens here, because the kernel never gets near the ceiling. Tuning
  transfers between consumer and datacenter Blackwell for these kernels in
  a way it does not for the prime-field ones.

## The `|F|²` pair table

`pairtable.cuh` builds the meet-in-the-middle table that
`RESEARCH_KOBLITZ_INDEX_CALCULUS.md` measures as the whole wall-clock
cost of a relation-collection worker:

> two workers run concurrently on the same four cores collected units 0
> and 1 in 23 s each (each builds its own `|F|²` pair table, **which is
> the whole cost**; the 64 probes are milliseconds)

The table is `|F|(|F|+1)/2` curve additions and nothing else.

**What makes it more than batch point addition.** An affine addition on
a binary curve costs one inversion, and an inversion is ~`m` squarings
and `m` multiplications — two orders of magnitude more than the three
multiplications the rest of the addition needs. So the table is not
built with `|F|²/2` inversions: fixing `P_i` and forming the whole row's
denominators `x(P_i) + x(P_j)` lets Montgomery's trick invert them with
**one** inversion and `3(k−1)` multiplications, exactly as
`PairSumTable::build_within` does on the CPU.

Montgomery's trick is a sequential prefix product, so it does not
parallelise across a row — it parallelises across *rows*, and there are
`|F|` of them. One thread owns one row. That leaves a triangular
imbalance (row 0 does `|F|` additions, the last does one), so the kernel
grid-strides and the host can pick a grid smaller than `|F|`.

**Packing is bit-for-bit `koblitz_fast::FastPoint::pack`**, so a table
built here and one built there sort and look up identically:
`infinity → 0`, `(x, y) → ((x+1) << 1) | (y > (x ^ y))`. The low bit is
the part that needed a test: on a binary curve `−P = (x, x+y)`, so `P`
and `−P` share an abscissa and a single bit of `y` separates them only
when `x` is odd. Comparing `y` against `x+y` separates them always.

`make test` runs `test_pt_*` for every curve: batch inversion against
one-at-a-time inversion (with a planted zero, which every real row has
at `j = i`), every row of a 48-point triangle against unbatched
`Koblitz::add`, every sum re-checked as on-curve, and the packing
checked for collisions and for keeping `P` from `−P`. The row test
deliberately includes negated points so that `P_i + P_j = O` occurs —
without them the infinity branch is never taken and the test would claim
coverage it does not have. `ecc2k95` skips it: at `m = 97` a point does
not fit in a `u64` and `FastCurve` refuses the same case.

## Files

| File | What it is |
|---|---|
| `ecref2k.py` | Python oracle: F(2^m), Koblitz curves, group order, the τ eigenvalue s, the normal-basis map, header and vector generation |
| `f2m.cuh` | Binary field: carry-less multiply, squaring, Itoh–Tsujii inversion, batch inversion, the class weight |
| `koblitz.cuh` | Points, Frobenius, scalar multiplication |
| `rho2k.cuh` | The Frobenius-class walk, distinguished points, canonicalisation, batched stepping |
| `scalar_ring.hpp` | Host-only Montgomery arithmetic mod r |
| `rho2k_host.hpp` | Walk replay, class relation, the solve |
| `kernels2k.cuh` | Kernels and launch structure |
| `bench2k.cu` | Device self-test and benchmarks |
| `test_cpu2k.cpp` | CPU verification harness |
| `ptx_stats2k.sh` | Static instruction and occupancy analysis, no GPU needed |

## Shared-memory detail worth keeping

The class-weight table is read `F2M_CB_WINDOWS` times per step by every
walk, at data-dependent indices — by far the hottest table in the kernel.
Its entries are padded to `F2M_WORDS + 1 = 5` words by the generator. An odd
stride is invertible mod 32, so the sixteen nibble values a warp can present
land in sixteen different banks; a power-of-two stride would serialise the
access four ways. This is the same reasoning that keeps the `inf` word in
the affine point struct over in `gpu/ecc`, and it is the kind of thing that
gets "cleaned up" later, so it is commented in place.

## Scope

This implements the attack's inner loop. It does not include the Certicom
challenge's actual base point and target, which would simply be plugged into
`ecref2k.py`; the generator here is derived deterministically. Nor does it
make ECC2K-95 newly breakable — it was broken twenty-eight years ago, and
`./bench2k rho` prints the device-years the full search would take at the
measured rate so the scale is explicit.
