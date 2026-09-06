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
j(P) = (g(x_P) mod 8) + 3
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
[solve] 21-bit Koblitz DLP, 128 walks, 24 trials
  mean 267 steps over 24 solves
  sqrt(pi/2 * r/(2m)) = 268  -> measured 1.00x   [class walk]
  sqrt(pi/2 * r/2)    = 1283 -> measured 0.21x   [negation only]

[solve] 40-bit Koblitz DLP, 128 walks, 8 trials, dp threshold g <= 11 (rate 1/463)
  mean 209744 steps over 8 solves
  sqrt(pi/2 * r/(2m)) = 145129 -> measured 1.45x   [class walk]
  sqrt(pi/2 * r/2)    = 929276 -> measured 0.23x   [negation only]
```

The 21-bit run lands on the class-walk prediction exactly. The 40-bit run
sits 45% above it, and that excess is accounted for: 128 walks each need
about 463 steps to reach a distinguished point, so 59k steps of the 210k
are the DP tail, and 145k + 59k ≈ 204k.

## What the field arithmetic actually costs

Measured statically (`./ptx_stats2k.sh`, clang for sm_90):

| Operation | PTX instructions | Notes |
|---|---|---|
| F(2^97) multiply | 405 | 9 carry-less 32×32 products, Karatsuba twice |
| F(2^97) squaring | 79 | 0.20 of a multiply — bit spreading, no multiplier |
| secp256k1 multiply | 468 | for comparison, from `gpu/ecc` |

**A 97-bit binary field multiply costs about as much as a 256-bit prime
field multiply.** That is the single most important fact about binary-field
ECC on a GPU, and it has one cause: there is no carry-less multiply
instruction. No PCLMULQDQ, no VMULL, nothing. The 32×32 carry-less product
is synthesised from sixteen ordinary widening multiplies with the
interleaved-mask trick — split each operand into four subsets by bit index
mod 4, multiply as integers, mask. Within a subset the per-bit partial sums
cannot exceed 8, which fits in the 4-bit gap between kept bits, so no carry
ever crosses into a bit that matters.

What Koblitz curves give back:

- **Squaring is a fifth of a multiply**, so τ is nearly free. This is the
  whole basis of the 13.9× class speedup.
- **Inversion is cheap.** Itoh–Tsujii computes a^(2^m−2) with about 7
  multiplications and 96 squarings, roughly 23 multiply-equivalents,
  against 270 for a Fermat inversion modulo a prime.

That last point changes the kernel's shape. Batching inversions across
walks is worth 27× over a prime field; here it is worth 3–4×, and W beyond
16 buys almost nothing:

| Batch W | Multiplies per rho step | |
|---|---|---|
| 1 | ~28 | one inversion per step |
| 8 | ~11 | |
| 16 | ~9.4 | |
| 32 | ~8.7 | diminishing |

(Per step: 2 for the addition, 1 squaring, 3 for Montgomery's trick, the
Frobenius applications, the class weight, and 23/W for the inversion.)

## Occupancy: measured

`./ptx_stats2k.sh`, `k2k_rho_walk_lowmem<8>`, block size 128:

| Arch | `R2K_MIN_BLOCKS` | Registers | Stack | Resident threads/SM |
|---|---|---|---|---|
| sm_90 | unset | 142 | 440 B | 384 |
| sm_100 | unset | 150 | 440 B | 384 |
| sm_120 | unset | 148 | 440 B | 384 |
| sm_90 / sm_100 / sm_120 | 4 | 128 | 448 B | 512 |

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
