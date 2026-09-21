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
3. **Every faster routine against the one it replaced.** `f2m_prod` against
   `clmul128` and `f2m_reduce` against `f2m_reduce_generic`, on 60,002
   products per curve: real multiplies, real squarings, random buffers of
   the widest degree a product can have, and both corners. Then each τ^k
   table against repeated squaring, the inversion that uses it against the
   inversion that does not, and a 300-step run of the real stepper with the
   tables against the same run without, which has to come out bit-identical
   down to the distinguished-point multiset.
4. **The identity the whole attack rests on**: τ(G) = s·G, checked both
   against Python and by computing s·G on the device path.
5. **Class equivariance**: f(τP) = τf(P) and f(−P) = −f(P) on sample points.
6. **Three steppers agree.** The batched, low-memory and unbatched
   implementations must produce bit-identical state over 300 iterations and
   the same multiset of distinguished points.
7. **End to end, and quantitatively.** The toy curves are solved through the
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
| F(2^97) multiply | 304 | 6 carry-less 32×32 products, three-way Karatsuba |
| F(2^97) squaring | 78 | 0.26 of a multiply — bit spreading, no multiplier |
| F(2^97) class weight | 165 | 25 windowed table reads, once per walk step |
| F(2^97) τ^k, windowed | 296 | 101 of them loads; worth it from k = 4 up |
| secp256k1 multiply | 472 | for comparison, from `gpu/ecc`; 210 with its inline-PTX carry chains |

The multiply was 405 and the squaring 79 before the widths below were cut
down; the whole before-and-after is in the next section.

**A 97-bit binary field multiply still costs most of what a 256-bit prime
field multiply costs** — because this backend multiplies in software. The
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
`ecc2k130/include/packed131.h`. This backend has not been ported, and it
remains the largest single thing left on the table: six software `clmul32`
leaves are still 96 of the multiply's 304 instructions in widening
multiplies alone, where `clmad` would do a 128-bit product in six
instructions and make the three-way split below unnecessary for the
multiply. Porting it needs a `ptxas` from CUDA 13.3 or newer, which is why
it is not done here — the toolchain `ptx_stats2k.sh` assembles from the pip
wheels tops out at 12.9, and nothing in this directory has ever run on a GPU.

What Koblitz curves give back:

- **Squaring is a quarter of a multiply**, so τ is nearly free. This is the
  whole basis of the 13.9× class speedup. It was a fifth before the multiply
  got cheaper; the ratio is a property of the pair, so it moves when either
  side does, and anything quoted in "multiply-equivalents" has to name the
  ratio it used.
- **Inversion is cheap.** Itoh–Tsujii computes a^(2^m−2) with about 7
  multiplications and 96 squarings, against 270 multiplications for a Fermat
  inversion modulo a prime.

That last point changes the kernel's shape. Batching inversions across
walks is worth 27× over a prime field; here it is worth 3–4×, and W beyond
16 buys little. The operation counts are counted, not derived: `make test`
builds with `F2M_COUNT_OPS` and `./test_ecc2k95` runs the real stepper with
the distinguished-point test disabled, so every walk pays exactly one step
and none reseeds. Multiplying them by the instruction counts above gives the
per-step cost:

| Batch W | mul | sqr | PTX instructions/step, before | after | ratio |
|---|---|---|---|---|---|
| 1 | 9.00 | 108.90 | 12425 | **11395** | 0.917 |
| 2 | 7.00 | 61.23 | 7844 | **7069** | 0.901 |
| 4 | 6.00 | 37.15 | 5534 | **4887** | 0.883 |
| 8 | 5.50 | 25.33 | 4397 | **3813** | 0.867 |
| 16 | 5.25 | 19.32 | 3820 | **3268** | 0.855 |
| 32 | 5.12 | 16.49 | 3543 | **3008** | 0.849 |

(Per step: 2 multiplies for the addition, 1 squaring for λ², 3 multiplies
for Montgomery's trick, 2j squarings for the two Frobenius chains, one class
weight, and one inversion — 96 squarings and 7 multiplies — split W ways.
The operation counts are the same in both columns: the arithmetic narrowing
below changes what a field operation costs, not what the walk computes.)

Where the W = 8 step goes: 44% multiplies, 27% the two Frobenius chains, 25%
the amortised inversion's squarings, 4% the class weight. Squarings are the
larger half, and every one of them is part of some τ^k — which is a fixed
linear map and need not cost k squarings. That is the "τ^k is not k
squarings" section below.

**This table counts arithmetic only, and arithmetic is not what picks W.**
`ptxas` puts `k2k_rho_walk<W>` at 344, 464, 656, 1048 and 1816 bytes of
per-thread local memory for W = 1, 2, 4, 8, 16 (sm_90, block 128,
`R2K_MIN_BLOCKS=4`), and `k2k_rho_walk_lowmem<W>` at 400, 464 and 592 for
W = 4, 8, 16. At 512 resident threads/SM, W = 16 on the batched variant is
930 KB of local footprint per SM. Since the inversion being amortised is
only 9,616 instructions (31.6 multiplies) to begin with, going from W = 8 to
16 buys 544 instructions a step, 14.3%, against nearly doubled spill traffic.
Nothing here has run on a GPU, so the crossover is unmeasured — pick W on
hardware, not from this table, and weigh `lowmem` seriously. It is counted
in the same run, and it costs this much more arithmetic for less than half
the local memory:

| Batch W | batched | lowmem | ratio |
|---|---|---|---|
| 1 | 11395 | 12467 | 1.094 |
| 2 | 7069 | 7850 | 1.110 |
| 4 | 4887 | 5513 | 1.128 |
| 8 | 3813 | 4369 | 1.146 |
| 16 | 3268 | 3789 | 1.159 |
| 32 | 3008 | 3516 | 1.169 |

**Nothing in the step is computed twice.** A step needs τ^j(P) and g(x_P),
and the obvious arrangement pays for each twice: both phases apply τ^j (and
to both coordinates, though phase A only reads x), and the class weight is
evaluated for the distinguished-point test and again for the next step's j.
Instead phase A walks only the x chain and hands phase B the denominator
x₁ + x₂ it already formed, phase B walks only the y chain, and the weight
is carried in the walk state. At W = 8 that was 12.91 → 10.62
multiply-equivalents and one class-weight evaluation per step instead of
two — and that table is the hottest in the kernel. The `lowmem` variant
takes the weight saving but still re-walks the x chain, 15.52 → 12.08.
(Those four figures are multiply-equivalents at the 1 squaring = 0.2
multiplies the arithmetic had when they were measured. The ratio is 0.26
now, so they are not comparable with the instruction counts above; the
operation counts they rest on are unchanged.)

Carrying the denominator is a trade, not a free win: it costs one extra
field element per walk, and that is why `k2k_rho_walk<8>` went from 872 to
1040 bytes of local memory while `lowmem<8>` went from 448 to 464.
`ecc2k130` makes both choices under a flag — `PACKED_CACHE_DENOM`, described
in its `DENOMINATOR-CACHE.md`, against the recompute path its `walk.h`
prices at "786 instructions for 1048 bytes of traffic per slot". Which side
wins is a memory-hierarchy question, so the two variants here keep opposite
answers until someone measures on a GPU.

## The container is not the field

`F2M_WORDS` is four words for every field here, and almost nothing needs all
four. An element is `ceil(m/32)` words, of which `m/32` are full and, when m
is not a multiple of 32, one holds `m mod 32` bits. A product of two of them
has degree at most 2m−2, and each fold of the reduction shortens what is
left again. Three places in the arithmetic could have been running at the
container's width rather than the field's. Two of them were:

- **The product.** Karatsuba over four words spends a full 32×32 leaf on a
  top word that at m = 97 holds **one bit**. Splitting it off —
  `A = A_lo + a_top·t^96`, so `A·B = A_lo·B_lo + (A_lo·b_top + a_top·B_lo)·t^96
  + a_top·b_top·t^192` — leaves a three-word Karatsuba, six leaves instead of
  nine, plus seven products by a single bit, which are masks. This is not a
  quirk of m = 97: the Koblitz challenge fields sit just above a word
  boundary, 131 = 4·32 + 3 as well.
- **The reduction.** Both folds ran the full eight words. At m = 97 the
  input is seven, `hi = T >> m` is three, the first fold reaches four, and
  the second fold's `hi2` is **five bits** — an extract and two XORs, not a
  second pass.
- **The spread.** Nothing to do: the 64-bit form already lets `clang` drop
  the top word's spread when it can see the operand is reduced. Rewriting it
  as two 32-bit spreads of the half-words looks like the natural shape for
  32-bit lanes and costs 104 instructions per chained squaring against 78,
  so it is deliberately still in 64-bit form.

Per operation, and per rho step at W = 8 (`./ptx_stats2k.sh`; the toy curves
are here because they show what the container was costing when it is most of
the element):

| | m = 97 multiply | squaring | step, W = 8 | m = 41 mul | m = 23 mul |
|---|---|---|---|---|---|
| before | 405.1 | 79.1 | 4397 | 468.0 | 461.0 |
| after | **304.0** | **78.0** | **3813** | **146.0** | **52.0** |
| ratio | 0.750 | 0.986 | 0.867 | 0.312 | 0.113 |

Widening multiplies per m = 97 multiply fall from 132 to 96, which is
exactly the six surviving leaves. Registers and occupancy do not move:
`k2k_rho_walk<8>` and `k2k_rho_walk_lowmem<8>` stay at 128 registers and 512
resident threads/SM on sm_90, sm_100 and sm_120, with at most 8 more bytes
of stack.

As a practicality note and not as the metric, the 40-bit end-to-end solve in
`make test` — eight real discrete logarithms on one core — drops from
2.1–2.7 s to 1.02–1.05 s over three runs each. It takes the same 191,840
steps on both sides, which is the point: same walk, cheaper arithmetic.

**This is engineering, not an advance.** It changes what a field operation
costs, not how many of them a walk does: the operation counts in the table
above are identical on both sides, the walk computes the same function, and
`S = operations/√n` — the unit `docs/index-calculus-scoreboard.html` is drawn
in — is exactly flat. There is no row to add there.

**How it is checked.** Both narrowings are the only way the arithmetic can
now be wrong, so `make test` checks each against the full-width version it
replaces, on every curve: `f2m_prod` against `clmul128` and `f2m_reduce`
against `f2m_reduce_generic`, which is kept in `f2m.cuh` for exactly this
purpose. 60,002 products per curve — from real multiplies, from real
squarings, from random buffers of the widest degree a product can have, and
at both corners — plus the Python field vectors and the end-to-end solves
that ran before.

## τ^k is not k squarings

Squaring is cheap, which is the whole reason to walk on Frobenius classes,
and it makes it easy to miss that **τ^k is a linear map and its cost need not
grow with k**. It is F2-linear, so it is fixed by where it sends each `t^j`,
and applying it is the same 4-bit-window table read the class weight already
does: 25 windows, one read and one word-XOR each, at any k.

Measured on sm_90: **296 instructions for a window application against 78 for
a squaring**, so a table pays from k = 4 up. That is the whole rule, and the
generator applies it — `ecref2k.py` reads the Frobenius exponents out of the
same walk of m−1 that `pow_2n_minus_1` does, keeps those ≥ 4, and emits a
table for each. At m = 97 the Itoh–Tsujii chain applies τ at k = 1, 3, 6, 12,
24, 48, so four tables come out: τ^6, τ^12, τ^24, τ^48. τ^48 alone is 3,744
instructions of squaring replaced by 296.

Counted by `./test_ecc2k95`, the W = 8 step goes from 25.33 squarings to
14.08 squarings and 0.50 table applications, with the multiply count
untouched:

| Batch W | batched | batched + τ^k | ratio | lowmem |
|---|---|---|---|---|
| 1 | 11395 | **5559** | 0.488 | 12467 |
| 2 | 7069 | **4151** | 0.587 | 7850 |
| 4 | 4887 | **3428** | 0.701 | 5513 |
| 8 | 3813 | **3083** | 0.809 | 4369 |
| 16 | 3268 | **2904** | 0.889 | 3789 |
| 32 | 3008 | **2824** | 0.939 | 3516 |

The tables help most where the inversion is least amortised, so they pull the
batching optimum down — and the interesting consequence is in the local
memory, not the arithmetic. **W = 4 with tables costs 3,428 instructions a
step against W = 8 without at 3,813, and `k2k_rho_walk<4>` is 656 bytes of
local memory against `<8>`'s 1,048.** Fewer instructions and a third less
spill traffic, from halving the batch.

**This is a trade, and only one side of it is measured.** The tables replace
arithmetic with table traffic: at m = 97 an inversion reads 4 × 25 × 4 = 400
words, which at W = 8 is 50 words per walk step against the class weight's
100 — the same kind of read at half the rate, on 32 KB more table. They stay
in global memory rather than being staged into shared like the class-weight
table, precisely because Montgomery's trick has already made them a W-times
colder read; giving up shared memory for them would be paying occupancy for
the coldest table in the kernel.

So they are **off by default**: `rho2k_ctx.ftb` is null, which is the
squaring chain this had before, and `./bench2k rho --frob 1` turns them on.
Whoever runs this on a GPU first should measure both; the instruction count
says −19% at W = 8 and cannot say what the reads cost.

`./bench2k selftest` takes its device inversion through the tables while the
host reference squares, so the walk comparison that was already there checks
the table path against the thing it replaces, on hardware, for free.

**The walk's own τ^j chain is the obvious next target and does not survive
the same arithmetic.** It is the larger half — about 12.5 squarings a step
against the inversion's 11.25 at W = 8 — so tabling it looks like the bigger
win. Price it the way the trade actually runs, in ALU instructions saved per
table read added:

| | ALU replaced | reads added | ALU per read | tables |
|---|---|---|---|---|
| inversion τ^k | 780 | 50 | **15.4** | 4 (k fixed) |
| walk τ^j chain | 584 | 202 | **2.9** | 8 (j data-dependent) |

The inversion wins because Montgomery's trick has already divided its reads
by W while the squarings it replaces are not divided by anything. The walk
chain gets no such discount: it runs twice per step per walk, so the same
substitution buys three ALU instructions per read instead of fifteen, on
eight tables rather than four because j is data-dependent. A table read has
to be cheaper than three ALU slots for that to pay, which on a GPU it
generally is not. Not implemented, and the number above is why.

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
`research/notes/ecc2k130/RESEARCH_KOBLITZ_INDEX_CALCULUS.md` measures as the whole wall-clock
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

**The folded key, `pt_canon`.** The CPU side stores one entry per
*signed Frobenius orbit* rather than one per pair, which at `n = 61` is
`2n` times fewer entries and buys a base `√(2n)` times wider at the same
memory. The key that names an orbit is `pt_canon`, and it is the same
function `koblitz_fast::FrobeniusCanon::canon` is: `π` is a squaring
only in a *polynomial* basis, and in a normal basis `{β, β², β⁴, …}` it
is a one-bit cyclic rotation of the coordinate word, so the orbit of `x`
is the set of rotations of that word and the least rotation names it.
Nothing is squared and nothing is reduced.

The basis change is **host data**, uploaded from
`FrobeniusCanon::tables()`, not rediscovered on the device. The normal
element comes from a randomised search, so a device that searched for
its own would find a different basis and name the same orbits
differently — valid on its own, and unable to read a table the CPU
built. `n ≤ 62` is the ceiling, so the fold does not reach `ecc2k95`.

**What is keyed is not yet what is stored.** `pairtable_kernel` takes
the tables and keys by orbit when given them, but it still emits one
entry per pair: the table is the same size and merely named
differently. The fold's actual saving needs the base sorted by orbit
host-side and one entry per signed orbit, and then the orbit tag that
makes summand recovery `O(n)` instead of `O(|F|)`. Both are the next
steps.

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
| `ecref2k.py` | Python oracle: F(2^m), Koblitz curves, group order, the τ eigenvalue s, the normal-basis map, the windowed linear-map tables, header and vector generation |
| `f2m.cuh` | Binary field: carry-less multiply, squaring, Itoh–Tsujii inversion, batch inversion, the class weight, τ^k as a linear map |
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
