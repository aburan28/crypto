# Sparse-basis compound route: frozen static screen

Source parent: `5434e208953527d2f37e14bde9f6c470542922cf`.

Reference: the independently replayed counter-free fused sigma path at
**15.436677 billion complete scalar updates/s** on one RTX PRO 6000.  This is
a hash-pinned result from commit
`e0b0858b179fe28f753353e7a7313cfe23d71884`: the headline result SHA-256 is
`d31f2d758d2e787abe8f043b1b2ae69d337ebad9a30403afe74d50d4a003dd8c`
and its independent-audit SHA-256 is
`1b1a46ca52aa14ef86b96fc9b765c9c10e2815c1fe9bcad867422785ab7976d1`.
This benchmark-only arithmetic engineering screen does not launch a GPU,
walk a challenge interval, collect distinguished points, solve a collision,
or recover a scalar.

## Hypothesis

Keep point coordinates, the Montgomery prefix/reverse products, and the
per-update lambda arithmetic in a sparse degree-131 polynomial basis.  The two
candidate moduli are fixed before synthesis:

```
f_1 = z^131 + z^8 + z^3 + z^2 + 1
f_2 = z^131 + z^8 + z^5 + z^2 + 1
```

For branch selection, convert sparse-polynomial X directly to the repository's
normal coordinates with a three-input-bit shared-memory lookup table.  For the
selected `j=3..10`, apply a direct sparse-polynomial linear map
`L_j(a) = a + a^(2^j)` independently to X and Y through the same table shape.
Convert only the completed batch product to the existing beta polynomial
basis, run the existing default eight-product Itoh--Tsujii inverse
(`PACKED_INV_POLY=0`), and convert the inverse back to the sparse basis.  That
inverse converts beta input to normal coordinates, keeps its accumulator in
the normal basis, and converts the unreduced products directly back through
`fromPolynomialProduct131`; its eight products do not call the direct
polynomial reducer.  At batch 16 the two new sparse/beta maps are charged at
two applications per batch.  The inverse's existing internal transforms stay
unchanged and common to the reference.

## Native construction and exact oracle

One C++17 program performs the complete screen.  Python and Python-generated
research artifacts are excluded.  For each candidate modulus the program:

1. proves the prime-degree Rabin irreducibility conditions;
2. constructs a deterministic order-263 root in the sparse base field,
   derives a root of the repository beta-basis modulus, and builds both
   directions of the exact field isomorphism;
3. checks both 131-by-131 conversion matrices have rank 131 and round-trip all
   131 basis vectors;
4. checks sparse multiplication/reduction against beta-basis multiplication
   on all 261 raw-product basis vectors and deterministic dense vectors;
5. derives sparse-to-normal and all eight sparse `L_j` matrices, proves the
   table evaluator on every input basis vector, and rechecks deterministic
   dense vectors against the repository's actual normal-basis conversion and
   Frobenius implementation; and
6. emits deterministic counts, table sizes, SHA-256-bound source/output, and a
   decision.  A failed assertion produces no admissible result.

All 131 Frobenius-conjugate embeddings are mathematically equivalent.  The
program reports the first deterministic embedding produced by the fixed
order-263 construction.  Table shape and evaluator counts do not depend on
matrix density, so no post-hoc embedding search can improve the admitted
three-bit route.

## Fixed circuit and accounting

Each linear map is split into 44 consecutive three-bit chunks.  Every table
entry is five 32-bit output words.  The low 128 bits are laid out as one
aligned vector and the top three bits in a separate word, matching the prior
direct-sigma screen.  Per map application the ledger charges:

- 7,040 table bytes (`44 * 8 * 5 * 4`);
- 44 aligned 128-bit shared loads plus 44 top-word shared loads, or 88 LDS
  instructions (220 scalar-word equivalents);
- 315 exact straight-line data ALU operations: 44 right shifts, 44 masks,
  three boundary left shifts and three ORs, 220 accumulator XORs, and one
  final top-word mask; the admission ledger retains the prior direct-sigma
  screen's conservative charge of 317; and
- five live output words, one extracted index, and the five input words.

Nine hot maps (normal conversion plus eight `L_j`) and two batch-boundary
basis maps are retained together.  Per complete scalar update at batch 16 the
candidate is therefore charged for `3 + 2/16 = 3.125` table applications.
The full table allocation, LDS/ALU counts, product/reduction counts, CLMAD
count, and minimum live storage are reported explicitly.

The arithmetic ledger inherits the fused reference's batch 16 defaults
`PACKED_INV_POLY=0`, `PACKED_ALU_SQUARE=0`, and `PACKED_ALU_SQR=0`: six
`CLMAD`s per field product, five per polynomial lambda square, and four per
normal-basis square.  These settings are stated because changing any of them
would change the pipe count and define a different compound route.

The sparse reducer is the fixed two-fold circuit obtained from
`z^131 = z^8 + z^b + z^2 + 1`, with `b=3` or `5`.  It is compared with a
literal native transcription of the shipping direct reducer.  Exact 32-bit
shift, OR, AND, and XOR counts are part of the result; compiler fusion or host
assembly is not substituted for that circuit count.

## Conservative RTX PRO 6000 gate

The static roofline uses the repository's measured RTX PRO 6000 rates and the
most favourable recorded clock: 188 SMs, 2.430 GHz, 64 ALU lanes/SM-clock,
9.2 random shared-memory lane-instructions/SM-clock, and 2.0 paired-CLMAD
lanes/SM-clock.  The rate source is
`benchmarks/fast-clmad/probe/summary.json`, SHA-256
`57d8541c28648911e3e59a9a6bdd34f114658f5df0b299bab651193293b4dc50`;
the SM register/shared limits come from
`benchmarks/hardware-limits/result.json`, SHA-256
`31d1c91ca70006086d0454f51be5a5f5866afdbab6c06adcda88026c88a54c1a`.
The lookup-only projection deliberately makes every carryless
product, reduction, conversion arithmetic operation, state access, and loop
instruction free.  The CLMAD row likewise makes every other instruction
free.  These are conditional one-pipe models, not predicted timings: the 9.2
rate was measured for random `LDS.U8`, not this exact aligned
`LDS.128`+`LDS.32` multicast layout, so the separate 16-lane sensitivity must
also be reported.

A bounded whole-walk GPU prototype is justified only if all exact checks pass,
the tables fit one block within the 101,376-byte opt-in limit, no static
storage estimate requires more than 255 registers/thread, and every frozen
one-pipe model exceeds the 15.436677 B/s reference by at least 5%.  Failure is a
static no-go for this compound route.  It does not rule out a different linear
circuit or a new measured shared-memory primitive.

## Pre-result correction

The first implementation attempt exposed an error in item 2 before it emitted
a result: `2^131 mod 263` is `1`, not `-1`.  The 263rd roots therefore lie in
the base field and are obtained with exponent `(2^131 - 1) / 263`; a quadratic
extension is neither required nor correct for this construction.  This
correction was committed before the admitted native run.  It changes no
candidate modulus, circuit family, counts, throughput reference, or decision
gate above.

The post-run source review found three, rather than four, chunk extractions
crossing a 32-bit word boundary.  The exact straight-line evaluator count is
therefore 315.  The frozen admission calculation deliberately keeps the
earlier, conservative 317-operation charge; the lookup-only wall decides the
result, so this accounting correction cannot change the gate.
