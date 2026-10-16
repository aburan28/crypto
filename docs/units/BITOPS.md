# `bitop/v1`: one unit of work for every ECDLP method

Every ECDLP cost in this repository, in `aburan28/cryptanalysis` and in cairn's
bound objectives is charged in one unit: the **bit operation**, one two-input
Boolean gate (AND, XOR, OR) in the textbook circuit for the primitive. A
method's cost is

```text
work = Σ count(primitive) · gates(primitive, width)
```

an exact integer that depends on the counts and this document, and on nothing
else. No host, compiler or run enters it. The implementation is
`src/cryptanalysis/bitops.rs`. `docs/units/bitop-v1.vectors.json` holds the
gate count of every primitive at a spread of widths. Any other implementation
must reproduce every row exactly.

## Why one unit

Before `bitop/v1` the repositories charged work in at least six units:

- group-addition equivalents, with non-group work converted by nanosecond
  ratios from one host;
- solver work priced by this host's wall time, or not at all;
- `full_prices.json` ratios;
- reference picoseconds (`cryptanalysis` ic-bench);
- Callgrind instructions;
- wall-clock online speedups.

You cannot divide a figure in one of these by a figure in another. An
exponent fitted in one is not an exponent in another.

The worst consequence was silent: work with no conversion ratio was charged
zero. Examples:

- linear algebra on any curve without a row in `docs/ic/calibration.json`;
- Buchberger, SAT and PDP solver work.

Every index-calculus total was therefore a lower bound, and the ratio of two
lower bounds bounds nothing.

## The conventions

Throughout, `b` is the bit length of `p`, `n` the degree of `F_{2^n}`, and
`w` the bit length of the subgroup order `r`.

### Integers

| primitive | gates | derivation |
|---|---|---|
| `w`-bit addition | `5w` | ripple carry, a full adder is 2 XOR + 2 AND + 1 OR |
| `w × w` product | `w² + 5w(w−1)` | `w²` partial-product ANDs, `w − 1` row additions |

### Prime field `F_p` (Montgomery representation)

| primitive | gates |
|---|---|
| `field_add` | `2·5b + 3b = 13b`: an addition, a conditional subtraction of `p`, and a `b`-bit select |
| `field_mul` | `2·(b² + 5b(b−1)) + 13b`: the product, the reduction's product, and the final subtraction |
| `field_sqr` | as `field_mul` |
| `field_inv` | `2b · 13b`: `2b` steps of the binary extended Euclid, each a subtraction, a comparison and a select |
| `field_sqrt` | `2b · field_mul`: exponentiation (`b` squarings, `b/2` products) plus a `b/2`-product allowance for the Tonelli–Shanks loop |
| `legendre` | as `field_inv` (the Jacobi algorithm is a gcd) |

### Binary field `F_{2^n}` (polynomial basis)

Reduction is by a pentanomial, i.e. four low terms. A trinomial is cheaper;
the convention takes the dearer modulus so that no curve is flattered by its
own.

| primitive | gates |
|---|---|
| `field_add` | `n` XORs |
| `field_mul` | `n² + (n−1)² + 4(n−1)`: schoolbook carry-less product plus reduction |
| `field_sqr` | `4(n−1)`: spreading the bits is wiring, then the reduction |
| `field_inv` | `2n · 5n`: polynomial extended Euclid, `2n` steps of a shift-XOR pair and a select |
| `as_solve` | `n² + n(n−1)`: half-trace as a precomputed `n × n` matrix over `F_2` times a vector |
| `canon` | `as_solve + 2n²`: normal-basis conversion, then the least of `n` rotations |

### `F_{p⁴}` (`ecbench::fp4`)

| primitive | gates |
|---|---|
| `field_mul`, `field_sqr` | `16 · field_mul(F_p)` (schoolbook) |
| `group_add`, `group_double` | `97 · field_mul(F_p)` (`ecbench::fp4`'s `C_ADD`, the extension inversion through the norm) |

### Group operations (affine, as every counted group in `ecbench` computes them)

| curve | addition | doubling |
|---|---|---|
| `y² = x³ + ax + b` over `F_p` | `I + 2M + S + 6A` | `I + 2M + 2S + 8A` |
| `y² + xy = x³ + ax² + b` over `F_{2^n}` | `I + 2M + S + 9A` | `I + 2M + 2S + 5A` |

A scalar multiplication counted as such (rather than as the additions and
doublings it performed) is charged `1.5·w` additions. That is the convention
of the strong rho reference, which is the only method that counts one this
way.

### Arithmetic modulo `r` and bit-packed work

| primitive | gates |
|---|---|
| `modr_mul` | `field_mul` of `F_p` at `w` bits |
| `modr_add` | `field_add` of `F_p` at `w` bits |
| `modr_inv` | `field_inv` of `F_p` at `w` bits (binary extended Euclid) |
| a multiply-subtract of elimination (`row_op`) | `modr_mul + modr_add` |
| `word_op` | 64: one 64-bit XOR, AND or compare of a bit-packed algorithm (F4 and Macaulay rows, Gray-code enumeration) |
| `hash_probe` | `4κ` with `κ` the bits of a field element: a fold-XOR hash, a key comparison and a select. Keys are field elements, already uniform, so no multiplier is needed. A table's real cost is memory, which is reported beside the unit and never folded in |

## What the unit leaves out, and how that is caught

The unit is a gate count. It deliberately leaves out several things, which are
reported beside it and never folded in:

- **Memory.** Table sizes and peak resident memory.
- **Control.** Branches, loop overhead and allocation.
- **Uncounted work.** Anything a method does that no counter sees.

The last of these is the dangerous one, and the **instruction audit** catches
it. `ecbench audit-ir` replays a record under Callgrind and divides each
phase's instructions by its counted bit operations. A phase whose ratio is far
above the reference walk's is doing work its counters miss: an F4 elimination
that counts its XORs but not its symbolic preprocessing, for example. A
record that fails the audit is not admitted as a bound.

Work that was uncounted when sessions were committed is counted in a phase's
`work_native` map, never in `native`, so committed records keep replaying
exactly. The code that counts an item of a family's uncounted list also writes
the marker `counted:<item>`; a record that carries the marker leaves the item
off its `known_uncounted` list, and a record without it (every record written
before the counter existed) keeps it.

A record whose counters include one with no `bitop/v1` rule has no total. Its
`work.complete` is `false` and the unpriced counters are named. An unknown
cost is never written as zero.

## Changing the convention

A committed figure must reprice to the same integer forever, so this version
is frozen. A different convention gets a new name:

- `bitop/v2` for Karatsuba products;
- a trinomial-aware reduction;
- a memory-cost model.

Every record keeps its primitive counts, not only its total, so any later
version re-prices the whole history without re-running it.
