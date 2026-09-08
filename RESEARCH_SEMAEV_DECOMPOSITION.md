# Fast factor-base decomposition for binary Semaev `S₄`

`src/cryptanalysis/semaev_decomp.rs`, `examples/semaev_decomp_bench.rs`

This is the sequel to [`RESEARCH_SAT_SEMAEV.md`](RESEARCH_SAT_SEMAEV.md),
which ended on a negative result: a CDCL solver on the Weil-descended
`S₄` system spends **one conflict per candidate triple** and never
prunes, so it is a strictly worse enumerator than evaluating the
polynomial directly, and the gap widens with the factor-base dimension
`l` — 5× at `l = 5`, 62× at `l = 8`.

That conclusion was about *solvers*.  The measurement underneath it
pointed somewhere else: the cost of an index-calculus run is the
**decomposition oracle** — given a target `x_R`, are there three
factor-base elements `X₁, X₂, X₃` with `S₄(X₁,X₂,X₃,x_R) = 0`? — and
every implementation on the table was answering it by walking all
`2^{3l}/3!` triples at ~8 µs each.  Both halves of that are attackable.

The result: **926× faster** at `l = 8` against the enumeration baseline
the SAT work was measured on, and `l = 12` — 11.5 billion triples,
previously a day and a half of compute — decided in 15 seconds.

And it changes nothing about whether the attack works.  [Why](#what-this-does-not-buy)
is the more interesting half.

## Don't enumerate triples

`S₄` is a **quartic in its last argument**.  Fix `X₁` and `X₂` and it
becomes a degree-4 polynomial in `X₃`; its four roots are the only
candidates worth considering.  So loop over **pairs** and solve:

```text
  for each pair (X₁ ≤ X₂):              2^{2l}/2 of them
      build the quartic q(t) in X₃
      find the roots of q that lie in the factor base
```

The whole thing turns on that last line being cheap.  It is, because
the factor base is an `F₂`-subspace `V = ⟨1, z, …, z^{l−1}⟩`, and the
polynomial vanishing exactly on a subspace,

```text
  L_V(t) = Π_{v ∈ V} (t + v) = Σ_{i=0}^{l} aᵢ · t^{2^i},
```

is **linearized**: degree `2^l`, but only `l + 1` non-zero
coefficients.  So `L_V mod q` costs `l` squarings of a degree-3
polynomial, and

```text
  gcd(q, L_V mod q)
```

has as its roots exactly the roots of `q` that lie in `V`.  No search
over the factor base at all.

`L_V` is built one basis vector at a time from
`L_{i+1}(t) = L_i(t)² + L_i(b)·L_i(t)`, which holds because `L_i` is
additive, so `Π_{v ∈ V_i ∪ (b+V_i)}(t+v) = L_i(t)·L_i(t+b)`.

Per pair that is `O(l)` field operations where evaluation over the
factor base would be `2^l`, so the algorithmic saving grows like
`2^l / l`.  Measured, it does exactly that — see the `speedup` column
below, which doubles per rung.

Note this is not root-finding in `F_{2ⁿ}`.  The textbook route,
`gcd(q, t^{2ⁿ} + t)`, needs `n` squarings; restricting to the subspace
up front needs `l ≈ n/3`.

## Making the field fast

The first working version ran at **28 µs per pair**, for a few hundred
field operations' worth of work.  Three things accounted for that, and
only the last is the one people reach for first.

**Inversions (≈10×).**  A textbook polynomial remainder divides by the
divisor's leading coefficient, and at `n = 24` one field inversion —
Fermat, `n` squarings and `n` multiplications — costs about fifty
multiplications.  A single pair ran *thirteen* of them, all inside
`rem` and `gcd`.  Two changes remove nearly all:

- the modulus is made monic **once per row** of the pair loop, by
  Montgomery batch inversion (`Gf2::batch_inv`): one inversion and
  `3(k−1)` multiplications for `k` quartics, so the per-pair share goes
  to zero;
- the gcd uses **pseudo-remainders**, which scale the dividend up by
  the divisor's leading coefficient rather than scaling the divisor
  down.  That changes the gcd by a non-zero scalar, which changes
  neither its degree nor its roots — the only two things the caller
  reads.

**The multiply itself (≈6×).**  The schoolbook shift-and-xor loop runs
`n` iterations with two unpredictable branches in each.  `Gf2` uses a
carry-less multiply (`pclmulqdq`, with a scalar carry-less loop where
the instruction is unavailable) and folds the 128-bit product down
through a byte-indexed table of `v · z^{n+8j} mod irr` — `⌈(n−1)/8⌉`
lookups instead of `n` conditional shifts.  Squaring skips the
multiply entirely: in characteristic 2 it is the `F₂`-linear
bit-spreading `Σ aᵢ zⁱ ↦ Σ aᵢ z^{2i}`, then the same reduction.

**Branches in straight-line code (≈1.4×).**  Two that looked free and
were not: a zero-operand early-out in `mul`, and a "while the high part
is non-zero" exit in the reduction loop.  Both save real work and both
cost more in mispredictions than the work they save; the reduction loop
now runs a fixed number of times and eats one or two redundant lookups.

Separately, the crate's general `F2mElement` is `Vec<u64>`-backed and
allocates twice per multiplication.  For a field that fits in one word
that costs about a hundred times the arithmetic, which is why this
module carries its own field rather than reusing it.

## Measured

`cargo run --release --example semaev_decomp_bench 12 8`

Moduli are found by search and checked irreducible with Rabin's test at
startup, not read from a table.  `n = 3l` throughout, which keeps the
decomposition probability near `1/3!` as the factor base grows.

Targets are chosen **not to decompose**, so every method must exhaust
its space: no trajectory luck, and rejection is ~85% of an attack's
work anyway.

| `l` | `n` | triples | triples, general field | triples, one word | pairs + solve | algorithmic speedup |
|---:|---:|---:|---:|---:|---:|---:|
| 5 | 15 | 5 984 | 0.060 s | 0.001 s | 0.001 s | 2.2× |
| 6 | 18 | 45 760 | 0.364 s | 0.007 s | 0.002 s | 3.4× |
| 7 | 21 | 357 760 | 3.032 s | 0.046 s | 0.007 s | 6.7× |
| 8 | 24 | 2 829 056 | 27.767 s | 0.385 s | **0.030 s** | 12.9× |
| 9 | 27 | 22 500 864 | — | 3.753 s | 0.171 s | 21.9× |
| 10 | 30 | 179 481 600 | — | — | 0.692 s | — |
| 11 | 33 | 1 433 753 600 | — | — | 3.187 s | — |
| 12 | 36 | 11 461 636 096 | — | — | **15.312 s** | — |

The `speedup` column is pairs-and-solve against triple enumeration
*using the same arithmetic*, so it isolates the algorithm from the
field work.  It doubles per rung, which is what `2^l / l` predicts.

At `l = 8`, where all three methods still run:

| | time | vs. previous | cumulative |
|---|---:|---:|---:|
| SAT solver (`RESEARCH_SAT_SEMAEV.md`) | 1 472 s | — | — |
| triples, general field | 27.8 s | 53× | 53× |
| triples, one-word field | 0.385 s | 72× | 3 820× |
| pairs + solve | **0.030 s** | 12.9× | **49 000×** |

Two thirds of that is arithmetic and one part in thirteen is the
algorithm — but the algorithmic part is the one that keeps growing.
Extrapolating the per-triple rates to `l = 12` (0.167 µs/triple
one-word, measured at `l = 9`; 9.8 µs/triple general, measured at
`l = 8`) against the measured 15.3 s: **125× algorithmic, ≈7 300×
overall**, or 31 hours down to 15 seconds.

## What this does not buy

Worth being explicit, because the numbers above are large enough to
mislead.  Take the whole relation-collection phase:

- it needs `2^l` relations, one per factor-base element;
- a random target decomposes with probability `≈ 2^{3l−n}/3!`, so it
  takes `3!·2^{n−3l}` targets per relation;
- each target costs `Θ(2^{2l})` with any oracle that enumerates the
  factor base, this one included.

Multiply:

```text
  2^l  ·  3!·2^{n−3l}  ·  Θ(2^{2l})  =  Θ(2^n)
```

**independent of `l`.**  A larger factor base needs fewer tries, and
makes each try proportionally more expensive, and the two cancel
exactly.  So index calculus built on an enumeration oracle costs `2^n`
however fast the inner loop is — against `2^{n/2}` for Pollard rho.

The benchmark measures this rather than asserting it.  Random targets,
`s/target` measured, projected collection cost `= 2^l / rate ×
s/target`, divided by `2^n`:

| `l` | `n` | samples | decomposes | `1/3!` | s/target | projected collection | `/ 2^n` |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 5 | 15 | 4 000 | 20.9% | 16.7% | 3.2e−4 | 4.9e−2 s | 1.5e−6 |
| 6 | 18 | 4 000 | 13.5% | 16.7% | 1.7e−3 | 7.9e−1 s | 3.0e−6 |
| 7 | 21 | 2 881 | 15.6% | 16.7% | 6.9e−3 | 5.7e0 s | 2.7e−6 |
| 8 | 24 | 638 | 14.4% | 16.7% | 3.1e−2 | 5.6e1 s | 3.3e−6 |
| 9 | 27 | 132 | 14.4% | 16.7% | 1.5e−1 | 5.4e2 s | 4.0e−6 |
| 10 | 30 | 32 | 15.6% | 16.7% | 6.4e−1 | 4.2e3 s | 3.9e−6 |
| 11 | 33 | 9 | 33.3% | 16.7% | 2.3e0 | 1.4e4 s | 1.7e−6 |
| 12 | 36 | 3 | 33.3% | 16.7% | 1.1e1 | 1.3e5 s | 1.9e−6 |

Across the well-sampled rungs (`l = 6…10`, 32 to 4 000 targets each)
the last column is **flat within 1.5×** while `2^n` grows 4 000-fold.
The predicted decomposition rate holds too: 13.5–15.6% measured against
16.7%.  The `l = 11, 12` rows rest on 9 and 3 samples and their rates
are noise — they are in the table for the timings, not the rates.

So: the oracle got roughly a thousand times faster and the attack's
exponent did not move. That is not a disappointment, it is the shape of
the problem.  The SAT route failed for the same reason at one further
remove — a CDCL solver on this encoding is an enumerator too, just a
slower one.

What the speed does buy is a usable experiment.  `l = 12` used to be a
day and a half per target and is now 15 seconds, which is the range
where the only thing that *would* change the conclusion — a
sub-`2^{2l}` decomposition oracle — could actually be tested.

## Correctness

The fast path is worth nothing if it answers a different question, so
every claim it makes is checked against something independent.

```bash
cargo test --release --lib cryptanalysis::semaev_decomp
```

- `decomposition_agrees_with_exhaustive_search` — the decisive one.
  On all 60 upstream corpus instances, pairs-and-solve must reach the
  same verdict as walking every sorted triple, and any witness it
  returns must genuinely satisfy the symmetrised `S₄`.  That includes
  `n19l6-19-U`, which is satisfiable despite its label (see
  `RESEARCH_SAT_SEMAEV.md`), so a naive implementation that trusted
  filenames would fail here.
- `decomposition_agrees_with_enumeration_on_random_targets` — the
  corpus only reaches `l = 6`, and the polynomial machinery is worked
  much harder above that.  40 random targets at `l = 7`, `n = 21`,
  every verdict compared against triple enumeration, and the test
  asserts that both verdicts actually occurred so it cannot pass
  vacuously.
- `quartic_matches_the_symmetrised_form` — the rearrangement into a
  quartic in `X₃`, against the twelve-term symmetrised `S₄` evaluated
  directly through the general field, over 64 000 triples.
- `field_arithmetic_matches_the_general_implementation` — `mul`, `sqr`,
  `inv` and `batch_inv` against `F2mElement` across `n = 15, 17, 18,
  21, 24, 27, 30, 33`, because the reduction table's width and the
  carry-less product's shape both depend on `n`.  Batch inversion is
  checked with zeros in the slice.
- `subspace_polynomial_vanishes_on_the_subspace` — `L_V` must vanish on
  `V` and nowhere else.  Swept completely for `n = 15, 21`; for
  `n = 36` all of `V` plus 20 000 random points outside it.

## What's still open

- **A sub-`2^{2l}` oracle** is the only thing that would matter, and
  nothing here is one.  The route worth trying is Gröbner basis (F4/F5)
  on the descended system for fixed `X₁`, which is a bivariate problem
  in `2l` unknowns — the premise being to exploit the algebra rather
  than search around it.  `src/cryptanalysis/groebner_f4.rs` exists.
- **Frobenius-stable factor bases.**  On the Koblitz curve
  `y² + xy = x³ + x² + 1`, `φ(x,y) = (x², y²)` is an endomorphism, so a
  relation for `x_R` gives relations for `x_R^{2^i}` free — a factor-`n`
  saving in relation collection.  It needs a factor base closed under
  squaring, which `⟨1, z, …, z^{l−1}⟩` is not; the subfield
  `F_{2^l} ⊂ F_{2^{3l}}` is (and `l | 3l` always holds on this ladder).
  That is a factor of `n`, not a change of exponent, but it is free and
  it changes the corpus comparison, which is why it is not done here.
- **Parallelism.**  The pair loop is embarrassingly parallel over `X₁`
  and nothing in it shares state.  Four cores would be ~2 bits of `l`.
  Recorded, not done: it does not interact with anything above.
- **`b = 1` only.**  Inherited from `binary_semaev_s4`: `S₄` does not
  depend on `a₂` but does depend on `b`, and the `b`-powers are folded
  into the twelve constants.
- **Constant factor.**  About 1 µs per pair at `l = 8`, roughly 190
  field multiplications.  Two thirds is the `l` squarings modulo the
  quartic; a `t⁴`/`t⁶` reduction table would cut that by a quarter, and
  Gray-coding the pair loop makes the quartic's constant and leading
  coefficients pure XOR updates.  Together maybe 25%, which is why
  neither is done.

## References

- **P. Gaudry**, *Index calculus for abelian varieties of small
  dimension and the elliptic curve discrete logarithm problem*,
  J. Symbolic Computation 44 (2009) — the pairs-and-solve structure.
- **C. Diem**, *On the discrete logarithm problem in elliptic curves*,
  Compositio Math. 147 (2011) — factor bases as subspaces.
- **R. Lidl, H. Niederreiter**, *Finite Fields*, §3.4 — linearized and
  subspace polynomials.
- **J.-C. Faugère, P. Gaudry, L. Huot, G. Renault**, *Using symmetries
  in the index calculus for elliptic curves discrete logarithm*,
  J. Cryptology 27 (2014) — the symmetrised `S₄` this consumes.
- **P. L. Montgomery**, *Speeding the Pollard and elliptic curve methods
  of factorization*, Math. Comp. 48 (1987) — batch inversion.
- **M. Trimoska et al.**, `EC-Index-Calculus-Benchmarks` — the corpus
  the correctness gate is built on.
