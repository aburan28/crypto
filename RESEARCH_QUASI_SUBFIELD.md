# Quasi-subfield polynomials over `F_{2^n}`: a census, a criterion, and where the barrier actually is

**Module:** `src/cryptanalysis/quasi_subfield.rs`
**Demos:** `cargo run --release --example quasi_subfield_census`,
`cargo run --release --example quasi_subfield_reach [n_max]`

**The paper.** M.-D. Huang, M. Kosters, C. Petit, S. L. Yeo, Y. Yun,
*Quasi-subfield polynomials and the elliptic curve discrete logarithm
problem*, J. Math. Cryptol. 14(1):25–38, 2020.  Open access,
`doi:10.1515/jmc-2015-0049`.

They build index-calculus factor bases from the roots of
`L(X) = X^{q^{n0}} − λ(X)` when `L` splits completely in `F_{q^n}` and
`deg λ` is small; the subfield case is `λ(X) = X`.  Definition 3.1 calls
`L` **quasi-subfield** when `log_q(deg λ) < n0²/n`.  Their paper closes
with:

> It remains an open problem to find (or rule out) the existence of
> quasi-subfield polynomials where `deg λ` is small enough to improve on
> the best (generic) algorithms for ECDLP.  A question of particular
> interest is whether the bound on `deg λ` provided by Lemma 4.1 is
> tight: in fact removing the term `n mod n0` in this bound would show
> that our approach cannot beat generic algorithms.

This note answers the second question, measures the first, and then
argues the first was the wrong place to look.

## 1. The search, made finite

Take `q = 2`.  A monic 2-linearised `L(X) = X^{2^{n0}} + Σ_{i<n0} c_i X^{2^i}`
acts `F_2`-linearly on `F_{2^n}`, and **splits completely there iff that
map has kernel of dimension exactly `n0`**.  So "splits completely" is a
rank condition on an `n × n` matrix over `F_2`, and the whole question is
a finite search over `(c_0, …, c_{n0−1})`.

`deg λ = 2^j` with `j = max{i : c_i ≠ 0}`, so Definition 3.1 reads:
**every coefficient from index `⌈n0²/n⌉` up must vanish.**

Counting heuristically: the space has `2^{t·n}` points with
`t = ⌈n0²/n⌉`, and a random `F_2`-linear map has an `n0`-dimensional
kernel with probability of order `2^{−n0²}`.  Since `t·n ≈ n0²`, the
expected number of solutions is of order **one** — the question sits
exactly on the first-moment threshold, which is why a statistical
argument does not settle it.

## 2. The `n mod n0` term cannot be removed

Lemma 4.1 gives `log₂ deg λ ≥ (n0 − (n mod n0)) / ⌊n/n0⌋`.  Deleting the
`n mod n0` term raises this to `n0 / ⌊n/n0⌋`.

Over `F_{2^7}` the exhaustive census finds quasi-subfield polynomials
that are **not** subfield polynomials, and they sit strictly below the
stripped bound:

| `n` | `n0` | witness `λ` | `j = log₂ deg λ` | Lemma 4.1 | stripped |
|---:|---:|---|---:|---:|---:|
| 7 | 3 | `X² + X` | 1 | 1.000 | 1.500 |
| 7 | 4 | `X⁴ + X² + X` | 2 | 1.000 | 4.000 |

The second is checked by hand as well as by machine.  As an operator
`L = σ⁴ + σ² + σ + 1` with `σ` the Frobenius, and
`t⁴ + t² + t + 1 = (t + 1)(t³ + t² + 1)`, both factors of `t⁷ + 1`, so
the kernel has dimension 4 and `L` splits.  `deg λ = 4`, and
`log₂ 4 = 2 < 16/7`.

**So the stripped bound is false, and the `n mod n0` term is needed.**
The paper's question of particular interest resolves in the direction
that keeps their approach alive rather than killing it.  Test:
`the_hand_checked_witness_is_a_genuine_quasi_subfield_polynomial`.

## 3. Existence is governed by `t^n − 1`, not by counting

An `F_2`-subspace `V ⊆ F_{2^n}` is Frobenius-stable exactly when
`V = ker g(σ)` for a divisor `g` of `t^n − 1`, and then the subspace
polynomial of `V` **is** `g(σ)`.  Its coefficients are those of `g`, so
they lie in `F_2`, and the quasi-subfield condition becomes a statement
about **gaps in divisors of `t^n − 1`**:

```text
    g = t^{n0} + (terms of degree < n0²/n).
```

Measured claim: on every cell the exhaustive census can reach, the
divisor criterion and the census **agree on existence** — every
quasi-subfield polynomial over `F_{2^n}` is accounted for by a
Frobenius-stable subspace.  Test:
`divisor_prediction_matches_the_exhaustive_census`.

The first-moment heuristic is badly wrong in the safe direction:

| `n` | `n0` | searched | first-moment expectation | actually found |
|---:|---:|---:|---:|---:|
| 9 | 4 | `2^18` | 4 | **0** |
| 10 | 4 | `2^20` | 16 | **0** |
| 11 | 4 | `2^22` | 64 | **0** |
| 13 | 4 | `2^26` | 1024 | **0** |
| 13 | 5 | `2^26` | 2 | **0** |

At `n = 11` and `n = 13` the value `2` is a primitive root, so `t^n − 1`
has no divisor of intermediate degree at all, and the count is zero
however large the heuristic says it should be.

This matters practically: testing divisors costs `2^{j_max+1}` trial
divisions instead of the census's `2^{(j_max+1)·n}` rank computations,
which is what lets the question be decided at cryptographic `n`.

## 4. The census at prime `n`

Scanning every odd `n ≤ 600` and every `n0 < 64`: 560 cells admit a
quasi-subfield polynomial, but 483 of those are the subfield case
`n0 | n`.  At **prime** `n`, which is the Koblitz setting, only 13
non-subfield cells exist in the whole range:

| `n` | `n0` | `j` | `n0/n` |
|---:|---:|---:|---:|
| 3, 5, 11, 13, 17, 19, 23 | `n − 1` | `n − 2` | 0.67 → 0.96 |
| 7 | 3, 4, 6 | 1, 2, 5 | 0.43, 0.57, 0.86 |
| 31 | 15, 16 | 7, 8 | 0.48, 0.52 |
| **73** | **9** | **1** | **0.123** |

The `n0 = n − 1` family is the trace hyperplane and is useless: the
factor base is half the field.  The only cell with a ratio small enough
to matter for index calculus is `n = 73`, `n0 = 9`, from
`g = t⁹ + t + 1`, a degree-9 factor of `t^73 + 1` — giving a 512-element
factor base in `F_{2^73}` with `deg λ = 2`.  It clears Definition 3.1
with almost no room: `1 < 81/73 = 1.1096`.

## 5. Where the barrier actually is

The paper's Remark 3.1 optimises its own complexity to
`Õ(q^{n(1 − α + 4.876 α²)})` and observes the exponent bottoms out near
`0.95`, "which beats brute force algorithms".  It never claims to beat
generic ones, and it cannot:

```text
    min_α (1 − α + 4.876 α²) = 1 − 1/(4 · 4.876) = 0.9487 ,
```

attained at `α ≈ 0.1025`, against the generic exponent `0.5`.  **The
minimum does not depend on `deg λ` at all.**  So no quasi-subfield
polynomial, however small its `deg λ`, beats generic algorithms inside
this framework.  Beating them needs the solving-exponent constant `c` in
`1 − α + c α²` to satisfy `min = 1 − 1/(4c) < 1/2`, i.e. **`c < 0.5`**,
against `c = 4.876` today — the `4.876` coming from Rojas's algorithm
for the polynomial system, not from the factor base.

That reframes the open problem.  The paper asks whether `deg λ` can be
made small enough; the arithmetic says `deg λ` is not the binding
constraint.  §4 adds that the supply of candidates is thin anyway: one
useful cell below `n = 600` at prime `n`.

## 6. What this does not settle

- Characteristic 2 only, and `n0 < 64`, which is where the `u128`
  polynomial arithmetic stays exact.  Odd characteristic is untouched.
- §3's agreement between the divisor criterion and the census is
  measured on the cells brute force can reach, not proved.  It is a
  conjecture with evidence, and the natural proof route — that a
  quasi-subfield gap forces Frobenius stability — is not attempted here.
- §5 takes the paper's Remark 3.1 at face value.  It is an arithmetic
  observation about their stated bound, not an independent lower bound,
  and it says nothing about approaches outside that framework — the
  paper's own list of generalisations (rational `λ`, an isogeny map for
  `L`, large-prime and unbalanced variations) is untested here.
- A regression worth recording: the first version of the reach sweep
  formed `t^n` as a `u128` bitmask, which wraps silently once `n ≥ 128`
  and manufactured "hits" at every large `n`.  Divisibility is now
  tested as `t^n ≡ 1 (mod g)` by square-and-multiply, so `t^n` is never
  represented; the regression test
  `divisibility_is_exact_past_the_word_width` pins it at
  `n = 163, 233, 283, 409, 571`.
