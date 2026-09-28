# Choosing the factor base for the solver, not for the yield

**Harness:** `cargo run --release --example groebner_base_sweep -- --out DIR`
**Frozen evidence:** `research/groebner_base_sweep_20260915/`
**Related:** [`research/notes/ecc2k130/RESEARCH_GROEBNER_STAGE.md`](../ecc2k130/RESEARCH_GROEBNER_STAGE.md) (the stage
being priced), [`research/notes/ecc2k130/RESEARCH_KOBLITZ_INDEX_CALCULUS.md`](../ecc2k130/RESEARCH_KOBLITZ_INDEX_CALCULUS.md)
(the pipeline), [`research/notes/index-calculus/RESEARCH_QUASI_SUBFIELD.md`](RESEARCH_QUASI_SUBFIELD.md) (which
subspaces exist), [`research/notes/ecc2k130/RESEARCH_ECC2K130_DECOMPOSITION.md`](../ecc2k130/RESEARCH_ECC2K130_DECOMPOSITION.md)
(the counting boundary, and the trace-zero yield factor this note confirms).

**The question.**  Summation polynomials and the Weil restriction are written
*over the factor-base subspace*: the oracle's system has `m·ℓ` Boolean unknowns
where `ℓ = dim V`, so the base is not only what the relations are written in, it
is what the solver has to solve.  Does the base the repository selects minimise
what the solver then does?

**Bottom line.  No — at `K_1/2^15` the selection objective ranks the candidate
bases almost exactly backwards.**  `koblitz_factor_base_search` minimises
expected *trials*; over five admissible invariant subspaces of that curve the
Spearman correlation between that ranking and measured whole-run cost is
`−0.90`, and the trials-optimal base costs **`12.12×`** the cheapest one over 60
complete, verified discrete logarithms.  At the other degrees measured the two
objectives agree or their difference is inside sampling error.  **And at
`n = 131` the question does not arise: the invariant-subspace lattice has two
members and neither is usable.**

## 0. Boundaries, stated before measuring

**The reference.**  Two of them, because there are two incumbent selectors:

- the **pipeline default** — `build_frobenius_factor_base(kc, 0)`, the first
  largest-degree factor of `x^n − 1`, used when no recipe is supplied;
- the **trials optimum** — the base minimising
  `T(F) = (U(F) + 1) / p_m(F)`, columns over coverage, which is what
  `koblitz_factor_base_search` scores.  Its own module says why this is only
  half the cost: *"solving cost per trial is a separate axis that depends on the
  oracle, and the report carries the standard proxies … alongside so a caller
  can trade them off."*  The proxies are `|F|^{m−1}` enumeration work and a SAT
  variable count; neither is what the Gröbner oracle costs.

**The floor.**  Unchanged and untouchable from here: the counting bound of
`research/notes/ecc2k130/RESEARCH_ECC2K130_DECOMPOSITION.md` §5, `m·2^131` group operations for the
whole attack *independently of `ℓ`*, because a larger base needs proportionally
fewer targets and makes each target proportionally dearer.  Choosing a better
base inside the family cannot move it, so every ratio to the floor below is
flat.

**The unit.**  64-bit word XORs in the Macaulay elimination, as in
`research/notes/ecc2k130/RESEARCH_GROEBNER_STAGE.md` — now also reported per whole run by
`ic run --json` as `counts.f4_word_ops`.

**Falsification target.**  The thread is a success if a base exists whose
measured collection cost is lower than the trials-optimal base's, at the same
`n` and `m`, on independent target samples, with the discrete logarithm actually
recovered and verified on every run.  It fails if the two objectives agree
wherever the lattice offers a real choice.  Inadmissible: comparing bases on
different curves, targets or summand counts; quoting stage cost where a complete
run was affordable; treating a coverage difference inside its own binomial
interval as a property of the base.

## 1. What choice exists at all

An `F_2`-subspace of `F_{2^n}` is Frobenius-stable exactly when it is
`ker g(σ)` for a divisor `g | t^n − 1` (`research/notes/index-calculus/RESEARCH_QUASI_SUBFIELD.md` §3), so the
candidates *are* the divisors and the available dimensions are the sums of
cyclotomic coset sizes:

| `n` | `ord₂(n)` | coset sizes | invariant dimensions `0 < ℓ < n` |
|---:|---:|:--|:--|
| 9 | 6 | 1×1, 2×1, 6×1 | 1, 2, 3, 6, 7, 8 |
| 13 | 12 | 1×1, 12×1 | 1, 12 |
| 15 | 4 | 1×1, 2×1, 4×3 | 1 … 14 (14 of them) |
| 17 | 8 | 1×1, 8×2 | 1, 8, 9, 16 |
| 23 | 11 | 1×1, 11×2 | 1, 11, 12, 22 |
| 31 | 5 | 1×1, 5×6 | 1, 5, 6, 10, … , 30 (12) |
| 73 | 9 | 1×1, 9×8 | 1, 9, 10, 18, … (16) |
| 127 | 7 | 1×1, 7×18 | 1, 7, 8, 14, … (36) |
| **131** | **130** | **1×1, 130×1** | **1, 130** |
| 163 | 162 | 1×1, 162×1 | 1, 162 |

**At the challenge degree there is nothing to choose.**  `2` is a primitive root
modulo `131`, so `t^131 − 1 = (t+1)·(irreducible of degree 130)` and the only
Frobenius-invariant subspaces are dimension `1` — two abscissae — and dimension
`130`, the trace hyperplane, half the field.  A factor base "optimised for
summation polynomials" in this family does not exist at `n = 131`; the families
that remain are non-invariant subspaces, which forfeit the `n`-fold Frobenius
collapse of the relation columns, and unions of Frobenius orbits, which
`research/notes/ecc2k130/RESEARCH_ECC2K130_DECOMPOSITION.md` §6 already prices at `2^125.55` and
`2^44.5` stored points.  `n = 127` and `n = 73`, one and two off the challenge
degree, have rich lattices; `131` does not.  That is a property of the
challenge's degree, not of the method.

**A free structural bit.**  The trace is `Tr = h(σ)` with `h = (t^n+1)/(t+1)`,
so `ker g(σ) ⊆ ker Tr` exactly when `g | h` — that is, **when `(x+1)` does not
divide the divisor**.  `research/notes/ecc2k130/RESEARCH_ECC2K130_DECOMPOSITION.md` §3.1 derives that a
base inside `ker Tr` yields twice, since its sums land in a subgroup of index 2.
That criterion is an algebraic fact, not a fit, and direct trace computation on
every candidate at `n = 9, 15, 17, 23, 31` agrees with it row for row.

The sweep shows the factor of two paying: adding the factor `(x+1)` to a divisor
doubles the subspace and takes it out of `ker Tr`, and at `K_1/2^15` the three
such pairs go `67 → 133` points for `36% → 48%` coverage, `77 → 123` for
`44% → 46%`, and `97 → 163` for `12% → 12%`.  Twice the base for nothing, half
of it, and nothing again — the trace-zero base is buying back most of what a
whole extra dimension would give, and it is decided by a divisibility with no
solving at all.

## 2. The table

`C = (U + 1)/p · E[word XORs per target]` — the expectation over *all* targets,
refutations included, because a trial is paid whether or not it succeeds.  `T`
is the same with the middle factor set to `1`.  Coverage is pooled over two
independent target samples (seeds 1 and 2) into one binomial with a Wilson 95%
interval; `C`'s band is `C` evaluated at the ends of it.

**`K_1/2^15`, `m = 2`, 128 targets per base** (the six cheapest of sixteen):

| divisor | `ℓ` | `\|F\|` | `U` | `Tr=0` | coverage (95% CI) | ops/target | `T` | `C` | `C`/best | ratio to floor |
|:--|--:|--:|--:|:-:|--:|--:|--:|--:|--:|:--|
| `1·4` | 6 | 77 | 5 | yes | 44% (36–53%) | 1.13e5 | 13.71 | **1.48e6** | 1.00× | flat |
| `1·3` | 6 | 67 | 4 | yes | 36% (28–45%) | 1.19e5 | 15.24 | 1.66e6 | 1.12× | flat |
| `0·1·3` | 7 | 133 | 7 | no | 48% (40–57%) | 3.28e5 | 16.52 | 5.41e6 | 3.66× | flat |
| `3·4` | 8 | 251 | 10 | yes | 99% (96–100%) | 9.42e5 | 11.09 | 1.04e7 | 7.07× | flat |
| `2·4` | 8 | 281 | 11 | yes | 99% (96–100%) | 1.15e6 | 12.09 | 1.39e7 | 9.39× | flat |
| `2·3` ← **trials pick** | 8 | 211 | 8 | yes | 85% (78–90%) | 1.45e6 | **10.57** | 1.53e7 | **10.35×** | flat |

**Class: engineering** by §3 of `AGENTS.md` — a cost fell, the ratio to the
floor did not move — with an **accounting** component: the objective being
corrected was under-counting, not mis-measured.

The mechanism is visible in the two middle columns.  Coverage saturates: past
`ℓ = 8` every target decomposes and there is nothing left to win.  Cost does
not: the system has `m·ℓ` unknowns, and ops/target rises `1.13e5 → 3.28e5 →
1.45e6` from `ℓ = 6` to `8`.  `T` sees only the columns and the coverage, so it
walks up the `ℓ` ladder collecting a factor of two in yield while paying an
order of magnitude in solving — and it cannot see the second number at all.

The other rungs, for completeness:

| curve | trials pick | cost pick | penalty | at the worst end of both CIs | verdict |
|:--|:--|:--|--:|--:|:--|
| `K_0/2^9`, `m = 2` | `2` | `2` | 1.00× | — | agree |
| `K_0/2^9`, `m = 3` | `2` | `2` | 1.00× | — | agree |
| `K_1/2^15`, `m = 2` | `2·3` | `1·4` | **10.35×** | **7.99×** | **disagree** |
| `K_1/2^17`, `m = 2` | `2` | `2` | 1.00× | — | agree |
| `K_1/2^23`, `m = 2` | `1` | `1` | 1.00× | — | agree |
| `K_0/2^31`, `m = 2` | `1·3·5` | `1·2·3` | 1.19× | 0.56× | **not established** |

At `n = 31` the six `ℓ = 15` subspaces differ by up to `1.8×` in `C`, but 24
targets per base put every coverage inside everyone else's interval, so **no
same-dimension preference is established there** and the row says so.  The `n =
15` disagreement is the one that survives its own error bars, and it survives
them by `7.99×`.

## 3. The whole pipeline, on the pair the objectives disagree about

Stage cost is a model.  The five admissible bases were therefore run end to end:
`ic run` at `K_1/2^15`, twelve random known-answer targets each (seeds 11…133),
identical curve, subgroup and solver, every logarithm recovered and verified.

| divisor | verified | `\|F\|` | columns | trials | F4 reductions | **word XORs** | ratio | wall |
|:--|:-:|--:|--:|--:|--:|--:|--:|--:|
| `1·3` | 12/12 | 67 | 2 | 53 | 365 | **6,226,452** | 1.00× | 0.093 s |
| `1·4` | 12/12 | 77 | 1 | 56 | 374 | 6,649,451 | 1.07× | 0.093 s |
| `3·4` | 12/12 | 251 | 5 | **42** | 1,053 | 37,392,007 | 6.01× | 0.289 s |
| `2·4` | 12/12 | 281 | 6 | 63 | 1,742 | 61,979,392 | 9.95× | 0.421 s |
| `2·3` | 12/12 | 211 | 4 | 52 | 2,157 | **75,493,337** | **12.12×** | 0.496 s |

Read the trials column and the word-XOR column together.  `3·4` needs the
**fewest trials of any base** and costs `6×` the cheapest; `2·3`, the base the
trials objective selects, costs `12.12×`.  Over these five bases the rank
correlation between `T` and measured whole-run cost is `ρ = −0.90`.  The stage
model predicted `10.35×` for `2·3` against `1·4`; the complete runs measured
`11.35×`, so the model is good to 10% on the quantity it is used for.

Two honest asymmetries, both against the finding's favour and neither large
enough to change it: the pipeline merges relation columns by cofactor
projection, so its column counts (2, 1, 5, 6, 4) are smaller than the
signed-orbit counts `U` the stage metric uses — `C` therefore overstates the
relation requirement, differently per base; and `1·3` and `1·4` are tied within
noise, so "the cheapest base" is that pair, not either one.

## 4. What this changes, and what it does not

- **Selection should score `C`, not `T`.**  The instrument exists now
  (`groebner_base_sweep` reuses the stage profiler), it costs one sweep per
  curve, and where it disagrees with `T` it disagrees by an order of magnitude.
  The default base is worse still: at `K_1/2^15` it is the `ℓ = 4` divisor,
  which decomposes **nothing** at `m = 2` (`0/32` in
  `research/notes/ecc2k130/RESEARCH_GROEBNER_STAGE.md`'s ladder), which is what the yield search was
  built to fix — this note is about the choice it makes *among* usable bases.
  **§6 does this**: the selector now scores `C` when asked to.
- **The free bit is worth taking.**  Prefer a divisor not divisible by `(x+1)`:
  it puts the base inside `ker Tr` and doubles the yield, decided without
  solving anything.
- **Nothing here bears on ECC2K-130.**  At `n = 131` the lattice is `{1, 130}`.
  The attack cost stays where `research/notes/ecc2k130/RESEARCH_ECC2K130_DECOMPOSITION.md` §5 put it,
  `m·2^131` and `2^71.77` times rho, and an order of magnitude in a toy-scale
  stage constant is `2^3.6` against `2^131`.
- **Not measured here:** base construction (`|F|` point enumeration), lifting,
  filtering and the relation linear algebra.  All of them grow with `|F|`, so
  pricing them would move every row in the same direction as the finding —
  toward the smaller base — but they are not in these numbers and the ratios
  above are the Gröbner stage's alone.
- **Not claimed:** any change to the exponent, the rho crossover, or the
  full-DLP complexity.  Per `AGENTS.md` §8 the end-to-end table is a matched
  toy-scale comparison with verified answers, not an ECDLP speed-up claim at
  deployed sizes.

## 5. Reproducing

```bash
cargo build --release --example groebner_base_sweep
./target/release/examples/groebner_base_sweep --degrees 9,15,17,23,31 --target-seed 1 \
    --out research/groebner_base_sweep_20260915/sample-a
./target/release/examples/groebner_base_sweep --degrees 9,15,17,23,31 --target-seed 2 \
    --out research/groebner_base_sweep_20260915/sample-b
python3 research/groebner_base_sweep_20260915/analyse.py \
    research/groebner_base_sweep_20260915/sample-a \
    research/groebner_base_sweep_20260915/sample-b --output /tmp/base-sweep.json

# rank the same lattice by trials, then by measured solving cost (§6)
./target/release/ic --json search --degree 15 --curve-a 1 --summands 2 --family divisor \
    --min-dimension 3 --max-dimension 8 --no-prune --no-saturate --validate-top 0 --seed 1
./target/release/ic --json search --degree 15 --curve-a 1 --summands 2 --family divisor \
    --min-dimension 3 --max-dimension 8 --no-prune --no-saturate --validate-top 0 \
    --solver groebner --solve-cost-targets 8 --seed 1

# one complete verified DLP on a chosen base
./target/release/ic run --degree 15 --curve-a 1 --solver groebner --batch 1 \
    --random-target --seed 11 \
    --factor-base research/groebner_base_sweep_20260915/end_to_end/base-1-4.json --json
```


## 6. The selector, corrected

§4 said selection should score `C`.  It now can:
`SearchOptions::solve_cost_targets` runs the Gröbner oracle on that many
census targets per candidate and ranks by
`expected_stage_ops = expected_trials × measured word XORs per target`,
exposed as `ic search --solve-cost-targets N` and as
`factor_base.solve_cost_targets` in a workflow parameter file.  `None`, the
default, leaves every earlier search scoring exactly as it did.

**What can be priced, and what cannot.**  The summation polynomial is
Weil-restricted over the subspace basis, so the system's solutions are the
whole span of that basis.  Only a base that *is* that span —
`FactorBaseDomain::LinearSubspace` — is described by its own system.  A
pruned subset, a 2-torsion saturation, a Frobenius union and an explicit
orbit set are all proper subsets carried by the SAT domain trie instead, and
asking the Gröbner oracle for them makes it enumerate the span and reject
non-members: it prices the span, not the base.  The number would be wrong
at any price, and it is not cheap either — the rejected enumeration runs to
the solver's node budget on every target, and a first cut of this selector,
which measured every candidate, did not finish a degree-15 sweep in ten
minutes where the linear candidates alone take `0.32 s`.  That comparison is
not matched and is quoted only as what the guard is for; it is also no
longer reproducible, because the guard is now in the binary.  Those
candidates
are therefore left unpriced, with the reason recorded per candidate, and
ranked below every priced one — the two scores are word XORs and trials, and
comparing them would let a trials count of `8` beat a word-XOR count of
`10^6`.  For the same reason the workflow refuses `solve_cost_targets`
together with `prune`, `saturate`, a non-linear family, or any collection
oracle but Gröbner: scoring one oracle and running another selects for the
wrong thing, and a workflow selects with no one reading the ranking.
Neither path measures under `IC_REDUCTION_CACHE` either — a memoised
reduction returns without running F4, so a warm cache would price whichever
candidates came later at a fraction of their cost — and both say so rather
than reporting a number they did not measure.  The other two cache layers
are harmless here: `word_ops` is counted inside the elimination, which a
preprocessing hit still performs.

**The table.**  Same unit as §3, 64-bit word XORs in the Macaulay
elimination, summed over twelve complete runs per base on the seeds §3 used
(`11 … 133`), same curve, subgroup and solver, every logarithm recovered and
verified — 72 runs, `72/72`.

| divisor | `ℓ` | ranked first by | verified | `\|F\|` | columns | trials | F4 reductions | **word XORs** | ratio | ratio to floor |
|:--|--:|:--|:-:|--:|--:|--:|--:|--:|--:|:--|
| `0·3` | 5 | `C` (`--solve-cost-targets`) | 12/12 | 31 | 1 | 152 | 232 | **2,405,397** | 1.00× | flat |
| `1·3` | 6 | — (§3's cheapest) | 12/12 | 67 | 2 | 53 | 365 | 9,543,849 | 3.97× | flat |
| `1·4` | 6 | — | 12/12 | 77 | 1 | 56 | 374 | 10,144,324 | 4.22× | flat |
| `3·4` | 8 | `T` (trials) | 12/12 | 251 | 5 | 42 | 1,053 | 53,895,954 | **22.41×** | flat |
| `2·4` | 8 | — | 12/12 | 281 | 6 | 63 | 1,742 | 88,192,793 | 36.66× | flat |
| `2·3` | 8 | — (§3's trials pick) | 12/12 | 211 | 4 | 52 | 2,157 | 109,566,018 | 45.55× | flat |

**Read against §3.**  Trials, F4 reductions, base sizes and column counts
are *identical* to §3's table row for row — the pipeline is the same one.
The word-XOR counts are not: every row is `1.42–1.53×` §3's, because the
decomposition system is now built through
`cryptanalysis::polynomial_reuse`, which changed what the elimination is
handed and so what it costs.  §3's absolute figures were measured before
that and are not reproducible from this tree; its ratios are — `2·3` against
`1·3` was `12.12×` there and is `11.48×` here, and the rank correlation
between `T` and measured cost over these six bases is `ρ = −0.83`.  The
table above is the
matched one for this code, and §2's stage model, which predates the same
change, should be re-run before its absolute `C` values are quoted again.

**Frozen evidence:** `research/factor_base_selector_20260915/` — both search
reports and all 72 runs.

Two things to read off it.  The corrected objective picks `0·3`, which
**§2's sweep never scored**: that sweep's dimension window started at
`ℓ = 6` and `0·3` is `ℓ = 5`, so the cheapest base at this curve was outside
the window the note measured in.  And the mechanism is the one §2 named,
running the other way: `0·3` pays **152 trials against `3·4`'s 42** — more
than three times as many — and still costs `22.41×` less, because each of
its trials solves a 10-unknown system instead of a 16-unknown one.  Trials
are the cheap axis.  That is the whole finding, and it is now what the
selector optimises.

**Unchanged.**  The floor, the ratio to it, and everything §4 says about
`n = 131`: the invariant-subspace lattice there is `{1, 130}` and an
objective cannot choose from a set with nothing usable in it.  Still not in
the score: base construction, lifting, filtering and the relation linear
algebra, all of which grow with `|F|` and so would favour the smaller base
further.  The report says so in its own `limitations` field.

## References

- P. Gaudry, *Index calculus for abelian varieties of small dimension and the
  elliptic curve discrete logarithm problem*, J. Symb. Comput. 44 (2009).
- J.-C. Faugère, L. Perret, C. Petit, G. Renault, *Improving the complexity of
  index calculus algorithms in elliptic curves over binary fields*, EUROCRYPT
  2012 — the Weil restriction over a subspace this note varies.
- S. Galbraith, S. Gebregiyorgis, *Summation polynomial algorithms for elliptic
  curves in characteristic two*, INDOCRYPT 2014 — on what multiplicative
  structure the symmetrised route would need from `V`.
- M.-D. Huang, M. Kosters, C. Petit, S. L. Yeo, Y. Yun, *Quasi-subfield
  polynomials and the elliptic curve discrete logarithm problem*, J. Math.
  Cryptol. 14 (2020).
