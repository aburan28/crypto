# Targets: six experiments on ECC2K-130 point decomposition

**Boundaries:** `scripts/ecc2k130_decomposition_targets.py`
**Frozen artefact:** `experiments/ecc2k130_decomposition_targets.json`
**Background:** [`RESEARCH_ECC2K130_DECOMPOSITION.md`](RESEARCH_ECC2K130_DECOMPOSITION.md)
**Pipeline under test:** `src/cryptanalysis/koblitz_index_calculus.rs`,
`src/cryptanalysis/semaev_decomp.rs`, the `ic` tool

Six scoped experiments, each with a boundary derived *before* anything is run, a
primary metric in one unit, and a falsifier specific enough that a run either
meets it or does not.  Written to be picked up one at a time.

The background note answered "can an ECC2K-130 point be decomposed?" with
*exists, admissible, not findable*, and closed the fix-some-summands oracle
family at `m·2^131`.  These are the six places that answer could still be wrong,
ordered by how much would change if it were.

**One of them already paid.**  E6 was written as a check on a line in the
background note and instead corrected it: the GGMP collapse needs `π(F) = F`,
not a subspace, so it *is* available at `n = 131` from a union of Frobenius
orbits, and the bottom of the family is `2^124.99` rather than `2^132.58`.  That
is `2^64.18×` rho instead of `2^71.77×` — an **accounting** correction by §3 of
`AGENTS.md`, recorded in §6.1 of the background note, and the reason the rest of
this list is worth running rather than assuming.

## The shared boundary

Everything is scored against Pollard rho with the `⟨−1⟩ × ⟨π⟩` speed-up on the
same subgroup:

```text
rho reference at n = 131:   2^60.8090        S = 0.077430
```

and the family's derived line, from the background note's §5.1:

```text
    2^l relations · 2^n/C(|F|,m) targets · C(|F|,m−1) oracle  =  m · 2^n
                                    ... and m · 2^n / n for a Frobenius-stable base.
```

Every experiment below reports **`Λ = total operations / 2^n`**, which is the
unit that makes the product law a constant (`Λ = m`, or `m/n`) and turns "did
anything move" into reading one column.  `S = ops/√r` is reported alongside so
the rows drop into `docs/index-calculus-scoreboard.html` unchanged.

---

## E1 — The scale model: a ladder of curves shaped like ECC2K-130

**Question.**  Everything at `n = 131` in the background note is derived.  Does
the product law actually hold on curves with the same structure?

**The ladder.**  Prime `n` where `2` is a primitive root, so `ord_n(2) = n−1`,
the 2-cyclotomic cosets are `{1, n−1}`, and the only invariant subspace
dimensions are `0, 1, n−1, n` — the exact obstruction that defines the ECC2K-130
instance.  Below 70 that is

```text
    n = 11, 13, 19, 29, 37, 53, 59, 61, 67          (then 83, 101, 107, 131)
```

on `K_0 : y² + xy = x³ + 1`, with a subspace factor base of dimension
`l = ⌈(n + log₂ m!)/m⌉` and `m = 3`.

**Feasibility, and the honesty line it forces.**  The predicted total is `3·2^n`,
so end-to-end runs stop at `n ≈ 37` (`2^38.6` group operations, hours) and the
rest must be **composed** from separately measured phase rates.  Report them in
different columns and never in the same one:

| class | rungs | what is measured |
|---|---|---|
| end-to-end | 11, 13, 19, 29, 37 | one planted logarithm recovered per run |
| composed | 53, 59, 61, 67 | measured yield, measured oracle rate, measured matrix, multiplied |

The composed rungs exist to extend the fit, and the end-to-end rungs exist to
validate the composition: at `n = 29` and `n = 37` the composed projection and
the real run must agree to within `1.3×` or the composition is not admissible
evidence.

**A distortion to control for.**  The odd part of `#E(F_2^n)` is not prime at
every rung — `n = 11` is `23²`, `n = 67` has a 26-bit largest prime factor.
Pose the discrete logarithm in the largest prime-order subgroup, report
`log₂` of it, and use `Λ` (which does not see the subgroup) as the primary
metric; `S` is the secondary and must be read with the subgroup column beside
it.  The faithful rungs, where the largest prime is within 4 bits of `2^n`, are
`n = 13, 19, 59`.

**Primary metric.**  `Λ = total operations / 2^n`, predicted flat at `m = 3`.

**Falsifier.**  A least-squares slope of `log₂(total)` against `n` below `0.95`
over four or more rungs, or any rung with `Λ < 0.5·m`.  Either says the product
law is not what governs these curves and the `n = 131` derivation does not
follow.

**Correctness gate.**  Every end-to-end rung recovers its planted logarithm;
every relation is re-verified in the group before entering the matrix; zero
trivial relations counted.  A run failing any of these is a bug report.

**What it costs.**  A few CPU-days for the end-to-end half.  This is the
experiment that most directly tests the background note, and the one to run
first.

---

## E2 — Is the total really flat in the factor-base dimension?

**Question.**  The product law's whole content is that `l` is not a lever.  The
note measures that as arithmetic; measure it as a run.

**Design.**  One curve from E1's ladder (`n = 59` for the faithful subgroup, or
`n = 29` if the budget is tight), `m = 3`, sweep `dim V = 8, 12, 16, 20, 24, 28`
against a saturating dimension of `20.53`.  Measure targets tried, oracle
seconds per target, relations, and linear algebra separately, and compose.

**Predicted**, from the frozen artefact: `Λ` constant up to the saturating
dimension, then rising as `2^{ml−n}`.

**Falsifier.**  Any dimension whose measured total is below half the flat line.
That would mean the cancellation between "fewer targets" and "dearer targets" is
not exact, and the exponent at 131 is not `n`.

**What it costs.**  Hours, on top of E1's harness.  Cheap, and it is the
experiment that would catch a modelling error in E1 before it propagates.

---

## E3 — Deciding versus localising: where the protection actually lives

**Question.**  This is the one with the most at stake, and the note only hints
at it.  The background note's free-oracle floor grants a hypothetical oracle
*everything* — it decides, and hands over the witness — and comes out at
`2^56.40` for `m = 4`, **below rho**.  That is what says the counting argument
does not protect this curve.  But an oracle that only answers **yes or no**, and
leaves you to find the summands, is a different and much weaker object.  How
much weaker is a derivable number, and it is large.

Turning a passed target into an explicit `m`-subset is a meet-in-the-middle, not
a search: tabulate every `⌈m/2⌉`-subset sum of the base **once**, then for each
*successful* target enumerate its `⌊m/2⌋`-subset complements and probe.  The
table is shared across every target ever tried; only the probes recur, and they
recur once per relation.  All costs are `log2` operations; the memory column is
the entries that table holds.

| `m` | free localising detector | free detector, witness by MITM | gap | witness table | beats rho? |
|---:|---:|---:|---:|---:|---|
| 2 | `2^89.25` | `2^89.45` | `2^0.20` | `2^43.15` | neither |
| 3 | `2^68.58` | `2^68.88` | `2^0.30` | `2^64.70` | neither |
| 4 | `2^56.41` | `2^68.29` | `2^11.88` | `2^44.50` | localising only |
| 5 | `2^48.44` | `2^59.75` | `2^11.31` | `2^56.97` | **both** |
| 6 | `2^42.85` | `2^62.00` | `2^19.15` | `2^45.27` | localising only |
| 7 | `2^38.74` | `2^56.93` | `2^18.19` | `2^53.82` | **both** |
| 8 | `2^35.61` | `2^59.27` | `2^23.66` | `2^45.82` | **both** |

*(The `m = 4` localising floor is the note's §5.3 free-oracle line; the two
sweeps differ by 0.01 bits of grid, `2^56.40` there against `2^56.41` here.)*

An earlier revision of this section charged the witness as a full `C(|F|, m−1)`
search and concluded that a free non-localising detector never reaches rho,
bottoming out at `2^73.29`.  **That was wrong by up to 16 bits, and the
conclusion it carried was wrong too.**  A free detector alone *does* dip below
rho — at `m = 5` (`2^−1.06`), `m = 7` (`2^−3.88`) and `m = 8` (`2^−1.54`) — so
deciding for free is, at some `m`, already enough.

The distinction survives, in a weaker and more precise form:

- **Localising crosses one summand earlier and never comes back up.**  It is
  under rho from `m = 4` onward, monotonically; the detector-only floor saws,
  because the MITM split `⌈m/2⌉ / ⌊m/2⌋` is better balanced at odd `m` than at
  the even `m` above it.  At `m = 6` a free detector is back *above* rho.
- **Localising is free of memory; deciding alone is not.**  Every detector-only
  cell that reaches rho needs a witness table of `2^45`–`2^57` entries — at
  `m = 5`, `2^56.97` entries against a `2^59.75` total, so the table *is* the
  algorithm.  Priced against BSGS at the same memory (§5.3 of the note), that is
  a much less interesting object than the row's total suggests.
- **The gap is still 11 to 24 bits from `m = 4` up**, and at the useful end of
  the range it is the difference between a floor near rho and one well under it.

So the honest statement is not "the whole protective margin sits in that gap".
It is: *localisation is worth 11–24 bits and removes a `2^45`-plus memory
requirement, and a free detector that cannot localise still gets within a couple
of bits of rho at the right `m`.*  No statement in the literature about
decomposition-oracle hardness distinguishes the two, and the second half of that
sentence is the part that should worry anyone proposing a detector.

Both lines are **floors**: derived from the product law, not measured, and not
attained by any algorithm in this repository or any other one known to it.

**Design.**  Two measurements, both at toy sizes on E1's harness:

1. For each oracle the repository has — pairs-and-solve, matrix-F4, SAT — measure
   the cost of a **sub-base query** `W ⊆ V` as a function of `dim W`, and the
   multiplier between "decide" and "produce a witness".  Predicted for
   pairs-and-solve: the query cost is `C(2^{dim W}, m−1)`, i.e. it localises for
   free, and bisection costs about `2×` a single full query.
2. Take any *proposed* cheap detector — a resultant-vanishing test, a trace
   condition, a partial Gröbner refutation at fixed degree — and measure the same
   two things.  A detector that does not localise is not thereby harmless — it
   buys `2^56.93` at `m = 7` — but it has to carry a `2^53.82`-entry table to do
   it, so the measurement to make is the *joint* one: cost **and** memory, scored
   against the BSGS line at that memory, not cost alone.

The Nagao/Riemann–Roch encoding landed on main by #316
(`research/nagao_relations/`) is the natural second case.  It is a different
*relation representation* rather than a subspace-membership oracle — it searches
for a function in `L(4O)` and reads the decomposition off its zeroes — so the
question E3 asks lands on it unchanged and is worth asking early: does a
coefficient search restricted to a sub-base cost less than the same search on
the whole base, in proportion, or not at all?  That note is explicit that it is
"implementation and toy correctness work, not an ECDLP speedup result", so
nothing here is scored against it yet; E3 is the frame it would be scored in.

**Falsifier / what would count.**  For the localising line: a detector whose
full-base query is cheaper than `C(|F|, m−1)` *and* whose sub-base query on `W`
costs less than the full-base query times `(|W|/|F|)^{m−1}`.  That combination
moves the `2^56.40` line into reach.  For the detector-only line, which this
revision showed is already under rho at `m ∈ {5, 7, 8}`: a cheap full-base
detector whose witness extraction beats the meet-in-the-middle — or, going the
other way, a proof that the MITM table cannot be shared across targets, which
would put the `2^73.29` line back.

**Why this one is worth designing carefully.**  Both hypothetical oracles look
identical when stated as "a fast decomposition oracle", and they differ by up to
`2^37.68`.  Anyone proposing one should be asked which it is, first.

---

## E4 — Do large primes move the product law, or relabel it?

**Question.**  Large primes are the one classical index-calculus lever the
background note does not price.  Allow the last summand to land anywhere in a
larger subspace `V' ⊃ V`, store the partial relation under that large prime, and
pair partials off by birthday.

**The model, and its guard.**  With `dim V = l`, `dim V' = l'`:

```text
    partial relations per target  =  C(2^l, m−1) · 2^{l'} / (m · 2^n)
    partials needed               =  2^{(l + l' + 1)/2}            (to yield 2^l full ones)
    collection                    =  m · 2^{n + (l − l' + 1)/2}
    memory                        =  2^{(l + l' + 1)/2}
```

The guard is **partial yield ≤ 1 per target**.  Without it `l' → n`, the oracle
stops filtering anything, and what is left is a birthday search on differences of
targets — a generic algorithm, which has to be scored as one and loses to rho.
Guarded, the derived optima are `2^70.50` at `m = 2` (store `2^66.50`), `2^72.31`
at `m = 3`, `2^73.73` at `m = 4`; memory-capped they run

| store | large-prime variant | BSGS, same memory |
|---:|---:|---:|
| `2^30` | `2^107.58` | `2^99.00` |
| `2^40` | `2^97.58` | `2^89.00` |
| `2^50` | `2^87.58` | `2^79.00` |
| `2^60` | `2^77.58` | `2^69.00` |

**Prediction.**  A genuine and large improvement on `m·2^131` — sixty bits at the
unguarded optimum — that nonetheless loses to plain baby-step giant-step at every
memory budget by a near-constant `2^8.6`, and to rho by more.  Class:
**relabelling**, by §3 of `AGENTS.md`.

**Falsifier.**  A guarded cell below the BSGS line at the same memory.

**Caveat stated in advance.**  The model above is crude on purpose: it does not
price collision bookkeeping, and it does not check that the full relations
produced by pairing are **independent** over the factor base.  Independence is
the most likely place for it to be wrong, and the experiment must measure the
rank of the relation matrix, not assume it.

---

## E5 — The yield distribution, not just its mean

**Question.**  Every extrapolation in the background note rests on
`λ = C(|F|,m)/#E` and a Poisson tail.  The twelve measured cells came in at
`0.95` to `1.38` times the prediction, eleven of them above `1.0`.  That is the
signature of **under-dispersion** — the decompositions spread more evenly across
targets than independence would give — and if it is real it is a small
security-negative correction that nobody has written down.

**Design.**  On the toy rungs where the full curve is enumerable (`n ≤ 19`),
compute the **exact** decomposition count for every target rather than a yes/no,
and test the index of dispersion `Var/mean` against its null value of `1`.  Its
standard error is `√(2/T)`, so resolving a 20% effect at 3σ needs **450 targets
per cell** — well inside budget.

**Primary metric.**  `Var/mean` per `(n, m, l)` cell.

**Falsifier.**  `|Var/mean − 1| > 0.2` at 3σ on three or more cells, with a
consistent sign.  That would mean the Poisson tail is the wrong model and the
`m·l ≥ n + log₂ m!` saturation threshold needs a correction term.

**Bound on what it can matter.**  A constant, not an exponent: the mean is what
the product law consumes, and E5 can only move the constant in front.  It is on
this list because it is cheap and because the background note's central table
quietly depends on it.

---

## E6 — Frobenius-stable bases without a subspace anywhere

**Status: this one already fired, and corrected the background note.**

**Question.**  The classification "the only invariant dimensions at `n = 131`
are `0, 1, 130, 131`" is a classification of invariant **subspaces**.  GGMP's
collapse asks only for `π(F) = F`.  Does the saving survive without the
subspace?

**Answer, derived.**  Yes.  `131` is prime, so every `x ∉ F_2` has an orbit of
size exactly `131`, and a union of `k` orbits is `π`-stable at any size.  The
collapse is therefore available at ECC2K-130, the product law becomes
`m·2^n/n`, and the bottom of the family moves from `2^132.58` to `2^124.99` —
`2^64.18×` rho instead of `2^71.77×`.  What it costs is the subspace root-find:
an orbit union has no low-degree membership polynomial, so the base must be
**materialised** (`2^44.53` points at `m = 3`, about 300 TB; `2^23.42`, about
90 MB, at `m = 6`).  **`131` being prime forbids having both the collapse and
the implicit base, not the collapse itself.**

**The concrete shape to build it out of.**  "A union of orbits" is not a
construction.  The natural one is a **Hamming-weight shell in a normal basis**:
in the basis `{β, β², β⁴, …}` the Frobenius is a cyclic shift of the 131
coordinates, so the weight-`w` shell is shift-invariant by definition, and
because 131 is prime every shell with `0 < w < 131` is exactly `C(131,w)/131`
orbits of size 131 — no fixed points to special-case.  The sizes land on what
§6 of the background note asks for almost exactly:

| `m` | orbit base wants | weight shell | shell size |
|---:|---:|---:|---:|
| 3 | `2^44.53` | `w = 9` | `2^44.43` |
| 4 | `2^33.90` | `w = 6` | `2^32.54` |
| 5 | `2^27.58` | `w = 5` | `2^28.15` |
| 6 | `2^23.42` | `w = 4` | `2^23.48` |

Shells jump three to five bits apart, so finer sizes need a union of adjacent
shells or a shift-invariant subset of one; the membership test is a popcount
either way.  **The representation is already in the repository and the factor
base is not**: `ecc2k130/FROBENIUS-NETWORK.md` applies Frobenius as a fixed
bit-permutation network over the permuted type-II ONB, and
`src/bin/ic/params.rs` carries ECC2K-130 only as an abstract normal-basis
profile with "point coordinates not imported".  So E6 needs the base built, not
merely selected.

*(A correction to something adjacent: `research/nagao_relations/README.md` on
main states that this repository "already uses normal-basis Hamming-weight
sets" as factor bases.  It does not — what exists is the ONB as a **field
representation** for the rho client.  The observation underneath that sentence
is right and is the same one as §6.1 of the background note; only the claim that
it is implemented is not.)*

**What is still to run.**  The derivation assumes `|F|/n` relations suffice,
which needs them to be **independent over the orbit unknowns**.  That is the
part most likely to be wrong, and it is measurable on E1's ladder:

- build an orbit-union base of `k` orbits at `n = 19, 29, 37`;
- collect relations, rewrite them onto orbit representatives with their `λ^{k_i}`
  weights, and measure the **rank** of the resulting matrix against `|F|/n`;
- report the realised saving as `rank-limited relations needed / (|F|/n)`.

**Falsifier.**  A realised saving differing from `n` by more than 20%, which
would mean orbit relations carry a dependency the unknown-count argument misses
— exactly the failure mode GGMP's own "multiple invariant factor bases" remark
worries about (`RESEARCH_KOBLITZ_INDEX_CALCULUS.md`, open problems: "a single
one need not yield `n` *independent* relations").

---

## Ordering, and what each one would change

| | experiment | cost | if it confirms | if it falsifies |
|---|---|---|---|---|
| 1 | **E1** ladder | CPU-days | the `n = 131` derivation is evidence, not arithmetic | the whole note's exponent is wrong |
| 2 | **E6** rank | hours on E1 | `m·2^n/n` stands as the family's floor | the seven-bit correction reverses |
| 3 | **E3** localisation | days | proposed oracles get a sharp question to answer, cost *and* memory | a cheap detector exists — localising it reaches rho from `m = 4`, and even a non-localising one reaches it at `m ∈ {5, 7, 8}` if it can carry a `2^45`-plus table |
| 4 | **E2** flatness | hours | `l` is confirmed not to be a lever | the cancellation is inexact and the exponent moves |
| 5 | **E4** large primes | days | one more classical lever priced and closed | the first sub-rho decomposition cell in this repository |
| 6 | **E5** dispersion | hours | the Poisson model is safe to extrapolate | a correction term to the saturation threshold |

E3 is the one with a real chance of surprise, because nothing in the literature
separates the two oracle strengths it separates — and because pricing the
witness correctly already moved its headline once, from "a free detector never
reaches rho" to "it reaches rho at three of seven summand counts, with a table
that is itself the size of the algorithm".  E1 is the one to run first,
because every other row on this page is scored against a law that only E1 can
confirm actually governs these curves.

## What is deliberately not on this list

- **Anything that changes the curve or the field.**  Isogeny class search (L5),
  extension fields, and quasi-subfield polynomials are each measured and closed
  in their own notes; re-opening them needs a new mechanism, not a new run.
- **A better Gröbner engine.**  The deciding degree grows as `l/2`
  (`RESEARCH_SEMAEV_DECOMPOSITION.md`), and F5's criteria remove zero reductions
  without lowering the degree of regularity, so a better implementation moves the
  constant and not the slope.
- **Factor bases that are neither subspaces nor Frobenius-stable.**  A base
  needs either a cheap algebraic membership test or materialisation; E6 covers
  the second, and the first is exactly the quasi-subfield census that came back
  empty at 131.  Note this excludes less than it first appears: normal-basis
  weight shells, which look like a third shape, are shift-invariant and so are
  orbit unions, already inside E6.  A base outside both would be a new idea, and
  this list is for experiments, not for ideas.
