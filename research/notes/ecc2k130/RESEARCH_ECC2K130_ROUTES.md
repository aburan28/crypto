# Routes: where to take ECC2K-130 index calculus next

Ranked. Each route says why it is open, what to run, the one number that
decides it, and what would close it. Companion to
[`research/notes/ecc2k130/RESEARCH_ECC2K130_IC_LITERATURE.md`](RESEARCH_ECC2K130_IC_LITERATURE.md)
(what the literature has and has not) and
[`research/notes/ecc2k130/RESEARCH_KOBLITZ_SCALING_TARGET.md`](RESEARCH_KOBLITZ_SCALING_TARGET.md)
(the measured ladder these routes feed).

## The frame

Three established facts set the ranking, and all three are measurements
or citations rather than expectations:

1. **Constants do not move the exponent.** Relation collection with any
   *enumeration* oracle costs `2^l · 3!·2^{n−3l} · Θ(2^{2l}) = Θ(2^n)`,
   independent of `l` — measured flat within 1.5× across `l = 6…10`
   ([`research/notes/index-calculus/RESEARCH_SEMAEV_DECOMPOSITION.md`](../index-calculus/RESEARCH_SEMAEV_DECOMPOSITION.md)).
   A faster oracle buys reach, not an exponent.
2. **Nothing here is an attack.** The algebraic line's best *conditional*
   estimate is `2^86` at `n = 131` against rho's `2^60.8` — about `2^25`
   short, and that is granting the suspect first fall degree assumption
   *and* the best proven bound. These routes are about finding where the
   crossover actually sits, not about threatening the curve.
3. **Our measurements are the state of the art, not a reproduction of
   it.** No rigorous bound on `D_reg` for these systems exists in either
   direction; the only rigorous Weil-restriction bound is numerically
   vacuous at `n = 131` (solving degree ≤ 263). So the FFD ladder is
   contributing evidence, not re-deriving known results.

---

## Route 1 — Crossbred on the `m = 3` systems, and then into the attack

**Why it is open.** The survey found **no publication** combining
Crossbred (or any modern MQ/hybrid solver) with binary ECDLP
decomposition systems — item 3 of the research brief returned zero
claims, and the structural constraints that kill much of the rest
(`131` prime, `2` primitive mod `131`) **do not touch it at all**,
because the choice of solver is orthogonal to the factor-base
construction.

And it is half-built. `src/cryptanalysis/crossbred.rs` implements
Joux–Vitse Crossbred and is correctness-gated against matrix-F4 at
`m = 2` and `m = 3`
(`crossbred_covers_the_f4_oracle_on_the_cubic_m3_system`, 27 unknowns,
`D = 3`, `k = 12`), and `examples/crossbred_bench.rs` already prints
`xb/brute` and `xb/F4` columns over a ladder. What is missing is that
**no research note records what that bench found**, and Crossbred is not
a `DecompositionStrategy` — so it has never reached relation collection.

**Run.** (a) Write the existing bench's ladder up as a note, with the
`(D, k)` sweep and where a crossbred space exists at all (`kernel_dim`).
(b) Add `DecompositionStrategy::Crossbred` and re-run the scaling ladder
end to end.

**Metric.** `xb/F4` wall-clock ratio at `m = 3` on the largest rung both
finish; then relations/second end to end.

**Falsifier.** `xb/F4 ≥ 1` across the ladder with no `(D, k)` choice
doing better → the route is closed and the note says so.

**Run on 2026-09-20; the route stays open, narrowly.**
[`RESEARCH_ECC2K130_CROSSBRED.md`](RESEARCH_ECC2K130_CROSSBRED.md) does both
halves. (a) The ladder is written up: `xb/F4 = 0.023` in bit operations at
`m = 3`, `0.696` in the wall clock this route asked for, correctness gate
`yes` on every row, so the falsifier does not fire. (b)
`DecompositionStrategy::Crossbred` exists and reaches relation collection.
Three findings temper it. The wall-clock margin **shrinks with size**
(`0.040 → 0.696` over one rung of `n` at `m = 3`), so one more rung decides
the route. **No `(D, k)` produced a single filter**, which is the structural
reason: without filters the method is Macaulay preprocessing in front of `2^k`
linear solves, and the preprocessing is paid per target. And end to end
Crossbred is the best algebraic arm at `m = 2` but the **worst** at `m = 3`,
where exhaustive enumeration beats every algebraic oracle by `1800×`. The GPU
bonus below survives, and is worth nothing while the margin is heading for
`1.0`.

**Bonus, and the reason this ranks first.** Crossbred's search phase is
`2^k` independent points with no shared state and a bitwise AND against
a precomputed table per point — by the module's own account, the one
part of this problem family that maps onto a GPU without an argument.
The repository already has CUDA work for ECC2K-130 rho to borrow from.

## Route 2 — Crossbred on the *symmetrised* system

**Why it is open.** A composition nobody has tried, on either side. The
symmetrised `u`-frame system is both smaller and lower-degree than the
chained `x`-system — 13 unknowns at Boolean degree 4 against 30 at
degree 3 for `n = 15`; 25 against 44 at `n = 17`
([`research/notes/index-calculus/RESEARCH_EXOTIC_COORDINATES.md`](../index-calculus/RESEARCH_EXOTIC_COORDINATES.md) §8.2)
— and Crossbred's trade (reach degree `D` by algebra, finish by
enumerating `k` variables) is most favourable exactly where the degree is
low and the variables few.

**Prerequisite is already met.** `build_symmetrised_system` returns the
same `F2BoolPoly` shape `extract_crossbred` consumes, so this is wiring,
not a port.

**Metric / falsifier.** As Route 1. Closed early if no crossbred space
exists at any `(D, k)` for the symmetrised systems (`kernel_dim = 0`
throughout) — which is itself worth recording, since it would say the
symmetrised systems are *too* small for the technique.

## Route 3 — Wire the symmetrised oracle into end-to-end collection

> **Closed at its gate; do not wire.**  On ECC2K-130's structure the only
> Frobenius-stable `V ∋ 1` are `F₂` and the field, so the symmetrised base
> loses its orbit collapse.  Per relation, the oracle then costs `3.1–14.7×`
> the enumeration it would replace (`n = 13–23`, priced from below).  See
> §X4′ of [`RESEARCH_ECC2K130_ROUTE_TARGETS.md`](RESEARCH_ECC2K130_ROUTE_TARGETS.md).

**Why it is open.** The symmetrised oracle is measured at ~350× the
chained `x`-system at `m = 3` (10 291 ms → 29 ms, `K₀/F₂¹⁵`) and **none
of it reaches the attack**: `DecompositionStrategy` has four variants and
none is symmetrised. This is §8.5 item 2 of the exotic-coordinates note,
still open.

**This PR adds the missing prerequisite.** `SymmetrisedFactorBase`
carried the bases the *system* needs but none of the orbit structure the
*attack* needs; `frobenius_view_of_symmetrised` supplies it, gated on
four instances.

**Run.** `DecompositionStrategy::Symmetrised` dispatching through the
bridge, then the scaling ladder.

**Metric.** End-to-end relations/second and the measured orbit collapse,
against the `x`-chained arm.

**Falsifier.** The 350× does not survive contact with relation
collection — eaten by the `via T` lifting, by the different found/refuted
mix, or by the orbit structure.

**Risk to watch.** `F_u` and `F_x` differ as sets, so only within-verdict
medians may be compared; the found/refuted counts must be printed beside
any ratio.

## Route 4 — Reach `m = 4`

**Why it is open.** The conditional theory's optimal summand count is
`m ≈ n^{1/3} ≈ 5.1` at `n = 131`; the harness reaches `m = 3` and
strains at `m = 4`. That gap — between the `m` we can run and the `m` the
theory needs — is the real practical barrier, and the exotic-coordinates
note names the missing arm itself (§8.5 item 1): chain the *symmetrised*
`S₃` for `m ≥ 4`, giving `4(ℓ−1) + 1 + n` unknowns with bilinear links.
Its §8.4 also concedes that the fairer production baseline — chaining the
symmetrised rather than the plain `S₃` — "does not exist yet".

**Metric.** The scaling target's primary metric, `m·ℓ + (m−2)·n` at
`m = ⌈n/ℓ⌉`, brought under 64 with a solve and a clean gate.

**Falsifier.** The first fall degree grows with `n` on the `m = 4`
systems. That would close the route *and* be a result in its own right —
it is H1 of the scaling target, and per the frame above, nobody has a
rigorous answer.

## Route 5 — A second literature pass on the under-searched items

**Why it is open.** Research items 4 and 5 returned **zero surviving
claims**, flagged honestly as absence of evidence in the collected corpus
rather than evidence of absence. The structural constraints pre-filter
much of item 4, but they do **not** obviously touch item 5 at all:
whether Koblitz structure beyond Frobenius — the τ-adic expansion, CM by
`(1 ± √−7)/2`, the class group of `Z[τ]` — has ever been used for *index
calculus* rather than for rho or scalar multiplication.

**Run.** A targeted second survey with different angles: CM-specific and
Jacobian-specific queries, theses and non-English venues, and the
citation graph around Koblitz-curve CM rather than around ECDLP.

**Metric.** Any source giving a concrete operation count at prime
extension degree.

**Falsifier.** A second independent pass also returns empty — at which
point item 5 should be treated as genuinely unexplored, and the decision
becomes whether to explore it ourselves rather than whether to read more.

---

## Closed — do not reopen

Recorded so the effort is not spent twice.

- **Joux–Vitse cover-and-decomposition.** Structurally inapplicable at
  prime extension degree, not merely inefficient: it needs a tower
  `F_{q^d}/F_q/F_p` with composite degree, and `131` is prime. The
  obvious evasion was tested — embedding into `F_2^{131m}` collapses both
  towers, with `d = 131` giving a GHS cover of genus `≈ 2^130`, exactly
  the constraint
  [`research/notes/ecc2k130/RESEARCH_ECC2K130_HYPERELLIPTIC.md`](RESEARCH_ECC2K130_HYPERELLIPTIC.md)
  derived independently. No 2011–2026 work removes the condition.
- **Yokoyama et al. (JMC 2020) as a shield.** It is not a lower bound
  that closes the question: conditional on unproven assumptions plus an
  unproved conjecture, scoped to Semaev's *naive* method, prime fields
  only, never instantiated at a cryptographic size — and its §5.4 carves
  out exactly Diem subspaces and quasi-subfield polynomials, both already
  ours.
- **Pairs-and-solve as an oracle speed-up for the Koblitz `m = 3`
  path.** It is `|F|²·ℓ`, against `Enumerate`'s `|F|²` (the last summand
  is already a hash lookup) and `PairTable`'s `|F|` probes. It wins only
  against a naive triple loop, which the Koblitz code never used. Its
  value is elsewhere — deciding a target without materialising the factor
  base at all, which is what
  [`research/notes/index-calculus/RESEARCH_SEMAEV_DECOMPOSITION.md`](../index-calculus/RESEARCH_SEMAEV_DECOMPOSITION.md)
  measures.
