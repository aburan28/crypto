# Note, 2026-09-30: Frobenius-orbit coordinates for the decomposition step

This closes, by structure, the "Frobenius-orbit coordinates" item that the 2026-09-30
deep-research shortlist and the stop decision's reopening list carry as a symmetry lever.
It runs nothing. It cites what exists, states the one reduction that is not already written
down, and names what would reopen it. It changes no decision.

## 1. The item, in its three readings

Thread 1 of [KOBLITZ_SUBFIELD_INDEX_CALCULUS_QUEUE.md](../../KOBLITZ_SUBFIELD_INDEX_CALCULUS_QUEUE.md)
asks to "solve directly in Frobenius representative + relative-shift coordinates; measure
whether the decomposition problem itself shrinks." Survey
[§3.2](DECOMPOSITION-SURVEY.md#32-symmetry-enlarged-descents-2-4-torsion-translations-frobenius-orbit-coordinates-as-a-lever-on-the-degree-slope-of-31)
lists it beside the 2-torsion and symmetric-group levers. There are three things it can mean.

### 1a. A per-target symmetry of the decomposition system: does not exist

A symmetry lever rewrites the system for one target `R` in invariants of a group acting on
the summands while fixing `R`. Frobenius does not fix `R`:

```text
   P₁ + … + P_m = R      ⇒      σP₁ + … + σP_m = σR,   and  σR ≠ R  unless  R ∈ E(F₂).
```

`E(F₂)` is `{O, T}` on `K₀` and `{O, T, (1, ·), (1, ·)}` on `K₁`, none in the prime-order
subgroup. So Frobenius carries the system for `R` to the system for `σR`; it is a map
**between** targets, not an automorphism of one target's system, and there is no invariant
rewriting of a single system. This is the sentence the RR encoding's own README already
wrote: "Frobenius transforms `R` as well as the summands; do not quotient individual point
orbits while silently holding `R` fixed" (`research/nagao_relations/README.md`, "Next
experiments" 3). The torsion translation of §3.2 is different precisely because `T` is
rational: `P ↦ P + T` on summands is absorbed by an even number of translations.

### 1b. Orbit representatives on a Frobenius-stable base: the fold, already priced

If `V = ker g(σ)` is Frobenius-stable then `F` is a union of `σ`-orbits and every relation
gives `n` relations. That is survey
[§3.4](DECOMPOSITION-SURVEY.md#34-frobenius-stable-and-quasi-subfield-factor-bases-a-polyn-lever-and-an-availability-table-not-an-exponent):
a factor of `n` on the collection phase, against the `√n` the matched rho already takes
from the same Frobenius, so at most `√n` relative — a constant in `n`, never an exponent.
Its availability is exact
([decomp-frobenius-stable.txt](decomp-frobenius-stable.txt)): **no** stable subspace exists
at 13, 19, 29, 37, 53, 59, 61, 131 or 163, where 2 is a primitive root; 17, 23, 41, 43, 47
have them only at about `n/2`; only `n = 31` has one near the `m = 4` optimum. The
per-target solve is unchanged under this reading: with `V` stable, the system for `R` still
has `mℓ + (m−2)n` unknowns, and 1a applies to it.

### 1c. Orbit coordinates on a non-stable base: the Frobenius-conjugate coset base

The remaining reading, and the one thread 1 most plausibly intends at a degree with no stable
subspace: keep an arbitrary `V`, enlarge the base to `F' = ∪_k σ^k(F)`, and write a summand as
`σ^{k_i}(P_i)` with `P_i ∈ F` and a shift `k_i ∈ Z/n`. Then

```text
   x(σ^k P) = x(P)^{2^k} ∈ V^{2^k},
```

so for a fixed shift vector `(k₁, …, k_m)` the decomposition is over the `m` subspaces
`V^{2^{k₁}}, …, V^{2^{k_m}}`, each of dimension `ℓ`. That is a **coset-typed factor base**
in the sense of Galbraith–Gebregiyorgis 2014 (`inputs/GG-2014-806/`) and Nagao 2015/984 §7
(`EQS3`/`EQS4`), with the cosets `V + v_i` replaced by the conjugate subspaces `V^{2^k}`.
It is not a new object:

- **Yield.** `|F'| = n·|F|` up to overlaps. The yield per target rises by about `n^m/m!`
  (the distinct-shift classes) and the number of typed systems to solve rises by the same
  factor, so the yield per unit of solving is unchanged. A single random subspace of
  dimension `ℓ + log₂ n` gives the same yield with one system, so the ordinary `ℓ` trade
  dominates this one unless the typed systems are individually cheaper.
- **Per-system cost.** Whether a typed system is cheaper than the untyped one — the
  removal of the `m!` symmetry — is exactly `H-SEMBIN-c59e50`'s claim
  (`crypto-autoresearcher/ledger/hypotheses/H-SEMBIN-c59e50.yaml`), on `EQS4`, in the
  SEMBIN lane. Survey §3.3 already says: do not duplicate it.
- **The fold.** `F'` is stable by construction, so relations fold by `n`; that is 1b again,
  bought at `n×` the base size, and it is the `√n` constant of 1b.

So 1c is "the coset-typed base with a specific choice of cosets". Nothing about the
conjugacy of the cosets enters the per-target algebra except through the field equations,
which are the same for every coset.

## 2. What this note claims and does not claim

- **Claims.** Under every reading, Frobenius-orbit coordinates change the per-target
  decomposition system either not at all (1a, 1b) or into a coset-typed system whose cost
  question is owned elsewhere (1c). The constant they buy on collection is `≤ √n` relative to
  the matched rho (1b). None of this can move an exponent, and none of it is available as a
  subspace fold at the tournament's prime degrees or at 131 and 163.
- **Does not claim.** Anything about the SEMBIN lane's typed-system cost, which is open and
  is measured there; anything about `n = 31`, 233 or 571, where a stable subspace of a useful
  dimension exists and the `√n` constant is real; anything beyond `m ≤ 5`, which is where
  the availability table stops.

## 3. What would reopen it

- A **rational** automorphism of the curve beyond `±1` and the torsion translations. Survey
  §4: `Aut(E) = {±1}` and `End(E) = Z[τ]`; there is none.
- A decomposition step whose target is Frobenius-fixed by construction, that is, an
  index-calculus variant that decomposes elements of `E(F₂)` only. There is no such variant:
  the targets are the walk's points.
- A measured typed-system cost below the untyped one in the SEMBIN lane. That would make 1c
  worth building at a degree where the conjugate cosets are the cheapest typing available,
  and it would be SEMBIN's result carried here, not a Frobenius result.

## 4. Where it sits in the records

- Survey §3.2 keeps the item listed as a lever; this note is the missing "measured here /
  cheapest falsification" pair for it: no measurement is needed, and the cheapest
  falsification is the one line in 1a.
- The stop decision's reopening item "a symmetry not tested here: … Frobenius-orbit
  coordinates" is met for readings 1a and 1b by structure, and routed for 1c.
- The symmetric-group action, the other half of that reopening item, is measured beside the
  Riemann–Roch norm form in
  [ic_rr_norm_ladder_20260930](../../ic_rr_norm_ladder_20260930/PREREGISTRATION.md), whose
  `rr` arm is fully symmetric in the summands.
