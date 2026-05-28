# Research Program v10: Moonshots

**Date:** 2026-05-28
**Status:** speculative; documented for completeness
**Motivation:** v7-v9 cover concrete engineering, plausible algorithmic
frontiers, and a multi-year roadmap. v10 documents ideas that are
*almost certainly wrong* but worth recording so they don't get
re-discovered and re-rejected later.

This is a record of *eliminated possibilities* and *implausible-but-
not-ruled-out* directions. Each entry has a reason it probably fails.

## M1: Solve ECDLP by guessing

**Idea:** Random guess `k`, check if `k*P = Q`. Expected `n/2` guesses.
For 80-bit: 2^79 group ops.

**Why probably wrong:** Worse than Pollard rho's 2^40.

**Worth recording?** Yes, as the absolute upper bound (Shoup 1997).

## M2: Brute-force "from both sides"

**Idea:** Meet-in-the-middle BSGS: compute `i*P` for i ∈ [0, √n)
and `Q - j*P` for j ∈ [0, √n). Memory and time both √n.

**Why probably wrong:** This IS Shanks BSGS; same as rho asymptotically.
Already implemented in Sage's `discrete_log`. No improvement.

**Worth recording?** Just as the canonical √n baseline.

## M3: Pollard kangaroo with adaptive jumps

**Idea:** Instead of fixed-stride kangaroo jumps, adaptively choose
jumps based on observed coordinate values. Maybe better expected
cycle length.

**Why probably wrong:** Pollard's lambda has been analyzed
extensively; adaptive variants don't significantly improve the
asymptotic constant (van Oorschot-Wiener 1999).

**Worth recording?** As a "tested and rejected" entry in the empirical
map. Pollard rho is already optimal-up-to-constants.

## M4: Use the Riemann Hypothesis

**Idea:** RH gives bounds on prime gaps and L-function values. If RH
is true, maybe explicit attack constants improve.

**Why probably wrong:** RH gives statistical bounds, not algorithmic
shortcuts. ECDLP attacks aren't bottlenecked by prime distribution.

**Worth recording?** As an example of "true mathematical theorem
without algorithmic consequence".

## M5: Encode ECDLP in 3-SAT and solve via modern SAT solvers

**Idea:** Modern SAT solvers (CaDiCaL, Kissat) routinely solve
millions-of-variables SAT instances. Encode ECDLP as a SAT problem
and solve.

**Why probably wrong:** ECDLP encoded in SAT has ~n variables and
~n^2 clauses. For 80-bit n = 2^80, that's 2^80 variables — way
beyond SAT solver capacity.

**Cross-reference:** `RESEARCH_SAT_SEMAEV.md` in this repo (already
exists). Likely confirms this direction is closed.

**Worth recording?** Yes, as the cleanest negative for the
"complexity-theoretic reduction" family of approaches.

## M6: Encode ECDLP in QUBO and solve via D-Wave

**Idea:** Convert to Quadratic Unconstrained Binary Optimization
and solve on a quantum annealer.

**Why probably wrong:** Same scaling issue as M5. Annealers handle
~10^4 variable problems, not 10^24.

**Worth recording?** Yes, parallels M5.

## M7: Train a neural network on ECDLP

**Idea:** Train a transformer on `(P, k*P)` → `k` pairs.

**Why probably wrong:** ECDLP is provably pseudorandom-permutation-
hard under standard assumptions; NN can't learn what's structurally
random.

**Tested?** Phase 28 in v8 plans this experiment to formally rule out.

**Worth recording?** Yes, as a popular naive proposal.

## M8: Apply Galois Theory directly

**Idea:** The Galois group `Gal(K̄ / K)` acts on `E(K̄)`. Maybe the
Galois action recovers `k`.

**Why probably wrong:** The Galois action on `E(F_p)` is trivial
(F_p is fixed). The interesting Galois action is over `Q(E[ℓ])`,
which is for the ℓ-torsion. We've already studied this (Phase 4,
Galois rep) and found no exploit.

**Worth recording?** Cross-reference Phase 4.

## M9: Use the Birch and Swinnerton-Dyer conjecture

**Idea:** BSD relates `L(E, 1)` to the rank and regulator. For
rank-0 curves, `L(E, 1) ≠ 0` is provable.

**Why probably wrong:** BSD tells us about the global curve over
`Q`, not the reduction mod `p`. The link to ECDLP is unclear.

**Cross-reference:** Phase 19.2 (p-adic L-function) tried this
angle and found no connection.

**Worth recording?** Yes, as the "deep math, no algorithmic
consequence" archetype.

## M10: Use Iwasawa theory

**Idea:** Iwasawa main conjecture relates `p`-adic L-functions
to Selmer groups. Maybe Selmer structure encodes ECDLP.

**Why probably wrong:** Selmer groups are global invariants of
`E/Q`; they don't tell you about `E(F_p)` ECDLP directly.

**Worth recording?** Yes, parallels M9. v9 Phase 29 plans a deeper
survey.

## M11: Topological methods

**Idea:** Treat `E(F_p)` as a topological space (it's finite, so this
is silly, but maybe use algebraic topology of `E(C)` and reduce).

**Why probably wrong:** Algebraic topology of `E(C) ≈ S^1 × S^1`
gives `π_1 = Z^2`. This is the fundamental group, not the discrete
log structure.

**Worth recording?** Yes, as an example of conceptual category error.

## M12: Modular forms direct attack

**Idea:** Modularity theorem (Wiles, Breuil-Conrad-Diamond-Taylor)
relates `E/Q` to a modular form `f`. Maybe Fourier coefficients of
`f` encode ECDLP info.

**Why probably wrong:** Fourier coefficients give `a_p = p + 1 -
|E(F_p)|`. We already use this (it's just the trace). No further
algorithmic info known.

**Cross-reference:** Phase 2 (`phase2_modular_form_survey.sage`)
already explored.

**Worth recording?** Yes.

## M13: Compute `Q ÷ P` symbolically

**Idea:** Just symbolic division! `Q = k*P` means `k = Q/P`.

**Why probably wrong:** Group inversion is the EC group operation
inverse. We can compute `-P`, but `Q/P` ≠ `Q + (-P)` (that's
subtraction in the group). There's no division operation; ECDLP
*is* the question of inverting scalar multiplication.

**Worth recording?** Yes, as an example of conflating group operations
with field operations.

## M14: Lookup tables for huge groups

**Idea:** Precompute the entire `E(F_p)` group as a table.

**Why probably wrong:** Storage 2^n elements, n=80 → 2^80 entries.
Beyond all available storage in the universe.

**Worth recording?** Yes, as the upper-bound space argument.

## M15: "Just try harder"

**Idea:** Pollard rho with 10x more constants, 10x more cores, ...

**Why probably wrong:** Doesn't change asymptotics. We've already
done this analysis (v6 Phase 21.1).

**Worth recording?** Yes, as a reminder that engineering has limits
without algorithmic novelty.

## How to use this document

When someone proposes a "moonshot" attack, check if it's documented
here. If yes, the entry has the reason it (probably) fails. If not,
add it after testing or analysis. Over time this becomes a
comprehensive catalog of negative results.

## Companion to other programs

- Many M1-M15 entries cross-reference v8 phases (formal investigation)
  or specific test scripts in `ecdlp_autolab/`
- v9 Phase 40 (publication) can cite this catalog for "what doesn't
  work" appendix

## Reminder

Documenting "what doesn't work" is as important as documenting "what
works". The literature is full of unrecorded null results that cause
re-investigation by every new researcher. v10 prevents that for
prime-field ECDLP.
