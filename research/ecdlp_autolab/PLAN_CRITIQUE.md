# Adversarial critique of the v1-v11 plan portfolio

**Date:** 2026-05-28
**Purpose:** Self-critique of the research program portfolio. Each plan
has weaknesses; this document surfaces them so future iterations can
patch them.

The critique is structured as: **assumption** → **why it might be wrong**
→ **what to do about it**.

## Critique of v7 (engineering to close C-extension gap)

### Assumption 1: Montgomery gives 2× speedup at 81-bit
**Why might be wrong:** At 81-bit, GMP's `mpz_mod` is already
optimized via Barrett internally. Our hand-rolled Barrett might
match GMP, not beat it. Real speedup might be 1.0-1.3×.

**Mitigation:** Day 1-2 should explicitly compare hand-rolled to GMP
on identical workloads, not project a 2× improvement.

### Assumption 2: SIMD multi-target gives 1.5× speedup
**Why might be wrong:** Modular arithmetic has data-dependent branches
(carry propagation in Barrett, conditional subtractions). SIMD lanes
must mask through these, possibly losing speedup. Real may be 1.0-1.2×.

**Mitigation:** Prototype Phase 22.2 on a tight kernel first; measure
before integrating into full rho.

### Assumption 3: Combined stack gives 5.15× (closes the 5.5× gap)
**Why might be wrong:** Multiplicative speedup assumes independent
optimizations. In reality, Barrett + SIMD + larger partition share
cache and instruction bandwidth; real combined speedup may be 3-4×.

**Mitigation:** Plan for 3-4× and explicitly accept that 80-bit may
remain infeasible in budget. Phase 22.6 fallback: relax success criterion.

### Assumption 4: Hyperthreading helps
**Why might be wrong:** Modular arithmetic is memory-bound with high
ILP; HT might not help and could hurt due to cache contention.

**Mitigation:** Phase 22.4 is explicitly "measure, don't assume". Already
in plan.

## Critique of v8 (algorithmic frontiers)

### Assumption 1: F5 specialized for Semaev escapes Yokoyama bound
**Why might be wrong:** The Yokoyama bound is about ideal regularity,
which is invariant under choice of Buchberger vs F4/F5 (these all
compute the same GB). Specialization might give constant-factor wins
but not asymptotic.

**Mitigation:** Phase 23 should explicitly state: even if F5 wins
asymptotically, it might be by a factor 10-100×, not the orders of
magnitude needed to beat rho. Discuss in literature review.

### Assumption 2: Drinfeld lift preserves DLP structure
**Why might be wrong:** Lifting `E/F_p` to a Drinfeld module setting
requires a specific construction. The lift might not exist or might
make ECDLP HARDER, not easier.

**Mitigation:** Phase 24's first week should be pure literature survey
to confirm the lift exists for typical primes.

### Assumption 3: Tropical reformulation preserves DLP
**Why might be wrong:** Tropical max-plus algebra doesn't preserve
group structure. The "discrete log" in tropical sense might be
totally different from classical.

**Mitigation:** Phase 25's first day should clarify what tropical-ECDLP
even means.

### Assumption 4: Phase 28 (ML representation) will fail cleanly
**Why might be wrong:** Could actually succeed at small bits (12-16),
giving a misleading positive that doesn't scale. The NN might "memorize"
small instances.

**Mitigation:** Phase 28 should test scaling from 8 to 16 to 24 bits,
showing accuracy drops to chance.

## Critique of v9 (long-horizon)

### Assumption 1: BSD-conjecture-conditional results matter
**Why might be wrong:** BSD might never be proven; even if proven,
the algorithmic content for ECDLP is unclear. v9 Phase 36 might be
purely vacuous.

**Mitigation:** Phase 36 should produce a CONDITIONAL result statement
even before BSD is proven — make the conditional explicit and useful.

### Assumption 2: Cross-disciplinary methods will work
**Why might be wrong:** Cross-disciplinary attacks rarely succeed
without prior anchor. Without specific prior research, v9 Phase 38
is speculation without scaffolding.

**Mitigation:** Phase 38 should explicitly start with literature
review, NOT new mathematics.

### Assumption 3: Distributed cryptanalysis is just engineering
**Why might be wrong:** It IS engineering, but coordinated distributed
attacks have specific bottlenecks (network latency, fault tolerance,
key distribution) that aren't trivial. Phase 39 underestimates this.

**Mitigation:** Phase 39 should reference specific prior work (Wenger-
Wolfger 2014 113-bit Koblitz, Bos et al. 2012 112-bit) for honest cost
estimates.

## Critique of v10 (moonshots)

### Assumption 1: All M1-M15 are correctly characterized as failing
**Why might be wrong:** Some of these are "almost certainly wrong"
but not formally disproven. Someone might find a clever angle.

**Mitigation:** Each M-entry should have an explicit `Status: tested`
vs `Status: speculation` flag. Many M-entries are speculation that we
haven't tested.

### Assumption 2: Documenting bad ideas helps
**Why might be wrong:** A researcher seeing M1-M15 might assume the
list is exhaustive ("if it's not here, it's worth trying") when it's
just our best guess.

**Mitigation:** v10's preamble should say "this is NOT exhaustive;
just our catalog. New M-entries should be added by experimenters."

## Critique of v11 (publication)

### Assumption 1: J. Math. Cryptol. will accept this paper
**Why might be wrong:** Negative results are notoriously hard to publish.
Most reviewers want a positive contribution. Our paper is mostly
negative with one engineering positive.

**Mitigation:** Frame the paper as "empirical methodology + map"
rather than "negative results catalog". The methodology contribution
(50-trial protocol, end-to-end scaling validation) is itself novel.

### Assumption 2: Open-science release is straightforward
**Why might be wrong:** Reproducible Sage builds are tricky (PARI/GP
version, Singular library, GMP version all matter). Docker container
helps but bloats releases.

**Mitigation:** Plan Phase 42 with explicit version pinning AND
docker container. Test on at least 2 different machines.

### Assumption 3: Conference talks will be accepted
**Why might be wrong:** Workshop talks are selective; without a positive
algorithmic result, getting accepted is harder.

**Mitigation:** Lead with the C-extension engineering result (positive)
and the methodology (novel), not the negative attack catalog.

## Critique of PROGRAMS.md and PROGRAMS portfolio

### Assumption 1: The numbering v1-v11 is informative
**Why might be wrong:** Numbers suggest chronology, but v10/v11 are
not "newer" than v8/v9; they're different *types* of programs. The
numbering creates false ordering.

**Mitigation:** Could rename to v6-engineering, v6-frontier, v6-roadmap,
etc. But that's churn; live with the current numbering.

### Assumption 2: Sequential recommended order works
**Why might be wrong:** v11 (paper) and v7 (engineering) are in
parallel, but the plan implies v11 first. In practice they're both
~3 months and can interleave. Implementation may differ.

**Mitigation:** PROGRAMS.md notes "parallel" but could be more explicit.

### Assumption 3: We will execute these plans
**Why might be wrong:** The user might decide ECDLP cryptanalysis isn't
where they want to invest more time. v6's negative result might suffice.

**Mitigation:** Treat v7-v11 as "ready to execute if/when desired".
Don't expect a specific schedule.

## Critique of overall portfolio

### Bias 1: Optimism on engineering
We assume engineering can give 5.15×. History shows engineering
projections fall short. Realistically expect 2-3×.

### Bias 2: Pessimism on algorithmic
We assume no algorithmic breakthrough is achievable. History shows
breakthroughs do happen (FPPR 2012 for binary fields, Castryck-Decru
2022 for SIDH). Maintain humility that v8 might surprise.

### Bias 3: Limited adversarial perspective
We didn't consider: what if someone has ALREADY found an attack but
hasn't published? Quantum + classical hybrid + side channel = future
attacks we can't predict.

**Mitigation:** Stay tuned to the literature; the resistance map is a
snapshot, not a permanent truth.

### Bias 4: Sample size of one
All work is on 4 LMFDB curves at 6 primes. The "uniform resistance"
claim might be specific to these instances. Other random curves might
behave differently.

**Mitigation:** v9 Phase 35 (curve-generation pipeline) could test
attack-resistance on hundreds of random curves.

## What this critique recommends

1. **Don't promise 80-bit feasibility from v7 alone.** Set expectation
   that 80-bit might remain infeasible, and we ship the engineering
   improvement regardless.

2. **Start v11 (paper) and v7 (engineering) in parallel.** Don't
   sequence them.

3. **Reframe v8 phases** to be explicit about their "research project
   not session experiment" nature. Don't expect quick wins.

4. **Audit v10** to mark which entries are tested vs speculated.

5. **Stay humble about v6's "exhaustive map"** claim. It's exhaustive
   for known classical attacks; new classes can emerge.

## Reading this in the future

If you're reading this and wondering "should I update this critique?"
— yes! Self-critique should evolve as new evidence arrives. Add
sections, contradict prior analysis, document what you learn.
