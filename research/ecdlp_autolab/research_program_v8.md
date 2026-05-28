# Research Program v8: Algorithmic frontiers

**Date:** 2026-05-28
**Status:** planned, not yet executed
**Motivation:** v6 closed every known cryptanalytic direction; v7
attacks the C-extension constant-factor gap. v8 asks: **are there
algorithmic ideas that have never been tried for prime-field ECDLP?**

This is the most speculative program. Each phase is a research-paper-
sized investigation, not a session experiment. We design feasibility
tests and small-scale validations that could indicate whether a
direction is worth pursuing.

## Phase 23: Signature-based Groebner (F5/F4*) specialized for Semaev

**Hypothesis:** Yokoyama et al. 2020's lower bound is for *naive*
Buchberger. F5 with signatures is empirically faster (sometimes by
super-polynomial factors) on systems with non-trivial syzygies.
The Semaev F_n polynomial system has rich syzygy structure (the
elimination ideal of n-1 dependent variables).

**Concrete experiment:**
1. Implement F5 for the Semaev ideal `⟨F_3(xR, X_2, X_3), Q(X_2), Q(X_3)⟩`
2. Specialized signature management for the polynomial structure
3. Measure first-fall-degree (FFD) — does it match the regularity?
4. Compare per-query cost to F4 baseline

**Risk:** F5 implementations are notoriously delicate; even Sage's
F5 is incomplete. May require months of effort.

**Reference:** Faugère 2002 (original F5), Eder-Faugère 2017 survey,
Huang-Kosters-Yeo 2015 (LFD vs FFD).

## Phase 24: Drinfeld-module-based ECDLP variants

**Hypothesis:** Drinfeld modules are function-field analogs of
elliptic curves. ECDLP has counterpart "Drinfeld module DLP" with
*sub-exponential* algorithms (Gekeler 2008, van der Heiden 2004).
If we can lift our prime-field ECDLP to a Drinfeld module setting
(e.g., via Artin-Schreier extension), the sub-exponential machinery
might apply.

**Concrete experiment:**
1. Choose a base curve (toy size) and construct its Drinfeld lift
2. Apply Gekeler's algorithm; measure cost
3. Compare to direct ECDLP

**Risk:** The lift might not exist or might preserve hardness.
This is essentially a literature survey + small-scale prototype.

**Reference:** Drinfeld 1974, Gekeler 2008, van der Heiden 2004.

## Phase 25: Tropical / max-plus algebra reformulation

**Hypothesis:** Tropical geometry (`+, max` instead of `+, *`) gives
linear algebra over the tropical semiring. ECDLP, being multiplicative
in nature, becomes "linear" in the tropical setting. Vaccon-Verron-
Yokoyama 2018 (ISSAC) showed tropical F5 has different complexity
than classical F5.

**Concrete experiment:**
1. Reformulate the Semaev system tropically
2. Apply Vaccon et al.'s tropical F5
3. Measure complexity at small bits

**Risk:** Tropical reformulation may not preserve the discrete-log
structure; correctness uncertain.

**Reference:** Vaccon-Verron-Yokoyama 2018 (ISSAC).

## Phase 26: Hyperelliptic descent

**Hypothesis:** For some elliptic curves, there exists a 2-isogeny
to a Jacobian of a genus-2 hyperelliptic curve. ECDLP on the elliptic
descends to DLP in the Jacobian. Genus-2 DLP has more sophisticated
algorithms (Gaudry 2007, Diem 2006, large-prime variant).

**Concrete experiment:**
1. For our LMFDB curves, search for explicit Jacobian descent maps
2. If found, measure genus-2 DLP cost vs. genus-1 cost
3. Cross-reference with `RESEARCH_DIEM_DESCENT.md` (already in the repo)

**Risk:** Hyperelliptic descent may not exist for our specific curves;
when it does, algorithms are still subexponential, not polynomial.

**Reference:** Gaudry 2007 ("Index calculus for abelian varieties of
small dimension"), already cited in our repo.

## Phase 27: Brandt matrix / quaternion-algebra-based attacks

**Hypothesis:** Supersingular curves have Brandt matrix structure
on their isogeny class. Ordinary curves are not directly amenable,
but the Eichler-Selberg trace formula relates ordinary and
supersingular curves via L-function values.

**Concrete experiment:**
1. Compute Eichler-Selberg trace formula at small prime p
2. Look for relationships between ordinary ECDLP and supersingular
   isogeny problems
3. Survey: does any literature suggest a path?

**Risk:** Highly speculative; probably nothing concrete to test.

**Reference:** Eichler 1955, Selberg 1956.

## Phase 28: Machine-learning representation of ECDLP

**Hypothesis:** Could a neural network learn the discrete log map
`P → log_G(P)` for an elliptic curve, given many (P, log) training
pairs? This would not "break" ECDLP in the formal sense, but might
reveal structure not visible to traditional algebraic methods.

**Concrete experiment:**
1. Train a small NN on (P, log_G(P)) pairs from a 20-bit curve
2. Measure generalization error
3. Compare to baseline random predictor

**Risk:** Almost certainly fails — ECDLP is provably random-oracle-like.
But worth verifying explicitly that no shortcut exists.

**Reference:** No prior work I'm aware of (worth checking SAT-based
papers like `RESEARCH_SAT_SEMAEV.md`).

## Phase 29: p-adic Newton-Iwasawa methods

**Hypothesis:** Iwasawa theory describes p-adic structure of L-functions
and Selmer groups. Could the Iwasawa main conjecture give us an
ECDLP shortcut at primes where the relevant Iwasawa module has
special structure?

**Concrete experiment:**
1. Survey: what p-adic invariants of `E/Q` are computable mod `p^k`
   for cryptographic-sized `p`?
2. Test for special structure on benchmark targets

**Risk:** Highly speculative; Iwasawa theory deals with global
fields, not directly with ECDLP. May not give attack at all.

**Reference:** Iwasawa 1973, modern surveys (Kato 2004).

## Phase 30: Multi-target relation across curves

**Hypothesis:** The 6 LMFDB benchmark targets share structure
(same originator curves, similar primes). A relation that involves
ALL 6 targets simultaneously might exist where individual relations
don't.

**Concrete experiment:**
1. Define multi-curve relations: `Σ_i a_i P_i + b_i Q_i = 0` (in some
   common ambient group)
2. Determine if such relations are findable in less time than 6×
   individual ECDLPs
3. Probably negative, but a clean negative result is valuable

**Risk:** No obvious common ambient group exists for distinct curves.

**Reference:** Faz-Hernandez et al. 2015 (cross-curve attacks for
specific construction).

## Priority order

For limited research bandwidth:
1. **Phase 23 (F5 for Semaev)** — best literature anchor; highest impact
2. **Phase 26 (Hyperelliptic descent)** — concrete computable test;
   complements existing `RESEARCH_DIEM_DESCENT.md`
3. **Phase 28 (ML representation)** — quick to test; even negative result
   is valuable
4. **Phase 30 (Cross-curve relations)** — quick negative result expected
5. **Phase 25 (Tropical)** — interesting but theoretically uncertain
6. **Phase 24 (Drinfeld lift)** — months of work, uncertain payoff
7. **Phase 29 (Iwasawa)** — most speculative
8. **Phase 27 (Brandt/quaternion)** — quickly verifiable null

## Success criteria

A "breakthrough" in v8 would mean:
- Phase 23: F5 specialized gives sub-`O(√n)` empirical scaling
- Phase 26: explicit genus-2 descent map exists for an LMFDB target
- Phase 28: NN learns the DL map with non-trivial accuracy (would be
  publishable if it works)

A "negative result" in v8 closes another door, extending the
empirical resistance map.

## What v8 does NOT pursue

- Quantum algorithms (out of classical scope)
- Side-channel attacks (require implementation access)
- Pure number-theoretic conjectures with no algorithmic content

## Cross-references

- Yokoyama lower bound is for *naive* IC; F5 specialized (Phase 23) is the
  most direct candidate to escape it
- `RESEARCH_DIEM_DESCENT.md` (this repo) is the foundation for Phase 26
- `RESEARCH_HIGHER_SEMAEV.md` (this repo) relates to Phase 23
- The Yokoyama paper itself acknowledges F5-based variants may have
  different complexity profile

## Estimated effort

- Phase 23: 2-4 weeks (F5 implementation) + 2-4 weeks (benchmarking)
- Phase 26: 1 week (descent maps) + 1 week (genus-2 DLP)
- Phase 28: 3-5 days (NN training + analysis)
- Phase 30: 2-3 days (formal definition + small test)

These are paper-sized investigations, not session experiments.

## What comes after v8

If any phase delivers, write a full research paper.

If v8 closes uniformly, we've mapped the cryptanalytic landscape
exhaustively for the foreseeable future. The next progress would
require either fundamentally new mathematics or dramatically more
compute (large-scale distributed, quantum).
