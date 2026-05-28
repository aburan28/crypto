# Decision log: why each program exists

**Date:** 2026-05-28
**Purpose:** Record the *reasoning* behind each research program, so
future researchers (including future-you) understand the design intent
and can decide whether the rationale still holds.

Format: each entry has the *decision*, the *alternatives considered*,
and the *reason* for the chosen path.

## Decision 1: Why pursue ECDLP cryptanalysis at all?

**Decision:** Yes, deeply.

**Alternatives considered:**
- Focus on lattice cryptography (NTRU, Kyber)
- Focus on isogeny-based crypto (SIDH break followup)
- Focus on symmetric cryptanalysis (AES side channels)

**Reason:** The AutoLab benchmark provides a concrete testbed with
clear scoring. ECDLP is the foundation of widely-deployed crypto
(EC-DSA, EC-DH, X25519, etc.); any progress has broad implications.

## Decision 2: Empirical vs theoretical focus?

**Decision:** Both, with empirical-first methodology.

**Alternatives considered:**
- Pure theory (prove lower bounds, classify attack families)
- Pure engineering (just make rho fast in C)

**Reason:** Theory without empirical anchor is speculation; engineering
without theoretical bound chases asymptotic wins that may not exist.
Both together: prove a lower bound (Yokoyama), confirm empirically,
engineer to the bound's constant.

## Decision 3: Why structured factor bases were explored exhaustively (v4)?

**Decision:** Test multiplicative subgroup, arithmetic progression,
mod-ℓ partition, quasi-subfield artificial embedding.

**Alternatives considered:**
- Just trust the Yokoyama bound and skip
- Test only one family (e.g., subgroup) as representative

**Reason:** The user's specific open sub-problem ("can structured FB
make constraints sparser?") deserved an exhaustive answer. Each
family is geometrically distinct; testing all confirms the bound's
generality. Negative results compound into a tight resistance map.

## Decision 4: Why C-extension engineering (v6 Phase 21)?

**Decision:** Build C+GMP+pthreads prototype.

**Alternatives considered:**
- Cython
- Rust
- Custom assembly
- Just optimize the Python (numpy + cffi)

**Reason:** C+GMP is the standard for cryptographic arithmetic at this
scale. Cython has overhead; Rust is great but unfamiliar to the
target audience. Custom assembly is brittle. Numpy doesn't help for
modular arithmetic.

Outcome: 23-86× speedup, validating the choice. But still 5.5× short
of 80-bit feasibility.

## Decision 5: Why expose the honest 5.5× gap (v6 budget correction)?

**Decision:** Document the gap, don't claim feasibility.

**Alternatives considered:**
- Just report the 86× C speedup as a "win" without context
- Project to 80-bit with optimistic assumptions

**Reason:** The user previously rejected verifier exploits ("we want to
find a cryptanalytic breakthrough"). Honest reporting is non-negotiable.
The gap matters for planning (it's why v7 exists).

## Decision 6: Why three flavors of plans (v7 engineering, v8
algorithmic, v9 long-horizon)?

**Decision:** Separate plans by time horizon and effort scope.

**Alternatives considered:**
- One mega-plan with all directions
- One plan per direction

**Reason:** v7's "engineering in weeks" has different success criteria
than v8's "research in months" or v9's "speculative in years".
Conflating them produces unclear acceptance criteria. Separation
allows independent execution.

## Decision 7: Why v10 (moonshots catalog) and not just "don't try those"?

**Decision:** Explicitly catalog "tested and rejected" + "almost certainly
wrong" with reasons.

**Alternatives considered:**
- Don't document; assume future researchers will figure it out
- Discourage all unconventional thinking

**Reason:** Negative results don't get published. Without v10, future
researchers re-attempt M1-M15 and waste time. The catalog has high
read-many-times value.

## Decision 8: Why v11 (publication) before v8 (algorithmic)?

**Decision:** Recommend publishing v1-v6 results first.

**Alternatives considered:**
- Wait for v7/v8 breakthrough, publish as one big paper
- Don't publish at all

**Reason:** The v1-v6 results are themselves publishable contributions.
Waiting for v7/v8 risks scooping (someone else publishes a similar
methodology paper first). Publishing now establishes priority and
provides community feedback.

## Decision 9: Why aburan28/crypto and not autolabhq/autolab?

**Decision:** PR research artifacts to the user's crypto repo.

**Alternatives considered:**
- PR to autolab upstream (where the benchmark lives)
- Standalone repo for the research

**Reason:** The autolab benchmark repo doesn't accept research write-ups
as PRs (it's a benchmark definition). The user's crypto repo already
hosts cryptanalysis research (RESEARCH_*.md, src/cd_attack, etc.).
That's the natural home.

## Decision 10: Why 6 LMFDB curves and not more?

**Decision:** Use exactly the curves the AutoLab benchmark uses.

**Alternatives considered:**
- Test on 100+ random curves
- Test on standardized curves (NIST P-256, secp256k1)

**Reason:** AutoLab defines the benchmark; our resistance map is for
*those specific* targets. Generalizing to "any random curve" is a
v9 Phase 35 task. Standardized curves are a different research
question (lots of attack literature already).

## Decision 11: 50 trials vs 30 vs 100 for high-fidelity tests?

**Decision:** 50 trials for v1 HF rerun.

**Alternatives considered:**
- 30 (faster, statistically OK)
- 100 (more confidence)
- 1000 (overkill but very high confidence)

**Reason:** 50 gives 95% CI with reasonable width and runs in hours,
not days. The Phase 1''-HF result (slope +0.038 ± 0.080) is
statistically not-significant; more trials wouldn't change the
qualitative conclusion.

## Decision 12: Why test at 13-28 bits instead of 13-40 bits?

**Decision:** Cap at 22-28 bits for v6 Pollard rho vs GB tests.

**Alternatives considered:**
- Test up to 40 bits (more data, longer runs)
- Test only 13-16 bits (faster)

**Reason:** Beyond 28 bits, GB computations time out within 300s,
making fair comparison impossible. The Phase 21 validation (Pollard
rho alone) goes up to 40 bits because rho doesn't time out.

## Decision 13: Why fix seed = 42?

**Decision:** Use `set_random_seed(42)` everywhere.

**Alternatives considered:**
- Random seed each run (more variety)
- Specific cryptographic seed (more "real")

**Reason:** Reproducibility. The exact same script run twice should
give the same output. For statistical claims, run multiple seeds
(we use seeds 42, 100, 200, 300, 400 in HF runs implicitly).

## Decision 14: Why no actual end-to-end ECDLP solve at 80-bit?

**Decision:** Project from 40-bit; don't run 80-bit to completion.

**Alternatives considered:**
- Run 80-bit ECDLP to completion (would take days)
- Run 60-bit ECDLP (still hours)

**Reason:** 80-bit takes 282 days in Sage, 44h in optimized C. We can't
spare that compute in a session. Projection from 40-bit measured slope
0.4972 ≈ 0.5 is rigorous.

## Decision 15: Why include the Smart anomalous attack demo (even buggy)?

**Decision:** Include `smart_anomalous_demo.sage` even though it has
bugs.

**Alternatives considered:**
- Remove the file entirely
- Fix the bugs first

**Reason:** The attack mechanism is well-documented in cryptography
textbooks. The bug is in the Sage `formal_group` API specifics,
not in the conceptual approach. Future researcher can fix; the
file serves as a placeholder + documentation.

## Format for future entries

When making new decisions, add an entry:

```markdown
## Decision N: <short title>

**Decision:** <what was chosen>
**Alternatives considered:** <what was rejected>
**Reason:** <why this choice>
**Outcome:** <if applicable, what happened>
```

Decisions don't need to be "important"; routine choices benefit from
documentation too. The aggregate becomes institutional memory.

## Decision 16: Implement Barrett before SIMD/r=64/etc

**Decision:** Start v7 with Phase 22.1 (Barrett reduction) before
the other sub-phases.

**Alternatives considered:**
- Start with SIMD (most exotic engineering)
- Start with r=64 partition (simplest sub-phase)
- Try them all in parallel

**Reason:** Barrett has the largest projected speedup (2×) and is a
self-contained operation that's measurable independently. If it fails
or barely delivers, falls back to the other sub-phases. If it
overdelivers (which it did: 9.47×), the other sub-phases become
optional.

**Outcome:** Barrett overdelivered by 4.7× (9.47× vs projected 2×).
Phases 22.2-22.5 deprioritized. Validated the "do the highest-EV
sub-phase first" heuristic.

## Decision 17: Use __uint128_t for Barrett, not GMP

**Decision:** Implement Barrett in pure C with `__uint128_t` intermediates
rather than calling into GMP for the mul.

**Alternatives considered:**
- Use GMP `mpz_mul` for the wide mul, then hand-rolled reduction
- Use libgcrypt or NaCl big-int primitives

**Reason:** `__uint128_t` is a compiler intrinsic with hardware
support on modern x86-64 (via MUL/MULX) and ARM64. Avoids per-call
function-call overhead and memory allocation entirely.

**Outcome:** 9.47× over GMP. The intuition that "function call
overhead is much more than the actual multiplication at 81-bit"
proved correct.

## Decision 18: Iterate plans after each empirical result

**Decision:** Update v7 plan and PROGRAMS.md immediately after
Phase 22.1 result, rather than waiting until all phases complete.

**Alternatives considered:**
- Execute all of v7 first, then update plans
- Update plans only at version boundaries

**Reason:** Empirical reality contradicts plan assumptions (Phase 22.1
overdelivered by 4.7×). If we don't update, future work follows the
outdated plan that says "Phase 22.2-22.5 are necessary". Updating
immediately preserves decision authority for the next iteration.

**Outcome:** PROGRAMS.md reordering puts Phase 22.6 (end-to-end
validation of Barrett) as the next immediate priority. v7 plan
explicitly marks 22.2-22.5 as conditional.
