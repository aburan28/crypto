---
name: rootfinding-autolab
description: Iterate on exact finite-field root-finding research with hypotheses, proofs, counterexamples, immutable bounded experiments, matched Gröbner baselines, and an active next-work queue. Use when developing or evaluating this project's algorithms, not for claiming an unmeasured breakthrough.
---

# Rootfinding AutoLab

Use this repository-scoped skill to keep the research program moving while
preserving its correctness and performance requirements. A failed candidate
ends that experiment; it supplies evidence for the next hypothesis.

## Start from evidence

Read `lab/METHODOLOGY.md`, the current program manifest and append-only ledger.
Consult `research/advance/ADVANCE.md`, `performance/GATE.md`, and
`research/advance_next/PREFIX_LIFT_NEXT_ROUTE.md` before extending their claims.
Use [the research queue](references/RESEARCH_QUEUE.md) to choose the next
highest-information task. Preserve old code, frozen manifests, raw runs, and
negative results; create new versioned variants and experiment IDs.

There are three independent assessments: **mathematical correctness and
completeness**, **novelty against primary literature**, and **measured effect
under qualified conditions**. A promising result in one does not imply another.
An LLM's novelty label and a smaller representation establish neither novelty
nor a speed gain.

## The iteration loop

1. **Hypothesis.** State the specific bottleneck, proposed mechanism, supported
   input family, output contract, and a falsifiable prediction. Identify the
   competing strong specialized method as well as Gröbner baselines.
2. **Derivation.** Prove necessary/sufficient conditions, coverage of every
   branch, multiplicities, exceptional fibers, and arithmetic assumptions.
   Give work, storage, and output bounds. Quantify surviving dimensions instead
   of assuming generic rank; count every conditional outer parameter.
3. **Counterexample/oracle.** Try degenerate curves, repeated roots, empty/full
   domains, torsion, zero coefficients, characteristic three, even binary
   extensions, and nonrational lifts. Compare independent exhaustive small
   oracles. A zero residual certifies membership, not completeness.
4. **Variant.** Implement the smallest component that tests the mechanism in
   a fresh directory. Run correctness checks before timing. An agent designs
   and reviews patches; the coordinator only runs declared commands.
5. **Preregistration.** Freeze source/input hashes, field basis, target rule,
   train/holdout split, command arrays, dependencies, budgets, seeds, and metrics
   before evaluating the holdout or comparison. A corrected implementation is
   a new variant and plan, never a silent edit to frozen evidence.
   Independent reviewers must use a unique author/version directory and check
   that it does not already exist before writing. A filename that looks like
   an appropriate review destination may already be registered. Never replace
   another review, receipt or log. If an accidental overwrite occurs, preserve
   the incident, restore exact registered bytes from the original author, and
   verify their hashes before continuing the checkpoint.
6. **Bounded experiment.** Validate hashes and run with process-group deadlines,
   RAM/work/state caps and one-core limits. Retain complete stdout/stderr and
   status. Do not co-run compilation or other experiments during performance
   isolation. Charge parsing, setup, extraction, verification and serialization.
   After building the complete result, serialize it once before accepting any
   exact/global completion flag. Demote late results to inconclusive and clear
   both original-source and auxiliary exact counts; retain verified tuples as
   lower bounds. The supervising process must additionally charge startup,
   final stdout and full child lifetime. A native or pre-evidence deadline
   check alone does not cover final reporting. Review this boundary before
   freezing every new worker or wrapper.
7. **Review.** Identify what was proved, disproved, estimated or unobserved.
   Keep timing/correctness evidence separate. Diagnose the failure mechanism;
   propose a new measurable consequence rather than merely raising a timeout.
8. **Next queue.** Append a decision and at least one concrete next action:
   repair an assumption, reduce width, derive a branch, build a missing baseline,
   improve representation, or seek prior art. Repeat within the session budget.

## Exact output and field scope

The main contract is every ordered original polynomial root over the declared
finite field and domain, with canonical coordinates, no duplicates, and an
exact completeness certificate. An exact count with truncated tuples is a
different contract. Unfinished search returns an unresolved count/lower bound;
it must never turn into an empty root set or an UNSAT claim.

Rational curve-point decompositions and unrestricted polynomial roots differ.
Keep nonrational-lift/ghost branches, intermediate coordinates, repeated sources,
and torsion unless a coverage proof permits their removal. Record rational-map
pole semantics. Separate characteristic-two and odd-characteristic proofs;
exercise prime and extension fields, including odd extensions and characteristic
three. Field arithmetic coverage does not establish algorithm or timing coverage.

The requested targets include the actual GF(2^83) family and the standard
NIST P-192, P-224, P-256, P-384 and P-521 prime fields and curves. Use the
source-identified presets and exact prime/group checks in
`lab/research/nist_fields.py`. Prime-field intervals and implicit domain maps
need separate algorithms: a nontrivial proper Fp-subspace does not exist in Fp.
Neither a field's bit length nor a small-domain calibration establishes hardness.

## Budgets and progression

Development defaults: 900 seconds per coordinator session, at most three ready
experiments, 4 GiB per worker, one core, and no network-dependent worker. Start
with a 90-second component experiment. Declare tighter algorithmic caps and
preflight large domains before materialization. These are configurable research
budgets, not release criteria or predictions of tractability.

Escalate from exhaustive small oracles to structural profiles, fixed calibration
ladder, unopened holdout, then matched hard comparisons. Reserve hard inputs
before seeing comparative timings. Calibration discoveries can motivate new
training variants, but they cannot select winners for the same holdout. Retire
and replace an exposed holdout under a recorded selection rule.

Initialize the declared manifest, then validate and run a bounded session:

```sh
python -m lab.runner init --manifest lab/research_program.json --ledger lab/ledger
python -m lab.runner validate --ledger lab/ledger
python -m lab.runner run-ready --ledger lab/ledger --max-experiments 3 --wall-budget 900
python -m lab.runner status --ledger lab/ledger
```

The finite-cycle `daemon` executes eligible plans; it does not invent or
patch algorithms. A worker returns exactly six `scientific_result_v1` fields:
`schema`, `experiment_id`, `status`, `summary`, `evidence`, and `next_actions`.
Place release-gate metadata inside `evidence`; extra root fields are rejected.
Retain refuted/inconclusive/blocked
attempts. A session budget pause or empty ready queue means **awaiting design**;
the program stays active unless the user ends it or its goal is satisfied.

## Frozen release gate

Never weaken `performance/GATE.md`: at least 20 distinct preregistered coupled
nonlinear hard inputs, fastest completed matched baseline median at least one
second, genuine source-identified F4/F5/F6 with the same all-roots contract,
every registered configuration complete, and independently identical root sets.
Require at least three cold repetitions, at least 2× geometric-mean gain against
the fastest baseline per input, no input slower than 1.10×, each A/A pair spread
at most 10%, and uncontended single-thread isolation. The actual GF(2^83) target
and separate requested field-family evidence remain part of the scope.

Timeouts, unsupported encoders, missing F6, truncated output and noisy hosts
make comparison evidence incomplete; they never become speed ratios. Missing
baselines create an independent engineering workstream. They do not stop proof,
oracle, circuit, conditional-splitting or prime-field experiments.

## Communicate progress

Report the current hypothesis, evidence gained, resulting bottleneck, and next
bounded experiment. Distinguish an experiment blocked on a dependency from a
research program that still has independent tasks. Work within existing user
authorization; this skill adds no approval step for local reversible research.
Do not publish, push, message other people, schedule hosted jobs, or install an
external AutoLab service merely because this local methodology uses that name.
