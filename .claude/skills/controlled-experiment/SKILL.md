---
name: controlled-experiment
description: Design, preregister, run and report a controlled experiment before making any empirical or comparative claim about curves, primes, groups, algorithms or their structure ("X is unusual", "X has no effect", "X is as random as", "method A beats B"). Requires classifying each claim, a negative-control population drawn from the right null distribution, a positive control proving the detector can see the effect, the same instrument on every arm, a frozen protocol before data, and reporting by rank/effect size with uncertainty. Use whenever a research note, PR or answer would state something an exact computation does not by itself establish.
---

# Controlled experiment

Read `AGENTS.md` and the `report-evidence` skill first. This skill decides
*whether a claim needs an experiment*, and if so, *what makes the experiment
able to fail*. `ecbench-measure` remains the harness for ECDLP cost
comparisons; use this skill to design what goes into it, and use this skill
alone for structural or statistical claims ecbench does not cover.

## 1. Classify every claim before writing it

Put each sentence you intend to publish in exactly one class, and label it
in the note.

| class | what establishes it | needs a control? |
|:--|:--|:--|
| **Theorem / identity** | proof, or a certificate (`identity.certificate/v1`) | no |
| **Exact computation** | a deterministic, verified computation (a factorisation checked by product and primality, an exact order) | no, but it must be reproducible and verified |
| **Comparative / statistical** | "unusual", "typical", "no effect", "as random as", "faster", "smaller", "independent" | **yes — this skill** |
| **Heuristic estimate** | a model (Dickman ρ, a square-probability argument, a cost formula) | it is a prediction; state it as one, and test it if the conclusion rests on it |
| **Out of scope** | a security or impossibility statement no experiment here can reach | say "not established" |

A comparative claim made from an exact computation alone ("P-256's n − 1 is
unusually smooth" from one factorisation) is not established. One number has
no reference. Do not write it until §2–§5 are done.

## 2. Build the control arms

Every comparative claim needs both.

- **Negative control (the null population).** A sample from the distribution
  the claim says the target is *not* different from, or *is* different
  from. Draw it from the closest valid null, and write down why it is the
  right one. For example: group orders of curves over a fixed `p` are
  well modelled by random primes in the Hasse interval `[p+1−2√p, p+1+2√p]`.
  Random 256-bit integers are a worse null, because they are not prime and
  not in the interval. Match everything except the property under test:
  size, primality, field and curve shape.
- **Positive control (detector power).** An input where the effect is known to
  be present, run through the *same* code path. A "found nothing" result is
  only evidence if the detector demonstrably finds the thing when it is
  there. Plant it: a curve with a constructed small point, a prime built to
  be smooth, an instance with a known answer.
- **Same instrument on every arm.** Same code, same budgets, same seeds policy,
  same metric, same stopping rule. If the instrument is incomplete (a
  factoring budget that leaves cofactors), that is acceptable *only* if it
  is identical across arms and the metric is defined on what it returns.
  Never compare a PARI factorisation of the target against a budgeted rho
  of the controls.

## 3. Freeze the protocol before running

Write `PROTOCOL.md` in the study directory and commit it (or at least hash it
into the results) before the first data point. It states:

1. the question and the exact claim each outcome would license;
2. the target(s), the null population and how it is sampled (generator,
   seed, sample size `N` and why `N` is enough for the resolution you need);
3. the positive control and the expected detection;
4. the metric(s), defined precisely, and the instrument with its budgets;
5. the decision rule, for example: "the target is *unusual* if its empirical
   rank is below the 2.5th or above the 97.5th percentile of the control,
   two-sided, with a Clopper–Pearson interval on the tail fraction";
6. how many metrics and targets are tested, and the multiple-comparison
   correction (Bonferroni over the family is fine);
7. what does *not* count: changing `N`, metric, null or threshold after
   seeing results. Any such change is a new, labelled protocol version, and
   both versions are reported.

**Power check before freezing (do the arithmetic).** An empirical rank test
cannot reject below its resolution. With the add-one estimate the smallest
attainable two-sided p is `2/(N + 1)`, so `N` must satisfy
`2/(N + 1) < α/m` for `m` corrected tests, with margin: `N + 1 ≥ 4m/α`
lets a target beyond every control reach `p = α/(2m)`. Check this for every
test family in the protocol and write the numbers into it. A protocol that
fails this check can only ever report "within the control distribution", so
its nulls are vacuous. (The first protocol written with this skill,
`research/nist_cyclotomic_controls_20261011/PROTOCOL.md`, failed it in all
three comparative families. The fix was a dated amendment, not an edit.)
**Pre-flight runtime estimate (write it into the protocol).** Before
freezing, time one instrument call per arm and per size class, and
separate the *per-sample* cost (scales with `N`) from the *fixed* cost
(target and positive controls, paid once per curve or target). Then
estimate `fixed + N · per_sample / threads` for every family. Fixed costs
are easy to miss. In the NIST re-test the dependent-pair positive control
(heights of `[3]G̃` and `[4]G̃`, about 16 times the generator's height and
super-linear in it) took longer than 100 random pairs did. Choose `N`
from the power check first. If the estimate is infeasible, change the
design (a cheaper positive control of equal strength, fewer families),
never the power. Record the estimate next to the actual run time, so the
next protocol calibrates against it.

Also do a dry run of the full analysis at small `N` on scratch output. It
catches power failures, runtime blow-ups and broken positive controls
before the real run. Disclose in the amendment anything that dry run
revealed about the target.

## 4. Run

- Native code per `AGENTS.md` (Rust), with a fixed seed recorded in the output.
- Emit raw per-sample rows (JSON lines or CSV), not just summaries, so the
  rank and intervals can be recomputed.
- Run the positive control first. If it fails, stop and fix the instrument.
  A negative result from a blind detector is reported as *invalid*, not as
  "no effect".

## 5. Report

Lead with the table. For each target × metric:

| target | metric | target value | control median [2.5%, 97.5%] | target rank / N | tail p (two-sided) | positive control detected? | verdict |

- **Verdict wording is fixed by the decision rule.** Use "within the control
  distribution (rank r/N)", "more extreme than q% of controls", "positive
  control failed — inconclusive". Never write "random", "as large as
  possible", "no effect" or "never" unless the claim class allows it.
- **Effect size before significance.** Give the difference in the metric's
  units next to the p-value.
- **Scope sentence.** State what the experiment does *not* show. "Within the
  null on these metrics" is not "secure", and "no small point below height H in
  this census" is not "no small point exists".
- Keep every earlier claim that this experiment weakens visible, with a
  dated correction: the original sentence, what was wrong, and the replacement.

## 6. Checklist (paste into the PR)

- [ ] every published sentence is classified (§1)
- [ ] comparative claims have a negative-control population with a stated null
- [ ] every "nothing found" has a passing positive control on the same path
- [ ] one instrument, identical budgets across arms
- [ ] power check passes: `2/(N + 1) < α/m` for every test family (§3)
- [ ] pre-flight runtime estimate (fixed and per-sample) in the protocol (§3)
- [ ] `PROTOCOL.md` frozen (hash recorded) before results; dry run on scratch output
- [ ] raw rows committed; summary recomputable
- [ ] decision rule and multiple-comparison correction applied as written
- [ ] verdicts use the fixed wording; scope sentence present
- [ ] superseded claims corrected in place, with date
