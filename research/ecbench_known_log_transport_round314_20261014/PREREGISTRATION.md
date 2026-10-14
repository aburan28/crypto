# Round 314 preregistration: exact-log transport versus rho

Date: 2026-10-14

## Open request

Measure and improve the end-to-end index-calculus line until it reaches
parity with a square-root Pollard-rho reference.  The motivating P-256
object remains
`icv1-fp256-t89188191154553853111372247798585809583-f188c491` with
`FB1h2f8621cda105`.

## Hypothesis

The strongest complete small-instance candidate is dominated by rebuilding
the factor-base logarithms.  If signed-Frobenius factor-base representatives
are selected as known scalar multiples of the generator, their logarithms
transport exactly through the orbit coefficients.  Removing the rank phase
may put the complete cold solve, including base construction, folded pair
table, target descent, witness replay, and scalar replay, within 1.30 times
the matched strong rho cost.

This is a mechanism screen, not a P-256 speedup claim.  A result promotes to
the P-256 Dickson program only if the transported points are members of
`FB1h2f8621cda105` and their logs are obtained without solving an ECDLP or
charging hidden precomputation.

## Frozen method boundary

Candidate method ID: `ic.transported_logs`.

For `K` columns, scan `[s]G` for increasing positive `s`, retain the first
`K` distinct signed-Frobenius orbits, and materialise every signed orbit
member.  The column log is the retained `s`; an orbit member has log
`s * (+/- lambda^j) mod r`.  Charge the scan, construction, pair-table
construction, target validation and shifts, failed decomposition probes,
witness replay, log combination, and final `[d]G = Q` replay.  Report all
native field, hash, and allocation work that remains unpriced.

The target is never supplied as a scalar.  The public-target law and the
runner's independent point check remain unchanged.  Every reported
decomposition and recovered scalar must replay as a full group equality.

The matched reference is `rho.signed_frobenius_strong` on the same one-target
workload.  The primary unit is cold end-to-end `ecbench.gae`; wall time and
instructions are diagnostics unless their required isolation is met.

## Development and freeze rule

1. Correctness tests use the registered `a=0, n=17` Koblitz instance and
   deterministic point-only targets.
2. A development-only K sweep may be run on `a=0, n=29` and `a=0, n=37`
   with target seed `20261014314001`.  It may change only the deterministic
   column-selection rule.
3. After that sweep, freeze one rule mapping `(r,n)` to `K`.  The initial
   rule is the nearest positive integer to
   `(3*r/(2*n^3))^(1/4)`.  Any replacement and its development evidence must
   be recorded before opening holdout results.
4. Holdout targets use seed `20261014314002`, are public hash-to-subgroup
   points, and may not influence the frozen rule.

## Holdout panel and controls

Use four registered one-word Koblitz curves (`a=0`) at degrees 19, 23, 29,
and 37.  Each curve has eight independent one-target workloads, five measured
rounds, and one warmup.  Each spec contains:

- the transported-log candidate;
- the matched strong signed-Frobenius rho reference;
- an identical transported-log A/A control.

Use alternating order, one thread, a fresh process per execution, the normal
resource lock, and the strongest isolation available in this environment.
No execution may be silently discarded.  Audit every record with
`ecbench verify --replay-all`.

## Promotion gates

All gates are conjunctive:

1. Every candidate and reference record verifies; no wrong answer, exhausted
   run, crash, timeout, false positive, or false negative occurs.
2. Candidate/control A/A is statistically consistent and has identical
   counted work on matched deterministic inputs.
3. At degree 37, mean cold candidate `S` divided by mean cold reference `S`
   is at most 1.30.  Online-only parity does not satisfy this gate.
4. A log-log fit across all four degrees has a candidate work exponent no
   greater than 0.55 in subgroup order.  The fit, residuals, and uncertainty
   are reported; four points are not called an asymptotic proof.
5. Peak measured RSS and the exact pair-table entry count are reported.  No
   omitted or reusable setup is removed from cold cost.
6. P-256 promotion requires an exact, reproducible membership/log-transport
   certificate for the Dickson base.  Selecting arbitrary known multiples of
   `G`, assuming their Dickson membership, or counting a DLP oracle as free
   fails this gate.

If gates 1--5 pass but gate 6 fails, report parity only for the Koblitz
analogue and a negative P-256 transfer result.  If the parity gates fail,
publish the measured ratio and dominant phase rather than relabelling a local
phase improvement as end-to-end parity.

## Iteration and stop condition

One post-holdout engineering iteration is allowed only for a defect or a
phase that the frozen records identify as at least 20% of candidate counted
cost.  It must get a new method identity and fresh target seed.  Stop Round
314 after that iteration or after the promotion gates pass.  Further rounds
require a new preregistration; “continue until parity” does not authorize
discarding failed branches or adapting to hidden holdouts.

## Amendment log

1. Before any development or holdout target was run, the first unit test
   showed that `koblitz_instance(0, 17)` returns no admissible subgroup in
   this harness.  The correctness-only fixture is therefore `a=1, n=17`,
   which is the instance already used by the shared-rank correctness tests.
   The development seeds, holdout curves, metrics, and gates are unchanged.
2. Before any development or holdout target was run, the registry check
   showed that degree 29 is registered only for `a=1`; there is no `a=0`
   degree-29 workload.  Both the development panel and the degree-29 holdout
   entry therefore use `a=1, n=29`.  Seeds, all other curves, metrics, and
   gates are unchanged.
3. After the declared development sessions and before opening the holdout
   seed, the curve panel is changed from degrees 19, 23, 29, 37 to the four
   registered one-word curves with strictly increasing subgroup sizes:
   `(a,n)=(0,37),(1,43),(1,47),(0,41)`, whose `log2(r)` values are 27.78,
   32.11, 36.64, and 39.00.  The original degree-based panel included the
   degree-29 curve with only a 15.37-bit subgroup and could not support the
   preregistered scaling fit.  This change uses registry facts only, not
   holdout results.  The fresh holdout seed, metrics, and gates are unchanged.

## Frozen rule before holdout

The 64-target development confirmation selected K=11 at degree 37.  The
frozen analytic rule is
`K = round(1.34 * (r / n^3)^(1/4))`, clamped to at least one, giving K values
11, 21, 43, and 71 on the four holdout curves in increasing subgroup order.
No holdout result was opened before fixing these values.
