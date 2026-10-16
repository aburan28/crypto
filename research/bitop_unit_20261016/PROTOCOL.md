# Protocol: ECDLP methods priced in `bitop/v1`

Frozen 2026-10-16, before any session of this study ran. Changes after the
first session are listed under "Deviations" with their reason.

## Question

What does each ECDLP method in the repository cost when every operation it
performs is charged in one machine-independent unit? How does that cost scale
with the subgroup order `r`?

The unit is `bitop/v1` (`docs/units/BITOPS.md`): two-input Boolean gates per
primitive, under frozen formulas.

The study answers with fitted exponents and constants, each carrying an
interval. It makes no claim that any method beats another beyond what those
intervals separate. It makes no security claim about any deployed curve.

## Instruments

- **Records.** Every record is priced by `ecbench::work` from its own counts.
- **Completeness.** A record's pricing is *complete* when it has:
  - no counter without a rule (`unpriced`);
  - no work its method's code is known not to count (`known_uncounted`).
- **Instruction audit.** `ecbench audit-ir` checks completeness independently.
  It replays a run under Callgrind and divides the solve's instructions by its
  counted bit operations. A method whose ratio exceeds twice the matched rho's
  on the same curve fails the audit.

## Sessions

1. **New prime ladder.** `specs/prime-ladder.json`:
   - one prime-order curve per bit size from 16 to 26, found by `prime_search`
     with seed 20261016;
   - eight public targets per curve;
   - six arms: `rho.negation` (reference), `bsgs.negation`, `kangaroo.vow`,
     and `ic.pipeline` three ways: MITM with bases of 32 and 128 abscissae,
     and `rootfind-subtract:m=3` with a base of 32.
2. **Re-runs of committed specs.** The binary is the instrumented one, so the
   new counters are present. Each spec is copied unchanged into `specs/`:
   - `ecbench_frobenius_orbit_vs_rho_20261014` a0-n13, a0-n19, a0-n23 and
     a0-n31: Koblitz toy tier, Frobenius-orbit and pdp3 index calculus,
     signed-Frobenius and strong rho;
   - `fiber_ic_scaling_20261011` scaling_h4: prime, 8 curves, fiber-closed
     MITM at five base sizes;
   - `ecbench_sat_pdp_20261009` count-qualified-main: SAT decomposition;
   - `f4_f5_f6_cross_family_20261011` K13: the F4 family through pdp3;
   - `ecbench_yield_sweep_20261004` koblitz: Buchberger, crossbred and F4
     descent.

Sessions run at isolation L0 (`--cpus none`). Bit-operation counts do not
depend on the host or on isolation, and wall time is not used.

## Measures (all from `ecbench`, none computed outside it)

- **Primary.** For each arm, `ecbench bound fit --unit bitop/v1`:
  - the exponent `α` of `bitops = C·r^α`, by least squares over all verified
    runs, with a two-stage cluster bootstrap 95% interval (2000 resamples,
    seed 20261005);
  - the same fit in addition-equivalents: `bitops / gates(group_add)` at each
    size, which takes out the field-arithmetic growth `~(log r)²`;
  - the constant as the ratio to the generic floor `√(πr/2A)` additions;
  - per-phase exponents where a phase has at least 4 sizes.
- **Secondary.**
  - `ecbench amortize`: per curve, the window of target counts in which an
    index-calculus arm costs no more than batch rho;
  - `ecbench work --summary` over every committed session: how much of the
    archive is completely priced, and what the rest misses.

## Decision rules, fixed now

- An exponent is reported only with its interval. Digits beyond the interval's
  width are not reported.
- **Positive control.** The generic methods (rho, BSGS, kangaroo) have the
  theoretical exponent ½ in addition-equivalents. If any generic arm's 95%
  interval in addition-equivalents excludes ½ by more than 0.05, the
  instrument is suspect, and the study reports that before any
  index-calculus figure.
- An arm whose pricing is incomplete, or that fails the audit, is reported as
  a lower bound and is not admitted as a bound.
- A size with any unverified run makes the arm inadmissible (ecbench's own
  rule). It is reported, not dropped.

## Deviations

1. **Sessions re-run after fixes.** Four fixes landed after the first
   sessions ran:
   - the fiber oracles' tables counted their insert attempts;
   - oracles with no table (descent, pdp-split, and pdp3's algebraic
     engines) declared their set-up counted;
   - pdp3's root-finding engine charged its table's inserts.

   Sessions that had already run on the earlier binary were discarded
   uncommitted and re-run unchanged from their frozen specs: the fiber
   ladder, the four Koblitz specs, SAT, F4 K13 and the yield sweep. The
   prime ladder ran once; its arms use none of the changed code.

   A later fix is to the reader only. The binary that ran these sessions
   writes `null` for the logarithm of a phase with zero gates, and `ecbench`
   now reads that as absent instead of failing to parse it. Records are
   unchanged.
2. **Instruction-audit panel.** The audit draws on the committed specs of
   the studies whose engines it checks, not on this study's sessions, so
   that each audited run has a committed record to compare with. A second
   audit covers the prime ladder: one run per arm on its largest curve.
3. **Koblitz index-calculus arms.** These arms are tuned per curve size,
   each with a different base, so they have one constant-only bound per
   size. Only the two rho arms are fitted across the four sizes.
