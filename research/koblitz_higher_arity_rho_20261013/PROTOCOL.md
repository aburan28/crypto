# Stride-ordered higher-arity extraction versus matched rho

**Status:** preregistered before implementation or measurement, 2026-10-13.

## Requested scope

| Requirement | Planned evidence | Acceptance state |
|---|---|---|
| Change the relation/extraction shape toward ECC2K-130 | Add a six-point, stride-ordered chain kernel that visits the same pair cells as the existing `S4Strided` kernel but obtains the complete S4 root set through two independently checked S3 chains | pending |
| Measure the rho baseline | Run `rho.signed_frobenius_strong` on the exact public target(s) used by the IC arms, in one native `ecbench` session | pending |
| Preserve full accounting | Record cold setup, the five exclusive IC online phases, native counters, failures, timeout/resource exits, and independent scalar replay | pending |
| Test the ECC2K-130 transfer claim | Run the registered n=83 confidence gate if the native harness and resource envelope admit it; otherwise preserve the exact implementation/resource blocker and leave transfer unestablished | pending |

## Boundary and reference

The primary unit is

`S = cold total group-addition equivalents / sqrt(r)`.

The generic signed-Frobenius floor is `sqrt(pi/(2A))`, with `A = 2n`.
The matched reference is `rho.signed_frobenius_strong`, one target at a
time.  An IC/rho ratio is admissible only when both arms recover and verify
the same public point and all costs represented by `S` are priced.  Any
unpriced S3/S4 or lookup work makes the ratio a lower-bound diagnostic.

## Hypothesis

For a fixed pair cell `(s1,s2)`, the existing quartic S4 polynomial has roots
`x(R +/- s1 +/- s2)`.  Computing the two points `R+s1` and `R-s1`, then solving
`S3(x(R+s1),x(s2),X)` and `S3(x(R-s1),x(s2),X)`, should recover the same root
set while preserving the diagonal scan order that de-clusters six-sum hits.
The candidate changes the extraction primitive; it does not change the factor
base, target, pair-cell order, relation verifier, or linear algebra.

## Frozen implementation and correctness gates

- Baseline revision: `759978250c1ce53d7b87bd8a1c2304781fcc4fc4`.
- Curve gate: registered Koblitz `a=0`, n=41
  (`icv1-f2m41-tm2308219-7f48b14a`), followed by n=83 only if supported.
- Factor-base gate: `compact-orbit-scan`, 20 columns, identical raw scan cap
  for baseline and candidate.
- Baseline extraction arms: `4sum` and `b2` (stride-ordered quartic S4).
- Candidate extraction arm: new stride-ordered two-S3-chain kernel.
- For every tested pair cell, candidate roots must equal the existing S4 root
  set after sort/dedup and every emitted relation must pass both fast and
  general curve re-addition.
- Planted, random-target, non-decomposable and budget-stop tests are retained.
- A disagreement, verification rejection, or changed pair visitation order is
  a correctness failure, not a performance result.

## Measurement plan

1. Record the host capsule, compiler, source revision and baseline binary hash.
2. Preserve an unmodified baseline probe of `b2` at n=41, K=20.
3. Build a native `ecbench.spec/v1` with public one-target workloads and arms:
   strong rho, 4-sum IC, quartic-strided six-sum IC, and chain-strided six-sum IC.
   Include a second identical strong-rho arm as the preregistered A/A wall-noise
   control; it does not change the matched operation-count reference.
4. Run the arms interleaved under the same one-thread resource envelope.  Keep
   every error, timeout, exhausted relation stream and isolation downgrade.
5. Audit with exact replay, save the one-unit table and pairwise comparisons,
   and report cold and online intervals separately.

The first session is a bounded discriminator, not a scaling claim: one n=41
size cannot establish an exponent and one host cannot supply independent
validation.  A wall-time statement requires admitted isolation and the
repository's A/A/noise gates; otherwise wall time is descriptive only.

## Success and stop conditions

The candidate passes its implementation gate if it is root-set identical to
`b2`, emits only independently verified relations, and completes the same
n=41 public workload.  It passes the local extraction-cost gate if its recorded
native work per verified relation is below the existing quartic-strided arm.

This round stops on the first of:

- completed audited session;
- a 900-second per-run timeout;
- a fail-closed correctness disagreement;
- memory pressure or another recorded resource refusal.

No outcome at n=41 establishes ECC2K-130 total-work parity.  The strategic
claim remains open unless the n=83 gate completes and the resulting full cold
`S` is competitive with matched rho; transfer from n=83 to n=131 remains an
extrapolation even then.

## Environment note

The required `conductor check` was attempted before editing with the source,
study and benchmark paths reserved above.  This runtime has no `conductor`
executable (exit 127) and the checkout has no `.conductor/WORKFLOW.md`; no
scope-conflict certificate is therefore available.
