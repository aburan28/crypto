# Follow-up: index calculus where the first session's base did not exist

Registered 2026-10-03, after the first Koblitz session completed and
before the follow-up session ran. The first session is kept as it is.

## What happened

In `sessions/koblitz` every `ic.pipeline` arm errored on two of the four
curves, before doing any work: `no invariant subspace for divisor [0, 1]
on this curve`. 192 records, status `error`, exit 0. The cause is in the
base builder, not the solver:

- `koblitz-orbit:divisor=0;1` names the factors of `x^n − 1` at indices 0
  (`x + 1`) and 1 of `all_factors_of_x_n_minus_1`, which enumerates
  irreducible factors of degree at most 24.
- `icv1-f2m19-t797-b6cf2467`: `ord_19(2) = 18`, so the only other factor
  has degree 18 and the two together span the whole field; no proper
  invariant subspace results.
- `icv1-f2m29-tm40309-30c52b96`: `ord_29(2) = 28 > 24`, so no non-trivial
  factor is enumerated and index 1 does not exist.

On `n = 17` (`ord = 8`) and `n = 23` (`ord = 11`) the divisor exists and
every IC run verified.

## What this session does

Measure index calculus on those two curves with the `binary-subspace`
base (abscissae in a polynomial-basis `F_2`-subspace), which needs no
Frobenius-invariant factor, under the two oracles that accept it:
`mitm:m=2` (unfolded; a folded table is forbidden on a base that is not
Frobenius-closed) and `subtract`. Dimensions 6 and 8. The reference
(`rho.signed_frobenius_strong`), the operations baseline
(`rho.signed_frobenius`) and the control are re-run in the same session
so every comparison is within one session. Targets, seeds and rounds are
the first session's.

## Predictions

- **P5:** every IC arm verifies on every workload of both curves.
- **P6:** every IC arm's `Σ S_B / Σ S_A` against the reference is above 1
  with its interval excluding 1, on both curves; the subspace base is
  expected to cost more than the Frobenius-invariant base did on `n = 17`
  and `n = 23`, because its columns are not merged by Frobenius.
- **P7:** dimension 8 costs more than dimension 6 in `S` at these sizes
  (the table dominates below the crossover).

## Inadmissible

The same list as `PROTOCOL.md`. The first session's error rows stay in
its table and are not replaced by these.
