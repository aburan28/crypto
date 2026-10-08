# Protocol D-3: can relation-matrix filtering lift the linear-algebra cap enough for a three-summand oracle?

Frozen 2026-10-07, before any instrument is built.  **Stage diagnostic.**
`S`, end-to-end cost and speedup are **unset**.  Status: **PENDING**.

## Derivation (stated before measuring)

The gap formula of
[`../geometric_v_linear_20261006/README.md`](../geometric_v_linear_20261006/README.md)
reads: at `n = 131`, a trial that resolves `k` summands at polynomial cost
leaves a bits-only gap to rho of about `70.19 − (k−1)·c`, where `2^c` is
the number of relation columns, and the linear algebra at `m·2^{2c}`
operations against rho's `2^{60.81}` caps `c` near 30.  Hence `k ≥ 4`.  The
cap is the only term in that formula that engineering could move.

Write the solve as block Wiedemann on an `N × N` matrix with `w` nonzeros
per row, cost about `2wN²`.  With `N = 2^c` and `w = 3` the cap is
`c ≤ 29.1`.  Filtering in the number-field-sieve sense, singleton and
clique removal followed by merging, replaces `N, w` by `N/f, w'`.  For a
three-summand oracle to suffice the formula needs `2c ≥ 70.19`, so
`c ≥ 35.1`, and the filtered solve must still fit:

    2 · w' · (2^{35.1} / f)² ≤ 2^{60.81}   ⇔   f ≥ 2^{5.2} · √w'.

At `w' = 10` that is `f ≥ 117`; at `w' = 4` it is `f ≥ 74`.  Published
NFS filtering gains are 5× to 20× in dimension at row weights that grow
into the tens.  So the arithmetic predicts that filtering **cannot** lift
the cap far enough, and this protocol measures the factor on the
repository's own relation matrices to replace that prediction with a
number.  The question matters because it is the last escape from `k ≥ 4`
inside the fix-some-summands family.

A separate observation fixes the instrument: the compact-orbit rank stage
collects exactly `K` rows for `K` columns (guided rank, zero surplus), so
there is no excess to filter.  Filtering needs excess rows, and excess
rows cost probes; the protocol charges them.

## Instrument (Rust, to build in the follow-on PR)

1. Extend `src/cryptanalysis/koblitz_relation_solver.rs` with a filtering
   pass over a relation set held as sparse rows over `Z/rZ`: singleton
   removal (a column in one row), clique removal with a target excess, and
   Markowitz-style merging up to a weight bound `w_max`.  Every removed or
   merged row is recorded so the final solution can be back-substituted
   and verified in the group as today.
2. Collect relation sets with excess `x ∈ {10%, 25%, 50%, 100%}` over `K`
   on `icv1-f2m41-tm2308219-7f48b14a` (`K = 255`),
   `icv1-f2m53-tm56619371-dac20a85` (`K = 440`) and
   `icv1-f2m61-t158598901-ab42b6c5` (`K = 600`), three seeds each, with the
   compact-orbit rank stage allowed to continue past full rank.  Charge the
   extra probes.
3. Report `f = N / N'`, `w'`, the implied Wiedemann cost `2w'N'²` against
   the dense and unfiltered sparse costs, and the back-substituted,
   group-verified solution.
4. Compute the implied admissible `c` at `n = 131` from the measured `f`
   and `w'` beside the `c = 35.1` requirement.

## Predictions (pass/fail)

- **F1 (singletons and cliques).**  At 50% excess, singleton plus clique
  removal cuts `N` by at most 3× on every curve.
- **F2 (merging).**  Merging to `w_max = 12` brings the total factor `f`
  to at most 10× on every curve, with `w'` in `[6, 12]`.
- **F3 (cap).**  The implied admissible `c` at `n = 131` rises by at most
  2.5 bits over 29.1, so stays below 32 and well below 35.1.
- **F4 (cost of excess).**  Collecting 50% excess costs at least 1.4× the
  probes of the exact-rank collection, since probes per relation are flat
  in the row count.
- **F5 (correctness).**  Every filtered system back-substitutes to the
  same column logs as the unfiltered solve and verifies in the group.

## Decision rule (registered)

If F1 to F3 pass, the `k ≥ 4` requirement stands with a measured margin,
and the open problem of the geometric-V note is confirmed as stated.
Class: **boundary**.  If `f ≥ 40` at `w' ≤ 10` on any curve, the
requirement weakens toward `k = 3` and D-2's three-summand reach becomes
the binding number; that outcome is recorded as **accounting**, since no
algorithm changed, and the gap table is recomputed.

## Stop condition and inadmissible moves

Bounded: three curves, four excess levels, three seeds.  Stops when F1 to
F5 are scored.

Inadmissible: dropping the excess probes from the charge; merging without
a weight bound; reporting `f` without `w'`; quoting the implied `c` as a
measurement at `n = 131`, which it is not; treating a filtered dimension as
a change to the product law.
