# Decision, 2026-09-29: stop the IC-versus-rho collector line and the `m ≥ 4` algebraic route at these sizes

This records a stop, not a negative theorem. Each item below names its evidence, the
boundary it was tested at, what stays uncertain, and what would reopen it.

## What stops

1. **More tournament rounds of the pair and triple collectors against rho.** This covers
   `pair_table`, `triple_counted` and the `pair_or_triple` switch, and any further constant
   tuning of them.
2. **The algebraic `m = 4` subspace decomposition** (chained `S₃`), on the frozen engine and
   on the current one, at `n ≤ 19`.
3. **An `m = 5` audit.** It was gated on a lever that lowers the degree's slope in `ℓ`, and
   the one lever the library can test lowers only its constant.

The campaign's own rule ([PLAN.md](PLAN.md)) stops an unsuccessful line after three rounds
without promotion. Rounds 0023, 0024 and 0025 produced no win against a matched rho.

## Evidence

| claim | measurement | where |
|:--|:--|:--|
| No IC collector is below rho in instructions | against the lean rho, the matched walk in the IC arm's own arithmetic: pair 1.36–2.56× on the panel and 2.1–13.7× at the crossover cells; switch 1.36–2.56× and 1.35–5.0×; both final stages, 12 cells, `n` 13–61 | [RESULTS.md, round 0025](RESULTS.md#round-0025-a-lean-rho-and-the-pairtriple-switch) |
| Round 0017's eight-cell win does not survive a matched rho | it fails at `n23a1` (1.126); past the crossing a matched rho wins by 1.9–9.3× | [ROUND24](ROUND24-pair-inverse.md), RESULTS round 0024 |
| Table decompositions cannot move the exponent | generic arity-`j` tables cost `r^{j/(2j−1)}` | [DECOMPOSITION-SURVEY.md](DECOMPOSITION-SURVEY.md) §3.5 |
| `m = 4` is closed on the frozen engine | `ĉ = 0.985` [0.983, 1.018] against `c* = 0.25`; controls pass | [ic_m4_exponent_audit_20260928](../../ic_m4_exponent_audit_20260928/RESULTS.md) |
| `m = 4` is closed on the current engine | `ĉ = 1.060` [1.058, 1.093]; 1.4–3.2× cheaper per cell, but steeper | [ic_m4_head_engine_20260929](../../ic_m4_head_engine_20260929/RESULTS.md) |
| Torsion symmetrisation is a constant lever | refutation degree rises about 1 per unit `ℓ`, `s̄ = 1.033` | [ic_symmetry_lever_slope_20260929](../../ic_symmetry_lever_slope_20260929/RESULTS.md) |
| The one IC advantage left is a constant in time | in-process 0.74–0.84 of the strongest rho at the four cells with `n ≤ 19`, parity at `n = 23`; it comes from higher instructions per second (1.65–2.03×), mostly because rho zero-fills a fixed 160 KB cache per solve | [walltime_isolation_20260929](walltime_isolation_20260929/RESULTS.md), [walltime_strong_rho_20260929](walltime_strong_rho_20260929/RESULTS.md), [throughput_gap_20260929](throughput_gap_20260929/RESULTS.md) |

**Classification (AGENTS.md §3).** The small-cell time lead is **engineering**. Nothing in
this line is an **advance**: no ratio to the generic floor fell with `n`.

## Scope, and what stays uncertain

- **Sizes.** Koblitz curves at `n` 13–61 for the tournament and `n ≤ 19` for the algebraic
  audits. There is no measurement near `m = 83` (AGENTS.md §8a) or at the challenge's 131.
- **Engines and solvers.** One IC implementation family. One Gröbner engine family
  (Macaulay/F4 degree growth), on two engine versions.
- **The throughput gap is attributed, but no fix has been measured.** The IC arm runs
  1.65–2.03× more instructions per second than rho, and a batched rho does not close the gap
  at these sizes.
  - A cache and branch simulation
    ([throughput_gap_20260929](throughput_gap_20260929/RESULTS.md)) traces most of it to rho
    zero-filling a fixed 163,840-byte cycle-detection cache (`RECENT_SLOTS`) on every solve.
  - A rho with the cache sized to the walk would probably remove IC's small-cell time lead.
    That rho has not been built or measured.
  - It is a constant either way.

## What would reopen it

Any one of these reopens the named item, and a proposal citing it should link this file.

- **A decomposition that is not a table and not the chained `S₃` system**, with a
  pre-registered exponent audit. Candidates are Nagao / Riemann–Roch decomposition
  (survey §3.3) and function-first relations.
- **A symmetry not tested here:** the symmetric-group action on summands, or Frobenius-orbit
  coordinates (survey §3.2). This reopens item 3 only if an `m = 3` ladder like
  `ic_symmetry_lever_slope_20260929` measures a slope `≤ 0.15`.
- **Any engine or solver family** that measures `ĉ < 0.25` at `m = 4` (it reopens item 2),
  or that runs beyond 64 unknowns, so that `n ≥ 23` can be audited.
- **A published degree bound** for Weil-descended or symmetrised summation systems, in
  either direction.
- **A rho stronger than the lean rho at small `n`.** It would shrink or remove the one
  constant-time lead. It cannot create a result.

## What is not stopped

- **The tooling and the rules** stay in force: the isolated benchmark runner, the AGENTS.md
  §10 isolation rule, and the evaluator-r25 checker.
- **Other lanes on the scoreboard** are not decided here: compact orbit, the Koblitz
  collection thread, and the F4/F5 backends.
- **The protocol for any reopened round** is in
  [next-proposal-single-target.json](next-proposal-single-target.json).

## Addendum, 2026-09-30: the sized-cache rho

This addendum supersedes one sentence above and changes no decision. The scope-and-uncertainty
bullet that says a rho with the cache sized to the walk "has not been built or measured" is
out of date. [sized_cache_rho_20260930](sized_cache_rho_20260930/RESULTS.md) built and measured it,
under the isolated runner and a registration made before the run:

- The sized rho takes 0.875–0.913 of the lean rho's wall time at the four small cells.
- IC is then 0.98–1.03 of the sized rho, so its small-cell time lead is gone. IC's one residual,
  2% at `n19a0`, is inside the A/A residue.
- At `n23a0` the sized rho is 4% faster than IC (`incumbent/sized` 1.043 [1.015, 1.074]).
- The walk is unchanged on all 60 timed cases, and every one of 1,536 recorded rho runs verifies.

The reopening item "a rho stronger than the lean rho at small `n`" is met for the time lead. The
sized rho is now the reference rho at `n ≤ 23` for any reopened round.
