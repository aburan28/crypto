# Every candidate, one harness: preregistration

Registered 2026-10-03, before any session ran. The two specs beside this
file are hashed by `ecbench` at run time; the session records carry the
spec id, so a spec edited after this commit would show.

## Question

What does every ECDLP method this repository implements cost, cold, to
solve one previously unseen target, in one unit, on the same workloads,
and where does each index-calculus candidate stand against the matched
rho on the same point? This is the cross-method question AGENTS.md §12
assigns to `ecbench`; it is not a tournament round and promotes nothing.

## Boundary

- **Floor:** the curve's generic floor `√(π/2A)` in `S`, recorded per run
  (`A = 2` on the prime curves, `A = 2n` on the Koblitz curves).
- **Reference:** `rho.negation` on the prime curves; on the Koblitz curves
  `rho.signed_frobenius_strong`, the strong single-target reference the IC
  measurement rules require since 2026-10-01. `rho.signed_frobenius` is
  carried as the operations-only baseline.
- **Unit:** `S = total group-addition equivalents / √r`, cold end to end,
  every phase charged, as `ecbench` defines it (README §7). BSGS's memory
  is not in `S` and is read as the time–memory trade, never a finding.

## Workloads

- Prime family: `prime_search` at 18, 20, 22 and 24 bits, seed 59297;
  Koblitz family: `koblitz` `a=1 n=17`, `a=0 n=19`, `a=0 n=23`,
  `a=1 n=29`. Four sizes per family (§5). Eight planted targets per curve,
  `target_seed 20261003`; three interleaved rounds, one warm-up; seed
  20261003. The `n=29` curve's prime subgroup is `2^15.37`, smaller than
  the others': sizes are read by `r`, not by `n`.
- Arms: every rho variant the harness has; every BSGS variant; the
  kangaroo; `ic.pipeline` with every factor base and oracle that applies
  to the family (`prime-abscissa` × {`mitm` folded, `subtract`};
  `koblitz-orbit` × {`mitm-frobenius` m=2, m=3, `subtract`}). The
  `descent-algebraic` oracle is excluded: its solver budget makes runs
  nondeterministic and it is measured by its own protocols. A control arm
  repeats the reference for the A/A noise floor.

## Predictions, scored at each family's largest `r`

- **P1 (calibration holds):** each generic method's mean `S` lands within
  its 95 % interval of its analysed constant, as in
  `research/ecbench_calibration_20261002`: plain rho `√(π/2)`, negation
  rho `√(π/4)`, signed-Frobenius rho `√(π/4n)`, BSGS 1.5 / 4/3 / 1,
  kangaroo 2. Rho's interval is expected to be wide (CV ≈ 0.5 at 24 runs).
- **P2 (no IC arm below the reference):** every `ic.pipeline` arm's
  `Σ S_B / Σ S_A` against the family's reference is above 1 with its
  interval excluding 1, at every curve. This is the repository's standing
  verdict; the run tests it on fresh targets with the harness's accounting.
- **P3 (direction with size):** each IC arm's ratio to the reference does
  not fall with `r` across the four sizes. A falling ratio would be
  reported as such with its per-curve rows, never pooled.
- **P4 (correctness):** every execution verifies (`[k]G = Q` and
  `k = planted`, checked by the runner and again by `ecbench verify`);
  every replay in the audit is identical.

## Success, stop and what is inadmissible

- The run is a success when P4 holds and the table is complete; the
  predictions are scored as written whether they hold or not. A failed
  prediction is a result, not a reason to rerun.
- A timed-out or exhausted run stays in the table as its own status and
  makes any comparison involving it `incomplete`; it is not dropped.
- Inadmissible after seeing results: changing curves, seeds, targets,
  rounds, parameters, the unit or the reference; dropping an arm;
  editing a session.
- **Wall time** is reported at whatever isolation level the host earns.
  This host is a four-vCPU cloud VM; L2 is not expected. Operation
  counts are the result; wall time is descriptive unless `admitted`.

## Classification

Whatever the table shows, the class under AGENTS.md §3 is fixed in
advance: no method is changed here, so nothing is an *advance* or
*engineering*; the session is a measurement of existing methods on fresh
targets. If an IC ratio reads below the reference the first question is
the reference, per the 2026-09-23 control.
