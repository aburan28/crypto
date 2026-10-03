# Native prepared SAT migration

This work continues the complete source-bound F4/F5 and SAT one-target goal.
The dashboard redesign and completed F5 v3 development control are accepted in
[PR #1255](https://github.com/aburan28/crypto/pull/1255). The underlying goal is
still active: no fresh paired family comparison has been completed.

The current repository rule forbids Python research execution, orchestration
and verification. The next prepared-SAT controller must therefore execute
entirely in Rust with the retained native exporter and CryptoMiniSat backend.
Historical Python sources, measurements and certificates stay immutable and
keep their provenance. The existing `icprog` harness is the native home; this
work does not introduce another tournament or scoreboard.

## Delivered verification path

`icprog sat-source-replay` is a bounded postexecution replay of the already
disclosed n17 fixture. It does not accept a new curve, public target, witness,
solver executable, seed or conflict budget. Its source byte pins are the copies
accepted in PR #1127, taken from trial 1 of the closed prepared SAT control.
The command checks those committed copies directly; it does not reopen the
original publication archive. The archive/execution hashes in its report are
historical associations, not a new archive custody claim.

Before substituting the source witness, native independent curve arithmetic
reconstructs the SAT preparation from its whole-certificate seal:

- Preserve all 63 geometric points, including the killed torsion point;
  cofactor projection has 62 usable points and 29 sign/Frobenius columns.
- Reconstruct orbit representatives and coefficients, relation rows, rank
  trajectory, duplicate/dependent counts and the exact saved matrix snapshot.
- Re-add all 37 retained ordinary witnesses. Independently check the 106
  `SOURCE_UNSAT` queries for complete geometric three-summand absence. Include
  repeated indices, torsion points and identity pair sums in that proof.
- Preserve the six `CONFLICT_BUDGET_INCONCLUSIVE` records without asserting
  whether they have a decomposition. They contribute no relation.
- Solve the reconstructed matrix over the prime subgroup scalar field and
  independently replay every one of the 29 column logarithms.

The source witness uses indices `[29,51,2]`, points
`[(33,16320),(57,96968),(3,95627)]` and query `(62577,27783)` on the retained
n17 curve. The verifier derives the 51 source bits from those x-coordinates and
the three elementary symmetric coefficients; it does not trust a saved model.
It reconstructs all 716 auxiliary AND values and checks all 50 ANF equations,
2,364 ordinary CNF clauses and 50 signed XOR rows. The complete source files
are SHA-256 pinned, and their export descriptors are checked with native
BLAKE3. ANF/CNF parsing rejects bad headers, counts, variable bounds,
unterminated/truncated rows and incomplete/duplicate auxiliary definitions.

Output creation uses `create_new`; existing evidence is never overwritten.
Source inputs must be bounded regular files; Unix leaf symlinks are refused.
The checker binds its executable before arithmetic, rejects a changed executable
at the end, and retains checker-source digests and host OS/CPU
architecture. These identify this replay, not a complete scientific runtime.
No elapsed time, speedup, candidate ID or family-promotion claim is produced.

The original eight-query SAT control remains exhausted and closed. Its
`CONFLICT_BUDGET_INCONCLUSIVE` outcome remains unchanged even though one of
its source queries has an independently verified model. Finding that model by
substitution does not demonstrate that CryptoMiniSat found it within a budget.
All three historical tournament confirmation sets remain closed.

## Native execution gates still pending

This document is a migration contract, **not an executable registration**.
No new solver budget or dispatch is frozen here.

1. Port prepared target recovery, source-model lifting, exporter/CMS process
   handling and progress receipts into Rust. Reuse existing native identity,
   query-law and independent curve code. Do not spawn or embed Python.
2. Retain the exact native asset bytes and sources for exporter, CryptoMiniSat,
   Cadical/Cadiback and their build dependencies. Bind the Rust controller and
   auditor's full executed sources, Cargo dependencies, compiler/flags and
   binary digests before any scientific dispatch. Validate custody before and
   after each child and retain the actual command/input/stdout/stderr/exit code.
3. Implement bounded child deadlines, process-group cancellation/drain and
   one-use execution claims. A timeout, OOM, budget exhaustion or interrupted
   execution remains a terminal row; resumption cannot retry that invocation.
4. Match ordinary-query and target-query laws. Import preparation before the
   online clock. Charge every target-dependent export, source gate, failed
   attempt, model check, descent and independent scalar replay exactly once.
   Retain the five exclusive online phases and their exact summed interval.
5. Test native transport/source gates, arithmetic recovery, solver status
   parsing, invalid/partial models and watchdog failure paths without a new
   scientific solver dispatch. Freeze a separate development registration
   with the disclosed point, explicit budget, resource envelope, stop rule and
   original known-input disclosure. Register the complete native auditor too.
6. Execute that new registration once and publish its original terminal
   result plus independent admission. A larger conflict budget is a separately
   named method parameter, never a revision or retry of the exhausted run.
7. Port the remaining comparison controller and F4/F5 execution path to native
   code. Reconcile the full exposure union, incumbent and strongest applicable
   one-target rho reference. Rebuild on the comparison hardware and freeze
   calibration, resources, paired public points, arm order, repetitions and
   acceptance/uncertainty rules before sampling fresh targets.
8. Execute the fresh paired comparison. Admit speedup only for independently
   verified complete one-target solves on the same point and resource envelope.
   Preserve incomplete cells; online speed is primary, with preparation/cold
   costs and stage diagnostics separate. Retain diverse whole pipelines for
   recombination; local stage winners do not establish a global optimum.

Rust-native retained-evidence replay tests run on the Linux and macOS CI hosts.
Those tests establish correctness of the replay on the reported architecture;
they are not timing benchmarks, device speedup validation, new natural-yield
measurements or authorization to solve a challenge/third-party target.

## Build serialization

`src/bin/isolated_bench/busy_unix.rs` can be bootstrapped with `rustc` before
Cargo. Its Unix `flock` uses `/tmp/crypto-bench.lock`, the same advisory lock as
the retained isolation tools. The parent holds the lock while the child runs;
the command's exit code is preserved. Concurrent-child, command-failure and
launch-error controls verify serialization and release.

The macOS `isolated_bench` supports **busy only**. It does not claim CPU
affinity, quiet-host sampling or calibrated timing. The Linux timed isolation
implementation is unchanged. Bootstrap compilation itself is a small startup
step outside every scientific interval; all subsequent heavy builds/tests
take the shared lock.
