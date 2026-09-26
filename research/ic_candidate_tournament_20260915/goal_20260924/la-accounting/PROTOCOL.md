# Failed relation-LA attempt accounting

This follows the PDP outcome correction in PR 789 and remains an admission
correctness change. No improvement round is consumed and no algorithmic speedup
is claimed. The first bounded round retained the incumbent; two rounds remain.

`FactorBaseLogSolver::try_solve` cloned its report, ran linear algebra into that
clone, then used `?` on an unsuccessful solve. Consequently, failed attempts,
their elapsed time and sparse-solver diagnostics vanished. On success only the
attempt counter was copied back, leaving the persistent report without its
successful time, verification status or sparse diagnostics. This defect is in
the generic path, not the archived optimized producer used in round one.

Before admitting that path, require the persistent report to retain every
actually attempted matrix solve. A check with too few rows does not count as a
solve. Dense and sparse modes obey the same contract. The fix updates the stored
report directly and records whether the latest attempt returned a certified table.

Frozen toy controls use `K_0/F_(2^9)`, the checked-in autolab dependency lock and
one release-test thread. A homogeneous singular matrix over a 36-point-requested
subgroup-orbit base isolates failure accounting: two attempted solves must remain
two attempts, cumulative time must retain a deliberately nonzero prior value,
and sparse diagnostics must survive. These manually supplied rows are a matrix
unit control, never measured IC relation evidence. A second control collects 600
ordinary toy queries over factor-base index zero with the existing collector
seed, then checks that repeated successful solves keep the same report in storage
and in the return value and that every returned column log verifies in the group.
Existing collected-relation/precompute equivalence controls remain required.

Stop on any failed control. The Linux integration workflow is the acceptance gate;
local timing in these tests is not a performance sample. Complete online/cold
cost, Ir, S, rho ratio and speedup remain unknown. The generic worker still needs
exclusive phases and full per-query statistics before tournament admission.
