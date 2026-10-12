# Round 0005 pre-registration: complete cold batches of 16

This is a separate workload panel. It does not establish single-target parity.
Each job constructs one curve and 16 independent public-hash targets, charges all
factor-base/table/precomputation work once, solves every target, performs every
internal check and general final scalar check, and verifies all outputs with the
independent Python checker. Rho solves the same 16 targets through the existing
per-target signed-Frobenius solver API on the same constructed curve. No cache or
precomputation outside the timed/profiled job is admitted for either method.

Parent: round-0004 selection `folded_lift_batch4`, conditional on full parent promotion.
Source: `/home/ubuntu/crypto/research/ic_candidate_tournament_20260915/runs/round-0004/source_candidates/folded_lift/source`. Fresh seed 2026091605; target count 16; five confirmation
cells, 60 fresh job fixtures (960 targets per repetition and arm), three
repetitions, independent replay. Budget 1800 paired jobs. Single logical CPU 7,
8 GiB address-space cap, 60-second watchdog per process. Builds avoid both threads
of the parent's measurement core, and measured execution waits for its completion.

Three source candidates: lazy first descent probe, cached fast verification of
recovered logs, and both. Base support, m=3, source accounting boundary and all
verification obligations stay fixed. The first-probe candidate keeps the same
probe order, avoids repeating the first failed probe, and respects the trial cap.

Reference is matched rho. The floor remains K instructions for K required
independent relation columns; it is weak and cannot establish a non-generic
advance. Class: engineering. Within this panel, promotion requires >=20% lower
instructions and native wall, upper paired 95% limits <1 and every cell <=1.10,
on confirmation and replay. Parity requires candidate/rho upper paired 95% limits
and every cell ratio <=1.10 in BOTH metrics on BOTH final stages.

The preceding 36-trial development screen completed correctly but did not meet
the parity rule: ratio estimates near one hid a regression in the smallest cell.
That screen selects hypotheses only; its cases are not reused for confirmation.
No mathematical or family-wide speedup is claimed. Comparisons use the shipped
rho implementation; additional cross-target rho optimizations are outside this
measured comparison and must be measured if introduced.
