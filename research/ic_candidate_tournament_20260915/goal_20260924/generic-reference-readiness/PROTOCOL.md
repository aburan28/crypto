# Generic reference-panel admission readiness

Depends on PR 822's generic scientific admission. That PR is not yet merged at
registration. This check is a prerequisite diagnostic for reference qualification,
not a speed comparison, candidate selection or improvement round. The incumbent
and the sealed first round remain unchanged, with two improvement rounds left.

Hypothesis: the generic source-bound worker can export independently admitted
complete pipelines on each of the five already registered development cells:
17a1, 19a0, 23a0, 23a1 and 31a0. Its present negative-PDP checker accepts only
two-summand proofs at n <= 13. Determine whether the actual three-summand reports
on this panel encounter that or another admission boundary. Do not weaken the
checker to make a report pass. No confirmation or replay targets are inspected.
The held-out 29a1 cell is excluded.

Freeze exactly fifteen jobs before execution: on each cell, pair-table PDP with
dense LA, pair-table PDP with sparse LA, and the existing signed-Frobenius rho
with requested width four. Each job takes the same one supplied public target
within its cell, generated outside measurement from seed 2026092561. The
algorithm seed is 2026092561. The IC factor-base recipe is subgroup_orbits with
seed 43 and requested point count 16n; retain actual usable counts and columns,
including any sampler overshoot. This recipe does not assert equality with an
archived incumbent's support. Summands=3, batch_trials=1, max_trials=65536,
collection_window=0; all remaining worker settings are their bound defaults.

Reuse the existing source-bound build and admission runner. Freeze its build
receipt, source manifest and executable digest in the output before execution.
The local build from scientific admission has source SHA-256
215f4c7fa338b1053481563a717c867f1b75b25d936eb3b8cafa7fd2af31e134 and executable
SHA-256 fefcc1bc9de9cc6314e2a36629a0b898685f24138962553813b65c5dd1bff215.
No algorithm or observer is changed for this diagnostic. Use one Rayon worker,
60-second child deadlines and a 20-minute total driver deadline. Linux uses one
pinned CPU and an 8-GiB address-space cap; macOS records unavailable caps as null.
Reserve run numbers 2000 through 2014. Output directories must be new.

Success requires all fifteen jobs to complete and pass exact source, dispatch,
query-law, base census, matrix/rank, phase-closure and scalar checks as applicable.
A worker completion rejected by the independent checker is an admission failure,
not an incorrect final scalar by implication and never a qualified timing result.
Preserve raw reports, unsuccessful attempts, timeouts and reasons. Retain the
eleven existing environment-override rejection controls. Do not retry a failed
measured execution. A subsequent checker correction replays retained reports;
a changed producer requires a separately registered pass.

The existing runner's PASS means that a report is consistently described; an
incomplete report can pass that audit. The readiness decision additionally
requires worker_status=complete for every job. Leave performance qualification,
speedup, S and promotion unset. Native online/cold clocks are raw diagnostics
only. This local check does not establish strong-reference ranking, observer
overhead, instruction calibration or a five-cell performance qualification.
It does not consume an improvement round.
