# Pre-outcome source/interface correction on draft PR #804

Original frozen head: `5f57144b2f259fcd080daec79137bc40689020b3`.
Original freeze SHA-256: `84c2be70ad91c5f5adb57b9cf74f9ae1719a931ef844f88072200dee4927fe44`.

A deliberate direct-producer release-gate check at that head failed before the
release check: `export.py` put the parent #802 source directory at the front of
`sys.path`, causing `produce.py` to import the parent `verify.py` instead of
this PR's verifier. The raw stderr, command, exit status, and confirmation that
no output directory was created are retained here. This was not a CNF, SAT,
proof, or n131 outcome.

Pre-outcome peer review also found that the first solver-output parser expected
every CaDiCaL `v` line to end in `0`, whereas the SHA-pinned archived CaDiCaL
3.0.1 SAT transcript continues the model across lines and terminates only at
the end. Another finding required nonempty ASCII DRAT before labeling an
UNSAT proof checked. The corrected freeze appends the parent module directory,
parses continued `v` lines with strict end-of-stream and complete-assignment
checks, removes `--quiet`, requires nonempty ASCII DRAT, and pins a complete
artifact manifest for any later capped run. Hash-only CI exercises the actual
archived solver-output interface. Both runner and direct producer remain
release-blocked until #802 merges and this PR is rebased/re-frozen.
