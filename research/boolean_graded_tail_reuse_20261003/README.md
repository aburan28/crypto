# Graded Boolean Macaulay batch reuse

`PROTOCOL.md` and `protocol.json` freeze an exact algebraic experiment before
implementation or timing. A fixed quadratic core makes the cubic Macaulay
block identical across changing affine tails. The proposed candidate reuses
the high-block elimination schedule, then recomputes the changed low-degree
tail and materializes every output. It must beat fresh matched packed controls
in complete cold batches to qualify.

The native Rust producer and verifier are implemented in the follow-on
branch. Status: **strongest-control discovery complete; universal gate
rejected before holdout**. [CONCLUSION.md](CONCLUSION.md) records the exact
2/8 discovery decision and n=24 finite stage gain. [ATTEMPTS.md](ATTEMPTS.md)
and [RUNS.json](RUNS.json) retain both earlier excluded raw artifacts: the
first failed post-seal float readback; the second replayed completely but
omitted an applicable 16-bit dense control at n=24. The admitted bundle is
in [qualified_discovery_01](qualified_discovery_01/manifest.json). Fresh
holdouts were not run. Full-solver cost, relation yield and rho ratio remain
null.
The new holdout seeds remain unused. Inputs are generated public Boolean
systems, with no curve or key interface. The previous full-trace screen in
PR #1202 is context, not a performance baseline.

The worker tests exhaustive small affine assignments, all five declared
sizes, both changing-affine families, repeat and support-escape guardrails,
ranked and dense coordinate maps, exact returned echelons, and native replay
rejection of altered output and contended resources. A qualified Linux ARM64
run uses `run.sh discovery NEW_DIRECTORY`. It retains all source, binary,
resource, raw-sample and verifier evidence under a manifest. The full phase
requires the verifier to accept an already sealed discovery that passed the
frozen 1.5x screen. Neither phase is launched by ordinary PR checks.
