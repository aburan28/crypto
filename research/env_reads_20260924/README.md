# Environment reads out of the Gröbner engine's inner loops

**Change.**
- `visit_macaulay_rows` and the two fused packers read the Macaulay row cap
  (`F4_F2_MAX_ROWS`, an environment lookup) once per build instead of once per
  row.
- `SolverEngine::degree_ladder` (`KIC_F4_MAX_DEGREE_ONLY`, read once per node)
  now reads its switch once per process, and so does the monomial-schedule
  cache (`KIC_F4_DISABLE_SCHEDULE_CACHE`).
- The row and column caps stay per-build reads, because tests
  (`f4_rref_tests.rs`) set them at run time.  `KIC_F4_NODE_DUMP` stays a
  per-node read, because a probe (`examples/inherited_f4_probe.rs`) sets it at
  run time.

**Why.**  A callgrind profile of `koblitz_decompose_bench 1 15 2 4` on main
(`232fceec`) spent 1.4% of its instructions in `getenv`.  There were 8,861
lookups, 6,062 of them from `max_f4_rows`, most inside the per-row loop.

**Class: engineering, on uncounted work.**  No counter moves.  This is a
stage diagnostic with no end-to-end claim (AGENTS.md §8).

**Evidence** (`run.sh`, `check.py` → `check.md`).  Main's binary and the
change's, the default policy on both, run interleaved with the arm order
rotating, three repetitions each, on the five stage ladders and the twenty
whole-logarithm seeds of
[`RESEARCH_SUPPORT_LOCAL_MULTIPLIERS.md`](../notes/ecc2k130/RESEARCH_SUPPORT_LOCAL_MULTIPLIERS.md) §4.

- **Identity.**  Every field but timings, on every rung of all 30 ladder runs,
  equals the merged default's registered run.  So do trials, relations and
  oracle word operations on every whole logarithm, and all 120 whole
  logarithms verify `[k]G = Q`.
- **Wall time**, a practicality note:
  - The ladders run `1.03–1.07×` faster.  The `0.98×` on the chain ladder is
    a `0.11 s` suite.
  - The whole logarithms are **flat**: `0.975× [0.954, 0.996]` on
    `K_0/2^13` and `1.000×` on `K_0/2^9`.  Those runs take 0.03–0.3 s, where
    process start-up and the host dominate.  Twelve more interleaved pairs on
    the heaviest seed (209) read `1.016×` in the change's favour.
  - Nothing is claimed beyond the ladders' few per cent.

Binaries: main's were copied from `target/release` built at `d69936e1`, whose
`src/` is identical to `232fceec`'s.  sha256:
- main `ic` 401caa87c3d13283…
- main `groebner_stage_bench` a915b32cb43efaf8…
- change `ic` fc351df7395b6d33…
- change `groebner_stage_bench` 72e57cda3e20e17d…

Host: Intel(R) Xeon(R) Processor @ 2.10GHz, 4 logical cores; a shared cloud host with no pinning.
